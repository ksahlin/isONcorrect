//! The structural-overcorrection guard and the per-read driver, matching
//! `fix_correction` and `correct_read`.
//!
//! This is the last thing that happens to a read, and it is easy to miss: after
//! the corrected intervals have been stitched back in, the result is aligned
//! against the **original** read and any indel run longer than 10 alignment
//! columns is **reverted to the original**. The assumption is that a difference
//! that large is real structure — exon variation between isoforms — rather than
//! a sequencing error, and correcting it would be flattening a genuine
//! transcript difference into the consensus. Runs of ten or fewer keep the
//! correction.
//!
//! # What decides the answer
//!
//! * the threshold is on **alignment columns**, not on inserted bases, and it is
//!   `> 10`, so a run of exactly 10 keeps the correction;
//! * a run mixes both directions. `l` counts columns of either kind while both
//!   `o_segm` and `c_segm` accumulate, so a deletion immediately followed by an
//!   insertion is **one** run. When it flushes, only one side survives: the
//!   original's bases if the run is long, the corrected read's if it is short.
//!   The other side's characters are dropped;
//! * the trailing run is flushed after the loop, so an alignment ending in an
//!   indel is handled — the reference repeats the same four lines to do it;
//! * `correct_read` aligns `(seq, corr)` in that order, so `orig` is the first
//!   row and the guard reverts *towards* it.

use crate::align;
use crate::parasail::{self, Scoring};
use crate::simd;

/// Align the read against its correction, for the guard.
///
/// **Default: [`simd::global_cigar`], which is SIMD and adaptively banded.** It finds
/// an optimal-scoring alignment — verified on 1 400 of 1 400 recorded
/// alignments, zero suboptimal — but reports a *different* equally-optimal path
/// than parasail in most cases, which shifts the guard's answer on roughly 1% of
/// reads. That is the one place the port deliberately diverges from the
/// reference; see PORTING.md.
///
/// Setting `ISONCORRECT_EXACT_GUARD=1` restores the exact parasail-compatible
/// DP, which is byte-identical to the reference. `bench/equivalence.sh` runs in
/// that mode, so the equivalence gate still covers the whole pipeline.
fn align(seq: &[u8], corr: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let _t = crate::profile::scope("guard alignment");
    let use_exact = exact_guard() || seq.len() < simd::MIN_LEN || corr.len() < simd::MIN_LEN;
    let cigar = if use_exact {
        parasail::semiglobal(seq, corr, Scoring::GUARD).cigar
    } else {
        match simd::global_cigar(seq, corr, Scoring::GUARD) {
            Some(cigar) => cigar,
            // Falling back rather than failing: an alignment that does not span
            // both sequences cannot be expanded, and the exact DP always can.
            None => parasail::semiglobal(seq, corr, Scoring::GUARD).cigar,
        }
    };
    align::cigar_to_seq(&cigar, seq, corr).expect("the aligner's own CIGAR expands over its inputs")
}

/// Whether `ISONCORRECT_EXACT_GUARD` asks for the parasail-compatible DP.
fn exact_guard() -> bool {
    use std::sync::atomic::{AtomicU8, Ordering};
    static MODE: AtomicU8 = AtomicU8::new(u8::MAX);
    let cached = MODE.load(Ordering::Relaxed);
    if cached != u8::MAX {
        return cached == 1;
    }
    let on = std::env::var("ISONCORRECT_EXACT_GUARD")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false);
    MODE.store(u8::from(on), Ordering::Relaxed);
    on
}

/// Longest indel run, in alignment columns, that keeps its correction.
pub const MAX_INDEL_RUN: usize = 10;

/// Revert over-long indel runs to the original read, matching
/// `fix_correction(orig, corr)`.
///
/// `orig` and `corr` are the two rows of an alignment, equal length, gaps as
/// `-`. Returns the adjusted corrected sequence, ungapped.
pub fn fix_correction(orig: &[u8], corr: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(corr.len());
    let mut o_segm: Vec<u8> = Vec::new();
    let mut c_segm: Vec<u8> = Vec::new();
    let mut run = 0usize;

    let flush = |out: &mut Vec<u8>, o_segm: &mut Vec<u8>, c_segm: &mut Vec<u8>, run: usize| {
        if run > MAX_INDEL_RUN {
            out.extend_from_slice(o_segm);
        } else if run > 0 {
            out.extend_from_slice(c_segm);
        }
        o_segm.clear();
        c_segm.clear();
    };

    for (&o, &c) in orig.iter().zip(corr) {
        if o != b'-' && c != b'-' {
            flush(&mut out, &mut o_segm, &mut c_segm, run);
            out.push(c);
            run = 0;
        } else if o == b'-' {
            c_segm.push(c);
            run += 1;
        } else {
            // Reached only when `o != '-'` and `c == '-'`. The reference has a
            // fourth branch here that raises, but it is dead: the `elif o ==
            // '-'` above already absorbs the gap-aligned-to-gap case.
            o_segm.push(o);
            run += 1;
        }
    }
    flush(&mut out, &mut o_segm, &mut c_segm, run);
    out
}

/// One interval to correct: the span it replaces in the read, and what to
/// replace it with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub stop: usize,
    pub corrected: Vec<u8>,
}

/// Splice the corrected intervals into the read, matching `correct_read`'s
/// stitching.
///
/// Intervals arrive in ascending order and do not overlap — `solve_WIS`
/// guarantees it, and `get_intervals_to_correct` reverses the indices to restore
/// ascending order. Everything between them is carried over from the original.
///
/// # Panics
///
/// On an empty `intervals`, which the reference also cannot handle
/// (`corr_seq[0][0]`). It never happens: `isoncorrect_main` skips
/// `correct_read` entirely when WIS selects nothing.
pub fn stitch(seq: &[u8], intervals: &[Interval]) -> Vec<u8> {
    assert!(!intervals.is_empty(), "correct_read with no intervals");
    let mut out = Vec::with_capacity(seq.len());
    out.extend_from_slice(&seq[..intervals[0].start.min(seq.len())]);
    for (cnt, iv) in intervals.iter().enumerate() {
        out.extend_from_slice(&iv.corrected);
        let from = iv.stop.min(seq.len());
        let to = match intervals.get(cnt + 1) {
            Some(next) => next.start.min(seq.len()).max(from),
            None => seq.len(),
        };
        out.extend_from_slice(&seq[from..to]);
    }
    out
}

/// Stitch, then guard: the whole tail of `correct_read`.
///
/// The reference realigns once more after `fix_correction` and throws the result
/// away; that call is not reproduced (see *Deferred improvements*).
pub fn correct_read(seq: &[u8], intervals: &[Interval]) -> Vec<u8> {
    let corr = stitch(seq, intervals);
    let (orig_aln, corr_aln) = align(seq, &corr);
    let _t = crate::profile::scope("fix_correction");
    fix_correction(&orig_aln, &corr_aln)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fix(orig: &str, corr: &str) -> String {
        String::from_utf8(fix_correction(orig.as_bytes(), corr.as_bytes())).unwrap()
    }

    #[test]
    fn a_gapless_alignment_returns_the_correction() {
        assert_eq!(fix("ACGTACGT", "ACGTTCGT"), "ACGTTCGT");
    }

    #[test]
    fn a_short_insertion_keeps_the_correction() {
        // Three inserted bases: a run of 3, well under the limit.
        assert_eq!(fix("ACGT---ACGT", "ACGTTTTACGT"), "ACGTTTTACGT");
    }

    #[test]
    fn a_short_deletion_keeps_the_correction() {
        assert_eq!(fix("ACGTTTTACGT", "ACGT---ACGT"), "ACGTACGT");
    }

    #[test]
    fn a_long_insertion_is_reverted_to_the_original() {
        // 11 inserted bases: over the limit, so the original's bases win --- and
        // the original has none there, so the insertion vanishes.
        let orig = format!("ACGT{}ACGT", "-".repeat(11));
        let corr = format!("ACGT{}ACGT", "T".repeat(11));
        assert_eq!(fix(&orig, &corr), "ACGTACGT");
    }

    #[test]
    fn a_long_deletion_is_reverted_to_the_original() {
        let orig = format!("ACGT{}ACGT", "T".repeat(11));
        let corr = format!("ACGT{}ACGT", "-".repeat(11));
        assert_eq!(fix(&orig, &corr), format!("ACGT{}ACGT", "T".repeat(11)));
    }

    #[test]
    fn a_run_of_exactly_ten_keeps_the_correction() {
        let orig = format!("ACGT{}ACGT", "-".repeat(10));
        let corr = format!("ACGT{}ACGT", "T".repeat(10));
        assert_eq!(fix(&orig, &corr), format!("ACGT{}ACGT", "T".repeat(10)));
        // Eleven is the first reverted length.
        let orig = format!("ACGT{}ACGT", "-".repeat(11));
        let corr = format!("ACGT{}ACGT", "T".repeat(11));
        assert_eq!(fix(&orig, &corr), "ACGTACGT");
    }

    #[test]
    fn a_trailing_indel_is_flushed_after_the_loop() {
        assert_eq!(fix("ACGT---", "ACGTTTT"), "ACGTTTT");
        let orig = format!("ACGT{}", "-".repeat(11));
        let corr = format!("ACGT{}", "T".repeat(11));
        assert_eq!(fix(&orig, &corr), "ACGT");
    }

    #[test]
    fn adjacent_deletion_and_insertion_count_as_one_run() {
        // 6 deleted then 6 inserted columns is a run of 12, so the original
        // side wins and the inserted bases are dropped entirely.
        let orig = "ACGTTTTTTT------ACGT";
        let corr = "ACGT------GGGGGGACGT";
        assert_eq!(fix(orig, corr), "ACGTTTTTTTACGT");
        // One matching column between them makes two runs of 6, each under the
        // limit, so each keeps the corrected side --- which for the deletion is
        // empty, so those bases go, and for the insertion is the inserted bases,
        // so those stay.
        let orig = "ACGTTTTTTTA------ACGT";
        let corr = "ACGT------AGGGGGGACGT";
        assert_eq!(fix(orig, corr), "ACGTAGGGGGGACGT");
    }

    fn iv(start: usize, stop: usize, corrected: &str) -> Interval {
        Interval {
            start,
            stop,
            corrected: corrected.as_bytes().to_vec(),
        }
    }

    #[test]
    fn stitching_replaces_each_span_and_keeps_the_rest() {
        let seq = b"AAAACCCCGGGGTTTT";
        let out = stitch(seq, &[iv(4, 8, "cccc"), iv(12, 16, "tttt")]);
        assert_eq!(out, b"AAAAccccGGGGtttt");
    }

    #[test]
    fn stitching_keeps_the_prefix_and_the_tail() {
        let seq = b"AAAACCCCGGGGTTTT";
        assert_eq!(stitch(seq, &[iv(4, 8, "xx")]), b"AAAAxxGGGGTTTT");
        assert_eq!(stitch(seq, &[iv(0, 4, "xx")]), b"xxCCCCGGGGTTTT");
        assert_eq!(stitch(seq, &[iv(12, 16, "xx")]), b"AAAACCCCGGGGxx");
    }

    #[test]
    fn a_correction_may_change_the_span_length() {
        let seq = b"AAAACCCCGGGG";
        assert_eq!(stitch(seq, &[iv(4, 8, "")]), b"AAAAGGGG");
        assert_eq!(stitch(seq, &[iv(4, 8, "cccccccc")]), b"AAAAccccccccGGGG");
    }

    #[test]
    fn an_unchanged_read_survives_the_whole_guard() {
        let seq = b"ACGTACGTTGCAACGTACGTTGCAACGT";
        // One interval whose correction is exactly the original span.
        let out = correct_read(seq, &[iv(8, 16, "TGCAACGT")]);
        assert_eq!(out, seq);
    }

    #[test]
    fn the_guard_reverts_a_long_insertion_end_to_end() {
        let seq = b"ACGTACGTTGCAACGTACGTTGCAACGTAAGGCCTTAACC";
        // Correction that inserts 15 bases into the middle of its span.
        let mut corrected = b"TGCAACGT".to_vec();
        corrected.extend_from_slice(b"GGGGGGGGGGGGGGG");
        let out = correct_read(
            seq,
            &[Interval {
                start: 8,
                stop: 16,
                corrected,
            }],
        );
        assert_eq!(out, seq, "a 15-base insertion should have been reverted");
    }
}

/// Differential oracles for the guard.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// FIX_CORRECTION_CASES=/tmp/d/fix_correction_cases.tsv \
///   CORRECT_READ_CASES=/tmp/d/correct_read_cases.tsv \
///   CORRECTION_CASES=/tmp/d/correction_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml guard::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    /// `fix_correction` is a pure function of the two aligned rows, so this
    /// replays it without parasail in the loop — a mismatch here is this
    /// function's bug, not the aligner's.
    #[test]
    fn fix_correction_matches_the_reference() {
        let Ok(path) = std::env::var("FIX_CORRECTION_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("FIX_CORRECTION_CASES={path} unreadable: {e}"));

        let (mut n, mut bad, mut shown) = (0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 3, "short row");
            let got = fix_correction(f[0].as_bytes(), f[1].as_bytes());
            if got != f[2].as_bytes() {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "fix_correction differs:\n  python: {}\n  rust  : {}",
                        f[2],
                        String::from_utf8_lossy(&got)
                    );
                }
                bad += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad, 0, "{bad} of {n} adjusted sequences differed");
        eprintln!("fix_correction matched the reference on all {n} calls");
    }

    /// `correct_read` end to end: stitching plus the guard.
    ///
    /// Intervals are recorded by source, so the corrected span of each comes
    /// either from a literal (a region carried over in
    /// `previously_corrected_regions`) or from the numbered
    /// `get_best_corrections` call in `correction_cases.tsv`. Reading it back
    /// that way keeps the harness from duplicating any reference logic.
    #[test]
    fn correct_read_matches_the_reference() {
        let (Ok(reads_path), Ok(corr_path)) = (
            std::env::var("CORRECT_READ_CASES"),
            std::env::var("CORRECTION_CASES"),
        ) else {
            return;
        };

        // call number -> the corrected span it returned for the current read
        let corr_data = std::fs::read_to_string(&corr_path).expect("correction cases");
        let mut by_call: Vec<String> = Vec::new();
        for line in corr_data.lines().filter(|l| l.starts_with("C\t")) {
            let f: Vec<&str> = line.split('\t').collect();
            by_call.push(f.get(5).unwrap_or(&"").to_string());
        }

        let data = std::fs::read_to_string(&reads_path).expect("correct_read cases");
        let (mut n, mut bad, mut shown) = (0usize, 0usize, 0usize);
        let mut seq = Vec::new();
        let mut intervals: Vec<Interval> = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "R" => {
                    seq = f[2].as_bytes().to_vec();
                    intervals.clear();
                }
                "I" => {
                    let corrected = match f[5].split_once(':') {
                        Some(("s", literal)) => literal.as_bytes().to_vec(),
                        Some(("c", call)) => {
                            by_call[call.parse::<usize>().unwrap()].as_bytes().to_vec()
                        }
                        other => panic!("bad interval source {other:?}"),
                    };
                    intervals.push(Interval {
                        start: f[3].parse().unwrap(),
                        stop: f[4].parse().unwrap(),
                        corrected,
                    });
                }
                "A" => {
                    let got = correct_read(&seq, &intervals);
                    if got != f[2].as_bytes() {
                        if shown < 3 {
                            shown += 1;
                            eprintln!(
                                "correct_read differs on a {}bp read:\n  python: {}\n  rust  : {}",
                                seq.len(),
                                f[2],
                                String::from_utf8_lossy(&got)
                            );
                        }
                        bad += 1;
                    }
                    n += 1;
                }
                other => panic!("unknown record {other}"),
            }
        }
        assert!(n > 0, "no cases in {reads_path}");
        assert_eq!(bad, 0, "{bad} of {n} corrected reads differed");
        eprintln!("correct_read matched the reference on all {n} reads");
    }
}
