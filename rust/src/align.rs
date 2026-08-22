//! Global alignment with an edlib-compatible CIGAR, matching
//! `edlib.align(query, target, task="path", mode="NW")`.
//!
//! `get_best_corrections` aligns every supporting segment back to the consensus
//! and consumes the **CIGAR**, not just the distance. That makes this the one
//! place where matching edlib's *choice* of alignment matters: many alignments
//! achieve the same optimal edit distance, and a different tie-break yields a
//! different-but-equally-optimal CIGAR — which changes the MSA and therefore the
//! corrected output.
//!
//! # The tie-break
//!
//! Determined empirically rather than guessed. Running the traceback with all
//! six preference orderings against 677 real `task="path"` calls recorded from
//! the reference:
//!
//! | preference (backwards from the end) | CIGARs matching |
//! | --- | --- |
//! | insert, delete, diagonal | **677 / 677** |
//! | delete, insert, diagonal | 637 / 677 |
//! | delete, diagonal, insert | 281 / 677 |
//! | insert, diagonal, delete | 131 / 677 |
//! | diagonal, insert, delete | 56 / 677 |
//! | diagonal, delete, insert | 55 / 677 |
//!
//! So edlib prefers **insertion, then deletion, then the diagonal**, walking
//! backwards from `(n, m)`. Note how badly the intuitive diagonal-first ordering
//! does — 8% — which is why this was measured.
//!
//! Every ordering reproduces the edit *distance*; only the path differs. That is
//! the whole point: the distance is unique, the alignment is not.
//!
//! CIGAR operations use edlib's extended alphabet: `=` match, `X` mismatch,
//! `I` insertion to target (consumes query), `D` deletion from target.

/// Alignment of `query` against `target`: edit distance and extended CIGAR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub edit_distance: usize,
    pub cigar: String,
}

/// Globally align `query` to `target`, returning an edlib-compatible CIGAR.
///
/// Straightforward `O(n*m)` Needleman-Wunsch with unit costs. The segments this
/// runs on are short — anchor spans are bounded by `--xmax`, 80 by default — so
/// the quadratic table is small. edlib is bit-parallel and would be faster on
/// long inputs; if that ever matters, the oracle test below is what makes a
/// faster implementation safe to swap in.
pub fn global(query: &[u8], target: &[u8]) -> Alignment {
    let (n, m) = (query.len(), target.len());

    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    let at = |i: usize, j: usize| i * (m + 1) + j;
    for i in 0..=n {
        dp[at(i, 0)] = i as u32;
    }
    for j in 0..=m {
        dp[at(0, j)] = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub = dp[at(i - 1, j - 1)] + u32::from(query[i - 1] != target[j - 1]);
            dp[at(i, j)] = sub.min(dp[at(i - 1, j)] + 1).min(dp[at(i, j - 1)] + 1);
        }
    }

    // Traceback, preferring insertion, then deletion, then the diagonal.
    let (mut i, mut j) = (n, m);
    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    while i > 0 || j > 0 {
        if i > 0 && dp[at(i - 1, j)] + 1 == dp[at(i, j)] {
            ops.push(b'I');
            i -= 1;
        } else if j > 0 && dp[at(i, j - 1)] + 1 == dp[at(i, j)] {
            ops.push(b'D');
            j -= 1;
        } else if i > 0 && j > 0 {
            let same = query[i - 1] == target[j - 1];
            ops.push(if same { b'=' } else { b'X' });
            i -= 1;
            j -= 1;
        } else {
            break;
        }
    }
    ops.reverse();

    Alignment {
        edit_distance: dp[at(n, m)] as usize,
        cigar: run_length_encode(&ops),
    }
}

/// One run of a CIGAR: `(length, operation)`.
pub type CigarOp = (usize, u8);

/// Parse an extended CIGAR into runs, matching `help_functions.cigar_to_seq`'s
/// `re.split(r'[=DXSMI]+', cigar)` walk.
///
/// The reference's regex tolerates `S` and `M` when splitting but its loop then
/// falls through to `print("error"); sys.exit()` for anything that is not
/// `=`, `X`, `I` or `D`. Nothing in the pipeline produces those operations —
/// edlib emits only the extended four — so this returns `None` instead of
/// aborting the process, and callers treat it as a bug rather than an input.
pub fn parse_cigar(cigar: &str) -> Option<Vec<CigarOp>> {
    let mut ops = Vec::new();
    let mut len = 0usize;
    let mut seen_digit = false;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            len = len * 10 + usize::from(c - b'0');
            seen_digit = true;
        } else {
            if !seen_digit || !matches!(c, b'=' | b'X' | b'I' | b'D') {
                return None;
            }
            ops.push((len, c));
            len = 0;
            seen_digit = false;
        }
    }
    if seen_digit {
        return None; // trailing length with no operation
    }
    Some(ops)
}

/// Expand a CIGAR into the two gapped strings, matching
/// `help_functions.cigar_to_seq(cigar, query, ref)`.
///
/// Returns `(query_aligned, ref_aligned)` — the reference's return order, which
/// `get_best_corrections` immediately unpacks as `read_alignment,
/// ref_alignment`. `I` is an insertion with respect to the reference (a gap in
/// `ref_aligned`); `D` is a deletion (a gap in `query_aligned`).
///
/// Returns `None` if the CIGAR does not parse or does not cover both sequences.
pub fn cigar_to_seq(cigar: &str, query: &[u8], reference: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    ops_to_seq(&parse_cigar(cigar)?, query, reference)
}

/// [`cigar_to_seq`] without the string.
///
/// The SIMD aligner reports operations directly, so on the hot path formatting a
/// CIGAR only to parse it straight back is pure overhead — 4.27 million times per
/// four real clusters. [`crate::simd::global_ops`] feeds this instead.
pub fn ops_to_seq(ops: &[CigarOp], query: &[u8], reference: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let (mut q_index, mut r_index) = (0usize, 0usize);
    let (mut q_aln, mut r_aln) = (Vec::new(), Vec::new());

    for &(len, op) in ops {
        match op {
            b'=' | b'X' => {
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_index += len;
                r_index += len;
            }
            b'I' => {
                r_aln.extend(std::iter::repeat_n(b'-', len));
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                q_index += len;
            }
            b'D' => {
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_aln.extend(std::iter::repeat_n(b'-', len));
                r_index += len;
            }
            _ => return None,
        }
    }
    // Python would silently accept a CIGAR that stops short; nothing produces
    // one, and a partial alignment would corrupt the MSA rather than fail.
    if q_index != query.len() || r_index != reference.len() {
        return None;
    }
    Some((q_aln, r_aln))
}

/// `[(2, b'='), (1, b'X')]` becomes `"2=1X"`. The inverse of [`parse_cigar`].
pub fn encode_cigar(ops: &[CigarOp]) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    for &(len, op) in ops {
        let _ = write!(out, "{len}{}", op as char);
    }
    out
}

/// `[b'=', b'=', b'X']` becomes `"2=1X"`.
fn run_length_encode(ops: &[u8]) -> String {
    let mut out = String::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(q: &str, t: &str) -> String {
        global(q.as_bytes(), t.as_bytes()).cigar
    }

    #[test]
    fn identical_sequences_are_all_matches() {
        let a = global(b"ACGT", b"ACGT");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "4=");
    }

    #[test]
    fn a_substitution_shows_as_x() {
        let a = global(b"ACGT", b"AGGT");
        assert_eq!(a.edit_distance, 1);
        assert_eq!(a.cigar, "1=1X2=");
    }

    #[test]
    fn empty_query_is_all_deletions() {
        let a = global(b"", b"ACGT");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4D");
    }

    #[test]
    fn empty_target_is_all_insertions() {
        let a = global(b"ACGT", b"");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4I");
    }

    #[test]
    fn both_empty_yields_an_empty_cigar() {
        let a = global(b"", b"");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "");
    }

    #[test]
    fn cigar_operation_counts_cover_both_sequences() {
        // Query length = = + X + I; target length = = + X + D.
        let (q, t) = ("ACGTACGTAA", "ACGTTACGT");
        let c = cigar(q, t);
        let (mut qn, mut tn, mut num) = (0usize, 0usize, String::new());
        for ch in c.chars() {
            if ch.is_ascii_digit() {
                num.push(ch);
            } else {
                let n: usize = num.parse().unwrap();
                num.clear();
                match ch {
                    '=' | 'X' => {
                        qn += n;
                        tn += n;
                    }
                    'I' => qn += n,
                    'D' => tn += n,
                    _ => panic!("unexpected op {ch}"),
                }
            }
        }
        assert_eq!(qn, q.len(), "query coverage");
        assert_eq!(tn, t.len(), "target coverage");
    }

    fn expand(cigar: &str, q: &str, r: &str) -> (String, String) {
        let (qa, ra) = cigar_to_seq(cigar, q.as_bytes(), r.as_bytes()).unwrap();
        (
            String::from_utf8(qa).unwrap(),
            String::from_utf8(ra).unwrap(),
        )
    }

    #[test]
    fn cigar_parses_into_runs() {
        assert_eq!(
            parse_cigar("2=1X10I").unwrap(),
            vec![(2, b'='), (1, b'X'), (10, b'I')]
        );
        assert_eq!(parse_cigar("").unwrap(), vec![]);
        assert!(parse_cigar("3M").is_none(), "M is not in the pipeline");
        assert!(parse_cigar("3").is_none(), "length with no operation");
        assert!(parse_cigar("=").is_none(), "operation with no length");
    }

    #[test]
    fn insertion_gaps_the_reference_and_deletion_gaps_the_query() {
        let (q, r) = expand("2=2I2=", "ACGTAC", "ACAC");
        assert_eq!((q.as_str(), r.as_str()), ("ACGTAC", "AC--AC"));
        let (q, r) = expand("2=2D2=", "ACAC", "ACGTAC");
        assert_eq!((q.as_str(), r.as_str()), ("AC--AC", "ACGTAC"));
    }

    #[test]
    fn an_expanded_cigar_round_trips_its_own_alignment() {
        let (q, r) = ("ACGTACGTAA", "ACGTTACGT");
        let a = global(q.as_bytes(), r.as_bytes());
        let (qa, ra) = expand(&a.cigar, q, r);
        assert_eq!(qa.len(), ra.len());
        assert_eq!(qa.replace('-', ""), q);
        assert_eq!(ra.replace('-', ""), r);
        // Gapped columns account for exactly the edit distance.
        let diffs = qa.bytes().zip(ra.bytes()).filter(|(a, b)| a != b).count();
        assert_eq!(diffs, a.edit_distance);
    }

    #[test]
    fn a_cigar_that_does_not_cover_both_sequences_is_rejected() {
        assert!(cigar_to_seq("2=", b"ACGT", b"ACGT").is_none());
        assert!(cigar_to_seq("4=", b"ACGT", b"ACG").is_none());
    }

    /// Tie cases verified directly against Python edlib, not assumed.
    ///
    /// The preference is applied walking *backwards*, so the gap lands at the
    /// end of the emitted CIGAR rather than the start --- which is the opposite
    /// of what "insertion first" suggests when read forwards.
    #[test]
    fn tie_cases_match_edlib_ground_truth() {
        assert_eq!(cigar("AA", "A"), "1=1I");
        assert_eq!(cigar("A", "AA"), "1=1D");
        assert_eq!(cigar("ACGT", "ACT"), "2=1I1=");
        assert_eq!(cigar("AAA", "A"), "1=2I");
    }
}

/// Differential oracle against Python's edlib.
///
/// ```bash
/// bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/cg
/// CIGAR_CASES=/tmp/cg/cigar_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml align::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    #[test]
    fn matches_edlib_cigars_exactly() {
        let Ok(path) = std::env::var("CIGAR_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CIGAR_CASES={path} unreadable: {e}"));

        let (mut n, mut bad_cigar, mut bad_ed, mut shown) = (0usize, 0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                continue;
            }
            let got = global(f[0].as_bytes(), f[1].as_bytes());
            let want_ed: usize = f[2].parse().unwrap();
            if got.edit_distance != want_ed {
                bad_ed += 1;
            }
            if got.cigar != f[3] {
                if shown < 3 {
                    shown += 1;
                    eprintln!("cigar mismatch:\n  edlib: {}\n  rust : {}", f[3], got.cigar);
                }
                bad_cigar += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad_ed, 0, "{bad_ed} of {n} edit distances differed");
        assert_eq!(bad_cigar, 0, "{bad_cigar} of {n} CIGARs differed");
        eprintln!("matched edlib on all {n} recorded path alignments");
    }

    /// Replay `help_functions.cigar_to_seq`: `#cigar query ref q_aln r_aln`.
    ///
    /// Separate from the CIGAR oracle above because the two can fail
    /// independently — expansion is deterministic given a CIGAR, so a mismatch
    /// here is this function's bug and not the aligner's.
    ///
    /// ```bash
    /// CIGAR_TO_SEQ_CASES=/tmp/cg/cigar_to_seq_cases.tsv \
    ///   cargo test --manifest-path rust/Cargo.toml align::oracle -- --nocapture
    /// ```
    #[test]
    fn expands_cigars_like_the_reference() {
        let Ok(path) = std::env::var("CIGAR_TO_SEQ_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CIGAR_TO_SEQ_CASES={path} unreadable: {e}"));

        let (mut n, mut bad) = (0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 5, "short row: {line}");
            let (q_aln, r_aln) = cigar_to_seq(f[0], f[1].as_bytes(), f[2].as_bytes())
                .unwrap_or_else(|| panic!("rejected a CIGAR the reference expanded: {}", f[0]));
            if q_aln != f[3].as_bytes() || r_aln != f[4].as_bytes() {
                if bad < 3 {
                    eprintln!(
                        "mismatch on {}:\n  python: {} / {}\n  rust  : {} / {}",
                        f[0],
                        f[3],
                        f[4],
                        String::from_utf8_lossy(&q_aln),
                        String::from_utf8_lossy(&r_aln)
                    );
                }
                bad += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad, 0, "{bad} of {n} expansions differed");
        eprintln!("expanded all {n} recorded CIGARs identically");
    }
}

/// Could `block-aligner` replace this? A measurement.
///
/// After block-aligner took the guard from 65% of runtime to 0.9%, this became
/// the largest single cost — 31.4% over ten real clusters, ~61 000 alignments
/// per cluster. The obvious move is to reuse the same SIMD aligner here.
///
/// Two things make it a real question rather than a formality.
///
/// **Length.** These are anchor-bounded segments against a consensus: *median
/// 55 x 57 bp, max 89 x 94*, where block-aligner's minimum block size is 32. It
/// won 99x on 1.4 kb guard alignments but only 8x on short ones, and it pays
/// per-call setup for `PaddedBytes` across tens of thousands of calls.
///
/// **Objective.** block-aligner asserts `gaps.open < gaps.extend`, so it cannot
/// express unit-cost edit distance at all. That is not the obstacle it first
/// appears: edit distance was the reference's choice *for speed in Python*, not
/// a modelling preference, and affine gaps describe indels better. So the
/// comparison worth making is against the guard's own scoring — `match 4,
/// mismatch -8, open 12, extend 1` — which block-aligner does express, with
/// [`crate::parasail::global_affine`] as the exact oracle for it.
///
/// The unit-cost CIGAR count below is therefore expected to be low. It measures
/// how far affine moves the answer, not an error.
///
/// ```bash
/// CIGAR_CASES=/tmp/d/cigar_cases.tsv \
///   cargo test --release --manifest-path rust/Cargo.toml align::block_option -- --nocapture
/// ```
#[cfg(test)]
mod block_option {
    use crate::parasail::{self, Scoring, TieBreak};

    /// `(query, target, edlib's CIGAR, edlib's edit distance)`.
    type Case = (Vec<u8>, Vec<u8>, String, usize);

    fn cases() -> Option<Vec<Case>> {
        let path = std::env::var("CIGAR_CASES").ok()?;
        let data = std::fs::read_to_string(path).ok()?;
        Some(
            data.lines()
                .filter(|l| !l.starts_with('#') && !l.is_empty())
                .filter_map(|l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    if f.len() < 4 {
                        return None;
                    }
                    // block-aligner asserts every character is A-Z, and the
                    // recorded corpus includes `min_ed`'s calls, whose query is
                    // a gap-padded insertion slot like "-TT-". Those would have
                    // to be handled separately; skip them here so the timing
                    // reflects the main segment-vs-consensus workload.
                    if !f[0]
                        .bytes()
                        .chain(f[1].bytes())
                        .all(|c| c.is_ascii_uppercase())
                    {
                        return None;
                    }
                    Some((
                        f[0].as_bytes().to_vec(),
                        f[1].as_bytes().to_vec(),
                        f[3].to_string(),
                        f[2].parse().ok()?,
                    ))
                })
                .collect(),
        )
    }

    #[test]
    fn speed_and_optimality_at_segment_lengths() {
        let Some(cases) = cases() else { return };
        if cases.is_empty() {
            return;
        }

        let sc = Scoring::GUARD;

        let t0 = std::time::Instant::now();
        for (q, r, _, _) in &cases {
            std::hint::black_box(super::global(q, r));
        }
        let scalar_ms = t0.elapsed().as_secs_f64() * 1000.0;

        // The shipped aligner, not a copy of it: `simd::global_ops` is exactly
        // what `corrections::segment_ops` calls, including the reused block and
        // the derived starting width.
        let mut ops = Vec::new();
        let t0 = std::time::Instant::now();
        let mut refused = 0usize;
        let mut same_as_unit_cost = 0usize;
        let mut got: Vec<Option<Vec<crate::align::CigarOp>>> = Vec::with_capacity(cases.len());
        for (q, r, want_cigar, _) in &cases {
            if !crate::simd::global_ops(q, r, sc, &mut ops) {
                refused += 1;
                got.push(None);
                continue;
            }
            if crate::align::encode_cigar(&ops) == *want_cigar {
                same_as_unit_cost += 1;
            }
            got.push(Some(ops.clone()));
        }
        let block_ms = t0.elapsed().as_secs_f64() * 1000.0;

        // Optimality against the exact affine DP. Not timed — it is the oracle.
        let mut suboptimal = 0usize;
        let mut unexpanded = 0usize;
        for ((q, r, _, _), ops) in cases.iter().zip(&got) {
            let Some(ops) = ops else { continue };
            let exact = parasail::global_affine(q, r, sc, TieBreak::PARASAIL);
            if score_of(ops, q, r, sc) != exact.score {
                suboptimal += 1;
            }
            if crate::align::ops_to_seq(ops, q, r).is_none() {
                unexpanded += 1;
            }
        }

        eprintln!(
            "{} recorded segment-vs-consensus alignments (median ~55x57 bp)",
            cases.len()
        );
        eprintln!("  scalar unit-cost DP     {scalar_ms:>9.1} ms");
        eprintln!(
            "  simd::global_ops        {block_ms:>9.1} ms   ({:.2}x)",
            scalar_ms / block_ms.max(0.001)
        );
        eprintln!(
            "  suboptimal {suboptimal}, refused {refused}, unexpandable {unexpanded}, \
             same path as unit cost {same_as_unit_cost}/{}",
            cases.len()
        );
        assert_eq!(suboptimal, 0, "missed the optimal affine score");
        assert_eq!(
            unexpanded, 0,
            "returned a CIGAR that does not span its inputs"
        );
        // Every case here is ACGTN by construction, so nothing should be refused.
        assert_eq!(refused, 0, "refused a plain ACGTN pair");
    }

    /// Score an alignment under `sc`, to compare against the exact DP's own.
    fn score_of(ops: &[crate::align::CigarOp], q: &[u8], r: &[u8], sc: Scoring) -> i32 {
        let (mut qi, mut ri, mut total) = (0usize, 0usize, 0i32);
        for &(len, op) in ops {
            match op {
                b'=' | b'X' => {
                    for _ in 0..len {
                        total += if q[qi] == r[ri] {
                            sc.match_score
                        } else {
                            sc.mismatch
                        };
                        qi += 1;
                        ri += 1;
                    }
                }
                b'I' | b'D' => {
                    total -= sc.open + sc.ext * (len as i32 - 1);
                    if op == b'I' {
                        qi += len;
                    } else {
                        ri += len;
                    }
                }
                _ => unreachable!("unexpected op"),
            }
        }
        total
    }
}
