//! Correcting one interval, matching `get_best_corrections`.
//!
//! Everything the earlier stages built comes together here: a consensus over
//! the supporting segments, each segment aligned back to it, those pairwise
//! alignments stacked into a matrix, and then — column by column — a decision
//! about what each read *should* have said.
//!
//! # The two decisions
//!
//! Each column takes one of two paths, and which one depends on whether the
//! column has [alternative references](crate::contexts):
//!
//! * **with alternatives** — the read is corrected towards whichever reference
//!   its own context is closest to. Its context is compared against the
//!   consensus and against every alternative; an exact match wins immediately,
//!   otherwise the strictly smallest edit distance does, and the consensus wins
//!   ties by being measured first. This is what keeps a real variant from being
//!   flattened into the consensus.
//! * **without alternatives** — the read is corrected to the consensus, unless
//!   its own base is common enough (`variant_threshold <= PFM[i][base]`) *and*
//!   the rest of its window matches the consensus exactly. That last condition
//!   is the point: a base that differs where everything around it agrees looks
//!   like a variant, whereas one surrounded by other differences looks like
//!   noise.
//!
//! # Details that decide the answer
//!
//! * **a consensus gap with no alternatives emits nothing at all.** The whole
//!   per-read loop sits under `if ref_nucl != "-"`, so an insertion relative to
//!   the consensus contributes no column — which is how insertions get dropped.
//!   Output-equivalent to emitting `-` and filtering it, but only because every
//!   row skips together;
//! * the two `k_size` trims at the end: the segment was taken from anchor to
//!   anchor *inclusive*, and both anchors are cut back off, so what is returned
//!   is strictly the span between them;
//! * rows are processed in partition order, and the returned per-read map keeps
//!   it. Nothing downstream depends on that today, but the reference's order is
//!   free to preserve here.
//!
//! # The alternatives loop used to be non-deterministic
//!
//! It iterates the alternatives with a `break` on an exact match and a strict
//! `<`, so **the order they arrive in decides the answer** whenever two of them
//! tie on edit distance. The reference supplied them as a Python `set`, whose
//! iteration order depends on `PYTHONHASHSEED` — measured, not deduced: on one
//! interval of a gene-level SIRV cluster, 51 of 267 corrected regions took
//! between two and five different values across 24 hash seeds.
//!
//! The reference now returns a list in insertion order, so there is a single
//! defined answer and this port reproduces it. `bench/check_seed_sensitivity.py`
//! is the regression check; PORTING.md has the measurements.

use indexmap::IndexMap;

use crate::align;
use crate::contexts;
use crate::editdist;
use crate::msa;
use crate::parasail::Scoring;
use crate::poa;
use crate::simd;

/// A corrected region of one read: `(pos1 + k_size, pos2, weight, sequence)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Correction {
    pub start: usize,
    pub stop: usize,
    /// Supporting segment count — the reference's `weight`.
    pub weight: usize,
    pub seq: Vec<u8>,
}

/// What `get_best_corrections` returns: the current read's corrected span, and
/// the same for every other read that supported the interval.
#[derive(Debug, Clone)]
pub struct Corrections {
    /// The corrected span of the read being corrected.
    pub curr_read: Vec<u8>,
    /// Corrected regions by `r_id`, in partition order.
    pub others: IndexMap<usize, Vec<Correction>>,
}

/// Correct one interval, matching
/// `get_best_corrections(curr_best_seqs, reads, k_size, ...)`.
///
/// `instance` is the interval's supporting segments as `(r_id, pos1, pos2)`,
/// **with the read being corrected first** — the reference reads
/// `curr_best_seqs[0]` for exactly that. `seqs` is indexed by `r_id - 1`.
///
/// Returns `None` only for an empty instance, which the reference never
/// produces.
pub fn best_corrections(
    instance: &[(usize, usize, usize)],
    seqs: &[Vec<u8>],
    k_size: usize,
    v_depth_ratio_threshold: f64,
    max_seqs_to_spoa: usize,
) -> Option<Corrections> {
    let (curr_read_id, _, _) = *instance.first()?;
    let weight = instance.len();

    // `reads[q_id][1][pos1 : pos2 + k_size]` --- Python slices clamp.
    let segments: Vec<&[u8]> = instance
        .iter()
        .map(|&(q_id, pos1, pos2)| {
            let seq = &seqs[q_id - 1];
            &seq[pos1.min(seq.len())..(pos2 + k_size).min(seq.len())]
        })
        .collect();

    let consensus = {
        let _t = crate::profile::scope("run_spoa");
        poa::consensus(&poa::apply_spoa_cutoff(&segments, max_seqs_to_spoa))?
    };
    let consensus = consensus.into_bytes();

    // The partition: the consensus aligned to itself, then every segment.
    let mut aligned: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(segments.len() + 1);
    aligned.push((consensus.clone(), consensus.clone()));
    let _align_t = crate::profile::scope("align (NW cigar)");
    let mut ops: Vec<align::CigarOp> = Vec::new();
    for seg in &segments {
        segment_ops(seg, &consensus, &mut ops);
        let (read_aln, ref_aln) =
            align::ops_to_seq(&ops, seg, &consensus).expect("an alignment's own CIGAR expands");
        aligned.push((read_aln, ref_aln));
    }

    let rows: Vec<(&[u8], &[u8])> = aligned
        .iter()
        .map(|(q, t)| (q.as_slice(), t.as_slice()))
        .collect();
    drop(_align_t);
    let matrix = {
        let _t = crate::profile::scope("create_multialignment_matrix");
        msa::multialignment_matrix(&rows)
    };
    let (ref_aln, supporting) = matrix.split_first().expect("the consensus row");

    let pfm = msa::pfm(supporting);
    let windows = {
        let _t = crate::profile::scope("get_contexts");
        contexts::contexts(ref_aln, k_size / 2)
    };

    // `max(3, weight * T)` --- the reference computes the same number twice,
    // as variant_threshold and context_threshold.
    let threshold = (weight as f64 * v_depth_ratio_threshold).max(3.0);
    let alternatives = {
        let _t = crate::profile::scope("get_alternative_ref_contexts");
        let fcm = contexts::frequency_context_matrix(&matrix, &windows, threshold);
        contexts::alternative_ref_contexts(&fcm, ref_aln, &windows, threshold)
    };

    let _loop_t = crate::profile::scope("correction loop over PFM");
    let mut corrected: Vec<Vec<u8>> = vec![Vec::with_capacity(ref_aln.len()); supporting.len()];

    for (i, &ref_nucl) in ref_aln.iter().enumerate() {
        let (start, stop) = windows[i];

        if !alternatives[i].is_empty() {
            let consensus_context = degap(&ref_aln[start..stop]);
            for (r, row) in supporting.iter().enumerate() {
                let read_nucl = row[i];
                if read_nucl == ref_nucl {
                    corrected[r].push(ref_nucl);
                    continue;
                }

                let read_context = &row[start..stop];
                let read_degapped = degap(read_context);
                // Seeded with the distance to the consensus, so an alternative
                // has to beat it strictly.
                let mut min_ed = editdist::bounded(&consensus_context, &read_degapped, 1000);
                let mut correct_to = ref_nucl;
                for alt in &alternatives[i] {
                    if read_context == alt.context {
                        correct_to = alt.variant;
                        break;
                    }
                    let ed = editdist::bounded(&degap(&alt.context), &read_degapped, 1000);
                    if ed < min_ed {
                        min_ed = ed;
                        correct_to = alt.variant;
                    }
                }
                corrected[r].push(correct_to);
            }
        } else {
            // A consensus gap with no alternatives emits no column at all.
            if ref_nucl == b'-' {
                continue;
            }
            for (r, row) in supporting.iter().enumerate() {
                let read_nucl = row[i];
                let keep = read_nucl != b'-'
                    && read_nucl != ref_nucl
                    && threshold <= f64::from(pfm[i].get(read_nucl))
                    && ref_aln[start..i] == row[start..i]
                    && ref_aln[i + 1..stop] == row[i + 1..stop];
                corrected[r].push(if keep { read_nucl } else { ref_nucl });
            }
        }
    }

    let mut others: IndexMap<usize, Vec<Correction>> = IndexMap::new();
    let mut curr_read = Vec::new();
    for (r, &(q_id, pos1, pos2)) in instance.iter().enumerate() {
        let seq = trim(&degap(&corrected[r]), k_size);
        // A read appearing twice would overwrite here, as it does in the
        // reference --- `find_most_supported_span` keys by r_id, so it cannot.
        if q_id == curr_read_id {
            curr_read = seq.clone();
        }
        others.entry(q_id).or_default().push(Correction {
            start: pos1 + k_size,
            stop: pos2,
            weight,
            seq,
        });
    }

    Some(Corrections { curr_read, others })
}

fn degap(v: &[u8]) -> Vec<u8> {
    v.iter().copied().filter(|&c| c != b'-').collect()
}

/// `s[k_size:-k_size]`, which Python leaves empty rather than panicking when
/// the string is shorter than the two trims.
fn trim(v: &[u8], k_size: usize) -> Vec<u8> {
    if v.len() <= 2 * k_size {
        return Vec::new();
    }
    v[k_size..v.len() - k_size].to_vec()
}

/// Align one supporting segment back to the consensus.
///
/// **Default: affine, via SIMD.** [`crate::simd`] explains why affine rather than
/// the reference's unit-cost edit distance — in short, edit distance was a
/// Python-era speed compromise, affine models indels better, and it is 2.7x
/// faster here besides. This stage was **33.7%** of the port's runtime on ten
/// real clusters, the single largest cost after the guard was fixed.
///
/// This is the port's deepest deliberate divergence from the reference: the CIGAR
/// chosen here builds the MSA, which builds the PFM, which decides every
/// corrected base. About 11% of alignments take a different path than edlib's
/// unit-cost one. `ISONCORRECT_EXACT_ALIGN=1` restores the exact edlib-compatible
/// path and `bench/equivalence.sh` sets it, so the byte-identity gate is intact.
fn segment_ops(seg: &[u8], consensus: &[u8], out: &mut Vec<align::CigarOp>) {
    if !exact_align()
        && seg.len() >= simd::MIN_LEN
        && consensus.len() >= simd::MIN_LEN
        && simd::global_ops(seg, consensus, Scoring::GUARD, out)
    {
        return;
    }
    let cigar = align::global(seg, consensus).cigar;
    *out = align::parse_cigar(&cigar).expect("the DP's own CIGAR parses");
}

/// Whether `ISONCORRECT_EXACT_ALIGN` asks for the edlib-compatible DP.
fn exact_align() -> bool {
    use std::sync::atomic::{AtomicU8, Ordering};
    static MODE: AtomicU8 = AtomicU8::new(u8::MAX);
    let cached = MODE.load(Ordering::Relaxed);
    if cached != u8::MAX {
        return cached == 1;
    }
    let on = std::env::var("ISONCORRECT_EXACT_ALIGN")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false);
    MODE.store(u8::from(on), Ordering::Relaxed);
    on
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Three identical reads: the consensus is the read, and correcting changes
    /// nothing but the two `k_size` trims.
    #[test]
    fn identical_reads_correct_to_themselves() {
        let read = b"ACGTTGCAACGTTGCAACGTACGT".to_vec();
        let seqs = vec![read.clone(), read.clone(), read.clone()];
        let instance = vec![(1, 0, 16), (2, 0, 16), (3, 0, 16)];
        let c = best_corrections(&instance, &seqs, 8, 0.1, 200).unwrap();
        assert_eq!(c.curr_read, &read[8..16]);
        assert_eq!(c.others.len(), 3);
    }

    #[test]
    fn a_lone_substitution_is_corrected_away() {
        let good = b"ACGTTGCAACGTTGCAACGTACGT".to_vec();
        let mut bad = good.clone();
        bad[12] = b'A'; // inside the span that survives the trim
        let seqs = vec![bad.clone(), good.clone(), good.clone(), good.clone()];
        let instance = vec![(1, 0, 16), (2, 0, 16), (3, 0, 16), (4, 0, 16)];
        let c = best_corrections(&instance, &seqs, 8, 0.1, 200).unwrap();
        assert_eq!(c.curr_read, &good[8..16], "the odd base out survived");
    }

    #[test]
    fn every_supporting_read_gets_its_own_region_back() {
        let read = b"ACGTTGCAACGTTGCAACGTACGT".to_vec();
        let seqs = vec![read.clone(), read.clone(), read.clone()];
        let instance = vec![(2, 0, 16), (1, 0, 16), (3, 2, 14)];
        let c = best_corrections(&instance, &seqs, 8, 0.1, 200).unwrap();
        // Keyed by r_id in partition order, and the current read is first.
        assert_eq!(c.others.keys().copied().collect::<Vec<_>>(), vec![2, 1, 3]);
        // start is pos1 + k_size, stop is pos2 untouched.
        assert_eq!((c.others[&3][0].start, c.others[&3][0].stop), (10, 14));
        assert!(c.others.values().all(|v| v[0].weight == 3));
    }

    #[test]
    fn an_empty_instance_has_nothing_to_correct() {
        assert!(best_corrections(&[], &[], 8, 0.1, 200).is_none());
    }

    #[test]
    fn a_span_shorter_than_the_two_trims_comes_back_empty() {
        assert_eq!(trim(b"ACGT", 8), b"");
        assert_eq!(trim(b"ACGTACGTACGTACGT", 8), b"");
        assert_eq!(trim(b"ACGTACGTXACGTACGT", 8), b"X");
    }
}

/// Differential oracle against `get_best_corrections` itself.
///
/// This is the end-to-end check the PFM never had: it runs the consensus, the
/// alignments, the matrix, the PFM, the contexts and the correction loop, and
/// compares the corrected strings the reference produced.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// CORRECTION_CASES=/tmp/d/correction_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml corrections::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    /// ```text
    /// C <case> <k_size> <T> <max_seqs_to_spoa> <curr_read_corr>
    /// S <case> <q_id> <pos1> <pos2> <segment>
    /// O <case> <q_id> <start> <stop> <weight> <corrected>
    /// ```
    ///
    /// Segments are recorded rather than whole reads, so a case is replayable
    /// without the fastq: the segment *is* `reads[q_id][1][pos1 : pos2 +
    /// k_size]`, and nothing downstream reads anything else.
    #[test]
    fn matches_the_reference_corrections() {
        let Ok(path) = std::env::var("CORRECTION_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CORRECTION_CASES={path} unreadable: {e}"));

        struct Case {
            k_size: usize,
            t: f64,
            max_seqs_to_spoa: usize,
            curr: String,
            instance: Vec<(usize, usize, usize)>,
            segments: Vec<Vec<u8>>,
            others: Vec<String>,
        }

        let mut cases: Vec<Case> = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "C" => cases.push(Case {
                    k_size: f[2].parse().unwrap(),
                    t: f[3].parse().unwrap(),
                    max_seqs_to_spoa: f[4].parse().unwrap(),
                    curr: f.get(5).unwrap_or(&"").to_string(),
                    instance: Vec::new(),
                    segments: Vec::new(),
                    others: Vec::new(),
                }),
                "S" => {
                    let c = cases.last_mut().unwrap();
                    c.instance.push((
                        f[2].parse().unwrap(),
                        f[3].parse().unwrap(),
                        f[4].parse().unwrap(),
                    ));
                    c.segments.push(f[5].as_bytes().to_vec());
                }
                "O" => cases.last_mut().unwrap().others.push(line[2..].to_string()),
                other => panic!("unknown record {other}"),
            }
        }
        assert!(!cases.is_empty(), "no cases in {path}");

        let (mut bad_curr, mut bad_other, mut shown) = (0usize, 0usize, 0usize);
        let mut mismatches: Vec<String> = Vec::new();
        for (i, case) in cases.iter().enumerate() {
            // Rebuild a read table the instance can index into: each segment is
            // placed at its own pos1, which is all the code ever reads.
            let mut seqs: Vec<Vec<u8>> =
                vec![Vec::new(); case.instance.iter().map(|t| t.0).max().unwrap()];
            for (&(q_id, pos1, _), seg) in case.instance.iter().zip(&case.segments) {
                let mut read = vec![b'A'; pos1];
                read.extend_from_slice(seg);
                seqs[q_id - 1] = read;
            }

            let got = best_corrections(
                &case.instance,
                &seqs,
                case.k_size,
                case.t,
                case.max_seqs_to_spoa,
            )
            .expect("a non-empty instance");

            if got.curr_read != case.curr.as_bytes() {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "case {i} current read:\n  python: {}\n  rust  : {}",
                        case.curr,
                        String::from_utf8_lossy(&got.curr_read)
                    );
                }
                mismatches.push(format!(
                    "C\t{i}\t{}",
                    String::from_utf8_lossy(&got.curr_read)
                ));
                bad_curr += 1;
            }

            let mine: Vec<String> = got
                .others
                .iter()
                .flat_map(|(q_id, regions)| {
                    regions.iter().map(move |r| {
                        format!(
                            "{q_id}\t{}\t{}\t{}\t{}",
                            r.start,
                            r.stop,
                            r.weight,
                            String::from_utf8_lossy(&r.seq)
                        )
                    })
                })
                .collect();
            let want: Vec<String> = case
                .others
                .iter()
                .map(|l| l.split_once('\t').map_or("", |x| x.1).to_string())
                .collect();
            if mine != want {
                if shown < 3 {
                    shown += 1;
                    eprintln!("case {i} regions:\n  python: {want:?}\n  rust  : {mine:?}");
                }
                for m in &mine {
                    mismatches.push(format!("O\t{i}\t{m}"));
                }
                bad_other += 1;
            }
        }

        // Mismatches are written out so they can be checked against the
        // reference under several PYTHONHASHSEEDs --- see the module docs on
        // set-iteration order, and `bench/check_seed_sensitivity.py`.
        if let Ok(out) = std::env::var("CORRECTION_MISMATCHES") {
            std::fs::write(&out, mismatches.join("\n") + "\n").expect("writable mismatch path");
            eprintln!("wrote {} mismatch records to {out}", mismatches.len());
        }

        assert_eq!(
            bad_curr,
            0,
            "{bad_curr} of {} corrected reads differed",
            cases.len()
        );
        assert_eq!(bad_other, 0, "{bad_other} region sets differed");
        eprintln!("matched the reference on all {} corrections", cases.len());
    }
}
