//! Collecting supporting segments per anchor pair, matching
//! `find_most_supported_span`.
//!
//! For one `m1` and each of its partners `m2`, gather the reads that also span
//! that pair and decide, for each, whether its segment is close enough to the
//! current read's to count as support. The result is a candidate interval that
//! WIS later chooses among.
//!
//! Details that decide the answer:
//!
//! * `to_add` is a **dict, iterated in insertion order** to build the payload.
//!   That order becomes the order supporting segments are written for the
//!   consensus, and POA is order-sensitive — so it reaches the output and is not
//!   free to change. Hence `IndexMap`;
//! * an anchor pair is only considered when at least **three** reads carry it
//!   (`len(relevant_reads)//3 >= 3`);
//! * `elif relevant_read_id in reads_visited: pass` does **nothing** and falls
//!   through to the alignment. A read already seen for this `m2` is therefore
//!   re-aligned rather than short-circuited — easy to misread as a `continue`;
//! * `already_computed` persists **across calls** for the read being corrected,
//!   so results depend on the order anchors are visited;
//! * the estimate `ed_est` is float arithmetic via `math.fabs`, compared against
//!   a float threshold, so it must be computed in the same order.
//!
//! # Why this is the port's hot spot, and what was done about it
//!
//! 28.4% of runtime after the MSA and `align` work, over 342 940 calls — so
//! per-call cost rather than one hot inner loop. Measured on three real clusters:
//! the bounded edit distance is only **27%** of the stage. The other 73% is the
//! scaffolding around **15.4 million** iterations of the inner loop (1 715 per
//! read), which was doing, per iteration, a SipHash of an ~80-byte key.
//!
//! Three changes, none of which touch what is computed:
//!
//! * **`FxHashMap` instead of SipHash** on `added_strings`, `already_computed`
//!   and `reads_visited`. All three are lookup-only — never iterated — so the
//!   hasher cannot reach output. `to_add` stays an `IndexMap` precisely because
//!   *its* order does;
//! * **the scratch containers are reused across partners and calls** rather than
//!   allocated fresh each time. `clear()` keeps the capacity;
//! * **`ref_seq` is a slice**, not a `to_vec()` per partner.
//!
//! The `Vec<u8>` keys of `added_strings` are still allocated on insert, but
//! inserts are bounded by the reads actually added, while the lookups were the
//! 15-million-iteration cost.

use indexmap::IndexMap;
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};

use crate::anchors::AnchorDb;
use crate::editdist;
use crate::quality;

/// A candidate interval, matching the tuple appended to `all_intervals`:
/// `(p1 + k_size, p2, len(seqs)//3, seqs)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Candidate {
    pub start: usize,
    pub stop: usize,
    /// Supporting segment count, including the read being corrected.
    pub support: usize,
    /// `(r_id, pos1, pos2)` per supporting segment, in `to_add` insertion order.
    pub instance: Vec<(usize, usize, usize)>,
}

/// Cached alignment for one read: `(ref_start, ref_end, read_start, read_end, ed)`.
type Computed = (usize, usize, usize, usize, f64);

/// Per-read state for span finding.
///
/// `already_computed` is deliberately owned here rather than rebuilt per call:
/// the reference threads one dict through every anchor of a read, and clearing
/// it would change results.
pub struct SpanFinder<'a> {
    /// Sequences by `r_id` (1-based; index with `r_id - 1`).
    pub seqs: &'a [Vec<u8>],
    /// Quality prefix sums by `r_id`, from [`quality::prefix_sums`].
    pub qvs: &'a [Vec<f64>],
    pub k: usize,
    already_computed: FxHashMap<usize, Computed>,
    /// Alignments performed; the reference reports this as `tmp_cnt`.
    pub edlib_calls: usize,
    /// Per-partner scratch, reused across partners and calls. Cleared at the top
    /// of each partner, so carrying it here is purely an allocation saving.
    scratch: Scratch,
}

/// The per-partner containers. `to_add` is insertion-ordered because that order
/// becomes the payload order — `IndexMap` keeps that regardless of hasher, so
/// only its *lookups* get the faster one. The rest are lookup-only.
#[derive(Default)]
struct Scratch {
    to_add: IndexMap<usize, (usize, usize, usize, f64), FxBuildHasher>,
    added_strings: FxHashMap<Vec<u8>, f64>,
    reads_visited: FxHashSet<usize>,
    /// Recycled `added_strings` keys. Each landed edit distance inserts the
    /// ~80-byte segment as a key, 15 million times per three clusters; the map is
    /// cleared per partner, so the buffers can simply be handed back.
    key_pool: Vec<Vec<u8>>,
}

/// An `added_strings` key, from the pool when one is free.
fn key_for(pool: &mut Vec<Vec<u8>>, bytes: &[u8]) -> Vec<u8> {
    let mut v = pool.pop().unwrap_or_default();
    v.clear();
    v.extend_from_slice(bytes);
    v
}

impl<'a> SpanFinder<'a> {
    pub fn new(seqs: &'a [Vec<u8>], qvs: &'a [Vec<f64>], k: usize) -> Self {
        Self {
            seqs,
            qvs,
            k,
            already_computed: FxHashMap::default(),
            edlib_calls: 0,
            scratch: Scratch::default(),
        }
    }

    /// Clear the cross-anchor cache, as `--exact` does per read.
    pub fn reset_cache(&mut self) {
        self.already_computed.clear();
    }

    /// Append candidate intervals for `m1` and each of its partners.
    pub fn find(
        &mut self,
        r_id: usize,
        m1: &[u8],
        p1: usize,
        partners: &[(&[u8], usize)],
        db: &AnchorDb,
        out: &mut Vec<Candidate>,
    ) {
        // Lift these out of `self` so the borrow checker does not tie the
        // sequence slices to the mutable borrows of the counters below. Both
        // are `&'a`, so copying the reference is free.
        let k = self.k;
        let seqs = self.seqs;
        let qvs = self.qvs;
        // One outer lookup for `m1`, not one per partner.
        let Some(m1_partners) = db.partners(m1) else {
            return;
        };
        for &(m2, p2) in partners {
            let relevant = m1_partners.get(m2).map(Vec::as_slice).unwrap_or(&[]);
            // The reference tests len(array)//3 >= 3 on the flat payload.
            if relevant.len() < 3 {
                continue;
            }

            let seq = &seqs[r_id - 1];
            let ref_end = p2 + k;
            if ref_end > seq.len() || p1 > ref_end {
                continue;
            }
            let ref_seq: &[u8] = &seq[p1..ref_end];
            let p_error_ref = quality::mean_error(&qvs[r_id - 1], p1, ref_end);

            let Scratch {
                to_add,
                added_strings,
                reads_visited,
                key_pool,
            } = &mut self.scratch;
            to_add.clear();
            key_pool.extend(added_strings.drain().map(|(k, _)| k));
            reads_visited.clear();

            to_add.insert(r_id, (r_id, p1, p2, 0.0));
            added_strings.insert(key_for(key_pool, ref_seq), 0.0);

            for &(other, pos1, pos2) in relevant {
                if other == r_id {
                    continue;
                }
                let other_seq_full = &seqs[other - 1];
                let read_end = pos2 + k;
                if read_end > other_seq_full.len() || pos1 > read_end {
                    continue;
                }
                let read_seq = &other_seq_full[pos1..read_end];

                if read_seq == ref_seq {
                    to_add.insert(other, (other, pos1, pos2, 0.0));
                    reads_visited.insert(other);
                    self.already_computed
                        .insert(other, (p1, p2, pos1, pos2, 0.0));
                    continue;
                }

                // NOTE: the reference's `elif ... in reads_visited: pass` is a
                // no-op that falls through to the alignment below. Only skip the
                // shortcuts when the read has *not* been visited for this m2.
                if !reads_visited.contains(&other) {
                    if let Some(&prev_ed) = added_strings.get(read_seq) {
                        if let Some(existing) = to_add.get(&other) {
                            if prev_ed >= existing.3 {
                                continue;
                            }
                        }
                        to_add.insert(other, (other, pos1, pos2, prev_ed));
                        reads_visited.insert(other);
                        self.already_computed
                            .insert(other, (p1, p2, pos1, pos2, prev_ed));
                        continue;
                    }

                    if let Some(&(c_ref_s, c_ref_e, c_read_s, c_read_e, c_ed)) =
                        self.already_computed.get(&other)
                    {
                        if c_read_s <= pos1 && pos2 <= c_read_e && c_ref_s <= p1 && p2 <= c_ref_e {
                            let p_error_read = quality::mean_error(&qvs[other - 1], pos1, read_end);
                            let thresh = (p_error_ref + p_error_read) * ref_seq.len() as f64;
                            let read_beg_diff = pos1 as i64 - c_read_s as i64;
                            let read_end_diff = pos2 as i64 - c_read_e as i64;
                            let ref_beg_diff = p1 as i64 - c_ref_s as i64;
                            let ref_end_diff = p2 as i64 - c_ref_e as i64;
                            let ed_est = c_ed
                                + ((ref_end_diff - read_end_diff) as f64).abs()
                                + ((read_beg_diff - ref_beg_diff) as f64).abs();
                            if (0.0..=thresh).contains(&ed_est) {
                                if let Some(existing) = to_add.get(&other) {
                                    if ed_est >= existing.3 {
                                        continue;
                                    }
                                }
                                to_add.insert(other, (other, pos1, pos2, ed_est));
                                added_strings.insert(key_for(key_pool, read_seq), ed_est);
                                reads_visited.insert(other);
                                // NOTE: already_computed is deliberately NOT
                                // refreshed here; the reference does not.
                                continue;
                            }
                        }
                    }
                }

                // Fall-through: align. The reference bounds by len(ref_seq).
                let editdist = editdist::bounded(ref_seq, read_seq, ref_seq.len());
                self.edlib_calls += 1;
                if editdist >= 0 {
                    let ed = editdist as f64;
                    if let Some(existing) = to_add.get(&other) {
                        if ed >= existing.3 {
                            continue;
                        }
                    }
                    to_add.insert(other, (other, pos1, pos2, ed));
                    added_strings.insert(key_for(key_pool, read_seq), ed);
                    reads_visited.insert(other);
                    self.already_computed
                        .insert(other, (p1, p2, pos1, pos2, ed));
                }
            }

            let instance: Vec<(usize, usize, usize)> =
                to_add.values().map(|&(r, a, b, _)| (r, a, b)).collect();
            out.push(Candidate {
                start: p1 + k,
                stop: p2,
                support: instance.len(),
                instance,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::{self, Minimizer};

    fn setup() -> (Vec<Vec<u8>>, Vec<Vec<f64>>) {
        // Four reads sharing an anchor pair, one with a substitution.
        let seqs: Vec<Vec<u8>> = [
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGACGT",
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGACGT",
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGACGT",
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTGGCACGT",
        ]
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
        let qvs = seqs
            .iter()
            .map(|s| quality::prefix_sums(&vec![b'+'; s.len()]))
            .collect();
        (seqs, qvs)
    }

    #[test]
    fn anchor_needs_three_supporting_reads() {
        let (seqs, qvs) = setup();
        let db = AnchorDb::default(); // empty: no pair reaches three
        let mut f = SpanFinder::new(&seqs, &qvs, 4);
        let mut out = Vec::new();
        f.find(1, b"ACGT", 0, &[(b"ACGT", 30)], &db, &mut out);
        assert!(out.is_empty(), "fewer than three reads must yield nothing");
    }

    fn build_db(seqs: &[Vec<u8>], k: usize) -> AnchorDb {
        // Every read carries ACGT at 0 and ACGT at 33.
        let by_read: Vec<(usize, Vec<Minimizer>)> = (1..=seqs.len())
            .map(|r| {
                (
                    r,
                    vec![
                        (b"ACGT".to_vec(), 0usize),
                        (b"GGG".to_vec(), 30usize),
                        (b"ACGT".to_vec(), 33usize),
                    ],
                )
            })
            .collect();
        anchors::build(&by_read, k, 10, 100, seqs.len())
    }

    #[test]
    fn identical_segments_are_supported_without_alignment() {
        let (seqs, qvs) = setup();
        let db = build_db(&seqs, 4);
        let mut f = SpanFinder::new(&seqs, &qvs, 4);
        let mut out = Vec::new();
        f.find(1, b"ACGT", 0, &[(b"ACGT".as_slice(), 33)], &db, &mut out);

        assert_eq!(out.len(), 1);
        let c = &out[0];
        // Reads 1..3 are identical; read 4 differs by one base and still aligns
        // within the bound, so all four support the interval.
        assert_eq!(c.support, 4);
        // The current read is always first: to_add[r_id] is inserted first.
        assert_eq!(c.instance[0].0, 1);
        // Exactly one alignment: the three identical reads short-circuit.
        assert_eq!(f.edlib_calls, 1, "identical segments must not be aligned");
    }

    #[test]
    fn interval_bounds_match_the_reference_tuple() {
        let (seqs, qvs) = setup();
        let db = build_db(&seqs, 4);
        let mut f = SpanFinder::new(&seqs, &qvs, 4);
        let mut out = Vec::new();
        f.find(1, b"ACGT", 0, &[(b"ACGT".as_slice(), 33)], &db, &mut out);
        // (p1 + k_size, p2, ...)
        assert_eq!((out[0].start, out[0].stop), (4, 33)); // p1 + k_size, p2
    }

    #[test]
    fn the_corrected_read_is_first_in_the_payload() {
        // POA is order-sensitive, so this ordering reaches the output.
        let (seqs, qvs) = setup();
        let db = build_db(&seqs, 4);
        let mut f = SpanFinder::new(&seqs, &qvs, 4);
        let mut out = Vec::new();
        f.find(3, b"ACGT", 0, &[(b"ACGT".as_slice(), 33)], &db, &mut out);
        assert_eq!(out[0].instance[0].0, 3);
    }

    #[test]
    fn cache_survives_across_calls_until_reset() {
        let (seqs, qvs) = setup();
        let db = build_db(&seqs, 4);
        let mut f = SpanFinder::new(&seqs, &qvs, 4);
        let mut out = Vec::new();
        f.find(1, b"ACGT", 0, &[(b"ACGT".as_slice(), 33)], &db, &mut out);
        let after_first = f.edlib_calls;
        assert!(after_first > 0);
        f.reset_cache();
        f.find(1, b"ACGT", 0, &[(b"ACGT".as_slice(), 33)], &db, &mut out);
        assert!(
            f.edlib_calls > after_first,
            "resetting the cache must force realignment"
        );
    }
}
