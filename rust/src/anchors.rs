//! The minimizer-pair anchor index, matching
//! `get_minimizer_combinations_database` and `minimizers_comb_iterator`.
//!
//! For every read, every ordered pair of minimizers `(m1, m2)` whose span falls
//! in `(xmin, xmax]` is recorded, so that other reads sharing the same pair can
//! be found later. Details that change results:
//!
//! * the span test is **exclusive at the low end and inclusive at the high
//!   end**: `x_low < p2 - p1 <= x_high`;
//! * the inner scan **breaks** as soon as the span exceeds `x_high`, relying on
//!   minimizer positions being increasing;
//! * `minimizers_comb_iterator` yields each partner list **reversed**
//!   (`m1_curr_spans[::-1]`), and that order reaches the payload, which decides
//!   the order supporting reads are collected downstream. It is not free to
//!   change;
//! * a pair where both k-mers are poly-A (`'A' * k`) is skipped entirely;
//! * after building, entries with a single occurrence are dropped, as are
//!   entries more abundant than the read count when the partner is poly-A or
//!   the count exceeds ten times the read count.
//!
//! Unlike the reference this stores only surviving entries. The reference's
//! `del` is immediately undone by a `defaultdict` read, leaving ~94% of keys
//! present but empty; a lookup on an empty entry and a lookup on an absent one
//! both mean "no support", so the data is the same. See PORTING.md.

use indexmap::IndexMap;
use rustc_hash::FxBuildHasher;

/// `(r_id, p1, p2)` --- a read supporting one anchor pair, and where.
pub type Occurrence = (usize, usize, usize);

/// A minimizer and its position in a read.
pub type Minimizer = (Vec<u8>, usize);

/// One read's minimizers, keyed by `r_id`, in `r_id` order.
pub type MinimizersByRead = [(usize, Vec<Minimizer>)];

/// A borrowed minimizer and its position.
pub type MinimizerRef<'a> = (&'a [u8], usize);

/// One `m1` and the partners it spans to, as `comb_iterator` yields them.
pub type Combination<'a> = (MinimizerRef<'a>, Vec<MinimizerRef<'a>>);

/// Anchor pair index: `m1 -> m2 -> occurrences`.
///
/// Insertion-ordered so dumps stay diffable against the reference. Only the
/// payload order actually feeds results.
#[derive(Debug, Default)]
pub struct AnchorDb {
    map: IndexMap<Vec<u8>, PartnerMap, FxBuildHasher>,
    /// Total pairs generated before filtering; the reference prints this.
    pub generated: usize,
    /// Entries dropped for having a single occurrence.
    pub singletons: usize,
    /// Entries dropped as poly-A or excessively abundant.
    pub high_abundance: usize,
}

/// The partners of one `m1`.
///
/// `FxBuildHasher` rather than SipHash on both levels: these maps are read by
/// lookup on the hot path, and `IndexMap` preserves *insertion* order regardless
/// of hasher, so [`AnchorDb::iter`] — the only place order is observed, when
/// dumping — is unaffected. Same for `retain`, which keeps relative order.
pub type PartnerMap = IndexMap<Vec<u8>, Vec<Occurrence>, FxBuildHasher>;

impl AnchorDb {
    /// Occurrences for a pair; empty when unknown.
    ///
    /// The reference indexes a `defaultdict`, so a miss yields an empty array
    /// rather than an error. Same contract here, without the insertion.
    pub fn get(&self, m1: &[u8], m2: &[u8]) -> &[Occurrence] {
        self.partners(m1)
            .and_then(|inner| inner.get(m2))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// All partners of `m1`, so a caller sweeping many `m2` for one `m1` pays
    /// the outer lookup once. `find_most_supported_span` does exactly that.
    pub fn partners(&self, m1: &[u8]) -> Option<&PartnerMap> {
        self.map.get(m1)
    }

    /// Number of surviving `(m1, m2)` entries.
    pub fn len(&self) -> usize {
        self.map.values().map(PartnerMap::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate surviving entries in insertion order, for dumping.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &[u8], &[Occurrence])> {
        self.map.iter().flat_map(|(m1, inner)| {
            inner
                .iter()
                .map(move |(m2, occ)| (m1.as_slice(), m2.as_slice(), occ.as_slice()))
        })
    }
}

/// Partner minimizers for each `m1`, in the reference's order.
///
/// Mirrors `minimizers_comb_iterator`, including the trailing reversal and the
/// fact that the last minimizer is never used as an `m1` (`minimizers[:-1]`).
pub fn comb_iterator(
    minimizers: &[Minimizer],
    x_low: usize,
    x_high: usize,
) -> Vec<Combination<'_>> {
    let mut out = Vec::new();
    if minimizers.is_empty() {
        return out;
    }
    for i in 0..minimizers.len() - 1 {
        let (ref m1, p1) = minimizers[i];
        let mut spans = Vec::new();
        for (m2, p2) in &minimizers[i + 1..] {
            let span = p2.saturating_sub(p1);
            if span > x_high {
                break;
            }
            if span > x_low {
                spans.push((m2.as_slice(), *p2));
            }
        }
        // The reference yields m1_curr_spans[::-1]; this order reaches the
        // payload and therefore the results.
        spans.reverse();
        out.push(((m1.as_slice(), p1), spans));
    }
    out
}

/// Build the anchor index over every read's minimizers.
///
/// `minimizers_by_read` must be in `r_id` order, matching the reference's
/// iteration over `M`.
pub fn build(
    minimizers_by_read: &MinimizersByRead,
    k: usize,
    x_low: usize,
    x_high: usize,
    n_reads: usize,
) -> AnchorDb {
    let forbidden = vec![b'A'; k];
    let mut db = AnchorDb::default();

    for (r_id, minimizers) in minimizers_by_read {
        for ((m1, p1), spans) in comb_iterator(minimizers, x_low, x_high) {
            for (m2, p2) in spans {
                if m2 == forbidden.as_slice() && m1 == forbidden.as_slice() {
                    continue;
                }
                db.generated += 1;
                // `entry(m1.to_vec())` would allocate both keys on every pair,
                // including the overwhelming majority that already exist. Look up
                // first and only allocate on a genuine miss. A new key still lands
                // at the end, so insertion order is unchanged.
                if !db.map.contains_key(m1) {
                    db.map.insert(m1.to_vec(), PartnerMap::default());
                }
                let inner = db.map.get_mut(m1).expect("just inserted");
                if !inner.contains_key(m2) {
                    inner.insert(m2.to_vec(), Vec::new());
                }
                inner
                    .get_mut(m2)
                    .expect("just inserted")
                    .push((*r_id, p1, p2));
            }
        }
    }

    // Filtering, matching the reference's two rules. Counters are local because
    // `retain` holds a mutable borrow of the map for its duration.
    let mut singletons = 0usize;
    let mut high_abundance = 0usize;
    for inner in db.map.values_mut() {
        inner.retain(|m2, occ| {
            // The reference tests `len(...) > 3` on a flat [r_id, p1, p2, ...]
            // array, i.e. "more than one occurrence".
            if occ.len() <= 1 {
                singletons += 1;
                return false;
            }
            if occ.len() > n_reads
                && (m2.as_slice() == forbidden.as_slice() || occ.len() > 10 * n_reads)
            {
                high_abundance += 1;
                return false;
            }
            true
        });
    }
    db.map.retain(|_, inner| !inner.is_empty());
    db.singletons = singletons;
    db.high_abundance = high_abundance;

    db
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mins(pairs: &[(&str, usize)]) -> Vec<Minimizer> {
        pairs
            .iter()
            .map(|(m, p)| (m.as_bytes().to_vec(), *p))
            .collect()
    }

    #[test]
    fn span_is_exclusive_low_inclusive_high() {
        // x_low < p2 - p1 <= x_high
        let m = mins(&[("AAA", 0), ("CCC", 10), ("GGG", 20), ("TTT", 30)]);
        let got = comb_iterator(&m, 10, 20);
        let (_, spans) = &got[0]; // m1 at p1 = 0
                                  // span 10 excluded (not > 10), span 20 included, span 30 excluded
        let positions: Vec<usize> = spans.iter().map(|(_, p)| *p).collect();
        assert_eq!(positions, vec![20]);
    }

    #[test]
    fn partner_list_is_reversed() {
        // The reference yields m1_curr_spans[::-1]; partners must come back in
        // descending position order.
        let m = mins(&[("AAA", 0), ("CCC", 15), ("GGG", 25), ("TTT", 35)]);
        let got = comb_iterator(&m, 10, 40);
        let positions: Vec<usize> = got[0].1.iter().map(|(_, p)| *p).collect();
        assert_eq!(positions, vec![35, 25, 15], "partners must be reversed");
    }

    #[test]
    fn last_minimizer_is_never_an_m1() {
        let m = mins(&[("AAA", 0), ("CCC", 30)]);
        let got = comb_iterator(&m, 10, 40);
        assert_eq!(got.len(), 1, "minimizers[:-1] excludes the final entry");
    }

    #[test]
    fn scan_breaks_past_the_upper_bound() {
        // Everything beyond x_high is unreachable, so a far-away partner that
        // would satisfy no test must not appear.
        let m = mins(&[("AAA", 0), ("CCC", 20), ("GGG", 1000), ("TTT", 1001)]);
        let got = comb_iterator(&m, 10, 40);
        assert_eq!(got[0].1.len(), 1);
    }

    #[test]
    fn singletons_are_dropped() {
        // One read, so every pair occurs exactly once and nothing survives.
        let m = vec![(1usize, mins(&[("AAA", 0), ("CCC", 30), ("GGG", 60)]))];
        let db = build(&m, 3, 10, 100, 1);
        assert!(db.is_empty());
        assert!(db.singletons > 0);
    }

    #[test]
    fn pairs_shared_by_two_reads_survive() {
        let one = mins(&[("ACG", 0), ("TTG", 30)]);
        let two = mins(&[("ACG", 5), ("TTG", 40)]);
        let db = build(&[(1, one), (2, two)], 3, 10, 100, 2);
        let occ = db.get(b"ACG", b"TTG");
        assert_eq!(occ, &[(1, 0, 30), (2, 5, 40)]);
    }

    #[test]
    fn poly_a_pairs_are_skipped() {
        // m1 == m2 == 'A'*k is skipped at generation time.
        let poly = mins(&[("AAA", 0), ("AAA", 30)]);
        let db = build(&[(1, poly.clone()), (2, poly)], 3, 10, 100, 2);
        assert!(db.get(b"AAA", b"AAA").is_empty());
    }

    #[test]
    fn missing_pairs_return_empty_not_an_error() {
        // The reference indexes a defaultdict, so a miss is "no support".
        let db = AnchorDb::default();
        assert!(db.get(b"NNN", b"NNN").is_empty());
    }

    #[test]
    fn occurrence_order_follows_read_order() {
        // Payload order decides the order supporting reads are collected, so it
        // must follow r_id order.
        let a = mins(&[("ACG", 0), ("TTG", 30)]);
        let db = build(&[(1, a.clone()), (2, a.clone()), (3, a)], 3, 10, 100, 3);
        let ids: Vec<usize> = db.get(b"ACG", b"TTG").iter().map(|(r, _, _)| *r).collect();
        assert_eq!(ids, vec![1, 2, 3]);
    }
}
