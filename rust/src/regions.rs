//! Skipping anchors whose span was already corrected, from `isoncorrect_main`.
//!
//! After a read is corrected, the regions it covered in *other* reads are
//! recorded in `previously_corrected_regions`, so later reads do not redo the
//! work. This is what makes correction order-dependent.
//!
//! The filtering has two parts, and the second is easy to miss:
//!
//! 1. keep an anchor span unless **both** ends fall inside already-considered
//!    positions;
//! 2. *re-admit* a span that rule 1 rejected when its two ends belong to
//!    **different contiguous groups** of considered positions — i.e. the span
//!    bridges a gap, so the middle has not actually been corrected.
//!
//! Under `--exact` the reference clears `previously_corrected_regions` at the
//! start of every read, so both rules become inert and every anchor is kept.

use std::collections::{BTreeSet, HashMap};

/// Positions covered by already-corrected regions: `range(p1, p2)` per region.
pub fn considered_positions(regions: &[(usize, usize)]) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    for &(p1, p2) in regions {
        // Python's range(p1, p2) is empty when p2 <= p1.
        for pos in p1..p2 {
            out.insert(pos);
        }
    }
    out
}

/// Group id per position, splitting at gaps.
///
/// Mirrors the reference loop over consecutive sorted positions: a jump of more
/// than one starts a new group. Positions not present get no entry.
pub fn pos_groups(positions: &BTreeSet<usize>) -> HashMap<usize, usize> {
    let mut groups = HashMap::new();
    if positions.len() <= 1 {
        // The reference guards with `if len(...) > 1`, leaving pos_group empty.
        return groups;
    }
    let sorted: Vec<usize> = positions.iter().copied().collect();
    let mut group_id = 0usize;
    let (mut last_p1, mut last_p2) = (sorted[0], sorted[0]);
    for w in sorted.windows(2) {
        let (p1, p2) = (w[0], w[1]);
        if p2 > p1 + 1 {
            groups.insert(p1, group_id);
            group_id += 1;
            groups.insert(p2, group_id);
        } else {
            groups.insert(p1, group_id);
        }
        last_p1 = p1;
        last_p2 = p2;
    }
    // Trailing `if p2 == p1 + 1: pos_group[p2] = group_id`, using the loop's
    // final pair.
    if last_p2 == last_p1 + 1 {
        groups.insert(last_p2, group_id);
    }
    groups
}

/// Anchor partners still worth examining for `p1`.
///
/// Returns the kept spans in the reference's order: rule-1 survivors first,
/// then the rule-2 re-admissions appended.
pub fn filter_spans<'a>(
    spans: &[(&'a [u8], usize)],
    p1: usize,
    k: usize,
    considered: &BTreeSet<usize>,
    groups: &HashMap<usize, usize>,
) -> Vec<(&'a [u8], usize)> {
    let both_ends_considered =
        |p2: usize| considered.contains(&(p1 + k)) && p2 > 0 && considered.contains(&(p2 - 1));

    let mut kept: Vec<(&[u8], usize)> = spans
        .iter()
        .copied()
        .filter(|&(_, p2)| !both_ends_considered(p2))
        .collect();

    // Rule 2: a rejected span whose ends sit in different groups spans a gap,
    // so it is put back.
    let readmitted: Vec<(&[u8], usize)> = spans
        .iter()
        .copied()
        .filter(|&(_, p2)| {
            if !both_ends_considered(p2) {
                return false;
            }
            match (
                groups.get(&(p1 + k)),
                p2.checked_sub(1).and_then(|x| groups.get(&x)),
            ) {
                (Some(a), Some(b)) => a != b,
                _ => false,
            }
        })
        .collect();

    kept.extend(readmitted);
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions_expand_to_half_open_ranges() {
        let p = considered_positions(&[(5, 8)]);
        assert_eq!(p.iter().copied().collect::<Vec<_>>(), vec![5, 6, 7]);
    }

    #[test]
    fn empty_or_inverted_regions_contribute_nothing() {
        assert!(considered_positions(&[]).is_empty());
        assert!(considered_positions(&[(9, 9)]).is_empty());
        assert!(considered_positions(&[(9, 4)]).is_empty());
    }

    #[test]
    fn a_single_position_yields_no_groups() {
        // The reference guards on `len(...) > 1`.
        let p = considered_positions(&[(5, 6)]);
        assert!(pos_groups(&p).is_empty());
    }

    #[test]
    fn contiguous_positions_form_one_group() {
        let p = considered_positions(&[(0, 5)]);
        let g = pos_groups(&p);
        let ids: BTreeSet<usize> = g.values().copied().collect();
        assert_eq!(ids.len(), 1, "contiguous run must be a single group");
    }

    #[test]
    fn a_gap_starts_a_new_group() {
        let p = considered_positions(&[(0, 3), (10, 13)]);
        let g = pos_groups(&p);
        assert_ne!(g[&0], g[&12], "positions across a gap differ");
        assert_eq!(g[&0], g[&2], "positions within a run agree");
    }

    #[test]
    fn nothing_is_filtered_when_no_region_was_corrected() {
        // The --exact case: every anchor survives.
        let spans: Vec<(&[u8], usize)> = vec![(b"AAA", 40), (b"CCC", 60)];
        let kept = filter_spans(&spans, 10, 9, &BTreeSet::new(), &HashMap::new());
        assert_eq!(kept.len(), 2);
    }

    #[test]
    fn spans_inside_one_corrected_run_are_dropped() {
        let considered = considered_positions(&[(0, 100)]);
        let groups = pos_groups(&considered);
        let spans: Vec<(&[u8], usize)> = vec![(b"AAA", 50)];
        // p1 + k = 19 and p2 - 1 = 49 are both inside the run, same group.
        let kept = filter_spans(&spans, 10, 9, &considered, &groups);
        assert!(kept.is_empty(), "fully covered span must be dropped");
    }

    #[test]
    fn spans_bridging_a_gap_are_readmitted() {
        // Two corrected runs with an uncorrected gap between them; a span whose
        // ends land in different runs still has work to do in the middle.
        let considered = considered_positions(&[(0, 30), (60, 100)]);
        let groups = pos_groups(&considered);
        let spans: Vec<(&[u8], usize)> = vec![(b"AAA", 80)];
        // p1 + k = 19 (first run), p2 - 1 = 79 (second run).
        let kept = filter_spans(&spans, 10, 9, &considered, &groups);
        assert_eq!(kept.len(), 1, "gap-bridging span must be kept");
    }

    #[test]
    fn readmitted_spans_come_after_the_survivors() {
        // Order matters: the reference appends the rule-2 list to the rule-1
        // list, and this order reaches find_most_supported_span.
        let considered = considered_positions(&[(0, 30), (60, 100)]);
        let groups = pos_groups(&considered);
        let spans: Vec<(&[u8], usize)> = vec![(b"BRIDGE", 80), (b"FREE", 200)];
        let kept = filter_spans(&spans, 10, 9, &considered, &groups);
        assert_eq!(kept[0].0, b"FREE", "rule-1 survivor first");
        assert_eq!(kept[1].0, b"BRIDGE", "rule-2 re-admission after");
    }
}
