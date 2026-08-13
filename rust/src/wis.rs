//! Weighted Interval Scheduling, matching `solve_WIS` and `fill_p2`.
//!
//! Given the candidate intervals found for a read, pick a maximum-weight set of
//! non-overlapping ones to correct over. Standard WIS dynamic programming, with
//! several reference-specific details that decide the answer:
//!
//! * the weight is `(w - 1) * (stop - start + epsilon)` with
//!   `epsilon = 0.0001` --- `w - 1` because the read's own interval is counted
//!   in the support, and the epsilon keeps zero-length intervals from being
//!   free. This is **floating point**, and the comparisons below are strict, so
//!   the arithmetic has to be performed in the same order;
//! * `fill_p2` builds `stop -> j` as a dict comprehension, so when several
//!   intervals share a `stop` the **last** one wins, and only intervals with
//!   `start < stop` are considered at all;
//! * the traceback takes an interval only on a strict `>`, so exact ties fall
//!   through to `j - 1`;
//! * the returned indices are in **decreasing** order; the caller reverses them
//!   (`opt_indicies[::-1]`).
//!
//! Input must already be sorted by `(stop, start)`, as `isoncorrect_main` does
//! with `all_intervals.sort(key=lambda x: (x[1], x[0]))`.

/// One candidate interval: `(start, stop, support)`.
///
/// `support` counts the supporting reads *including* the read being corrected,
/// which is why the weight subtracts one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub stop: usize,
    pub support: usize,
}

const EPSILON: f64 = 0.0001;

/// `p[j]`: for each interval, the index one past the last interval that ends at
/// or before its start. Mirrors `fill_p2`.
fn fill_p(intervals: &[Interval]) -> Vec<usize> {
    // p is 1-based with a leading placeholder, matching `p = [None]`.
    let mut p = vec![0usize; intervals.len() + 1];
    if intervals.is_empty() {
        return p;
    }

    // Python: {stop: j for j, ... if start < stop}. Later entries overwrite
    // earlier ones, so the *last* interval with a given stop wins.
    let max_stop = intervals[intervals.len() - 1].stop;
    let mut coord_to_max_j = vec![0usize; max_stop + 1];
    let mut stop_to_max_j: Vec<Option<usize>> = vec![None; max_stop + 1];
    for (j, iv) in intervals.iter().enumerate() {
        if iv.start < iv.stop && iv.stop <= max_stop {
            stop_to_max_j[iv.stop] = Some(j);
        }
    }

    let mut j_max = 0usize;
    for i in 0..=max_stop {
        if let Some(j) = stop_to_max_j[i] {
            j_max = j;
        }
        coord_to_max_j[i] = j_max;
    }

    for (j, iv) in intervals.iter().enumerate() {
        p[j + 1] = coord_to_max_j[iv.start.min(max_stop)];
    }
    p
}

/// Weight of interval `j`, exactly as the reference computes it.
fn weight(iv: &Interval) -> f64 {
    // (w - 1) * (stop - start + epsilon). `w - 1` can be 0 but never negative
    // in practice: the read itself always contributes one unit of support.
    let w = iv.support as f64 - 1.0;
    let span = (iv.stop as f64 - iv.start as f64) + EPSILON;
    w * span
}

/// Chosen interval indices, in the reference's order (decreasing).
///
/// Returns an empty vector when nothing is worth correcting --- which the
/// caller treats as "leave the read alone".
pub fn solve(intervals: &[Interval]) -> Vec<usize> {
    let mut opt_indices = Vec::new();
    if intervals.is_empty() {
        return opt_indices;
    }

    let p = fill_p(intervals);

    // v and OPT are 1-based with a leading placeholder.
    let mut v = vec![0.0f64; intervals.len() + 1];
    for (j, iv) in intervals.iter().enumerate() {
        v[j + 1] = weight(iv);
    }

    let mut opt = vec![0.0f64; intervals.len() + 1];
    for j in 1..=intervals.len() {
        opt[j] = (v[j] + opt[p[j]]).max(opt[j - 1]);
    }

    let mut j = intervals.len();
    while j > 0 {
        // Strict `>`: an exact tie falls through to j - 1.
        if v[j] + opt[p[j]] > opt[j - 1] {
            opt_indices.push(j - 1);
            j = p[j];
        } else {
            j -= 1;
        }
    }
    opt_indices
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: usize, stop: usize, support: usize) -> Interval {
        Interval {
            start,
            stop,
            support,
        }
    }

    #[test]
    fn empty_input_selects_nothing() {
        assert!(solve(&[]).is_empty());
    }

    #[test]
    fn support_of_one_has_zero_weight_and_is_not_selected() {
        // w - 1 == 0, so the interval is worth nothing and the strict `>` in the
        // traceback rejects it. This is the "all intervals had only the read
        // itself as support" case the caller checks for.
        assert!(solve(&[iv(0, 100, 1)]).is_empty());
        assert!(solve(&[iv(0, 50, 1), iv(60, 100, 1)]).is_empty());
    }

    #[test]
    fn a_single_supported_interval_is_selected() {
        assert_eq!(solve(&[iv(0, 100, 5)]), vec![0]);
    }

    /// The reference's WIS is **suboptimal**, and the port reproduces that.
    ///
    /// `fill_p2` stores 0-based interval indices in `p`, but `p[j]` is used to
    /// index the 1-based `OPT` array (`OPT[p[j]]`). Every predecessor is
    /// therefore off by one, so a compatible earlier interval is often treated
    /// as incompatible. Three pairwise-disjoint intervals — which a correct WIS
    /// would all select — come back as just two.
    ///
    /// Verified directly against `solve_WIS`: this returns `[2, 0]`, not
    /// `[2, 1, 0]`. Fixing it would change corrected output on real data, so it
    /// is a reference bug to preserve, not to repair. See PORTING.md.
    #[test]
    fn disjoint_intervals_reproduce_the_reference_off_by_one() {
        assert_eq!(
            solve(&[iv(0, 10, 3), iv(20, 30, 3), iv(40, 50, 3)]),
            vec![2, 0]
        );
        assert_eq!(
            solve(&[iv(0, 10, 4), iv(20, 30, 4), iv(40, 50, 4)]),
            vec![2, 0]
        );
    }

    #[test]
    fn five_overlapping_intervals_match_the_reference() {
        // Also verified directly against solve_WIS.
        let got = solve(&[
            iv(0, 30, 5),
            iv(10, 40, 8),
            iv(35, 60, 4),
            iv(50, 90, 9),
            iv(80, 120, 6),
        ]);
        assert_eq!(got, vec![3, 0]);
    }

    #[test]
    fn overlapping_intervals_pick_the_heavier_set() {
        // One wide well-supported interval versus a narrow one it covers.
        let got = solve(&[iv(10, 20, 2), iv(0, 100, 10)]);
        assert_eq!(got, vec![1], "the heavier interval wins");
    }

    #[test]
    fn indices_are_returned_in_decreasing_order() {
        let got = solve(&[iv(0, 10, 4), iv(20, 30, 4), iv(40, 50, 4)]);
        for pair in got.windows(2) {
            assert!(pair[0] > pair[1], "expected decreasing: {got:?}");
        }
    }

    #[test]
    fn selected_intervals_never_overlap() {
        // The off-by-one makes the selection conservative, never invalid: it can
        // drop compatible intervals but must not return overlapping ones.
        let intervals = vec![
            iv(0, 30, 5),
            iv(10, 40, 8),
            iv(35, 60, 4),
            iv(50, 90, 9),
            iv(80, 120, 6),
        ];
        let chosen = solve(&intervals);
        let mut picked: Vec<Interval> = chosen.iter().map(|&j| intervals[j]).collect();
        picked.sort_by_key(|i| i.start);
        for pair in picked.windows(2) {
            assert!(
                pair[0].stop <= pair[1].start,
                "overlap between {:?} and {:?}",
                pair[0],
                pair[1]
            );
        }
    }

    #[test]
    fn zero_length_intervals_are_handled() {
        // start == stop is excluded from stop_to_max_j by `if start < stop`,
        // but the interval still carries weight epsilon * (w - 1).
        let got = solve(&[iv(5, 5, 3), iv(0, 10, 3)]);
        assert!(!got.is_empty());
    }

    #[test]
    fn epsilon_makes_zero_length_intervals_worth_something() {
        // Without the epsilon a zero-length interval would weigh exactly 0 and
        // the strict `>` would drop it.
        assert_eq!(solve(&[iv(7, 7, 4)]), vec![0]);
    }
}

/// Replay recorded `solve_WIS` calls from `bench/dump_reference.py`.
///
/// Hand-written cases cannot cover the float weights and strict-`>` tie-breaks
/// that real data produces; this feeds the reference's own inputs back through
/// and compares the chosen indices exactly.
#[cfg(test)]
mod replay {
    use super::*;
    use std::collections::BTreeMap;

    /// Recorded calls: inputs per call id, and the indices the reference chose.
    type Recorded = (BTreeMap<usize, Vec<Interval>>, BTreeMap<usize, Vec<usize>>);

    fn load(dir: &str) -> Option<Recorded> {
        let inp = std::fs::read_to_string(format!("{dir}/wis_input.tsv")).ok()?;
        let out = std::fs::read_to_string(format!("{dir}/wis_output.tsv")).ok()?;

        let mut calls: BTreeMap<usize, Vec<Interval>> = BTreeMap::new();
        for line in inp.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
            let f: Vec<&str> = line.split('\t').collect();
            calls
                .entry(f[0].parse().unwrap())
                .or_default()
                .push(Interval {
                    start: f[2].parse().unwrap(),
                    stop: f[3].parse().unwrap(),
                    support: f[4].parse().unwrap(),
                });
        }
        let mut expected: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for line in out.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
            let f: Vec<&str> = line.split('\t').collect();
            expected
                .entry(f[0].parse().unwrap())
                .or_default()
                .push(f[1].parse().unwrap());
        }
        Some((calls, expected))
    }

    #[test]
    fn matches_recorded_reference_calls() {
        let dir = match std::env::var("WIS_DUMP") {
            Ok(d) => d,
            // Not every checkout has a dump; skip rather than fail spuriously.
            Err(_) => return,
        };
        let Some((calls, expected)) = load(&dir) else {
            panic!("WIS_DUMP={dir} set but no wis_input.tsv / wis_output.tsv there");
        };
        assert!(!calls.is_empty(), "no recorded calls in {dir}");

        let mut checked = 0usize;
        for (call, intervals) in &calls {
            let got = solve(intervals);
            let want = expected.get(call).cloned().unwrap_or_default();
            assert_eq!(
                got,
                want,
                "call {call} with {} intervals disagreed",
                intervals.len()
            );
            checked += 1;
        }
        eprintln!("replayed {checked} solve_WIS calls from {dir}");
    }
}
