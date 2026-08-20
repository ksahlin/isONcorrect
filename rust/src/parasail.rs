//! Semi-global affine alignment matching `parasail.sg_trace_scan_16`.
//!
//! This is the alignment behind the structural-overcorrection guard: the
//! corrected read is aligned back against the original, and [`crate::guard`]
//! reverts any indel run longer than 10 columns on the assumption it is real
//! structure (an exon difference) rather than an error.
//!
//! `correct_read` calls it as
//!
//! ```text
//! parasail_alignment(seq, corr, match_score=4, mismatch_penalty=-8,
//!                    opening_penalty=12, gap_ext=1)
//! ```
//!
//! — note 12, not the function's own default of 24.
//!
//! # What "sg" means, measured
//!
//! parasail's plain `sg` leaves gaps at **both ends of both sequences** free.
//! Verified against the library rather than read off documentation:
//!
//! ```text
//! sg("ACGTACGT", "TTTACGTACGTTTT") -> 3D8=3D, score 32 == 8 matches * 4
//! sg("TTTACGTACGTTTT", "ACGTACGT") -> 3I8=3I, score 32
//! ```
//!
//! Both end gaps contribute nothing to the score, and the CIGAR still spans
//! both sequences in full. A gap of length `L` inside the alignment costs
//! `open + (L - 1) * ext`, so a single-base gap costs `open` — also measured:
//! `sg("ACGTACGT", "ACGACGT")` scores 16, which is the free-end-gap alignment
//! tying exactly with `7 * 4 - 12` for the one-base-deletion alignment.
//!
//! That tie is the whole problem. Free end gaps plus affine costs leave many
//! alignments at the optimal score, and parasail returns one of them. Which one
//! changes where indel runs fall, which changes what the guard reverts, which
//! changes the output — so the choice has to be reproduced, not merely matched
//! on score.
//!
//! # CIGAR convention
//!
//! `I` consumes `s1` (the query), `D` consumes `s2` (the database), matching
//! what `cigar_to_seq(cigar, s1, s2)` then expects — the same convention as
//! edlib's, so [`crate::align::cigar_to_seq`] expands these unchanged.

use crate::align::CigarOp;

/// Scoring, in parasail's sign convention: `match_score` positive, `mismatch`
/// negative, and both gap penalties positive numbers that get subtracted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub open: i32,
    pub ext: i32,
}

impl Scoring {
    /// What `correct_read` passes: `match=4, mismatch=-8, open=12, ext=1`.
    pub const GUARD: Scoring = Scoring {
        match_score: 4,
        mismatch: -8,
        open: 12,
        ext: 1,
    };
}

/// Which way a tie is broken. Every field was determined by sweeping against
/// recorded parasail output; see the `sweep` test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TieBreak {
    /// Order to try the three predecessors of an `H` cell.
    pub h_order: [Step; 3],
    /// On a tie inside a gap chain, prefer opening a new gap over extending.
    pub prefer_open: bool,
    /// When several cells on the last row/column share the best score, scan the
    /// last column before the last row.
    pub column_first: bool,
    /// Keep the last equally-good end cell rather than the first.
    pub last_max: bool,
}

/// A predecessor of an `H` cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// Diagonal: a match or mismatch.
    Diagonal,
    /// From the gap-in-`s1` chain, emitting `D`.
    Delete,
    /// From the gap-in-`s2` chain, emitting `I`.
    Insert,
}

impl TieBreak {
    /// The ordering that reproduces `sg_trace_scan_16`, measured rather than
    /// guessed — see the `sweep` test for the numbers.
    ///
    /// `column_first` and `last_max` only matter when several cells on the last
    /// row or column tie for the best score, which never happens on the recorded
    /// corpus: the two sequences are a read and its own correction, so a
    /// full-length alignment always beats ending early. Both values are
    /// therefore **unpinned by evidence**, and the sweep reports all four
    /// combinations as equally perfect. If a corpus ever reaches such a tie,
    /// re-run the sweep before trusting these two.
    pub const PARASAIL: TieBreak = TieBreak {
        h_order: [Step::Diagonal, Step::Insert, Step::Delete],
        prefer_open: false,
        column_first: false,
        last_max: false,
    };

    /// Every combination, for the sweep that picks the real one.
    #[cfg(test)]
    fn all() -> Vec<TieBreak> {
        use Step::*;
        let orders = [
            [Diagonal, Delete, Insert],
            [Diagonal, Insert, Delete],
            [Delete, Diagonal, Insert],
            [Delete, Insert, Diagonal],
            [Insert, Diagonal, Delete],
            [Insert, Delete, Diagonal],
        ];
        let mut out = Vec::new();
        for h_order in orders {
            for prefer_open in [true, false] {
                for column_first in [true, false] {
                    for last_max in [true, false] {
                        out.push(TieBreak {
                            h_order,
                            prefer_open,
                            column_first,
                            last_max,
                        });
                    }
                }
            }
        }
        out
    }
}

/// Result of a semi-global alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub cigar: String,
    pub ops: Vec<CigarOp>,
}

const NEG: i32 = i32::MIN / 4;

/// Align `s1` against `s2`, reproducing `parasail.sg_trace_scan_16`.
///
/// `O(n*m)` time and memory. The guard runs this once per read against a
/// sequence of nearly the same length, so the table is the read length squared —
/// the one place in the port where that is a real cost, and the reason
/// *Deferred improvements* notes banding it.
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    semiglobal_with(s1, s2, sc, TieBreak::PARASAIL)
}

/// [`semiglobal`] with an explicit tie-break, for the sweep.
pub fn semiglobal_with(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> Alignment {
    let (n, m) = (s1.len(), s2.len());
    let w = m + 1;
    let at = |i: usize, j: usize| i * w + j;

    // H: best overall. E: ends with a gap in s1 (emits D). F: ends with a gap
    // in s2 (emits I).
    let mut h = vec![NEG; (n + 1) * w];
    let mut e = vec![NEG; (n + 1) * w];
    let mut f = vec![NEG; (n + 1) * w];

    // Free end gaps: an unaligned prefix of either sequence costs nothing, so
    // the whole first row and column start at zero — including the gap chains,
    // which is what keeps a *leading* gap free rather than charging `open`.
    for i in 0..=n {
        h[at(i, 0)] = 0;
        f[at(i, 0)] = 0;
    }
    for j in 0..=m {
        h[at(0, j)] = 0;
        e[at(0, j)] = 0;
    }

    for i in 1..=n {
        for j in 1..=m {
            let open_d = h[at(i, j - 1)] - sc.open;
            let ext_d = e[at(i, j - 1)] - sc.ext;
            e[at(i, j)] = open_d.max(ext_d);

            let open_i = h[at(i - 1, j)] - sc.open;
            let ext_i = f[at(i - 1, j)] - sc.ext;
            f[at(i, j)] = open_i.max(ext_i);

            let diag = h[at(i - 1, j - 1)]
                + if s1[i - 1] == s2[j - 1] {
                    sc.match_score
                } else {
                    sc.mismatch
                };
            h[at(i, j)] = diag.max(e[at(i, j)]).max(f[at(i, j)]);
        }
    }

    // Trailing gaps are free too, so the alignment may end anywhere on the last
    // row or the last column.
    let mut best: (i32, usize, usize) = (NEG, n, m);
    let mut ends: Vec<(i32, usize, usize)> = Vec::with_capacity(n + m + 2);
    let row = (0..=m).map(|j| (h[at(n, j)], n, j));
    let col = (0..=n).map(|i| (h[at(i, m)], i, m));
    if tb.column_first {
        ends.extend(col);
        ends.extend(row);
    } else {
        ends.extend(row);
        ends.extend(col);
    }
    for cand in ends {
        if cand.0 > best.0 || (tb.last_max && cand.0 == best.0) {
            best = cand;
        }
    }
    let (score, mut i, mut j) = best;

    // Traceback, emitting ops backwards.
    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    // Free trailing gaps first: whatever is left of either sequence.
    ops.extend(std::iter::repeat_n(b'D', m - j));
    ops.extend(std::iter::repeat_n(b'I', n - i));

    // `Some(Step::Delete)` means we are inside a D chain and must decide
    // open-vs-extend before leaving it; `None` means we are on H.
    let mut chain: Option<Step> = None;
    while i > 0 && j > 0 {
        match chain {
            Some(Step::Delete) => {
                let open_d = h[at(i, j - 1)] - sc.open;
                let ext_d = e[at(i, j - 1)] - sc.ext;
                let take_open = if tb.prefer_open {
                    open_d >= ext_d
                } else {
                    open_d > ext_d
                };
                ops.push(b'D');
                j -= 1;
                chain = if take_open { None } else { Some(Step::Delete) };
            }
            Some(Step::Insert) => {
                let open_i = h[at(i - 1, j)] - sc.open;
                let ext_i = f[at(i - 1, j)] - sc.ext;
                let take_open = if tb.prefer_open {
                    open_i >= ext_i
                } else {
                    open_i > ext_i
                };
                ops.push(b'I');
                i -= 1;
                chain = if take_open { None } else { Some(Step::Insert) };
            }
            Some(Step::Diagonal) => unreachable!("the diagonal is not a chain"),
            None => {
                let cur = h[at(i, j)];
                let same = s1[i - 1] == s2[j - 1];
                let diag = h[at(i - 1, j - 1)] + if same { sc.match_score } else { sc.mismatch };
                let mut taken = None;
                for step in tb.h_order {
                    let hit = match step {
                        Step::Diagonal => diag == cur,
                        Step::Delete => e[at(i, j)] == cur,
                        Step::Insert => f[at(i, j)] == cur,
                    };
                    if hit {
                        taken = Some(step);
                        break;
                    }
                }
                match taken.expect("some predecessor achieves the cell score") {
                    Step::Diagonal => {
                        ops.push(if same { b'=' } else { b'X' });
                        i -= 1;
                        j -= 1;
                    }
                    step => chain = Some(step),
                }
            }
        }
    }
    // Free leading gaps.
    ops.extend(std::iter::repeat_n(b'D', j));
    ops.extend(std::iter::repeat_n(b'I', i));

    ops.reverse();
    let (cigar, runs) = encode(&ops);
    Alignment {
        score,
        cigar,
        ops: runs,
    }
}

fn encode(ops: &[u8]) -> (String, Vec<CigarOp>) {
    let mut out = String::new();
    let mut runs = Vec::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
        runs.push((k - start, c));
    }
    (out, runs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sg(a: &str, b: &str) -> Alignment {
        semiglobal(a.as_bytes(), b.as_bytes(), Scoring::GUARD)
    }

    /// Ground truth read off the parasail library itself, not assumed.
    #[test]
    fn end_gaps_are_free_on_both_sequences() {
        let a = sg("ACGTACGT", "TTTACGTACGTTTT");
        assert_eq!(a.cigar, "3D8=3D");
        assert_eq!(a.score, 32);

        let b = sg("TTTACGTACGTTTT", "ACGTACGT");
        assert_eq!(b.cigar, "3I8=3I");
        assert_eq!(b.score, 32);
    }

    #[test]
    fn identical_sequences_score_every_match() {
        let a = sg("ACGTACGT", "ACGTACGT");
        assert_eq!(a.cigar, "8=");
        assert_eq!(a.score, 32);
    }

    #[test]
    fn a_single_base_gap_costs_the_opening_penalty() {
        // Both sides long enough that end gaps cannot win.
        let s1 = "ACGTACGTACGTTGCAACGT";
        let s2 = "ACGTACGTACGTGCAACGT"; // one T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 19 * 4 - 12);
        assert_eq!(a.cigar.matches('I').count(), 1);
    }

    #[test]
    fn a_longer_gap_costs_open_plus_extensions() {
        let s1 = "ACGTACGTACGTTTTTGCAACGTACGT";
        let s2 = "ACGTACGTACGTGCAACGTACGT"; // four T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 23 * 4 - (12 + 3));
    }

    #[test]
    fn the_cigar_spans_both_sequences() {
        for (s1, s2) in [
            ("ACGTACGT", "TTTACGTACGTTTT"),
            ("ACGTACGTACGTTGCAACGT", "ACGTACGTACGTGCAACGT"),
            ("ACGT", "ACGT"),
            ("", "ACGT"),
            ("ACGT", ""),
        ] {
            let a = semiglobal(s1.as_bytes(), s2.as_bytes(), Scoring::GUARD);
            let (mut q, mut r) = (0usize, 0usize);
            for (len, op) in &a.ops {
                match op {
                    b'=' | b'X' => {
                        q += len;
                        r += len;
                    }
                    b'I' => q += len,
                    b'D' => r += len,
                    _ => panic!("unexpected op"),
                }
            }
            assert_eq!((q, r), (s1.len(), s2.len()), "on {s1:?} vs {s2:?}");
        }
    }

    #[test]
    fn an_empty_sequence_is_all_free_gaps() {
        let a = sg("", "ACGT");
        assert_eq!((a.cigar.as_str(), a.score), ("4D", 0));
        let b = sg("ACGT", "");
        assert_eq!((b.cigar.as_str(), b.score), ("4I", 0));
        let c = sg("", "");
        assert_eq!((c.cigar.as_str(), c.score), ("", 0));
    }
}

/// Differential oracle against the parasail library.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml parasail::oracle -- --nocapture
/// ```
///
/// `sweep` is the tool that *chose* [`TieBreak::PARASAIL`]: it runs every
/// combination against the recorded cases and prints the score. Re-run it if a
/// mismatch ever appears — a new corpus reaching a tie the current ordering gets
/// wrong would show up as a different winner.
#[cfg(test)]
mod oracle {
    use super::*;

    struct Case {
        s1: Vec<u8>,
        s2: Vec<u8>,
        cigar: String,
        score: i32,
        scoring: Scoring,
    }

    fn load() -> Option<Vec<Case>> {
        let path = std::env::var("PARASAIL_CASES").ok()?;
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("PARASAIL_CASES={path} unreadable: {e}"));
        let mut cases = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 8, "short row");
            cases.push(Case {
                s1: f[0].as_bytes().to_vec(),
                s2: f[1].as_bytes().to_vec(),
                cigar: f[2].to_string(),
                score: f[3].parse().unwrap(),
                scoring: Scoring {
                    match_score: f[4].parse().unwrap(),
                    mismatch: f[5].parse().unwrap(),
                    open: f[6].parse().unwrap(),
                    ext: f[7].parse().unwrap(),
                },
            });
        }
        Some(cases)
    }

    #[test]
    fn matches_parasail_exactly() {
        let Some(cases) = load() else { return };
        assert!(!cases.is_empty(), "no cases");

        let (mut bad_score, mut bad_cigar, mut shown) = (0usize, 0usize, 0usize);
        for c in &cases {
            let got = semiglobal(&c.s1, &c.s2, c.scoring);
            if got.score != c.score {
                bad_score += 1;
            }
            if got.cigar != c.cigar {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "cigar mismatch on {}x{}:\n  parasail: {}\n  rust    : {}",
                        c.s1.len(),
                        c.s2.len(),
                        c.cigar,
                        got.cigar
                    );
                }
                bad_cigar += 1;
            }
        }
        eprintln!(
            "parasail oracle: {} cases, {bad_score} score mismatches, {bad_cigar} cigar mismatches",
            cases.len()
        );
        assert_eq!(bad_score, 0, "scores differed");
        assert_eq!(bad_cigar, 0, "CIGARs differed");
    }

    /// Not an assertion --- a measurement. Prints how each tie-break ordering
    /// scores so the right one can be picked from evidence.
    #[test]
    fn sweep() {
        let Some(mut cases) = load() else { return };
        let Ok(setting) = std::env::var("PARASAIL_SWEEP") else {
            return;
        };
        // 96 configurations times full-length reads is hours of quadratic DP, so
        // the variable doubles as a case cap: PARASAIL_SWEEP=200 sweeps the
        // first 200. Ties are what matter here, not volume.
        if let Ok(cap) = setting.parse::<usize>() {
            if cap > 1 {
                cases.truncate(cap);
            }
        }
        let mut results: Vec<(usize, TieBreak)> = TieBreak::all()
            .into_iter()
            .map(|tb| {
                let hits = cases
                    .iter()
                    .filter(|c| semiglobal_with(&c.s1, &c.s2, c.scoring, tb).cigar == c.cigar)
                    .count();
                (hits, tb)
            })
            .collect();
        results.sort_by_key(|(hits, _)| std::cmp::Reverse(*hits));
        eprintln!("tie-break sweep over {} cases:", cases.len());
        for (hits, tb) in results.iter().take(8) {
            eprintln!(
                "  {hits:>6} / {}  h_order={:?} open={} col_first={} last_max={}",
                cases.len(),
                tb.h_order,
                tb.prefer_open,
                tb.column_first,
                tb.last_max
            );
        }
    }
}
