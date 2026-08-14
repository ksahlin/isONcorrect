//! Sequence contexts and alternative references, matching `get_contexts`,
//! `test_numba` and `get_alternative_ref_contexts`.
//!
//! The PFM alone cannot tell a sequencing error from a real variant: both show
//! up as a minority base. What separates them is the **surrounding sequence**.
//! An error is independent of its neighbours, so the reads carrying it differ
//! from the consensus at that one column and agree elsewhere; a variant travels
//! with the rest of its allele, so the reads carrying it share a whole context.
//!
//! So each column gets a window of roughly `k/2` non-gap characters either
//! side, the distinct context strings across all rows are counted, and a
//! context that is both frequent enough and far enough from the consensus
//! becomes an **alternative reference** — reads matching it are corrected
//! towards it instead of towards the consensus.
//!
//! # What decides the answer
//!
//! * the window is measured in **non-gap characters**, so it widens across
//!   insertion columns. It is also asymmetric: the left half stops on the
//!   `k/2 + 1`-th character, the right half goes one further (`+ 1`);
//! * contexts are counted over **every row including the consensus**
//!   (`A = np.array([l for l in alignment_matrix.values()])` keeps `"ref"`);
//! * the eligible contexts are ordered **by depth descending, ties by the
//!   context string ascending** — the tie-break comes from `np.unique`, which
//!   sorts the raw bytes, and `sorted(..., reverse=True)` being stable. Verified
//!   against numpy, not assumed;
//! * `eligible_contexts_hcomp` is **iterated while being appended to**, and the
//!   nearest previous entry wins on a strict `<`, so its insertion order
//!   selects which alternative a context is filed under. Hence [`IndexMap`];
//! * both the raw context and its homopolymer-compressed form compete for that
//!   minimum, and whichever wins becomes the key.
//!
//! # The reference's own set-iteration order
//!
//! `get_alternative_ref_contexts` returns a **`set`** per column, and
//! `get_best_corrections` iterates it with a `break` on an exact match and a
//! strict `<` on edit distance. Set iteration order over tuples of strings
//! depends on `PYTHONHASHSEED`, so the reference is order-dependent here by
//! construction. Measured: columns with two and three alternatives do occur
//! (15 and 2 respectively on a 342-read SIRV cluster), but corrected output was
//! byte-identical across three hash seeds — the alternatives have to *tie* on
//! edit distance before the order can bite. This port emits insertion order,
//! which is deterministic; see PORTING.md.

use indexmap::IndexMap;

use crate::editdist;

/// Offset of the `k + 1`-th non-gap character, matching `get_context_offset`.
///
/// Returns the index at which the count exceeds `k`, or the last index if it
/// never does — and `0` for an empty input, which falls out of the loop rather
/// than needing the reference's explicit guard.
pub fn context_offset<I: Iterator<Item = u8>>(vector: I, k: usize) -> usize {
    let mut nuc_obs = 0usize;
    let mut offset = 0usize;
    for (i, n) in vector.enumerate() {
        offset = i;
        if n != b'-' {
            nuc_obs += 1;
        }
        if nuc_obs > k {
            break;
        }
    }
    offset
}

/// Context window `[start, stop)` per alignment column, matching
/// `get_contexts(alignment_matrix, k_size)`.
///
/// `k_size` is `int(k / 2)` at the call site. The window is **not symmetric**:
/// the right edge takes `get_context_offset(...) + 1`, so it reaches one column
/// further than the left. `stop` can equal `ref_aln.len()`; `start` is always a
/// valid index.
pub fn contexts(ref_aln: &[u8], k_size: usize) -> Vec<(usize, usize)> {
    (0..ref_aln.len())
        .map(|i| {
            let p1 = context_offset(ref_aln[..i].iter().rev().copied(), k_size);
            let p2 = context_offset(ref_aln[i + 1..].iter().copied(), k_size) + 1;
            (i - p1, i + p2)
        })
        .collect()
}

/// Collapse homopolymer runs: `AACCGT` becomes `ACGT`.
///
/// `''.join(ch for ch, _ in itertools.groupby(seq))`.
pub fn homopolymer_compress(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if out.last() != Some(&c) {
            out.push(c);
        }
    }
    out
}

/// One distinct context in a window, with the number of rows carrying it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Eligible<'a> {
    /// The character this column contributes, i.e. `context[i - start]`.
    pub variant: u8,
    /// The context as it appears in the matrix — gaps included.
    pub context: &'a [u8],
    pub depth: u32,
}

/// Distinct contexts of one window, ordered as the reference produces them.
#[derive(Debug, Clone)]
struct Window {
    start: usize,
    contexts: Vec<(Vec<u8>, u32)>,
}

/// Frequency context matrix: the eligible contexts of every column.
///
/// The reference's `test_numba` reuses the previous column's entries whenever
/// the window is unchanged, only shifting which character counts as the
/// variant. This stores that sharing directly — one [`Window`] per distinct
/// window, with columns pointing at it — rather than copying the tuples, which
/// is the same data in far less space. Consecutive columns share a window
/// often, since the window only moves when a non-gap character enters or leaves
/// it.
#[derive(Debug, Clone)]
pub struct Fcm {
    windows: Vec<Window>,
    per_column: Vec<usize>,
}

impl Fcm {
    /// Eligible contexts of column `i`, in the reference's order.
    pub fn column(&self, i: usize) -> impl Iterator<Item = Eligible<'_>> + '_ {
        let w = &self.windows[self.per_column[i]];
        w.contexts.iter().map(move |(ctx, depth)| Eligible {
            variant: ctx[i - w.start],
            context: ctx,
            depth: *depth,
        })
    }

    pub fn len(&self) -> usize {
        self.per_column.len()
    }

    pub fn is_empty(&self) -> bool {
        self.per_column.is_empty()
    }
}

/// Build the frequency context matrix, matching `test_numba`.
///
/// `rows` is the whole multialignment matrix **including the consensus row** —
/// the reference builds `A` from `alignment_matrix.values()`, which keeps
/// `"ref"`, so the consensus contributes one count to its own context.
///
/// A context is eligible when its depth is strictly greater than
/// `max(context_threshold / 10, 5)`. The reference expresses that as a `break`
/// out of a depth-descending walk, which is the same thing.
pub fn frequency_context_matrix(
    rows: &[Vec<u8>],
    contexts: &[(usize, usize)],
    context_threshold: f64,
) -> Fcm {
    let cutoff = (context_threshold / 10.0).max(5.0);
    let mut windows: Vec<Window> = Vec::new();
    let mut per_column = Vec::with_capacity(contexts.len());
    let mut prev: Option<(usize, usize)> = None;

    for &(start, stop) in contexts {
        if prev == Some((start, stop)) {
            per_column.push(windows.len() - 1);
            continue;
        }
        prev = Some((start, stop));

        // `np.unique` sorts the raw bytes, which for these strings is plain
        // lexicographic order --- so a BTreeMap gives the tie-break for free.
        let mut counts: std::collections::BTreeMap<&[u8], u32> = std::collections::BTreeMap::new();
        for row in rows {
            *counts.entry(&row[start..stop]).or_insert(0) += 1;
        }

        let mut contexts: Vec<(Vec<u8>, u32)> = counts
            .into_iter()
            .filter(|&(_, n)| f64::from(n) > cutoff)
            .map(|(ctx, n)| (ctx.to_vec(), n))
            .collect();
        // Stable, so equal depths keep the lexicographic order above.
        contexts.sort_by_key(|&(_, n)| std::cmp::Reverse(n));

        windows.push(Window { start, contexts });
        per_column.push(windows.len() - 1);
    }

    Fcm {
        windows,
        per_column,
    }
}

/// An alternative reference for one column: correct reads matching `context`
/// towards `variant` rather than towards the consensus.
#[derive(Debug, Clone, PartialEq)]
pub struct Alternative {
    pub variant: u8,
    /// The context as it appears in the matrix — gaps included, so it compares
    /// directly against a row's slice.
    pub context: Vec<u8>,
    pub depth: u32,
    pub threshold: f64,
}

/// Alternative references per column, matching `get_alternative_ref_contexts`.
///
/// A column is skipped when it has no eligible contexts, or when every eligible
/// context contributes the same character — nothing to choose between.
///
/// Otherwise each context is filed under the nearest existing entry, starting
/// from the consensus context and its homopolymer-compressed form. A context
/// whose compressed form is already filed is dropped as a duplicate; one that
/// survives is admitted only if its depth clears `context_threshold / min_ed`,
/// so a context far from everything seen needs less support than a near one.
///
/// The two consensus entries are removed at the end: they are seeded to catch
/// contexts that are really the consensus, not to be reported as alternatives.
pub fn alternative_ref_contexts(
    fcm: &Fcm,
    ref_aln: &[u8],
    contexts: &[(usize, usize)],
    context_threshold: f64,
) -> Vec<Vec<Alternative>> {
    let mut out = vec![Vec::new(); ref_aln.len()];

    for j in 0..ref_aln.len() {
        let eligible: Vec<Eligible> = fcm.column(j).collect();
        // `len({v for ...}) == 1` --- one variant means nothing to choose.
        if eligible.is_empty() || eligible.iter().all(|e| e.variant == eligible[0].variant) {
            continue;
        }

        let (start, stop) = contexts[j];
        let consensus_context: Vec<u8> = degap(&ref_aln[start..stop]);
        let consensus_hpol = homopolymer_compress(&consensus_context);

        // `None` marks the two seeded consensus entries, which are removed
        // below; the reference seeds them with an empty `set()`.
        let mut filed: IndexMap<Vec<u8>, Option<Alternative>> = IndexMap::new();
        filed.insert(consensus_hpol.clone(), None);
        filed.insert(consensus_context.clone(), None);

        // Already depth-descending; the reference's re-sort is stable and so a
        // no-op, but the order is what decides which entries exist when a later
        // context looks for its nearest neighbour.
        for e in &eligible {
            let seq = degap(e.context);
            let seq_hpol = homopolymer_compress(&seq);
            if filed.contains_key(&seq_hpol) {
                continue;
            }

            // Nearest already-filed entry, with the raw and compressed forms
            // competing. Strict `<`, so the earliest-filed winner keeps it.
            let mut min_ed = 1000i64;
            let mut new_alt_ref: Vec<u8> = Vec::new();
            for prev in filed.keys() {
                let ed = editdist::bounded(prev, &seq, 1000);
                if ed < min_ed {
                    min_ed = ed;
                    new_alt_ref = seq.clone();
                }
                let ed = editdist::bounded(prev, &seq_hpol, 1000);
                if ed < min_ed {
                    min_ed = ed;
                    new_alt_ref = seq_hpol.clone();
                }
            }

            let threshold = context_threshold / (min_ed as f64).max(1.0);
            if f64::from(e.depth) >= threshold {
                // Insert, not replace-in-place: a repeated key keeps its
                // original position, which is Python dict semantics too.
                filed.insert(
                    new_alt_ref,
                    Some(Alternative {
                        variant: e.variant,
                        context: e.context.to_vec(),
                        depth: e.depth,
                        threshold,
                    }),
                );
            }
        }

        filed.shift_remove(&consensus_hpol);
        if consensus_hpol != consensus_context {
            filed.shift_remove(&consensus_context);
        }
        out[j] = filed.into_values().flatten().collect();
    }

    out
}

fn degap(v: &[u8]) -> Vec<u8> {
    v.iter().copied().filter(|&c| c != b'-').collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx(aln: &str, k: usize) -> Vec<(usize, usize)> {
        contexts(aln.as_bytes(), k)
    }

    fn rows(rs: &[&str]) -> Vec<Vec<u8>> {
        rs.iter().map(|r| r.as_bytes().to_vec()).collect()
    }

    #[test]
    fn the_offset_counts_non_gap_characters_only() {
        // k = 2: stop once a third character is seen, at index 2.
        assert_eq!(context_offset(b"ACGT".iter().copied(), 2), 2);
        // Gaps do not count, so the same k reaches further.
        assert_eq!(context_offset(b"A-C-G-T".iter().copied(), 2), 4);
    }

    #[test]
    fn a_short_vector_gives_its_last_index() {
        assert_eq!(context_offset(b"AC".iter().copied(), 9), 1);
        assert_eq!(context_offset(b"".iter().copied(), 9), 0);
    }

    #[test]
    fn the_window_is_asymmetric_by_one() {
        // Away from the edges the right side reaches one column further.
        let c = ctx("ACGTACGTACGT", 2);
        let (start, stop) = c[6];
        assert_eq!(6 - start, 2, "left reaches the third character back");
        assert_eq!(stop - 6, 3, "right reaches one further");
    }

    #[test]
    fn windows_stay_inside_the_alignment() {
        let aln = "AC-GTA-CGT";
        for k in [0, 1, 4, 100] {
            for (i, &(start, stop)) in ctx(aln, k).iter().enumerate() {
                assert!(start <= i && i < stop, "column {i} outside its own window");
                assert!(stop <= aln.len(), "window past the end at {i}");
            }
        }
    }

    #[test]
    fn a_gap_run_widens_the_window() {
        let plain = ctx("ACGTACGT", 1);
        let gapped = ctx("AC--GTACGT", 1);
        let (s, e) = gapped[5]; // the G after the gaps
        assert!(e - s > plain[3].1 - plain[3].0);
    }

    #[test]
    fn homopolymers_collapse_to_one_character_each() {
        assert_eq!(homopolymer_compress(b"AACCCGT"), b"ACGT");
        assert_eq!(homopolymer_compress(b"AAAA"), b"A");
        assert_eq!(homopolymer_compress(b""), b"");
        assert_eq!(homopolymer_compress(b"ACGT"), b"ACGT");
    }

    #[test]
    fn contexts_below_the_depth_cutoff_are_dropped() {
        // Cutoff is max(threshold/10, 5), so six identical rows clear it and
        // five do not.
        let m = rows(&["ACGTACG"; 6]);
        let c = ctx("ACGTACG", 1);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        assert_eq!(fcm.column(3).count(), 1);

        let m = rows(&["ACGTACG"; 5]);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        assert_eq!(fcm.column(3).count(), 0);
    }

    #[test]
    fn eligible_contexts_come_back_deepest_first() {
        let mut m = rows(&["ACGTACG"; 10]);
        m.extend(rows(&["ACGAACG"; 6]));
        let c = ctx("ACGTACG", 1);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        let got: Vec<(u32, u8)> = fcm.column(3).map(|e| (e.depth, e.variant)).collect();
        assert_eq!(got, vec![(10, b'T'), (6, b'A')]);
    }

    #[test]
    fn equal_depths_break_ties_lexicographically() {
        let mut m = rows(&["ACGTACG"; 8]);
        m.extend(rows(&["ACGAACG"; 8]));
        let c = ctx("ACGTACG", 1);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        let got: Vec<u8> = fcm.column(3).map(|e| e.variant).collect();
        assert_eq!(got, vec![b'A', b'T'], "ACGA sorts before ACGT");
    }

    #[test]
    fn the_consensus_row_is_counted_like_any_other() {
        // Nine supporting rows plus the consensus makes ten.
        let mut m = rows(&["ACGTACG"]); // the consensus row
        m.extend(rows(&["ACGTACG"; 9]));
        let c = ctx("ACGTACG", 1);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        assert_eq!(fcm.column(3).next().unwrap().depth, 10);
    }

    #[test]
    fn a_column_with_one_variant_has_no_alternatives() {
        let m = rows(&["ACGTACG"; 10]);
        let c = ctx("ACGTACG", 1);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        let alts = alternative_ref_contexts(&fcm, b"ACGTACG", &c, 3.0);
        assert!(alts.iter().all(Vec::is_empty));
    }

    #[test]
    fn a_deep_divergent_context_becomes_an_alternative() {
        // Half the rows carry a different allele across the whole window.
        let mut m = rows(&["ACGTACGTACGT"; 10]);
        m.extend(rows(&["ACGTTGATACGT"; 10]));
        let c = ctx("ACGTACGTACGT", 2);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        let alts = alternative_ref_contexts(&fcm, b"ACGTACGTACGT", &c, 3.0);
        let with_alts: Vec<usize> = alts
            .iter()
            .enumerate()
            .filter(|(_, a)| !a.is_empty())
            .map(|(i, _)| i)
            .collect();
        assert!(!with_alts.is_empty(), "the divergent allele went unnoticed");
        for i in with_alts {
            assert_eq!(alts[i].len(), 1);
            assert_eq!(alts[i][0].depth, 10);
            // The alternative reports the column's own character, not the
            // consensus one.
            assert_eq!(alts[i][0].variant, b"ACGTTGATACGT"[i]);
        }
    }

    #[test]
    fn the_consensus_context_is_never_itself_an_alternative() {
        let mut m = rows(&["ACGTACGTACGT"; 10]);
        m.extend(rows(&["ACGTTGATACGT"; 10]));
        let c = ctx("ACGTACGTACGT", 2);
        let fcm = frequency_context_matrix(&m, &c, 3.0);
        let alts = alternative_ref_contexts(&fcm, b"ACGTACGTACGT", &c, 3.0);
        for (j, column) in alts.iter().enumerate() {
            let (start, stop) = c[j];
            let consensus = &b"ACGTACGTACGT"[start..stop];
            assert!(
                column.iter().all(|a| a.context != consensus),
                "column {j} reported the consensus as an alternative"
            );
        }
    }

    #[test]
    fn a_shallow_divergent_context_is_rejected_by_the_threshold() {
        // One row differs: depth 1 against a threshold of 1000/min_ed.
        let mut m = rows(&["ACGTACGTACGT"; 10]);
        m.extend(rows(&["ACGTTGATACGT"; 6]));
        let c = ctx("ACGTACGTACGT", 2);
        let fcm = frequency_context_matrix(&m, &c, 1000.0);
        let alts = alternative_ref_contexts(&fcm, b"ACGTACGTACGT", &c, 1000.0);
        assert!(alts.iter().all(Vec::is_empty));
    }
}

/// Differential oracle against `get_contexts` and `get_alternative_ref_contexts`.
///
/// ```bash
/// bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/d
/// CONTEXT_CASES=/tmp/d/context_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml contexts::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    /// One record per `get_alternative_ref_contexts` call:
    ///
    /// ```text
    /// C <case> <k_size> <context_threshold> <ref_aln>
    /// R <case> <row>                                -- matrix rows, ref first
    /// W <case> <start:stop,start:stop,...>          -- contexts_per_pos
    /// A <case> <col> <variant> <context> <depth> <threshold>
    /// ```
    ///
    /// Alternatives are compared as a **set**: the reference returns a Python
    /// `set` per column, whose iteration order depends on `PYTHONHASHSEED`. The
    /// port's order is its own choice, so ordering it here would test the
    /// harness's dump order rather than the algorithm. Everything else — which
    /// alternatives exist, their depths and thresholds — is compared exactly.
    ///
    /// Thresholds are compared as `f64`, after both sides go through Rust's
    /// formatter, so Python's `repr` conventions do not leak into the check.
    #[test]
    fn matches_the_reference_contexts_and_alternatives() {
        let Ok(path) = std::env::var("CONTEXT_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CONTEXT_CASES={path} unreadable: {e}"));

        #[derive(Default)]
        struct Case {
            k_size: usize,
            threshold: f64,
            ref_aln: Vec<u8>,
            rows: Vec<Vec<u8>>,
            windows: Vec<(usize, usize)>,
            alts: Vec<Vec<String>>,
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
                    threshold: f[3].parse().unwrap(),
                    ref_aln: f[4].as_bytes().to_vec(),
                    ..Default::default()
                }),
                "R" => cases
                    .last_mut()
                    .unwrap()
                    .rows
                    .push(f[2].as_bytes().to_vec()),
                "W" => {
                    let c = cases.last_mut().unwrap();
                    c.windows = f[2]
                        .split(',')
                        .map(|p| {
                            let (a, b) = p.split_once(':').expect("start:stop");
                            (a.parse().unwrap(), b.parse().unwrap())
                        })
                        .collect();
                }
                "A" => {
                    let c = cases.last_mut().unwrap();
                    let col: usize = f[2].parse().unwrap();
                    if c.alts.len() <= col {
                        c.alts.resize(col + 1, Vec::new());
                    }
                    let threshold: f64 = f[6].parse().unwrap();
                    c.alts[col].push(format!("{}\t{}\t{}\t{}", f[3], f[4], f[5], threshold));
                }
                other => panic!("unknown record {other}"),
            }
        }
        assert!(!cases.is_empty(), "no cases in {path}");

        let (mut bad_win, mut bad_alt, mut n_alt, mut shown) = (0usize, 0usize, 0usize, 0usize);
        for (i, case) in cases.iter().enumerate() {
            let got_windows = contexts(&case.ref_aln, case.k_size);
            if got_windows != case.windows {
                bad_win += 1;
                continue;
            }

            let fcm = frequency_context_matrix(&case.rows, &got_windows, case.threshold);
            let got = alternative_ref_contexts(&fcm, &case.ref_aln, &got_windows, case.threshold);

            for (col, column) in got.iter().enumerate() {
                let want: std::collections::BTreeSet<String> = case
                    .alts
                    .get(col)
                    .map(|v| v.iter().cloned().collect())
                    .unwrap_or_default();
                let mine: std::collections::BTreeSet<String> = column
                    .iter()
                    .map(|a| {
                        format!(
                            "{}\t{}\t{}\t{}",
                            a.variant as char,
                            String::from_utf8_lossy(&a.context),
                            a.depth,
                            a.threshold,
                        )
                    })
                    .collect();
                n_alt += want.len();
                if want != mine {
                    if shown < 3 {
                        shown += 1;
                        eprintln!("case {i} col {col}:\n  python: {want:?}\n  rust  : {mine:?}");
                    }
                    bad_alt += 1;
                }
            }
        }
        assert_eq!(
            bad_win,
            0,
            "{bad_win} of {} context vectors differed",
            cases.len()
        );
        assert_eq!(bad_alt, 0, "{bad_alt} columns of alternatives differed");
        eprintln!(
            "matched the reference on {} context vectors and {n_alt} alternatives",
            cases.len()
        );
    }
}
