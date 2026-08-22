//! Bounded global edit distance, matching `edlib_alignment`.
//!
//! The reference calls
//! `edlib.align(x, y, "NW", 'dist', k)` and takes `result["editDistance"]`,
//! which is the **global (Needleman-Wunsch) Levenshtein distance**, or `-1`
//! when it exceeds `k`.
//!
//! # Why this can be native Rust, and the CIGAR path cannot
//!
//! isONcorrect uses edlib in two ways, and they have different requirements:
//!
//! * here, only the **distance** is used. Edit distance is a uniquely defined
//!   number, so any correct implementation agrees with edlib by definition.
//!   Nothing about edlib's internals leaks into the result.
//! * in `get_best_corrections`, edlib is called with `task="path"` and the
//!   **CIGAR** is consumed. That is *not* unique --- many alignments achieve the
//!   same optimal distance and edlib picks one by its own tie-breaking. A
//!   different implementation would return a different, equally optimal
//!   alignment, changing the MSA and therefore the corrected output.
//!
//! So this module is safe to implement natively; the consensus stage is not.
//! See PORTING.md.
//!
//! # Why `rapidfuzz` and not `triple_accel`
//!
//! This was the port's single hottest function --- **67% of runtime** on a
//! profiled 200-read cluster --- and the reason was not the algorithm but the
//! architecture: `triple_accel` guards its SIMD paths with `target_arch = "x86"`
//! and `"x86_64"` and has no `aarch64` equivalent, so on Apple Silicon or any
//! ARM host every call fell back to naive `O(n*m)` scalar DP.
//!
//! `rapidfuzz` implements the same bit-parallel Myers algorithm edlib uses, in
//! plain `u64` arithmetic with **no architecture guards at all**, so it is fast
//! everywhere rather than fast on x86. `O(n * m / 64)` instead of `O(n * m)`.
//!
//! # Measured and rejected: a cached pattern
//!
//! `find_most_supported_span` compares one `ref_seq` against every read carrying
//! an anchor pair, so hoisting the bit-parallel pattern table out of the loop
//! with `levenshtein::BatchComparator` looks free. It is a **loss**: 4.11 s to
//! 7.51 s over 15.4 million calls on three real clusters.
//!
//! The reason is that `distance_with_args` has a specialised path for patterns
//! that fit one 64-bit block, using a stack `PatternMatchVector`, while
//! `BatchComparator` always builds the heap `BlockPatternMatchVector`. Segments
//! here are ~80 bp, so the plain call is on the fast path and the "optimisation"
//! moves it off. Do not re-try this without re-measuring.
//!
//! Swapping the implementation is safe precisely because of the argument above:
//! the answer is a uniquely defined number, so there is no tie-break to
//! preserve. The differential tests below check it anyway --- against a naive DP
//! over exhaustive short inputs and random longer ones, and against Python
//! edlib over recorded real pairs.

use rapidfuzz::distance::levenshtein;

/// Sentinel for "distance exceeds the bound", matching edlib's `-1`.
pub const OVER_BOUND: i64 = -1;

/// Global edit distance between `x` and `y`, or [`OVER_BOUND`] if it exceeds `k`.
///
/// `k` is a maximum edit distance, not a band width. The reference always passes
/// `len(ref_seq)`, so the bound only bites when the two segments differ wildly
/// in length.
pub fn bounded(x: &[u8], y: &[u8], k: usize) -> i64 {
    // `score_cutoff` yields None when the distance exceeds k, which is exactly
    // edlib's -1 contract.
    let args = levenshtein::Args::default().score_cutoff(k);
    match levenshtein::distance_with_args(x.iter().copied(), y.iter().copied(), &args) {
        Some(d) => d as i64,
        None => OVER_BOUND,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_have_distance_zero() {
        assert_eq!(bounded(b"ACGTACGT", b"ACGTACGT", 8), 0);
    }

    #[test]
    fn empty_against_empty_is_zero() {
        assert_eq!(bounded(b"", b"", 0), 0);
    }

    #[test]
    fn substitutions_insertions_and_deletions_all_cost_one() {
        assert_eq!(bounded(b"ACGT", b"AGGT", 4), 1, "substitution");
        assert_eq!(bounded(b"ACGT", b"ACGGT", 5), 1, "insertion");
        assert_eq!(bounded(b"ACGT", b"AGT", 4), 1, "deletion");
    }

    #[test]
    fn distance_is_symmetric() {
        let (a, b) = (b"ACGTACGTAGGCTA".as_slice(), b"ACGTCGTAGGCTAA".as_slice());
        assert_eq!(bounded(a, b, 20), bounded(b, a, 20));
    }

    #[test]
    fn exceeding_the_bound_returns_the_edlib_sentinel() {
        // Distance is 4; a bound of 2 must report "over" rather than a number.
        assert_eq!(bounded(b"AAAA", b"CCCC", 2), OVER_BOUND);
        // At the bound it is reported normally.
        assert_eq!(bounded(b"AAAA", b"CCCC", 4), 4);
    }

    #[test]
    fn empty_against_nonempty_is_the_length() {
        assert_eq!(bounded(b"", b"ACGT", 4), 4);
        assert_eq!(bounded(b"ACGT", b"", 4), 4);
        // ...and is bounded like anything else.
        assert_eq!(bounded(b"", b"ACGT", 3), OVER_BOUND);
    }

    /// Textbook full-matrix Levenshtein. Slow and obviously correct, which is
    /// the point: it is the reference the fast implementation is checked against.
    fn naive(x: &[u8], y: &[u8]) -> usize {
        let mut prev: Vec<usize> = (0..=y.len()).collect();
        let mut cur = vec![0usize; y.len() + 1];
        for i in 1..=x.len() {
            cur[0] = i;
            for j in 1..=y.len() {
                let sub = prev[j - 1] + usize::from(x[i - 1] != y[j - 1]);
                cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        prev[y.len()]
    }

    /// Every pair of strings over `AC` up to length 6, against the naive DP.
    ///
    /// Exhaustive rather than sampled, because bit-parallel edit distance fails
    /// in narrow structural cases --- word boundaries, all-match rows, empty
    /// inputs --- not in "typical" ones.
    #[test]
    fn exhaustive_short_inputs_match_a_naive_dp() {
        let mut words: Vec<Vec<u8>> = Vec::new();
        for len in 0..=6usize {
            for bits in 0..(1u32 << len) {
                words.push(
                    (0..len)
                        .map(|i| if bits >> i & 1 == 1 { b'A' } else { b'C' })
                        .collect(),
                );
            }
        }
        let mut checked = 0usize;
        for a in &words {
            for b in &words {
                let want = naive(a, b);
                assert_eq!(
                    bounded(a, b, a.len().max(b.len())),
                    want as i64,
                    "on {:?} vs {:?}",
                    String::from_utf8_lossy(a),
                    String::from_utf8_lossy(b)
                );
                checked += 1;
            }
        }
        assert!(checked > 15_000, "expected a real sweep, got {checked}");
    }

    /// Longer random sequences, including lengths that straddle the 64-bit word
    /// boundary the bit-parallel algorithm blocks on.
    #[test]
    fn random_inputs_across_the_word_boundary_match_a_naive_dp() {
        // Deterministic xorshift: a fixed corpus beats a random one that cannot
        // be re-run.
        let mut state = 0x243f_6a88_85a3_08d3u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let alphabet = b"ACGT";
        for _ in 0..3000 {
            // Lengths around 64 and 128 are where blocking bugs live.
            let pick = |n: u64| (n % 200) as usize;
            let (la, lb) = (pick(next()), pick(next()));
            let a: Vec<u8> = (0..la).map(|_| alphabet[(next() % 4) as usize]).collect();
            let b: Vec<u8> = (0..lb).map(|_| alphabet[(next() % 4) as usize]).collect();
            let want = naive(&a, &b) as i64;
            assert_eq!(bounded(&a, &b, la.max(lb)), want, "len {la} vs {lb}");
            // And the bound behaves: one below the answer must report "over".
            if want > 0 {
                assert_eq!(bounded(&a, &b, (want - 1) as usize), OVER_BOUND);
                assert_eq!(bounded(&a, &b, want as usize), want);
            }
        }
    }

    #[test]
    fn lengths_at_and_beyond_one_word_are_exact() {
        for len in [63usize, 64, 65, 127, 128, 129] {
            let a = vec![b'A'; len];
            let mut b = a.clone();
            // One substitution in the middle, one at each end.
            b[0] = b'C';
            b[len / 2] = b'C';
            b[len - 1] = b'C';
            assert_eq!(bounded(&a, &b, len), 3, "at length {len}");
        }
    }

    #[test]
    fn is_global_not_local_or_infix() {
        // A local aligner would score this 0 by matching the shared core; a
        // global one must pay for the flanks. The reference uses mode "NW".
        assert_eq!(bounded(b"ACGT", b"TTTTACGTTTTT", 20), 8);
    }
}

/// Parity against the reference's edlib, replayed from recorded pairs.
///
/// Generate the cases inside the reference environment, then run the test:
///
/// ```bash
/// bench/gen_editdist_cases.py --fastq test_data/isoncorrect/0.fastq --out /tmp/ed.tsv
/// EDLIB_CASES=/tmp/ed.tsv cargo test --manifest-path rust/Cargo.toml editdist::parity -- --nocapture
/// ```
#[cfg(test)]
mod parity {
    use super::*;

    #[test]
    fn agrees_with_python_edlib_on_real_pairs() {
        let Ok(path) = std::env::var("EDLIB_CASES") else {
            return; // no recorded cases in this checkout
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("EDLIB_CASES={path} unreadable: {e}"));

        let (mut n, mut bad) = (0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                continue;
            }
            let got = bounded(f[0].as_bytes(), f[1].as_bytes(), f[2].parse().unwrap());
            let want: i64 = f[3].parse().unwrap();
            if got != want {
                if bad < 3 {
                    eprintln!(
                        "mismatch: len({}, {}) k={} want={want} got={got}",
                        f[0].len(),
                        f[1].len(),
                        f[2]
                    );
                }
                bad += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad, 0, "{bad} of {n} pairs disagreed with edlib");
        eprintln!("compared {n} real pairs against Python edlib, all equal");
    }
}
