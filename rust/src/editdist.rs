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

use triple_accel::levenshtein::levenshtein_simd_k;

/// Sentinel for "distance exceeds the bound", matching edlib's `-1`.
pub const OVER_BOUND: i64 = -1;

/// Global edit distance between `x` and `y`, or [`OVER_BOUND`] if it exceeds `k`.
///
/// `k` is a maximum edit distance, not a band width. The reference always passes
/// `len(ref_seq)`, so the bound only bites when the two segments differ wildly
/// in length.
pub fn bounded(x: &[u8], y: &[u8], k: usize) -> i64 {
    // triple_accel returns None when the distance exceeds k, which is exactly
    // edlib's -1 contract.
    match levenshtein_simd_k(x, y, k as u32) {
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
