//! Affine global alignment via `block-aligner`, shared by the two call sites
//! that were the port's top two costs.
//!
//! # Why affine, and why this is not "approximate"
//!
//! block-aligner asserts `gaps.open < gaps.extend`, so it cannot express
//! unit-cost edit distance. That looks like a blocker for [`crate::align`]'s call
//! site, which reproduces edlib. It is not: **edit distance was the reference's
//! choice for speed in Python, not a modelling preference.** Affine gaps describe
//! sequencing indels better, and the reference already uses them where it can
//! afford to — the guard's `match 4, mismatch -8, open 12, extend 1`. Using that
//! scoring here is a modelling improvement that happens to also be faster.
//!
//! It is a *divergence from byte-identity*, though, and a much deeper one than
//! the guard's: this alignment builds the MSA, which builds the PFM, which
//! decides every corrected base. `ISONCORRECT_EXACT_ALIGN=1` restores the exact
//! edlib-compatible path, and `bench/equivalence.sh` sets it, so the 29-case gate
//! still covers the whole pipeline.
//!
//! # The block width is derived, not tuned
//!
//! block-aligner grows its block greedily, following the best score. On a *global*
//! alignment of unequal-length sequences that heuristic can fail outright: the
//! path must deviate from the diagonal by at least `|len(a) - len(b)|`, and if the
//! block starts narrower than that the gap is never found. Measured on one real
//! cluster at a fixed block of 32: **10 of 14 831 alignments came back
//! suboptimal**, every one of them a ~43 bp segment against an ~86 bp consensus.
//!
//! So [`min_block`] sizes the starting block from the length difference. That is
//! a *lower bound on the deviation*, not a threshold read off a corpus — which is
//! precisely what the rejected banding scheme (see [`crate::parasail`]) lacked.
//! It costs about 15% of the raw speedup and removes the silent failure mode.
//!
//! Verified against [`crate::parasail::global_affine`] on **158 559 recorded
//! segment-vs-consensus alignments across five corpora: zero suboptimal.**
//! `align::block_option` is that test.

use crate::parasail::Scoring;

/// Largest block the aligner may grow to. Larger costs nothing when the
/// alignment stays near the diagonal.
const MAX_BLOCK: usize = 256;

/// Smallest block for a pair of these lengths: the length difference is a hard
/// lower bound on how far the path must leave the diagonal. See the module docs.
fn min_block(a: usize, b: usize) -> usize {
    let need = a.abs_diff(b) + 32;
    need.next_power_of_two().clamp(32, MAX_BLOCK)
}

/// Shortest sequence worth handing to SIMD.
///
/// Below this the exact DP is already trivial and block-aligner's minimum block
/// size would dominate, so there is nothing to gain.
pub const MIN_LEN: usize = 64;

/// Affine global alignment, returning an extended CIGAR.
///
/// Gap penalties translate directly: block-aligner charges
/// `open + extend * (n - 1)`, the convention parasail uses, and its `I`/`D` mean
/// the same thing (`I` consumes the query, `D` the reference).
///
/// Returns `None` when the aligner cannot be used or its result does not span
/// both sequences — a non-ACGTN character (block-aligner's `NucMatrix` has no
/// entry, and `convert_char` asserts `A..=Z`), a penalty outside `i8`, or a
/// clipped result. Every caller falls back to an exact DP, so `None` is "use the
/// slow path", never an error.
pub fn global_cigar(query: &[u8], reference: &[u8], sc: Scoring) -> Option<String> {
    use block_aligner::{cigar::*, scan_block::*, scores::*};

    if !usable(query) || !usable(reference) {
        return None;
    }
    let matrix = NucMatrix::new_simple(
        i8::try_from(sc.match_score).ok()?,
        i8::try_from(sc.mismatch).ok()?,
    );
    let gaps = Gaps {
        open: i8::try_from(-sc.open).ok()?,
        extend: i8::try_from(-sc.ext).ok()?,
    };

    let q = PaddedBytes::from_bytes::<NucMatrix>(query, MAX_BLOCK);
    let r = PaddedBytes::from_bytes::<NucMatrix>(reference, MAX_BLOCK);
    let mut block = Block::<true, false>::new(query.len(), reference.len(), MAX_BLOCK);
    let lo = min_block(query.len(), reference.len());
    block.align(&q, &r, &matrix, gaps, lo..=MAX_BLOCK, 0);

    let res = block.res();
    if res.query_idx != query.len() || res.reference_idx != reference.len() {
        return None;
    }
    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    block
        .trace()
        .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
    Some(cigar.to_string())
}

/// Whether `block-aligner`'s nucleotide matrix covers every character.
///
/// `msa::min_ed` aligns a gap-padded insertion slot like `"-TT-"`, and
/// `convert_char` asserts `A..=Z`, so that call site has to stay on the exact DP.
/// It is a small minority — 27 of 677 `task="path"` calls at default parameters.
fn usable(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'N'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_starting_block_covers_the_length_difference() {
        // Equal lengths need nothing beyond the floor.
        assert_eq!(min_block(80, 80), 32);
        // 43 against 86 is the shape that failed at a fixed 32.
        assert!(min_block(43, 86) >= 43);
        // And it never exceeds the ceiling.
        assert_eq!(min_block(10, 5000), MAX_BLOCK);
    }

    #[test]
    fn a_gap_padded_slot_is_rejected_rather_than_asserting() {
        // `NucMatrix` has no entry for `-`, and block-aligner would assert
        // inside `convert_char`. Returning None sends the caller to the exact DP.
        assert!(global_cigar(b"-TT-", b"ATTA", Scoring::GUARD).is_none());
        assert!(global_cigar(b"ACGT", b"ACRT", Scoring::GUARD).is_none());
    }

    #[test]
    fn an_identical_pair_aligns_end_to_end() {
        let s = b"ACGTACGTACGTTTTTACGTACGTACGTGGGGACGTACGTACGTACGT";
        assert_eq!(
            global_cigar(s, s, Scoring::GUARD).as_deref(),
            Some(format!("{}=", s.len()).as_str())
        );
    }
}
