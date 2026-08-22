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
/// **Measured, not assumed.** At 64 this sent 39% of segment-vs-consensus
/// alignments (1.68 million of 4.27 million on four real clusters) to the scalar
/// DP, since those segments are median ~55 bp. Lowering it to 32 moved them onto
/// the SIMD path and cut the fallback to 75 587 calls — the gap-padded `min_ed`
/// queries, which cannot use it at all.
///
/// Below 32 there is genuinely nothing to gain: that is block-aligner's minimum
/// block size, so the block already spans the whole alignment.
pub const MIN_LEN: usize = 32;

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
    let mut ops = Vec::new();
    global_ops(query, reference, sc, &mut ops).then(|| crate::align::encode_cigar(&ops))
}

/// The same alignment, reported as operations rather than a CIGAR string.
///
/// This is the form the aligner produces natively, and what
/// [`crate::align::ops_to_seq`] consumes, so the hot path never builds a string.
/// `out` is cleared first and is a caller-owned buffer so it can be reused.
/// Returns `false` in exactly the cases [`global_cigar`] returns `None`.
pub fn global_ops(
    query: &[u8],
    reference: &[u8],
    sc: Scoring,
    out: &mut Vec<crate::align::CigarOp>,
) -> bool {
    SCRATCH.with_borrow_mut(|scratch| scratch.align_into(query, reference, sc, out))
}

thread_local! {
    /// See [`Scratch`]. Thread-local rather than threaded through every caller:
    /// the aligner is called from deep inside `get_best_corrections`, and
    /// clusters are processed one per thread.
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::new());
}

/// Reusable allocations for the aligner.
///
/// **This is not a micro-optimisation.** The DP over a ~55x57 segment is a few
/// thousand cells, and there are **4.27 million calls per four real clusters**, so
/// per-call setup competes with the arithmetic. A fresh `Block`, two
/// `PaddedBytes` and a `Cigar` (which allocates and zeroes `q + r + 5` 16-byte
/// entries, ~2 KB here) per call cost **18.84 s of 4 clusters; reused, 11.54 s**.
///
/// `Block::new` takes *maximum* lengths and `align` then takes the actual
/// sequences, so one block sized to the largest pair seen so far serves every
/// smaller one. It only ever grows.
struct Scratch {
    block: block_aligner::scan_block::Block<true, false>,
    q: block_aligner::scan_block::PaddedBytes,
    r: block_aligner::scan_block::PaddedBytes,
    /// Reusable because `cigar_eq` clears it itself — `Cigar::clear` is
    /// `pub(crate)`, but the traceback calls it, so a long-lived `Cigar` is safe.
    /// `Cigar::new` allocates and zeroes `q + r + 5` 16-byte entries, which at
    /// these lengths is ~2 KB per call.
    cigar: block_aligner::cigar::Cigar,
    cap: usize,
}

/// Starting capacity, and the granularity it grows in. Round up so a run of
/// slightly-longer sequences does not reallocate on every call.
const CAP_STEP: usize = 256;

impl Scratch {
    fn new() -> Self {
        use block_aligner::{scan_block::*, scores::*};
        Self {
            block: Block::new(CAP_STEP, CAP_STEP, MAX_BLOCK),
            q: PaddedBytes::new::<NucMatrix>(CAP_STEP, MAX_BLOCK),
            r: PaddedBytes::new::<NucMatrix>(CAP_STEP, MAX_BLOCK),
            cigar: block_aligner::cigar::Cigar::new(CAP_STEP, CAP_STEP),
            cap: CAP_STEP,
        }
    }

    fn reserve(&mut self, len: usize) {
        use block_aligner::{scan_block::*, scores::*};
        if len <= self.cap {
            return;
        }
        self.cap = len.next_multiple_of(CAP_STEP);
        self.block = Block::new(self.cap, self.cap, MAX_BLOCK);
        self.q = PaddedBytes::new::<NucMatrix>(self.cap, MAX_BLOCK);
        self.r = PaddedBytes::new::<NucMatrix>(self.cap, MAX_BLOCK);
        self.cigar = block_aligner::cigar::Cigar::new(self.cap, self.cap);
    }

    fn align_into(
        &mut self,
        query: &[u8],
        reference: &[u8],
        sc: Scoring,
        out: &mut Vec<crate::align::CigarOp>,
    ) -> bool {
        use block_aligner::{cigar::Operation, scores::*};

        out.clear();
        if !usable(query) || !usable(reference) {
            return false;
        }
        let (Ok(m), Ok(x), Ok(open), Ok(extend)) = (
            i8::try_from(sc.match_score),
            i8::try_from(sc.mismatch),
            i8::try_from(-sc.open),
            i8::try_from(-sc.ext),
        ) else {
            return false;
        };
        let matrix = NucMatrix::new_simple(m, x);
        let gaps = Gaps { open, extend };

        self.reserve(query.len().max(reference.len()));
        self.q.set_bytes::<NucMatrix>(query, MAX_BLOCK);
        self.r.set_bytes::<NucMatrix>(reference, MAX_BLOCK);
        let lo = min_block(query.len(), reference.len());
        self.block
            .align(&self.q, &self.r, &matrix, gaps, lo..=MAX_BLOCK, 0);

        let res = self.block.res();
        if res.query_idx != query.len() || res.reference_idx != reference.len() {
            return false;
        }
        self.block.trace().cigar_eq(
            &self.q,
            &self.r,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        out.reserve(self.cigar.len());
        for i in 0..self.cigar.len() {
            let op_len = self.cigar.get(i);
            let op = match op_len.op {
                Operation::Eq => b'=',
                Operation::X => b'X',
                Operation::I => b'I',
                Operation::D => b'D',
                // `cigar_eq` never emits M or Sentinel; treat either as "cannot
                // use this alignment" rather than guessing.
                _ => return false,
            };
            out.push((op_len.len, op));
        }
        true
    }
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
