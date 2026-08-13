//! Minimizer extraction, matching `get_kmer_minimizers_lex`.
//!
//! `hash_fcn` is hardcoded to `"lex"` in `isoncorrect_main`, so minimizers are
//! selected by **lexicographic order on the k-mer string**, not by a hash. The
//! random-hash and maximizer variants are dead in the default path.
//!
//! Three details in the reference are easy to get wrong and all of them change
//! output:
//!
//! 1. **Ties resolve to the last position in the window** (`rindex`), not the
//!    first.
//! 2. The rescan only happens when the outgoing k-mer *was* the minimizer *and*
//!    the minimizer has actually left the window (`minimizer_pos < i - w`).
//!    Comparing by value alone would rescan too eagerly.
//! 3. The loop is `range(w+1, len(seq) - k)`, whose upper bound is **exclusive**
//!    --- so the final k-mer at `len(seq) - k` is never examined. That is an
//!    off-by-one in the reference, and it is part of the contract.

use std::collections::VecDeque;

/// `seq[i:i+k]` with Python's clamping slice semantics.
///
/// Python never panics on an out-of-range slice, it truncates. Short reads
/// therefore compare truncated k-mers, and reproducing that matters more than
/// it looks: the comparison is lexicographic, and a truncated k-mer sorts
/// before any string it is a prefix of.
fn kmer(seq: &[u8], i: usize, k: usize) -> &[u8] {
    let start = i.min(seq.len());
    let end = i.saturating_add(k).min(seq.len());
    &seq[start..end]
}

/// Lexicographic minimizers as `(kmer, position)`, in the reference's order.
pub fn minimizers_lex(seq: &[u8], k: usize, w_size: usize) -> Vec<(Vec<u8>, usize)> {
    let mut out = Vec::new();
    if k == 0 || w_size < k {
        return out;
    }
    let w = w_size - k;

    // Initial window covers k-mer start positions 0..=w.
    let mut window: VecDeque<&[u8]> = (0..=w).map(|i| kmer(seq, i, k)).collect();

    let Some(curr) = window.iter().min() else {
        return out;
    };
    let mut curr_min: Vec<u8> = curr.to_vec();
    // rindex: last occurrence, so ties go to the rightmost position.
    let mut minimizer_pos = window
        .iter()
        .rposition(|x| *x == curr_min.as_slice())
        .expect("min came from this window");

    out.push((kmer(seq, minimizer_pos, k).to_vec(), minimizer_pos));

    // Python: range(w + 1, len(seq) - k). Upper bound exclusive; empty when the
    // sequence is shorter than k.
    let upper = seq.len().saturating_sub(k);
    for i in (w + 1)..upper {
        let new_kmer = kmer(seq, i, k);
        let discarded = window.pop_front().expect("window is never empty");
        window.push_back(new_kmer);

        // `i - w` is the leftmost start position now in the window.
        if discarded == curr_min.as_slice() && minimizer_pos < i - w {
            let m = window.iter().min().expect("window is never empty").to_vec();
            let rel = window
                .iter()
                .rposition(|x| *x == m.as_slice())
                .expect("min came from this window");
            curr_min = m;
            minimizer_pos = rel + i - w;
            out.push((kmer(seq, minimizer_pos, k).to_vec(), minimizer_pos));
        } else if new_kmer < curr_min.as_slice() {
            curr_min = new_kmer.to_vec();
            minimizer_pos = i;
            out.push((kmer(seq, i, k).to_vec(), i));
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(seq: &str, k: usize, w: usize) -> Vec<(String, usize)> {
        minimizers_lex(seq.as_bytes(), k, w)
            .into_iter()
            .map(|(m, p)| (String::from_utf8(m).unwrap(), p))
            .collect()
    }

    #[test]
    fn ties_resolve_to_the_last_position() {
        // Three identical k-mers in the opening window; rindex picks the last.
        let got = run("AAAAAAAAAAAA", 3, 5);
        assert_eq!(got[0].1, 2, "expected the rightmost tied position");
    }

    #[test]
    fn first_minimizer_comes_from_the_opening_window() {
        // Window covers start positions 0..=(w-k) = 0..=2 for k=3, w=5.
        // k-mers: TTT, TTA, TAA -> min is TAA at position 2.
        let got = run("TTTAACCCGGG", 3, 5);
        assert_eq!(got[0], ("TAA".to_string(), 2));
    }

    #[test]
    fn final_kmer_is_never_examined() {
        // The reference's loop upper bound is exclusive, so a uniquely smallest
        // k-mer sitting at the very last start position is skipped.
        let seq = "TTTTTTTTTTAAA"; // len 13, k=3 -> last start position is 10
        let got = run(seq, 3, 5);
        assert!(
            !got.iter().any(|(_, p)| *p == 10),
            "position 10 is the final k-mer and must not appear: {got:?}"
        );
    }

    #[test]
    fn short_sequences_yield_the_reference_s_degenerate_result() {
        // Verified directly against get_kmer_minimizers_lex: any sequence too
        // short to fill the opening window yields exactly [("", 11)] at k=9,
        // w=20. Python clamps out-of-range slices to "", empty sorts first, and
        // rindex picks the last of the w+1 tied entries -- so position 11, which
        // does not even exist in the sequence.
        //
        // This is nonsense, but it is the reference's nonsense, and reads this
        // short would have to be handled identically to stay byte-compatible.
        for s in ["", "A", "ACGT", "ACGTACGTA"] {
            assert_eq!(
                run(s, 9, 20),
                vec![(String::new(), 11)],
                "short sequence {s:?} must match the reference"
            );
        }
        // One k-mer longer than the window and it behaves sensibly again.
        assert_eq!(
            run("ACGTACGTACGTACGTACGTACGT", 9, 20),
            vec![("ACGTACGTA".to_string(), 8)]
        );
    }

    #[test]
    fn window_equal_to_k_emits_every_position() {
        // w == k means w_size - k == 0, so the window holds a single k-mer and
        // every new k-mer smaller than the current one is emitted.
        let got = run("GGGTTTAAA", 3, 3);
        assert!(!got.is_empty());
        assert_eq!(got[0].1, 0);
    }

    #[test]
    fn positions_are_strictly_increasing() {
        let seq = "GGGTCAGATGCCTGTAATGAGCAACACAACTGGGCCCATGGTAGGGTCAACAATTCTCTCC";
        let got = run(seq, 9, 20);
        assert!(got.len() > 1);
        for pair in got.windows(2) {
            assert!(
                pair[0].1 < pair[1].1,
                "positions must increase: {:?} then {:?}",
                pair[0],
                pair[1]
            );
        }
    }

    #[test]
    fn reported_kmer_matches_the_sequence_at_that_position() {
        let seq = "GGGTCAGATGCCTGTAATGAGCAACACAACTGGGCCCATGGTAGGGTCAACAATTCTCTCC";
        for (m, p) in run(seq, 9, 20) {
            assert_eq!(&seq[p..p + m.len()], m);
        }
    }
}
