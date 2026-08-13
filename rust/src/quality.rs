//! Per-base error probabilities and their prefix sums, matching `get_qvs`.
//!
//! The reference builds, for every read, a running total of per-base error
//! probabilities:
//!
//! ```python
//! D = {chr(i): min(10**(-(ord(chr(i)) - 33)/10.0), 0.79433) for i in range(128)}
//! quality_values_database[r_id] = [0]
//! tmp_tot_sum = 0
//! for char_ in qual:
//!     quality_values_database[r_id].append(tmp_tot_sum + qv)
//!     tmp_tot_sum += qv
//! ```
//!
//! Two details matter:
//!
//! * the array has a **leading zero**, so it is `len(qual) + 1` long and
//!   `qvs[i]` is the summed error probability of the first `i` bases. Spans are
//!   then `qvs[end] - qvs[start]`, which is why the indices in
//!   `find_most_supported_span` can reach `pos2 + k_size`;
//! * probabilities are **clamped at 0.79433**, the value for Phred 1. Phred 0
//!   would otherwise give 1.0. Anything below `!` (33) clamps to the same
//!   ceiling.

/// Phred 0 maps to 10^0 = 1.0, which the reference caps at this value.
const MAX_ERROR_PROB: f64 = 0.79433;

/// Per-base error probability for a Phred+33 quality character.
pub fn error_probability(c: u8) -> f64 {
    let phred = c as f64 - 33.0;
    let p = 10f64.powf(-phred / 10.0);
    p.min(MAX_ERROR_PROB)
}

/// Prefix sums of per-base error probability, with a leading zero.
///
/// `out[i] - out[j]` is the summed error probability of bases `j..i`.
pub fn prefix_sums(qual: &[u8]) -> Vec<f64> {
    let mut out = Vec::with_capacity(qual.len() + 1);
    out.push(0.0);
    let mut total = 0.0f64;
    for &c in qual {
        // The reference appends `tmp_tot_sum + qv` *then* adds to the running
        // total, so entry i is the sum over the first i bases.
        let qv = error_probability(c);
        out.push(total + qv);
        total += qv;
    }
    out
}

/// Mean error probability over `start..end`, as the reference computes it.
///
/// Panics on an empty span, matching Python's `ZeroDivisionError`; the callers
/// only ever pass spans of at least `k`.
pub fn mean_error(prefix: &[f64], start: usize, end: usize) -> f64 {
    (prefix[end] - prefix[start]) / (end - start) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_values_match_the_reference_table() {
        // 10^(-(q-33)/10), capped at 0.79433.
        assert!((error_probability(b'+') - 0.1).abs() < 1e-12, "Phred 10");
        assert!((error_probability(b'5') - 0.01).abs() < 1e-12, "Phred 20");
        assert!((error_probability(b'?') - 0.001).abs() < 1e-12, "Phred 30");
    }

    #[test]
    fn low_quality_is_clamped() {
        // Phred 0 would be 1.0; the reference caps it.
        assert_eq!(error_probability(b'!'), MAX_ERROR_PROB);
        // Phred 1 is exactly the cap.
        assert!((error_probability(b'"') - MAX_ERROR_PROB).abs() < 1e-5);
        // Anything below '!' clamps too.
        assert_eq!(error_probability(0), MAX_ERROR_PROB);
    }

    #[test]
    fn prefix_has_a_leading_zero_and_is_one_longer() {
        let p = prefix_sums(b"+++++");
        assert_eq!(p.len(), 6);
        assert_eq!(p[0], 0.0);
    }

    #[test]
    fn prefix_entry_i_covers_the_first_i_bases() {
        let p = prefix_sums(b"++++");
        // Each '+' is 0.1, so entry i is i * 0.1.
        for (i, entry) in p.iter().enumerate().take(5) {
            assert!((entry - 0.1 * i as f64).abs() < 1e-12, "entry {i}");
        }
    }

    #[test]
    fn spans_are_differences_of_prefix_entries() {
        let p = prefix_sums(b"++++++++");
        assert!((p[6] - p[2] - 0.4).abs() < 1e-12);
    }

    #[test]
    fn mean_error_over_a_uniform_read_is_the_per_base_value() {
        let p = prefix_sums(b"+++++++++++");
        assert!((mean_error(&p, 0, 10) - 0.1).abs() < 1e-12);
        assert!((mean_error(&p, 3, 7) - 0.1).abs() < 1e-12);
    }

    #[test]
    fn empty_quality_yields_just_the_leading_zero() {
        assert_eq!(prefix_sums(b""), vec![0.0]);
    }
}
