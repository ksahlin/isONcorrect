//! Argument validation, in the reference's order.
//!
//! Order matters and is observable. In `isONcorrect.py::main` the `--xmin`
//! clamp happens *before* the window check, and the output folder is created
//! *before* the window check too --- so `--w 150` still leaves an empty
//! outfolder behind on the way out. These functions are separate so the
//! sequencing stays visible at the call site in each binary.

/// Clamp `--xmin` up to `2 * k`, as `main` does.
///
/// Returns the effective value and, when it changed, the notice the reference
/// prints to stderr.
pub fn clamp_xmin(k: usize, xmin: usize) -> (usize, Option<String>) {
    if xmin < 2 * k {
        let clamped = 2 * k;
        (clamped, Some(format!("xmin set to {clamped}")))
    } else {
        (xmin, None)
    }
}

/// The window-size check. `Err` carries the reference's exact message.
pub fn check_window(k: usize, w: usize) -> Result<(), String> {
    if 100 < w || w < k {
        Err(
            "Please specify a window of size larger or equal to k, and smaller than 100."
                .to_string(),
        )
    } else {
        Ok(())
    }
}

/// `run_isoncorrect`'s `--split_mod` / `--residual` check.
///
/// The reference uses a bare `assert`, so it dies with `AssertionError` and
/// exit 1. We keep the exit code and say something useful instead.
pub fn check_split_mod(split_mod: usize, residual: usize) -> Result<(), String> {
    if split_mod > 1 && residual >= split_mod {
        Err(format!(
            "--residual must be less than --split_mod (got --split_mod {split_mod} --residual {residual})"
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xmin_is_clamped_to_twice_k() {
        // Observed from the reference: --k 20 --xmin 10 prints "xmin set to 40".
        let (xmin, notice) = clamp_xmin(20, 10);
        assert_eq!(xmin, 40);
        assert_eq!(notice.as_deref(), Some("xmin set to 40"));
    }

    #[test]
    fn xmin_at_or_above_twice_k_is_untouched() {
        let (xmin, notice) = clamp_xmin(9, 18);
        assert_eq!(xmin, 18);
        assert!(notice.is_none());

        let (xmin, notice) = clamp_xmin(9, 80);
        assert_eq!(xmin, 80);
        assert!(notice.is_none());
    }

    #[test]
    fn default_k_and_xmin_are_exactly_at_the_boundary() {
        // Defaults are k=9, xmin=18 == 2k, so the default path never clamps.
        let (xmin, notice) = clamp_xmin(9, 18);
        assert_eq!(xmin, 18);
        assert!(notice.is_none());
    }

    #[test]
    fn window_must_be_between_k_and_100() {
        assert!(check_window(9, 20).is_ok());
        assert!(check_window(9, 9).is_ok());
        assert!(check_window(9, 100).is_ok());
        // w < k
        assert!(check_window(15, 10).is_err());
        // w > 100
        assert!(check_window(9, 150).is_err());
        assert!(check_window(9, 101).is_err());
    }

    #[test]
    fn window_error_matches_reference_text() {
        let err = check_window(9, 150).unwrap_err();
        assert_eq!(
            err,
            "Please specify a window of size larger or equal to k, and smaller than 100."
        );
    }

    #[test]
    fn split_mod_requires_smaller_residual() {
        assert!(check_split_mod(1, 0).is_ok());
        // split_mod == 1 disables the check entirely in the reference, so even a
        // nonsensical residual passes.
        assert!(check_split_mod(1, 7).is_ok());
        assert!(check_split_mod(2, 0).is_ok());
        assert!(check_split_mod(2, 1).is_ok());
        assert!(check_split_mod(2, 2).is_err());
        assert!(check_split_mod(2, 5).is_err());
    }
}
