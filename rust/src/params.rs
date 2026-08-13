//! Resolved algorithm parameters.
//!
//! `run_isoncorrect` forwards exactly this set down to each per-cluster run, so
//! keeping it in one place stops the two binaries drifting apart. The Python
//! reference threads the same values through a dict in
//! `run_isoncorrect.py::run_isoncorrect_parallel`.

/// Parameters that affect corrected output.
///
/// Anything here changes results and therefore needs a case in
/// `bench/equivalence.sh`. Anything that does not affect output (`--verbose`,
/// `--t`) belongs elsewhere.
#[derive(Debug, Clone, PartialEq)]
pub struct Params {
    pub k: usize,
    pub w: usize,
    pub xmin: usize,
    pub xmax: usize,
    /// `--T`: minimum fraction for keeping a substitution.
    pub t_threshold: f64,
    pub max_seqs: usize,
    pub max_seqs_to_spoa: usize,
    pub exact: bool,
    pub exact_instance_limit: usize,
    pub set_w_dynamically: bool,
}

impl Params {
    /// Window used for a batch of `n_reads`.
    ///
    /// Mirrors `isoncorrect_main`: with `--set_w_dynamically`, batches of 100
    /// or more reads get `min(w, k + (n/100 + 4))`, smaller ones
    /// `k + 1 + n/30`. Note the help text for the flag describes a different
    /// formula (`k + max(2k, n/1000)`) than the code implements --- the code is
    /// the reference.
    pub fn window_for_batch(&self, n_reads: usize) -> usize {
        if !self.set_w_dynamically {
            return self.w;
        }
        if n_reads >= 100 {
            self.w.min(self.k + (n_reads / 100 + 4))
        } else {
            self.k + 1 + n_reads / 30
        }
    }

    /// Whether this batch runs the exact path.
    ///
    /// `--exact` forces it, and `isoncorrect_main` additionally sets it for any
    /// batch of at most `--exact_instance_limit` reads. Since
    /// `run_isoncorrect` defaults that limit to 50, every small cluster takes
    /// the exact path by default --- it is not an opt-in.
    pub fn exact_for_batch(&self, n_reads: usize) -> bool {
        self.exact || n_reads <= self.exact_instance_limit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base() -> Params {
        Params {
            k: 9,
            w: 20,
            xmin: 18,
            xmax: 80,
            t_threshold: 0.1,
            max_seqs: 2000,
            max_seqs_to_spoa: 200,
            exact: false,
            exact_instance_limit: 0,
            set_w_dynamically: false,
        }
    }

    #[test]
    fn static_window_is_just_w() {
        let p = base();
        assert_eq!(p.window_for_batch(5), 20);
        assert_eq!(p.window_for_batch(5000), 20);
    }

    #[test]
    fn dynamic_window_matches_reference_branches() {
        let p = Params {
            set_w_dynamically: true,
            ..base()
        };
        // n >= 100: min(w, k + (n/100 + 4))
        assert_eq!(p.window_for_batch(100), 9 + (1 + 4));
        assert_eq!(p.window_for_batch(250), 9 + (2 + 4));
        // capped by w once k + n/100 + 4 exceeds it
        assert_eq!(p.window_for_batch(100_000), 20);
        // n < 100: k + 1 + n/30
        assert_eq!(p.window_for_batch(0), 10);
        assert_eq!(p.window_for_batch(59), 9 + 1 + 1);
        assert_eq!(p.window_for_batch(90), 9 + 1 + 3);
    }

    #[test]
    fn exact_limit_makes_small_batches_exact() {
        // run_isoncorrect's default: every cluster of <=50 reads is exact.
        let p = Params {
            exact_instance_limit: 50,
            ..base()
        };
        assert!(p.exact_for_batch(50));
        assert!(p.exact_for_batch(1));
        assert!(!p.exact_for_batch(51));

        // isONcorrect's own default of 0 still catches empty batches, matching
        // the reference's `<=` comparison.
        let p = base();
        assert!(!p.exact_for_batch(1));
        assert!(p.exact_for_batch(0));

        // --exact forces it regardless of size.
        let p = Params {
            exact: true,
            ..base()
        };
        assert!(p.exact_for_batch(10_000));
    }
}
