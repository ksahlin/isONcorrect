//! Command-line parity with the Python implementation.
//!
//! Flag names, defaults, validation order and exit codes are copied from
//! `src/isoncorrect/isONcorrect.py::main` and `run_isoncorrect.py::main`. Where
//! this file and the Python reference disagree, the reference is right --- see
//! `PORTING.md`.
//!
//! Two groups of flags are deliberately absent or rejected; see "Scope" in
//! `PORTING.md` for why. The observable contract, which `bench/equivalence.sh`
//! asserts, is the same for both groups: **non-zero exit, and the flag named in
//! the output**.

use std::path::PathBuf;

use clap::{ArgAction, Parser};

/// The version string the CLI advertises. Matches argparse's
/// `%(prog)s 0.2.0`, so `--version` prints e.g. `isONcorrect 0.2.0`.
pub const VERSION: &str = "0.2.0";

/// Exit code used when a deliberately unported flag is supplied.
pub const EXIT_UNSUPPORTED: i32 = 2;

/// `isONcorrect` --- correct a single cluster.
#[derive(Debug, Parser)]
#[command(
    name = "isONcorrect",
    version = VERSION,
    about = "De novo error correction of long-read transcriptome reads",
    disable_help_subcommand = true
)]
pub struct CorrectArgs {
    /// Path to input fastq file with reads
    #[arg(long)]
    pub fastq: Option<PathBuf>,

    /// Kmer size
    #[arg(long, default_value_t = 9)]
    pub k: usize,

    /// Window size
    #[arg(long, default_value_t = 20)]
    pub w: usize,

    /// Lower interval length
    #[arg(long, default_value_t = 18)]
    pub xmin: usize,

    /// Upper interval length
    #[arg(long, default_value_t = 80)]
    pub xmax: usize,

    /// Minimum fraction keeping substitution
    #[arg(long = "T", default_value_t = 0.1)]
    pub t_threshold: f64,

    /// Get exact solution for WIS for every read (much slower)
    #[arg(long, action = ArgAction::SetTrue)]
    pub exact: bool,

    /// Maximum number of seqs to spoa
    #[arg(long = "max_seqs_to_spoa", default_value_t = 200)]
    pub max_seqs_to_spoa: usize,

    /// Maximum number of seqs to correct at a time (in case of large clusters)
    #[arg(long = "max_seqs", default_value_t = 2000)]
    pub max_seqs: usize,

    /// Activates slower exact mode for instances smaller than this limit
    #[arg(long = "exact_instance_limit", default_value_t = 0)]
    pub exact_instance_limit: usize,

    /// Set w = k + max(2*k, floor(cluster_size/1000))
    #[arg(long = "set_w_dynamically", action = ArgAction::SetTrue)]
    pub set_w_dynamically: bool,

    /// Print various developer stats
    #[arg(long, action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Output folder
    #[arg(long)]
    pub outfolder: Option<PathBuf>,

    // --- Tier 2: recognised so the error can name them, then rejected. ---
    // These work in the Python implementation, so they may appear in existing
    // pipeline scripts. Hidden from --help because they are not supported.
    #[arg(long, action = ArgAction::SetTrue, hide = true)]
    pub randstrobes: bool,

    #[arg(long, hide = true)]
    pub layers: Option<i64>,

    #[arg(long = "set_layers_manually", action = ArgAction::SetTrue, hide = true)]
    pub set_layers_manually: bool,

    #[arg(long = "use_racon", action = ArgAction::SetTrue, hide = true)]
    pub use_racon: bool,
    // --- Tier 1: --disable_numpy and --compression are absent entirely. ---
    // clap rejects them with "unexpected argument '--disable_numpy' found",
    // which already satisfies the contract (non-zero exit, flag named).
}

/// `run_isoncorrect` --- fan `isONcorrect` out over a folder of clusters.
#[derive(Debug, Parser)]
#[command(
    name = "run_isoncorrect",
    version = VERSION,
    about = "De novo clustering of long-read transcriptome reads",
    disable_help_subcommand = true
)]
pub struct RunArgs {
    /// Path to input fastq folder with reads in clusters
    #[arg(long = "fastq_folder")]
    pub fastq_folder: Option<PathBuf>,

    /// Number of cores allocated for clustering
    #[arg(long = "t", default_value_t = 8)]
    pub nr_cores: usize,

    /// Kmer size
    #[arg(long, default_value_t = 9)]
    pub k: usize,

    /// Window size
    #[arg(long, default_value_t = 20)]
    pub w: usize,

    /// Lower interval length
    #[arg(long, default_value_t = 18)]
    pub xmin: usize,

    /// Upper interval length
    #[arg(long, default_value_t = 80)]
    pub xmax: usize,

    /// Minimum fraction keeping substitution
    #[arg(long = "T", default_value_t = 0.1)]
    pub t_threshold: f64,

    /// Do exact correction for clusters under this threshold
    #[arg(long = "exact_instance_limit", default_value_t = 50)]
    pub exact_instance_limit: usize,

    /// Do not recompute previous results if corrected_reads.fastq is complete
    #[arg(long = "keep_old", action = ArgAction::SetTrue)]
    pub keep_old: bool,

    /// Set w = k + max(2*k, floor(cluster_size/1000))
    #[arg(long = "set_w_dynamically", action = ArgAction::SetTrue)]
    pub set_w_dynamically: bool,

    /// Maximum number of seqs to correct at a time (in case of large clusters)
    #[arg(long = "max_seqs", default_value_t = 2000)]
    pub max_seqs: usize,

    /// Split cluster ids into n partitions by residual of cluster_id / n
    #[arg(long = "split_mod", default_value_t = 1)]
    pub split_mod: usize,

    /// Run on cluster ids with this residual of cluster_id / --split_mod
    #[arg(long, default_value_t = 0)]
    pub residual: usize,

    /// Process reads per batch instead of per cluster
    #[arg(long = "split_wrt_batches", action = ArgAction::SetTrue)]
    pub split_wrt_batches: bool,

    /// Outfolder with all corrected reads
    #[arg(long)]
    pub outfolder: Option<PathBuf>,

    // --- Tier 2, as above ---
    #[arg(long, action = ArgAction::SetTrue, hide = true)]
    pub randstrobes: bool,

    #[arg(long, hide = true)]
    pub layers: Option<i64>,

    #[arg(long = "set_layers_manually", action = ArgAction::SetTrue, hide = true)]
    pub set_layers_manually: bool,

    #[arg(long = "use_racon", action = ArgAction::SetTrue, hide = true)]
    pub use_racon: bool,
}

/// Build the rejection message for an unported flag.
///
/// The message must name the flag: `bench/equivalence.sh` greps for it, and a
/// user whose pipeline silently lost `--use_racon` polishing is worse off than
/// one whose run stopped.
fn unsupported_message(flag: &str, reason: &str) -> String {
    format!(
        "error: {flag} is not supported by the Rust port of isONcorrect.\n\
         \n\
         {reason}\n\
         \n\
         See the \"Scope\" section of PORTING.md. Use the Python implementation\n\
         if you need this flag."
    )
}

/// First unported flag present, if any, as a ready-to-print message.
pub fn unsupported_flag(
    randstrobes: bool,
    layers: Option<i64>,
    set_layers_manually: bool,
    use_racon: bool,
) -> Option<String> {
    if randstrobes {
        return Some(unsupported_message(
            "--randstrobes",
            "It is experimental, did not improve accuracy, and is not reproducible:\n\
             it depends on Python's randomised string hash, so its output changes\n\
             between runs. There is no fixed behaviour to reproduce.",
        ));
    }
    if layers.is_some() {
        return Some(unsupported_message(
            "--layers",
            "It only takes effect together with --randstrobes, which is not supported.",
        ));
    }
    if set_layers_manually {
        return Some(unsupported_message(
            "--set_layers_manually",
            "It only takes effect together with --randstrobes, which is not supported.",
        ));
    }
    if use_racon {
        return Some(unsupported_message(
            "--use_racon",
            "It shells out to the `racon` binary once per correction interval,\n\
             reintroducing the per-interval process spawn this port exists to remove.",
        ));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn correct_defaults_match_python() {
        let a = CorrectArgs::parse_from(["isONcorrect"]);
        assert_eq!(a.k, 9);
        assert_eq!(a.w, 20);
        assert_eq!(a.xmin, 18);
        assert_eq!(a.xmax, 80);
        assert_eq!(a.t_threshold, 0.1);
        assert_eq!(a.max_seqs_to_spoa, 200);
        assert_eq!(a.max_seqs, 2000);
        // Note: 0 here, but 50 in run_isoncorrect. The asymmetry is in the
        // reference and is load-bearing --- see PORTING.md on --exact.
        assert_eq!(a.exact_instance_limit, 0);
        assert!(!a.exact);
        assert!(!a.set_w_dynamically);
        assert!(!a.verbose);
        assert!(a.fastq.is_none());
        assert!(a.outfolder.is_none());
    }

    #[test]
    fn run_defaults_match_python() {
        let a = RunArgs::parse_from(["run_isoncorrect"]);
        assert_eq!(a.nr_cores, 8);
        assert_eq!(a.k, 9);
        assert_eq!(a.w, 20);
        assert_eq!(a.xmin, 18);
        assert_eq!(a.xmax, 80);
        assert_eq!(a.t_threshold, 0.1);
        assert_eq!(a.exact_instance_limit, 50);
        assert_eq!(a.max_seqs, 2000);
        assert_eq!(a.split_mod, 1);
        assert_eq!(a.residual, 0);
        assert!(!a.keep_old);
        assert!(!a.split_wrt_batches);
    }

    #[test]
    fn paper_settings_parse() {
        let a =
            CorrectArgs::parse_from(["isONcorrect", "--k", "9", "--w", "10", "--max_seqs", "1000"]);
        assert_eq!((a.k, a.w, a.max_seqs), (9, 10, 1000));
    }

    #[test]
    fn t_uses_capital_flag() {
        let a = CorrectArgs::parse_from(["isONcorrect", "--T", "0.05"]);
        assert_eq!(a.t_threshold, 0.05);
        // `--t` is run_isoncorrect's core count and must not collide with --T.
        let r = RunArgs::parse_from(["run_isoncorrect", "--t", "4", "--T", "0.2"]);
        assert_eq!(r.nr_cores, 4);
        assert_eq!(r.t_threshold, 0.2);
    }

    #[test]
    fn tier1_flags_are_absent_and_named_in_the_error() {
        for flag in ["--disable_numpy", "--compression"] {
            let err = CorrectArgs::try_parse_from(["isONcorrect", flag])
                .expect_err("tier 1 flag must be rejected");
            assert!(
                err.to_string().contains(flag),
                "error must name {flag}, got: {err}"
            );
        }
    }

    #[test]
    fn tier2_flags_parse_then_report_unsupported() {
        let a = CorrectArgs::parse_from(["isONcorrect", "--randstrobes"]);
        let msg = unsupported_flag(a.randstrobes, a.layers, a.set_layers_manually, a.use_racon)
            .expect("must be reported unsupported");
        assert!(msg.contains("--randstrobes"));

        let a = CorrectArgs::parse_from(["isONcorrect", "--use_racon"]);
        let msg = unsupported_flag(a.randstrobes, a.layers, a.set_layers_manually, a.use_racon)
            .expect("must be reported unsupported");
        assert!(msg.contains("--use_racon"));

        let a = CorrectArgs::parse_from(["isONcorrect", "--layers", "3"]);
        assert!(
            unsupported_flag(a.randstrobes, a.layers, a.set_layers_manually, a.use_racon).is_some()
        );
    }

    #[test]
    fn supported_flags_are_not_reported_unsupported() {
        let a = CorrectArgs::parse_from(["isONcorrect", "--exact", "--set_w_dynamically"]);
        assert!(
            unsupported_flag(a.randstrobes, a.layers, a.set_layers_manually, a.use_racon).is_none()
        );
    }

    #[test]
    fn version_string_matches_python() {
        assert_eq!(VERSION, "0.2.0");
        // argparse prints "<prog> <version>"; clap does the same.
        assert_eq!(CorrectArgs::command().get_name(), "isONcorrect");
        assert_eq!(RunArgs::command().get_name(), "run_isoncorrect");
    }

    #[test]
    fn clap_definitions_are_well_formed() {
        CorrectArgs::command().debug_assert();
        RunArgs::command().debug_assert();
    }

    fn long_names(mut cmd: clap::Command) -> std::collections::BTreeSet<String> {
        // clap materialises the built-in --help/--version only on build(),
        // so without this they are missing from the surface we assert on.
        cmd.build();
        cmd.get_arguments()
            .filter_map(|a| a.get_long().map(str::to_string))
            .collect()
    }

    /// The exact flag set, locked against the Python reference.
    ///
    /// This exists because clap derives `--field-name` from `field_name` by
    /// default, which silently renamed `--max_seqs` to `--max-seqs` and would
    /// have broken every multi-word flag. Any change to the surface has to be
    /// made here deliberately.
    #[test]
    fn correct_flag_set_matches_python() {
        let expected: std::collections::BTreeSet<String> = [
            // argparse built-ins
            "help",
            "version",
            // ported
            "fastq",
            "k",
            "w",
            "xmin",
            "xmax",
            "T",
            "exact",
            "max_seqs_to_spoa",
            "max_seqs",
            "exact_instance_limit",
            "set_w_dynamically",
            "verbose",
            "outfolder",
            // tier 2: parsed so the rejection can name them
            "randstrobes",
            "layers",
            "set_layers_manually",
            "use_racon",
            // tier 1 (--disable_numpy, --compression) are deliberately absent
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();

        assert_eq!(long_names(CorrectArgs::command()), expected);
    }

    #[test]
    fn run_flag_set_matches_python() {
        let expected: std::collections::BTreeSet<String> = [
            "help",
            "version",
            "fastq_folder",
            "t",
            "k",
            "w",
            "xmin",
            "xmax",
            "T",
            "exact_instance_limit",
            "keep_old",
            "set_w_dynamically",
            "max_seqs",
            "split_mod",
            "residual",
            "split_wrt_batches",
            "outfolder",
            "randstrobes",
            "layers",
            "set_layers_manually",
            "use_racon",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();

        assert_eq!(long_names(RunArgs::command()), expected);
    }

    /// No flag may contain a hyphen: the reference uses underscores throughout,
    /// and a hyphenated alias would be a silently different CLI.
    #[test]
    fn no_flag_is_hyphenated() {
        for (what, names) in [
            ("isONcorrect", long_names(CorrectArgs::command())),
            ("run_isoncorrect", long_names(RunArgs::command())),
        ] {
            for name in names {
                assert!(
                    !name.contains('-'),
                    "{what}: --{name} is hyphenated; the reference uses underscores"
                );
            }
        }
    }

    /// Every underscore-bearing flag actually parses under that exact spelling.
    #[test]
    fn underscore_flags_parse() {
        let a = CorrectArgs::parse_from([
            "isONcorrect",
            "--max_seqs",
            "40",
            "--max_seqs_to_spoa",
            "20",
            "--exact_instance_limit",
            "200",
            "--set_w_dynamically",
        ]);
        assert_eq!(a.max_seqs, 40);
        assert_eq!(a.max_seqs_to_spoa, 20);
        assert_eq!(a.exact_instance_limit, 200);
        assert!(a.set_w_dynamically);

        let r = RunArgs::parse_from([
            "run_isoncorrect",
            "--fastq_folder",
            "clusters",
            "--split_wrt_batches",
            "--split_mod",
            "2",
            "--residual",
            "1",
            "--keep_old",
        ]);
        assert_eq!(r.fastq_folder.unwrap().to_str().unwrap(), "clusters");
        assert!(r.split_wrt_batches);
        assert_eq!((r.split_mod, r.residual), (2, 1));
        assert!(r.keep_old);
    }
}
