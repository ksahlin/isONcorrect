//! `run_isoncorrect` --- fan `isONcorrect` out over a folder of clusters.
//!
//! CLI parity with `src/isoncorrect/run_isoncorrect.py`. The work distribution
//! is not ported yet; this validates arguments exactly as the reference does
//! and then reports that clearly.
//!
//! When it is ported, note that the reference spawns `python isONcorrect.py`
//! per batch via `subprocess`. In Rust this should be in-process work
//! distribution --- see "Deferred improvements" in PORTING.md.

use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isoncorrect::cli::{unsupported_flag, RunArgs, EXIT_UNSUPPORTED};
use isoncorrect::params::Params;
use isoncorrect::validate::check_split_mod;

fn main() -> ExitCode {
    if std::env::args().len() == 1 {
        RunArgs::command()
            .print_help()
            .expect("failed to write help");
        println!();
        return ExitCode::SUCCESS;
    }

    let args = RunArgs::parse();

    if let Some(msg) = unsupported_flag(
        args.randstrobes,
        args.layers,
        args.set_layers_manually,
        args.use_racon,
    ) {
        eprintln!("{msg}");
        return ExitCode::from(EXIT_UNSUPPORTED as u8);
    }

    // The reference asserts this, dying with AssertionError and exit 1.
    if let Err(msg) = check_split_mod(args.split_mod, args.residual) {
        eprintln!("error: {msg}");
        return ExitCode::FAILURE;
    }

    if let Some(outfolder) = args.outfolder.as_ref() {
        if !outfolder.exists() {
            if let Err(e) = std::fs::create_dir_all(outfolder) {
                eprintln!(
                    "error: could not create outfolder {}: {e}",
                    outfolder.display()
                );
                return ExitCode::FAILURE;
            }
        }
    }

    let Some(fastq_folder) = args.fastq_folder.as_ref() else {
        RunArgs::command()
            .print_help()
            .expect("failed to write help");
        println!();
        return ExitCode::SUCCESS;
    };

    // Note --xmin is NOT clamped here: run_isoncorrect passes it through and the
    // per-cluster isONcorrect process does the clamping. Preserve that.
    let params = Params {
        k: args.k,
        w: args.w,
        xmin: args.xmin,
        xmax: args.xmax,
        t_threshold: args.t_threshold,
        max_seqs: args.max_seqs,
        // run_isoncorrect does not expose --max_seqs_to_spoa, so every cluster
        // gets isONcorrect's own default.
        max_seqs_to_spoa: 200,
        exact: false,
        exact_instance_limit: args.exact_instance_limit,
        set_w_dynamically: args.set_w_dynamically,
    };

    eprintln!(
        "run_isoncorrect (Rust port): argument handling is in place, but the work\n\
         distribution and correction algorithm are not ported yet, so no output\n\
         was written.\n\
         \n\
         Would have corrected clusters in: {}\n\
         Cores: {}  split_mod: {}  residual: {}  split_wrt_batches: {}  keep_old: {}\n\
         Parameters: {params:?}\n\
         \n\
         Use the Python implementation for real runs. Progress: PORTING.md",
        fastq_folder.display(),
        args.nr_cores,
        args.split_mod,
        args.residual,
        args.split_wrt_batches,
        args.keep_old,
    );
    ExitCode::FAILURE
}
