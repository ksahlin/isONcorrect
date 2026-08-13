//! `isONcorrect` --- correct a single cluster.
//!
//! CLI parity with `src/isoncorrect/isONcorrect.py`. The correction algorithm
//! is not ported yet; this validates arguments exactly as the reference does
//! and then reports that clearly.

use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isoncorrect::cli::{unsupported_flag, CorrectArgs, EXIT_UNSUPPORTED};
use isoncorrect::params::Params;
use isoncorrect::validate::{check_window, clamp_xmin};

fn main() -> ExitCode {
    // The reference prints help and exits 0 when invoked with no arguments at
    // all, so this has to be decided before parsing (nothing is required).
    if std::env::args().len() == 1 {
        CorrectArgs::command()
            .print_help()
            .expect("failed to write help");
        println!();
        return ExitCode::SUCCESS;
    }

    let args = CorrectArgs::parse();

    // Reject unported flags before doing anything else, so a pipeline that
    // passes --use_racon stops rather than quietly producing unpolished output.
    if let Some(msg) = unsupported_flag(
        args.randstrobes,
        args.layers,
        args.set_layers_manually,
        args.use_racon,
    ) {
        eprintln!("{msg}");
        return ExitCode::from(EXIT_UNSUPPORTED as u8);
    }

    // Order below mirrors isONcorrect.py::main exactly. It is observable: the
    // xmin notice precedes the window check, and the outfolder is created
    // before the window check, so an invalid --w still leaves the folder behind.
    let (xmin, notice) = clamp_xmin(args.k, args.xmin);
    if let Some(notice) = notice {
        eprintln!("{notice}");
    }

    // The reference raises AttributeError here (it tests args.flnc / args.ccs,
    // which argparse never defines). Printing help is the evident intent; this
    // is an error-path divergence, noted in PORTING.md.
    let Some(fastq) = args.fastq.as_ref() else {
        CorrectArgs::command()
            .print_help()
            .expect("failed to write help");
        println!();
        return ExitCode::SUCCESS;
    };

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

    if let Err(msg) = check_window(args.k, args.w) {
        eprintln!("{msg}");
        return ExitCode::FAILURE;
    }

    let params = Params {
        k: args.k,
        w: args.w,
        xmin,
        xmax: args.xmax,
        t_threshold: args.t_threshold,
        max_seqs: args.max_seqs,
        max_seqs_to_spoa: args.max_seqs_to_spoa,
        exact: args.exact,
        exact_instance_limit: args.exact_instance_limit,
        set_w_dynamically: args.set_w_dynamically,
    };

    eprintln!(
        "isONcorrect (Rust port): argument handling is in place, but the correction\n\
         algorithm is not ported yet, so no output was written.\n\
         \n\
         Would have corrected: {}\n\
         Parameters: {params:?}\n\
         \n\
         Use the Python implementation for real runs. Progress: PORTING.md",
        fastq.display()
    );
    ExitCode::FAILURE
}
