//! `isONcorrect` --- correct a single cluster.
//!
//! CLI parity with `src/isoncorrect/isONcorrect.py`, and the correction itself
//! via [`isoncorrect::driver`]. Argument handling mirrors the reference's
//! `main` exactly, including the order side effects happen in.

use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isoncorrect::cli::{unsupported_flag, CorrectArgs, EXIT_UNSUPPORTED};
use isoncorrect::driver;
use isoncorrect::fastq;
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

    let reads = match std::fs::File::open(fastq)
        .and_then(|f| fastq::read_fastq(std::io::BufReader::new(f)))
    {
        Ok(reads) => reads,
        Err(e) => {
            eprintln!("error: could not read {}: {e}", fastq.display());
            return ExitCode::FAILURE;
        }
    };
    eprintln!("Total cluster of {} reads.", reads.len());

    let (corrected, _stats) = driver::correct_cluster(&reads, &params);

    // The reference always writes into --outfolder; argparse defaults it, so a
    // missing one here means the caller passed an empty value.
    let Some(outfolder) = args.outfolder.as_ref() else {
        eprintln!("error: --outfolder is required");
        return ExitCode::FAILURE;
    };
    let path = outfolder.join("corrected_reads.fastq");
    let file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: could not write {}: {e}", path.display());
            return ExitCode::FAILURE;
        }
    };
    let mut writer = std::io::BufWriter::new(file);
    if let Err(e) = driver::write_fastq(&mut writer, &corrected) {
        eprintln!("error: could not write {}: {e}", path.display());
        return ExitCode::FAILURE;
    }
    if let Err(e) = writer.into_inner().map(|f| f.sync_all()) {
        eprintln!("error: could not flush {}: {e:?}", path.display());
        return ExitCode::FAILURE;
    }

    let (checked, differed) = isoncorrect::parasail::band_check_report();
    if checked > 0 {
        eprintln!("band-check: {checked} guard alignments, {differed} differed from exact");
    }
    let (checked, cigars, outputs) = isoncorrect::parasail::global_check_report();
    if checked > 0 {
        eprintln!(
            "global-check: {checked} alignments, {cigars} different CIGARs, \
             {outputs} different guard outputs"
        );
    }

    ExitCode::SUCCESS
}
