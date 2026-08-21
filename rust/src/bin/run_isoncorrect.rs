//! `run_isoncorrect` --- fan `isONcorrect` out over a folder of clusters.
//!
//! CLI parity with `src/isoncorrect/run_isoncorrect.py`, and the work
//! distribution via [`isoncorrect::fanout`].
//!
//! The reference spawns `python isONcorrect.py` per cluster through a
//! `multiprocessing.Pool`; this runs the same work on an in-process thread
//! pool, which removes one interpreter startup per cluster.

use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isoncorrect::cli::{unsupported_flag, RunArgs, EXIT_UNSUPPORTED};
use isoncorrect::fanout::{self, Fanout};
use isoncorrect::params::Params;
use isoncorrect::validate::{check_split_mod, clamp_xmin};

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

    // The reference passes --xmin through untouched and the *per-cluster*
    // isONcorrect process clamps it to 2k. Running in-process, that clamp has to
    // happen here instead, or the two disagree on every --xmin below 2k.
    let (xmin, notice) = clamp_xmin(args.k, args.xmin);
    if let Some(notice) = notice {
        eprintln!("{notice}");
    }
    let params = Params {
        k: args.k,
        w: args.w,
        xmin,
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

    let Some(outfolder) = args.outfolder.as_ref() else {
        eprintln!("error: --outfolder is required");
        return ExitCode::FAILURE;
    };

    let cfg = Fanout {
        fastq_folder: fastq_folder.clone(),
        outfolder: outfolder.clone(),
        threads: args.nr_cores,
        split_mod: args.split_mod,
        residual: args.residual,
        split_wrt_batches: args.split_wrt_batches,
        keep_old: args.keep_old,
    };

    // Under --split_wrt_batches the work list is built from a temporary
    // directory of per-batch files rather than from the input folder.
    let work_dir = if cfg.split_wrt_batches {
        match tempdir() {
            Ok(dir) => {
                println!("Temporary workdirektory: {}", dir.display());
                Some(dir)
            }
            Err(e) => {
                eprintln!("error: could not create a temporary directory: {e}");
                return ExitCode::FAILURE;
            }
        }
    } else {
        None
    };
    let source = match work_dir.as_ref() {
        Some(dir) => match fanout::split_clusters(fastq_folder, dir, args.max_seqs) {
            Ok(split) => split,
            Err(e) => {
                eprintln!("error: {e}");
                return ExitCode::FAILURE;
            }
        },
        None => fastq_folder.clone(),
    };

    let instances = match fanout::plan(&source, &cfg) {
        Ok(instances) => instances,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    println!("Using {} cores.", args.nr_cores);
    if let Err(e) = fanout::run_all(&instances, &params, args.nr_cores) {
        eprintln!("error: {e}");
        return ExitCode::FAILURE;
    }

    if cfg.split_wrt_batches {
        if let Err(e) = fanout::join_batches(&source, outfolder, args.split_mod, args.residual) {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
        if let Some(dir) = work_dir.as_ref() {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    ExitCode::SUCCESS
}

/// A private temporary directory, matching `tempfile.mkdtemp()`.
fn tempdir() -> std::io::Result<std::path::PathBuf> {
    let base = std::env::temp_dir();
    for attempt in 0..64 {
        let candidate = base.join(format!("run_isoncorrect_{}_{attempt}", std::process::id()));
        match std::fs::create_dir(&candidate) {
            Ok(()) => return Ok(candidate),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        "could not find an unused temporary directory name",
    ))
}
