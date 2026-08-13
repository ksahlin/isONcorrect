//! Development tool: dump the port's intermediate structures.
//!
//! Emits exactly the format `bench/dump_reference.py` writes, so each porting
//! stage can be diffed against the reference in isolation:
//!
//! ```text
//! bench/dump_reference.py --fastq c.fastq --outdir /tmp/py
//! isoncorrect-dump --fastq c.fastq --outdir /tmp/rs
//! diff -r /tmp/py/batch_000 /tmp/rs/batch_000
//! ```
//!
//! Not part of the shipped CLI contract; it exists to make failures locatable.

use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use isoncorrect::fastq::read_fastq;
use isoncorrect::minimizers::minimizers_lex;
use isoncorrect::params::Params;

#[derive(Parser)]
#[command(
    name = "isoncorrect-dump",
    about = "Dump the Rust port's intermediate stages"
)]
struct Args {
    #[arg(long)]
    fastq: PathBuf,
    #[arg(long)]
    outdir: PathBuf,
    #[arg(long, default_value_t = 9)]
    k: usize,
    #[arg(long, default_value_t = 20)]
    w: usize,
    #[arg(long = "max_seqs", default_value_t = 2000)]
    max_seqs: usize,
    #[arg(long = "set_w_dynamically", default_value_t = false)]
    set_w_dynamically: bool,
}

fn main() -> ExitCode {
    let args = Args::parse();

    let file = match fs::File::open(&args.fastq) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: {}: {e}", args.fastq.display());
            return ExitCode::FAILURE;
        }
    };
    let reads = match read_fastq(std::io::BufReader::new(file)) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error reading fastq: {e}");
            return ExitCode::FAILURE;
        }
    };

    let params = Params {
        k: args.k,
        w: args.w,
        xmin: 18,
        xmax: 80,
        t_threshold: 0.1,
        max_seqs: args.max_seqs,
        max_seqs_to_spoa: 200,
        exact: false,
        exact_instance_limit: 0,
        set_w_dynamically: args.set_w_dynamically,
    };

    if let Err(e) = fs::create_dir_all(&args.outdir) {
        eprintln!("error: {e}");
        return ExitCode::FAILURE;
    }

    // The reference batches with `batch(all_reads, max_seqs)`; one chunk per
    // dump directory, numbered the same way.
    for (batch_id, chunk) in reads.chunks(args.max_seqs).enumerate() {
        let dir = args.outdir.join(format!("batch_{batch_id:03}"));
        if let Err(e) = fs::create_dir_all(&dir) {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }

        let mut rf = BufWriter::new(fs::File::create(dir.join("reads.tsv")).unwrap());
        writeln!(rf, "#r_id\tacc\tseq\tqual").unwrap();
        for r in chunk {
            writeln!(rf, "{}\t{}\t{}\t{}", r.r_id, r.acc, r.seq, r.qual).unwrap();
        }
        rf.flush().unwrap();

        let w = params.window_for_batch(chunk.len());
        let mut mf = BufWriter::new(fs::File::create(dir.join("minimizers.tsv")).unwrap());
        writeln!(mf, "#r_id\tidx\tkmer\tpos").unwrap();
        let mut total = 0usize;
        for r in chunk {
            for (i, (m, p)) in minimizers_lex(r.seq.as_bytes(), params.k, w)
                .into_iter()
                .enumerate()
            {
                writeln!(
                    mf,
                    "{}\t{}\t{}\t{}",
                    r.r_id,
                    i,
                    String::from_utf8_lossy(&m),
                    p
                )
                .unwrap();
                total += 1;
            }
        }
        mf.flush().unwrap();

        eprintln!(
            "batch {batch_id:03}: reads={} w={} k={} minimizers={}",
            chunk.len(),
            w,
            params.k,
            total
        );
    }

    ExitCode::SUCCESS
}
