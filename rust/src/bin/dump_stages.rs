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
use isoncorrect::anchors::{self, Minimizer};
use isoncorrect::fastq::read_fastq;
use isoncorrect::minimizers::minimizers_lex;
use isoncorrect::params::Params;
use isoncorrect::quality;
use isoncorrect::regions;
use isoncorrect::support::SpanFinder;
use isoncorrect::validate::clamp_xmin;
use std::collections::{BTreeSet, HashMap};

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
    #[arg(long, default_value_t = 18)]
    xmin: usize,
    #[arg(long, default_value_t = 80)]
    xmax: usize,
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

    // isONcorrect.py::main clamps --xmin up to 2*k before anything else runs,
    // so the dump has to do it too or the anchor spans diverge.
    let (xmin, notice) = clamp_xmin(args.k, args.xmin);
    if let Some(n) = notice {
        eprintln!("{n}");
    }

    let params = Params {
        k: args.k,
        w: args.w,
        xmin,
        xmax: args.xmax,
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
        let mut by_read: Vec<(usize, Vec<Minimizer>)> = Vec::new();
        for r in chunk {
            let mins = minimizers_lex(r.seq.as_bytes(), params.k, w);
            by_read.push((r.r_id, mins.clone()));
            for (i, (m, p)) in mins.into_iter().enumerate() {
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

        // Anchor index, dumped in the same format as bench/dump_reference.py.
        let db = anchors::build(&by_read, params.k, params.xmin, params.xmax, chunk.len());
        let mut af = BufWriter::new(fs::File::create(dir.join("anchors.tsv")).unwrap());
        let mut kf = BufWriter::new(fs::File::create(dir.join("anchor_keys.tsv")).unwrap());
        writeln!(af, "#m1\tm2\tr_id\tp1\tp2").unwrap();
        writeln!(kf, "#m1\tm2\tn_entries").unwrap();
        let mut n_entries = 0usize;
        for (m1, m2, occ) in db.iter() {
            let (m1s, m2s) = (String::from_utf8_lossy(m1), String::from_utf8_lossy(m2));
            writeln!(kf, "{}\t{}\t{}", m1s, m2s, occ.len()).unwrap();
            for (r_id, p1, p2) in occ {
                writeln!(af, "{}\t{}\t{}\t{}\t{}", m1s, m2s, r_id, p1, p2).unwrap();
                n_entries += 1;
            }
        }
        af.flush().unwrap();
        kf.flush().unwrap();
        // Per-read span collection, mirroring isoncorrect_main's inner loop.
        //
        // previously_corrected_regions is left empty: populating it needs the
        // correction stage. Under --exact the reference clears it per read, so
        // an --exact run is directly comparable to this.
        let seq_by_id: Vec<Vec<u8>> = chunk.iter().map(|r| r.seq.as_bytes().to_vec()).collect();
        let qvs: Vec<Vec<f64>> = chunk
            .iter()
            .map(|r| quality::prefix_sums(r.qual.as_bytes()))
            .collect();
        let mut finder = SpanFinder::new(&seq_by_id, &qvs, params.k);
        let no_regions = BTreeSet::new();
        let no_groups = HashMap::new();

        let mut sf = BufWriter::new(fs::File::create(dir.join("spans.tsv")).unwrap());
        writeln!(sf, "#r_id\tm1\tp1\tstart\tstop\tweight\tinstance").unwrap();
        let mut n_spans = 0usize;
        for (r_id, mins) in &by_read {
            // --exact resets the cache per read.
            finder.reset_cache();
            for ((m1, p1), spans) in anchors::comb_iterator(mins, params.xmin, params.xmax) {
                let kept = regions::filter_spans(&spans, p1, params.k, &no_regions, &no_groups);
                if kept.is_empty() {
                    continue;
                }
                let mut out_iv = Vec::new();
                finder.find(*r_id, m1, p1, &kept, &db, &mut out_iv);
                for c in out_iv {
                    let payload: Vec<String> = c
                        .instance
                        .iter()
                        .flat_map(|&(a, b, cc)| [a.to_string(), b.to_string(), cc.to_string()])
                        .collect();
                    writeln!(
                        sf,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        r_id,
                        String::from_utf8_lossy(m1),
                        p1,
                        c.start,
                        c.stop,
                        c.support,
                        payload.join(",")
                    )
                    .unwrap();
                    n_spans += 1;
                }
            }
        }
        sf.flush().unwrap();
        eprintln!(
            "batch {batch_id:03}: spans={n_spans} edlib_calls={}",
            finder.edlib_calls
        );

        eprintln!(
            "batch {batch_id:03}: anchor_keys={} entries={} generated={} singletons={} high_abundance={}",
            db.len(),
            n_entries,
            db.generated,
            db.singletons,
            db.high_abundance
        );

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
