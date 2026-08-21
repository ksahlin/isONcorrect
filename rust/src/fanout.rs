//! Fanning out over a folder of clusters, matching `run_isoncorrect.py`.
//!
//! The reference builds a work list from the input directory and hands it to a
//! `multiprocessing.Pool`, where each worker **shells out to
//! `python isONcorrect.py`**. Here the same work runs in-process on a thread
//! pool, which removes one interpreter startup per cluster.
//!
//! Clusters are independent — each one reads a single fastq and writes a single
//! output folder, with no shared state — so scheduling order cannot affect
//! output. That is what makes the parallelism safe, and it is worth stating
//! because it is the only reason this can be reordered at all.
//!
//! # The observable contract
//!
//! ```text
//! <outfolder>/<batch_id>/corrected_reads.fastq
//! <outfolder>/<batch_id>/{stdout,stderr}.txt
//! ```
//!
//! `bench/equivalence.sh` compares every `corrected_reads.fastq` by relative
//! path, so **the folder layout is part of the contract** while the two log
//! files are diagnostics.
//!
//! # Quirks that have to be reproduced
//!
//! * `batch_id` is the filename before the **first** `.`, and `cl_id` is the
//!   part of that before the first `_`. Only `.fastq` files are considered;
//! * `int(cl_id) % split_mod != residual` skips a cluster, so cluster names must
//!   parse as integers;
//! * under `--split_wrt_batches`, joining concatenates batches in
//!   `sorted(glob(...))` order, which is **lexicographic on the path**: a
//!   cluster with ten or more batches concatenates 0, 1, 10, 11, 2, … That
//!   reaches the output file, so it is reproduced;
//! * `split_cluster_in_batches` **latches**. Once one cluster is small enough to
//!   symlink rather than split, every later cluster is symlinked without being
//!   measured. It only makes sense because isONclust numbers clusters
//!   largest-first, so sizes descend — see *Known bugs* in PORTING.md;
//! * `--keep_old` compares **line counts** (`wc -l`) of the output and input
//!   files, not read counts.

use std::collections::BTreeSet;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::driver;
use crate::fastq;
use crate::params::Params;

/// One unit of work: a fastq to correct and where its output folder goes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Instance {
    /// Filename before the first `.` — `12` normally, `12_0` when split.
    pub batch_id: String,
    pub fastq: PathBuf,
    pub outfolder: PathBuf,
}

/// Everything `run_isoncorrect` needs beyond [`Params`].
#[derive(Debug, Clone)]
pub struct Fanout {
    pub fastq_folder: PathBuf,
    pub outfolder: PathBuf,
    pub threads: usize,
    pub split_mod: usize,
    pub residual: usize,
    pub split_wrt_batches: bool,
    pub keep_old: bool,
}

/// `batch_id` and `cl_id` for a fastq filename, or `None` if it is not a fastq.
///
/// `read_fastq_file.split(".")[0]`, then `batch_id.split("_")[0]`.
fn ids(filename: &str) -> Option<(String, String)> {
    if !filename.ends_with(".fastq") {
        return None;
    }
    let batch_id = filename.split('.').next().unwrap_or_default().to_string();
    let cl_id = batch_id.split('_').next().unwrap_or_default().to_string();
    Some((batch_id, cl_id))
}

/// Whether `cl_id` survives `--split_mod` / `--residual`.
///
/// The reference does `int(cl_id) % split_mod != residual`, so a
/// non-numeric cluster name raises. Here it is reported instead of panicking.
fn selected(cl_id: &str, split_mod: usize, residual: usize) -> Result<bool, String> {
    let n: usize = cl_id
        .parse()
        .map_err(|_| format!("cluster id {cl_id:?} is not an integer"))?;
    Ok(n % split_mod.max(1) == residual)
}

/// Line count, matching the reference's `wc -l` shell-out.
fn line_count(path: &Path) -> io::Result<usize> {
    let data = fs::read(path)?;
    Ok(data.iter().filter(|&&b| b == b'\n').count())
}

/// Build the work list, matching `run_isoncorrect_parallel`'s first loop.
///
/// Sorted by `batch_id` **as a string**, which is the reference's scheduling
/// order (`instances.sort(key=lambda x: x[3])`). Output does not depend on it,
/// but matching it keeps progress reporting comparable.
pub fn plan(dir: &Path, cfg: &Fanout) -> Result<Vec<Instance>, String> {
    let entries = fs::read_dir(dir).map_err(|e| format!("cannot read {}: {e}", dir.display()))?;
    let mut out = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|e| format!("cannot read {}: {e}", dir.display()))?;
        let name = entry.file_name().to_string_lossy().to_string();
        let Some((batch_id, cl_id)) = ids(&name) else {
            continue;
        };
        if !selected(&cl_id, cfg.split_mod, cfg.residual)? {
            println!(
                "skipping {batch_id} because args.split_mod:{} and args.residual:{} set.",
                cfg.split_mod, cfg.residual
            );
            continue;
        }
        let fastq = dir.join(&name);
        let outfolder = cfg.outfolder.join(&batch_id);

        if cfg.keep_old {
            let candidate = outfolder.join("corrected_reads.fastq");
            if candidate.is_file() {
                if let (Ok(a), Ok(b)) = (line_count(&candidate), line_count(&fastq)) {
                    if a == b {
                        println!("already computed cluster and complete file {batch_id}");
                        continue;
                    }
                }
            }
        }

        out.push(Instance {
            batch_id,
            fastq,
            outfolder,
        });
    }
    out.sort_by(|a, b| a.batch_id.cmp(&b.batch_id));
    Ok(out)
}

/// Split one fastq into `chunk_lines`-line pieces named `<cl_id>_<i>.fastq`,
/// matching `splitfile`.
///
/// An input whose length is an exact multiple of the chunk size would leave a
/// final empty file; the reference removes it, and so does this.
fn split_file(src: &Path, out_dir: &Path, cl_id: &str, chunk_lines: usize) -> io::Result<usize> {
    let data = fs::read(src)?;
    let mut written = 0usize;
    let mut start = 0usize;
    let mut i = 0usize;
    while start < data.len() {
        let mut end = start;
        for _ in 0..chunk_lines {
            match data[end..].iter().position(|&b| b == b'\n') {
                Some(p) => end += p + 1,
                None => {
                    end = data.len();
                    break;
                }
            }
            if end >= data.len() {
                break;
            }
        }
        let path = out_dir.join(format!("{cl_id}_{i}.fastq"));
        fs::write(&path, &data[start..end])?;
        written += 1;
        start = end;
        i += 1;
    }
    Ok(written)
}

/// Prepare the `--split_wrt_batches` input directory, matching
/// `split_cluster_in_batches`.
///
/// Large clusters are split into `4 * max_seqs`-line files; small ones are
/// symlinked as `<cl_id>_0.fastq` so the rest of the pipeline sees a uniform
/// naming scheme.
///
/// **The size test latches**, exactly as the reference's does: after the first
/// cluster that fits, later clusters are symlinked without being measured.
pub fn split_clusters(indir: &Path, work_dir: &Path, max_seqs: usize) -> Result<PathBuf, String> {
    let split_dir = work_dir.join("split_in_batches");
    fs::create_dir_all(&split_dir).map_err(|e| format!("cannot create split dir: {e}"))?;

    // `sorted(os.listdir(indir), key=lambda x: int(x.split('.')[0]))` --- the
    // key is applied to *every* entry, so a non-numeric name raises there. Here
    // such entries sort last instead of aborting, and are skipped below.
    let mut names: Vec<String> = fs::read_dir(indir)
        .map_err(|e| format!("cannot read {}: {e}", indir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    names.sort_by_key(|n| {
        n.split('.')
            .next()
            .and_then(|s| s.parse::<i64>().ok())
            .unwrap_or(i64::MAX)
    });

    let chunk_lines = 4 * max_seqs;
    let mut latched_small = false;
    for name in names {
        if !name.ends_with(".fastq") {
            continue;
        }
        let src = indir.join(&name);
        let cl_id = name.split('.').next().unwrap_or_default().to_string();

        if !latched_small {
            let lines = line_count(&src).map_err(|e| format!("cannot read {name}: {e}"))?;
            latched_small = lines <= chunk_lines;
        }

        if latched_small {
            // A symlink is enough --- nothing writes to the input.
            let link = split_dir.join(format!("{cl_id}_0.fastq"));
            let _ = fs::remove_file(&link);
            let target = src
                .canonicalize()
                .map_err(|e| format!("cannot resolve {name}: {e}"))?;
            symlink(&target, &link).map_err(|e| format!("cannot link {name}: {e}"))?;
        } else {
            split_file(&src, &split_dir, &cl_id, chunk_lines)
                .map_err(|e| format!("cannot split {name}: {e}"))?;
        }
    }
    Ok(split_dir)
}

#[cfg(unix)]
fn symlink(target: &Path, link: &Path) -> io::Result<()> {
    std::os::unix::fs::symlink(target, link)
}

#[cfg(not(unix))]
fn symlink(target: &Path, link: &Path) -> io::Result<()> {
    // No symlinks to rely on; a copy is observationally identical here, since
    // nothing writes through the link.
    fs::copy(target, link).map(|_| ())
}

/// Concatenate per-batch outputs back into one folder per cluster, matching
/// `join_back_corrected_batches_into_cluster`.
///
/// Batch folders are concatenated in **lexicographic path order** and then
/// removed, so a cluster with ten or more batches joins 0, 1, 10, 11, 2, …
/// That is the reference's `sorted(glob.glob(...))` and it reaches the bytes of
/// the output file.
pub fn join_batches(
    split_dir: &Path,
    outfolder: &Path,
    split_mod: usize,
    residual: usize,
) -> Result<(), String> {
    // Cluster ids come from the *split* directory's filenames, keeping only
    // names of the form `<cl>_<batch>.fastq`.
    let mut cl_ids: BTreeSet<String> = BTreeSet::new();
    let entries =
        fs::read_dir(split_dir).map_err(|e| format!("cannot read {}: {e}", split_dir.display()))?;
    for entry in entries.filter_map(|e| e.ok()) {
        let name = entry.file_name().to_string_lossy().to_string();
        let parts: Vec<&str> = name.split('_').collect();
        if parts.len() != 2 {
            continue;
        }
        let cl_id = parts[0].to_string();
        if selected(&cl_id, split_mod, residual)? {
            cl_ids.insert(cl_id);
        }
    }

    for cl_id in cl_ids {
        let target_dir = outfolder.join(&cl_id);
        fs::create_dir_all(&target_dir).map_err(|e| format!("cannot create outfolder: {e}"))?;
        let outfile = target_dir.join("corrected_reads.fastq");
        fs::write(target_dir.join("cat.stderr"), b"")
            .map_err(|e| format!("cannot write cat.stderr: {e}"))?;

        // `glob(outdir/<cl_id>_*)`, sorted lexicographically.
        let mut batches: Vec<PathBuf> = fs::read_dir(outfolder)
            .map_err(|e| format!("cannot read {}: {e}", outfolder.display()))?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().starts_with(&format!("{cl_id}_")))
                    .unwrap_or(false)
            })
            .collect();
        batches.sort();

        let mut joined: Vec<u8> = Vec::new();
        for batch in &batches {
            let file = batch.join("corrected_reads.fastq");
            if file == outfile {
                continue;
            }
            match fs::read(&file) {
                Ok(bytes) => joined.extend_from_slice(&bytes),
                Err(e) if e.kind() == io::ErrorKind::NotFound => continue,
                Err(e) => return Err(format!("cannot read {}: {e}", file.display())),
            }
        }
        fs::write(&outfile, &joined)
            .map_err(|e| format!("cannot write {}: {e}", outfile.display()))?;
        for batch in &batches {
            let _ = fs::remove_dir_all(batch);
        }
    }
    Ok(())
}

/// Correct one instance and write its output folder.
fn run_one(inst: &Instance, p: &Params) -> Result<usize, String> {
    fs::create_dir_all(&inst.outfolder)
        .map_err(|e| format!("cannot create {}: {e}", inst.outfolder.display()))?;

    let file = fs::File::open(&inst.fastq)
        .map_err(|e| format!("cannot open {}: {e}", inst.fastq.display()))?;
    let reads = fastq::read_fastq(io::BufReader::new(file))
        .map_err(|e| format!("cannot parse {}: {e}", inst.fastq.display()))?;

    let (corrected, stats) = driver::correct_cluster(&reads, p);

    let out = fs::File::create(inst.outfolder.join("corrected_reads.fastq"))
        .map_err(|e| format!("cannot write output for {}: {e}", inst.batch_id))?;
    let mut writer = io::BufWriter::new(out);
    driver::write_fastq(&mut writer, &corrected)
        .map_err(|e| format!("cannot write output for {}: {e}", inst.batch_id))?;
    writer
        .flush()
        .map_err(|e| format!("cannot flush output for {}: {e}", inst.batch_id))?;

    // The reference redirects the child's streams to these two files. They are
    // diagnostics, not part of the compared output, but pipelines look for them.
    let _ = fs::write(
        inst.outfolder.join("stderr.txt"),
        format!("Total cluster of {} reads.\n", reads.len()),
    );
    let _ = fs::write(
        inst.outfolder.join("stdout.txt"),
        format!(
            "corrected {} reads, {} left uncorrected, {} intervals\n",
            corrected.len(),
            stats.uncorrected,
            stats.intervals
        ),
    );
    Ok(reads.len())
}

/// Run every instance across `threads` worker threads.
///
/// Returns the first error encountered, after the running work drains.
pub fn run_all(instances: &[Instance], p: &Params, threads: usize) -> Result<(), String> {
    if instances.is_empty() {
        return Ok(());
    }
    let next = AtomicUsize::new(0);
    let errors: Mutex<Vec<String>> = Mutex::new(Vec::new());

    std::thread::scope(|scope| {
        for _ in 0..threads.max(1).min(instances.len()) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= instances.len() {
                    break;
                }
                let inst = &instances[i];
                match run_one(inst, p) {
                    Ok(n) => println!("Done with batch_id:{} ({n} reads).", inst.batch_id),
                    Err(e) => errors.lock().expect("poisoned").push(e),
                }
            });
        }
    });

    let errors = errors.into_inner().expect("poisoned");
    match errors.into_iter().next() {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("fanout_test_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_cluster(dir: &Path, name: &str, reads: usize) {
        let mut s = String::new();
        for i in 0..reads {
            s.push_str(&format!("@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"));
        }
        fs::write(dir.join(name), s).unwrap();
    }

    fn cfg(indir: &Path, outdir: &Path) -> Fanout {
        Fanout {
            fastq_folder: indir.to_path_buf(),
            outfolder: outdir.to_path_buf(),
            threads: 1,
            split_mod: 1,
            residual: 0,
            split_wrt_batches: false,
            keep_old: false,
        }
    }

    #[test]
    fn only_fastq_files_become_instances() {
        assert_eq!(ids("12.fastq"), Some(("12".to_string(), "12".to_string())));
        assert_eq!(
            ids("12_3.fastq"),
            Some(("12_3".to_string(), "12".to_string()))
        );
        assert_eq!(ids("clusters.tsv"), None);
        assert_eq!(ids("0.fastq.gz"), None);
    }

    #[test]
    fn the_plan_is_sorted_by_batch_id_as_a_string() {
        let indir = tmp("plan");
        let outdir = tmp("plan_out");
        for n in ["0", "1", "2", "10"] {
            write_cluster(&indir, &format!("{n}.fastq"), 1);
        }
        fs::write(indir.join("clusters.tsv"), "ignored").unwrap();
        let plan = plan(&indir, &cfg(&indir, &outdir)).unwrap();
        // String order, so 10 sorts between 1 and 2 --- the reference's own
        // `sort(key=lambda x: x[3])` on batch ids.
        assert_eq!(
            plan.iter().map(|i| i.batch_id.as_str()).collect::<Vec<_>>(),
            vec!["0", "1", "10", "2"]
        );
    }

    #[test]
    fn split_mod_and_residual_select_clusters() {
        let indir = tmp("mod");
        let outdir = tmp("mod_out");
        for n in 0..4 {
            write_cluster(&indir, &format!("{n}.fastq"), 1);
        }
        let mut c = cfg(&indir, &outdir);
        c.split_mod = 2;
        c.residual = 1;
        let plan = plan(&indir, &c).unwrap();
        assert_eq!(
            plan.iter().map(|i| i.batch_id.as_str()).collect::<Vec<_>>(),
            vec!["1", "3"]
        );
    }

    #[test]
    fn a_non_numeric_cluster_name_is_reported_not_panicked() {
        let indir = tmp("badname");
        let outdir = tmp("badname_out");
        write_cluster(&indir, "abc.fastq", 1);
        let mut c = cfg(&indir, &outdir);
        c.split_mod = 2;
        assert!(plan(&indir, &c).is_err());
    }

    #[test]
    fn keep_old_skips_a_complete_output() {
        let indir = tmp("keep");
        let outdir = tmp("keep_out");
        write_cluster(&indir, "0.fastq", 3); // 12 lines
        let done = outdir.join("0");
        fs::create_dir_all(&done).unwrap();
        fs::write(done.join("corrected_reads.fastq"), "a\nb\nc\nd\n".repeat(3)).unwrap();

        let mut c = cfg(&indir, &outdir);
        c.keep_old = true;
        assert!(plan(&indir, &c).unwrap().is_empty(), "should have skipped");

        // A short (incomplete) output is recomputed.
        fs::write(done.join("corrected_reads.fastq"), "a\n").unwrap();
        assert_eq!(plan(&indir, &c).unwrap().len(), 1);
    }

    #[test]
    fn splitting_produces_chunks_of_whole_records() {
        let indir = tmp("split");
        let work = tmp("split_work");
        write_cluster(&indir, "0.fastq", 10); // 40 lines
        let split = split_clusters(&indir, &work, 2).unwrap(); // 8 lines per chunk
        let mut names: Vec<String> = fs::read_dir(&split)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect();
        names.sort();
        assert_eq!(
            names,
            vec![
                "0_0.fastq",
                "0_1.fastq",
                "0_2.fastq",
                "0_3.fastq",
                "0_4.fastq"
            ]
        );
        for n in &names {
            assert_eq!(line_count(&split.join(n)).unwrap(), 8);
        }
    }

    #[test]
    fn an_exact_multiple_leaves_no_empty_chunk() {
        let indir = tmp("exact");
        let work = tmp("exact_work");
        write_cluster(&indir, "0.fastq", 4); // 16 lines, chunk = 16
        let split = split_clusters(&indir, &work, 4).unwrap();
        let names: Vec<_> = fs::read_dir(&split).unwrap().map(|e| e.unwrap()).collect();
        assert_eq!(names.len(), 1, "no trailing empty file");
    }

    #[test]
    fn a_small_cluster_is_symlinked_as_batch_zero() {
        let indir = tmp("small");
        let work = tmp("small_work");
        write_cluster(&indir, "0.fastq", 2); // 8 lines, chunk = 4000
        let split = split_clusters(&indir, &work, 1000).unwrap();
        assert!(split.join("0_0.fastq").exists());
        assert_eq!(line_count(&split.join("0_0.fastq")).unwrap(), 8);
    }

    #[test]
    fn the_size_test_latches_after_the_first_small_cluster() {
        let indir = tmp("latch");
        let work = tmp("latch_work");
        // Cluster 0 is small, cluster 1 is large. Ascending id order means 0 is
        // seen first, latches, and 1 is symlinked rather than split --- which is
        // the reference's behaviour, bug and all.
        write_cluster(&indir, "0.fastq", 1); // 4 lines
        write_cluster(&indir, "1.fastq", 20); // 80 lines
        let split = split_clusters(&indir, &work, 2).unwrap(); // chunk = 8 lines
        assert!(split.join("1_0.fastq").exists());
        assert!(
            !split.join("1_1.fastq").exists(),
            "cluster 1 should have been symlinked whole, not split"
        );
    }

    #[test]
    fn joining_concatenates_in_lexicographic_batch_order() {
        let split = tmp("join_split");
        let outdir = tmp("join_out");
        // Eleven batches, so lexicographic order differs from numeric.
        for i in 0..11 {
            fs::write(split.join(format!("0_{i}.fastq")), "x").unwrap();
            let d = outdir.join(format!("0_{i}"));
            fs::create_dir_all(&d).unwrap();
            fs::write(d.join("corrected_reads.fastq"), format!("{i}\n")).unwrap();
        }
        join_batches(&split, &outdir, 1, 0).unwrap();
        let joined = fs::read_to_string(outdir.join("0/corrected_reads.fastq")).unwrap();
        assert_eq!(joined, "0\n1\n10\n2\n3\n4\n5\n6\n7\n8\n9\n");
        // The per-batch folders are gone.
        assert!(!outdir.join("0_0").exists());
        assert!(outdir.join("0/cat.stderr").exists());
    }
}
