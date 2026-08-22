//! Correcting one cluster, matching `isoncorrect_main`.
//!
//! This is the loop every other module hangs off: batch the reads, build the
//! minimizer and anchor indexes per batch, then walk the reads in `r_id` order
//! correcting each one.
//!
//! # The feedback loop, and why it makes reads order-dependent
//!
//! Correcting read *i* also produces corrected spans for every *other* read
//! that supported its intervals. Those land in `previously_corrected_regions`,
//! and when read *j* comes up later they do two things: they suppress anchors
//! whose whole span is already covered ([`crate::regions::filter_spans`]), and
//! they enter WIS as ready-made intervals whose payload is a literal string
//! rather than a set of segments to align. So the answer for read *j* depends on
//! which reads were corrected before it — `sorted(reads)`, i.e. ascending
//! `r_id`, is part of the contract.
//!
//! **`--exact` switches the whole mechanism off** by clearing the dict at the
//! top of every read. Since `isoncorrect_main` forces `--exact` for any cluster
//! at or below `--exact_instance_limit`, and `run_isoncorrect` defaults that to
//! 50, small clusters never exercise the feedback path at all.
//!
//! # Details that decide the answer
//!
//! * batches are contiguous `--max_seqs` chunks in `r_id` order, and
//!   `previously_corrected_regions` is **per batch** — nothing carries across;
//! * the window can be recomputed per batch under `--set_w_dynamically`, from
//!   the *batch* size, while the `--exact` decision uses the *cluster* size;
//! * `all_intervals` is sorted by `(stop, start)` with a **stable** sort, so
//!   equal keys keep the order they were appended in: span candidates first, in
//!   anchor-visit order, then carried-over regions;
//! * WIS indices come back descending and are reversed, giving intervals in
//!   ascending order for stitching;
//! * a read with no intervals, or whose intervals WIS all rejects, is emitted
//!   **uncorrected**;
//! * the emitted quality string is `'+'` repeated, not the input quality.

use std::collections::HashMap;

use crate::anchors;
use crate::corrections::{self, Correction};
use crate::fastq::Read;
use crate::guard;
use crate::minimizers;
use crate::params::Params;
use crate::quality;
use crate::regions;
use crate::support::SpanFinder;
use crate::wis;

/// One corrected read, ready to write.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Corrected {
    pub r_id: usize,
    pub acc: String,
    pub seq: Vec<u8>,
}

impl Corrected {
    /// The quality string the reference writes: `'+'` per base, not the input
    /// qualities. Deliberate in the reference, and part of the output contract.
    pub fn qual(&self) -> String {
        "+".repeat(self.seq.len())
    }
}

/// What an interval carries into the correction step.
enum Payload {
    /// Segments to align and correct, from `find_most_supported_span`.
    Instance(Vec<(usize, usize, usize)>),
    /// A span already corrected while another read was being corrected.
    Literal(Vec<u8>),
}

struct Pending {
    start: usize,
    stop: usize,
    support: usize,
    payload: Payload,
}

/// Counters the reference prints under `--verbose`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Stats {
    /// `tmp_cnt`: bounded edit distances computed.
    pub edlib_calls: usize,
    /// Reads emitted without any correction.
    pub uncorrected: usize,
    /// Correction intervals actually run.
    pub intervals: usize,
}

/// Correct a whole cluster, matching `isoncorrect_main`.
///
/// Returns one entry per input read, in `r_id` order — the order the reference
/// writes `corrected_reads.fastq` in.
pub fn correct_cluster(reads: &[Read], p: &Params) -> (Vec<Corrected>, Stats) {
    let mut out: Vec<Corrected> = Vec::with_capacity(reads.len());
    let mut stats = Stats::default();

    // The exact decision is made once, from the cluster size, before batching.
    let exact = p.exact_for_cluster(reads.len());

    for batch in reads.chunks(p.max_seqs.max(1)) {
        let w = p.window_for_batch(batch.len());

        // Sequences and quality prefix sums are indexed by `r_id - 1`, so they
        // are built over the whole cluster even though only this batch's reads
        // are corrected --- `find_most_supported_span` only ever looks up reads
        // that share an anchor, and anchors are built per batch.
        let seqs: Vec<Vec<u8>> = reads.iter().map(|r| r.seq.as_bytes().to_vec()).collect();
        let qvs: Vec<Vec<f64>> = {
            let _t = crate::profile::scope("get_qvs");
            reads
                .iter()
                .map(|r| quality::prefix_sums(r.qual.as_bytes()))
                .collect()
        };

        // `(r_id, minimizers)` in r_id order: `anchors::build` iterates it in
        // order and that order reaches the payload.
        let mins: Vec<(usize, Vec<anchors::Minimizer>)> = batch
            .iter()
            .map(|r| {
                let _t = crate::profile::scope("get_minimizers_and_positions");
                (r.r_id, minimizers::minimizers_lex(r.seq.as_bytes(), p.k, w))
            })
            .collect();
        let db = {
            let _t = crate::profile::scope("get_minimizer_combinations_database");
            // `batch.len()`, NOT `reads.len()`. The reference rebinds `reads` to
            // the current batch at the top of its batch loop, so the abundance
            // threshold in `get_minimizer_combinations_database` -- "more
            // occurrences than there are reads" -- counts the *batch*. Passing
            // the cluster size raises the threshold, so highly abundant anchors
            // that the reference drops survive, and the read gets extra
            // high-support intervals.
            //
            // Only a cluster larger than `--max_seqs` can show this, which is
            // why it took a drosophila cluster of 3 528 reads to find: it
            // changed 587 of that cluster's reads.
            anchors::build(&mins, p.k, p.xmin, p.xmax, batch.len())
        };

        // Per batch, exactly as the reference rebuilds it.
        let mut previously_corrected: HashMap<usize, Vec<Correction>> = HashMap::new();

        // See [`SpansDump`]. Off unless ISONCORRECT_SPANS_DUMP is set.
        let mut spans_dump = SpansDump::from_env();

        // Hoisted out of the read loop so its scratch buffers survive, with
        // `reset_cache` per read reproducing the reference's `already_computed =
        // {}`. Building a `SpanFinder` per read would throw the buffers away and
        // give back most of what reusing them bought.
        let mut finder = SpanFinder::new(&seqs, &qvs, p.k);

        for read in batch {
            let r_id = read.r_id;
            let seq = read.seq.as_bytes();

            if exact {
                previously_corrected.clear();
            }

            let carried = previously_corrected.remove(&r_id).unwrap_or_default();
            let spans_covered: Vec<(usize, usize)> =
                carried.iter().map(|c| (c.start, c.stop)).collect();
            let considered = regions::considered_positions(&spans_covered);
            let groups = if carried.is_empty() {
                HashMap::new()
            } else {
                regions::pos_groups(&considered)
            };

            let mut all: Vec<Pending> = Vec::new();
            // The reference resets `already_computed` per read, unconditionally.
            finder.reset_cache();
            let minimizers_of_read = &mins
                .iter()
                .find(|(id, _)| *id == r_id)
                .expect("every batch read has minimizers")
                .1;
            for ((m1, p1), partners) in anchors::comb_iterator(minimizers_of_read, p.xmin, p.xmax) {
                let kept = regions::filter_spans(&partners, p1, p.k, &considered, &groups);
                if kept.is_empty() {
                    continue;
                }
                let mut found = Vec::new();
                {
                    let _t = crate::profile::scope("find_most_supported_span");
                    finder.find(r_id, m1, p1, &kept, &db, &mut found);
                }
                spans_dump.record(r_id, m1, p1, &found);
                all.extend(found.into_iter().map(|c| Pending {
                    start: c.start,
                    stop: c.stop,
                    support: c.support,
                    payload: Payload::Instance(c.instance),
                }));
            }

            // Carried-over regions join the solver after the fresh candidates.
            all.extend(carried.into_iter().map(|c| Pending {
                start: c.start,
                stop: c.stop,
                support: c.weight,
                payload: Payload::Literal(c.seq),
            }));

            if all.is_empty() {
                stats.uncorrected += 1;
                out.push(Corrected {
                    r_id,
                    acc: read.acc.clone(),
                    seq: seq.to_vec(),
                });
                continue;
            }

            // `sort(key=lambda x: (x[1], x[0]))`, stable.
            all.sort_by_key(|iv| (iv.stop, iv.start));
            let intervals: Vec<wis::Interval> = all
                .iter()
                .map(|iv| wis::Interval {
                    start: iv.start,
                    stop: iv.stop,
                    support: iv.support,
                })
                .collect();
            let chosen = {
                let _t = crate::profile::scope("solve_WIS");
                wis::solve(&intervals)
            };
            if chosen.is_empty() {
                stats.uncorrected += 1;
                out.push(Corrected {
                    r_id,
                    acc: read.acc.clone(),
                    seq: seq.to_vec(),
                });
                continue;
            }

            // `opt_indicies[::-1]` --- WIS returns them latest-first.
            let mut to_correct: Vec<guard::Interval> = Vec::with_capacity(chosen.len());
            for &idx in chosen.iter().rev() {
                let iv = &all[idx];
                let corrected = match &iv.payload {
                    Payload::Literal(s) => s.clone(),
                    Payload::Instance(instance) => {
                        stats.intervals += 1;
                        let _t = crate::profile::scope("get_best_corrections");
                        let Some(c) = corrections::best_corrections(
                            instance,
                            &seqs,
                            p.k,
                            p.t_threshold,
                            p.max_seqs_to_spoa,
                        ) else {
                            continue;
                        };
                        // Every other supporting read gets its span recorded for
                        // later, including this read's own --- the reference does
                        // not filter it out either.
                        for (other, regions_of) in c.others {
                            previously_corrected
                                .entry(other)
                                .or_default()
                                .extend(regions_of);
                        }
                        c.curr_read
                    }
                };
                to_correct.push(guard::Interval {
                    start: iv.start,
                    stop: iv.stop,
                    corrected,
                });
            }

            let corrected = if to_correct.is_empty() {
                stats.uncorrected += 1;
                seq.to_vec()
            } else {
                let _t = crate::profile::scope("correct_read");
                guard::correct_read(seq, &to_correct)
            };
            out.push(Corrected {
                r_id,
                acc: read.acc.clone(),
                seq: corrected,
            });
        }
        // Once per batch, not once per read: `finder` now spans the batch, so its
        // counter is cumulative and adding it inside the loop would double-count.
        stats.edlib_calls += finder.edlib_calls;
        spans_dump.flush();
    }

    (out, stats)
}

/// The `spans.tsv` oracle, recorded from the **live driver** rather than the
/// dump binary.
///
/// This closes the one real hole in the port's stage coverage.
/// `isoncorrect-dump` can emit spans, but it does not correct, so
/// `previously_corrected_regions` is always empty there and only the `--exact`
/// trajectory is comparable. The anchor filtering that consumes those regions
/// (`regions::filter_spans`, `considered_positions`, `pos_groups`) is therefore
/// *inert* in every dump ever taken, and had no oracle at all — it was covered
/// only by end-to-end output on a 100-read fixture where it barely fires.
///
/// Recording here instead makes the non-exact trajectory directly comparable to
/// `bench/dump_reference.py`'s `spans.tsv`, which has always recorded the live
/// reference. Format is identical, so the two files diff.
///
/// ```bash
/// ISONCORRECT_SPANS_DUMP=/tmp/rs_spans.tsv isONcorrect --fastq c.fastq --outfolder /tmp/x
/// diff /tmp/py/spans.tsv /tmp/rs_spans.tsv
/// ```
struct SpansDump {
    path: Option<String>,
    rows: Vec<String>,
}

impl SpansDump {
    fn from_env() -> Self {
        Self {
            path: std::env::var("ISONCORRECT_SPANS_DUMP")
                .ok()
                .filter(|p| !p.is_empty()),
            rows: Vec::new(),
        }
    }

    fn record(&mut self, r_id: usize, m1: &[u8], p1: usize, found: &[crate::support::Candidate]) {
        if self.path.is_none() {
            return;
        }
        for c in found {
            let payload: Vec<String> = c
                .instance
                .iter()
                .flat_map(|&(r, a, b)| [r.to_string(), a.to_string(), b.to_string()])
                .collect();
            self.rows.push(format!(
                "{r_id}\t{}\t{p1}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(m1),
                c.start,
                c.stop,
                c.support,
                payload.join(",")
            ));
        }
    }

    /// Appended per batch, so a multi-batch cluster produces one file in read
    /// order — the order the reference's list is built in.
    fn flush(&mut self) {
        let Some(path) = self.path.clone() else {
            return;
        };
        use std::io::Write;
        let existed = std::path::Path::new(&path).exists();
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .expect("spans dump path is writable");
        if !existed {
            writeln!(f, "#r_id\tm1\tp1\tstart\tstop\tweight\tinstance").expect("spans dump header");
        }
        for row in self.rows.drain(..) {
            writeln!(f, "{row}").expect("spans dump row");
        }
    }
}

/// Write `corrected_reads.fastq`, matching the reference's format exactly.
pub fn write_fastq<W: std::io::Write>(out: &mut W, corrected: &[Corrected]) -> std::io::Result<()> {
    for c in corrected {
        writeln!(out, "@{}", c.acc)?;
        out.write_all(&c.seq)?;
        writeln!(out)?;
        writeln!(out, "+")?;
        writeln!(out, "{}", c.qual())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> Params {
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

    fn read(r_id: usize, seq: &str) -> Read {
        Read {
            r_id,
            acc: format!("read{r_id}"),
            seq: seq.to_string(),
            qual: "I".repeat(seq.len()),
        }
    }

    #[test]
    fn a_cluster_with_nothing_to_correct_is_returned_verbatim() {
        // Three reads too short to share anchors: no intervals, no corrections.
        let reads = vec![read(1, "ACGTACGTAC"), read(2, "ACGTACGTAC")];
        let (out, stats) = correct_cluster(&reads, &params());
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].seq, b"ACGTACGTAC");
        assert_eq!(stats.uncorrected, 2);
    }

    #[test]
    fn output_is_in_r_id_order_across_batches() {
        let seq = "ACGTTGCAACGTTGCAACGTACGTTGCAACGTTGCAACGT";
        let reads: Vec<Read> = (1..=5).map(|i| read(i, seq)).collect();
        let p = Params {
            max_seqs: 2,
            ..params()
        };
        let (out, _) = correct_cluster(&reads, &p);
        assert_eq!(
            out.iter().map(|c| c.r_id).collect::<Vec<_>>(),
            vec![1, 2, 3, 4, 5]
        );
    }

    #[test]
    fn identical_reads_come_back_unchanged() {
        // A real cluster of identical reads: the consensus is the read, so
        // whatever gets corrected must reproduce it exactly.
        let seq = "ACGTTGCAACCTGGATCAGGACTTGCAACGTACGTTGCAACGTTGCAACGTAGGCTTAACCGGTT";
        let reads: Vec<Read> = (1..=6).map(|i| read(i, seq)).collect();
        let (out, _) = correct_cluster(&reads, &params());
        for c in &out {
            assert_eq!(c.seq, seq.as_bytes(), "read {} changed", c.r_id);
        }
    }

    #[test]
    fn the_quality_string_is_plus_per_base() {
        let c = Corrected {
            r_id: 1,
            acc: "x".into(),
            seq: b"ACGT".to_vec(),
        };
        assert_eq!(c.qual(), "++++");
        let mut buf = Vec::new();
        write_fastq(&mut buf, &[c]).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "@x\nACGT\n+\n++++\n");
    }

    #[test]
    fn an_empty_cluster_produces_no_output() {
        let (out, _) = correct_cluster(&[], &params());
        assert!(out.is_empty());
    }
}

/// End-to-end oracle: the reference's own `corrected_reads.fastq`.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// CLUSTER_FASTQ=cluster.fastq EXPECTED_FASTQ=/tmp/d/_run/corrected_reads.fastq \
///   cargo test --manifest-path rust/Cargo.toml driver::oracle -- --nocapture
/// ```
///
/// Parameters come from `CLUSTER_ARGS`, a space-separated `--k 9 --w 10`-style
/// string, so a case can be replayed at the settings it was recorded with.
#[cfg(test)]
mod oracle {
    use super::*;
    use crate::fastq;

    fn params_from_env() -> Params {
        let mut p = Params {
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
        };
        let args = std::env::var("CLUSTER_ARGS").unwrap_or_default();
        let f: Vec<&str> = args.split_whitespace().collect();
        let mut i = 0;
        while i < f.len() {
            let value = |i: usize| f[i + 1].parse().expect("numeric flag value");
            match f[i] {
                "--k" => {
                    p.k = value(i);
                    i += 2;
                }
                "--w" => {
                    p.w = value(i);
                    i += 2;
                }
                "--xmin" => {
                    p.xmin = value(i);
                    i += 2;
                }
                "--xmax" => {
                    p.xmax = value(i);
                    i += 2;
                }
                "--T" => {
                    p.t_threshold = f[i + 1].parse().unwrap();
                    i += 2;
                }
                "--max_seqs" => {
                    p.max_seqs = value(i);
                    i += 2;
                }
                "--max_seqs_to_spoa" => {
                    p.max_seqs_to_spoa = value(i);
                    i += 2;
                }
                "--exact_instance_limit" => {
                    p.exact_instance_limit = value(i);
                    i += 2;
                }
                "--exact" => {
                    p.exact = true;
                    i += 1;
                }
                "--set_w_dynamically" => {
                    p.set_w_dynamically = true;
                    i += 1;
                }
                other => panic!("unhandled flag {other} in CLUSTER_ARGS"),
            }
        }
        // `main` clamps xmin up to 2k before anything runs; the dump binary has
        // to mirror this and so does the oracle, or the comparison lies.
        if p.xmin < 2 * p.k {
            p.xmin = 2 * p.k;
        }
        p
    }

    #[test]
    fn matches_the_reference_output_byte_for_byte() {
        let (Ok(cluster), Ok(expected)) = (
            std::env::var("CLUSTER_FASTQ"),
            std::env::var("EXPECTED_FASTQ"),
        ) else {
            return;
        };

        let file = std::fs::File::open(&cluster).expect("cluster fastq");
        let reads = fastq::read_fastq(std::io::BufReader::new(file)).expect("parse");
        let (corrected, stats) = correct_cluster(&reads, &params_from_env());

        let mut got = Vec::new();
        write_fastq(&mut got, &corrected).unwrap();
        let want = std::fs::read(&expected).expect("expected fastq");

        if got != want {
            // Locate the first differing record rather than dumping megabytes.
            let g: Vec<&[u8]> = got.split(|&b| b == b'\n').collect();
            let wv: Vec<&[u8]> = want.split(|&b| b == b'\n').collect();
            let mut shown = 0;
            for (i, (a, b)) in g.iter().zip(&wv).enumerate() {
                if a != b {
                    eprintln!(
                        "line {i} (record {}):\n  python: {}\n  rust  : {}",
                        i / 4,
                        String::from_utf8_lossy(&b[..b.len().min(120)]),
                        String::from_utf8_lossy(&a[..a.len().min(120)])
                    );
                    shown += 1;
                    if shown == 3 {
                        break;
                    }
                }
            }
            eprintln!("{} vs {} lines, stats {stats:?}", g.len(), wv.len());
        }
        assert_eq!(got, want, "corrected_reads.fastq differed");
        eprintln!(
            "byte-identical on {} reads ({} uncorrected, {} intervals corrected)",
            corrected.len(),
            stats.uncorrected,
            stats.intervals
        );
    }
}
