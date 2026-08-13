#!/usr/bin/env python3
"""Dump the reference's intermediate structures, stage by stage.

End-to-end equivalence tells you *that* the port is wrong, never *where*. This
writes the structures each porting stage produces, in a stable line-oriented
format, so a Rust stage can be diffed in isolation long before the whole
pipeline runs.

It does not modify `src/isoncorrect/` --- the reference stays normative. Like
`profile_python.py` it wraps the functions of interest at import time.

Must run inside the pinned reference environment:

    conda activate isoncorrect-ref
    bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/dump

Extra isONcorrect flags pass straight through:

    bench/dump_reference.py --fastq c.fastq --outdir /tmp/dump -- --k 9 --w 10

Output, one directory per batch (`--max_seqs` chunk):

    batch_000/reads.tsv       r_id, acc, seq, qual   -- parsing and r_id order
    batch_000/minimizers.tsv  r_id, idx, kmer, pos   -- lex minimizers, ties last
    batch_000/anchors.tsv     m1, m2, r_id, p1, p2   -- the minimizer-pair index
    batch_000/anchor_keys.tsv m1, m2, n_entries      -- key set and multiplicity
    meta.txt                  parameters and totals

Every file is emitted in the reference's own iteration order, NOT sorted:
Python dicts are insertion-ordered and that order feeds later stages, so the
port has to reproduce it. Diff these files directly.
"""

from __future__ import annotations

import argparse
import os
import sys

BATCH = {"n": 0, "outdir": None, "meta": []}


def _batch_dir():
    d = os.path.join(BATCH["outdir"], f"batch_{BATCH['n']:03d}")
    os.makedirs(d, exist_ok=True)
    return d


def dump_reads(reads):
    """reads: {r_id: (acc, seq, qual)} in insertion order (1-based r_id)."""
    path = os.path.join(_batch_dir(), "reads.tsv")
    with open(path, "w") as fh:
        fh.write("#r_id\tacc\tseq\tqual\n")
        for r_id, (acc, seq, qual) in reads.items():
            fh.write(f"{r_id}\t{acc}\t{seq}\t{qual}\n")
    return len(reads)


def dump_minimizers(m):
    """M: {r_id: [(kmer, pos), ...]}."""
    path = os.path.join(_batch_dir(), "minimizers.tsv")
    total = 0
    with open(path, "w") as fh:
        fh.write("#r_id\tidx\tkmer\tpos\n")
        for r_id, mins in m.items():
            for i, (kmer, pos) in enumerate(mins):
                fh.write(f"{r_id}\t{i}\t{kmer}\t{pos}\n")
                total += 1
    return total


def dump_anchors(m2):
    """M2: {m1: {m2: array('I', [r_id, p1, p2, ...])}} after filtering.

    Two files. `anchors.tsv` is the payload; `anchor_keys.tsv` records the key
    set with its entry count.

    Empty keys are skipped in both files. They exist only as an artefact: the
    filtering loop does `del M2[m1][m2]` and then immediately reads
    `M2[m1][m2]`, which on a defaultdict *recreates* the key with an empty
    array, so deleted singletons come back as zero-entry keys. A lookup on one
    is indistinguishable from a lookup on an absent key --- both mean "no
    support" --- so the port stores only non-empty entries, and dumping them
    here would make every diff fail for no reason. The count is reported in
    `meta.txt` instead, since it is the memory story.
    """
    d = _batch_dir()
    n_triples = 0
    n_keys = 0
    n_empty = 0
    with open(os.path.join(d, "anchors.tsv"), "w") as fh, open(
        os.path.join(d, "anchor_keys.tsv"), "w"
    ) as kh:
        fh.write("#m1\tm2\tr_id\tp1\tp2\n")
        kh.write("#m1\tm2\tn_entries\n")
        for m1 in m2:
            for m2k in m2[m1]:
                payload = m2[m1][m2k]
                n = len(payload) // 3
                if n == 0:
                    n_empty += 1
                    continue
                kh.write(f"{m1}\t{m2k}\t{n}\n")
                n_keys += 1
                for i in range(0, len(payload), 3):
                    fh.write(
                        f"{m1}\t{m2k}\t{payload[i]}\t{payload[i + 1]}\t{payload[i + 2]}\n"
                    )
                    n_triples += 1
    return n_keys, n_triples, n_empty


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Dump the isONcorrect reference's intermediate structures",
        epilog="Flags after `--` are forwarded verbatim to isONcorrect.",
    )
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("passthrough", nargs="*")
    args = ap.parse_args()

    # The default path is deterministic; pin the seed anyway so dumps taken at
    # different times are comparable. See PORTING.md.
    os.environ.setdefault("PYTHONHASHSEED", "0")

    try:
        from isoncorrect import isONcorrect as ion
    except ImportError as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("Run inside the reference env: conda activate isoncorrect-ref", file=sys.stderr)
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    BATCH["outdir"] = args.outdir

    orig_minpos = ion.get_minimizers_and_positions
    orig_db = ion.get_minimizer_combinations_database

    def wrapped_minpos(reads, w, k, hash_fcn):
        # First call of a batch: `reads` is the batch dict, so this is also the
        # only place the per-batch read set is visible.
        n_reads = dump_reads(reads)
        result = orig_minpos(reads, w, k, hash_fcn)
        n_mins = dump_minimizers(result)
        BATCH["meta"].append(
            f"batch {BATCH['n']:03d}: reads={n_reads} w={w} k={k} "
            f"hash_fcn={hash_fcn} minimizers={n_mins}"
        )
        return result

    def wrapped_db(reads, m, k, x_low, x_high):
        result = orig_db(reads, m, k, x_low, x_high)
        n_keys, n_triples, n_empty = dump_anchors(result)
        BATCH["meta"].append(
            f"batch {BATCH['n']:03d}: anchor_keys={n_keys} entries={n_triples} "
            f"empty_keys={n_empty} xmin={x_low} xmax={x_high}"
        )
        BATCH["n"] += 1
        return result

    # solve_WIS is pure: intervals in, chosen indices out. Capturing both sides
    # of every call lets the Rust implementation be replayed against the exact
    # same inputs, which is far more searching than any synthetic test --- the
    # weights are floats and the tie-breaks are strict `>`.
    orig_wis = ion.solve_WIS
    wis_calls = {"n": 0, "inp": [], "out": []}

    def wrapped_wis(intervals):
        call = wis_calls["n"]
        for j, (start, stop, w, _instance) in enumerate(intervals):
            wis_calls["inp"].append(f"{call}\t{j}\t{start}\t{stop}\t{w}")
        result = orig_wis(intervals)
        for idx in result:
            wis_calls["out"].append(f"{call}\t{idx}")
        wis_calls["n"] += 1
        return result

    # find_most_supported_span appends to `all_intervals` in place. Recording
    # what each call appends, tagged with the anchor that produced it, pins both
    # the candidate intervals and the order supporting segments were collected
    # in --- the latter reaches the consensus, which is order-sensitive.
    orig_span = ion.find_most_supported_span
    spans = {"rows": []}

    def wrapped_span(
        r_id, m1, p1, m1_curr_spans, db, reads_, all_intervals, k_size, *rest
    ):
        before = len(all_intervals)
        result = orig_span(
            r_id, m1, p1, m1_curr_spans, db, reads_, all_intervals, k_size, *rest
        )
        for start, stop, weight, instance in all_intervals[before:]:
            payload = ",".join(str(x) for x in instance)
            spans["rows"].append(
                f"{r_id}\t{m1}\t{p1}\t{start}\t{stop}\t{weight}\t{payload}"
            )
        return result

    ion.get_minimizers_and_positions = wrapped_minpos
    ion.get_minimizer_combinations_database = wrapped_db
    ion.solve_WIS = wrapped_wis
    ion.find_most_supported_span = wrapped_span

    argv = ["isONcorrect", "--fastq", args.fastq, "--outfolder", os.path.join(args.outdir, "_run")]
    argv += [a for a in args.passthrough if a != "--"]
    sys.argv = argv
    print(f"==> dumping: {' '.join(argv)}")

    try:
        ion.main()
    except SystemExit as exc:
        if exc.code not in (0, None):
            print(f"isONcorrect exited with {exc.code}", file=sys.stderr)
            return int(exc.code)

    with open(os.path.join(args.outdir, "spans.tsv"), "w") as fh:
        fh.write("#r_id\tm1\tp1\tstart\tstop\tweight\tinstance\n")
        fh.write("\n".join(spans["rows"]))
        fh.write("\n" if spans["rows"] else "")
    BATCH["meta"].append(f"find_most_supported_span intervals={len(spans['rows'])}")

    with open(os.path.join(args.outdir, "wis_input.tsv"), "w") as fh:
        fh.write("#call\tj\tstart\tstop\tw\n")
        fh.write("\n".join(wis_calls["inp"]))
        fh.write("\n" if wis_calls["inp"] else "")
    with open(os.path.join(args.outdir, "wis_output.tsv"), "w") as fh:
        fh.write("#call\topt_index\n")
        fh.write("\n".join(wis_calls["out"]))
        fh.write("\n" if wis_calls["out"] else "")
    BATCH["meta"].append(
        f"solve_WIS calls={wis_calls['n']} intervals={len(wis_calls['inp'])} "
        f"selected={len(wis_calls['out'])}"
    )

    with open(os.path.join(args.outdir, "meta.txt"), "w") as fh:
        fh.write(f"argv: {' '.join(argv)}\n")
        fh.write(f"batches: {BATCH['n']}\n")
        for line in BATCH["meta"]:
            fh.write(line + "\n")

    print()
    print(f"==> wrote {BATCH['n']} batch dump(s) to {args.outdir}")
    for line in BATCH["meta"]:
        print("   ", line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
