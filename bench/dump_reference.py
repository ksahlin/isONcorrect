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

Whole-run oracles, written once at the end rather than per batch:

    cigar_cases.tsv           edlib task="path" query, target, distance, CIGAR
    cigar_to_seq_cases.tsv    CIGAR expanded into its two gapped strings
    spoa_cases.tsv            consensus and the sequences it came from
    msa_cases.tsv             pairwise alignments in, multialignment rows out
    spans.tsv                 what find_most_supported_span appended
    wis_input.tsv/.tsv        solve_WIS inputs and chosen indices

Every file is emitted in the reference's own iteration order, NOT sorted:
Python dicts are insertion-ordered and that order feeds later stages, so the
port has to reproduce it. Diff these files directly.
"""

from __future__ import annotations

import argparse
import os
import sys

# Pin the hash seed before anything imports the reference.
#
# `os.environ["PYTHONHASHSEED"] = ...` inside a running interpreter does
# NOTHING --- CPython reads the variable once at startup and the string hash key
# is already fixed by the time any Python code runs. This file used to do
# exactly that, which made the dumps look pinned while they were in fact
# randomly seeded on every run.
#
# That matters more than it sounds: `get_alternative_ref_contexts` returns a
# `set`, and `get_best_corrections` iterates it, so the reference's *output*
# depends on the seed. See PORTING.md and `check_seed_sensitivity.py`.
# Re-exec is the only way to fix it from inside the script.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

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

    # Capture edlib's task="path" calls: query, target, and the CIGAR it chose.
    # The CIGAR is NOT unique -- many alignments share the optimal score -- so
    # this is the oracle for deciding whether any Rust aligner can replace it.
    import edlib as _edlib

    orig_align = _edlib.align
    cigar_cases = {"rows": []}

    def wrapped_align(query, target, *a, **kw):
        # The reference calls edlib both positionally --- align(x, y, "NW",
        # 'dist', k) --- and by keyword, so accept either shape.
        res = orig_align(query, target, *a, **kw)
        task = kw.get("task", a[1] if len(a) > 1 else None)
        if task == "path" and res.get("cigar"):
            cigar_cases["rows"].append(
                f"{query}\t{target}\t{res['editDistance']}\t{res['cigar']}"
            )
        return res

    _edlib.align = wrapped_align

    # Capture what actually goes into spoa and what comes back: the oracle for
    # deciding whether a Rust POA reproduces spoa's consensus exactly.
    from isoncorrect import create_augmented_reference as car

    orig_spoa = car.run_spoa
    spoa_cases = {"rows": []}

    def wrapped_spoa(reads_path, spoa_out_file, spoa_path):
        with open(reads_path) as fh:
            seqs = [l.strip() for l in fh if not l.startswith(">")]
        consensus = orig_spoa(reads_path, spoa_out_file, spoa_path)
        spoa_cases["rows"].append("\t".join([consensus, *seqs]))
        return consensus

    car.run_spoa = wrapped_spoa
    ion.create_augmented_reference.run_spoa = wrapped_spoa

    # cigar_to_seq expands a CIGAR into the two gapped strings. Deterministic
    # given its inputs, so recording both sides makes it directly replayable.
    from isoncorrect import help_functions as hf

    orig_cigar_to_seq = hf.cigar_to_seq
    cigar_seq_cases = {"rows": []}

    def wrapped_cigar_to_seq(cigar, query, ref):
        q_aln, r_aln = orig_cigar_to_seq(cigar, query, ref)
        cigar_seq_cases["rows"].append(f"{cigar}\t{query}\t{ref}\t{q_aln}\t{r_aln}")
        return q_aln, r_aln

    hf.cigar_to_seq = wrapped_cigar_to_seq
    ion.help_functions.cigar_to_seq = wrapped_cigar_to_seq

    # The multialignment matrix stacks the pairwise alignments into shared
    # columns. Several non-unique choices land here at once --- edlib's CIGAR
    # tie-break, which insertion is widest, and get_best_solution's four-way
    # fallthrough --- so real cases are worth much more than constructed ones.
    #
    # The partition's tuple keys are not recorded: the matrix depends only on
    # the pairwise alignments and their order. Note create_multialignment_matrix
    # calls position_query_to_alignment(s_alignment, m_alignment), i.e. query
    # first, target second, which is the reverse of the tuple's own order.
    from isoncorrect import correct_seqs as cs

    orig_mam = cs.create_multialignment_matrix
    msa_cases = {"rows": [], "n": 0}

    def wrapped_mam(partition):
        matrix = orig_mam(partition)
        case = msa_cases["n"]
        msa_cases["n"] += 1
        for i, q_acc in enumerate(partition):
            (_ed, m_alignment, s_alignment, _degree) = partition[q_acc]
            row = "".join(matrix[q_acc])
            msa_cases["rows"].append(
                f"{case}\t{i}\t{s_alignment}\t{m_alignment}\t{row}"
            )
        return matrix

    cs.create_multialignment_matrix = wrapped_mam
    ion.correct_seqs.create_multialignment_matrix = wrapped_mam

    # Context windows and alternative references. get_contexts is called
    # immediately before get_alternative_ref_contexts and its k_size is not
    # passed on, so it is picked up here and paired by call order.
    orig_contexts = ion.get_contexts
    orig_alt = ion.get_alternative_ref_contexts
    ctx_cases = {"rows": [], "n": 0, "k_size": None, "alts": 0}

    def wrapped_contexts(alignment_matrix, k_size):
        ctx_cases["k_size"] = k_size
        return orig_contexts(alignment_matrix, k_size)

    def wrapped_alt(alignment_matrix, contexts_per_pos, context_threshold):
        result = orig_alt(alignment_matrix, contexts_per_pos, context_threshold)
        case = ctx_cases["n"]
        ctx_cases["n"] += 1
        rows = ctx_cases["rows"]
        ref_aln = "".join(alignment_matrix["ref"])
        rows.append(f"C\t{case}\t{ctx_cases['k_size']}\t{context_threshold!r}\t{ref_aln}")
        # Every row, consensus included: test_numba counts A over all of them.
        for aln in alignment_matrix.values():
            rows.append(f"R\t{case}\t{''.join(aln)}")
        windows = ",".join(f"{a}:{b}" for a, b in contexts_per_pos)
        rows.append(f"W\t{case}\t{windows}")
        for col, alts in enumerate(result):
            for variant, aln_seq, depth, threshold in alts:
                rows.append(
                    f"A\t{case}\t{col}\t{variant}\t{''.join(aln_seq)}\t{depth}\t{threshold!r}"
                )
                ctx_cases["alts"] += 1
        return result

    ion.get_contexts = wrapped_contexts
    ion.get_alternative_ref_contexts = wrapped_alt

    # get_best_corrections end to end: the interval's segments in, the corrected
    # spans out. This is the only oracle that covers the PFM, which the
    # reference never returns and which cannot be captured without
    # reimplementing it here.
    #
    # Segments are recorded rather than whole reads: `seq` below is exactly what
    # the function reads out of `reads`, so a case replays without the fastq.
    orig_best = ion.get_best_corrections
    corr_cases = {"rows": [], "n": 0}

    def wrapped_best(curr_best_seqs, reads_, k_size, work_dir, *rest, **kw):
        curr, others = orig_best(curr_best_seqs, reads_, k_size, work_dir, *rest, **kw)
        v_depth = rest[0] if rest else kw.get("v_depth_ratio_threshold", 0.1)
        max_spoa = rest[1] if len(rest) > 1 else kw.get("max_seqs_to_spoa", 200)
        case = corr_cases["n"]
        corr_cases["n"] += 1
        rows = corr_cases["rows"]
        rows.append(f"C\t{case}\t{k_size}\t{v_depth!r}\t{max_spoa}\t{curr}")
        for q_id, pos1, pos2 in ion.grouper(curr_best_seqs, 3):
            seq = reads_[q_id][1][pos1 : pos2 + k_size]
            rows.append(f"S\t{case}\t{q_id}\t{pos1}\t{pos2}\t{seq}")
        for q_id, regions in others.items():
            for start, stop, weight, corrected in regions:
                rows.append(f"O\t{case}\t{q_id}\t{start}\t{stop}\t{weight}\t{corrected}")
        return curr, others

    ion.get_best_corrections = wrapped_best

    # The structural-overcorrection guard, in three pieces.
    #
    # parasail's alignment is the one part of this that is NOT uniquely
    # defined --- affine gaps with free end gaps leave many optimal alignments,
    # and parasail picks one --- so it gets its own oracle, exactly like edlib's
    # task="path" CIGAR did.
    orig_parasail = ion.parasail_alignment
    parasail_cases = {"rows": []}

    def wrapped_parasail(s1, s2, **kw):
        s1_aln, s2_aln, cigar, tuples, score = orig_parasail(s1, s2, **kw)
        parasail_cases["rows"].append(
            "\t".join(
                [
                    s1,
                    s2,
                    cigar,
                    str(score),
                    str(kw.get("match_score", 2)),
                    str(kw.get("mismatch_penalty", -2)),
                    str(kw.get("opening_penalty", 24)),
                    str(kw.get("gap_ext", 1)),
                ]
            )
        )
        return s1_aln, s2_aln, cigar, tuples, score

    # fix_correction is a pure function of the two aligned strings, so recording
    # both sides makes it replayable without parasail in the loop.
    orig_fix = ion.fix_correction
    fix_cases = {"rows": []}

    def wrapped_fix(orig, corr):
        adjusted = orig_fix(orig, corr)
        fix_cases["rows"].append(f"{orig}\t{corr}\t{adjusted}")
        return adjusted

    # correct_read stitches the corrected intervals back into the read and then
    # runs the guard. Intervals are recorded by *source* rather than by value:
    # `c:<n>` means the n-th get_best_corrections call (whose output is in
    # correction_cases.tsv), `s:<literal>` a region carried over from
    # previously_corrected_regions. That keeps the harness from reimplementing
    # any of the reference's logic.
    orig_correct_read = ion.correct_read
    read_cases = {"rows": [], "n": 0}

    def wrapped_correct_read(seq, reads_, intervals, *rest):
        case = read_cases["n"]
        read_cases["n"] += 1
        next_call = corr_cases["n"]
        rows = read_cases["rows"]
        rows.append(f"R\t{case}\t{seq}")
        for idx, (start, stop, _weights, instance) in enumerate(intervals):
            if isinstance(instance, str):
                source = f"s:{instance}"
            else:
                source = f"c:{next_call}"
                next_call += 1
            rows.append(f"I\t{case}\t{idx}\t{start}\t{stop}\t{source}")
        adjusted, others = orig_correct_read(seq, reads_, intervals, *rest)
        rows.append(f"A\t{case}\t{adjusted}")
        return adjusted, others

    ion.parasail_alignment = wrapped_parasail
    ion.fix_correction = wrapped_fix
    ion.correct_read = wrapped_correct_read

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

    with open(os.path.join(args.outdir, "cigar_cases.tsv"), "w") as fh:
        fh.write("#query\ttarget\ted\tcigar\n")
        fh.write("\n".join(cigar_cases["rows"]))
        fh.write("\n" if cigar_cases["rows"] else "")
    BATCH["meta"].append(f"edlib path-calls={len(cigar_cases['rows'])}")

    with open(os.path.join(args.outdir, "spoa_cases.tsv"), "w") as fh:
        fh.write("#consensus\tseq...\n")
        fh.write("\n".join(spoa_cases["rows"]))
        fh.write("\n" if spoa_cases["rows"] else "")
    BATCH["meta"].append(f"spoa calls={len(spoa_cases['rows'])}")

    with open(os.path.join(args.outdir, "cigar_to_seq_cases.tsv"), "w") as fh:
        fh.write("#cigar\tquery\tref\tquery_aligned\tref_aligned\n")
        fh.write("\n".join(cigar_seq_cases["rows"]))
        fh.write("\n" if cigar_seq_cases["rows"] else "")
    BATCH["meta"].append(f"cigar_to_seq calls={len(cigar_seq_cases['rows'])}")

    with open(os.path.join(args.outdir, "msa_cases.tsv"), "w") as fh:
        fh.write("#case\trow\tquery_aligned\ttarget_aligned\taln_row\n")
        fh.write("\n".join(msa_cases["rows"]))
        fh.write("\n" if msa_cases["rows"] else "")
    BATCH["meta"].append(
        f"multialignment matrices={msa_cases['n']} rows={len(msa_cases['rows'])}"
    )

    with open(os.path.join(args.outdir, "context_cases.tsv"), "w") as fh:
        fh.write("#C case k_size threshold ref_aln / R case row / W case start:stop,... / A case col variant context depth threshold\n")
        fh.write("\n".join(ctx_cases["rows"]))
        fh.write("\n" if ctx_cases["rows"] else "")
    BATCH["meta"].append(
        f"alternative-context calls={ctx_cases['n']} alternatives={ctx_cases['alts']}"
    )

    with open(os.path.join(args.outdir, "correction_cases.tsv"), "w") as fh:
        fh.write("#C case k_size T max_seqs_to_spoa curr / S case q_id pos1 pos2 segment / O case q_id start stop weight corrected\n")
        fh.write("\n".join(corr_cases["rows"]))
        fh.write("\n" if corr_cases["rows"] else "")
    BATCH["meta"].append(f"get_best_corrections calls={corr_cases['n']}")

    with open(os.path.join(args.outdir, "parasail_cases.tsv"), "w") as fh:
        fh.write("#s1\ts2\tcigar\tscore\tmatch\tmismatch\topen\text\n")
        fh.write("\n".join(parasail_cases["rows"]))
        fh.write("\n" if parasail_cases["rows"] else "")
    BATCH["meta"].append(f"parasail calls={len(parasail_cases['rows'])}")

    with open(os.path.join(args.outdir, "fix_correction_cases.tsv"), "w") as fh:
        fh.write("#orig_aligned\tcorr_aligned\tadjusted\n")
        fh.write("\n".join(fix_cases["rows"]))
        fh.write("\n" if fix_cases["rows"] else "")
    BATCH["meta"].append(f"fix_correction calls={len(fix_cases['rows'])}")

    with open(os.path.join(args.outdir, "correct_read_cases.tsv"), "w") as fh:
        fh.write("#R case seq / I case idx start stop c:<call>|s:<literal> / A case adjusted\n")
        fh.write("\n".join(read_cases["rows"]))
        fh.write("\n" if read_cases["rows"] else "")
    BATCH["meta"].append(f"correct_read calls={read_cases['n']}")

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
