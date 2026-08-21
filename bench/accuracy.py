#!/usr/bin/env python3
"""Per-read accuracy against a known reference transcriptome.

Equivalence testing answers "do two implementations agree". It cannot answer
"is the correction any good", and once the port deliberately diverges from the
reference — the guard's aligner reports a different equally-optimal alignment,
changing ~0.8% of reads — that second question is the one that matters. This
measures it.

Truth comes from one of two places:

* **simulated reads** carry their source transcript in the header
  (`@read_1_from_SIRV612`), so no search is needed;
* **real reads** are assigned by aligning against every transcript and taking
  the best match, ties broken by transcript name so the choice is deterministic.

Either way the read is then scored against its assigned transcript with edlib's
infix mode (`HW`), which finds the read's best placement inside it, and the error
rate is `edit_distance / len(read)`.

**The assignment is made once and reused for every set** — from the first set, or
whichever `--assign-from` names. That matters: if each implementation picked its
own best-matching transcript, every one would be flattered by construction, and
a correction that dragged a read toward the wrong isoform would score *better*
rather than worse. **Assign from the uncorrected reads**, which are neutral
between implementations.

`--pair-against` chooses the baseline for the per-read comparison separately, so
truth can stay neutral while the comparison is made against the reference
implementation.

Corrected reads keep their input headers, so the same mapping works for every
implementation's output.

    conda activate isoncorrect-ref
    bench/accuracy.py --transcriptome sirv_transcriptome.fasta \\
        --reads uncorrected=/tmp/folder_corpus \\
        --reads python=/tmp/fold_py \\
        --reads rust-exact=/tmp/fold_ex \\
        --reads rust-block=/tmp/fold_ba

Each `--reads NAME=PATH` is either a folder (searched for `corrected_reads.fastq`
and `*.fastq`) or a single fastq. Reads are matched across sets by header, and
only headers present in *every* set are scored, so the comparison is paired.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import statistics
import sys

HEADER_TRANSCRIPT = re.compile(r"_from_(\S+)")


def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_fastqs(path: str) -> dict[str, str]:
    """Header -> sequence, over a fastq or a folder of them."""
    if os.path.isfile(path):
        files = [path]
    else:
        files = sorted(glob.glob(os.path.join(path, "**", "corrected_reads.fastq"), recursive=True))
        if not files:
            files = sorted(glob.glob(os.path.join(path, "*.fastq")))
    out: dict[str, str] = {}
    for f in files:
        with open(f) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                out[header[1:].strip()] = seq
    return out


# Set once per worker process by `_init`. macOS spawns rather than forks, so the
# transcriptome is pickled to each worker once instead of being inherited.
_TRANSCRIPTS: dict[str, str] = {}
_NAMES: list[str] = []


def _init(transcripts: dict[str, str]) -> None:
    global _TRANSCRIPTS, _NAMES
    _TRANSCRIPTS = transcripts
    _NAMES = sorted(transcripts)


def _score_chunk(chunk):
    """One chunk of `(header, forced_truth_or_None, [seq_per_set])`.

    Assignment and scoring happen in the same worker so each read's sequences
    cross the process boundary exactly once. Returns
    `(header, truth_name, assign_err, [error_per_set])`.
    """
    import edlib

    out = []
    for header, forced, seqs in chunk:
        assign_seq = seqs[0]
        if forced is not None:
            truth_name, assign_err = forced, 0.0
        elif not assign_seq:
            continue
        else:
            best_name, best_ed = None, None
            for tname in _NAMES:
                ed = edlib.align(
                    assign_seq, _TRANSCRIPTS[tname], mode="HW", task="distance"
                )["editDistance"]
                if best_ed is None or ed < best_ed:
                    best_name, best_ed = tname, ed
            if best_name is None:
                continue
            truth_name = best_name
            assign_err = 100.0 * best_ed / len(assign_seq)

        truth = _TRANSCRIPTS[truth_name]
        errors = []
        for seq in seqs:
            if not seq:
                errors.append(None)
                continue
            ed = edlib.align(seq, truth, mode="HW", task="distance")["editDistance"]
            errors.append((ed, len(seq)))
        out.append((header, truth_name, assign_err, errors))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--transcriptome", required=True)
    ap.add_argument(
        "--reads",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="a labelled read set; repeat for each implementation",
    )
    ap.add_argument("--max-reads", type=int, default=0, help="sample this many (0 = all)")
    ap.add_argument(
        "--assign-from",
        default=None,
        metavar="NAME",
        help="read set used to pick each read's transcript (default: the first). "
        "Should be the UNCORRECTED reads: assigning from a corrected set lets that "
        "correction choose its own target, so a correction that pulled a read "
        "toward the wrong isoform would score better rather than worse.",
    )
    ap.add_argument(
        "--pair-against",
        default=None,
        metavar="NAME",
        help="baseline for the per-read paired comparison (default: the first "
        "set). Independent of --assign-from, so truth can stay neutral while the "
        "comparison is made against whichever implementation is the reference.",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=0,
        help="worker processes (0 = all cores). The best-of-N assignment is the "
        "cost here and is embarrassingly parallel per read.",
    )
    ap.add_argument(
        "--max-assign-error",
        type=float,
        default=0.0,
        help="drop reads whose best transcript match is worse than this percent "
        "(0 = keep all). Real data contains chimeras and contamination with no "
        "true SIRV of origin; the filter uses the *assignment* set, so it "
        "removes the same reads from every implementation.",
    )
    args = ap.parse_args()

    try:
        import edlib
    except ImportError:
        print("error: edlib missing. conda activate isoncorrect-ref", file=sys.stderr)
        return 1

    transcripts = read_fasta(args.transcriptome)
    print(f"transcriptome: {len(transcripts)} sequences")

    sets: dict[str, dict[str, str]] = {}
    for spec in args.reads:
        name, _, path = spec.partition("=")
        sets[name] = read_fastqs(path)
        print(f"  {name:<14} {len(sets[name]):>7} reads  ({path})")

    # Paired comparison: only reads every set has.
    common = set.intersection(*(set(v) for v in sets.values()))
    ordered = sorted(common)
    if args.max_reads:
        ordered = ordered[: args.max_reads]

    # Assign each read a transcript, once, from the first set, and score every
    # set against that assignment --- both in the same worker pass.
    assign_from = args.assign_from or next(iter(sets))
    if assign_from not in sets:
        print(f"error: --assign-from {assign_from!r} is not one of {list(sets)}", file=sys.stderr)
        return 1
    set_names = list(sets)
    workers = args.threads or (os.cpu_count() or 1)

    work = []
    from_header = 0
    for header in ordered:
        m = HEADER_TRANSCRIPT.search(header)
        forced = m.group(1) if m and m.group(1) in transcripts else None
        if forced is not None:
            from_header += 1
        # The assignment set must come first; `_score_chunk` reads seqs[0].
        seqs = [sets[assign_from][header]] + [
            sets[n][header] for n in set_names if n != assign_from
        ]
        work.append((header, forced, seqs))
    order = [assign_from] + [n for n in set_names if n != assign_from]

    chunk = 200
    chunks = [work[i : i + chunk] for i in range(0, len(work), chunk)]
    print(f"scoring on {workers} worker(s), {len(chunks)} chunks of <= {chunk}")

    records = []
    if workers > 1 and len(chunks) > 1:
        import multiprocessing as mp

        with mp.Pool(workers, initializer=_init, initargs=(transcripts,)) as pool:
            done = 0
            for got in pool.imap_unordered(_score_chunk, chunks):
                records.extend(got)
                done += 1
                if done % 50 == 0 or done == len(chunks):
                    print(f"\r  {done}/{len(chunks)} chunks", end="", flush=True)
        print()
    else:
        _init(transcripts)
        for c in chunks:
            records.extend(_score_chunk(c))

    if from_header:
        print(f"truth from read headers for {from_header} reads")
    searched = len(records) - from_header
    if searched:
        errs = sorted(r[2] for r in records if r[2] > 0)
        med = errs[len(errs) // 2] if errs else 0.0
        print(
            f"truth by best-of-{len(transcripts)} alignment for {searched} reads "
            f"(assigned from {assign_from!r}, median best-match error {med:.2f}%)"
        )
    if args.max_assign_error:
        before = len(records)
        records = [r for r in records if r[2] <= args.max_assign_error]
        print(
            f"dropped {before - len(records)} reads whose best match exceeded "
            f"{args.max_assign_error}% --- no plausible SIRV of origin"
        )
    print(f"scored {len(records)} reads present in every set\n")

    # Keyed by header, not appended positionally: `imap_unordered` returns chunks
    # in completion order, and a set missing a read would shift its list against
    # the others. The paired comparison below depends on this being right.
    per_read: dict[str, dict[str, float]] = {name: {} for name in sets}
    totals: dict[str, tuple[int, int]] = {name: (0, 0) for name in sets}
    for header, _truth, _aerr, errors in records:
        for name, cell in zip(order, errors):
            if cell is None:
                continue
            ed, length = cell
            per_read[name][header] = 100.0 * ed / length
            errs, bases = totals[name]
            totals[name] = (errs + ed, bases + length)
    results: dict[str, list[float]] = {
        name: [per_read[name][h] for h in sorted(per_read[name])] for name in sets
    }

    width = max(len(n) for n in sets)
    print(f"{'set':<{width}}  {'mean%':>7} {'median%':>8} {'p90%':>7} "
          f"{'total%':>7} {'perfect':>8} {'n':>7}")
    for name, vals in results.items():
        if not vals:
            print(f"{name:<{width}}  (no reads scored)")
            continue
        errs, bases = totals[name]
        perfect = sum(1 for v in vals if v == 0.0)
        vals_sorted = sorted(vals)
        print(
            f"{name:<{width}}  {statistics.mean(vals):>7.3f} "
            f"{statistics.median(vals):>8.3f} "
            f"{vals_sorted[int(0.9 * len(vals))]:>7.3f} "
            f"{100.0 * errs / bases:>7.3f} "
            f"{perfect:>8} {len(vals):>7}"
        )

    # Paired differences against the first set, which is the interesting number
    # when two implementations disagree on only a fraction of reads.
    names = list(sets)
    if len(names) > 1:
        base = args.pair_against or names[0]
        if base not in sets:
            print(f"error: --pair-against {base!r} is not one of {names}", file=sys.stderr)
            return 1
        print(f"\npaired against {base!r}, per read:")
        for name in [n for n in names if n != base]:
            shared = sorted(set(per_read[base]) & set(per_read[name]))
            deltas = [per_read[name][h] - per_read[base][h] for h in shared]
            better = sum(1 for d in deltas if d < 0)
            worse = sum(1 for d in deltas if d > 0)
            print(
                f"  {name:<{width}}  mean delta {statistics.mean(deltas):+.4f} pp  "
                f"better on {better}, worse on {worse}, equal on {len(deltas) - better - worse}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
