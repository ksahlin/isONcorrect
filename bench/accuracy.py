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

**The assignment is made once, from the first read set, and reused for all of
them.** That matters: if each implementation picked its own best-matching
transcript, every one would be flattered by construction, and a correction that
dragged a read toward the wrong isoform would score *better* rather than worse.
Pass the uncorrected reads first so the assignment is independent of any
correction.

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

    # Assign each read a transcript, once, from the first set.
    assign_from = next(iter(sets))
    names_sorted = sorted(transcripts)
    assigned: dict[str, str] = {}
    assign_err: dict[str, float] = {}
    from_header = 0
    for header in ordered:
        m = HEADER_TRANSCRIPT.search(header)
        if m and m.group(1) in transcripts:
            assigned[header] = m.group(1)
            assign_err[header] = 0.0
            from_header += 1
            continue
        seq = sets[assign_from][header]
        if not seq:
            continue
        best_name, best_ed = None, None
        for tname in names_sorted:
            ed = edlib.align(seq, transcripts[tname], mode="HW", task="distance")[
                "editDistance"
            ]
            if best_ed is None or ed < best_ed:
                best_name, best_ed = tname, ed
        if best_name is not None:
            assigned[header] = best_name
            assign_err[header] = 100.0 * best_ed / len(seq)

    if from_header:
        print(f"truth from read headers for {from_header} reads")
    searched = len(assigned) - from_header
    if searched:
        errs = sorted(assign_err[h] for h in assigned if assign_err[h] > 0)
        med = errs[len(errs) // 2] if errs else 0.0
        print(
            f"truth by best-of-{len(transcripts)} alignment for {searched} reads "
            f"(assigned from {assign_from!r}, median best-match error {med:.2f}%)"
        )
    if args.max_assign_error:
        before = len(assigned)
        assigned = {
            h: t for h, t in assigned.items() if assign_err[h] <= args.max_assign_error
        }
        print(
            f"dropped {before - len(assigned)} reads whose best match exceeded "
            f"{args.max_assign_error}% --- no plausible SIRV of origin"
        )
    ordered = [h for h in ordered if h in assigned]
    print(f"scoring {len(ordered)} reads present in every set\n")

    results: dict[str, list[float]] = {name: [] for name in sets}
    totals: dict[str, tuple[int, int]] = {name: (0, 0) for name in sets}
    for header in ordered:
        truth = transcripts[assigned[header]]
        for name, reads in sets.items():
            seq = reads[header]
            if not seq:
                continue
            # "HW" = infix: the read is aligned in full, the transcript's ends are
            # free. That is the right model for a read drawn from a transcript.
            ed = edlib.align(seq, truth, mode="HW", task="distance")["editDistance"]
            results[name].append(100.0 * ed / len(seq))
            errs, bases = totals[name]
            totals[name] = (errs + ed, bases + len(seq))

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
        base = names[0]
        print(f"\npaired against {base!r}, per read:")
        for name in names[1:]:
            deltas = [b - a for a, b in zip(results[base], results[name])]
            better = sum(1 for d in deltas if d < 0)
            worse = sum(1 for d in deltas if d > 0)
            print(
                f"  {name:<{width}}  mean delta {statistics.mean(deltas):+.4f} pp  "
                f"better on {better}, worse on {worse}, equal on {len(deltas) - better - worse}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
