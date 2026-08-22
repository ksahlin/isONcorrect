#!/usr/bin/env python3
"""Per-read accuracy against a genome, via spliced alignment.

The companion to `accuracy.py`. That one needs a transcriptome and scores each
read against one transcript with edlib; this one needs only a **genome** and lets
`minimap2 -ax splice` place the read across exons. Use it when there is no
trustworthy transcriptome for the organism — real ONT cDNA from a species whose
annotation you do not want to depend on.

    bench/accuracy_genome.py --genome fruitfly.fa \\
        --reads raw=/tmp/droso_clusters \\
        --reads python=/tmp/droso_py \\
        --reads rust=/tmp/droso_rs

# What is measured

Per read, from the primary alignment's CIGAR and `NM` tag:

* **error rate** = `NM / aligned_read_bases`, where `NM` is edit distance to the
  reference over the aligned block. This is the headline number.
* **aligned fraction** = `aligned_read_bases / read_length`, i.e. how much of the
  read minimap2 could place at all. Soft clips do not count as aligned.

**Both are reported, and the second is not optional.** Error rate alone is
trivially gamed by aligning less of the read: a correction that mangled a read's
ends would soft-clip them away and *improve* its error rate. The pair together
cannot be gamed in that direction.

`NM` over a genome includes real biological difference — SNPs against the
reference assembly, unannotated splice sites minimap2 mis-places — so the
absolute rate is an overestimate of sequencing error. That inflation is identical
across implementations, so comparisons hold even where the absolute number does
not. Do not quote the absolute figure as "the error rate" without that caveat.

# Why per-set alignment is fair here, unlike per-set assignment there

`accuracy.py` has to assign each read to one of many similar transcripts, and
letting each implementation pick its own would flatter every one of them — a
correction dragging a read toward the wrong isoform would score *better*. Here
each set is aligned independently, which is fine for a different reason: a genome
locus is essentially unambiguous, so there is no menu of near-identical targets
to be flattered by. What *is* still a risk is reads dropping out, so:

* only reads with a primary alignment **in every set** are scored, making the
  comparison paired;
* the number each set failed to align is reported, because a set that silently
  lost reads would otherwise look good.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def read_fastqs(path: str) -> dict[str, str]:
    """Header -> sequence, over a fastq or a folder of them.

    Same layout rules as `accuracy.py`, so the two take identical --reads paths.
    """
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


def write_fasta(seqs: dict[str, str], path: str) -> None:
    """One record per read, named by an index so minimap2 never splits a header.

    Read names here contain spaces and pipes; SAM's QNAME stops at the first
    whitespace, so the index is the join key and the mapping is kept in memory.
    """
    with open(path, "w") as fh:
        for i, (header, seq) in enumerate(seqs.items()):
            if seq:
                fh.write(f">r{i}\n{seq}\n")


def align(genome: str, fasta: str, threads: int) -> str:
    """Spliced alignment, primary hits only, as SAM on a temp path."""
    sam = fasta + ".sam"
    cmd = [
        "minimap2",
        "-ax",
        "splice",
        "-uf",  # cDNA is stranded: only consider the transcript strand
        "--secondary=no",
        "-t",
        str(threads),
        genome,
        fasta,
    ]
    with open(sam, "w") as out, open(os.devnull, "w") as null:
        rc = subprocess.call(cmd, stdout=out, stderr=null)
    if rc != 0:
        raise SystemExit(f"minimap2 failed ({rc}) on {fasta}")
    return sam


def score_sam(sam: str) -> dict[str, tuple[int, int, int]]:
    """`r<i>` -> (NM, aligned read bases, read length).

    Aligned read bases count `M`/`=`/`X`/`I` — everything the read contributes
    inside the alignment. `S`/`H` are the clipped remainder and `N` is the intron,
    neither of which the read spends bases on within the block.
    """
    out: dict[str, tuple[int, int, int]] = {}
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            name, flag, cigar = f[0], int(f[1]), f[5]
            # 0x4 unmapped, 0x100 secondary, 0x800 supplementary. Supplementary
            # rows are dropped rather than merged: a chimeric read has no single
            # error rate, and dropping it affects every set alike.
            if flag & 0x904 or cigar == "*":
                continue
            nm = None
            for tag in f[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break
            if nm is None:
                continue
            aligned = 0
            read_len = 0
            for num, op in CIGAR_OP.findall(cigar):
                n = int(num)
                if op in "M=XI":
                    aligned += n
                    read_len += n
                elif op in "SH":
                    read_len += n
            if aligned == 0:
                continue
            # First primary wins; with --secondary=no there should only be one.
            out.setdefault(name, (nm, aligned, read_len))
    return out


def pct(values: list[float]) -> tuple[float, float, float]:
    ordered = sorted(values)
    return (
        statistics.fmean(ordered),
        statistics.median(ordered),
        ordered[min(len(ordered) - 1, int(len(ordered) * 0.9))],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genome", required=True)
    ap.add_argument(
        "--reads",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="a labelled read set; repeat for each implementation",
    )
    ap.add_argument(
        "--pair-against",
        help="baseline for the per-read paired comparison (default: the first set)",
    )
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument(
        "--keep",
        help="directory to keep the SAM files in, instead of a temp dir that is deleted",
    )
    args = ap.parse_args()

    if not shutil.which("minimap2"):
        print("error: minimap2 not on PATH", file=sys.stderr)
        return 1

    sets: list[tuple[str, str]] = []
    for spec in args.reads:
        if "=" not in spec:
            print(f"error: --reads wants NAME=PATH, got {spec!r}", file=sys.stderr)
            return 1
        name, path = spec.split("=", 1)
        sets.append((name, path))

    workdir = args.keep or tempfile.mkdtemp(prefix="accgenome_")
    os.makedirs(workdir, exist_ok=True)

    # Read every set first, so the index -> header mapping is shared. Headers
    # survive correction unchanged, which is what makes the sets comparable.
    loaded: dict[str, dict[str, str]] = {}
    for name, path in sets:
        loaded[name] = read_fastqs(path)
        print(f"  {name:<16}{len(loaded[name]):>7} reads  ({path})")

    common = set.intersection(*(set(d) for d in loaded.values()))
    order = [h for h in loaded[sets[0][0]] if h in common]
    print(f"{len(order)} reads present in every set")

    scored: dict[str, dict[int, tuple[int, int, int]]] = {}
    for name, _ in sets:
        seqs = {h: loaded[name][h] for h in order}
        fasta = os.path.join(workdir, f"{name}.fa")
        write_fasta(seqs, fasta)
        sam = align(args.genome, fasta, args.threads)
        raw = score_sam(sam)
        scored[name] = {int(k[1:]): v for k, v in raw.items()}
        missing = len(order) - len(scored[name])
        print(f"  {name:<16}aligned {len(scored[name]):>7}, unaligned {missing}")

    # Paired: only indices with a primary alignment everywhere.
    usable = sorted(set.intersection(*(set(d) for d in scored.values())))
    print(f"\n{len(usable)} reads aligned in every set (paired comparison)\n")
    if not usable:
        return 1

    header = f"{'set':<16}{'err_mean%':>11}{'err_med%':>10}{'err_p90%':>10}"
    header += f"{'err_total%':>12}{'aln_mean%':>11}{'aln_med%':>10}"
    print(header)
    print("-" * len(header))
    per_set_err: dict[str, dict[int, float]] = {}
    for name, _ in sets:
        errs, alns = [], []
        nm_sum, aln_sum, len_sum = 0, 0, 0
        per_read: dict[int, float] = {}
        for i in usable:
            nm, aligned, read_len = scored[name][i]
            e = 100.0 * nm / aligned
            errs.append(e)
            per_read[i] = e
            alns.append(100.0 * aligned / read_len if read_len else 0.0)
            nm_sum += nm
            aln_sum += aligned
            len_sum += read_len
        per_set_err[name] = per_read
        em, emed, ep90 = pct(errs)
        am, amed, _ = pct(alns)
        total = 100.0 * nm_sum / aln_sum if aln_sum else 0.0
        print(
            f"{name:<16}{em:>11.3f}{emed:>10.3f}{ep90:>10.3f}"
            f"{total:>12.3f}{am:>11.3f}{amed:>10.3f}"
        )

    base = args.pair_against or sets[0][0]
    if base not in per_set_err:
        print(f"error: --pair-against {base!r} is not one of the sets", file=sys.stderr)
        return 1
    print(f"\npaired against {base!r}, per read:")
    for name, _ in sets:
        if name == base:
            continue
        deltas = [per_set_err[name][i] - per_set_err[base][i] for i in usable]
        better = sum(1 for d in deltas if d < 0)
        worse = sum(1 for d in deltas if d > 0)
        print(
            f"  {name:<16}mean delta {statistics.fmean(deltas):+.4f} pp  "
            f"better on {better}, worse on {worse}, equal on {len(deltas) - better - worse}"
        )

    if not args.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
