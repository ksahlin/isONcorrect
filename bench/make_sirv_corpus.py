#!/usr/bin/env python3
"""Split a simulated SIRV fastq into per-cluster fastqs, by read header.

The `test_data/isoncorrect/` fixtures are one 100-read cluster wearing two
filenames (see PORTING.md), which is a smoke test and nothing else. This builds
a corpus with the properties that fixture lacks: clusters large enough to leave
the `--exact` path, isoform variation inside a cluster, and enough reads to
cross `--max_seqs` and exercise batching.

The reads are simulated with the source transcript in the header:

    @read_1_from_SIRV612

so grouping needs no aligner and no clustering tool --- and, unlike an
isONclust run, it is reproducible from the fastq alone, which is what a
verification corpus wants.

Two groupings, and they stress different things:

    --group gene        SIRV6, SIRV5, ...   isoforms of one gene share a
                        cluster, so reads differ by whole exons. This is the
                        case the structural-overcorrection guard exists for.
    --group transcript  SIRV612, SIRV601... one transcript per cluster: high
                        coverage, differences are sequencing error only.

Output is isONclust's layout, so `run_isoncorrect --fastq_folder` reads it
directly:

    <outdir>/0.fastq, 1.fastq, ...   largest cluster first, matching isONclust
    <outdir>/clusters.tsv            cluster id, group name, reads

Clusters below --min-reads are dropped: isONcorrect needs three reads carrying
an anchor before it will correct anything.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

# "read_1_from_SIRV612" -> "SIRV612", and SIRV6 at gene level.
HEADER = re.compile(r"_from_(?P<transcript>\S+)")
GENE = re.compile(r"^(?P<gene>[A-Za-z]+\d)")


def group_of(header: str, level: str) -> str | None:
    m = HEADER.search(header)
    if not m:
        return None
    transcript = m.group("transcript")
    if level == "transcript":
        return transcript
    g = GENE.match(transcript)
    return g.group("gene") if g else transcript


def read_fastq(path):
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
            if not qual:
                print(f"error: truncated record at end of {path}", file=sys.stderr)
                return
            yield header.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--group", choices=("gene", "transcript"), default="gene")
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="truncate each cluster to this many reads (0 = keep all)",
    )
    args = ap.parse_args()

    clusters: dict[str, list] = {}
    skipped = 0
    for header, seq, qual in read_fastq(args.fastq):
        name = group_of(header, args.group)
        if name is None:
            skipped += 1
            continue
        clusters.setdefault(name, []).append((header, seq, qual))

    if skipped:
        print(f"warning: {skipped} reads had no _from_ tag in the header", file=sys.stderr)

    # isONclust numbers clusters largest first; ties by name so a rerun on the
    # same input always produces the same file names.
    ordered = sorted(clusters.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    ordered = [(n, r) for n, r in ordered if len(r) >= args.min_reads]

    os.makedirs(args.outdir, exist_ok=True)
    sizes = []
    with open(os.path.join(args.outdir, "clusters.tsv"), "w") as index:
        index.write("#cluster\tgroup\treads\n")
        for cluster_id, (name, reads) in enumerate(ordered):
            if args.max_reads:
                reads = reads[: args.max_reads]
            with open(os.path.join(args.outdir, f"{cluster_id}.fastq"), "w") as fh:
                for header, seq, qual in reads:
                    fh.write(f"{header}\n{seq}\n+\n{qual}\n")
            index.write(f"{cluster_id}\t{name}\t{len(reads)}\n")
            sizes.append(len(reads))

    total = sum(sizes)
    print(f"==> {len(ordered)} clusters, {total} reads -> {args.outdir}")
    if sizes:
        print(f"    sizes: min {min(sizes)}, median {sorted(sizes)[len(sizes) // 2]}, max {max(sizes)}")
        print(f"    above the --exact_instance_limit of 50: {sum(1 for s in sizes if s > 50)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
