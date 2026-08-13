#!/usr/bin/env python3
"""Record real (x, y, k, edit_distance) cases from the reference's edlib.

The Rust port computes this distance natively rather than binding edlib, which
is safe because edit distance is a uniquely defined number (see
`rust/src/editdist.rs`). "Safe in principle" still needs checking, so this
samples real substrings from a cluster and records what edlib actually returns.

    conda activate isoncorrect-ref
    bench/gen_editdist_cases.py --fastq test_data/isoncorrect/0.fastq --out /tmp/ed.tsv
    EDLIB_CASES=/tmp/ed.tsv cargo test --manifest-path rust/Cargo.toml editdist -- --nocapture
"""
import argparse, random, sys

ap = argparse.ArgumentParser()
ap.add_argument("--fastq", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--n", type=int, default=4000)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

try:
    from isoncorrect.help_functions import readfq
    from isoncorrect.isONcorrect import edlib_alignment
except ImportError as exc:
    sys.exit(f"error: {exc}\nRun inside the reference env: conda activate isoncorrect-ref")

reads = [s for _, (s, _) in readfq(open(args.fastq))]
if not reads:
    sys.exit("no reads")

random.seed(args.seed)
written = 0
with open(args.out, "w") as fh:
    fh.write("#x\ty\tk\ted\n")
    while written < args.n:
        a, b = random.choice(reads), random.choice(reads)
        i = random.randrange(0, max(1, len(a) - 60))
        j = random.randrange(0, max(1, len(b) - 60))
        x = a[i : i + random.randrange(10, 90)]
        y = b[j : j + random.randrange(10, 90)]
        if not x or not y:
            continue
        k = len(x)
        fh.write(f"{x}\t{y}\t{k}\t{edlib_alignment(x, y, k)}\n")
        written += 1
print(f"wrote {written} cases to {args.out}")
