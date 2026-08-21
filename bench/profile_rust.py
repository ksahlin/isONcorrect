#!/usr/bin/env python3
"""Stage attribution for the Rust port, summed over several clusters.

The counterpart to `profile_python.py`, printing the same columns so the two can
be read side by side — which is the only way to answer "which stage is the port
slower at than the reference", as opposed to "where does the port spend time".

Runs the binary once per cluster with `ISONCORRECT_PROFILE=1` and sums the
per-cluster tables. Single-threaded on purpose: one cluster per process keeps the
stage accumulators clean, and the totals are what matter, not wall clock.

    bench/profile_rust.py --clusters /tmp/real_trunc --n 10

`--binary` defaults to the release build. `--exact-guard` profiles the exact
parasail-compatible DP instead of block-aligner, which is how the two are
compared.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

ROW = re.compile(r"^(.+?)\s{2,}(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$")
UNATTRIB = re.compile(r"^unattributed.*?([\d.]+)\s+([\d.]+)\s*$")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--clusters", required=True, help="folder of per-cluster fastqs")
    ap.add_argument("--n", type=int, default=10, help="how many clusters")
    ap.add_argument(
        "--binary", default="rust/target/release/isONcorrect", help="path to isONcorrect"
    )
    ap.add_argument("--outdir", default="/tmp/profile_rust")
    ap.add_argument(
        "--exact-guard",
        action="store_true",
        help="profile the exact parasail-compatible guard instead of block-aligner",
    )
    ap.add_argument("--args", default="", help="extra flags forwarded to isONcorrect")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.clusters, "*.fastq")))[: args.n]
    if not files:
        print(f"error: no fastq files in {args.clusters}", file=sys.stderr)
        return 1

    env = dict(os.environ, ISONCORRECT_PROFILE="1")
    if args.exact_guard:
        env["ISONCORRECT_EXACT_GUARD"] = "1"

    calls: dict[str, int] = defaultdict(int)
    self_s: dict[str, float] = defaultdict(float)
    incl_s: dict[str, float] = defaultdict(float)
    unattributed = 0.0
    reads = 0

    for path in files:
        name = os.path.basename(path)[: -len(".fastq")]
        with open(path) as fh:
            reads += sum(1 for _ in fh) // 4
        out = subprocess.run(
            [args.binary, "--fastq", path, "--outfolder", os.path.join(args.outdir, name)]
            + args.args.split(),
            env=env,
            capture_output=True,
            text=True,
        )
        if out.returncode != 0:
            print(f"error: {args.binary} failed on {path}", file=sys.stderr)
            print(out.stderr[-2000:], file=sys.stderr)
            return 1
        for line in out.stdout.splitlines():
            if line.startswith(("stage", "-", "self column", "accounted")):
                continue
            m = UNATTRIB.match(line.strip())
            if m:
                unattributed += float(m.group(1))
                continue
            m = ROW.match(line.rstrip())
            if not m:
                continue
            label = m.group(1).strip()
            calls[label] += int(m.group(2))
            self_s[label] += float(m.group(3))
            incl_s[label] += float(m.group(4))

    guard = "exact DP" if args.exact_guard else "block-aligner"
    total = sum(self_s.values()) + unattributed
    print(f"\n{len(files)} clusters, {reads} reads, guard = {guard}")
    print("self column sums to the total. 'inclusive' double-counts nesting.\n")
    header = f"{'stage':<40}{'calls':>10}{'self_s':>10}{'self_%':>8}{'incl_s':>10}"
    print(header)
    print("-" * len(header))
    for label in sorted(self_s, key=lambda k: -self_s[k]):
        pct = 100.0 * self_s[label] / total if total else 0.0
        print(
            f"{label:<40}{calls[label]:>10}{self_s[label]:>10.2f}"
            f"{pct:>8.1f}{incl_s[label]:>10.2f}"
        )
    print("-" * len(header))
    pct = 100.0 * unattributed / total if total else 0.0
    print(f"{'unattributed (I/O, parsing, glue)':<40}{'':>10}{unattributed:>10.2f}{pct:>8.1f}")
    print(f"{'total':<40}{'':>10}{total:>10.2f}{100.0:>8.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
