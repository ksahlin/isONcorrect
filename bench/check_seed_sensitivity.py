#!/usr/bin/env python3
"""Decide whether a `get_best_corrections` mismatch is a port bug or the
reference disagreeing with itself.

`get_alternative_ref_contexts` returns a **`set`** per column, and
`get_best_corrections` iterates it with a `break` on an exact context match and
a strict `<` on edit distance. Set iteration order over tuples of strings
depends on `PYTHONHASHSEED`, so on a column with two or more alternatives that
tie, the reference's answer depends on the hash seed — there is no single
correct output to port to.

This replays the recorded cases the Rust oracle flagged, once per seed, and
reports for each whether the reference is stable and whether the port's answer
is one the reference produces:

    bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
    CORRECTION_CASES=/tmp/d/correction_cases.tsv \
      CORRECTION_MISMATCHES=/tmp/d/mismatches.tsv \
      cargo test --manifest-path rust/Cargo.toml corrections::oracle
    bench/check_seed_sensitivity.py --cases /tmp/d/correction_cases.tsv \
      --mismatches /tmp/d/mismatches.tsv

Each seed needs its own interpreter: `PYTHONHASHSEED` is read at startup, so
setting it from inside a running process does nothing. That is worth knowing —
it is exactly the mistake that made this non-determinism invisible for so long.
"""

from __future__ import annotations

import argparse
import array
import os
import subprocess
import sys
import tempfile

REPLAY = r"""
import sys, os, array, json
sys.argv = ["x"]
from isoncorrect import isONcorrect as ion

case_path, wanted = sys.argv[1] if False else os.environ["CASE_FILE"], os.environ["CASE_IDS"]
wanted = set(wanted.split(","))
cases = {}
for line in open(case_path):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[1] not in wanted:
        continue
    cases.setdefault(f[1], {"S": []})
    if f[0] == "C":
        cases[f[1]]["C"] = f
    elif f[0] == "S":
        cases[f[1]]["S"].append(f)

out = {}
workdir = os.environ["WORKDIR"]
for cid, case in cases.items():
    C = case["C"]
    k, T, ms = int(C[2]), float(C[3]), int(C[4])
    reads, flat = {}, array.array("I")
    for _, _, q, p1, p2, seg in case["S"]:
        q, p1, p2 = int(q), int(p1), int(p2)
        reads[q] = ("acc", "A" * p1 + seg, "I" * (p1 + len(seg)))
        flat.extend([q, p1, p2])
    curr, others = ion.get_best_corrections(flat, reads, k, workdir, T, ms, False, False)
    rows = [f"C\t{cid}\t{curr}"]
    for q_id, regions in others.items():
        for start, stop, weight, corrected in regions:
            rows.append(f"O\t{cid}\t{q_id}\t{start}\t{stop}\t{weight}\t{corrected}")
    out[cid] = rows
print("---JSON---")
print(json.dumps(out))
"""


def run_seed(seed: int, case_file: str, case_ids: list[str], workdir: str) -> dict:
    import json

    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(REPLAY)
        script = fh.name
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = str(seed)
    env["CASE_FILE"] = case_file
    env["CASE_IDS"] = ",".join(case_ids)
    env["WORKDIR"] = workdir
    proc = subprocess.run(
        [sys.executable, script], capture_output=True, text=True, env=env
    )
    os.unlink(script)
    if "---JSON---" not in proc.stdout:
        print(proc.stdout[-2000:], file=sys.stderr)
        print(proc.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"replay failed at seed {seed}")
    return json.loads(proc.stdout.split("---JSON---", 1)[1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cases", required=True, help="correction_cases.tsv")
    ap.add_argument("--mismatches", required=True, help="the Rust oracle's output")
    ap.add_argument("--seeds", type=int, default=8)
    args = ap.parse_args()

    port: dict[str, list[str]] = {}
    for line in open(args.mismatches):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        port.setdefault(f[1], []).append(line.rstrip("\n"))
    if not port:
        print("no mismatches to check")
        return 0

    case_ids = sorted(port, key=int)
    print(f"checking {len(case_ids)} mismatching case(s) across {args.seeds} hash seeds")

    workdir = tempfile.mkdtemp(prefix="seedcheck")
    by_seed = {s: run_seed(s, args.cases, case_ids, workdir) for s in range(args.seeds)}

    # Compare record by record, not answer by answer: one seed-sensitive record
    # would otherwise mask every stable one beside it, and a stable record the
    # port gets wrong is the only thing here that is actually a bug.
    unstable, stable_ok, stable_bad = 0, 0, 0
    for cid in case_ids:
        per_record: dict[str, set[str]] = {}
        for result in by_seed.values():
            for line in result.get(cid, []):
                key, value = _split(line)
                per_record.setdefault(key, set()).add(value)

        case_unstable, case_bad = [], []
        for line in port[cid]:
            key, value = _split(line)
            seen = per_record.get(key)
            if seen is None:
                case_bad.append((key, value, "<the reference emits no such record>"))
            elif len(seen) > 1:
                case_unstable.append((key, value, value in seen))
            elif value not in seen:
                case_bad.append((key, value, next(iter(seen))))

        unstable += len(case_unstable)
        stable_bad += len(case_bad)
        stable_ok += len(port[cid]) - len(case_unstable) - len(case_bad)

        print(
            f"  case {cid}: {len(port[cid])} record(s) --- "
            f"{len(case_unstable)} seed-dependent, {len(case_bad)} stable but wrong"
        )
        for key, value, agrees in case_unstable:
            mark = "port is one of them" if agrees else "PORT NOT AMONG THEM"
            print(f"      {key}: {len(per_record[key])} reference answers, {mark}")
        for key, value, want in case_bad:
            print(f"      {key}: STABLE MISMATCH\n        reference: {want}\n        port     : {value}")

    print()
    print(f"{stable_ok} record(s) agree, {unstable} seed-dependent, {stable_bad} stable mismatches")
    # A stable reference the port disagrees with is a real bug; the rest is the
    # reference's own non-determinism.
    return 1 if stable_bad else 0


def _split(line: str) -> tuple[str, str]:
    """`O <case> <q_id> <start> <stop> <weight> <seq>` keyed by everything but
    the sequence, so the same region compares across seeds."""
    f = line.split("\t")
    return "\t".join(f[:-1]), f[-1]


if __name__ == "__main__":
    sys.exit(main())
