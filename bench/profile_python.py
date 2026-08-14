#!/usr/bin/env python3
"""Stage-level profiler for the isONcorrect Python reference.

Where the time goes has never been measured, and the port should be aimed by
measurement rather than by intuition. This wraps the interesting functions and
reports both inclusive and exclusive (self) wall time, so nested costs -- e.g.
the spoa subprocess inside get_best_corrections -- are attributed to the right
place instead of being counted twice.

Must run inside the pinned reference environment:

    conda activate isoncorrect-ref
    bench/profile_python.py --fastq test_data/isoncorrect/0.fastq --outfolder /tmp/prof

Extra isONcorrect flags pass straight through:

    bench/profile_python.py --fastq cluster.fastq --outfolder /tmp/prof -- --k 9 --w 10

Add --cprofile to also dump a .prof for snakeviz/pstats.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import defaultdict

# ---------------------------------------------------------------------------
# A small nesting-aware timer.
#
# cProfile alone answers "which function", but it does not separate the spoa
# fork/exec from the Python around it in a way that reads naturally, and its
# overhead distorts the many-small-calls stages. This is coarse and cheap.
# ---------------------------------------------------------------------------

inclusive: dict[str, float] = defaultdict(float)
exclusive: dict[str, float] = defaultdict(float)
counts: dict[str, int] = defaultdict(int)
_stack: list[list[float]] = []  # each entry: [start_time, accumulated_child_time]


def instrument(module, name, label=None):
    """Wrap module.name with a timer. Returns silently if the attr is missing."""
    label = label or f"{module.__name__.split('.')[-1]}.{name}"
    original = getattr(module, name, None)
    if original is None:
        print(f"warning: {label} not found, skipping", file=sys.stderr)
        return

    def wrapper(*args, **kwargs):
        _stack.append([time.perf_counter(), 0.0])
        try:
            return original(*args, **kwargs)
        finally:
            start, child = _stack.pop()
            elapsed = time.perf_counter() - start
            inclusive[label] += elapsed
            exclusive[label] += elapsed - child
            counts[label] += 1
            if _stack:
                _stack[-1][1] += elapsed

    wrapper.__name__ = getattr(original, "__name__", name)
    wrapper.__doc__ = getattr(original, "__doc__", None)
    setattr(module, name, wrapper)


def report(total_wall: float) -> None:
    print()
    print("=" * 78)
    print("STAGE ATTRIBUTION")
    print("=" * 78)
    print(f"total wall time: {total_wall:.2f}s")
    print()
    print("'self' excludes time spent inside other instrumented functions, so the")
    print("self column sums to at most the total. 'inclusive' double-counts nesting.")
    print()
    header = f"{'stage':<46}{'calls':>10}{'self_s':>10}{'self_%':>8}{'incl_s':>10}"
    print(header)
    print("-" * len(header))
    for label in sorted(exclusive, key=lambda k: -exclusive[k]):
        self_s = exclusive[label]
        pct = 100.0 * self_s / total_wall if total_wall else 0.0
        print(f"{label:<46}{counts[label]:>10}{self_s:>10.2f}{pct:>8.1f}{inclusive[label]:>10.2f}")

    accounted = sum(exclusive.values())
    print("-" * len(header))
    print(f"{'accounted for':<46}{'':>10}{accounted:>10.2f}{100.0 * accounted / total_wall if total_wall else 0:>8.1f}")
    print(f"{'unattributed (I/O, parsing, interpreter, glue)':<46}{'':>10}{total_wall - accounted:>10.2f}"
          f"{100.0 * (total_wall - accounted) / total_wall if total_wall else 0:>8.1f}")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Stage-level profiler for the isONcorrect Python reference",
        epilog="Flags after `--` are forwarded verbatim to isONcorrect.",
    )
    ap.add_argument("--fastq", required=True, help="input cluster fastq")
    ap.add_argument("--outfolder", required=True, help="output folder")
    ap.add_argument("--cprofile", action="store_true", help="also write a cProfile .prof dump")
    ap.add_argument("--prof-out", default=None, help="path for the .prof dump")
    ap.add_argument("passthrough", nargs="*", help="extra isONcorrect flags")
    args = ap.parse_args()

    # No hash-seed pin here, deliberately. Setting PYTHONHASHSEED from inside a
    # running interpreter does nothing --- CPython reads it at startup --- so the
    # line that used to sit here was decorative. Timings do not depend on it;
    # `dump_reference.py` re-execs itself because its *output* does. See
    # PORTING.md.

    try:
        import edlib
        from isoncorrect import correct_seqs, create_augmented_reference
        from isoncorrect import isONcorrect as ion
    except ImportError as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("Run inside the reference env: conda activate isoncorrect-ref", file=sys.stderr)
        print("Create it with: bench/setup_reference_env.sh", file=sys.stderr)
        return 1

    # Anchor/indexing stages
    instrument(ion, "get_minimizers_and_positions")
    instrument(ion, "get_minimizer_combinations_database")
    instrument(ion, "get_qvs")
    # Per-read search and scheduling
    instrument(ion, "find_most_supported_span")
    instrument(ion, "solve_WIS")
    instrument(ion, "get_intervals_to_correct")
    # Correction
    instrument(ion, "correct_read")
    instrument(ion, "get_best_corrections")
    instrument(ion, "get_contexts")
    instrument(ion, "get_alternative_ref_contexts")
    # The structural-overcorrection guard at the end of correct_read
    instrument(ion, "parasail_alignment")
    instrument(ion, "fix_correction")
    # Consensus and MSA
    instrument(create_augmented_reference, "run_spoa", label="spoa subprocess")
    instrument(correct_seqs, "create_multialignment_matrix")
    # edlib is imported by both modules but patching the module attribute catches
    # every call site, since they all go through `edlib.align`.
    instrument(edlib, "align", label="edlib.align")

    argv = ["isONcorrect", "--fastq", args.fastq, "--outfolder", args.outfolder]
    argv += [a for a in args.passthrough if a != "--"]
    sys.argv = argv

    print(f"==> profiling: {' '.join(argv)}")

    profiler = None
    if args.cprofile:
        import cProfile

        profiler = cProfile.Profile()
        profiler.enable()

    start = time.perf_counter()
    try:
        ion.main()
    except SystemExit as exc:
        if exc.code not in (0, None):
            print(f"isONcorrect exited with {exc.code}", file=sys.stderr)
            return int(exc.code)
    total = time.perf_counter() - start

    if profiler is not None:
        profiler.disable()
        out = args.prof_out or os.path.join(args.outfolder, "isoncorrect.prof")
        profiler.dump_stats(out)
        print(f"==> cProfile dump: {out}")
        print(f"    inspect with: python -m pstats {out}")

    report(total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
