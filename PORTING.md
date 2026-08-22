# isONcorrect — Rust rewrite

## Goal

Port isONcorrect from Python to Rust. **Identical CLI**, primarily faster and lower memory. The
Python implementation in `src/isoncorrect/` is the **normative reference**: when Rust and Python
disagree, Python is right until a human decides otherwise.

**Output is byte-identical except at two alignment call sites, both of which now use affine SIMD
alignment rather than reproducing the reference's aligner.** Both are deliberate, measured
divergences taken for large speedups on the port's two dominant costs; each has an env var that
restores the exact reference-compatible path, `bench/equivalence.sh` sets both, and the 29-case
gate still covers the whole pipeline.

| Call site | Reference uses | Port defaults to | Restore with |
| --- | --- | --- | --- |
| the structural-overcorrection guard | parasail, affine semi-global | `block-aligner`, same scoring, own tie-break — changes ~0.8% of reads | `ISONCORRECT_EXACT_GUARD=1` |
| segment against consensus, which builds the MSA | edlib, **unit-cost** global | `block-aligner`, **affine** (`match 4, mismatch -8, open 12, extend 1`) — changes ~11% of alignments | `ISONCORRECT_EXACT_ALIGN=1` |

See *The guard's aligner: a deliberate divergence* and *The MSA's aligner: affine instead of edit
distance*. The second is the deeper of the two: that CIGAR builds the MSA, which builds the PFM,
which decides every corrected base.

Two binaries must keep their exact current names, flags, and defaults:

- `isONcorrect` — corrects one cluster (a single fastq), writes `<outfolder>/corrected_reads.fastq`
- `run_isoncorrect` — fans out `isONcorrect` over a folder of per-cluster fastqs (isONclust output)

## Layout

| Path | Role |
| --- | --- |
| `src/isoncorrect/isONcorrect.py` | Reference: minimizers, anchor spans, WIS, correction driver |
| `src/isoncorrect/correct_seqs.py` | Reference: MSA matrix, PFM, consensus correction |
| `src/isoncorrect/create_augmented_reference.py` | Reference: spoa/racon subprocess wrappers |
| `src/isoncorrect/run_isoncorrect.py` | Reference: batch splitting + multiprocessing driver |
| `rust/` | The port. `src/cli.rs` CLI parity, `src/params.rs` resolved parameters, `src/validate.rs` argument validation, `src/bin/` the two binaries |
| `bench/` | Equivalence + performance harness |
| `test_data/isoncorrect/{0,1}.fastq` | The CI / harness fixture. **These two files are byte-identical** — see below |
| `tools/repo-slim/` | Staged tooling to strip committed data from history |
| `paper/` | The snakemake evaluations behind the paper (`evaluation*`, `ont_error_rates`). Not part of the port. Do not modify. |
| `scripts/` | Paper experiment scripts. Not part of the port. `scripts/correction_pipeline.sh` is linked from README and must keep its path. |

## Port status

| Stage | State |
| --- | --- |
| CLI parity (both binaries) | done, locked by unit tests |
| fastq reading + `r_id` order | done, byte-identical dumps |
| minimizers (`lex`) | done, byte-identical dumps across 6 parameter settings |
| minimizer-pair database | done, byte-identical dumps across 6 parameter settings; stores 1 459 keys where the reference holds 23 260 |
| `solve_WIS` + `fill_p2` | done; replayed against every recorded reference call |
| bounded edit distance (edlib `'dist'`) | done natively, no C++; 4 000 real pairs verified against Python edlib |
| quality-value prefix sums (`get_qvs`) | done; `D` table cross-checked against the reference |
| `find_most_supported_span` | done; 52 558 intervals byte-identical to the reference |
| anchor filtering (`previously_corrected_regions`, `pos_group`) | done; live on the 100-read fixture and verified through end-to-end output |
| correction driver (`correct_read`) | done, see below |
| consensus POA (`run_spoa`) | done via `spoars`; oracle green on 5 922 real intervals |
| NW CIGAR (edlib `task="path"`) | done natively, no C++; 187 744 real alignments byte-identical |
| `cigar_to_seq` | done; 186 330 real expansions byte-identical |
| MSA matrix (`create_multialignment_matrix`) | done; 192 252 rows across 5 922 real matrices byte-identical |
| PFM (`get_best_corrections`) | done; covered by the `get_best_corrections` oracle below |
| context windows (`get_contexts`) | done; 5 922 real vectors byte-identical |
| frequency context matrix (`test_numba`) | done, no numpy; verified through the alternatives it feeds |
| alternative references (`get_alternative_ref_contexts`) | done; 636 real alternatives replayed, order included |
| the correction loop over `PFM` (`get_best_corrections`) | done; 5 922 of 5 922 real calls byte-identical |
| semi-global alignment (parasail `sg_trace_scan_16`) | done natively, no C++; tie-break measured, 2 408 real alignments byte-identical. Kept behind `ISONCORRECT_EXACT_GUARD=1` |
| the guard's production aligner | `block-aligner` (SIMD, banded). Optimal score on 1 400/1 400 recorded alignments; differs from parasail only in which equally-optimal path it reports, changing ~0.8% of reads |
| `fix_correction` | done; replayed against every recorded call |
| `correct_read` (stitching + guard) | done; replayed end to end against every recorded read |
| per-read driver (`isoncorrect_main`) | done; `corrected_reads.fastq` byte-identical on 1 204 reads |
| `run_isoncorrect` (batch driver) | done; in-process threads, all 29 equivalence cases green |

Argument names, defaults, validation order, stderr text and exit codes match the reference.
**Both binaries now do real work**, and `bench/equivalence.sh verify` is fully green.

### Stage-by-stage verification

End-to-end equivalence says *that* the port is wrong, never *where*. `bench/dump_reference.py`
wraps the reference (without modifying it) and writes each stage's structures in a stable,
line-oriented format; `isoncorrect-dump` writes the same format from the Rust side. Diff them
directly:

```bash
bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/py
rust/target/release/isoncorrect-dump --fastq test_data/isoncorrect/0.fastq --outdir /tmp/rs
diff -r /tmp/py/batch_000 /tmp/rs/batch_000
```

Verified identical so far, on both test clusters:

- `reads.tsv` — parsing, header rewriting, `r_id` order
- `minimizers.tsv` at `--k/--w` of 9/20, 9/10, 11/20, 7/15, 15/40, 9/9 (3 230 to 19 008 minimizers)
- `anchors.tsv` and `anchor_keys.tsv` at six `--k/--w/--xmin/--xmax` combinations, including the
  paper settings (90 036 entries) — and the Rust `generated` counter matches the reference's own
  printed "25953 MINIMIZER COMBINATIONS GENERATED"

Pure functions get a stronger check than diffing: `bench/dump_reference.py` records every
`solve_WIS` call's inputs *and* outputs, and a Rust test replays them. Hand-written cases cannot
cover float weights and strict-`>` tie-breaks; real ones do.

```bash
bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/wd -- --k 9 --w 10
WIS_DUMP=/tmp/wd cargo test --manifest-path rust/Cargo.toml wis::replay -- --nocapture
```

Replayed clean at `--k/--w` of 9/20 (100 calls, 1 342 intervals), 9/10 (100 calls, 20 689) and
11/25 (78 calls, 313). The test skips silently when `WIS_DUMP` is unset.

`spans.tsv` records what every `find_most_supported_span` call appended to `all_intervals`, tagged
with the anchor that produced it and the full payload.

**Verifying it needs `--exact`, and that is not a workaround — it is the only way to break the
circularity.** `previously_corrected_regions` is populated by the *correction* stage, so ordinarily
read 2 onwards depends on output that is not ported yet. But `isoncorrect_main` clears that dict at
the start of every read under `--exact`, making span-finding independent of correction. An
`--exact` reference run is therefore directly comparable, for all reads, today.

Byte-identical on both clusters at `--k/--w` of 9/20 (2 208 intervals, 6 859 alignments), 9/10
(**52 558** intervals, 128 062 alignments), 11/25 (429) and 7/15 (10 199). That covers the payload
ordering, the `already_computed` cache across anchors, the `reads_visited` fall-through and the
float thresholds.

~~Still unverified: the `previously_corrected_regions` / `pos_group` filtering.~~ **Now covered.**
It is inert under `--exact`, but `--exact_instance_limit` defaults to **0** in `isONcorrect` itself,
so the 100-read fixture runs the *non-exact* path — the feedback loop is live there, and the
end-to-end output is byte-identical. The circularity that made this hard to check disappeared once
the correction stage existed to close the loop.

`cigar_to_seq_cases.tsv` and `msa_cases.tsv` do the same for the two stages that turn pairwise
alignments into columns. Both are pure functions of their recorded inputs, so the Rust side replays
them rather than diffing a dump:

```bash
bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/d -- --k 9 --w 10
CIGAR_TO_SEQ_CASES=/tmp/d/cigar_to_seq_cases.tsv MSA_CASES=/tmp/d/msa_cases.tsv \
  cargo test --manifest-path rust/Cargo.toml oracle -- --nocapture
```

Green over a 16-run sweep (both clusters × `--k/--w` of 9/20, 9/10, 11/25, 7/15, 15/40, plus
`--max_seqs_to_spoa 3`, `--T 0.05` and `--xmin 30 --xmax 120`): **6 269 CIGAR expansions** and
**16 572 matrix rows across 1 622 multialignments**, byte-identical. The same sweep re-ran the two
older oracles on a wider corpus than they were recorded against — 6 378 edlib CIGARs and 674 spoa
consensus sequences, both still clean.

Three facts that sweep settled, each of which would otherwise have been a guess:

- **`get_best_solution`'s fourth branch is never reached.** Counting branch hits inside the
  reference at `--k 9 --w 10`: 165 "no insertion", 238 substring, 31 threaded through `min_ed`,
  and **0** for the shift-left fallback. So the corpus exercises three of four strategies, and the
  fallback rests on unit tests alone. It is also the branch whose Python is loosest — a `max_p > 0`
  guard whose false case falls through to code producing the same string — so a corpus that reaches
  it is worth having.
- **`min_ed` calls edlib again, and it shows up in the counts.** At default parameters the
  reference makes 677 `task="path"` calls but only 650 `cigar_to_seq` calls; the 27-call gap is
  `min_ed` fitting insertions into a widened slot. The CIGAR tie-break therefore reaches the matrix
  through two independent paths.
- **Partition keys cannot collide.** `partition[(q_id, pos1, pos2)]` looks like it could silently
  drop a duplicate segment, but `find_most_supported_span` builds `to_add` keyed by `r_id`, so an
  instance holds each read at most once. Confirmed on the corpus: matrix rows equal
  `cigar_to_seq` calls plus one consensus row per matrix, in every run. No dedup ever happens, and
  the port does not need an order-preserving map here.

~~**The PFM has no reference oracle.**~~ **It does now.** `correction_cases.tsv` records every
`get_best_corrections` call — the interval's segments in, the corrected spans out — so the PFM is
verified through the only thing that consumes it. Segments are recorded rather than whole reads,
which is exactly what the function reads, so a case replays without the fastq:

```bash
CORRECTION_CASES=/tmp/d/correction_cases.tsv \
  cargo test --manifest-path rust/Cargo.toml corrections::oracle -- --nocapture
```

**All 5 922 recorded calls are byte-identical.** They were 5 916 of 5 922 until the reference's own
non-determinism was fixed — see the next section.

Because this oracle runs the whole of `get_best_corrections`, one run of it re-checks every stage
underneath at once. Over the 12-run corpus (the fixture, 4 SIRV transcript clusters, 3 gene
clusters, and `--k 9 --w 10`, `--T 0.05`, `--max_seqs_to_spoa 5`, `--exact`):

| Stage | Cases, all byte-identical |
| --- | --- |
| spoa consensus | 5 922 |
| edlib `task="path"` CIGARs | 187 744 |
| `cigar_to_seq` expansions | 186 330 |
| multialignment rows | 192 252 across 5 922 matrices |
| context vectors | 5 922, carrying 636 alternatives |
| `get_best_corrections` | 5 922 |

`context_cases.tsv` records every `get_alternative_ref_contexts` call: the matrix it was given, the
context windows, and the alternatives it returned. The Rust side rebuilds the windows and the
frequency context matrix from the same rows and replays it.

```bash
CONTEXT_CASES=/tmp/d/context_cases.tsv \
  cargo test --manifest-path rust/Cargo.toml contexts::oracle -- --nocapture
```

Green on **5 922 context vectors and 636 alternatives** across 12 runs. Two things that corpus
made obvious, and neither would have shown up on the fixture:

- **Alternatives need isoform variation to appear at all.** One transcript per cluster at 7% error
  gives 68 alternatives across 1 987 calls; putting the isoforms of one gene together gives 568
  across 3 398. The stage exists for allelic difference, and a corpus without any barely tests it.
- **Error-free reads produce exactly zero.** All three `err0` clusters returned no alternatives from
  125 calls: with no errors every row carries the same context, so every column fails the
  "more than one variant" guard immediately. A clean corpus is the wrong corpus here.

The frequency context matrix has no oracle of its own — it is an intermediate the reference never
returns. It is verified transitively: a wrong count, a wrong eligibility cutoff or a wrong order
changes which alternatives come out, and those are compared exactly.

**Alternatives are compared in order**, and that only became meaningful once the reference was
fixed. It used to return a Python `set` per column, whose iteration order depends on
`PYTHONHASHSEED` — so there was no defined order to compare against, and the oracle compared the
alternatives as a set. `get_alternative_ref_contexts` now returns a list in insertion order, the
order is part of the answer, and the oracle checks it as such. That matters because the order is
exactly what `get_best_corrections` uses to pick a winner. See *The reference used to disagree with
itself* below.

### The reference used to disagree with itself, and that is now fixed

**`get_best_corrections` was not deterministic.** Same input, same flags, different corrected
regions, depending on `PYTHONHASHSEED`. That was a defect in the reference, not a porting problem,
and it accounted for every mismatch the correction oracle ever reported.

The mechanism was one line. `get_alternative_ref_contexts` returned a **`set`** per column, and the
correction loop iterates it:

```python
for (variant, ref_context, depth, thresh_) in alternative_refs[i]:
    if read_context == ref_context: ... break
    ed = edlib_alignment(...)
    if ed < min_ed: ...          # strict: the first minimum wins
```

Set iteration order over tuples of strings depends on the per-process hash seed, so when two
alternatives tied on edit distance against a read's context, which one won was decided by the seed.

Measured before the fix, on one interval of a gene-level SIRV cluster, across 24 hash seeds:

| | |
| --- | --- |
| corrected regions the interval returns | 267 |
| regions that are **identical under every seed** | 216 |
| regions that take **2 to 5 different values** | 51 |

**The fix** (its own commit, as *Working agreements* requires): `alternative_contexts` is a list in
insertion order rather than a set. `eligible_contexts_hcomp` is a dict, so that order is
deterministic and highest-depth-first — which is the order the *intent* of the code implies, and the
order the port already produced.

Two measurements make this a safe change rather than a hopeful one:

- **Every one of the 43 recorded golden outputs is byte-identical after the fix.** The goldens did
  not need re-recording. Re-check with `bench/equivalence.sh record --golden /tmp/after` and diff.
- **The final `corrected_reads.fastq` was already seed-independent** on every cluster tried, before
  the fix (3–4 seeds each, including the gene cluster with 51 unstable regions). The unstable values
  land in `previously_corrected_regions` and were never reused on this corpus. That was never a
  guarantee, though — `correct_read`'s `if isinstance(instance, str): best_corr = instance` path
  feeds a stored region straight into the output, so a propagation route existed. The fix closes it.

Also worth recording, since it bounded the blast radius: **the read's own correction was never
affected**, only `other_corrections`, the spans computed for the *other* reads supporting an
interval. Across 5 922 calls, `curr_read` matched under every seed.

`bench/check_seed_sensitivity.py` is the regression check. It replays recorded cases under N seeds
and classifies every record as *stable and agreeing*, *seed-dependent*, or *stable and wrong*:

```bash
CORRECTION_CASES=/tmp/d/correction_cases.tsv CORRECTION_MISMATCHES=/tmp/d/mismatches.tsv \
  cargo test --manifest-path rust/Cargo.toml corrections::oracle
bench/check_seed_sensitivity.py --cases /tmp/d/correction_cases.tsv \
  --mismatches /tmp/d/mismatches.tsv --seeds 24
```

**The seed count matters** — a lesson worth keeping. Before the fix, at 6 seeds, 6 records looked
like the port had invented an answer; at 24 all of them turned out to be reference answers the
smaller sample had missed. Do not draw conclusions from a short run.

The two options not taken, for the record: reproducing CPython's set order in Rust (SipHash-1-3 with
the `PYTHONHASHSEED=0` key, tuple hashing, CPython's set probing — exact, and permanently strange),
or excluding these regions from the equivalence contract as `--randstrobes` is. Fixing the reference
was cheapest, removed a real defect from a published tool, and cost nothing in output.

**The harness's hash-seed pin was fake, which is how this stayed hidden.** `dump_reference.py` did
`os.environ.setdefault("PYTHONHASHSEED", "0")` *inside* the running interpreter — CPython reads that
variable once at startup, so the call did nothing and every dump was randomly seeded. It now
re-execs itself, and repeat dumps are byte-identical. `profile_python.py` had the same line;
timings do not depend on it, so it is simply gone. `equivalence.sh` and `benchmark.sh` export the
variable in the shell before launching Python, which does work.

### The structural-overcorrection guard

Three oracles, because the three pieces fail independently:

```bash
bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
  FIX_CORRECTION_CASES=/tmp/d/fix_correction_cases.tsv \
  CORRECT_READ_CASES=/tmp/d/correct_read_cases.tsv \
  CORRECTION_CASES=/tmp/d/correction_cases.tsv \
  cargo test --manifest-path rust/Cargo.toml oracle -- --nocapture
```

Over the fixture plus 3 gene clusters, 3 transcript clusters and two parameter settings:

| Stage | Cases, all byte-identical |
| --- | --- |
| parasail semi-global alignment | 2 408 (scores **and** CIGARs) |
| `fix_correction` | 1 204 |
| `correct_read`, stitching and guard end to end | 1 204 reads |

`correct_read_cases.tsv` records intervals by **source** rather than by value — `c:<n>` for the
n-th `get_best_corrections` call, `s:<literal>` for a region carried over in
`previously_corrected_regions` — so the harness never has to reimplement any of the reference's
logic to say what a corrected span should be.

**parasail is reproduced natively, and both facts it rests on were measured.** First, what `sg`
means: plain `sg` leaves gaps free at **both ends of both sequences**, confirmed by running the
library — `sg("ACGTACGT", "TTTACGTACGTTTT")` returns `3D8=3D` scoring 32, exactly 8 matches at +4
with nothing charged for the end gaps. Second, the tie-break. Affine costs plus free end gaps leave
many alignments at the optimal score, so the same problem as edlib's CIGAR appears again — and the
same method settles it. Sweeping all 96 combinations of four tie-break parameters:

| | CIGARs matching |
| --- | --- |
| diagonal, insert, delete + **extend** over open | **200 / 200** |
| diagonal, insert, delete + **open** over extend | 182 / 200 |
| everything else | worse |

So parasail prefers the diagonal, then a gap in `s2`, then a gap in `s1`, and on a tie inside a gap
chain it **extends rather than opens**. The naive open-first choice scores 91%, which is exactly the
kind of near-miss that end-to-end testing would have found late and localised badly.

Two of the four parameters — which of the last row and last column to scan first, and whether to
keep the first or last of several equally-good end cells — are **unpinned by the corpus**, and the
sweep reports all four combinations as perfect. That is not laziness: the two sequences are a read
and its own correction, so a full-length alignment always beats stopping early and the end-cell tie
never arises. The values are recorded as arbitrary in `parasail.rs`, and the sweep is the tool to
re-run if a corpus ever reaches one.

`PARASAIL_SWEEP` doubles as a case cap (`PARASAIL_SWEEP=40`), because 96 configurations over
1 400 bp reads is hours of quadratic DP and ties are what matter, not volume.

**The port is now free of C++ entirely.** edlib, spoa and parasail were the three native libraries
the reference depends on; the first and third are reimplemented natively and the second goes through
`spoars`, which is pure Rust. No CMake, no vendored source, no toolchain beyond cargo.

Two things about `fix_correction` worth knowing before touching it:

- **an indel run mixes both directions.** `l` counts columns of either kind while both `o_segm` and
  `c_segm` accumulate, so a deletion immediately followed by an insertion is *one* run — and when it
  flushes, only one side survives. Reading it as two independent runs gives the wrong answer, and a
  unit test pins the case.
- **the fourth branch is dead.** `else: raise("Parsing alignment error...")` cannot be reached: the
  `elif o == '-'` above already absorbs the gap-aligned-to-gap case. It would also raise a
  `TypeError` rather than the intended error, since `raise` needs an exception, not a string.

### The batch driver

`run_isoncorrect` is orchestration rather than algorithm: enumerate the cluster fastqs, optionally
split large ones, correct each, optionally join the pieces back. The per-cluster result is already
byte-identical, so what is left to get right is *which reads end up in which file*.

**Clusters are independent**, which is what makes the parallelism safe to reorder: each one reads a
single fastq and writes a single output folder, with no shared state. The reference hands them to a
`multiprocessing.Pool` that spawns `python isONcorrect.py` per cluster; the port runs the same work
on an in-process thread pool, so one interpreter startup per cluster disappears.

`bench/equivalence.sh` compares **every `corrected_reads.fastq` by relative path**, so the folder
layout is as much part of the contract as the bytes. Four quirks decide it, and all four are
reproduced rather than tidied up — see *Known bugs* for why each is a defect:

- the size test in `split_cluster_in_batches` **latches** after the first cluster that fits;
- `--split_wrt_batches` builds `<cl_id>_<i>.fastq`, and *unsplit* clusters are symlinked as
  `<cl_id>_0.fastq`, so they too go through the join;
- joining concatenates in `sorted(glob(...))` order, which is lexicographic — batch 10 before batch
  2;
- `--keep_old` compares `wc -l` of output against input, so a truncated output with a coincidentally
  matching line count would be kept.

**The dump tool must apply the same argument massaging `main` does.** A mismatch at
`--xmin 14 --k 9` turned out to be the dump binary, not the anchor logic: `main` clamps `--xmin` up
to `2*k` before anything runs, so the reference was building spans at 18 while the port used 14.
Anything added to `main`'s preamble has to be mirrored in `dump_stages.rs` or the comparison lies.

Current `bench/equivalence.sh verify`: **29 passed, 0 failed.** All 16 single-cluster cases, all 6
`run_isoncorrect` folder cases — including `--split_wrt_batches`, `--split_mod`/`--residual` and
`--t 1` — and the 7 unsupported-flag cases. That is the whole sweep, byte-identical.

Verified equal to the reference by hand: `--version` on both binaries, exit codes for
`--w 150` (1), `--split_mod 2 --residual 5` (1), `--version` (0), no-args (0), and the exact stderr
text of `xmin set to 40` and the window-range message. Help *layout* differs — clap leads with the
description, argparse with usage — which is cosmetic and outside the contract.

```bash
cargo build --release --manifest-path rust/Cargo.toml
cargo test --manifest-path rust/Cargo.toml
```

### edlib: native for distance, native for the CIGAR too

isONcorrect calls edlib two ways, and they start from **different portability**:

| Call site | Uses | Uniquely defined? |
| --- | --- | --- |
| `edlib_alignment` → `edlib.align(x, y, "NW", 'dist', k)` | the **distance** only | **Yes.** Edit distance is a uniquely defined number, so any correct implementation agrees with edlib by definition. Done: `editdist.rs` on `triple_accel`, verified against Python edlib over 4 000 real substring pairs. |
| `get_best_corrections` → `edlib.align(seq, spoa_ref, task="path", mode="NW")` and `correct_seqs.min_ed` | the **CIGAR** | **No.** Many alignments achieve the same optimal distance and edlib picks one by its own tie-breaking. A different tie-break returns a different, equally optimal alignment — changing the MSA and the corrected output. |

**Resolved: the CIGAR is reproduced natively**, by `align.rs`, and the tie-break was *measured*
rather than reasoned about. Running the traceback under all six preference orderings against the
recorded corpus, only one — insertion, then deletion, then the diagonal, walking backwards from
`(n, m)` — reproduces edlib; the intuitive diagonal-first ordering gets 8%. Byte-identical on
**6 378 real `task="path"` calls**. The rest of this section records how that decision was reached,
since the alternatives are what to fall back on if the tie-break ever proves incomplete.

There is no native Rust *port* of edlib; `edlib_rs` and `rsedlib` are both C++ bindings.
`edlib_rs` additionally fails to build against CMake 4 (its vendored `CMakeLists.txt` requires
`cmake_minimum_required` compatibility that CMake 4 removed) — it only builds with
`CMAKE_POLICY_VERSION_MINIMUM=3.5` set, and pulls in `xz2`, `winapi` and `buf_redux`. Avoided for
the distance path.

**Measured.** `bench/dump_reference.py` now captures every `task="path"` call as
`cigar_cases.tsv` (query, target, edit distance, CIGAR) — 677 calls on cluster 0 at defaults. Tested
against that oracle:

| Crate | Result |
| --- | --- |
| `edlib_rs` 0.1.2 | **677/677 CIGARs identical**, edit distances too. But only builds with `CMAKE_POLICY_VERSION_MINIMUM=3.5` set, because its vendored `CMakeLists.txt` requires compatibility CMake 4 removed. Pulls in `xz2`, `winapi`, `buf_redux`. |
| `rsedlib` 0.1.1 | Compiles, then **fails to link** — its build script does not produce the vendored library (`ld: library 'edlib' not found`). Unusable here. |
| native Rust (`align.rs`) | **Chosen.** 6 378/6 378 CIGARs identical once the tie-break was measured. Pure Rust, no toolchain, no vendored source. |

So three options existed. The one taken was **(2)**, because the oracle made the tie-break an
empirical question rather than a guess:

1. **`edlib_rs` + documented CMake workaround.** Correct, at the cost of a C++ toolchain and an env
   var users must know about (or a `build.rs` that sets it).
2. **Reimplement edlib's NW traceback natively**, validated against `cigar_cases.tsv`. Keeps the
   build pure Rust — the port is otherwise C++-free — but the tie-breaking must be replicated
   exactly, not merely optimally.
3. **Vendor edlib's C source directly** with a `cc`-based build, skipping CMake entirely. edlib is a
   single self-contained `.cpp`, so this sidesteps the toolchain fragility while keeping exactness.

**The corpus is the guarantee, and it is still small.** `align.rs` is plain `O(n*m)`
Needleman-Wunsch, not edlib's bit-parallel Myers, so agreement rests entirely on the recorded cases
— all of which come from one 100-read cluster and segments bounded by `--xmax`. If a mismatch ever
appears, options 1 and 3 above are the fallback, and the `cigar_cases.tsv` oracle is what would
catch it.

### Porting gotchas found so far

- **clap rewrites `field_name` to `--field-name`.** The reference uses underscores throughout, so
  every multi-word flag silently became a different flag (`--max_seqs` → `--max-seqs`) until each
  was pinned with an explicit `long = "..."`. `cli.rs` has tests asserting the exact flag set and
  that no flag contains a hyphen; keep them.
- `--T` and `--t` are different flags, and `--t` exists only on `run_isoncorrect`. Case matters.
- `Cargo.toml` cannot hold the four-component version `0.1.3.5`; the crate is `0.1.3` and the CLI
  string lives in `cli::VERSION`. Keep them in step.
- **Python slices clamp, Rust slices panic.** `seq[i:i+k]` past the end silently truncates, and the
  truncated k-mer then takes part in lexicographic comparison — an empty string sorts before
  everything. Consequence, verified against the reference: **any read too short to fill the opening
  window returns exactly `[("", 11)]` at `--k 9 --w 20`** — one minimizer, empty k-mer, at a
  position that does not exist in the read. Nonsense, but reproducing it is required. `minimizers.rs`
  has a test pinning it.
- **`del M2[m1][m2]` does not delete.** `get_minimizer_combinations_database` deletes singleton
  anchors, then immediately evaluates `M2[m1][m2]` again on the next line — and because `M2` is a
  `defaultdict`, that *recreates* the key with an empty array. Measured on cluster 0 at default
  parameters: **23 260 keys, of which 21 801 are empty** — 94% of the index is dead weight.

  **The port does not reproduce this**, and does not need to: the only access anywhere is a lookup
  (`isONcorrect.py:1022`), never an iteration, and a lookup on an empty entry is indistinguishable
  from a lookup on an absent one — both mean "no support". Rust stores only non-empty entries and
  returns an empty slice on a miss. Same data, ~6% of the keys. This is a data-structure choice, not
  a behaviour change; see *Working agreements*.

  Consequence for the dumps: `anchor_keys.tsv` lists only non-empty keys **on both sides**, so the
  files diff cleanly, and the empty-key count is recorded in `meta.txt` instead.

## The pipeline, in one pass

Per cluster fastq, reads are batched into `--max_seqs` (default 2000) chunks. For each batch:

1. **Minimizers** — `get_minimizers_and_positions`, `hash_fcn` hardcoded to `"lex"`, so minimizers
   are chosen by **lexicographic order on the k-mer string**, not a hash. Ties within a window
   resolve to the **last** occurrence (`rindex`).
2. **Anchor pairs** — `get_minimizer_combinations_database` builds, for every read, pairs of
   minimizers `(m1, m2)` whose span is in `[--xmin, --xmax]`, indexed so other reads sharing the
   same pair can be found.
3. **Per read**, for each anchor pair, `find_most_supported_span` collects the supporting segments
   from other reads (edlib is used only when an anchor repeats within a read).
4. **Weighted Interval Scheduling** (`solve_WIS`) picks a non-overlapping set of high-support
   intervals to correct over.
5. **Per interval** (`get_best_corrections`) — build a consensus with **spoa**, align every
   supporting segment back to it with edlib (`task="path", mode="NW"`), build an MSA matrix and a
   position frequency matrix, then emit the corrected substring subject to `--T`.
6. **Structural-overcorrection guard** (`correct_read`, end). The corrected read is aligned back
   against the *original* read with parasail semi-global
   (`match=4, mismatch=-8, gap_open=12, gap_ext=1` via `sg_trace_scan_16`, falling back to `_32` —
   note `parasail_alignment` *defaults* to `opening_penalty=24`, but `correct_read` passes 12),
   and `fix_correction` walks that alignment: **any indel run longer than 10 alignment columns is
   reverted to the original read's sequence**, on the assumption it is a structural difference
   (e.g. exon variation) rather than an error. Runs of ≤10 keep the correction. This is the last
   thing that happens to a read and it is easy to miss — it is not optional.
7. Corrected regions are recorded in `previously_corrected_regions` so later reads skip work
   already done. **This makes reads order-dependent** — reads are processed in `sorted(reads)`
   order (r_id = 1-based input order). Preserve this exactly.

Output fastq preserves input headers, sequences are corrected, and the **quality string is
`"+" * len(seq)`** — not real qualities. Records are written in r_id order.

## Scope: what gets ported

Identical output is promised for the **supported set** below, not for every flag the Python CLI
accepts. Several options are experimental, deprecated, or broken, and reproducing them bit for bit
would cost more than they are worth. Dropping them is a deliberate decision, recorded here.

### Ported — inside the equivalence contract

`--fastq`, `--outfolder`, `--version`, `-h`, `--k`, `--w`, `--xmin`, `--xmax`, `--T`,
`--max_seqs`, `--max_seqs_to_spoa`, `--exact`, `--exact_instance_limit`, `--set_w_dynamically`

`run_isoncorrect`: `--fastq_folder`, `--t`, `--split_wrt_batches`, `--split_mod`, `--residual`,
`--keep_old`

> **`--exact` is not optional.** `run_isoncorrect` defaults `--exact_instance_limit` to 50, and
> `isoncorrect_main` sets `args.exact = True` whenever the batch has ≤ that many reads. Every
> cluster of 50 reads or fewer therefore runs the exact path by default. It is a default path
> wearing an optional flag's clothing.

### Dropped — not implemented

Two tiers, split by whether the flag has any working users to protect.

**Tier 1 — removed entirely.** The flag does not exist in the Rust CLI. The argument parser's own
unknown-flag error is the correct and sufficient response; do not carry a special case for these.

| Flag | Why it has no users to protect |
| --- | --- |
| `--disable_numpy` | Broken in the reference — it crashes (see *Known bugs*). Nobody can be running it successfully today. |
| `--compression` | Deprecated in its own help text, and **not exposed by `run_isoncorrect` at all**, so it is unreachable from the batch driver that real pipelines use. Only settable on a direct single-cluster invocation. |

Both are also absent from `run_isoncorrect`, so removing them cannot break a batch pipeline.

**Tier 2 — recognised and rejected.** These *work today* in Python, so they may appear in existing
pipeline scripts. The flag must parse, then exit non-zero with a message naming it and saying it is
not supported in the Rust port.

| Flag | Why dropped |
| --- | --- |
| `--randstrobes`, `--layers`, `--set_layers_manually` | Experimental; did not improve accuracy. Also **provably not reproducible**: output changes with `PYTHONHASHSEED` (measured — seeds 0/1/2 all differ), so there is no fixed reference behaviour to match. |
| `--use_racon` | Shells out to the `racon` binary, reintroducing exactly the subprocess-per-interval cost the port exists to remove. The code path also contains unconditional debug `print`s. |

Do not silently accept and ignore a Tier 2 flag: a pipeline that passes `--use_racon` and gets
unpolished output with a zero exit status is worse off than one that fails loudly. A generic
"unrecognised argument" error is not good enough here either — the message should be specific
enough to act on.

Either way the observable contract is the same and is what `bench/equivalence.sh` asserts:
**non-zero exit, and the flag named in the output.**

### Diagnostic only

`--verbose` prints developer statistics to stdout/stderr. Measured: it does **not** change
`corrected_reads.fastq`. Port it as best-effort; matching the Python diagnostic text word for word
is not required and is not part of the contract.

## spoa

spoa has to be replaced regardless of profiling, because a Rust port cannot shell out to a Python
subprocess wrapper. Today, **every correction interval** does: write a temp fasta → `fork/exec` the
`spoa` binary → write stdout to a temp file → read it back, so process spawn and disk I/O sit in
the innermost loop.

**Do not assume this is the dominant cost.** It is one plausible hot spot among several — building
the minimizer-pair database, `find_most_supported_span`, the MSA/PFM construction, and WIS are all
candidates. Profile first (see *Profiling*) and let measurements set the priority.

First measurements (`bench/PROFILE.md`, 100-read cluster) put spoa at 69% of self time at default
parameters but only 51% at paper parameters, where `get_minimizer_combinations_database` rises from
1.7% to 13%. **The profile is a function of the parameters and the cluster size, not a property of
the program**, and it has not yet been measured on a cluster large enough to matter.

The only invocation on the hot path is `create_augmented_reference.run_spoa`:

```
spoa <reads.fa> -l 0 -r 0 -g -2
```

Resolving that against spoa's CLI defaults (`m=5, n=-4, g=-8, e=-6, q=-10, c=-4`) and
`AlignmentEngine::Create`'s subtype dispatch (`g >= e` ⇒ linear, and then `e := g`), the effective
configuration is:

- **Alignment type:** `kSW` (local / Smith-Waterman)
- **Gap model:** **linear**, `gap = -2` per gap position
- **Match `m` = +5, mismatch `n` = -4**
- **Result:** consensus only (heaviest-bundle traversal)

Linear-gap local POA is a much smaller target than affine or convex. Only the consensus is used —
the MSA output path is not on the hot path.

### Measured: `spoars` reproduces spoa exactly on this configuration

`bench/dump_reference.py` captures the sequences handed to `run_spoa` and the consensus it returns,
giving a differential oracle. Against the real `spoa` binary (bioconda 4.1.5):

**505 / 505 unique correction intervals produced an identical consensus**, gathered across 8
parameter combinations (`--k/--w` of 9/20, 9/10, 11/25, 7/15, plus `--max_seqs_to_spoa 3` and
`--T 0.05`) on both test clusters, with up to 28 sequences per POA and a mean of 8.3.

The mapping from isONcorrect's invocation to the spoars API:

```rust
// spoa <fa> -l 0 -r 0 -g -2  ->  kSW (local), linear gap -2, m=+5, n=-4
let scoring = Scoring::new(5, -4, -2, -2, -2, -2)?;
let mut engine = SimdEngine::new(AlignmentType::Local, scoring);
// then per sequence: engine.align(seq, &graph) -> graph.add_alignment(&aln, seq, &weights)
```

`spoars` is pure Rust, builds in seconds, and needs no C++ toolchain — which also keeps the port
free of the CMake fragility that made `edlib_rs` unusable.

Caveats that keep this a decision rather than a formality: `spoars` is v0.1.3 with low usage, and the
corpus above comes from a single distinct 100-read cluster (the two fixtures are the same file), so
sequence diversity is limited. The differential oracle is cheap to re-run and should stay a
permanent gate, not a one-off.

**Decision: adopt `spoars` now, write a native POA later.** It gets the port to end-to-end
equivalence sooner without a C++ toolchain, and the oracle proves any future replacement is
faithful when it lands. The native implementation is a separate piece of work, tracked under
*Deferred improvements* — not a prerequisite.

Costs accepted with this choice:

- **MSRV rose from 1.74 to 1.85**, which `spoars` requires.
- A v0.1.3 dependency with low usage sits on the critical path. The oracle
  (`cargo test poa::oracle` with `SPOA_CASES` set) is the mitigation and must stay green on every
  upgrade.

Remaining options if `spoars` ever fails that gate:

1. **Reimplement linear-gap kSW POA + heaviest-bundle consensus natively.** The oracle above makes
   this a tractable, checkable project rather than a leap.
2. **`spoa` / `spoa-sys` bindings.** Identical by construction, at the cost of a C++ dependency.

`abpoa-rs` / `rsabpoa` (adaptive banded) and `poasta` (optimal gap-affine) implement *different*
algorithms and will not reproduce spoa's output, so they were never candidates while byte-identity
was the specification. It no longer is, and they have since been measured and rejected on their own
merits — abPOA aborts the process, `poa-consensus` is 11.6x slower, `poasta` has no consensus at all.
See *The POA bake-off*.

**Sequence insertion order into the POA graph changes the consensus.** Keep the exact order in
which `get_best_corrections` writes segments to the temp fasta, including the `i > max_seqs_to_spoa`
cutoff — note it is `>`, not `>=`, so it admits `max_seqs_to_spoa + 1` sequences.

## Determinism rules

The default code path is deterministic and must be reproduced exactly.

- **Do not use `hash()` semantics from Python.** Python's string `hash()` is randomized per process
  (`PYTHONHASHSEED`). It appears in `get_kmer_minimizers` and in the randstrobe functions, but the
  default path uses `hash_fcn = "lex"` and never reaches them.
- **`--randstrobes` is not reproducible in Python across runs.** It depends on randomized `hash()`.
  A Rust port cannot be bit-identical to an arbitrary Python run here. Either exclude it from the
  equivalence contract, or pin the Python reference with `PYTHONHASHSEED=0` and port SipHash-1-3
  with Python's exact string-hashing. Decide this explicitly; do not paper over it.
- **The default path had one leak of its own, now closed.** `get_alternative_ref_contexts` returned
  a `set` and `get_best_corrections` iterates it, so corrected regions depended on `PYTHONHASHSEED`
  — measured, 51 of 267 regions on one real interval. It returns a list in insertion order now. See
  *The reference used to disagree with itself* under *Stage-by-stage verification*.
- Ties in minimizer selection go to the **last** position in the window.
- `sorted(reads)` is a numeric sort on 1-based r_id.
- Dict iteration in the reference is insertion-ordered (Python 3.7+). Where iteration order feeds
  output, use an order-preserving map in Rust — not `HashMap`.

## Verification

Equivalence is the acceptance criterion, and it is checked, not assumed.

- `bench/equivalence.sh` runs Python and Rust on the same inputs and diffs
  `corrected_reads.fastq` byte for byte. A non-empty diff is a failure.
- Cover both binaries: `isONcorrect` on a single cluster, and `run_isoncorrect` over a folder.
- Sweep every ported parameter that changes behaviour: `--k`, `--w`, `--xmin`, `--xmax`, `--T`,
  `--max_seqs`, `--max_seqs_to_spoa`, `--exact`, `--exact_instance_limit`, `--set_w_dynamically`.
  Also cover the paper settings `--k 9 --w 10 --max_seqs 1000`.
- Dropped flags (see *Scope*) are not swept for output. They need one test each asserting a
  non-zero exit and a message naming the flag.
- Corpus beyond the 100-read test clusters is required before trusting the port.

> **The "two test clusters" are one cluster.** `test_data/isoncorrect/0.fastq`,
> `test_data/isoncorrect/1.fastq` and `test_data/100reads.fq` all share git blob `64cae773` —
> one blob, three paths, identical MD5. CI, the README's install check, and the default
> `bench/` corpus therefore all exercise **a single 100-read input**, and the recorded goldens for
> cluster 0 and cluster 1 are identical for that reason. This is a smoke test and nothing more.
> Building a real corpus (large clusters, exon variation, low coverage, repeated anchors) is a
> prerequisite for trusting the port, not a nice-to-have. See `bench/README.md`.

### The SIRV corpus

`bench/make_sirv_corpus.py` builds a real one from simulated SIRV reads, splitting a fastq by the
source transcript in each header (`@read_1_from_SIRV612`). No aligner, no isONclust run, and
reproducible from the fastq alone.

Two groupings, and they are not interchangeable: `--group gene` puts isoforms of one gene in one
cluster, so reads differ by **whole exons**, and the largest cluster crosses `--max_seqs`;
`--group transcript` gives one transcript per cluster, where every difference is sequencing error.
On the 10k-read 7%-error set that is 7 clusters of 713–2242 reads, or 68 of 20–342.

What that buys, beyond size: **62 of the 68 transcript clusters are above the
`--exact_instance_limit` of 50**, so they run the non-`--exact` path — which the 100-read fixture
never does, because `isoncorrect_main` forces `--exact` at or below the limit. The
`previously_corrected_regions` filtering is live on this corpus for the first time.

Every oracle was re-run against it, and all four stages held on the first try:

| Run | matrices | matrix rows |
| --- | --- | --- |
| 9 transcript clusters at defaults | 5 210 | 377 704 |
| `--k 9 --w 10`, `--k 11 --w 25`, `--exact` on one cluster | 2 390 | 97 850 |
| 3 gene clusters, 200 reads each | 3 398 | 103 747 |
| 3 error-free transcript clusters | 125 | 32 956 |

That is well over 100× the previous corpus, on a different organism, at 7% error rather than the
fixture's much cleaner reads — and `cigar_to_seq`, the CIGAR tie-break, the multialignment matrix
and the spoa consensus were all byte-identical throughout, every run, first try.

The error-free set (`reads_10k_err0.fastq`) is worth keeping as a *contrast*, not as the main
corpus: it drives the consensus and MSA hard but produces zero alternative contexts, because with
no errors every read carries the same context. Corpora have to be chosen per stage.

Simulated reads are still not real ones: uniform error, no chimeras, no adapters. Real ONT cDNA
through isONclust remains worth adding.

Performance is measured on the same corpus: wall clock and peak RSS, Rust vs Python, single-core
and at `--t N`.

## Performance, measured at each step

### Where it started

On one 200-read SIRV gene cluster (~1.4 kb reads, default parameters, single core):

| | wall | user | sys | peak RSS |
| --- | --- | --- | --- | --- |
| Python reference | 16.1 s | 9.0 s | 6.5 s | 339 MB |
| Rust, first working version | 10.7 s | 10.6 s | 0.1 s | 420 MB |

1.5x faster and 24% *more* memory — well short of the goal at the top of this file. The shape said
why: Python's 6.5 s of `sys` time is one `spoa` `fork`/`exec` per correction interval, which the
port removes outright, but the port *lost* on `user` time because the reference delegates its inner
loops to C (bit-parallel edlib, SIMD spoa and parasail) where the port had scalar `O(n*m)` DP.

And on a realistic workload — 68 SIRV transcript clusters (20–342 reads, ~1.4 kb), `--t 8`, both
implementations parallel — it was **dead even**: Python 83.7 s, Rust 84.1 s. That was worth stating
plainly, because it corrected the assumption behind porting the batch driver for speed: the
reference is *already* parallel via `multiprocessing`, so removing the subprocess and interpreter
startup only buys back what per-cluster inefficiency gives away. **Per-cluster speed is the lever**,
and it multiplies through the pipeline.

### Bounded edit distance: 67% of runtime, and an architecture trap

Profiling the port (`sample`, 5 064 samples, the 200-read gene cluster) put two thirds of runtime in
one function:

| component | share |
| --- | --- |
| `find_most_supported_span` → **bounded edit distance** | **67%** |
| the parasail guard | 19% |
| `get_best_corrections` in total (spoa 2%, NW CIGAR 3%, MSA <1%) | ~10% |
| everything else (minimizers, anchors, WIS, I/O) | <1% |

The cause was not the algorithm but the architecture. **`triple_accel` guards its SIMD paths with
`target_arch = "x86"` / `"x86_64"` and has no `aarch64` equivalent** — 67 occurrences of those two,
zero for ARM — so on Apple Silicon, Graviton, or any ARM host every call fell back to naive scalar
DP. A benchmark on an x86 laptop would never have shown this.

`rapidfuzz` implements the same bit-parallel Myers algorithm edlib uses, in plain `u64` arithmetic
with **no architecture guards at all**. Swapping it in is safe for the reason at the top of
`editdist.rs`: the answer is a uniquely defined number, so there is no tie-break to preserve.

| | wall | peak RSS |
| --- | --- | --- |
| before (`triple_accel`) | 10.7 s | 420 MB |
| after (`rapidfuzz`) | **3.9 s** | 420 MB |

**2.7x on one cluster, and 4.1x against the Python reference.** On the 68-cluster folder at `--t 8`:
**84.1 s → 24.2 s, a 3.5x speedup**, with all 68 outputs still byte-identical and the full
29-case sweep green. Memory is untouched, because the memory peak is the guard, not this.

### The guard: 3x less memory

After the edit-distance fix the parasail guard was **65%** of runtime and the whole of the memory
peak: three `i32` matrices over a read against its own correction, ~24 MB per read at 1.4 kb.

`parasail.rs` now keeps **two rows** of `H`/`E`/`F` and stores, per cell, only the four bits the
traceback cannot re-derive: which predecessor won the `H` cell (2 bits, from `TieBreak::h_order`),
and whether each gap chain leaves or extends (1 bit each). Same computation, so the recorded CIGARs
do not move.

| | wall | peak RSS |
| --- | --- | --- |
| three `i32` matrices | 3.9 s | 440 MB |
| two rows + a decision byte per cell | **3.5 s** | **147 MB** |

Two things that measurement settled, both the opposite of the obvious guess:

- **Four bits per cell is slower than eight.** Packing two cells per byte halves the array but adds
  a read-modify-write and a shift per cell, and the DP already streams row by row so the cache never
  benefited: 4.3 s against 3.5 s. One byte per cell it is.
- **Computing the tie-break winner by walking `h_order` cost ~10% of total runtime.** It is a
  data-dependent branch chain over a runtime array, executed once per cell, and the predictor cannot
  learn it. An 8-entry lookup table indexed by a 3-bit "which predecessors achieved this cell" mask
  removes the branches entirely.

This path is still what `ISONCORRECT_EXACT_GUARD=1` runs, so it remains the byte-identical
reference implementation and the thing the equivalence gate exercises.

### The guard's aligner: a deliberate divergence

The guard aligns each read against its own correction, and after the edit-distance fix that single
`O(n*m)` DP was **65% of the port's runtime**. It is now `block-aligner` by default: SIMD (Neon,
AVX2, SSE2) and adaptively banded.

**It is not an approximation on this data, and that is the point.** The first measurement of it was
misleading because it compared block-aligner's *global* alignment against parasail's *semi-global*
CIGARs — two different problems. Comparing like with like, against exact affine global:

| | |
| --- | --- |
| optimal score found | **1 400 / 1 400** recorded alignments |
| suboptimal | **0**, worst loss 0, at `max_block` 256, 1 024 and 4 096 |
| identical CIGAR | almost never |
| guard output differs | 16 / 1 400 (**1.1%**) |

So block-aligner finds an optimal-scoring alignment every time and simply reports a *different*
equally-optimal path. Roughly 1% of the time that shifts an indel run across the guard's 10-column
threshold and changes the corrected read. The penalties map exactly, which the matching scores
prove: block-aligner charges `open + extend * (n - 1)`, the same convention as parasail, and its
`I`/`D` mean the same thing.

End to end against the Python reference:

| | | |
| --- | --- | --- |
| one 200-read gene cluster | 16.1 s -> **1.4 s** | **11.6x** |
| 68-cluster folder, `--t 8` | 83.7 s -> **13.3 s** | **6.3x** |
| reads differing from Python, simulated | **78 / 10 000** | 0.78% |
| reads differing, **real isONclust data** | **2 364 / 75 442** | **3.13%** |

**The divergence is four times larger on real data than on simulated data** — 3.13% against 0.78%.
That is the same lesson the banding experiment taught: real ONT reads have variable lengths and
messier indels, so the guard's alignment has more ties to break differently. Any statement about
this divergence has to be sourced from the real corpus, not the simulated one.

`ISONCORRECT_EXACT_GUARD=1` selects the exact parasail-compatible DP instead, which is
byte-identical. `bench/equivalence.sh` exports it, so **the 29-case gate still passes 29/29** and
still covers every other stage byte-for-byte. Sequences shorter than 64 bp also take the exact path,
where the SIMD minimum block size would dominate anyway.

**Why the reference uses semi-global at all, and why that matters less than it looks.** The helper
is copied: `parasail_alignment(s1, s2, match_score=2, mismatch_penalty=-2, opening_penalty=..., gap_ext=1)`
calling `sg_trace_scan_16` exists verbatim in isONclust's `modules/consensus.py`, where it compares
*different reads* and semi-global is exactly right. In `correct_read` it is applied to a read against
its own correction — which shares an exact prefix and suffix with it, so free end gaps can almost
never pay — and none of the helper's defaults are used at the call site. So global is the natural
choice here and semi-global looks inherited rather than chosen. It is **not inert**, though:
measured on 75 442 real alignments, exact global and exact semi-global disagree on 414 CIGARs and
**145 guard outputs (0.19%)**. Most of the 1.1% above is therefore tie-break, not this.

### Banding our own DP: measured and rejected

Before adopting block-aligner, banding the exact DP was tried. It is the cautionary tale of this
whole exercise, so the numbers stay:

| corpus | half-band 32 | half-band 64 | half-band 128 |
| --- | --- | --- | --- |
| simulated SIRV, 2 408 alignments | 12 differ | **0 differ** | 0 differ |
| real isONclust, 13 963 alignments | — | **23 differ** | — |
| real isONclust, 75 442 alignments | — | — | **34 differ** |

The simulated corpus endorsed a band width that is wrong on real data, and widening only lowers the
rate rather than reaching zero. Two further findings worth keeping:

- **An edge-touch check is not a sufficient validator.** At half-band 32, twelve alignments returned
  a different CIGAR while the path never touched the band edge.
- **The provable band is worth almost nothing.** Any path reaching deviation `W` needs `g >= W` gap
  columns, so it scores at most `2(n+m) - 2W`; a banded score above that is provably optimal. But
  `4*min(n,m) - score` has median 753 here, so the provable half-band is `~slack/2 ~= 380`.
  Measured: full DP 3.5 s, provable band 3.1 s, unproven +-64 band 2.1 s. **The rigorous version
  buys 9%.**

`ISONCORRECT_BAND_CHECK=<half>` re-runs the comparison on any dataset, and the `parasail::banding`
tests keep the measurement. block-aligner made this moot: it is banded *and* SIMD *and* optimal.

### The MSA's aligner: affine instead of edit distance

Once block-aligner took the guard from 65% of runtime to 0.9%, the largest single cost became
`align::global` — the scalar `O(n*m)` Needleman-Wunsch that aligns every supporting segment back to
the spoa consensus. **33.7% of runtime over ten real clusters**, roughly 61 000 alignments of ~80 bp
per cluster.

Reusing block-aligner here looks blocked: it asserts `gaps.open < gaps.extend`, i.e. **strictly
affine gaps**, so it cannot express the unit-cost edit distance edlib computes. The way past that is
not a workaround but a correction of the premise — **edit distance was the reference's choice for
speed in Python, not a modelling preference.** Affine gaps describe sequencing indels better, and
the reference already uses them where it can afford to, in the guard. So the port uses the guard's
own scoring here too: `match 4, mismatch -8, open 12, extend 1`.

Measured on the same ten real isONclust clusters, 30 000 reads:

| | `align` stage | total |
| --- | --- | --- |
| scalar unit-cost DP (the reference's objective) | 62.52 s, 33.7% | 185.5 s |
| block-aligner, affine | **38.77 s, 22.6%** | **171.9 s** |

The stage gains 1.6x rather than the 2.7x the isolated benchmark shows, because the real mix
includes the two fall-through cases below — sub-64 bp segments and `min_ed`'s gap-padded slots —
which stay on the scalar DP.

**The block width is derived, not tuned, and that is the whole difference from the rejected banding
scheme.** block-aligner grows its block greedily by score, and on a *global* alignment of
unequal-length sequences that heuristic can fail outright: the path must leave the diagonal by at
least `|len(a) - len(b)|`, and a block narrower than that never finds the gap. Measured at a fixed
block of 32:

| | |
| --- | --- |
| suboptimal alignments, one real cluster | **10 of 14 831** |
| every one of them | a ~43 bp segment against an ~86 bp consensus |

`simd::min_block` therefore sizes the starting block from the length difference — a *lower bound on
the required deviation*, which is exactly the guarantee the ±64 band never had. It costs ~15% of the
raw speedup. After it:

| | |
| --- | --- |
| recorded segment-vs-consensus alignments checked against the exact affine DP | **158 559 across five corpora** |
| suboptimal | **0** |
| clipped (did not span both sequences) | 0 |
| taking the *same path* as edlib's unit-cost CIGAR | ~89% |

That last row is the divergence, and it is the intended one: ~11% of alignments choose a different
path because affine gaps prefer one long gap where unit cost is indifferent between one long gap and
several short ones. `align::block_option` is the test; `parasail::global_affine` is its oracle.

**This is a deeper divergence than the guard's**, and worth being explicit about. The guard runs once
per read and its output is either "keep the correction" or "revert this run", so a different path
changes little. This CIGAR builds the MSA, which builds the PFM, which decides *every corrected
base*. `ISONCORRECT_EXACT_ALIGN=1` restores `align::global` and `bench/equivalence.sh` sets it.

#### Then the per-call cost, which was most of it

Profiling *inside* the stage found three things, and the first two were assumptions worth checking:

- **the 64 bp threshold was excluding the majority of the work.** These segments are median ~55 bp,
  so `simd::MIN_LEN = 64` sent **1 683 397 of 4 272 008 alignments (39%)** to the scalar DP. Below
  32 there is genuinely nothing to gain — that is block-aligner's minimum block size, so the block
  already spans the whole alignment — but between 32 and 64 there was. Lowering it leaves only
  **75 587** fallbacks, which are the gap-padded `min_ed` queries that cannot use SIMD at all;
- **per-call setup competed with the arithmetic.** A ~55x57 DP is a few thousand cells, and there are
  4.27 million calls per four real clusters, against which a fresh `Block`, two `PaddedBytes` and a
  `Cigar` (which allocates *and zeroes* `q + r + 5` 16-byte entries, ~2 KB here) are not free.
  `simd::Scratch` reuses all four. `Block::new` takes *maximum* lengths and `align` takes the actual
  sequences, so one block serves every smaller pair; and `cigar_eq` calls `Cigar::clear` itself, so a
  long-lived `Cigar` is safe even though `clear` is `pub(crate)`;
- **the CIGAR was round-tripping through a `String`.** The aligner reports operations natively and
  `cigar_to_seq` immediately parsed them back. `simd::global_ops` and `align::ops_to_seq` skip it.

| | four real clusters |
| --- | --- |
| `MIN_LEN` 64, fresh allocations, string CIGAR | 12.93 s + 7.75 s fallback + 1.13 s expand = **21.81 s** |
| `MIN_LEN` 32, reused scratch, ops | 11.54 s + 0.11 s fallback + 0.88 s expand = **12.53 s** |

**1.74x on the stage**, and 58.0 s -> **48.1 s** on the 32 real clusters at `--t 8`.

Lowering `MIN_LEN` changes output — it extends the same divergence to more alignments — so it was
checked the same way. `align::block_option` now exercises `simd::global_ops` itself rather than its
own copy of the aligner: **159 209 recorded alignments across six corpora, 0 suboptimal, 0 refused,
0 that fail to span their inputs.** And on accuracy it is better again:

| set | mean % | median % | p90 % | total % | perfect | paired vs Python |
| --- | --- | --- | --- | --- | --- | --- |
| affine, `MIN_LEN` 64 | 1.702 | 1.154 | **3.571** | 1.243 | 156 | -0.0129 pp |
| affine, `MIN_LEN` 32 | **1.696** | **1.143** | 3.608 | **1.232** | **158** | **-0.0185 pp** |

Better on 11 891 reads and worse on 6 614 — the same 64% ratio on a larger changed set. Four of the
five summary statistics improve; **p90 is marginally worse** (3.608 against 3.571), so the tail of
hard reads gives back a little of what the bulk gains. Worth recording rather than rounding away.

Two call sites stay on the exact DP regardless:

- **sequences shorter than 32 bp**, where the block already spans everything;
- **`msa::min_ed`**, whose query is a gap-padded insertion slot like `"-TT-"`. block-aligner's
  `NucMatrix` has no entry for `-` and `convert_char` asserts `A..=Z`, so `simd::global_ops`
  refuses non-ACGTN input and the caller falls through. 75 587 of 4 272 008 calls.

### The MSA matrix: one buffer instead of `2t + 1` allocations

`create_multialignment_matrix` was **14.3%** of runtime, most of it allocator traffic rather than
computation. The reference holds each positioned vector as a list of strings, and the port copied
that shape as `Vec<Vec<u8>>`: `2 * len(consensus) + 1` tiny allocations per supporting segment per
correction interval, nearly all of them holding the single byte `"-"`.

Every slot is either one character or a **contiguous run of the aligned query** — an insertion is
exactly the query characters opposite a run of consensus gaps — so all of them fit in one flat buffer
with an offset table. Two allocations per segment instead of `2t + 1`, and nothing is copied twice.

| | `create_multialignment_matrix` | total, ten real clusters |
| --- | --- | --- |
| `Vec<Vec<u8>>` per slot | 29.73 s (14.3%) | 207.9 s |
| one buffer + offsets (`msa::Positioned`) | **14.94 s (8.1%)** | **185.5 s** |

A representation change only, so the oracle must not move, and it does not: **181 205 rows across
5 484 matrices byte-identical** over the nine-run corpus, with the `get_best_corrections` oracle
above it green on all 5 484 calls.

### `find_most_supported_span`: 2.25x, and none of it was the alignment

After the MSA and `align` work this was the top cost — **28.4% over ten real clusters, 342 940
calls**, so per-call work rather than one hot inner loop. The obvious suspect was the bounded edit
distance. It was not: a temporary scope around it attributed only **27%** of the stage there, over
**15.4 million calls** on three clusters — 1 715 per read.

The other 73% was scaffolding around those 15 million iterations. Three changes, none of which touch
what is computed:

| | |
| --- | --- |
| **`FxHashMap`/`FxHashSet` instead of SipHash** on `added_strings`, `already_computed`, `reads_visited`, and on `to_add`'s lookups | 20.3 s -> 11.2 s |
| **the scratch containers reused across partners and calls** instead of three fresh allocations per anchor partner | included above |
| **recycling `added_strings`'s keys** rather than a fresh `Vec<u8>` per landed edit distance | 11.2 s -> 9.3 s |

All three are safe for the same reason: **those maps are looked up, never iterated.** `to_add` stays
an `IndexMap` precisely because *its* order becomes the payload order, and `IndexMap` keeps insertion
order regardless of hasher, so even its hasher is free to change. `ref_seq` also became a slice
rather than a `to_vec()` per partner.

| | ten real clusters | share |
| --- | --- | --- |
| before | 48.86 s | 28.4% |
| after | **21.70 s** | **14.8%** |

Total 171.9 s -> 147.0 s, and the stage is no longer the top cost — `align` is again, at 27%.

**Measured and rejected: caching the edit-distance pattern.** One `ref_seq` is compared against every
read carrying the anchor pair, so hoisting the bit-parallel pattern table out with
`levenshtein::BatchComparator` looks free. It is a **loss** — 4.11 s to 7.51 s over the same 15.4
million calls — because `distance_with_args` has a specialised path for patterns fitting one 64-bit
block while `BatchComparator` always builds the heap multi-block one. The note is in `editdist.rs`;
do not retry it without re-measuring.

**Measured and nearly neutral: the anchor database.** `AnchorDb` also moved to `FxBuildHasher` on both
levels and gained `partners(m1)`, so `find` pays the outer lookup once per call rather than once per
partner. That is **21.70 s to 21.78 s — no change**; the outer lookup was never hot. What the same
commit *did* buy is in the database build: `entry(m1.to_vec())` allocated both keys on every
generated pair, including the overwhelming majority that already existed, and looking up first took
`get_minimizer_combinations_database` from 5.43 s to **4.35 s**.

Verified by the stage oracle rather than by end-to-end output alone: `spans.tsv` byte-identical to
the reference on both clusters at `--k/--w` of 9/20, 9/10 (**52 558 intervals**), 11/25 and 7/15,
with `anchors.tsv`, `anchor_keys.tsv` and `minimizers.tsv` identical too, and 29/29 equivalence.

### Accuracy against ground truth

Equivalence testing answers "do two implementations agree"; it cannot answer "is the correction any
good". Once the guard's aligner diverges deliberately, the second question is the one that matters.
`bench/accuracy.py` measures it.

Truth comes from one of two places. Simulated reads carry their source transcript in the header, so
no search is needed. **Real reads are assigned by aligning against every transcript and taking the
best match**, ties broken by name so the choice is deterministic. Either way the read is scored
against its assigned transcript with edlib's infix mode (`HW`) and the error rate is
`edit_distance / len(read)`.

**The assignment is made once, from the first read set, and reused for all of them.** That is
load-bearing: if each implementation picked its own best-matching transcript, every one would be
flattered by construction, and a correction that dragged a read toward the wrong isoform would score
*better* rather than worse. Pass the uncorrected reads first.

```bash
bench/accuracy.py --transcriptome sirv_transcriptome.fasta \
  --reads uncorrected=/tmp/real_trunc --reads python=/tmp/real_py \
  --reads rust-exact=/tmp/real_ex --reads rust-block=/tmp/real_ba
```

**Real isONclust data — 32 clusters, 75 442 reads, the corpus that matters:**

| set | mean % | median % | p90 % | total % | perfect |
| --- | --- | --- | --- | --- | --- |
| uncorrected | 6.245 | 5.580 | 9.655 | 6.112 | 1 |
| Python reference | 1.787 | 1.168 | 3.672 | 1.310 | 156 |
| Rust, `ISONCORRECT_EXACT_GUARD=1` | 1.787 | 1.168 | 3.672 | 1.310 | 156 |
| Rust, block-aligner (default) | **1.782** | **1.163** | **3.665** | **1.303** | 156 |

Restricted to the 2 364 reads where block-aligner's answer differs at all:

| | |
| --- | --- |
| block-aligner **better** | **1 373** |
| block-aligner worse | 794 |
| equal | 197 |
| mean delta | **-0.157 pp** (negative is better) |
| median delta | **-0.203 pp** |

**block-aligner is measurably better on real data**: better on 63% of the reads it changes, and
better in aggregate (1.303% against 1.310% total error). Simulated data said the opposite — 30
better, 41 worse of 71 — which is a reminder that the simulated corpus is not a substitute for real
reads on any question about this guard.

Simulated corpus, for comparison (68 transcript clusters, 10 000 reads): uncorrected 6.866% total,
all three implementations 0.532%, and block-aligner differing on only 71 reads.

Two further conclusions:

- **The exact path is confirmed byte-identical through a second, independent lens** — identical error
  rate on all 75 442 real reads and all 10 000 simulated ones, measured against ground truth rather
  than against the reference's bytes.
- **Correction works**, which is worth stating since nothing else in this file checks it: 6.2% -> 1.8%
  mean error on real reads, improving 75 254 of 75 442 and worsening 116.

Caveats, since they are inferences rather than measurements. Assignment uses the *uncorrected* read,
so at ~6% error a read may be assigned to the wrong isoform of the right gene — SIRV isoforms are
deliberately similar. That inflates absolute error rates for every implementation equally, so the
comparison holds even where the absolute numbers are pessimistic. And `bench/accuracy.py` is
single-threaded: the best-of-68 search over 75 442 reads takes roughly 20 minutes, which is worth
parallelising before this is run routinely.

#### Is affine actually more accurate? Measured: yes

Byte-identity cannot answer this, which is the point — the whole reason to make this change is that
the reference's objective is not the better one. `bench/accuracy.py` on the 32 real isONclust
clusters, 75 156 reads with a plausible SIRV of origin, truth assigned from the **uncorrected** reads:

| set | mean % | median % | p90 % | total % |
| --- | --- | --- | --- | --- |
| raw | 6.171 | 5.573 | 9.557 | 6.063 |
| Python reference | 1.714 | 1.165 | 3.601 | 1.259 |
| Rust, exact mode | 1.714 | 1.165 | 3.601 | 1.259 |
| Rust, block-aligner guard only | 1.710 | 1.160 | 3.590 | 1.253 |
| Rust, **+ affine MSA aligner** | **1.702** | **1.154** | **3.571** | **1.243** |

Paired per read against the Python reference:

| | mean delta | better | worse | equal |
| --- | --- | --- | --- | --- |
| exact mode | +0.0000 pp | 0 | 0 | 75 156 |
| guard only | -0.0048 pp | 1 355 | 783 | 73 018 |
| **+ affine MSA aligner** | **-0.0129 pp** | **9 740** | **5 631** | 59 785 |

**Affine is better, and by 2.7x more than the guard divergence was.** It touches twenty times as
many reads (15 371 against 2 138) and still wins 63% of them — the same ratio the guard showed, on a
much larger sample, which is what makes it a signal rather than a coincidence. Every summary statistic
moves the same direction.

So this change is faster *and* more accurate, and the reference's use of edit distance here is
confirmed as the speed compromise it was. The effect is small in absolute terms (0.013 pp on a 1.7%
error rate); the claim is the sign, not the magnitude.

### Where it stands

All rows below are **after** the `fill_p2` fix, on both sides, so Python and Rust are doing the same
(larger) amount of correction.

| | Python | Rust default | |
| --- | --- | --- | --- |
| one 200-read gene cluster | 16.1 s, 339 MB | **1.4 s, ~147 MB** | 11.6x faster |
| 32 real SIRV clusters, 75 442 reads, `--t 8` | 432.2 s | **43.5 s** | **9.9x faster** |
| **266 real drosophila clusters, 54 721 reads, `--t 8`** | **181.2 s** | **18.5 s** | **9.8x faster** |

**Byte-identity in exact mode**, measured against the fixed reference:

| corpus | reads | clusters differing | reads differing |
| --- | --- | --- | --- |
| the 100-read fixture, 29 equivalence cases | — | 0 | 0 |
| 32 real SIRV clusters | 75 442 | **0 of 32** | **0** |
| 266 real drosophila clusters | 54 721 | 2 of 266 | **2** |

The two drosophila reads are the parasail semi-global end-cell tie-break, where the scores match and
only the reported path differs; both happen to be *better* by the accuracy measure. See *The parasail
tie-break the corpus finally reached*.

The two deliberate divergences that the default path carries on top are the guard's aligner (~0.8% of
reads) and the MSA's aligner (~11% of alignments); see *Goal*.

Gates: **29/29** equivalence cases (exact mode), stage oracles green on every recorded corpus
including drosophila, 195 unit tests, and `spans.tsv` byte-identical over eight dump comparisons.

Where the time goes now, ten real clusters and 30 000 reads, after the MSA and `align` work:

| stage | calls | self | share |
| --- | --- | --- | --- |
| `align` (segment vs consensus) | 11 775 | 22.11 s | 21.1% |
| `run_spoa` | 11 775 | 21.21 s | 20.2% |
| `find_most_supported_span` | 342 940 | 20.55 s | 19.6% |
| `create_multialignment_matrix` | 11 775 | 9.66 s | 9.2% |
| `get_alternative_ref_contexts` | 11 775 | 6.99 s | 6.7% |
| `get_minimizer_combinations_database` | 20 | 4.15 s | 4.0% |
| the guard | 30 000 | ~2.4 s | ~2% |
| I/O, parsing, glue | | 10.96 s | 10.5% |
| **total** | | **104.85 s** | |

Cumulatively **207.9 s -> 104.9 s** across five rounds: the MSA representation, the affine aligner,
`find_most_supported_span`, the segment aligner's per-call cost, and the contexts/MSA trees.

**Three stages now share the top and the rest are small** — `align`, `run_spoa` and
`find_most_supported_span` at 21/20/20%, then a gap to 9%. All three are close to their floor under
the current constraints: `align` is SIMD with a derived band, `run_spoa` is spoa's own algorithm and
its per-call setup measures at zero (see `poa.rs`), and `find_most_supported_span` is 343 000 calls of
already-lean per-call work. Further gains here need a *different algorithm*, not tuning — which for
`run_spoa` means giving up byte-identity with spoa.

Two notes on reading these numbers. Absolute totals move with machine load (the same build measured
132 s and 121 s in two runs while other jobs came and went); *shares* are stable to a point or so and
are what the table is for. And instrumenting inside a stage distorts it badly at these call counts —
a temporary scope on the 7 million per-call alignments added ~10 s of mutex traffic on its own, so
sub-stage measurements have to be taken and then removed, never left in.

### What is not done

- ~~**The accuracy benchmark.**~~ **Done** — see *Accuracy against ground truth*. Both divergences
  *improve* accuracy: the guard by 0.005 pp, the affine MSA aligner by 0.013 pp.
- ~~**Re-profiling.**~~ **Done** — see *Where it stands*. `find_most_supported_span` is now the top
  cost at 28.4%, and it is 343 000 calls rather than one hot loop, so the next win is per-call work.
- **A native POA** to replace `spoars` — still deferred, and *Deferred improvements* now records why
  the case for it is weaker: allocation reuse measures at zero, leaving only the dependency argument.
- ~~**Real long-read data.**~~ **Done** — real ONT drosophila cDNA through isONclust, scored against
  the genome by spliced alignment. See *The drosophila corpus, and the bug it found*; it found one.
  Reads are still not multi-kb (median 553 bp, p99 2 705), so a transcriptome with genuinely long
  transcripts would still stress the guard and the block-width bound harder than anything tested.
- **Decide on `ISONCORRECT_FIX_WIS`.** It is measured, provably optimal and a large accuracy win, and
  it is still off by default. Turning it on — or better, fixing `fill_p2` in the reference — is a
  deliberate decision that re-records every golden and moves every accuracy number in this file.
- **`tools/repo-slim/`** is staged and unrun. It rewrites history and force-pushes, so it needs a
  human at the keyboard.

### The drosophila corpus, and the bug it found

The SIRV corpora are all short (median ~600 bp) and, more importantly, all *one organism's spike-in
transcripts*. Real ONT cDNA from drosophila adds genuine isoform complexity and a longer tail
(p99 2 705 bp, max 7 228 against SIRV's 2 898). 100 000 reads through isONclust gives **266 clusters,
54 721 reads**, the largest 3 528 — comfortably over `--max_seqs`.

**It found a real bug in the port on the first run.** In exact mode, 2 of 266 clusters and 588 of
54 721 reads differed from the reference.

`get_minimizer_combinations_database` drops an anchor pair whose occurrences outnumber the reads
("highly abundant"). The reference's batch loop is `for batch_id, reads in enumerate(batch(all_reads,
args.max_seqs))`, which **rebinds `reads` to the batch**, so `len(reads)` in that filter is the
*batch* size. The port passed the *cluster* size. A higher threshold means abundant anchors the
reference drops survive, and the read picks up extra very-high-support intervals — on the offending
anchor, three extra partners at support 1070 where the reference kept four at 715, 17, 1123 and 55.

Two things make this worth recording beyond the one-line fix:

- **Only a cluster larger than `--max_seqs` can expose it**, and even then only when some anchor's
  occurrence count falls between the batch size and the cluster size. Every SIRV real cluster is
  3 000 reads, so batching *was* exercised there — and the output was byte-identical both before and
  after the fix. The bug was latent on 75 442 reads across two corpora.
- **No stage oracle could see it.** Every one of them passed on the drosophila clusters: 845 810
  edlib CIGARs, 803 268 MSA rows, 800 586 `cigar_to_seq` expansions, 2 682 corrections, 11 526
  alternatives, 3 522 `solve_WIS` calls, 3 522 `fix_correction` calls, 3 522 `correct_read` reads. They
  replay *recorded reference inputs*, so they cannot detect the port following a different trajectory.
  What found it was comparing `spans.tsv` — see below.

#### `spans.tsv` from the live driver closes the last hole

`isoncorrect-dump` can emit spans, but it does not correct, so `previously_corrected_regions` is
always empty there and only the `--exact` trajectory was ever comparable. The anchor filtering that
consumes those regions was therefore **inert in every dump ever taken**, and the non-exact trajectory
had no stage-level check at all.

`driver::SpansDump` records the same format from the real driver, where the feedback loop is live:

```bash
ISONCORRECT_SPANS_DUMP=/tmp/rs_spans.tsv isONcorrect --fastq cluster.fastq --outfolder /tmp/x
diff /tmp/py/spans.tsv /tmp/rs_spans.tsv
```

That diff pointed at row 2 580 of 160 797 — one anchor, on read 1, where the port emitted three
intervals the reference did not. Read 1 is the *first* read, where the regions map is empty and the
filtering is inert, which immediately ruled out the filtering and pointed at the database instead.

After the fix: **1 of 266 clusters and 1 of 54 721 reads** differ, and that one is the parasail
tie-break below.

#### Accuracy against a genome, by spliced alignment

There is no transcriptome here worth depending on, so `bench/accuracy_genome.py` uses the genome:
`minimap2 -ax splice -uf --secondary=no`, error rate `NM / aligned_read_bases` from the primary
alignment.

**It reports the aligned fraction too, and that is not optional.** Error rate alone is gamed by
aligning less of the read — a correction that mangled the ends would soft-clip them away and score
*better*. The pair cannot be gamed in that direction. `NM` over a genome also includes real
biological difference (SNPs against the assembly, mis-placed unannotated splice sites), so the
absolute rate overstates sequencing error; the inflation is identical across implementations, so
comparisons hold where the absolute number does not.

| set | err mean % | err med % | err p90 % | err total % | aligned mean % | unaligned |
| --- | --- | --- | --- | --- | --- | --- |
| raw | 8.445 | 7.407 | 14.472 | 8.507 | 94.03 | 1 853 |
| Python reference | 2.117 | 1.438 | 4.360 | 2.102 | 94.19 | 1 127 |
| Rust, exact mode | 2.117 | 1.438 | 4.360 | 2.102 | 94.19 | 1 127 |
| Rust, default | **2.092** | **1.416** | **4.314** | **2.076** | **94.31** | **1 120** |

Paired against the reference: exact mode differs on **1 read of 52 848**; default is -0.0244 pp,
better on 11 417 and worse on 9 901.

**That is a 53.6% win rate, against 63% on SIRV** — the divergences are still net better on every
statistic, and correction still makes more reads alignable (1 853 unaligned to 1 120), but on real
drosophila the margin is much closer to a coin flip than the simulated and spike-in corpora implied.
Any claim about the size of the divergence's benefit should cite this number, not the SIRV one.

Wall clock, 266 clusters at `--t 8`: Python 192 s, Rust exact ~40 s, **Rust default 19.2 s — 10x**.

#### The parasail tie-break the corpus finally reached

One divergence survives in exact mode, and it is the one this file predicted. *The structural-
overcorrection guard* records that two of the four tie-break parameters — which of the last row and
column to scan first, and whether to keep the first or last equally-good end cell — were **unpinned
by the corpus**, because a read against its own correction never produces an end-cell tie. Drosophila
produces them: 1 of 420 CIGARs on one cluster, 20 of 3 522 on another, with **0 score mismatches** in
both. Purely which equally-optimal path parasail reports.

Re-running the sweep on the new cases **confirms the current configuration rather than overturning
it** — it is still the best at 419/420, and no combination in the 96-parameter space reaches 420. So
parasail's choice here is not expressible as a fixed preference order over predecessors, and
reproducing it would mean reproducing its vectorised scan. Not worth it for one read in 54 721, and
the divergence is *better* on that read by the accuracy measure.

Note also that of the 20 differing CIGARs on cluster 0, **none changed a corrected read**:
`fix_correction` collapsed all of them to the same sequence, which is the point of recording guard
output differences separately from CIGAR differences.

### `solve_WIS` was suboptimal, and fixing it is the largest accuracy win in the project

*Known bugs* #2 recorded from early in the port that `fill_p2` is off by one and that the port
reproduced it, with a note to measure what fixing it does to accuracy first, "since it may be a free
improvement". Measured, then fixed — **in the reference**, because that is where the defect is. It is
**13x larger than both aligner divergences combined**, and the only cost is about 7% in wall clock.

The defect is two defects in one line (`isONcorrect.py:952`):

1. it stores a **0-based** interval index in `p` and then indexes the **1-based** `OPT` array with it,
   so every predecessor is shifted down by one and compatible earlier intervals look incompatible;
2. `j_max` starts at `0`, which means both "interval 0" and "no compatible predecessor", so an
   interval with nothing before it is credited with interval 0's optimum.

**Both are now fixed, in the reference and in the port.** `fill_p2` stores a 1-based index and
reserves `0` for "none"; `wis::fill_p` does the same, `wis::fill_p_reference` keeps the old behaviour
for the test below, and there is no env var — the reference is correct now, so there is nothing to
switch back to. Goldens re-recorded, and fixed Python and fixed Rust are byte-identical on the
fixture.

**It is provably the fix, not merely a different answer.** `wis::fixed_vs_reference` compares both
tables against exhaustive maximum-weight-independent-set search over 3 000 random instances of up to
9 intervals:

| | |
| --- | --- |
| fixed table suboptimal | **0 of 3 000** |
| reference suboptimal | **2 040 of 3 000** |
| reference *better* than the optimum | 0 — impossible, and confirmed |

Accuracy on the 266 real drosophila clusters, 52 849 reads with a primary spliced alignment, with the
fix applied to **the reference** as well as the port:

| set | err mean % | err med % | err p90 % | err total % | aligned mean % |
| --- | --- | --- | --- | --- | --- |
| raw | 8.445 | 7.407 | 14.472 | 8.507 | 94.03 |
| Python, before the fix | 2.117 | 1.438 | 4.360 | 2.102 | 94.19 |
| **Python, after the fix** | **1.800** | **1.188** | **3.714** | **1.755** | 94.21 |
| Rust, exact mode | 1.800 | 1.188 | 3.714 | 1.755 | 94.21 |
| Rust, default (both aligner divergences) | **1.773** | **1.164** | **3.659** | **1.728** | **94.35** |

Paired per read against the fixed reference, the *unfixed* reference is **+0.3162 pp — worse on
26 166 reads and better on 8 065**, a 76% loss rate and a **15% relative increase in error**. Every
summary statistic moves the same way, including the aligned fraction, so this is not trading coverage
for accuracy. For scale, the guard's aligner is worth -0.005 pp and the affine MSA aligner -0.013 pp.

Rust exact mode against fixed Python: **2 reads of 52 849**, both the parasail tie-break, both
*better*. Confirmed by the oracles on the two clusters involved — 0 score mismatches, 2 CIGAR
mismatches each, and every other stage clean (78 288 edlib CIGARs, 76 540 MSA rows, 838 corrections,
838 POA consensuses, 1 130 alternatives).

**Confirmed independently on the SIRV corpus**, which matters because it is a different organism, a
different truth source (edlib against a transcriptome rather than spliced alignment against a genome)
and a different error profile — 32 real isONclust clusters, 75 156 reads. Measured before the fix
landed, comparing the port with and without it: mean 1.696% to 1.632%, total 1.232% to 1.157%,
perfect reads 158 to 169, **-0.0638 pp, better on 20 198 and worse on 9 753** — a 67% win rate.
Smaller in absolute terms than on drosophila (-0.316 pp), which is consistent with everything else in
this file: the spike-in and simulated corpora understate how much these decisions matter on real data.

Re-measured after the fix landed in the reference, comparing the reference against itself on the same
75 156 reads:

| set | mean % | median % | p90 % | total % | perfect |
| --- | --- | --- | --- | --- | --- |
| Python, before the fix | 1.714 | 1.165 | 3.601 | 1.259 | 156 |
| **Python, after the fix** | **1.649** | **1.104** | **3.448** | **1.183** | **163** |
| Rust default | 1.632 | 1.086 | 3.456 | 1.157 | 169 |

The unfixed reference is +0.0652 pp: worse on 20 384 reads and better on 10 119, a 67% loss rate —
the same direction and the same ratio as drosophila, on a different organism and a different truth
source. **And the fixed port is byte-identical to the fixed reference here** — all 32 clusters,
75 442 reads, in exact mode. The two drosophila reads that still differ are the parasail tie-break,
which this corpus does not reach.

It costs about 7% in wall clock (48.1 s to 51.8 s on the 32 SIRV clusters at `--t 8`), which is
exactly what one expects from correcting more regions — the defect's only virtue was doing less work.

This was a defect in the published tool rather than a porting artefact, which is why it is fixed in
`isONcorrect.py` and not merely in the port. Two consequences worth stating:

- **every golden was re-recorded**, and the two unit tests that pinned the old selection
  (`[2, 0]` for three disjoint intervals, `[3, 0]` for the five-interval case) now assert the optimal
  answers `[2, 1, 0]` and `[3, 1]`. The five-interval case is worth reading: weights are 120, 210, 75,
  320, 200, and `{1, 3}` at 530 beats `{0, 3}` at 440, which the old table could not see because it
  did not know 1 precedes 3;
- **the accuracy figures elsewhere in this file predate the fix.** The tables under *Accuracy against
  ground truth* and *Is affine actually more accurate?* compare implementations that all shared this
  defect, so their *relative* conclusions stand — the aligner divergences are still worth what they
  were worth — but every absolute error rate there is now pessimistic by roughly the 0.3 pp above.

### The POA bake-off: `spoars` wins, and that closes the native-POA question

`run_spoa` is 20% of runtime, so it is the obvious next target, and the *reason* it was previously
off-limits no longer holds. This file used to say of the banded POA crates: "implement **different**
algorithms and will not reproduce spoa's output. Do not use them for the identical-output port."
That prohibition was correct then and is void now — byte-identity with spoa is not the specification,
accuracy is. So the alternatives were actually tried.

Every candidate on crates.io, on real recorded intervals (drosophila cluster 0: 2 682 intervals,
287 830 sequences, 17.6 Mbases, mean 107 sequences per interval at 61 bp), configured with
isONcorrect's own scoring:

| candidate | time | agrees with spoa | verdict |
| --- | --- | --- | --- |
| **`spoars` (incumbent, pure Rust)** | **2.83 s** | **2 682 / 2 682** | keep |
| `abpoa-rs` -> abPOA (C, adaptive banded) | — | — | **aborts the process** |
| `poa-consensus` (pure Rust, banded) | 32.7 s | 1 654 / 2 682 (62%) | **11.6x slower** |
| `poasta` (pure Rust, optimal gap-affine) | — | — | **no consensus API** |

Three findings, none of which was predictable from the crate descriptions:

- **abPOA calls `exit(1)` from inside the C library** on real isONcorrect intervals —
  `[simd_abpoa_align_sequence_to_subgraph] Error in lg_backtrack.` and the process is gone. It does
  this in **both** local and global mode. For something invoked once per correction interval, with
  no way to catch it, that is disqualifying regardless of speed. It would also have reintroduced a C
  toolchain, ending the port's freedom from one;
- **`poa-consensus` is an order of magnitude slower**, not faster, despite being banded. Confirmed on
  a second corpus (1 735 intervals: 831 ms against 7 548 ms). It is also a different consensus by
  construction — it picks the *median-length* read as its backbone and ignores the caller's seed, so
  insertion order, which `get_best_corrections` controls deliberately, stops mattering;
- **`poasta` has no consensus function at all.** It is an aligner; heaviest-bundle traversal would
  have to be written on top.

So `spoars` is simultaneously the fastest option available and the only one byte-identical to spoa.
`run_spoa`'s 20% is therefore **not addressable by swapping libraries**, and writing a native POA
would mean beating both spoa and abPOA at their own game — a research project, not a porting task.
The dependency risk (v0.1.3, low usage) stands, and `poa::oracle` is the mitigation: it has held on
**5 705 recorded intervals in a single session** across the fixture, SIRV and drosophila, and must
stay green on every upgrade.

The earlier measurement in `poa.rs` — that hoisting the engine into a thread-local is worth nothing —
fits the same picture: the time is in the DP, and the DP is already vectorised. `spoars` does have
real NEON kernels with spoa's own int16/int32 escalation ladder, so this is *not* the `triple_accel`
situation where SIMD was silently x86-only.

### The real-data corpus

`isONclust` is **not** in the reference environment; it is installed in a separate `isonclust-tool`
conda env, deliberately, so that adding it cannot shift the pinned edlib/parasail/numpy versions the
goldens depend on. `pip install isONclust` fails on osx-arm64 (parasail-python has no wheel), so:
`conda create -n isonclust-tool python=3.10`, `conda install -c bioconda -c conda-forge
parasail-python pysam`, `pip install --no-deps isONclust`.

```bash
isONclust --fastq SIRV_real_full.fastq --outfolder /tmp/isonclust_full --t 8 --ont
isONclust write_fastq --clusters /tmp/isonclust_full/final_clusters.tsv \
  --fastq SIRV_real_full.fastq --outfolder /tmp/real_clusters_full --N 50
```

1.3 M reads -> 490 clusters, 32 with >= 50 reads, largest 555 403. Read lengths 71-2 898, median
~609 — much more variable than the simulated set, which is exactly why it caught the banding
failure the simulated corpus missed. Clustering the full 1.8 GB file takes about ten minutes.

## Profiling

Where the time goes is an open question — it has not been measured. Before optimising anything,
profile the Python reference on a realistic cluster and attribute cost across at least: minimizer
extraction, minimizer-pair database construction, `find_most_supported_span`, `solve_WIS`, the spoa
subprocess, edlib realignment, MSA/PFM construction, and the final parasail guard. Record the
numbers in `bench/` so later claims about bottlenecks are checkable rather than folklore.

Re-profile the Rust port on the same corpus. Bottlenecks will move, and the second profile is the
one that should drive optimisation work.

## Deferred improvements

Ideas discovered while porting that are **out of scope for the identical-output milestone** go
here rather than into the code. Anything that changes output must not ride along with the port.
Keep each entry short: what, why it is deferred, and where it was noticed.

### Known bugs in the reference

These are **defects in the Python implementation**, not porting notes. Fixing them changes output
or behaviour, so they must not ride along with the port — each needs its own commit, and the
goldens must be re-recorded afterwards.

Eight found so far. "Reproduced" means the port deliberately behaves the same way, because the
behaviour reaches the output:

| # | Defect | Affects | Port's stance |
| --- | --- | --- | --- |
| 1 | `get_alternative_ref_contexts` returned a `set`, so corrected regions depended on `PYTHONHASHSEED` | corrected output | **fixed in the reference** (271be1b); output unchanged on all 43 goldens |
| 1b | `create_multialignment_format_NEW`'s `start`/`stop` were always constants, and `correct_read` realigned after `fix_correction` and discarded the result | nothing — both dead | **deleted from the reference**; all 43 goldens byte-identical |
| 2 | `solve_WIS`'s `fill_p2` was off by one, so compatible intervals looked incompatible | corrected fewer regions than intended | **fixed in the reference**; the largest accuracy win in the port — see below |
| 3 | `split_cluster_in_batches` latches its size test after the first small cluster | which files `--split_wrt_batches` creates | reproduced |
| 4 | Joining batched clusters concatenates in lexicographic order (batch 10 before batch 2) | read order in the joined fastq | reproduced |
| 5 | `--split_wrt_batches` raises `ValueError` on any non-numeric filename in the input folder | crashes on real isONclust folders | **diverged**: skipped instead |
| 6 | `--disable_numpy` raises `ValueError` on any eligible context — the flag has never worked | the flag is unusable | flag removed from the Rust CLI |
| 7 | `batch()` raises `UnboundLocalError` on an empty fastq | error path | **diverged**: writes an empty output |
| 8 | `args.flnc` / `args.ccs` are read but never defined by argparse, so no-`--fastq` raises | error path | **diverged**: prints help, exits 0 |

**And one bug in the port, for symmetry**, since it is the only one a corpus has ever found rather
than inspection: `anchors::build` was given the cluster size where the reference passes the batch
size, so the high-abundance anchor filter used the wrong threshold on any cluster over `--max_seqs`.
Fixed; see *The drosophila corpus, and the bug it found*. The lesson is not the bug but that **no
stage oracle could see it** — they replay recorded reference inputs, so a port on a different
trajectory still passes. Comparing `spans.tsv` from the live driver is what caught it.

- ~~**`get_best_corrections` is non-deterministic on the default path.**~~ **Fixed.**
  `get_alternative_ref_contexts` returned a `set`, and the correction loop's `break` and strict `<`
  made the answer depend on its iteration order, hence on `PYTHONHASHSEED` — 51 of 267 corrected
  regions on one real interval took 2–5 different values across 24 seeds. It now returns a list in
  insertion order. All 43 golden outputs are byte-identical after the change, so this was a defect
  removed at no cost in output. See *The reference used to disagree with itself*;
  `bench/check_seed_sensitivity.py` is the regression check.
- **`split_cluster_in_batches` latches its size test.** `smaller_than_max_seqs` is set once a
  cluster fits in one batch and then never re-measured, so every *later* cluster is symlinked whole
  instead of being split. It only works because isONclust numbers clusters largest-first, making
  sizes descend in numeric-id order. A folder not in that order silently under-splits. The port
  reproduces it (`fanout::split_clusters`); fixing it changes which files exist and therefore the
  output layout.
- **`--split_wrt_batches` crashes on any non-numeric filename.** The sort key
  `int(x.split('.')[0])` is applied to *every* directory entry before the `.fastq` filter, so
  isONclust's own `clusters.tsv` — or any stray file — raises `ValueError`. Confirmed by running it:

  ```text
  ValueError: invalid literal for int() with base 10: 'clusters'
  ```

  This matters practically: the reference cannot run `--split_wrt_batches` on a directory that has
  anything but numerically-named fastqs in it. **Intentional divergence:** the port sorts such
  entries last and skips them, so the same command succeeds. Error path only, and the alternative
  would be reproducing a crash.
- **Joining batched clusters concatenates in lexicographic order.**
  `sorted(glob.glob(cl_id + '_*'))` puts batch 10 before batch 2, so a cluster split into ten or
  more batches comes back with its reads in the order 0, 1, 10, 11, 2, … The file is still a valid
  fastq and holds every read, but the order is neither the input's nor numeric. Reproduced, since it
  is the bytes of the output file.
- **`batch()` crashes on an empty fastq.** The loop variable `i` is used after the loop
  (`if i/size != 0`), so a zero-read input raises `UnboundLocalError` rather than writing an empty
  output. Error path only. The port writes an empty output file instead.
- **`--disable_numpy` is broken and always has been.** `test_numba` (numpy path) builds `FCM`
  entries as 3-tuples `(variant, context, depth)`; `sep_function` (non-numpy path) builds a
  `defaultdict(int)` keyed by 2-tuples. The unpack at `src/isoncorrect/isONcorrect.py:693` then
  raises `ValueError: not enough values to unpack (expected 3, got 2)` on any position with an
  eligible context. The flag is documented as "1.5x slower" but it does not run at all. Marked
  `xfail` in `bench/equivalence.sh`. **There is no reference behaviour for the port to match here** —
  decide whether the Rust port rejects the flag, accepts it as a no-op, or waits for a fix.
- ~~**`solve_WIS` is suboptimal: `fill_p2` is off by one.**~~ **Fixed, in the reference and the port.**
  It stored 0-based interval indices in `p` while `p[j]` indexed the **1-based** `OPT` array, and
  `j_max` started at `0`, which meant both "interval 0" and "no compatible predecessor". Three
  pairwise-disjoint intervals came back as two (`[2, 0]` rather than `[2, 1, 0]`).

  The measurement said to fix it: **-0.319 pp mean read error on real drosophila, a 15% relative
  reduction and 13x both aligner divergences combined**, for about 7% in wall clock. See
  *`solve_WIS` is suboptimal, and fixing it is the largest accuracy win found so far*. Goldens
  re-recorded; fixed Python and fixed Rust are byte-identical.

- **`args.flnc` and `args.ccs` are referenced but never defined by argparse**
  (`src/isoncorrect/isONcorrect.py:1679`). Running `isONcorrect` with no `--fastq` raises
  `AttributeError` instead of printing help. Error path only; no effect on correct runs.
  **Intentional divergence:** the Rust port prints help and exits 0 there, matching the evident
  intent and its own no-args behaviour. Reproducing a traceback is not worth it, and no correct run
  reaches this path.

### Performance and structure



- ~~**Write a native linear-gap kSW POA + heaviest-bundle consensus**, replacing `spoars`.~~
  **Closed, by measurement rather than by deferral.** See *The POA bake-off* below.
- ~~**The guard's semi-global alignment is `O(n*m)` in time and memory.**~~ **Done both ways.**
  `block-aligner` replaced it (SIMD + banded, ~0.8% of reads change), and the exact DP behind
  `ISONCORRECT_EXACT_GUARD=1` got the two-row rewrite. Original note follows for context:
- ~~**Memory.**~~
  `parasail.rs` keeps two rows of `H`/`E`/`F` plus one decision byte per cell instead of three
  `i32` matrices: 440 MB peak down to 147 MB, and ~13% faster. Banding was investigated and
  **rejected as unprovable** — see *The guard: 3x less memory, and why it is not banded*. The
  remaining time win is vectorising the DP, which is sound (the recurrence does not change) but a
  substantially larger piece of work; `parasail_cases.tsv` is what makes it checkable. It is still
  ~65% of runtime, so it is the top of the list.
- ~~**The reference realigns after `fix_correction` and discards the result.**~~ **Done** — deleted.
  It was output-neutral by inspection (the return value fed nothing), and all 43 goldens are
  byte-identical after the deletion, which is the check that matters.
- ~~**The positioned vector allocates one `Vec<u8>` per slot.**~~ **Done** — `msa::Positioned` is one
  flat buffer plus an offset table, halving the stage. See *The MSA matrix: one buffer instead of
  `2t + 1` allocations*.
- ~~**`create_multialignment_format_NEW` takes `start`/`stop` that are always constants.**~~ **Done** —
  the only caller passed `0` and `2*len(repr_seq)`, and `position_query_to_alignment` always returns
  `t_vector_start = 0`, `t_vector_end = 2*t` (it asserts as much), so the filter admitted every row
  and the slice was the whole vector. Parameters and filter deleted; all 43 goldens byte-identical.
- The `hash_fcn` parameter is threaded through but hardcoded to `"lex"`. `get_kmer_minimizers`
  (random) and `get_kmer_maximizers` were dead in the default path and are now deleted; the
  parameter itself remains, threaded and constant.
- ~~**Dropping the flags listed under *Scope* leaves a large amount of unreachable Python.**~~
  **Done.** The six out-of-scope flags are gone from the Python CLI too, along with the code only
  they reached: 15 functions in `isONcorrect.py` (`randstrobe*`, `seq_to_strobes_iter*`,
  `get_randstrobes_*`, `get_kmer_minimizers`/`maximizers`,
  `get_minimizers_and_positions_compressed`, `sep_function*`, `get_primes`, `argmin`) and 10 in
  `create_augmented_reference.py` (`run_racon`, five unused spoa variants, `longest_path_botond`,
  `kmer_counter`). **3 286 lines to 2 634, a 20% cut, with all 43 goldens byte-identical.**

  This makes the reference's supported set equal the port's, which is the point — previously the two
  disagreed about which flags exist. It does mean the *reference* now rejects `--randstrobes` and
  friends as unrecognised arguments, where the Rust port names them specifically, so the port's error
  message is the more helpful of the two. The 7 unsupported-flag equivalence cases test the Rust
  binary only and are unaffected.

## Repo hygiene

The repo is ~2.4 GB to clone. **2.920 GB of the 2.935 GB of blob content in history is `data/` and
large files under `test_data/`**; all source and scripts together are ~15 MB.

`tools/repo-slim/` implements the fix in three staged steps (archive → analyze → rewrite), none of
which push or upload on their own. Measured result: `.git` drops to ~900 KB, a full working clone
to ~2.9 MB. See `tools/repo-slim/README.md`, which also documents the two traps that make this
harder than it looks (blob deduplication hiding paths, and `git rev-parse` printing to stdout on
failure).

- **Never commit data files.** No `.csv`, `.fa`, `.fasta`, `.fastq`, `.tar.bz2`, `.bam`, or
  `.part*` beyond the small fixtures under `test_data/`. The root `.gitignore` now enforces this.
- Large fixtures and paper data belong in a GitHub Release or Zenodo, fetched by a script.
- `test_data/chr6_ensemble.fa` (18 MB) is unused by any test and is stripped. Its blob is
  byte-identical to `data/chr6_transcripts.fa`, which is exactly how it evaded a first attempt at
  building the removal list.

## Working agreements

- **Never run `git commit`, `git push`, or any other history-changing git command. Ever.** Stage
  nothing, commit nothing, rewrite nothing. Leave finished work in the working tree, say what
  changed and why, and let a human commit it. This holds even when the work is complete, verified,
  and obviously correct, and even when a note below says a change "needs its own commit" — those
  notes describe how commits should be *organised* when a human makes them, not permission to make
  one. If a commit seems warranted, propose the message and stop.
- Don't touch `paper/` or `scripts/` — those reproduce
  published figures.
- Don't "improve" the algorithm while porting. Behaviour changes and the port must not land in the
  same commit; an intentional divergence needs its own commit and a note here.
- **"Behaviour" means observable output, not internal representation.** Data structures are free:
  different containers, dropping provably-dead entries, avoiding an accidental `defaultdict`
  insertion, arena allocation — none of that is a behaviour change if the emitted bytes are
  identical and the harness stays green. Reproducing the reference's *data* is the contract;
  reproducing its *memory layout* is not, and doing so would forfeit most of the memory win the
  port exists to deliver. What is **not** free: anything that feeds output — iteration order where
  it reaches results, tie-breaking, the order supporting reads are collected, arithmetic and
  rounding.
- **When you spot a possible improvement, write it into *Deferred improvements* and move on.**
  That section is the backlog — losing these ideas is the failure mode, and acting on them
  mid-port is the other one.
- Keep the Python implementation in the tree until the Rust port passes the full equivalence sweep.
- Rust: `cargo fmt`, `cargo clippy -- -D warnings`, and `cargo test` clean before any commit.
