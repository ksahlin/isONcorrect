# isONcorrect — Rust rewrite

## Goal

Port isONcorrect from Python to Rust. **Identical CLI, identical output**, primarily faster and
lower memory. The Python implementation in `src/isoncorrect/` is the **normative reference**: when
Rust and Python disagree, Python is right until a human decides otherwise.

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
| semi-global alignment (parasail `sg_trace_scan_16`) | done natively, no C++; tie-break measured, 200 real alignments byte-identical |
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
algorithms and will not reproduce spoa's output. Do not use them for the identical-output port.

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

## First end-to-end performance measurement, and it is below target

Now that `isONcorrect` runs, the port can be timed. On one 200-read SIRV gene cluster
(~1.4 kb reads, default parameters, single core):

| | wall | user | sys | peak RSS |
| --- | --- | --- | --- | --- |
| Python reference | 16.1 s | 9.0 s | 6.5 s | 339 MB |
| Rust port | 10.7 s | 10.6 s | 0.1 s | 420 MB |

**1.5x faster, and 24% *more* memory.** That is well short of the goal at the top of this file, and
the shape of the numbers says exactly why:

- **The win is all in `sys` time.** Python spends 6.5 s in the kernel, almost entirely
  `fork`/`exec` plus temp-file I/O for one `spoa` subprocess *per correction interval*. The port's
  sys time is 0.1 s. Removing that was the easy, structural win.
- **The port loses on `user` time** — 10.6 s against 9.0 s. The reference delegates its inner loops
  to C: edlib is bit-parallel Myers, spoa and parasail are SIMD. The port answers with plain
  `O(n*m)` scalar DP in `align.rs` and `parasail.rs`. Correctness first was the right order, but the
  arithmetic has to get faster before the port is worth switching to.
- **Memory is worse for the same reason.** `parasail.rs` allocates three `i32` tables over the full
  read against its own correction — ~24 MB per read at 1.4 kb — and `align.rs` allocates a table per
  segment alignment. Both are transient, but they set the peak.

Now that `run_isoncorrect` works, the same question can be asked of a realistic workload: 68 SIRV
transcript clusters (20–342 reads each, ~1.4 kb), `--t 8`, both implementations parallel:

| | wall |
| --- | --- |
| Python `run_isoncorrect --t 8` | 83.7 s |
| Rust `run_isoncorrect --t 8` | 84.1 s |

**Dead even** — and all 68 `corrected_reads.fastq` byte-identical. That result is worth stating
plainly, because it corrects the assumption that motivated porting the batch driver first: the
reference is *already* parallel via `multiprocessing`, so removing the subprocess and interpreter
startup buys back only what per-cluster inefficiency gives away. Per-cluster speed is the lever, and
it multiplies through the whole pipeline.

A profile of the port on one cluster says where that speed is (5 064 samples, 200-read gene
cluster):

| component | share |
| --- | --- |
| `find_most_supported_span` → **bounded edit distance** | **67%** |
| the parasail guard | 19% |
| `get_best_corrections` in total (spoa 2%, NW CIGAR 3%, MSA <1%) | ~10% |
| everything else (minimizers, anchors, WIS, I/O) | <1% |

And the reason the top entry is so large is worth knowing before optimising anything else:
**`triple_accel` ships SIMD for `x86`/`x86_64` only** — 67 `target_arch` guards for those two, none
for `aarch64` — so on Apple Silicon or any ARM host every call falls back to naive scalar DP. Edit
distance is *uniquely defined*, so replacing it carries no tie-break risk at all, and
`bench/gen_editdist_cases.py` already provides the oracle. That makes it both the largest and the
safest of the remaining wins.

None of this is new information about *what* to do; every item is already in *Deferred improvements*
(band or linear-space the guard alignment, replace the scalar NW, drop the per-slot `Vec<u8>` in the
MSA). What is new is that they are now measured rather than suspected, and that the equivalence
oracles make each one safe to attempt. Do them in profile order, not in the order they were noticed.

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
| 2 | `solve_WIS`'s `fill_p2` is off by one, so compatible intervals look incompatible | corrects fewer regions than intended | reproduced; conservative, but worth measuring before fixing |
| 3 | `split_cluster_in_batches` latches its size test after the first small cluster | which files `--split_wrt_batches` creates | reproduced |
| 4 | Joining batched clusters concatenates in lexicographic order (batch 10 before batch 2) | read order in the joined fastq | reproduced |
| 5 | `--split_wrt_batches` raises `ValueError` on any non-numeric filename in the input folder | crashes on real isONclust folders | **diverged**: skipped instead |
| 6 | `--disable_numpy` raises `ValueError` on any eligible context — the flag has never worked | the flag is unusable | flag removed from the Rust CLI |
| 7 | `batch()` raises `UnboundLocalError` on an empty fastq | error path | **diverged**: writes an empty output |
| 8 | `args.flnc` / `args.ccs` are read but never defined by argparse, so no-`--fastq` raises | error path | **diverged**: prints help, exits 0 |

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
- **`solve_WIS` is suboptimal: `fill_p2` is off by one.** It stores 0-based interval indices in
  `p`, but `p[j]` is then used to index the **1-based** `OPT` array (`OPT[p[j]]`). Every predecessor
  is therefore shifted, and compatible earlier intervals get treated as incompatible. Three
  pairwise-disjoint intervals, which a correct WIS would all select, come back as two — verified
  directly against `solve_WIS`, which returns `[2, 0]` rather than `[2, 1, 0]`.

  The effect is to correct over **fewer** intervals than intended, so it is conservative rather than
  wrong-in-a-dangerous-way, but it does mean isONcorrect leaves correctable regions untouched.
  **The port reproduces it**, because fixing it changes corrected output on real data. If it is ever
  fixed, that is a behaviour change needing its own commit, a note here, and re-recorded goldens —
  and it is worth measuring what it does to accuracy first, since it may be a free improvement.

- **`args.flnc` and `args.ccs` are referenced but never defined by argparse**
  (`src/isoncorrect/isONcorrect.py:1679`). Running `isONcorrect` with no `--fastq` raises
  `AttributeError` instead of printing help. Error path only; no effect on correct runs.
  **Intentional divergence:** the Rust port prints help and exits 0 there, matching the evident
  intent and its own no-args behaviour. Reproducing a traceback is not worth it, and no correct run
  reaches this path.

### Performance and structure

- **Write a native linear-gap kSW POA + heaviest-bundle consensus**, replacing `spoars`. Removes the
  only low-usage dependency on the critical path and allows reusing graph allocations across
  intervals, which matters because POA runs once per correction interval. `poa.rs` already isolates
  the surface (`consensus`), and `poa::oracle` is the acceptance test — 505 real intervals against
  the `spoa` binary. Deferred because `spoars` is already byte-identical on that corpus.
- **The guard's semi-global alignment is `O(n*m)` in time and memory**, on the full read against its
  own correction — roughly 2 M cells and 24 MB of scratch for a 1 400 bp read, once per read.
  `parasail.rs` is a plain three-matrix Gotoh; parasail itself is striped SIMD. Two ways out, both
  output-neutral and both checkable against `parasail_cases.tsv`: keep only two rows plus a
  traceback-direction byte per cell, which drops the 24 MB to ~2 MB, or band the DP, since the two
  sequences differ by very little. Deferred because it changes nothing observable and the oracle
  makes it safe to do later. Profile first — the guard runs once per read, while spoa runs once per
  correction interval.
- **The reference realigns after `fix_correction` and discards the result**
  (`isONcorrect.py:1312`). The port omits that second call. It is output-neutral by inspection — the
  return value feeds nothing — and the `correct_read` oracle agrees on 1 204 reads, but deleting it
  from the Python is a separate commit.
- **The positioned vector allocates one `Vec<u8>` per slot**, i.e. `2 * len(consensus) + 1` small
  allocations per supporting segment per correction interval — the shape the reference's list of
  strings has. Nearly all of them are the one-byte `"-"`. A slot is either a single character or a
  contiguous range of the query, so an enum of `(char)` / `(start, len)` would remove essentially
  all of it. Deferred because it is an inner-loop rewrite with no output change, and `msa.rs` has
  the oracle that makes it safe to do later.
- `create_multialignment_format_NEW` takes `start`/`stop` and filters to rows covering that window,
  but the only caller passes the full vector, so the filter always admits everything and
  `t_vector_start`/`t_vector_end` are constants. The port drops the parameters. Deleting them from
  the Python is a separate commit.
- The `hash_fcn` parameter is threaded through but hardcoded to `"lex"`; `get_kmer_minimizers`
  (random) and `get_kmer_maximizers` are dead in the default path.
- Dropping the flags listed under *Scope* leaves a large amount of now-unreachable Python
  (`randstrobe*`, `seq_to_strobes_iter*`, `run_racon`, `get_minimizers_and_positions_compressed`,
  and most of `create_augmented_reference`). Deleting it from the Python tree is tempting but is a
  separate commit, and only after the port has passed the full sweep.

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
