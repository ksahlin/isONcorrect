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
| `find_most_supported_span`, WIS | next |
| consensus / POA, MSA, PFM | not started |
| structural-overcorrection guard | not started |

Argument names, defaults, validation order, stderr text and exit codes match the reference. The
correction algorithm is not ported: both binaries validate their arguments and then exit non-zero
saying so.

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

**The dump tool must apply the same argument massaging `main` does.** A mismatch at
`--xmin 14 --k 9` turned out to be the dump binary, not the anchor logic: `main` clamps `--xmin` up
to `2*k` before anything runs, so the reference was building spans at 18 while the port used 14.
Anything added to `main`'s preamble has to be mirrored in `dump_stages.rs` or the comparison lies.

Current `bench/equivalence.sh verify`: **7 passed, 22 failed** — the 7 are the unsupported-flag
cases, and the 22 failures are the supported cases with no algorithm behind them yet. That is the
expected state; it becomes the progress bar as the algorithm lands.

Verified equal to the reference by hand: `--version` on both binaries, exit codes for
`--w 150` (1), `--split_mod 2 --residual 5` (1), `--version` (0), no-args (0), and the exact stderr
text of `xmin set to 40` and the window-range message. Help *layout* differs — clap leads with the
description, argparse with usage — which is cosmetic and outside the contract.

```bash
cargo build --release --manifest-path rust/Cargo.toml
cargo test --manifest-path rust/Cargo.toml
```

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
   (`match=4, mismatch=-8, gap_open=24, gap_ext=1` via `sg_trace_scan_16`, falling back to `_32`),
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

Candidate approaches, in order of preference:

1. **Reimplement linear-gap kSW POA + heaviest-bundle consensus natively in Rust.** Small, no FFI,
   no C++ toolchain in the build, and it is the only option that also removes the process spawn.
   Must be validated to be **bit-identical** against the `spoa` binary over a large corpus of real
   intervals (see `bench/`).
2. **`spoars` crate** — advertises itself as a faithful native-Rust reimplementation of spoa. If it
   is genuinely faithful this is option 1 for free. It is very new and low-usage; treat "faithful"
   as a claim to verify, not a fact.
3. **`spoa` / `spoa-sys` crate** (bindings to the C++ library). Guaranteed-identical output and
   removes the fork+temp-file overhead, at the cost of a C++ dependency in the build.

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

Performance is measured on the same corpus: wall clock and peak RSS, Rust vs Python, single-core
and at `--t N`.

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

- **`--disable_numpy` is broken and always has been.** `test_numba` (numpy path) builds `FCM`
  entries as 3-tuples `(variant, context, depth)`; `sep_function` (non-numpy path) builds a
  `defaultdict(int)` keyed by 2-tuples. The unpack at `src/isoncorrect/isONcorrect.py:693` then
  raises `ValueError: not enough values to unpack (expected 3, got 2)` on any position with an
  eligible context. The flag is documented as "1.5x slower" but it does not run at all. Marked
  `xfail` in `bench/equivalence.sh`. **There is no reference behaviour for the port to match here** —
  decide whether the Rust port rejects the flag, accepts it as a no-op, or waits for a fix.
- **`args.flnc` and `args.ccs` are referenced but never defined by argparse**
  (`src/isoncorrect/isONcorrect.py:1679`). Running `isONcorrect` with no `--fastq` raises
  `AttributeError` instead of printing help. Error path only; no effect on correct runs.
  **Intentional divergence:** the Rust port prints help and exits 0 there, matching the evident
  intent and its own no-args behaviour. Reproducing a traceback is not worth it, and no correct run
  reaches this path.

### Performance and structure

- `correct_read` recomputes the parasail alignment after `fix_correction`
  (`src/isoncorrect/isONcorrect.py:1312`) and discards the result — it feeds nothing. Removing it
  is output-neutral, but confirm that against the equivalence harness before dropping it.
- `run_isoncorrect` shells out to `python isONcorrect.py` per batch via `subprocess`, with four
  near-duplicate `check_call` branches for the flag combinations. In Rust this should be in-process
  work distribution, which also removes the per-batch interpreter startup.
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
