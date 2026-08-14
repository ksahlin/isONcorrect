# bench/ — equivalence and performance harness

The Rust port's acceptance criterion is **byte-identical output** against the Python reference.
This directory is how that gets checked, plus the profiling and timing that tell us whether the
port is actually worth having.

## Quick start

```bash
bench/setup_reference_env.sh                  # pinned conda env with spoa + parasail + edlib
bench/equivalence.sh record                   # run Python, store goldens
bench/equivalence.sh verify                   # run Rust, diff against goldens
```

Before the Rust binaries exist, `record` and the profiler are still useful on their own.

## Files

| File | Purpose |
| --- | --- |
| `env/reference.yml` | Pinned reference environment. spoa decides every consensus — treat version bumps as breaking. |
| `env/resolved-*.txt` | Exact versions that got resolved, written at setup. Makes a golden run auditable. |
| `setup_reference_env.sh` | Creates the env, installs isONcorrect editable, verifies, records versions. |
| `equivalence.sh` | `record` / `verify` / `both` / `list`. The acceptance gate. |
| `benchmark.sh` | Wall clock and peak RSS, Python vs Rust, with a thread sweep. |
| `profile_python.py` | Stage attribution for the reference: self vs inclusive time per stage. |
| `corpus/` | Test clusters. Gitignored — never commit fastq here. |
| `golden/` | Recorded reference outputs. Gitignored. |
| `results/` | Benchmark CSVs and logs. Gitignored. |

## Corpus

`bench/corpus/` is empty by default and the harness falls back to `test_data/isoncorrect`
(two 100-read clusters, which are in fact **one** cluster under two names — see `PORTING.md`).
**That fallback is a smoke test, not proof of equivalence** — the scripts say so on every run.

A real corpus needs clusters that exercise the paths the small fixtures don't:

- large clusters (>2000 reads) so the `--max_seqs` batching path runs for real
- clusters above the `--exact_instance_limit` of 50, since at or below it `isoncorrect_main`
  forces `--exact` and the `previously_corrected_regions` filtering never runs
- clusters with genuine exon variation, which is what the structural-overcorrection guard exists for
- low-coverage clusters where interval support is thin
- clusters with repeated anchors within a read, the only thing that triggers edlib in
  `find_most_supported_span`

### Building one from simulated SIRV reads

`make_sirv_corpus.py` splits a simulated fastq into per-cluster files by the source transcript in
the read header (`@read_1_from_SIRV612`), so it needs no aligner and no clustering run — and unlike
an isONclust pass it is reproducible from the fastq alone, which is what a verification corpus
wants.

```bash
bench/make_sirv_corpus.py --fastq reads_10k_err7%.fastq --outdir $CORPUS/gene       --group gene
bench/make_sirv_corpus.py --fastq reads_10k_err7%.fastq --outdir $CORPUS/transcript --group transcript
```

The two groupings stress different things, and both are worth having:

| Grouping | On the 10k 7%-error SIRV set | What it exercises |
| --- | --- | --- |
| `gene` | 7 clusters, 713–2242 reads | isoforms of one gene in one cluster, so reads differ by whole exons; the largest crosses `--max_seqs`, so batching runs |
| `transcript` | 68 clusters, 20–342 reads | one transcript each: deep coverage where every difference is sequencing error. 62 of the 68 are above the `--exact` limit |

Output is isONclust's layout (`0.fastq`, `1.fastq`, … largest first, plus `clusters.tsv`), so
`run_isoncorrect --fastq_folder` reads it directly.

Write it **outside the repo** — anything committed to git is a mistake, and `.gitignore` here is
deliberately aggressive.

Real ONT cDNA through isONclust is still worth adding: simulated reads have uniform error and no
chimeras, so they under-represent exactly the messy cases the guard exists for.

## Interpreting the profiler

```bash
conda activate isoncorrect-ref
bench/profile_python.py --fastq test_data/isoncorrect/0.fastq --outfolder /tmp/prof
```

Two columns matter:

- **self** — time in that stage excluding other instrumented stages nested inside it. This is the
  column that tells you what to optimise. It sums to at most the total.
- **inclusive** — total time under that stage. Double-counts nesting, useful for seeing that e.g.
  `get_best_corrections` dominates while its own self time is small because the cost is in spoa
  and edlib beneath it.

The gap between "accounted for" and total wall time is fastq parsing, output writing, interpreter
startup, and glue.

Record the numbers here so later claims about bottlenecks can be checked rather than argued.

## Case types

`bench/equivalence.sh list` shows the matrix. Cases carry an optional note:

- **(no note)** — the normal case. Python output is recorded as golden; Rust must match it byte for
  byte.
- **`unsupported:<flag>`** — a flag deliberately not ported (see *Scope* in `PORTING.md`). No golden
  is recorded. `verify` instead asserts that the Rust binary **exits non-zero and names the flag**
  in its message. Exiting 0, or failing with a message that doesn't mention the flag, is a failure —
  a pipeline that passes `--use_racon` and gets a clean exit with no polishing is worse off than one
  that stops.

  This covers both scope tiers, because the observable contract is identical.
  `--disable_numpy` and `--compression` are removed from the Rust CLI entirely and satisfy it via
  the argument parser's unknown-flag error; `--randstrobes` and `--use_racon` parse and are
  rejected with a specific message. The harness does not care which, only that the binary refuses
  and says why.
- **`xfail:<reason>`** — the *reference* cannot run this at all, so there is nothing to compare.
  Currently unused.

## Rules

- A non-empty diff in `verify` is a failure. There is no tolerance threshold.
- Goldens are only comparable to the environment that produced them; `PROVENANCE.txt` is written
  alongside them for exactly this reason. Re-record after any dependency change.
- `PYTHONHASHSEED=0` is exported by the harness so the `--randstrobes` path is stable run to run.
  It is *not* stable in ordinary Python usage — see the determinism section in `PORTING.md`.
- Performance numbers are meaningless without a passing equivalence run alongside them.
