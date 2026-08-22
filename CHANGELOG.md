Changelog
=========

## v0.2.0 — the Rust rewrite

isONcorrect is now a Rust program. The command-line interface is unchanged — same flags, same
defaults, same output format — so existing pipelines do not need editing. The Python implementation
is still in the repository and still runs, but it is **deprecated**: it is kept as the reference
implementation that the Rust version is verified against.

Install with `cd rust && cargo build --release`. The only build dependency is a Rust toolchain
(1.85+); there is no C or C++ toolchain, no CMake and no vendored native source.

### Corrected output has changed

**`solve_WIS` was leaving correctable regions uncorrected, and that is fixed** in both
implementations. Its predecessor table stored 0-based indices into a 1-based array and used `0` to
mean both "the first interval" and "no previous interval", so compatible regions looked incompatible.
Verified against exhaustive maximum-weight-independent-set search: the fixed table is optimal on 3000
random instances, the old one suboptimal on 2040 and never better.

On 52 849 real drosophila reads, mean per-read error goes from **2.117% to 1.800%** — a 15% relative
reduction, improving 26 166 reads and worsening 8 065. Confirmed independently on 75 156 real SIRV
reads (1.714% to 1.649%). Every summary statistic moves the same way, including the fraction of each
read that aligns, so this is not trading coverage for accuracy.

Anyone comparing against output from an earlier version should expect it to differ, and to be better.

Two further changes are deliberate, each measured to be faster *and* more accurate on real data:

| Where | Instead of | Effect | Restore with |
| --- | --- | --- | --- |
| the structural-overcorrection guard | parasail semi-global | SIMD banded; ~0.8% of reads change | `ISONCORRECT_EXACT_GUARD=1` |
| segments against the consensus | edlib, unit-cost edit distance | SIMD affine; ~11% of alignments change | `ISONCORRECT_EXACT_ALIGN=1` |

With both variables set, output is byte-identical to the fixed Python reference on 75 442 real SIRV
reads and all 43 recorded test cases. On 54 721 real drosophila reads it differs on 2, where two
equally-optimal alignments exist and the two implementations report different ones.

### Faster and smaller

Real ONT data, 8 threads, against the Python implementation:

| | Python | Rust | |
| --- | --- | --- | --- |
| 32 SIRV clusters, 75 442 reads | 432 s | **43.5 s** | 9.9x |
| 266 drosophila clusters, 54 721 reads | 181 s | **18.5 s** | 9.8x |
| one 200-read cluster, peak memory | 339 MB | **~147 MB** | 2.3x less |

### Removed flags

Six flags are gone from both implementations: `--randstrobes`, `--layers`, `--set_layers_manually`,
`--use_racon`, `--disable_numpy`, `--compression`. They were experimental, deprecated, unreachable
from `run_isoncorrect`, or — in the case of `--disable_numpy` — broken: it raised `ValueError` on any
input that reached it, so it had never worked. The Rust binaries name the specific flag and exit
non-zero; the Python ones now reject them as unrecognised arguments.

### Repository

The paper result data (2.9 GB) has been removed from git history and archived on Zenodo at
<https://zenodo.org/records/21920617>; `tools/fetch_data.sh` restores and verifies it. A full clone is
now **3.2 MB** rather than ~2.4 GB. **Every commit SHA before the rewrite has changed**, so
commit-pinned links will not resolve; the `v0.1.0` and `v0.1.3.5` tags were rewritten with the rest
and still point into the current history.

### Also

- Two unused parameters and a discarded realignment call removed from the Python reference, along with
  25 functions left unreachable by the flag removals — 3 286 lines to 2 634.
- A non-determinism defect fixed in the reference: `get_alternative_ref_contexts` returned a `set`, so
  corrected regions depended on `PYTHONHASHSEED`.
- `bench/accuracy_genome.py` scores accuracy against a genome by spliced alignment, for organisms with
  no transcriptome worth depending on.
- CI now runs on GitHub Actions (the old Travis badge had been dead for years).
