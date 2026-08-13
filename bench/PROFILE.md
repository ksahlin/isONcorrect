# Profiling record — Python reference

Numbers produced by `bench/profile_python.py`. Recorded so later claims about bottlenecks are
checkable rather than folklore. **Append to this file; don't overwrite it.** Each entry must say
what was profiled and on what, because the answer changes with both.

---

## 2026-08-13 — first measurement

**Caveat that matters more than the numbers:** this is a single 100-read cluster
(`test_data/isoncorrect/0.fastq`). It is a smoke test. The stages that are cheap here are exactly
the ones expected to grow with cluster size, so **do not treat this as the profile for real data.**
Re-run on clusters of 1000+ reads from real isONclust output before optimising anything.

- Machine: Apple Silicon (arm64), macOS
- Env: python 3.12.13, numpy 2.1.3, parasail 1.3.4, edlib 1.3.9.post1, **spoa 4.1.5**
- Input: `test_data/isoncorrect/0.fastq` (100 reads)

### Default parameters (`--k 9 --w 20 --max_seqs 2000`)

Total wall: 0.78 s

| stage | calls | self s | self % |
| --- | ---: | ---: | ---: |
| spoa subprocess | 106 | 0.54 | **68.9** |
| `get_alternative_ref_contexts` | 106 | 0.12 | 14.7 |
| `get_best_corrections` | 106 | 0.03 | 4.0 |
| `get_contexts` | 106 | 0.02 | 1.9 |
| `create_multialignment_matrix` | 106 | 0.01 | 1.8 |
| `find_most_supported_span` | 2 533 | 0.01 | 1.7 |
| `get_minimizer_combinations_database` | 1 | 0.01 | 1.7 |
| `parasail_alignment` (structural guard) | 200 | 0.01 | 1.6 |
| `edlib.align` | 3 506 | 0.01 | 1.3 |
| everything else | — | <0.01 | <1 |

### Paper parameters (`--k 9 --w 10 --max_seqs 1000`)

Total wall: 1.99 s

| stage | calls | self s | self % |
| --- | ---: | ---: | ---: |
| spoa subprocess | 168 | 1.01 | **50.8** |
| `get_minimizer_combinations_database` | 1 | 0.26 | **13.0** |
| `get_alternative_ref_contexts` | 168 | 0.19 | 9.6 |
| `find_most_supported_span` | 9 629 | 0.13 | 6.4 |
| `get_best_corrections` | 168 | 0.08 | 4.2 |
| `edlib.align` | 48 584 | 0.08 | 3.8 |
| `create_multialignment_matrix` | 168 | 0.03 | 1.6 |
| `get_contexts` | 168 | 0.02 | 1.1 |
| `parasail_alignment` (structural guard) | 200 | 0.01 | 0.7 |

### What this does and does not show

- **spoa is the largest single line item on this input, but it is not the whole story and its
  share is not stable.** Halving `w` (default → paper) drops spoa from 69% to 51% while
  `get_minimizer_combinations_database` goes from 1.7% to 13% and `edlib.align` call count goes
  from 3.5k to 48.6k. The profile is a function of the parameters, not a property of the program.
- ~5 ms per spoa call for alignments this small is almost entirely `fork`/`exec` plus temp-file
  I/O, not alignment. That overhead disappears in any in-process implementation regardless of how
  fast the POA itself is.
- `get_minimizer_combinations_database` is a **single call** already costing 13% at paper settings
  on 100 reads. It builds an index over minimizer pairs shared across reads, so it is the stage
  most likely to dominate both time and memory as cluster size grows. It is the top candidate to
  overtake spoa on realistic input, and it is untested at scale here.
- `get_alternative_ref_contexts` at 10–15% is consistently third and is pure Python list/dict
  churn — likely a large constant-factor win in Rust.
- The structural-overcorrection guard (`parasail_alignment` + `fix_correction`) is cheap: 200 calls
  for 100 reads, under 2%. Note the 2:1 ratio — `correct_read` aligns twice per read and discards
  the second result (see *Deferred improvements* in `PORTING.md`).

### Open questions for the next profiling round

1. How do these shares move at 1 000 / 2 000 / 5 000 reads per cluster?
2. Is `get_minimizer_combinations_database` memory-bound? Peak RSS is not captured by this tool —
   use `bench/benchmark.sh` alongside it.
3. Does `find_most_supported_span` grow super-linearly with cluster size? Its call count already
   scales with anchor density.
