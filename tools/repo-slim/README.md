# tools/repo-slim — shrinking the repository

The repo is ~2.4 GB to clone. **99.5% of that is committed data**: `data/` and large files under
`test_data/`. All source, scripts, and evaluation code together are ~15 MB.

Measured on the current history:

| | |
| --- | --- |
| total distinct blob bytes in history | 2.935 GB |
| stripped by these tools | 2.920 GB (99.5%) |
| retained source, scripts, docs | 0.015 GB |
| retained fixtures under `test_data/` | 0.267 MB |
| **rewritten `.git`** | **~900 KB** |
| **full working clone after rewrite** | **~2.9 MB** |

## The three steps

Run them in this order. Each stops before doing anything irreversible.

```bash
tools/repo-slim/archive_data.sh     # 0. get the data out of git, into a Release
tools/repo-slim/analyze.sh          # 1. compute + review the removal list
tools/repo-slim/slim.sh             # 2. rewrite history in a scratch clone, verify
```

Only after all three pass do you force-push, and `slim.sh` prints the exact commands. Nothing here
pushes or uploads on its own.

| Script | Does | Never does |
| --- | --- | --- |
| `archive_data.sh` | Builds tarballs of `data/` + `chr6_ensemble.fa` with a checksummed manifest | Upload, unless you pass `--upload` |
| `analyze.sh` | Writes `removal-paths.txt` and `analysis.txt` | Modify the repository at all |
| `slim.sh` | Clones a scratch mirror, rewrites it, verifies it | Touch your working repo, or push |

`removal-paths.txt` is generated, reviewable, and the single input to the rewrite. Read it before
running `slim.sh`.

## Requirements

`git-filter-repo` — `pipx install git-filter-repo` or `pip install git-filter-repo`. If it is
installed somewhere unusual, point at it with `GIT_FILTER_REPO=/path/to/git-filter-repo`.

## Two traps this tooling exists to avoid

**1. Deduplicated blobs.** `git rev-list --objects` emits each blob exactly *once*, with a single
path. Any file whose content is byte-identical to another file is therefore invisible in that
listing. Building the removal list from it silently missed `test_data/chr6_ensemble.fa`, because
its 18 MB blob is identical to `data/chr6_transcripts.fa` — the rewrite ran, reported success, and
left the file in place. The removal list is now built from a tree walk
(`git log --all --name-only`), which enumerates every path regardless of content sharing.

The same trap applies to `git rev-parse HEAD:<missing-path>`: it prints the unresolved string to
stdout *and* exits non-zero, so `$(git rev-parse ... || echo gone)` yields both. Existence checks
use `git cat-file -e`.

**2. `test_data/isoncorrect/0.fastq` and `1.fastq` are the same file.** They share blob
`64cae773` with each other *and* with `test_data/100reads.fq` — one blob, three paths, identical
MD5. The "two test clusters" used by CI and by `bench/` are **one cluster twice**. Stripping is
path-based so keeping all three is correct, but do not mistake them for independent coverage. This
is why `bench/README.md` insists on a real corpus.

## What the rewrite costs

- **Every commit SHA changes.** All 9 forks diverge permanently and cannot be fast-forwarded;
  anyone holding a clone must re-clone. None of the forks are ahead of `master` today, so the
  practical cost is low — but it is not zero.
- **Commit-pinned links break**, including any in the paper or in Zenodo records.
- **GitHub keeps old objects reachable for a while.** Stripped data may stay downloadable via old
  SHAs until GitHub garbage-collects; ask GitHub Support to run `gc` if that matters. This is
  public paper data, so it is a tidiness question rather than a disclosure one.

Mitigation `slim.sh` recommends: tag the pre-rewrite state and push that tag first, so the old
history stays reachable by name.

## Verification

`slim.sh` refuses to declare success unless all of these hold:

- commits touching source are preserved exactly (the total count legitimately falls — data-only
  commits become empty and are pruned; that number is reported as info)
- named source files have byte-identical blobs before and after
- every retained fixture has a byte-identical blob before and after
- known-stripped paths are absent from `HEAD`
- **no** path under `data/` or `test_data/` survives anywhere in history except the keep list

The end-to-end check worth repeating by hand: clone the rewritten repo and run
`bench/equivalence.sh record` against its fixtures. The output must match the goldens byte for
byte. It did when this tooling was built.
