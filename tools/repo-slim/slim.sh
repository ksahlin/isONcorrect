#!/usr/bin/env bash
#
# Step 2 of repo slimming: rewrite history in a throwaway clone and verify it.
#
# This script NEVER touches your working repository and NEVER pushes. It produces
# a rewritten mirror in a scratch directory, checks it, and prints the commands
# you would run to publish it. Pushing is your call and yours alone -- it is the
# irreversible part.
#
#   tools/repo-slim/analyze.sh      # first: produce and review removal-paths.txt
#   tools/repo-slim/slim.sh         # then: rewrite + verify
#
# Options:
#   --workdir DIR   where to build the rewritten mirror (default: a temp dir)
#   --keep-workdir  don't delete the workdir on failure
#
# Requires git-filter-repo:
#   pipx install git-filter-repo      (or)   pip install git-filter-repo
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOVAL="$HERE/removal-paths.txt"

WORKDIR=""
KEEP_WORKDIR=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)      WORKDIR="$2"; shift 2 ;;
    --keep-workdir) KEEP_WORKDIR=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# --- preconditions -----------------------------------------------------------

if [[ ! -f "$REMOVAL" ]]; then
  echo "error: $REMOVAL not found. Run tools/repo-slim/analyze.sh first." >&2
  exit 1
fi

FILTER_REPO="${GIT_FILTER_REPO:-}"
if [[ -z "$FILTER_REPO" ]]; then
  if command -v git-filter-repo >/dev/null 2>&1; then
    FILTER_REPO="$(command -v git-filter-repo)"
  else
    echo "error: git-filter-repo not found." >&2
    echo "  install with:  pipx install git-filter-repo" >&2
    echo "            or:  pip install git-filter-repo" >&2
    echo "  or set GIT_FILTER_REPO=/path/to/git-filter-repo" >&2
    exit 1
  fi
fi

ORIGIN="$(git -C "$REPO_ROOT" remote get-url origin)"
echo "==> source repo:  $REPO_ROOT"
echo "==> origin:       $ORIGIN"
echo "==> filter-repo:  $FILTER_REPO"
echo "==> removing:     $(grep -vc '^#' "$REMOVAL") paths"
echo

# --- build a throwaway mirror ------------------------------------------------

CLEANUP_WORKDIR=0
if [[ -z "$WORKDIR" ]]; then
  WORKDIR="$(mktemp -d)"
  CLEANUP_WORKDIR=1
fi
mkdir -p "$WORKDIR"
MIRROR="$WORKDIR/isONcorrect-slim.git"

cleanup() {
  if [[ $? -ne 0 && $KEEP_WORKDIR -eq 0 && $CLEANUP_WORKDIR -eq 1 ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

echo "==> cloning a fresh mirror into $MIRROR"
echo "    (filter-repo requires a pristine clone; your working repo is untouched)"
rm -rf "$MIRROR"
git clone --mirror --no-local "$REPO_ROOT" "$MIRROR" 2>&1 | sed 's/^/    /'

size_of() { du -sk "$1" | awk '{printf "%.0f", $1/1024}'; }
BEFORE_MB="$(size_of "$MIRROR")"
echo "==> size before: ${BEFORE_MB} MB"

# --- rewrite -----------------------------------------------------------------

echo
echo "==> rewriting history"
# --invert-paths: the file lists what to REMOVE.
# --force: the mirror is freshly made by us, so filter-repo's freshness guard is
#          satisfied but the mirror carries origin refs it wants confirmation on.
(
  cd "$MIRROR"
  "$FILTER_REPO" --paths-from-file "$REMOVAL" --invert-paths --force 2>&1 | sed 's/^/    /'
  git reflog expire --expire=now --all
  git gc --prune=now --aggressive --quiet
)

AFTER_MB="$(size_of "$MIRROR")"
echo "==> size after:  ${AFTER_MB} MB"

# --- verify ------------------------------------------------------------------

echo
echo "==> verifying"
FAIL=0
check() { # check <description> <expected> <actual>
  if [[ "$2" == "$3" ]]; then
    printf '    ok    %-52s %s\n' "$1" "$3"
  else
    printf '    FAIL  %-52s expected %s, got %s\n' "$1" "$2" "$3"
    FAIL=1
  fi
}

# blob_at <repo> <path> -> blob sha, or the literal "absent". `git rev-parse`
# prints the unresolved string on stdout when a path is missing, so it cannot be
# used here; cat-file -e is the reliable existence test.
blob_at() {
  if git -C "$1" cat-file -e "HEAD:$2" 2>/dev/null; then
    git -C "$1" rev-parse "HEAD:$2"
  else
    echo absent
  fi
}

# filter-repo drops commits that become empty once their only content is
# stripped, so the total count legitimately falls. What must NOT change is the
# number of commits touching source.
ORIG_COMMITS="$(git -C "$REPO_ROOT" rev-list --count --all)"
NEW_COMMITS="$(git -C "$MIRROR" rev-list --count --all)"
printf '    info  %-52s %s -> %s (%s data-only commits pruned)\n' \
  "commit count" "$ORIG_COMMITS" "$NEW_COMMITS" "$((ORIG_COMMITS - NEW_COMMITS))"

ORIG_SRC="$(git -C "$REPO_ROOT" rev-list --count --all -- src scripts evaluation README.md setup.py pyproject.toml)"
NEW_SRC="$(git -C "$MIRROR" rev-list --count --all -- src scripts evaluation README.md setup.py pyproject.toml)"
check "commits touching source preserved" "$ORIG_SRC" "$NEW_SRC"

# Source must survive untouched.
for f in src/isoncorrect/isONcorrect.py src/isoncorrect/correct_seqs.py \
         src/isoncorrect/run_isoncorrect.py setup.py pyproject.toml README.md; do
  check "unchanged: $f" "$(blob_at "$REPO_ROOT" "$f")" "$(blob_at "$MIRROR" "$f")"
done

# Test fixtures must survive.
for f in test_data/isoncorrect/0.fastq test_data/isoncorrect/1.fastq \
         test_data/100reads.fq test_data/six_exons.fa \
         test_data/sirv_transcriptome.fasta; do
  check "fixture kept: $f" "$(blob_at "$REPO_ROOT" "$f")" "$(blob_at "$MIRROR" "$f")"
done

# Stripped paths must be gone from HEAD.
for f in data/chr6_transcripts.fa test_data/chr6_ensemble.fa data/12_results.csv; do
  check "stripped from HEAD: $f" "absent" "$(blob_at "$MIRROR" "$f")"
done

# ...and from all of history. Walk the trees, not the object listing, for the
# same dedup reason that made the removal list wrong when built from objects.
KEEP_RE='^test_data/(isoncorrect/[01]\.fastq|100reads\.fq|six_exons\.fa|sirv_transcriptome.*\.fasta)$'
LEFTOVER="$(git -C "$MIRROR" log --all --pretty=format: --name-only --no-renames \
            | sed '/^$/d' | LC_ALL=C sort -u \
            | grep -E '^(data/|test_data/)' \
            | grep -vE "$KEEP_RE" || true)"
if [[ -z "$LEFTOVER" ]]; then
  printf '    ok    %-52s %s\n' "no stripped paths anywhere in history" "clean"
else
  printf '    FAIL  %-52s %s remain\n' "stripped paths still in history" "$(wc -l <<<"$LEFTOVER" | tr -d ' ')"
  head -5 <<<"$LEFTOVER" | sed 's/^/          /'
  FAIL=1
fi

# The whole point.
SHRINK="$(python3 -c "print(f'{$BEFORE_MB/max($AFTER_MB,1):.0f}x')")"
printf '    info  %-52s %s MB -> %s MB (%s smaller)\n' "repository size" "$BEFORE_MB" "$AFTER_MB" "$SHRINK"

echo
if [[ $FAIL -ne 0 ]]; then
  echo "VERIFICATION FAILED — do not push this. Workdir kept at:" >&2
  echo "  $MIRROR" >&2
  KEEP_WORKDIR=1
  exit 1
fi
echo "==> verification passed"

# --- hand over ---------------------------------------------------------------

cat <<EOF

================================================================================
Rewritten mirror is ready. NOTHING HAS BEEN PUSHED.

  $MIRROR

Before you push, understand what this does:

  * Every commit SHA changes. All 9 forks diverge permanently and cannot be
    fast-forwarded. Anyone with a clone must re-clone.
  * Commit-pinned links break, including any in the paper or in Zenodo records.
  * GitHub keeps the old objects reachable for a while; stripped data may remain
    downloadable via old SHAs until GitHub garbage-collects. Ask GitHub Support
    to run gc if the data is sensitive. It is public paper data here, so this is
    a tidiness issue rather than a disclosure one.

Recommended order:

  1. Archive the data first, so it is not lost:
       tools/repo-slim/archive_data.sh

  2. Tag the current state so the pre-rewrite history stays reachable:
       git -C "$REPO_ROOT" tag pre-slim-\$(date +%Y%m%d)
       git -C "$REPO_ROOT" push origin pre-slim-\$(date +%Y%m%d)

  3. Push the rewritten history:
       git -C "$MIRROR" push --force --mirror origin

  4. Re-clone fresh and confirm:
       git clone $ORIGIN /tmp/isONcorrect-verify
       du -sh /tmp/isONcorrect-verify

  5. Re-run the equivalence harness against the fresh clone to confirm the
     fixtures still work:
       cd /tmp/isONcorrect-verify && bench/equivalence.sh record

Step 3 is the irreversible one.
================================================================================
EOF

CLEANUP_WORKDIR=0
