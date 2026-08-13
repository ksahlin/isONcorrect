#!/usr/bin/env bash
#
# Step 0 of repo slimming: get the data out of git before it is stripped.
#
# Packages data/ and the large test_data/ files into tarballs suitable for a
# GitHub Release, and writes a manifest with checksums so the archive can be
# verified later.
#
# By DEFAULT this only builds archives locally and prints what it would upload.
# Uploading publishes ~1.1 GB under your account and is not something this script
# does on its own -- pass --upload once you have looked at what was built.
#
#   tools/repo-slim/archive_data.sh                 # build locally, report
#   tools/repo-slim/archive_data.sh --upload        # build, then upload to a Release
#
# Options:
#   --outdir DIR   where to build archives (default: ./repo-slim-archive)
#   --tag TAG      release tag (default: data-archive-YYYYMMDD)
#   --upload       actually create the Release and upload
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUTDIR="$REPO_ROOT/repo-slim-archive"
TAG="data-archive-$(date +%Y%m%d)"
UPLOAD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tag)    TAG="$2"; shift 2 ;;
    --upload) UPLOAD=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

cd "$REPO_ROOT"

if [[ ! -d data ]]; then
  echo "error: no data/ directory in the working tree." >&2
  echo "       Run this BEFORE slimming, on a checkout that still has the data." >&2
  exit 1
fi

mkdir -p "$OUTDIR"
echo "==> building archives in $OUTDIR"

# GitHub Releases cap individual assets at 2 GB. data/ is ~1.1 GB, so one
# tarball is fine, but split anyway to keep individual downloads manageable and
# resumable on poor connections.
DATA_TAR="$OUTDIR/isONcorrect-paper-data.tar"

echo "    packing data/ ..."
tar -cf "$DATA_TAR" data
echo "    compressing (this takes a while on ~1.1 GB) ..."
if command -v pigz >/dev/null 2>&1; then
  pigz -f "$DATA_TAR"
else
  gzip -f "$DATA_TAR"
fi
DATA_TGZ="${DATA_TAR}.gz"

# The one large test_data file, kept separate: it is a reference genome slice,
# not paper results, and someone may want it without the 1.1 GB.
EXTRA=()
if [[ -f test_data/chr6_ensemble.fa ]]; then
  echo "    packing test_data/chr6_ensemble.fa ..."
  tar -czf "$OUTDIR/chr6_ensemble.fa.tar.gz" test_data/chr6_ensemble.fa
  EXTRA+=("$OUTDIR/chr6_ensemble.fa.tar.gz")
fi

echo "==> writing manifest"
MANIFEST="$OUTDIR/MANIFEST.txt"
{
  echo "isONcorrect data archive"
  echo "Created: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "From commit: $(git rev-parse HEAD)"
  echo
  echo "This archive holds the paper data and large fixtures that were removed"
  echo "from git history to make the repository cloneable. The analysis code that"
  echo "consumes it lives in evaluation*/ and scripts/ in the repository."
  echo
  echo "## Contents"
  for f in "$DATA_TGZ" "${EXTRA[@]:-}"; do
    [[ -f "$f" ]] || continue
    echo "  $(basename "$f")  $(du -h "$f" | cut -f1)"
  done
  echo
  echo "## SHA256"
  for f in "$DATA_TGZ" "${EXTRA[@]:-}"; do
    [[ -f "$f" ]] || continue
    if command -v shasum >/dev/null 2>&1; then
      shasum -a 256 "$f" | sed "s#$OUTDIR/##"
    else
      sha256sum "$f" | sed "s#$OUTDIR/##"
    fi
  done
  echo
  echo "## Original contents of data/"
  git ls-tree -r --name-only HEAD -- data | sed 's/^/  /'
} > "$MANIFEST"

echo
cat "$MANIFEST"
echo

ASSETS=("$DATA_TGZ" "${EXTRA[@]:-}" "$MANIFEST")

if [[ $UPLOAD -eq 0 ]]; then
  cat <<EOF
================================================================================
Archives built locally. NOTHING HAS BEEN UPLOADED.

Review the files in $OUTDIR (~$(du -sh "$OUTDIR" | cut -f1)), then pick a home.

RECOMMENDED: Zenodo. This is data behind a published paper, so it should be
citable and archival, which a GitHub Release is not.

  1. https://zenodo.org  ->  New upload
  2. Drag in every file from $OUTDIR
  3. Upload type "Dataset"; title e.g.
       "isONcorrect: paper result data and large fixtures"
  4. Link it to the publication DOI (10.1038/s41467-020-20340-8) under
     "Related/alternate identifiers"
  5. Publish, then note the record id from the URL
     (https://zenodo.org/records/<RECORD_ID>)
  6. Put that id into tools/fetch_data.sh, and cite the DOI in README.md

  Zenodo allows 50 GB per record and versions immutably.

ALTERNATIVE: a GitHub Release, if you want it beside the code and do not need a
DOI. Re-run with --upload, or by hand:

  gh release create "$TAG" \\
    --title "Paper data and large fixtures" \\
    --notes "Data removed from git history to make the repository cloneable. See MANIFEST.txt." \\
$(printf '    %s \\\n' "${ASSETS[@]}" | sed '$ s/ \\$//')

Either way: confirm the files download cleanly BEFORE running slim.sh. Once
history is rewritten, this working tree is the only remaining copy.
================================================================================
EOF
  exit 0
fi

if ! command -v gh >/dev/null 2>&1; then
  echo "error: gh CLI not found; cannot upload." >&2
  exit 1
fi

echo "==> creating release $TAG and uploading $(du -sh "$OUTDIR" | cut -f1)"
gh release create "$TAG" \
  --title "Paper data and large fixtures" \
  --notes "Data removed from git history to make the repository cloneable. See MANIFEST.txt for contents and checksums." \
  "${ASSETS[@]}"

echo
echo "==> done. Verify the assets are downloadable BEFORE running slim.sh:"
echo "    gh release view $TAG"
