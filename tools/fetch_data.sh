#!/usr/bin/env bash
#
# Restore the paper data that was removed from git history.
#
# data/ holds the result CSVs and figure inputs behind the isONcorrect paper. It
# is ~1.1 GB and is no longer committed, because carrying it made the repository
# 2.4 GB to clone. It is archived on Zenodo instead, with a DOI.
#
# Nothing in the repository imports data/ at run time -- it is consumed by hand
# by the plotting scripts under scripts/ and evaluation*/. You only need this if
# you are regenerating paper figures.
#
#   tools/fetch_data.sh                 # download and extract into data/
#   tools/fetch_data.sh --check         # verify an existing data/ against checksums
#   tools/fetch_data.sh --record 12345  # override the Zenodo record id
#
set -euo pipefail

# Zenodo record for the isONcorrect paper data.
# Version-specific record id, so this always fetches the same bytes.
ZENODO_RECORD="${ZENODO_RECORD:-21920617}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$REPO_ROOT"
CHECK_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --record) ZENODO_RECORD="$2"; shift 2 ;;
    --dest)   DEST="$2"; shift 2 ;;
    --check)  CHECK_ONLY=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$ZENODO_RECORD" ]]; then
  cat >&2 <<'EOF'
error: no Zenodo record configured.

The paper data has not been archived yet, or this script has not been updated
with the record id. Set it at the top of this file, or pass one:

  tools/fetch_data.sh --record 1234567

To create the archive in the first place, see tools/repo-slim/README.md.
EOF
  exit 1
fi

API="https://zenodo.org/api/records/${ZENODO_RECORD}"

echo "==> querying Zenodo record ${ZENODO_RECORD}"
META="$(mktemp)"
trap 'rm -f "$META"' EXIT
if ! curl -fsSL "$API" -o "$META"; then
  echo "error: could not fetch $API" >&2
  exit 1
fi

python3 - "$META" <<'PY'
import json, sys
meta = json.load(open(sys.argv[1]))
print(f"    title:   {meta.get('title','?')}")
doi = meta.get('doi') or meta.get('metadata', {}).get('doi')
print(f"    doi:     {doi}")
files = meta.get('files', [])
total = sum(f.get('size', 0) for f in files)
print(f"    files:   {len(files)}  ({total/1e9:.2f} GB)")
PY

# A plain file rather than `mapfile`: macOS ships bash 3.2, which has neither
# mapfile nor readarray.
ENTRIES="$(mktemp)"
trap 'rm -f "$META" "$ENTRIES"' EXIT

python3 - "$META" > "$ENTRIES" <<'PY'
import json, sys
meta = json.load(open(sys.argv[1]))
for f in meta.get('files', []):
    links = f.get('links', {})
    link = links.get('self') or links.get('download') or links.get('content')
    key = f.get('key') or f.get('filename')
    chk = (f.get('checksum') or '').replace('md5:', '')
    if key and link:
        print(f"{key}\t{link}\t{chk}")
PY

if [[ ! -s "$ENTRIES" ]]; then
  echo "error: record has no downloadable files" >&2
  exit 1
fi

DL="$DEST/.zenodo-download"
mkdir -p "$DL"

md5_of() {
  if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi
}

MISSING=0
while IFS=$'\t' read -r key link chk; do
  [[ -z "$key" ]] && continue
  target="$DL/$key"

  if [[ -f "$target" && -n "$chk" ]] && [[ "$(md5_of "$target")" == "$chk" ]]; then
    echo "    have   $key (checksum ok)"
    continue
  fi

  if [[ $CHECK_ONLY -eq 1 ]]; then
    echo "    MISSING or CHANGED: $key"
    MISSING=1
    continue
  fi

  echo "==> downloading $key"
  # Zenodo file keys may contain slashes (they mirror the uploaded path), so the
  # parent directory is not guaranteed to exist.
  mkdir -p "$(dirname "$target")"
  curl -fL --progress-bar "$link" -o "$target"

  if [[ -n "$chk" ]]; then
    got="$(md5_of "$target")"
    if [[ "$got" != "$chk" ]]; then
      echo "error: checksum mismatch for $key" >&2
      echo "  expected $chk" >&2
      echo "  got      $got" >&2
      exit 1
    fi
    echo "    checksum ok"
  fi
done < "$ENTRIES"

if [[ $CHECK_ONLY -eq 1 ]]; then
  echo "==> check complete"
  exit $MISSING
fi

echo "==> extracting"
while IFS=$'\t' read -r key _ _; do
  [[ -z "$key" ]] && continue
  case "$key" in
    *.tar.gz|*.tgz) tar -xzf "$DL/$key" -C "$DEST" && echo "    extracted $key" ;;
    *.tar.bz2)      tar -xjf "$DL/$key" -C "$DEST" && echo "    extracted $key" ;;
    *.tar)          tar -xf  "$DL/$key" -C "$DEST" && echo "    extracted $key" ;;
    MANIFEST.txt)   cp "$DL/$key" "$DEST/data-MANIFEST.txt" 2>/dev/null || true ;;
    *)              echo "    (leaving $key in $DL)" ;;
  esac
done < "$ENTRIES"

echo
echo "==> done. data/ restored under $DEST"
echo "    downloads cached in $DL (safe to delete)"
