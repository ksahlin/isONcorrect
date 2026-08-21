#!/usr/bin/env bash
#
# Equivalence harness: the acceptance criterion for the Rust port.
#
# The Python reference is slow, so its output is recorded once as a set of
# goldens and then compared against cheaply and repeatedly. A non-empty diff on
# any corrected_reads.fastq is a failure — there is no "close enough" here.
#
#   bench/equivalence.sh record          # run Python, store goldens
#   bench/equivalence.sh verify          # run Rust, diff against goldens
#   bench/equivalence.sh both            # record then verify
#   bench/equivalence.sh list            # print the case matrix and exit
#
# Options:
#   --corpus DIR    directory of per-cluster .fastq files
#                   (default: bench/corpus if non-empty, else test_data/isoncorrect)
#   --golden DIR    where goldens live (default: bench/golden)
#   --case NAME     run only this case (repeatable)
#   --jobs N        parallel cases (default 1; raise once you trust the harness)
#   --keep-going    do not stop at the first failing case
#
# Environment overrides:
#   ISONCORRECT_REF_ENV   conda env holding the Python reference
#   PY_ISONCORRECT        path to the Python `isONcorrect` binary
#   PY_RUN_ISONCORRECT    path to the Python `run_isoncorrect` binary
#   RS_ISONCORRECT        path to the Rust `isONcorrect` binary
#   RS_RUN_ISONCORRECT    path to the Rust `run_isoncorrect` binary
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_ENV="${ISONCORRECT_REF_ENV:-isoncorrect-ref}"

MODE="${1:-both}"
[[ $# -gt 0 ]] && shift || true

CORPUS=""
GOLDEN="$REPO_ROOT/bench/golden"
ONLY_CASES=()
JOBS=1
KEEP_GOING=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --corpus)     CORPUS="$2"; shift 2 ;;
    --golden)     GOLDEN="$2"; shift 2 ;;
    --case)       ONLY_CASES+=("$2"); shift 2 ;;
    --jobs)       JOBS="$2"; shift 2 ;;
    --keep-going) KEEP_GOING=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# The case matrix.
#
# Format: name|kind|args|note
#   kind=single -> `isONcorrect --fastq <cluster>` once per corpus cluster
#   kind=folder -> `run_isoncorrect --fastq_folder <corpus>` once
#   note        -> optional:
#                  "xfail:<reason>"       the *reference* cannot run this, so there
#                                         is no golden and nothing to verify.
#                  "unsupported:<flag>"   deliberately not ported (see Scope in
#                                         PORTING.md). No golden. Verified by
#                                         asserting the Rust binary exits non-zero
#                                         and names the flag.
#
# Every parameter here changes behaviour. Adding a flag to the CLI means adding
# a case, otherwise the flag is untested and the port is free to get it wrong.
# ---------------------------------------------------------------------------
CASES=(
  "default|single||"
  "paper|single|--k 9 --w 10 --max_seqs 1000|"
  "k11|single|--k 11|"
  "k7|single|--k 7|"
  "w15|single|--w 15|"
  "w_eq_k_plus1|single|--k 9 --w 10|"
  "xspan_narrow|single|--xmin 14 --xmax 40|"
  "xspan_wide|single|--xmin 25 --xmax 120|"
  "T_low|single|--T 0.05|"
  "T_high|single|--T 0.3|"
  # max_seqs below the cluster size forces the multi-batch path, which is where
  # batch-boundary and read-ordering bugs surface.
  "batched|single|--max_seqs 40|"
  "spoa_cap|single|--max_seqs_to_spoa 20|"
  "spoa_cap_tiny|single|--max_seqs_to_spoa 3|"
  "exact|single|--exact|"
  "exact_limit|single|--exact_instance_limit 200|"
  "dynamic_w|single|--set_w_dynamically|"

  # --- Deliberately not ported. See "Scope" in PORTING.md. ---
  #
  # Whichever tier a flag is in, the assertion is the same: non-zero exit and
  # the flag named in the output. Tier 1 flags are absent from the Rust CLI
  # entirely and satisfy this via the argument parser's unknown-flag error;
  # Tier 2 flags parse and are rejected with a specific message.
  #
  # Tier 1 (removed entirely): --disable_numpy crashes in the reference
  # (isONcorrect.py:693 unpacks 3-tuples from a dict keyed by 2-tuples), and
  # --compression is deprecated and not exposed by run_isoncorrect at all.
  # Neither has working users to protect.
  "no_numpy|single|--disable_numpy|unsupported:--disable_numpy"
  "compression|single|--compression|unsupported:--compression"
  # Tier 2 (recognised and rejected): these work today in Python, so they may
  # appear in existing pipeline scripts and deserve a specific message.
  "randstrobes|single|--randstrobes|unsupported:--randstrobes"
  "randstrobes_layers|single|--randstrobes --set_layers_manually --layers 3|unsupported:--randstrobes"
  "racon|single|--use_racon|unsupported:--use_racon"
  "folder_racon|folder|--t 2 --use_racon|unsupported:--use_racon"
  "folder_randstrobes|folder|--t 2 --randstrobes|unsupported:--randstrobes"

  "folder_default|folder|--t 3|"
  "folder_paper|folder|--t 3 --k 9 --w 10 --max_seqs 1000|"
  "folder_split_batches|folder|--t 3 --split_wrt_batches --max_seqs 40|"
  "folder_dynamic_w|folder|--t 3 --set_w_dynamically|"
  "folder_single_core|folder|--t 1|"
  "folder_split_mod|folder|--t 2 --split_mod 2 --residual 0|"
)

if [[ "$MODE" == "list" ]]; then
  printf '%-24s %-8s %-46s %s\n' "CASE" "KIND" "ARGS" "NOTE"
  for c in "${CASES[@]}"; do
    IFS='|' read -r n k a note <<<"$c"
    printf '%-24s %-8s %-46s %s\n' "$n" "$k" "${a:-<defaults>}" "$note"
  done
  exit 0
fi

# ---------------------------------------------------------------------------
# Resolve corpus and binaries
# ---------------------------------------------------------------------------
if [[ -z "$CORPUS" ]]; then
  if compgen -G "$REPO_ROOT/bench/corpus/*.fastq" >/dev/null; then
    CORPUS="$REPO_ROOT/bench/corpus"
  else
    CORPUS="$REPO_ROOT/test_data/isoncorrect"
    echo "note: bench/corpus is empty, falling back to test_data/isoncorrect" >&2
    echo "      (100-read clusters are a smoke test, not proof of equivalence)" >&2
  fi
fi
CORPUS="$(cd "$CORPUS" && pwd)"

conda_prefix_for() {
  conda env list | awk -v n="$1" '$1==n {print $NF; exit}'
}

resolve_python_bins() {
  if [[ -n "${PY_ISONCORRECT:-}" && -n "${PY_RUN_ISONCORRECT:-}" ]]; then
    return
  fi
  local prefix
  prefix="$(conda_prefix_for "$REF_ENV")"
  if [[ -z "$prefix" || ! -x "$prefix/bin/isONcorrect" ]]; then
    echo "error: reference env '$REF_ENV' not found or incomplete." >&2
    echo "       run: bench/setup_reference_env.sh" >&2
    exit 1
  fi
  PY_ISONCORRECT="$prefix/bin/isONcorrect"
  PY_RUN_ISONCORRECT="$prefix/bin/run_isoncorrect"
  # run_isoncorrect re-invokes `python isONcorrect.py`, so the env's bin
  # directory must lead on PATH or it will pick up the wrong interpreter.
  export PATH="$prefix/bin:$PATH"
}

resolve_rust_bins() {
  RS_ISONCORRECT="${RS_ISONCORRECT:-$REPO_ROOT/rust/target/release/isONcorrect}"
  RS_RUN_ISONCORRECT="${RS_RUN_ISONCORRECT:-$REPO_ROOT/rust/target/release/run_isoncorrect}"
  if [[ ! -x "$RS_ISONCORRECT" || ! -x "$RS_RUN_ISONCORRECT" ]]; then
    echo "error: Rust binaries not found." >&2
    echo "       looked for: $RS_ISONCORRECT" >&2
    echo "                   $RS_RUN_ISONCORRECT" >&2
    echo "       build with: cargo build --release --manifest-path rust/Cargo.toml" >&2
    exit 1
  fi
}

# The reference is only deterministic on the default path; --randstrobes depends
# on Python's randomised string hash. Pinning the seed makes even that path
# reproducible run to run, so goldens stay stable. See PORTING.md.
export PYTHONHASHSEED=0

# Two stages of the port default to affine SIMD alignment rather than
# reproducing the reference's aligner: the guard (parasail) and the
# segment-vs-consensus alignment that builds the MSA (edlib, unit cost). Both
# find an optimal-scoring alignment, but of a different objective and with their
# own tie-break, so both change output. Those are the port's two deliberate
# divergences (see PORTING.md); this gate is about everything else, so it runs
# the exact reference-compatible DP at both. Unset these to measure the
# divergence instead of asserting byte-identity.
export ISONCORRECT_EXACT_GUARD=1
export ISONCORRECT_EXACT_ALIGN=1

# ---------------------------------------------------------------------------
# Running one case
# ---------------------------------------------------------------------------

# run_case <impl:py|rs> <name> <kind> <args> <outdir>
run_case() {
  local impl="$1" name="$2" kind="$3" args="$4" outdir="$5"
  local single run_folder log
  if [[ "$impl" == "py" ]]; then
    single="$PY_ISONCORRECT"; run_folder="$PY_RUN_ISONCORRECT"
  else
    single="$RS_ISONCORRECT"; run_folder="$RS_RUN_ISONCORRECT"
  fi

  rm -rf "$outdir"
  mkdir -p "$outdir"
  log="$outdir/.harness.log"

  if [[ "$kind" == "single" ]]; then
    local f stem
    for f in "$CORPUS"/*.fastq; do
      stem="$(basename "$f" .fastq)"
      mkdir -p "$outdir/$stem"
      # shellcheck disable=SC2086  # args are intentionally word-split
      if ! "$single" --fastq "$f" --outfolder "$outdir/$stem" $args >>"$log" 2>&1; then
        echo "    ! $impl/$name failed on $stem (see $log)" >&2
        return 1
      fi
    done
  else
    # shellcheck disable=SC2086
    if ! "$run_folder" --fastq_folder "$CORPUS" --outfolder "$outdir" $args >>"$log" 2>&1; then
      echo "    ! $impl/$name failed (see $log)" >&2
      return 1
    fi
  fi
  return 0
}

# List every corrected_reads.fastq under a dir, as paths relative to that dir,
# sorted so the two sides line up.
list_outputs() {
  ( cd "$1" && find . -name corrected_reads.fastq -type f | LC_ALL=C sort )
}

# ---------------------------------------------------------------------------
# record / verify
# ---------------------------------------------------------------------------

selected_cases() {
  local c n
  for c in "${CASES[@]}"; do
    IFS='|' read -r n _ _ <<<"$c"
    if [[ ${#ONLY_CASES[@]} -eq 0 ]]; then
      echo "$c"
    else
      local want
      for want in "${ONLY_CASES[@]}"; do
        [[ "$want" == "$n" ]] && echo "$c"
      done
    fi
  done
}

do_record() {
  resolve_python_bins
  echo "==> recording goldens"
  echo "    corpus: $CORPUS ($(ls "$CORPUS"/*.fastq 2>/dev/null | wc -l | tr -d ' ') clusters)"
  echo "    golden: $GOLDEN"
  echo "    python: $PY_ISONCORRECT"
  mkdir -p "$GOLDEN"

  # Stamp the goldens with the environment that produced them.
  {
    echo "corpus: $CORPUS"
    echo "recorded_by: $PY_ISONCORRECT"
    for f in "$REPO_ROOT"/bench/env/resolved-*.txt; do
      [[ -e "$f" ]] && { echo "--- $(basename "$f") ---"; sed -n '1,12p' "$f"; }
    done
  } > "$GOLDEN/PROVENANCE.txt"

  local c name kind args note rc=0
  while IFS= read -r c; do
    IFS='|' read -r name kind args note <<<"$c"
    printf '  %-24s ' "$name"
    if [[ "$note" == xfail:* ]]; then
      echo "xfail (${note#xfail:}) — no golden recorded"
      continue
    fi
    if [[ "$note" == unsupported:* ]]; then
      echo "unsupported (${note#unsupported:}) — no golden recorded"
      continue
    fi
    local t0 t1
    t0=$(date +%s)
    if run_case py "$name" "$kind" "$args" "$GOLDEN/$name"; then
      t1=$(date +%s)
      local n
      n=$(list_outputs "$GOLDEN/$name" | wc -l | tr -d ' ')
      echo "ok  ($n outputs, $((t1 - t0))s)"
    else
      echo "FAILED"
      rc=1
      [[ $KEEP_GOING -eq 1 ]] || return 1
    fi
  done < <(selected_cases)
  return $rc
}

do_verify() {
  resolve_rust_bins
  if [[ ! -d "$GOLDEN" ]]; then
    echo "error: no goldens at $GOLDEN — run 'bench/equivalence.sh record' first" >&2
    exit 1
  fi
  echo "==> verifying Rust against goldens"
  echo "    corpus: $CORPUS"
  echo "    golden: $GOLDEN"
  echo "    rust:   $RS_ISONCORRECT"

  local workdir
  workdir="$(mktemp -d)"
  trap 'rm -rf "$workdir"' RETURN

  local c name kind args note pass=0 fail=0 xfail=0
  while IFS= read -r c; do
    IFS='|' read -r name kind args note <<<"$c"
    printf '  %-24s ' "$name"

    if [[ "$note" == xfail:* ]]; then
      echo "xfail (${note#xfail:})"
      xfail=$((xfail + 1))
      continue
    fi

    # Dropped flags: assert the port refuses them clearly instead of quietly
    # doing something else. A zero exit here is a failure.
    if [[ "$note" == unsupported:* ]]; then
      local flag="${note#unsupported:}" out rc
      if [[ "$kind" == "single" ]]; then
        # shellcheck disable=SC2086
        out="$("$RS_ISONCORRECT" --fastq "$CORPUS"/$(ls "$CORPUS" | head -1) \
               --outfolder "$workdir/$name" $args 2>&1)" && rc=0 || rc=$?
      else
        # shellcheck disable=SC2086
        out="$("$RS_RUN_ISONCORRECT" --fastq_folder "$CORPUS" \
               --outfolder "$workdir/$name" $args 2>&1)" && rc=0 || rc=$?
      fi
      if [[ "${rc:-0}" -eq 0 ]]; then
        echo "FAILED (exited 0; must reject $flag)"
        fail=$((fail + 1))
        [[ $KEEP_GOING -eq 1 ]] || { echo; echo "stopping (use --keep-going to continue)"; return 1; }
      elif ! grep -qF -- "$flag" <<<"$out"; then
        echo "FAILED (rejected, but message does not name $flag)"
        echo "$out" | tail -3 | sed 's/^/        /'
        fail=$((fail + 1))
        [[ $KEEP_GOING -eq 1 ]] || { echo; echo "stopping (use --keep-going to continue)"; return 1; }
      else
        echo "ok  (rejected $flag, exit $rc)"
        pass=$((pass + 1))
      fi
      continue
    fi

    if [[ ! -d "$GOLDEN/$name" ]]; then
      echo "SKIP (no golden)"
      continue
    fi

    if ! run_case rs "$name" "$kind" "$args" "$workdir/$name"; then
      echo "FAILED (rust exited non-zero)"
      fail=$((fail + 1))
      [[ $KEEP_GOING -eq 1 ]] || { echo; echo "stopping (use --keep-going to continue)"; return 1; }
      continue
    fi

    # Compare the *set* of outputs first: a missing or extra cluster is just as
    # much a failure as wrong bases, and diffing file-by-file would hide it.
    local g_list r_list
    g_list="$(list_outputs "$GOLDEN/$name")"
    r_list="$(list_outputs "$workdir/$name")"
    if [[ "$g_list" != "$r_list" ]]; then
      echo "FAILED (output set differs)"
      diff <(echo "$g_list") <(echo "$r_list") | sed 's/^/      /' || true
      fail=$((fail + 1))
      [[ $KEEP_GOING -eq 1 ]] || { echo; echo "stopping (use --keep-going to continue)"; return 1; }
      continue
    fi

    local rel bad=0
    while IFS= read -r rel; do
      [[ -z "$rel" ]] && continue
      if ! cmp -s "$GOLDEN/$name/$rel" "$workdir/$name/$rel"; then
        if [[ $bad -eq 0 ]]; then echo "FAILED"; fi
        bad=$((bad + 1))
        echo "      differs: $rel"
        if [[ $bad -le 2 ]]; then
          diff "$GOLDEN/$name/$rel" "$workdir/$name/$rel" | head -20 | sed 's/^/        /' || true
        fi
      fi
    done <<<"$r_list"

    if [[ $bad -eq 0 ]]; then
      echo "ok  ($(echo "$r_list" | wc -l | tr -d ' ') outputs identical)"
      pass=$((pass + 1))
    else
      echo "      ($bad file(s) differ)"
      fail=$((fail + 1))
      [[ $KEEP_GOING -eq 1 ]] || { echo; echo "stopping (use --keep-going to continue)"; return 1; }
    fi
  done < <(selected_cases)

  echo
  echo "==> $pass passed, $fail failed, $xfail xfail"
  [[ $fail -eq 0 ]]
}

case "$MODE" in
  record) do_record ;;
  verify) do_verify ;;
  both)   do_record && do_verify ;;
  *) echo "usage: $0 {record|verify|both|list} [options]" >&2; exit 2 ;;
esac
