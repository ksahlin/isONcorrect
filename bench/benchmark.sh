#!/usr/bin/env bash
#
# Performance harness: wall clock and peak RSS, Python reference vs Rust port,
# over the same corpus the equivalence harness uses.
#
# Numbers only mean something next to a passing equivalence run. A port that is
# 40x faster and produces different bases is not a port.
#
#   bench/benchmark.sh                       # both impls, default + paper params
#   bench/benchmark.sh --impl py             # reference only (use before the port exists)
#   bench/benchmark.sh --threads 1,4,8       # scaling sweep for run_isoncorrect
#   bench/benchmark.sh --reps 3              # repeat and keep the fastest
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_ENV="${ISONCORRECT_REF_ENV:-isoncorrect-ref}"

CORPUS=""
IMPLS="py,rs"
THREADS="1,4"
REPS=1
OUT="$REPO_ROOT/bench/results"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --corpus)  CORPUS="$2"; shift 2 ;;
    --impl)    IMPLS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --reps)    REPS="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$CORPUS" ]]; then
  if compgen -G "$REPO_ROOT/bench/corpus/*.fastq" >/dev/null; then
    CORPUS="$REPO_ROOT/bench/corpus"
  else
    CORPUS="$REPO_ROOT/test_data/isoncorrect"
    echo "note: using test_data/isoncorrect; too small for meaningful timings" >&2
  fi
fi
CORPUS="$(cd "$CORPUS" && pwd)"

conda_prefix_for() { conda env list | awk -v n="$1" '$1==n {print $NF; exit}'; }

PY_PREFIX="$(conda_prefix_for "$REF_ENV" || true)"
PY_ISONCORRECT="${PY_ISONCORRECT:-$PY_PREFIX/bin/isONcorrect}"
PY_RUN_ISONCORRECT="${PY_RUN_ISONCORRECT:-$PY_PREFIX/bin/run_isoncorrect}"
RS_ISONCORRECT="${RS_ISONCORRECT:-$REPO_ROOT/rust/target/release/isONcorrect}"
RS_RUN_ISONCORRECT="${RS_RUN_ISONCORRECT:-$REPO_ROOT/rust/target/release/run_isoncorrect}"

[[ -n "$PY_PREFIX" ]] && export PATH="$PY_PREFIX/bin:$PATH"
export PYTHONHASHSEED=0

# GNU time on Linux, BSD time on macOS; both report peak RSS but in different
# units and under different labels.
TIME_BIN="/usr/bin/time"
[[ -x "$TIME_BIN" ]] || { echo "error: /usr/bin/time not found" >&2; exit 1; }

case "$(uname -s)" in
  Darwin) TIME_FLAG="-l"; RSS_UNIT="bytes" ;;
  *)      TIME_FLAG="-v"; RSS_UNIT="kbytes" ;;
esac

# measure <logfile> <cmd...>  -> echoes "<wall_seconds> <peak_rss_mb>"
measure() {
  local log="$1"; shift
  local tf; tf="$(mktemp)"
  local start end wall rss

  start=$(python3 -c 'import time; print(time.time())')
  "$TIME_BIN" "$TIME_FLAG" "$@" >>"$log" 2>"$tf" || {
    echo "    ! command failed, see $log" >&2
    cat "$tf" >> "$log"
    rm -f "$tf"
    echo "NA NA"
    return 1
  }
  end=$(python3 -c 'import time; print(time.time())')
  wall=$(python3 -c "print(f'{$end - $start:.2f}')")

  if [[ "$RSS_UNIT" == "bytes" ]]; then
    rss=$(awk '/maximum resident set size/ {print $1/1048576; exit}' "$tf")
  else
    rss=$(awk -F': ' '/Maximum resident set size/ {print $2/1024; exit}' "$tf")
  fi
  [[ -z "$rss" ]] && rss="NA"
  cat "$tf" >> "$log"
  rm -f "$tf"
  printf '%s %.1f\n' "$wall" "${rss:-0}"
}

mkdir -p "$OUT"
STAMP="$(date +%Y%m%d-%H%M%S)"
CSV="$OUT/benchmark-$STAMP.csv"
LOGDIR="$OUT/logs-$STAMP"
mkdir -p "$LOGDIR"

echo "impl,mode,params,threads,rep,wall_s,peak_rss_mb" > "$CSV"

echo "==> corpus: $CORPUS ($(ls "$CORPUS"/*.fastq 2>/dev/null | wc -l | tr -d ' ') clusters)"
echo "==> results: $CSV"
echo

PARAM_SETS=(
  "default|"
  "paper|--k 9 --w 10 --max_seqs 1000"
)

printf '%-6s %-8s %-10s %-8s %10s %12s\n' IMPL MODE PARAMS THREADS WALL_S RSS_MB

for impl in ${IMPLS//,/ }; do
  if [[ "$impl" == "py" ]]; then
    single="$PY_ISONCORRECT"; folder="$PY_RUN_ISONCORRECT"
  else
    single="$RS_ISONCORRECT"; folder="$RS_RUN_ISONCORRECT"
  fi
  if [[ ! -x "$single" ]]; then
    echo "  (skipping $impl: $single not executable)"
    continue
  fi

  for ps in "${PARAM_SETS[@]}"; do
    IFS='|' read -r pname pargs <<<"$ps"

    # --- single cluster, one core ---
    for rep in $(seq 1 "$REPS"); do
      work="$(mktemp -d)"; log="$LOGDIR/$impl-single-$pname-r$rep.log"
      total_wall=0; max_rss=0
      ok=1
      for f in "$CORPUS"/*.fastq; do
        stem="$(basename "$f" .fastq)"
        # shellcheck disable=SC2086
        read -r w r < <(measure "$log" "$single" --fastq "$f" --outfolder "$work/$stem" $pargs) || { ok=0; break; }
        total_wall=$(python3 -c "print(f'{$total_wall + $w:.2f}')")
        max_rss=$(python3 -c "print(f'{max($max_rss, $r):.1f}')")
      done
      rm -rf "$work"
      if [[ $ok -eq 1 ]]; then
        printf '%-6s %-8s %-10s %-8s %10s %12s\n' "$impl" single "$pname" - "$total_wall" "$max_rss"
        echo "$impl,single,$pname,1,$rep,$total_wall,$max_rss" >> "$CSV"
      fi
    done

    # --- whole folder, thread sweep ---
    if [[ -x "$folder" ]]; then
      for t in ${THREADS//,/ }; do
        for rep in $(seq 1 "$REPS"); do
          work="$(mktemp -d)"; log="$LOGDIR/$impl-folder-$pname-t$t-r$rep.log"
          # shellcheck disable=SC2086
          if read -r w r < <(measure "$log" "$folder" --fastq_folder "$CORPUS" --outfolder "$work" --t "$t" $pargs); then
            printf '%-6s %-8s %-10s %-8s %10s %12s\n' "$impl" folder "$pname" "$t" "$w" "$r"
            echo "$impl,folder,$pname,$t,$rep,$w,$r" >> "$CSV"
          fi
          rm -rf "$work"
        done
      done
    fi
  done
done

echo
echo "==> wrote $CSV"
if command -v python3 >/dev/null 2>&1; then
  python3 - "$CSV" <<'PY'
import csv, sys
from collections import defaultdict

rows = list(csv.DictReader(open(sys.argv[1])))
if not rows:
    sys.exit()
best = defaultdict(lambda: (float("inf"), None))
for r in rows:
    try:
        w = float(r["wall_s"])
    except ValueError:
        continue
    key = (r["mode"], r["params"], r["threads"], r["impl"])
    if w < best[key][0]:
        best[key] = (w, r["peak_rss_mb"])

seen = {(m, p, t) for (m, p, t, _) in best}
print("\n==> speedup (best of reps, Python / Rust)")
for m, p, t in sorted(seen):
    py = best.get((m, p, t, "py"), (None,))[0]
    rs = best.get((m, p, t, "rs"), (None,))[0]
    if py and rs and rs != float("inf") and py != float("inf"):
        print(f"  {m:8s} {p:10s} t={t:3s}  {py:8.2f}s -> {rs:8.2f}s   {py/rs:6.1f}x")
PY
fi
