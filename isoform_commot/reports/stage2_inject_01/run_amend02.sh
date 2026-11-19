#!/usr/bin/env bash
# Finish injection 01 under amendment_02.md, from the existing main run (its OT runs are not repeated).
#
#   cd /home/liu/Research/isoform_commot
#   mv reports/stage2_inject_01/summary analysis_inject_01/summary_partial_unused   # unopened leftover of 11-19
#   # smoke (plumbing only, on the smoke run; not interpreted)
#   bash reports/stage2_inject_01/run_amend02.sh --smoke --runs analysis_inject_smoke --check-runs analysis_amend02_smoke \
#       --check-nitermax 40000 --jobs 8
#   # main (plan.md, amendment_01.md and amendment_02.md committed and unchanged)
#   bash reports/stage2_inject_01/run_amend02.sh --runs analysis_inject_01 --check-runs analysis_inject_01_n400k \
#       --results reports/stage2_inject_01/summary --jobs 8 --v12-report reports/stage2_ppi_01/scores/score_report.json
#
# Steps: unit tests; check-run commands = the 8 STABILITY_RUNS lines (stage2_inject.py) of <runs>/commands.txt with
# --nitermax 400000 and outputs under --check-runs (same env, device, truth, scores); run them; summaries (main:
# --convergence checks, sens: plan rule); stability = used runs vs the longer-cap runs (amendment_02 section 3); provenance.
# --runs must be spelled exactly as in the original run. Environment: commot-research for tests and summaries; the
# check runs use the env and device written in the original commands (inject-gpu / cuda for analysis_inject_01).
set -euo pipefail

PLAN=reports/stage2_inject_01/plan.md
AMEND1=reports/stage2_inject_01/amendment_01.md
AMEND2=reports/stage2_inject_01/amendment_02.md
BENCH_ENV=commot-research
RUNS= ; CHECK= ; RES= ; JOBS=8 ; SMOKE=0 ; V12= ; CHECK_N=400000 ; CONDA=${CONDA_EXE:-conda}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;
    --check-runs) CHECK=${2%/}; shift 2;;
    --results) RES=$2; shift 2;;
    --jobs) JOBS=$2; shift 2;;
    --v12-report) V12=$2; shift 2;;
    --check-nitermax) CHECK_N=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$RUNS" && -n "$CHECK" ]] || { echo "--runs and --check-runs are required" >&2; exit 2; }
[[ -f "$RUNS/commands.txt" ]] || { echo "$RUNS/commands.txt not found" >&2; exit 1; }
[[ -e "$CHECK" ]] && { echo "$CHECK exists; use a new --check-runs" >&2; exit 1; }
committed() {
  git ls-files --error-unmatch "$1" >/dev/null 2>&1 || { echo "$1 is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$1" || { echo "$1 has uncommitted changes" >&2; exit 1; }
}
if [[ $SMOKE -eq 1 ]]; then
  RES="$CHECK/summary"; SMOKE_FLAG=(--allow-unconverged); REQUIRE=()
else
  [[ -n "$RES" ]] || { echo "--results is required" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists (move the unopened leftover into $RUNS first)" >&2; exit 1; }
  [[ "$CHECK_N" == 400000 ]] || { echo "amendment_02.md fixes --check-nitermax 400000" >&2; exit 2; }
  committed "$PLAN"; committed "$AMEND1"; committed "$AMEND2"
  SMOKE_FLAG=(); REQUIRE=(--require-all)
fi
commit_of() { local c; c=$(git log -1 --format=%H -- "$1" 2>/dev/null || true); echo "${c:-none}"; }
mkdir -p "$CHECK/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[amend02 $(date -u +%H:%M:%S)] $*" | tee -a "$CHECK/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)

log "1/6 unit tests"
run "$BENCH_ENV" python -m unittest test_stage2_inject 2>&1 | tee "$CHECK/unittest.log" | tail -3

log "2/6 check-run commands (stability runs of the main setting, --nitermax $CHECK_N)"
run "$BENCH_ENV" python stage2_inject.py recheck-commands --commands "$RUNS/commands.txt" --runs "$RUNS" \
    --check-runs "$CHECK" --nitermax "$CHECK_N" --output "$CHECK/commands.txt" ${REQUIRE[@]+"${REQUIRE[@]}"} \
    | tee "$CHECK/recheck.json"
mkdir -p "$CHECK/main" "$CHECK/cmds"
i=0
while IFS= read -r line; do
  i=$((i + 1)); printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\n' "$line" > "$CHECK/cmds/$(printf %04d $i).sh"
done < "$CHECK/commands.txt"

log "3/6 check runs ($i commands, ${JOBS} in parallel)"
OT_START=$(date +%s)
ls "$CHECK"/cmds/*.sh | xargs -P "$JOBS" -n 1 bash || { echo "some check runs failed; see $CHECK/logs" >&2; exit 1; }
CHECK_WALL=$(( $(date +%s) - OT_START ))
log "   check runs took ${CHECK_WALL} s wall time"

log "4/6 summaries"
mkdir -p "$RES/stability"
run "$BENCH_ENV" python stage2_inject.py summarise --runs "$RUNS/main" --output "$RES/main" --convergence checks \
    ${V12:+--v12-report "$V12"} ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"} | tee "$CHECK/logs/summary_main.log"
run "$BENCH_ENV" python stage2_inject.py summarise --runs "$RUNS/sens" --output "$RES/sens" \
    ${V12:+--v12-report "$V12"} ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"} | tee "$CHECK/logs/summary_sens.log"

log "5/6 stability (amendment_02 section 3)"
run "$BENCH_ENV" python stage2_inject.py stability --used "$RUNS/main" --long "$CHECK/main" \
    --output "$RES/stability/stability_main.json" | tee "$CHECK/logs/stability.log"
cp "$CHECK/recheck.json" "$RES/stability/recheck.json"

log "6/6 provenance"
mkdir -p "$RES/seeds"
for f in "$RUNS"/truth_s*/truth.json; do
  s=$(basename "$(dirname "$f")"); s=${s#truth_}
  cp "$f" "$RES/seeds/truth_$s.json"; cp "$RUNS/truth_$s/pairs.csv" "$RES/seeds/pairs_$s.csv"
  cp "$RUNS/scores_$s/scores.json" "$RES/seeds/scores_$s.json"
done
if [[ -d "$RUNS/equivalence" ]]; then mkdir -p "$RES/equivalence"; cp "$RUNS"/equivalence/*.json "$RES/equivalence/"; fi
if [[ -f "$RUNS/gpu_env.txt" ]]; then cp "$RUNS/gpu_env.txt" "$RES/gpu_env.txt"; fi
MAIN_WALL=$(sed -n 's/.*OT runs took \([0-9]*\) s wall time.*/\1/p' "$RUNS/run.log" | tail -1)
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$(commit_of "$PLAN")",
  "amendment_01_commit": "$(commit_of "$AMEND1")",
  "amendment_02_commit": "$(commit_of "$AMEND2")",
  "head_commit_at_finish_start": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "runs_dir": "$RUNS",
  "check_runs_dir": "$CHECK",
  "recheck": $(cat "$CHECK/recheck.json"),
  "main_convergence_rule": "checks (amendment_02)",
  "sens_convergence_rule": "cap (plan.md)",
  "v12_report": "$V12",
  "main_run_ot_wall_seconds": ${MAIN_WALL:-null},
  "check_ot_wall_seconds": $CHECK_WALL,
  "parallel_jobs": $JOBS,
  "finish_started_utc": "$START",
  "finish_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES ($(grep -o '"verdict": "[^"]*"' "$RES/stability/stability_main.json"))"
