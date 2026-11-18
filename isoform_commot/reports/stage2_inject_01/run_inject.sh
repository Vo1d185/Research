#!/usr/bin/env bash
# Semi-synthetic injection test, reports/stage2_inject_01/plan.md.
#
#   cd /home/liu/Research/isoform_commot
#   # 0) region input (once; same pair list as stage2_input_post is checked below)
#   python stage2.py prepare --sample posterior --data data/visium_mouse_brain_sagittal_posterior \
#       --max-spots 400 --output stage2_input_post_roi400
#   # 1) smoke (2 seeds, oracle + permutation, raw, beta 0.1 and 1, 20000 iterations; not interpreted)
#   bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 \
#       --check-lr stage2_input_post --runs analysis_inject_smoke --jobs 8
#   # 2) main (plan.md must be committed and unchanged)
#   bash reports/stage2_inject_01/run_inject.sh --input stage2_input_post_roi400 --check-lr stage2_input_post \
#       --runs analysis_inject_01 --results reports/stage2_inject_01/summary --jobs 24 \
#       --v12-report reports/stage2_ppi_01/scores/score_report.json
#
# Environment: commot-research for every step (no GPU). --runs is git-ignored by the analysis*/ rule;
# --results is committed (do not name it results*/, which .gitignore ignores).
set -euo pipefail

PLAN=reports/stage2_inject_01/plan.md
BENCH_ENV=commot-research
INPUT= ; CHECK_LR= ; RUNS= ; RES= ; JOBS=8 ; SMOKE=0 ; V12= ; CONDA=${CONDA_EXE:-conda}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT=$2; shift 2;;
    --check-lr) CHECK_LR=$2; shift 2;;
    --runs) RUNS=$2; shift 2;;
    --results) RES=$2; shift 2;;
    --jobs) JOBS=$2; shift 2;;
    --v12-report) V12=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$INPUT" ]] || { echo "--input is required" >&2; exit 2; }
[[ -n "$RUNS" ]] || { echo "--runs is required" >&2; exit 2; }
[[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs" >&2; exit 1; }
if [[ -n "$CHECK_LR" ]]; then
  cmp -s "$INPUT/lr.csv" "$CHECK_LR/lr.csv" || { echo "$INPUT/lr.csv differs from $CHECK_LR/lr.csv" >&2; exit 1; }
fi
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; SEEDS=(1 2); VARIANTS=(q100 q100_perm); NORMS=(raw); BETAS="0.1 1"; NITER=20000
  SMOKE_FLAG=(--allow-unconverged)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  SEEDS=(1 2 3); VARIANTS=(q100 q100_perm q90 q90_perm q80 q80_perm q70 q70_perm q60 q60_perm)
  NORMS=(raw ligand); BETAS="0.05 0.1 0.2 0.5 1"; NITER=200000
  SMOKE_FLAG=()
fi
setting() { case "$1" in main) echo "--eps 0.1 --rho 10";; sens) echo "--eps 0.2 --rho 1";; esac; }
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || echo none)
mkdir -p "$RUNS/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[inject $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)

log "1/5 unit tests"
run "$BENCH_ENV" python -m unittest test_stage2_inject 2>&1 | tee "$RUNS/unittest.log" | tail -3

log "2/5 truth and scores per seed"
for s in "${SEEDS[@]}"; do
  run "$BENCH_ENV" python stage2_inject.py truth --input "$INPUT" --output "$RUNS/truth_s$s" --seed "$s"
  run "$BENCH_ENV" python stage2_inject.py scores --truth "$RUNS/truth_s$s" --output "$RUNS/scores_s$s" --seed "$s" \
      > "$RUNS/logs/scores_s$s.log"
done

log "3/5 OT runs (${JOBS} in parallel)"
CMDS="$RUNS/commands.txt"; : > "$CMDS"
for st in main sens; do
  for s in "${SEEDS[@]}"; do
    common="python stage2_inject.py run --input $INPUT --truth $RUNS/truth_s$s --scores $RUNS/scores_s$s --seed $s $(setting $st) --nitermax $NITER"
    echo "$CONDA run --no-capture-output -n $BENCH_ENV $common --variant none --output $RUNS/$st/s${s}_none > $RUNS/logs/${st}_s${s}_none.log 2>&1" >> "$CMDS"
    echo "$CONDA run --no-capture-output -n $BENCH_ENV $common --variant uniform --beta $BETAS --output $RUNS/$st/s${s}_uniform > $RUNS/logs/${st}_s${s}_uniform.log 2>&1" >> "$CMDS"
    for v in "${VARIANTS[@]}"; do
      for n in "${NORMS[@]}"; do
        echo "$CONDA run --no-capture-output -n $BENCH_ENV $common --variant $v --norm $n --beta $BETAS --output $RUNS/$st/s${s}_${v}_${n} > $RUNS/logs/${st}_s${s}_${v}_${n}.log 2>&1" >> "$CMDS"
      done
    done
  done
done
log "   $(wc -l < "$CMDS") run commands"
mkdir -p "$RUNS/main" "$RUNS/sens" "$RUNS/cmds"
i=0
while IFS= read -r line; do
  i=$((i + 1)); printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\n' "$line" > "$RUNS/cmds/$(printf %04d $i).sh"
done < "$CMDS"
ls "$RUNS"/cmds/*.sh | xargs -P "$JOBS" -n 1 bash || { echo "some runs failed; see $RUNS/logs" >&2; exit 1; }

log "4/5 summaries"
mkdir -p "$RES"
for st in main sens; do
  run "$BENCH_ENV" python stage2_inject.py summarise --runs "$RUNS/$st" --output "$RES/$st" \
      ${V12:+--v12-report "$V12"} ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"} | tee "$RUNS/logs/summary_$st.log"
done

log "5/5 provenance"
mkdir -p "$RES/seeds"
for s in "${SEEDS[@]}"; do
  cp "$RUNS/truth_s$s/truth.json" "$RES/seeds/truth_s$s.json"
  cp "$RUNS/truth_s$s/pairs.csv" "$RES/seeds/pairs_s$s.csv"
  cp "$RUNS/scores_s$s/scores.json" "$RES/seeds/scores_s$s.json"
done
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "head_commit_at_start": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "input": "$INPUT",
  "check_lr": "$CHECK_LR",
  "v12_report": "$V12",
  "seeds": "${SEEDS[*]}",
  "betas": "$BETAS",
  "nitermax": $NITER,
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "runs_dir": "$RUNS",
  "parallel_jobs": $JOBS
}
EOF
log "done -> $RES"
