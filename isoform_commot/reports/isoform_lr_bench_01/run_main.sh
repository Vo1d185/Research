#!/usr/bin/env bash
# Main (pre-registered) run of the isoform ranking benchmark, plan.md section 7.
# Refuses to start unless plan.md is committed and unchanged in the working tree.
#
#   cd /home/liu/Research/isoform_commot
#   bash reports/isoform_lr_bench_01/run_main.sh \
#       --runs isoform_lr_bench_runs/main_01 \
#       --results reports/isoform_lr_bench_01/results_humanV12 \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
#       --base-model models/esm2_t33_650M_UR50D \
#       --training-root data/lr_benchmark/plm_training
#
# Environments: commot-research (jobs, merge, evaluate); plminteract (plm_score.py, GPU).
# --runs holds jobs and raw scores (git-ignored); --results is committed.
set -euo pipefail

PLAN=reports/isoform_lr_bench_01/plan.md
DATASET=reports/isoform_lr_bench_01/dataset/benchmark.csv.gz
BENCH_ENV=commot-research
PLM_ENV=plminteract
CKPT_NAME=danliu1226/PLM-interact-650M-humanV12
CKPT_REV=2d0e7fb3deef94efa3f217daccc79092a6a42ceb
MODEL=humanV12
RUNS= ; RES= ; CKPT= ; BASE= ; TRAIN= ; CONDA=${CONDA_EXE:-conda}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=$2; shift 2;;
    --results) RES=$2; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;
    --base-model) BASE=$2; shift 2;;
    --training-root) TRAIN=$2; shift 2;;
    --dataset) DATASET=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --plm-env) PLM_ENV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$RUNS" ]] || { echo "--runs is required" >&2; exit 2; }
[[ -n "$RES" ]] || { echo "--results is required" >&2; exit 2; }
[[ -n "$CKPT" ]] || { echo "--checkpoint is required" >&2; exit 2; }
[[ -n "$BASE" ]] || { echo "--base-model is required" >&2; exit 2; }
[[ -n "$TRAIN" ]] || { echo "--training-root is required" >&2; exit 2; }
[[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs" >&2; exit 1; }
[[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }

# --- pre-registration gate
git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN")
PLAN_TIME=$(git log -1 --format=%cI -- "$PLAN")
HEAD_COMMIT=$(git rev-parse HEAD)
echo "plan.md committed in $PLAN_COMMIT at $PLAN_TIME; HEAD $HEAD_COMMIT"

mkdir -p "$RUNS"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[main $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)

log "1/4 jobs (all contrasts, both orientations)"
run "$BENCH_ENV" python isoform_lr_bench.py jobs --dataset "$DATASET" --output-dir "$RUNS/jobs" > "$RUNS/jobs.log"

log "2/4 score ($PLM_ENV, full length, no truncation)"
run "$PLM_ENV" python plm_score.py score --jobs "$RUNS/jobs/plm_input.csv" \
    --checkpoint "$CKPT" --base-model "$BASE" \
    --checkpoint-name "$CKPT_NAME" --checkpoint-revision "$CKPT_REV" \
    --output "$RUNS/scores.csv" --allow-long 2>&1 | tee "$RUNS/score.log" | grep -v "pairs$" || true
[[ -f "$RUNS/scores.csv.meta.json" ]] || { echo "scoring did not complete" >&2; exit 1; }

log "3/4 merge"
run "$BENCH_ENV" python isoform_lr_bench.py merge --dataset "$DATASET" --jobs-dir "$RUNS/jobs" \
    --scores "$RUNS/scores.csv" --output "$RUNS/merged.csv"

log "4/4 evaluate (pre-registered; training overlap $MODEL)"
run "$BENCH_ENV" python isoform_lr_bench.py evaluate --scored "$RUNS/merged.csv" --output-dir "$RES" \
    --training-root "$TRAIN" --model "$MODEL" | tee "$RUNS/evaluate.log"

cp "$RUNS/jobs/jobs.json" "$RES/jobs.json"
cp "$RUNS/scores.csv.meta.json" "$RES/scores.meta.json"
cp "$RUNS/merged.csv.meta.json" "$RES/merged.meta.json"
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "plan_commit_time": "$PLAN_TIME",
  "head_commit_at_start": "$HEAD_COMMIT",
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "runs_dir": "$RUNS",
  "bench_env": "$BENCH_ENV",
  "plm_env": "$PLM_ENV"
}
EOF
log "done -> $RES"
