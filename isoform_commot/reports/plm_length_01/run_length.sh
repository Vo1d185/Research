#!/usr/bin/env bash
# Experiment 1 (reports/plm_length_01/plan.md): humanV12 on long vs short pairs.
#
#   cd /home/liu/Research/isoform_commot
#   # 0) test split, once (commot-research)
#   python plm_length.py fetch --output data/lr_benchmark/plm_training/Bernett_benchmarking/pairs_uniprot_seqs_test.csv
#   # 1) smoke (200 pairs + 10 E2 pairs per label; not interpreted)
#   bash reports/plm_length_01/run_length.sh --smoke --runs analysis_length_smoke \
#       --test-file data/lr_benchmark/plm_training/Bernett_benchmarking/pairs_uniprot_seqs_test.csv \
#       --training-root data/lr_benchmark/plm_training --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
#       --base-model models/esm2_t33_650M_UR50D --plm-python /home/liu/.conda/envs/plminteract/bin/python --workers 2
#   # 2) main (plan.md committed and unchanged): same arguments without --smoke, plus
#   #    --runs analysis_length_main --results reports/plm_length_01/summary
#   # interrupted main run: add --resume (keeps the prepared data; scored chunks are skipped)
#
# Environments: commot-research runs this script's python steps; scoring runs in the environment of --plm-python.
# --runs is git-ignored by the analysis*/ rule; --results is committed (do not name it results*/, which .gitignore ignores).
set -euo pipefail

PLAN=reports/plm_length_01/plan.md
BENCH_ENV=commot-research
LR_SCORES=reports/plm_lr_benchmark/results/pair_scores/H-A__humanV12.csv.gz
RUNS= ; RES= ; TEST= ; TRAIN= ; CKPT= ; BASE= ; PLM_PY= ; WORKERS=2 ; SMOKE=0 ; RESUME=0 ; CONDA=${CONDA_EXE:-conda}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;
    --results) RES=$2; shift 2;;
    --test-file) TEST=$2; shift 2;;
    --training-root) TRAIN=$2; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;
    --base-model) BASE=$2; shift 2;;
    --plm-python) PLM_PY=$2; shift 2;;
    --workers) WORKERS=$2; shift 2;;
    --lr-scores) LR_SCORES=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --resume) RESUME=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
need() { [[ -n "$1" ]] || { echo "missing $2" >&2; exit 2; }; }
need "$RUNS" --runs; need "$TEST" --test-file; need "$TRAIN" --training-root
need "$CKPT" --checkpoint; need "$BASE" --base-model; need "$PLM_PY" --plm-python
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; LIMIT=(--limit 200 --e2-n 10); NBOOT=200
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  LIMIT=(); NBOOT=2000
fi
if [[ $RESUME -eq 1 ]]; then
  [[ -f "$RUNS/data/dataset.json" ]] || { echo "--resume: $RUNS/data/dataset.json not found" >&2; exit 1; }
else
  [[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs or --resume" >&2; exit 1; }
fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[length $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE resume=$RESUME workers=$WORKERS plan_commit=$PLAN_COMMIT"

log "1/5 unit tests"
run "$BENCH_ENV" python -m unittest test_plm_length 2>&1 | tee "$RUNS/logs/unittest.log" | tail -3

if [[ $RESUME -eq 0 ]]; then
  log "2/5 prepare"
  run "$BENCH_ENV" python plm_length.py prepare --test-file "$TEST" --training-root "$TRAIN" \
      --output-dir "$RUNS/data" ${LIMIT[@]+"${LIMIT[@]}"} | tee "$RUNS/logs/prepare.log"
else
  log "2/5 prepare: kept (--resume)"
fi

log "3/5 score ($WORKERS workers)"
SCORE_START=$(date +%s)
pids=()
for ((w = 0; w < WORKERS; w++)); do
  run "$BENCH_ENV" python plm_length.py work --runs "$RUNS" --worker $w --of "$WORKERS" --plm-python "$PLM_PY" \
      --checkpoint "$CKPT" --base-model "$BASE" > "$RUNS/logs/worker_$w.log" 2>&1 &
  pids+=($!)
done
fail=0; for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
[[ $fail -eq 0 ]] || { echo "a scoring worker failed; see $RUNS/logs and $RUNS/scores/*.log (rerun with --resume)" >&2; exit 1; }
SCORE_WALL=$(( $(date +%s) - SCORE_START ))
log "   scoring took ${SCORE_WALL} s"

log "4/5 merge"
run "$BENCH_ENV" python plm_length.py merge --runs "$RUNS" --output "$RUNS/scores.csv" | tee "$RUNS/logs/merge.log"

log "5/5 evaluate"
LR=(); [[ -f "$LR_SCORES" ]] && LR=(--lr-scores "$LR_SCORES")
run "$BENCH_ENV" python plm_length.py evaluate --runs "$RUNS" --scores "$RUNS/scores.csv" --output-dir "$RES" \
    --n-boot "$NBOOT" ${LR[@]+"${LR[@]}"} | tee "$RUNS/logs/evaluate.log"

cp "$RUNS/data/dataset.json" "$RES/dataset.json"; cp "$RUNS/scores.csv.json" "$RES/scores.json"
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "head_commit": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "resumed": $RESUME,
  "runs_dir": "$RUNS",
  "test_file": "$TEST",
  "training_root": "$TRAIN",
  "checkpoint": "$CKPT",
  "plm_python": "$PLM_PY",
  "workers": $WORKERS,
  "n_boot": $NBOOT,
  "lr_scores": "${LR[1]:-}",
  "scoring_wall_seconds": $SCORE_WALL,
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES"
