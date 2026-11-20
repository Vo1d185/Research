#!/usr/bin/env bash
# Experiment 3 (reports/affinity_01/plan.md): can BindPred rank measured binding affinities?
#
#   cd /home/liu/Research/isoform_commot
#   # 0) once: BindPred weights + training table, and catboost in the predict env
#   #    (see reports/affinity_01/README.md)
#   # 1) smoke (3 receptor groups; results are not interpreted)
#   bash reports/affinity_01/run_affinity.sh --smoke --runs analysis_affinity_smoke \
#       --bindpred-dir models/bindpred --esm2 models/esm2_t33_650M_UR50D \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
#       --plm-python /home/liu/.conda/envs/plminteract/bin/python
#   # 2) main (plan.md committed and unchanged): same arguments without --smoke, plus
#   #    --runs analysis_affinity_main --results reports/affinity_01/summary
#   # interrupted: rerun with --resume (finished steps are skipped)
#
# Environments: commot-research for prepare/evaluate (needs biopython), --plm-python for the
# embeddings and the V12 scores, --predict-env for catboost. --runs is git-ignored (analysis*/);
# --results is committed (do not name it results*/).
set -euo pipefail

PLAN=reports/affinity_01/plan.md
INPUTS=reports/affinity_01/inputs
BENCH_ENV=commot-research ; PREDICT_ENV=commot-research
RUNS= ; RES= ; BPDIR= ; ESM2= ; CKPT= ; PLM_PY= ; DEVICE= ; SMOKE=0 ; RESUME=0 ; CONDA=${CONDA_EXE:-conda}
V12_NAME=danliu1226/PLM-interact-650M-humanV12
V12_REV=2d0e7fb3deef94efa3f217daccc79092a6a42ceb
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;
    --results) RES=$2; shift 2;;
    --bindpred-dir) BPDIR=${2%/}; shift 2;;
    --esm2) ESM2=$2; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;
    --plm-python) PLM_PY=$2; shift 2;;
    --inputs) INPUTS=${2%/}; shift 2;;
    --device) DEVICE=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --predict-env) PREDICT_ENV=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --resume) RESUME=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
need() { [[ -n "$1" ]] || { echo "missing $2" >&2; exit 2; }; }
need "$RUNS" --runs; need "$BPDIR" --bindpred-dir; need "$ESM2" --esm2
need "$CKPT" --checkpoint; need "$PLM_PY" --plm-python
MODEL="$BPDIR/ESM2_BindPred.cbm" ; TRAIN_CSV="$BPDIR/embeddings_only.csv"
for f in "$MODEL" "$TRAIN_CSV" "$INPUTS/gtopdb_interactions.csv" "$INPUTS/gtopdb_peptides.csv" \
         "$INPUTS/uniprot_sequences.tsv"; do
  [[ -f "$f" ]] || { echo "missing input $f" >&2; exit 2; }
done
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; LIMIT=(--limit 3)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  LIMIT=()
fi
if [[ $RESUME -eq 1 ]]; then
  [[ -f "$RUNS/data/dataset.json" ]] || { echo "--resume: $RUNS/data/dataset.json not found" >&2; exit 1; }
elif [[ -e "$RUNS" ]]; then
  echo "$RUNS exists; use a new --runs or --resume" >&2; exit 1
fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[affinity $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE resume=$RESUME plan_commit=$PLAN_COMMIT"

log "1/6 unit tests"
run "$BENCH_ENV" python -m unittest test_affinity 2>&1 | tee "$RUNS/logs/unittest.log" | tail -3

log "2/6 prepare (GtoPdb -> pairs, leakage excluded)"
if [[ ! -f "$RUNS/data/dataset.json" ]]; then
  rm -rf "$RUNS/data"
  run "$BENCH_ENV" python affinity.py prepare --gtopdb-dir "$INPUTS" \
      --uniprot "$INPUTS/uniprot_sequences.tsv" --training-csv "$TRAIN_CSV" \
      --output-dir "$RUNS/data" ${LIMIT[@]+"${LIMIT[@]}"} | tee "$RUNS/logs/prepare.log"
fi

log "3/6 embeddings (frozen ESM-2 650M, BindPred recipe)"
[[ -f "$RUNS/embeddings.npz" ]] || \
  run "$BENCH_ENV" env PLM_PY="$PLM_PY" sh -c '"$PLM_PY" affinity.py embed --pairs '"$RUNS"'/data/pairs.csv \
      --esm2 '"$ESM2"' --output '"$RUNS"'/embeddings.npz'"${DEVICE:+ --device $DEVICE}" \
      > "$RUNS/logs/embed.log" 2>&1

log "4/6 BindPred predictions"
[[ -f "$RUNS/bindpred.csv" ]] || \
  run "$PREDICT_ENV" python affinity.py predict --pairs "$RUNS/data/pairs.csv" \
      --embeddings "$RUNS/embeddings.npz" --model "$MODEL" --output "$RUNS/bindpred.csv" \
      | tee "$RUNS/logs/predict.log"

log "5/6 V12 baseline scores"
[[ -f "$RUNS/v12.csv" ]] || \
  run "$BENCH_ENV" env PLM_PY="$PLM_PY" sh -c '"$PLM_PY" plm_score.py score --jobs '"$RUNS"'/data/plm_input.csv \
      --checkpoint '"$CKPT"' --base-model '"$ESM2"' --checkpoint-name '"$V12_NAME"' \
      --checkpoint-revision '"$V12_REV"' --allow-long --output '"$RUNS"'/v12.csv' \
      > "$RUNS/logs/v12.log" 2>&1

log "6/6 evaluate"
run "$BENCH_ENV" python affinity.py evaluate --pairs "$RUNS/data/pairs.csv" --bindpred "$RUNS/bindpred.csv" \
    --v12 "$RUNS/v12.csv" --output "$RES" | tee "$RUNS/logs/evaluate.log"

cat > "$RES/run_info.json" <<JSON
{
  "status": "complete", "started_utc": "$START", "finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "smoke": $SMOKE, "resume": $RESUME, "plan": "$PLAN", "plan_commit": "$PLAN_COMMIT",
  "runs": "$RUNS", "results": "$RES", "inputs": "$INPUTS", "bindpred_dir": "$BPDIR",
  "esm2": "$ESM2", "v12_checkpoint": "$CKPT", "predict_env": "$PREDICT_ENV",
  "git_commit": "$(git rev-parse HEAD)", "git_dirty": $(git diff --quiet && echo false || echo true)
}
JSON
log "done -> $RES"
