#!/usr/bin/env bash
# Step 1 of reports/stage2_ppi_01/plan.md: PLM-interact humanV12 scores for the stage-2 LR panel. No OT.
#
#   cd /home/liu/Research/isoform_commot
#   # smoke (6 entity pairs + the longest one; not interpreted)
#   bash reports/stage2_ppi_01/run_scoring.sh --smoke \
#       --input stage2_input_post --input stage2_input_ant \
#       --runs plm_jobs_stage2_ppi_smoke \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D
#   # main (plan.md must be committed and unchanged)
#   bash reports/stage2_ppi_01/run_scoring.sh \
#       --input stage2_input_post --input stage2_input_ant \
#       --runs plm_jobs_stage2_ppi_01 --results reports/stage2_ppi_01/scores \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D
#   (both: --cellchat data/cellchatdb is the default CellChatDB folder used for the unannotated sample)
#
# Environments: commot-research (tests, genes, proteins.py map, jobs, import, designs);
#               plminteract (plm_score.py, GPU). Override with --bench-env / --plm-env.
# --runs: working files, git-ignored by the existing plm_jobs*/ rule. --results: committed (main only);
# do not name it results*/ (that pattern is in .gitignore).
# --proteins DIR reuses an existing proteins.py map output instead of mapping again.
set -euo pipefail

PLAN=reports/stage2_ppi_01/plan.md
BENCH_ENV=commot-research
PLM_ENV=plminteract
CKPT_NAME=danliu1226/PLM-interact-650M-humanV12
CKPT_REV=2d0e7fb3deef94efa3f217daccc79092a6a42ceb
INPUTS=() ; RUNS= ; RES= ; CKPT= ; BASE= ; PROTEINS= ; SMOKE=0 ; CONDA=${CONDA_EXE:-conda}
CELLCHAT=data/cellchatdb ; CROSS_RATIO=3

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUTS+=("$2"); shift 2;;
    --runs) RUNS=$2; shift 2;;
    --results) RES=$2; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;
    --base-model) BASE=$2; shift 2;;
    --proteins) PROTEINS=$2; shift 2;;
    --cellchat) CELLCHAT=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --plm-env) PLM_ENV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ ${#INPUTS[@]} -ge 1 ]] || { echo "--input is required (one per slice)" >&2; exit 2; }
[[ -n "$RUNS" ]] || { echo "--runs is required" >&2; exit 2; }
[[ -n "$CKPT" ]] || { echo "--checkpoint is required" >&2; exit 2; }
[[ -n "$BASE" ]] || { echo "--base-model is required" >&2; exit 2; }
[[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs" >&2; exit 1; }
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/results"; LIMIT=(--limit-pairs 6 --include-longest --seed 0); EXPL=(--exploratory)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  LIMIT=(); EXPL=()
fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || echo none)
mkdir -p "$RUNS"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[stage2-ppi $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)

log "1/7 unit tests"
run "$BENCH_ENV" python -m unittest test_stage2_ppi 2>&1 | tee "$RUNS/unittest.log" | tail -3

log "2/7 subunit genes"
run "$BENCH_ENV" python stage2_ppi.py genes --input "${INPUTS[@]}" --output "$RUNS/genes.csv"

if [[ -z "$PROTEINS" ]]; then
  log "3/7 proteins.py map (mouse Swiss-Prot canonical; unresolved genes listed, never guessed)"
  PROTEINS="$RUNS/proteins"
  run "$BENCH_ENV" python proteins.py map --genes "$RUNS/genes.csv" --species mouse --synonyms --allow-unresolved \
      --output "$PROTEINS"
else
  log "3/7 reusing protein map $PROTEINS"
fi

log "4/7 jobs (both input orders, full length; plus ${CROSS_RATIO}x unannotated combinations)"
run "$BENCH_ENV" python stage2_ppi.py jobs --input "${INPUTS[@]}" --proteins "$PROTEINS" --output "$RUNS/jobs" \
    --cross-ratio "$CROSS_RATIO" --cross-seed 1 --cellchat "$CELLCHAT" \
    ${LIMIT[@]+"${LIMIT[@]}"} | tee "$RUNS/jobs.log"

log "5/7 score ($PLM_ENV, no truncation)"
run "$PLM_ENV" python plm_score.py score --jobs "$RUNS/jobs/plm_input.csv" \
    --checkpoint "$CKPT" --base-model "$BASE" \
    --checkpoint-name "$CKPT_NAME" --checkpoint-revision "$CKPT_REV" \
    --output "$RUNS/scores.csv" --allow-long 2>&1 | tee "$RUNS/score.log" | grep -v "pairs$" || true
[[ -f "$RUNS/scores.csv.meta.json" ]] || { echo "scoring did not complete" >&2; exit 1; }

log "6/7 import (orders averaged; complexes: max over subunit pairs)"
mkdir -p "$(dirname "$RES")"
run "$BENCH_ENV" python stage2_ppi.py import --jobs "$RUNS/jobs" --scores "$RUNS/scores.csv" \
    --output "$RES" ${EXPL[@]+"${EXPL[@]}"}

log "7/7 designs per slice (descriptive only)"
for IN in "${INPUTS[@]}"; do
  NAME=$(basename "$IN")
  run "$BENCH_ENV" python stage2_ppi.py designs --input "$IN" --ppi "$RES/ppi_scores.csv" \
      --output "$RES/designs_${NAME}.csv"
done

cp "$PROTEINS/protein_map.csv" "$RES/protein_map.csv"
cp "$PROTEINS/proteins.json" "$RES/proteins.json"
if [[ -f "${PROTEINS}_candidates.csv" ]]; then cp "${PROTEINS}_candidates.csv" "$RES/unresolved_candidates.csv"; fi
cp "$RUNS/jobs/jobs.json" "$RES/jobs.json"
cp "$RUNS/jobs/entity_pairs.csv" "$RES/entity_pairs.csv"
cp "$RUNS/scores.csv.meta.json" "$RES/scores.meta.json"
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "head_commit_at_start": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "inputs": "${INPUTS[*]}",
  "runs_dir": "$RUNS",
  "bench_env": "$BENCH_ENV",
  "plm_env": "$PLM_ENV"
}
EOF
log "done -> $RES"
