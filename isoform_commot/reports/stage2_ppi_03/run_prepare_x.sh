#!/usr/bin/env bash
# Step 1 of reports/stage2_ppi_03/plan.md (section 6): data + PLM-interact humanV12 scores. No OT, no evaluation.
#
#   cd /home/liu/Research/isoform_commot
#   # smoke (6 entity pairs per scoring group + the longest; inputs are prepared in full; not interpreted)
#   bash reports/stage2_ppi_03/run_prepare_x.sh --smoke --runs analysis_x_prep_smoke \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D
#   # main (plan.md committed and unchanged)
#   bash reports/stage2_ppi_03/run_prepare_x.sh --runs analysis_x_prep --results reports/stage2_ppi_03/scores \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D
#
# Environments: commot-research (tests, downloads, preparation, proteins.py, stage2_ppi.py); plminteract (plm_score.py).
# Every path has a default below and can be changed with the option of the same name. Prepared inputs
# (stage2_input_x*/, stage2_genes_x*/) and downloads (data/) are git-ignored; an existing prepared folder is reused
# (its preparation.json checksums are verified by every later step), never overwritten.
set -euo pipefail

PLAN=reports/stage2_ppi_03/plan.md
BENCH_ENV=commot-research ; PLM_ENV=plminteract ; CONDA=${CONDA_EXE:-conda}
CKPT_NAME=danliu1226/PLM-interact-650M-humanV12 ; CKPT_REV=2d0e7fb3deef94efa3f217daccc79092a6a42ceb
RUNS= ; RES= ; CKPT= ; BASE= ; SMOKE=0
ZEN=data/zenodo_7272562 ; BREAST_DATA=data/visium_human_breast_cancer ; POST_DATA=data/visium_mouse_brain_sagittal_posterior
CC_MOUSE=data/cellchatdb ; CC_HUMAN=data/cellchatdb_human ; DB_HUMAN=data/scseqcomm_human ; WORK=data/paper_outputs
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;          --results) RES=${2%/}; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;        --base-model) BASE=$2; shift 2;;
    --zenodo-dir) ZEN=$2; shift 2;;         --breast-data) BREAST_DATA=$2; shift 2;;
    --post-data) POST_DATA=$2; shift 2;;    --cellchat-mouse) CC_MOUSE=$2; shift 2;;
    --cellchat-human) CC_HUMAN=$2; shift 2;; --db-human) DB_HUMAN=$2; shift 2;;
    --work) WORK=$2; shift 2;;              --bench-env) BENCH_ENV=$2; shift 2;;
    --plm-env) PLM_ENV=$2; shift 2;;        --smoke) SMOKE=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$RUNS" && -n "$CKPT" && -n "$BASE" ]] || { echo "--runs, --checkpoint and --base-model are required" >&2; exit 2; }
[[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs" >&2; exit 1; }
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/scores"; LIMIT=(--limit-pairs 6 --include-longest --seed 0)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  LIMIT=()
fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[exp2x-prep $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE plan_commit=$PLAN_COMMIT"

log "1/5 unit tests"
run "$BENCH_ENV" python -m unittest test_exp2x test_stage2_s test_stage2_ppi 2>&1 | tee "$RUNS/logs/unittest.log" | tail -3

log "2/5 downloads (verified; existing files are checked, not replaced)"
run "$BENCH_ENV" python exp2x_data.py fetch-zenodo --output "$ZEN" --name downstream_gene_correlation.zip \
    visium-mouse-brain.zip visium-human-breast_cancer.zip seqfishplus-mouse-sscortex.zip | tee "$RUNS/logs/zenodo.log"
[[ -f "$CC_HUMAN/db.json" ]] || run "$BENCH_ENV" python exp2x_data.py fetch-cellchat --species human --output "$CC_HUMAN"
[[ -f "$DB_HUMAN/db.json" ]] || run "$BENCH_ENV" python exp2x_data.py fetch-db-human --output "$DB_HUMAN"
[[ -f "$POST_DATA/source.json" ]] || run "$BENCH_ENV" python stage2.py fetch --sample posterior --output "$POST_DATA"

log "3/5 inputs (A line: breast, seqfish; B line: post, breast, seqfish)"
prep() { local dir=$1; shift; if [[ -f "$dir/preparation.json" ]]; then log "   reuse $dir"; else run "$BENCH_ENV" python exp2x_data.py "$@" 2>&1 | grep -v -i warning | tee -a "$RUNS/logs/prepare.log"; fi; }
prep stage2_input_x_breast prepare-a --dataset breast --zenodo-dir "$ZEN" --data "$BREAST_DATA" --cellchat "$CC_HUMAN" \
    --output stage2_input_x_breast --genes-output stage2_genes_x_breast
prep stage2_input_x_seqfish prepare-a --dataset seqfish --zenodo-dir "$ZEN" --cellchat "$CC_MOUSE" \
    --output stage2_input_x_seqfish --genes-output stage2_genes_x_seqfish
prep stage2_input_xb_post prepare-b --dataset post --zenodo-dir "$ZEN" --data "$POST_DATA" --work "$WORK" --output stage2_input_xb_post
prep stage2_input_xb_breast prepare-b --dataset breast --zenodo-dir "$ZEN" --data "$BREAST_DATA" --work "$WORK" --output stage2_input_xb_breast
prep stage2_input_xb_seqfish prepare-b --dataset seqfish --zenodo-dir "$ZEN" --work "$WORK" --output stage2_input_xb_seqfish
for d in stage2_input_x_breast stage2_input_x_seqfish stage2_input_xb_post stage2_input_xb_breast stage2_input_xb_seqfish; do
  [[ -f "$d/preparation.json" ]] || { echo "preparation of $d failed; see $RUNS/logs/prepare.log" >&2; exit 1; }
done

log "4/5 humanV12 scores per group (both input orders, full length; complexes: max over subunit pairs)"
score_group() {   # name species input...
  local name=$1 species=$2; shift 2
  local d="$RUNS/$name" out="$RES/$name"
  mkdir -p "$d"
  log "   $name ($species): $*"
  run "$BENCH_ENV" python stage2_ppi.py genes --input "$@" --output "$d/genes.csv"
  run "$BENCH_ENV" python proteins.py map --genes "$d/genes.csv" --species "$species" --synonyms --allow-unresolved \
      --output "$d/proteins" > "$d/proteins.log" 2>&1
  run "$BENCH_ENV" python stage2_ppi.py jobs --input "$@" --proteins "$d/proteins" --output "$d/jobs" \
      ${LIMIT[@]+"${LIMIT[@]}"} | tee "$d/jobs.log"
  run "$PLM_ENV" python plm_score.py score --jobs "$d/jobs/plm_input.csv" --checkpoint "$CKPT" --base-model "$BASE" \
      --checkpoint-name "$CKPT_NAME" --checkpoint-revision "$CKPT_REV" --output "$d/scores.csv" --allow-long \
      > "$d/score.log" 2>&1
  [[ -f "$d/scores.csv.meta.json" ]] || { echo "scoring of $name did not complete; see $d/score.log" >&2; exit 1; }
  mkdir -p "$(dirname "$out")"
  run "$BENCH_ENV" python stage2_ppi.py import --jobs "$d/jobs" --scores "$d/scores.csv" --output "$out" | tee "$d/import.log"
  local inp
  for inp in "$@"; do
    run "$BENCH_ENV" python stage2_ppi.py designs --input "$inp" --ppi "$out/ppi_scores.csv" \
        --output "$out/designs_$(basename "$inp").csv" > /dev/null
  done
  cp "$d/proteins/protein_map.csv" "$d/proteins/proteins.json" "$out/"
  if [[ -f "$d/proteins_candidates.csv" ]]; then cp "$d/proteins_candidates.csv" "$out/unresolved_candidates.csv"; fi
  cp "$d/jobs/jobs.json" "$d/jobs/entity_pairs.csv" "$out/"; cp "$d/scores.csv.meta.json" "$out/scores.meta.json"
}
score_group mouse_A mouse stage2_input_x_seqfish
score_group human_A human stage2_input_x_breast
score_group mouse_B mouse stage2_input_xb_post stage2_input_xb_seqfish
score_group human_B human stage2_input_xb_breast

log "5/5 data report"
run "$BENCH_ENV" python - "$RES" <<'EOF' | tee "$RUNS/logs/data_report.log"
import json, sys
from pathlib import Path
res = Path(sys.argv[1]); out = {}
for d in ['stage2_input_x_breast', 'stage2_input_x_seqfish', 'stage2_input_xb_post', 'stage2_input_xb_breast', 'stage2_input_xb_seqfish']:
    p = json.loads((Path(d)/'preparation.json').read_text())
    keep = ['dataset', 'line', 'species', 'spots', 'genes', 'neighbour_spacing_um', 'um_per_unit', 'pairs_selected', 'pairs',
            'ligands', 'receptors', 'contact_neighbours_per_cell', 'diffusion_cutoff_um', 'paper_eval', 'paper_input_x_max_abs_diff',
            'normalization', 'prepared_sha256']
    out[d] = {k: p[k] for k in keep if k in p}
for g in ['mouse_A', 'human_A', 'mouse_B', 'human_B']:
    r = json.loads((res/g/'score_report.json').read_text())
    out[g] = {k: r[k] for k in ['n_pairs', 'n_scored', 'n_missing', 'by_type', 'score_distribution', 'order_asymmetry']}
    out[g]['unresolved_genes'] = json.loads((res/g/'proteins.json').read_text()).get('unresolved', [])
(res/'data_report.json').write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
print(json.dumps({k: {kk: v.get(kk) for kk in ['pairs_selected', 'pairs', 'receptors', 'n_pairs', 'n_scored', 'n_missing']} for k, v in out.items()}, indent=1))
EOF
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN", "plan_commit": "$PLAN_COMMIT", "head_commit": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE, "runs_dir": "$RUNS", "zenodo_dir": "$ZEN", "work": "$WORK",
  "checkpoint_name": "$CKPT_NAME", "checkpoint_revision": "$CKPT_REV",
  "run_started_utc": "$START", "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES"
