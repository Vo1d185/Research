#!/usr/bin/env bash
# Experiment 4 (reports/bindpred_ot_01/plan.md): BindPred S' in the OT cost on the two Visium slices.
#
#   cd /home/liu/Research/isoform_commot
#   # smoke (2 permutations per slice; results are not interpreted)
#   bash reports/bindpred_ot_01/run_bindpred_ot.sh --smoke --runs analysis_bpot_smoke \
#       --post-input stage2_input_post --ant-input stage2_input_ant \
#       --post-genes stage2_genes_post_v2 --ant-genes stage2_genes_ant --db data/scseqcomm \
#       --bindpred-dir models/bindpred --esm2 models/esm2_t33_650M_UR50D \
#       --plm-python /home/liu/.conda/envs/plminteract/bin/python --jobs 8
#   # main (plan.md committed and unchanged): same arguments without --smoke, plus
#   #   --runs analysis_bpot_main --results reports/bindpred_ot_01/summary
#   # interrupted: add --resume (finished cases and steps are skipped)
#
# Environments: commot-research for everything except the ESM-2 embeddings (--plm-python).
# Needs catboost + biopython in commot-research (installed for experiment 3) and models/bindpred (experiment 3).
set -euo pipefail

PLAN=reports/bindpred_ot_01/plan.md
SEQS=reports/bindpred_ot_01/inputs/mouse_panel_sequences.tsv
SCORES=reports/stage2_ppi_01/scores          # V12 scores and protein map (experiment-2 inputs)
EXP2=reports/stage2_ppi_02/summary           # committed experiment-2 receptor values (reproducibility gate)
BENCH_ENV=commot-research
RUNS= ; RES= ; POST= ; ANT= ; PGENES= ; AGENES= ; DB=data/scseqcomm ; BPDIR= ; ESM2= ; PLM_PY= ; DEVICE=
JOBS=8 ; BATCH=5 ; SMOKE=0 ; RESUME=0 ; CONDA=${CONDA_EXE:-conda}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;
    --results) RES=$2; shift 2;;
    --post-input) POST=${2%/}; shift 2;;
    --ant-input) ANT=${2%/}; shift 2;;
    --post-genes) PGENES=$2; shift 2;;
    --ant-genes) AGENES=$2; shift 2;;
    --db) DB=$2; shift 2;;
    --bindpred-dir) BPDIR=${2%/}; shift 2;;
    --esm2) ESM2=$2; shift 2;;
    --plm-python) PLM_PY=$2; shift 2;;
    --device) DEVICE=$2; shift 2;;
    --scores-dir) SCORES=${2%/}; shift 2;;
    --exp2-summary) EXP2=${2%/}; shift 2;;
    --sequences) SEQS=$2; shift 2;;
    --jobs) JOBS=$2; shift 2;;
    --batch) BATCH=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --resume) RESUME=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
need() { [[ -n "$1" ]] || { echo "missing $2" >&2; exit 2; }; }
need "$RUNS" --runs; need "$POST" --post-input; need "$ANT" --ant-input; need "$PGENES" --post-genes
need "$AGENES" --ant-genes; need "$BPDIR" --bindpred-dir; need "$ESM2" --esm2; need "$PLM_PY" --plm-python
for f in "$BPDIR/ESM2_BindPred.cbm" "$SEQS" "$SCORES/protein_map.csv" "$SCORES/ppi_scores.csv" \
         "$EXP2/post_main_receptor_partial_spearman.csv" "$EXP2/ant_main_receptor_partial_spearman.csv"; do
  [[ -f "$f" ]] || { echo "missing input $f" >&2; exit 2; }
done
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; SMOKE_FLAG=(--smoke)
else
  need "$RES" --results
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  SMOKE_FLAG=()
fi
if [[ $RESUME -eq 0 && -e "$RUNS" ]]; then echo "$RUNS exists; use a new --runs or --resume" >&2; exit 1; fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs" "$RUNS/cmds"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[bpot $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
input_of() { if [[ $1 == post ]]; then echo "$POST"; else echo "$ANT"; fi; }
genes_of() { if [[ $1 == post ]]; then echo "$PGENES"; else echo "$AGENES"; fi; }
args_of() {   # slice tag -> stage2_s.py run arguments (without --case)
  local inp b src; inp=$(input_of "$1"); b=$(basename "$inp")
  case "$2" in
    bp)    src="--designs $RUNS/scores/designs_bp_$b.csv --ppi-scores $RUNS/scores/ppi_bindpred.csv";;
    bprec) src="--designs $RUNS/scores/designs_bprec_$b.csv --ppi-scores $RUNS/scores/ppi_bindpred_receptor_only.csv";;
    v12)   src="--designs $SCORES/designs_$b.csv --ppi-scores $SCORES/ppi_scores.csv";;
  esac
  echo "$src --input $inp --protein-map $SCORES/protein_map.csv --genes $(genes_of "$1") --db $DB" \
       "--cases $RUNS/$1_$2/cases.json --output $RUNS/$1_$2 --setting main --slice $1 --device cpu"   # one line
}
names_of() { run "$BENCH_ENV" python -c "import json,sys; print(' '.join(c['name'] for c in json.load(open(sys.argv[1]))['cases']))" "$1"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE resume=$RESUME jobs=$JOBS plan_commit=$PLAN_COMMIT"

log "1/8 unit tests"
run "$BENCH_ENV" python -m unittest test_bindpred_ot test_stage2_s test_affinity > "$RUNS/logs/unittest.log" 2>&1 || true
tail -3 "$RUNS/logs/unittest.log"
grep -q '^OK' "$RUNS/logs/unittest.log" || { echo "unit tests failed; see $RUNS/logs/unittest.log" >&2; exit 1; }
if grep -q 'skipped' "$RUNS/logs/unittest.log"; then echo "some tests were skipped (catboost missing?)" >&2; exit 1; fi

log "2/8 subunit pairs of the panel (ligand first)"
if [[ ! -f "$RUNS/jobs/jobs.json" ]]; then
  rm -rf "$RUNS/jobs"
  run "$BENCH_ENV" python bindpred_ot.py jobs --input "$POST" --input "$ANT" --protein-map "$SCORES/protein_map.csv" \
      --sequences "$SEQS" --output-dir "$RUNS/jobs" | tee "$RUNS/logs/jobs.log"
fi

log "3/8 ESM-2 embeddings"
T0=$(date +%s)
if [[ ! -f "$RUNS/embeddings.npz" ]]; then
  "$PLM_PY" affinity.py embed --pairs "$RUNS/jobs/pairs.csv" --esm2 "$ESM2" --output "$RUNS/embeddings.npz" \
      ${DEVICE:+--device "$DEVICE"} > "$RUNS/logs/embed.log" 2>&1
fi
EMBED_S=$(( $(date +%s) - T0 ))

log "4/8 BindPred predictions and designs"
if [[ ! -f "$RUNS/bindpred.csv" ]]; then
  run "$BENCH_ENV" python affinity.py predict --pairs "$RUNS/jobs/pairs.csv" --embeddings "$RUNS/embeddings.npz" \
      --model "$BPDIR/ESM2_BindPred.cbm" --output "$RUNS/bindpred.csv" | tee "$RUNS/logs/predict.log"
fi
if [[ ! -f "$RUNS/scores/import.json" ]]; then
  rm -rf "$RUNS/scores"
  run "$BENCH_ENV" python bindpred_ot.py import --jobs-dir "$RUNS/jobs" --bindpred "$RUNS/bindpred.csv" --output-dir "$RUNS/scores"
  for sl in post ant; do
    inp=$(input_of $sl); b=$(basename "$inp")
    run "$BENCH_ENV" python stage2_ppi.py designs --input "$inp" --ppi "$RUNS/scores/ppi_bindpred.csv" \
        --output "$RUNS/scores/designs_bp_$b.csv" > "$RUNS/logs/designs_bp_$sl.log"
    run "$BENCH_ENV" python stage2_ppi.py designs --input "$inp" --ppi "$RUNS/scores/ppi_bindpred_receptor_only.csv" \
        --output "$RUNS/scores/designs_bprec_$b.csv" > "$RUNS/logs/designs_bprec_$sl.log"
  done
fi

log "5/8 case lists"
if [[ ! -f "$RUNS/post_bp/cases.json" ]]; then
  run "$BENCH_ENV" python bindpred_ot.py cases --runs "$RUNS" ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"}
fi

log "6/8 reproducibility gate (no_ppi and V12 S' beta 1 against experiment 2)"
if [[ ! -f "$RUNS/gate.json" ]]; then
  for sl in post ant; do
    run "$BENCH_ENV" python stage2_s.py run $(args_of $sl bp) --case no_ppi >> "$RUNS/logs/gate.log" 2>&1
    run "$BENCH_ENV" python stage2_s.py run $(args_of $sl v12) --case v12_snorm >> "$RUNS/logs/gate.log" 2>&1
  done
  run "$BENCH_ENV" python bindpred_ot.py check --runs "$RUNS" --exp2-summary "$EXP2" --output "$RUNS/gate.json" \
      | tee -a "$RUNS/run.log"
  [[ -f "$RUNS/gate.json" ]] && grep -q '"passed": true' "$RUNS/gate.json" \
      || { echo "reproducibility gate failed; see $RUNS/gate.json" >&2; exit 1; }
fi

log "7/8 OT cases (CPU, $JOBS in parallel)"
T0=$(date +%s); : > "$RUNS/cmds.txt"
for sl in post ant; do
  for tag in bp bprec; do
    arr=($(names_of "$RUNS/${sl}_$tag/cases.json"))
    for ((i = 0; i < ${#arr[@]}; i += BATCH)); do
      echo "$CONDA run --no-capture-output -n $BENCH_ENV python stage2_s.py run $(args_of $sl $tag) --case ${arr[*]:i:BATCH} >> $RUNS/logs/${sl}_$tag.log 2>&1" >> "$RUNS/cmds.txt"
    done
  done
done
rm -f "$RUNS"/cmds/run_*.sh; k=0
while IFS= read -r line; do
  k=$((k + 1)); printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\n' "$line" > "$RUNS/cmds/run_$(printf %04d $k).sh"
done < "$RUNS/cmds.txt"
ls "$RUNS"/cmds/run_*.sh | xargs -P "$JOBS" -n 1 bash || { echo "some runs failed; see $RUNS/logs" >&2; exit 1; }
OT_S=$(( $(date +%s) - T0 ))

log "8/8 summary"
run "$BENCH_ENV" python bindpred_ot.py summarise --runs "$RUNS" --output "$RES" | tee "$RUNS/logs/summary.log"
mkdir -p "$RES/scores" "$RES/cases"
cp "$RUNS/gate.json" "$RUNS/scores/import.json" "$RUNS/scores/ppi_bindpred.csv" "$RUNS/jobs/jobs.json" "$RES/scores/"
cp "$RUNS/bindpred.csv" "$RUNS/bindpred.csv.meta.json" "$RES/scores/"
for d in "$RUNS"/post_* "$RUNS"/ant_*; do cp "$d/cases.json" "$RES/cases/$(basename "$d").json"; done
DIRTY=true; git diff --quiet 2>/dev/null && DIRTY=false
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN", "plan_commit": "$PLAN_COMMIT", "head_commit": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "git_dirty": $DIRTY, "smoke": $SMOKE, "resumed": $RESUME,
  "runs_dir": "$RUNS", "post_input": "$POST", "ant_input": "$ANT", "post_genes": "$PGENES", "ant_genes": "$AGENES",
  "db": "$DB", "scores_dir": "$SCORES", "exp2_summary": "$EXP2", "bindpred_dir": "$BPDIR", "esm2": "$ESM2",
  "jobs": $JOBS, "batch": $BATCH, "embed_seconds": $EMBED_S, "ot_wall_seconds": $OT_S,
  "run_started_utc": "$START", "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES"
