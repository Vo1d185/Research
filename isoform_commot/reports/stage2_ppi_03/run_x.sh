#!/usr/bin/env bash
# Step 2 of reports/stage2_ppi_03/plan.md (section 6): A line (decisions), references, B line (paper protocol).
# Needs step 1 (run_prepare_x.sh, main) committed: reports/stage2_ppi_03/scores/.
#
#   cd /home/liu/Research/isoform_commot
#   bash reports/stage2_ppi_03/run_x.sh --smoke --runs analysis_x_smoke --jobs 8 --gpu-env inject-gpu --gpu-jobs 4
#   bash reports/stage2_ppi_03/run_x.sh --runs analysis_x_main --results reports/stage2_ppi_03/summary \
#       --jobs 8 --gpu-env inject-gpu --gpu-jobs 4          # plan.md committed and unchanged
#   # interrupted: add --resume (finished cases are skipped)
#
# Environments: commot-research (CPU: tests, A line, COMMOT runs, summaries); --gpu-env (torch + CUDA; inject-gpu
# of stage2_inject_01 amendment_01) for the B-line solves of our cost (eps 0.1, rho 10). All paths are options with
# the defaults below. --runs is git-ignored (analysis*/); --results is committed.
set -euo pipefail

PLAN=reports/stage2_ppi_03/plan.md
BENCH_ENV=commot-research ; GPU_ENV=inject-gpu ; CONDA=${CONDA_EXE:-conda}
RUNS= ; RES= ; JOBS=8 ; GPU_JOBS=4 ; BATCH=5 ; SMOKE=0 ; RESUME=0
SCORES=reports/stage2_ppi_03/scores ; EXP2_SCORES=reports/stage2_ppi_01/scores
EXP2_REF=reports/stage2_ppi_02/summary/post_main_receptor_partial_spearman.csv
ZIP=data/zenodo_7272562/downstream_gene_correlation.zip ; DB_MOUSE=data/scseqcomm ; DB_HUMAN=data/scseqcomm_human
POST=stage2_input_post ; ANT=stage2_input_ant ; PGENES=stage2_genes_post_v2 ; AGENES=stage2_genes_ant
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;            --results) RES=${2%/}; shift 2;;
    --jobs) JOBS=$2; shift 2;;                --gpu-jobs) GPU_JOBS=$2; shift 2;;
    --batch) BATCH=$2; shift 2;;              --bench-env) BENCH_ENV=$2; shift 2;;
    --gpu-env) GPU_ENV=$2; shift 2;;          --scores) SCORES=${2%/}; shift 2;;
    --exp2-scores) EXP2_SCORES=${2%/}; shift 2;; --exp2-reference) EXP2_REF=$2; shift 2;;
    --zenodo-zip) ZIP=$2; shift 2;;           --db-mouse) DB_MOUSE=$2; shift 2;;
    --db-human) DB_HUMAN=$2; shift 2;;        --post-input) POST=${2%/}; shift 2;;
    --ant-input) ANT=${2%/}; shift 2;;        --post-genes) PGENES=$2; shift 2;;
    --ant-genes) AGENES=$2; shift 2;;         --smoke) SMOKE=1; shift;;
    --resume) RESUME=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$RUNS" ]] || { echo "--runs is required" >&2; exit 2; }
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; A_SETS=(seqfish); B_SETS=(seqfish); SMOKE_FLAG=(--smoke); B_CASES=(ours_no_ppi)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  A_SETS=(breast seqfish post ant); B_SETS=(post breast seqfish); SMOKE_FLAG=()
  B_CASES=(ours_no_ppi ours_raw ours_snorm ours_uniform)
fi
[[ -f "$SCORES/data_report.json" ]] || { echo "$SCORES/data_report.json missing: run step 1 (run_prepare_x.sh) first" >&2; exit 1; }
if [[ $RESUME -eq 0 ]]; then [[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs or --resume" >&2; exit 1; }; fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs" "$RUNS/cmds"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[exp2x $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
a_args() {   # dataset -> run arguments
  case $1 in
    breast)  echo "--input stage2_input_x_breast --genes stage2_genes_x_breast --db $DB_HUMAN --designs $SCORES/human_A/designs_stage2_input_x_breast.csv --ppi-scores $SCORES/human_A/ppi_scores.csv --protein-map $SCORES/human_A/protein_map.csv";;
    seqfish) echo "--input stage2_input_x_seqfish --genes stage2_genes_x_seqfish --db $DB_MOUSE --designs $SCORES/mouse_A/designs_stage2_input_x_seqfish.csv --ppi-scores $SCORES/mouse_A/ppi_scores.csv --protein-map $SCORES/mouse_A/protein_map.csv";;
    post)    echo "--input $POST --genes $PGENES --db $DB_MOUSE --designs $EXP2_SCORES/designs_$(basename "$POST").csv --ppi-scores $EXP2_SCORES/ppi_scores.csv --protein-map $EXP2_SCORES/protein_map.csv";;
    ant)     echo "--input $ANT --genes $AGENES --db $DB_MOUSE --designs $EXP2_SCORES/designs_$(basename "$ANT").csv --ppi-scores $EXP2_SCORES/ppi_scores.csv --protein-map $EXP2_SCORES/protein_map.csv";;
  esac
}
a_input() { a_args "$1" | awk '{print $2}'; }
a_genes() { a_args "$1" | awk '{print $4}'; }
a_db() { a_args "$1" | awk '{print $6}'; }
b_group() { case $1 in breast) echo human_B;; *) echo mouse_B;; esac; }
b_args() { local g; g=$(b_group "$1"); echo "--input stage2_input_xb_$1 --designs $SCORES/$g/designs_stage2_input_xb_$1.csv --ppi-scores $SCORES/$g/ppi_scores.csv --protein-map $SCORES/$g/protein_map.csv --zenodo-zip $ZIP"; }
cutoff_of() { run "$BENCH_ENV" python -c "import json,sys; print(json.load(open(sys.argv[1]))['diffusion_cutoff_um'])" "$1/config.json"; }
execute() {   # command file, parallel
  local f=$1 k=0 tag; tag=$(basename "$f" .txt)
  rm -f "$RUNS"/cmds/"$tag"_*.sh
  while IFS= read -r line; do
    k=$((k + 1)); printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\n' "$line" > "$RUNS/cmds/${tag}_$(printf %04d $k).sh"
  done < "$f"
  [[ $k -eq 0 ]] && return 0
  ls "$RUNS"/cmds/"$tag"_*.sh | xargs -P "$2" -n 1 bash || { echo "some runs failed; see $RUNS/logs" >&2; exit 1; }
}
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE resume=$RESUME jobs=$JOBS gpu_jobs=$GPU_JOBS plan_commit=$PLAN_COMMIT"

log "1/9 unit tests"
run "$BENCH_ENV" python -m unittest test_exp2x test_stage2_s 2>&1 | tee "$RUNS/logs/unittest.log" | tail -3

log "2/9 step 0a: the paper's metric on the paper's own output"
for ds in "${B_SETS[@]}"; do
  [[ -f "$RUNS/${ds}_B_paper_ref/paper_ref.json" ]] || run "$BENCH_ENV" python exp2x.py paper-ref --input "stage2_input_xb_$ds" \
      --zenodo-zip "$ZIP" --output "$RUNS/${ds}_B_paper_ref" 2>&1 | grep -v -i warning | tee -a "$RUNS/run.log"
done

log "3/9 A-line case lists"
for ds in "${A_SETS[@]}" post; do
  [[ -f "$RUNS/${ds}_A/cases.json" ]] || run "$BENCH_ENV" python exp2x.py cases --dataset "$ds" \
      --output "$RUNS/${ds}_A/cases.json" ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"}
done

log "4/9 gates: regression on the Exp2 posterior (no_ppi, raw); equivalence with evaluate_targets.py (breast, seqfish)"
if [[ ! -f "$RUNS/regression.json" ]]; then
  run "$BENCH_ENV" python exp2x.py run --dataset post $(a_args post) --cases "$RUNS/post_A/cases.json" --output "$RUNS/post_A" \
      --case no_ppi raw > "$RUNS/logs/regression_run.log" 2>&1
  run "$BENCH_ENV" python exp2x.py regression --runs "$RUNS/post_A" --reference "$EXP2_REF" --case no_ppi raw \
      --output "$RUNS/regression.json" | tee -a "$RUNS/run.log" || { echo "regression gate failed" >&2; exit 1; }
fi
for ds in breast seqfish; do
  [[ -f "$RUNS/equivalence_$ds.json" ]] && continue
  [[ -f "$RUNS/${ds}_A/cases.json" ]] || run "$BENCH_ENV" python exp2x.py cases --dataset "$ds" --output "$RUNS/${ds}_A/cases.json" \
      ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"}
  rm -rf "$RUNS/full_$ds" "$RUNS/eval_gate_$ds"
  run "$BENCH_ENV" python exp2x.py run --dataset "$ds" $(a_args "$ds") --cases "$RUNS/${ds}_A/cases.json" --output "$RUNS/${ds}_A" \
      --case no_ppi raw --save-case "$RUNS/full_$ds" > "$RUNS/logs/gate_run_$ds.log" 2>&1
  run "$BENCH_ENV" python evaluate_targets.py run --input "$(a_input "$ds")" --genes "$(a_genes "$ds")" --db "$(a_db "$ds")" \
      --results "$RUNS/full_$ds/no_ppi" "$RUNS/full_$ds/raw" --output "$RUNS/eval_gate_$ds" > "$RUNS/logs/gate_eval_$ds.log" 2>&1
  run "$BENCH_ENV" python stage2_s.py equivalence --runs "$RUNS/${ds}_A" --eval-dir "$RUNS/eval_gate_$ds" --case no_ppi raw \
      --output "$RUNS/equivalence_$ds.json" | tee -a "$RUNS/run.log" || { echo "equivalence gate ($ds) failed" >&2; exit 1; }
  rm -rf "$RUNS/full_$ds"
done

log "5/9 A line: all cases (CPU, $JOBS in parallel, batches of $BATCH)"
T0=$(date +%s); : > "$RUNS/a_cmds.txt"
for ds in "${A_SETS[@]}"; do
  names=$(run "$BENCH_ENV" python -c "import json,sys; print(' '.join(c['name'] for c in json.load(open(sys.argv[1]))['cases']))" "$RUNS/${ds}_A/cases.json")
  arr=($names)
  for ((i = 0; i < ${#arr[@]}; i += BATCH)); do
    echo "$CONDA run --no-capture-output -n $BENCH_ENV python exp2x.py run --dataset $ds $(a_args "$ds") --cases $RUNS/${ds}_A/cases.json --output $RUNS/${ds}_A --case ${arr[*]:i:BATCH} >> $RUNS/logs/${ds}_A.log 2>&1" >> "$RUNS/a_cmds.txt"
  done
done
execute "$RUNS/a_cmds.txt" "$JOBS"
A_WALL=$(( $(date +%s) - T0 )); log "   A line took ${A_WALL} s"

log "6/9 A-line reference: COMMOT pure collective, eps 0.2, rho 1, 363.0757 um (as Exp2's Beps)"
: > "$RUNS/ref_cmds.txt"
for ds in "${A_SETS[@]}"; do
  [[ -f "$RUNS/${ds}_A/ref_commot_B/case.json" ]] && continue
  echo "rm -rf $RUNS/commot_A_$ds && $CONDA run --no-capture-output -n $BENCH_ENV python commot_baseline.py --input $(a_input "$ds") --output $RUNS/commot_A_$ds --dis-thr-um 363.0757 --weights 1 0 0 0 --nitermax 50000 --stopthr 1e-11 --eps 0.2 --rho 1 > $RUNS/logs/commot_A_$ds.log 2>&1 && $CONDA run --no-capture-output -n $BENCH_ENV python exp2x.py eval-commot --input $(a_input "$ds") --genes $(a_genes "$ds") --db $(a_db "$ds") --results $RUNS/commot_A_$ds --output $RUNS/${ds}_A/ref_commot_B --name ref_commot_B >> $RUNS/logs/commot_A_$ds.log 2>&1" >> "$RUNS/ref_cmds.txt"
done
execute "$RUNS/ref_cmds.txt" "$JOBS"

log "7/9 B line: COMMOT rerun (default mixture, eps 0.1, rho 10, paper cutoff) and step 0b"
: > "$RUNS/b_commot_cmds.txt"
for ds in "${B_SETS[@]}"; do
  [[ -f "$RUNS/${ds}_B_eval/eval_b.json" ]] && continue
  cut=$(cutoff_of "stage2_input_xb_$ds")
  echo "rm -rf $RUNS/commot_B_$ds && $CONDA run --no-capture-output -n $BENCH_ENV python commot_baseline.py --input stage2_input_xb_$ds --output $RUNS/commot_B_$ds --dis-thr-um $cut > $RUNS/logs/commot_B_$ds.log 2>&1 && $CONDA run --no-capture-output -n $BENCH_ENV python exp2x.py eval-b --input stage2_input_xb_$ds --zenodo-zip $ZIP --results $RUNS/commot_B_$ds --output $RUNS/${ds}_B_eval >> $RUNS/logs/commot_B_$ds.log 2>&1" >> "$RUNS/b_commot_cmds.txt"
done
execute "$RUNS/b_commot_cmds.txt" "$JOBS"

log "8/9 B line: our cost on GPU (gate first: gpu_unot.py equivalence, posterior B input, 200 iterations)"
if [[ ! -f "$RUNS/gpu_equivalence.json" ]]; then
  run "$GPU_ENV" python -c "import torch; assert torch.cuda.is_available(); print(torch.__version__, torch.cuda.get_device_name(0))" | tee "$RUNS/gpu_env.txt"
  run "$GPU_ENV" python gpu_unot.py equivalence --input stage2_input_xb_post --variant none --eps 0.1 --rho 10 --iters 200 \
      --device cuda --output "$RUNS/gpu_equivalence.json" > "$RUNS/logs/gpu_gate.log" 2>&1 \
    || { echo "GPU equivalence gate failed; see $RUNS/logs/gpu_gate.log" >&2; exit 1; }
fi
T0=$(date +%s); : > "$RUNS/b_cmds.txt"
for ds in "${B_SETS[@]}"; do for c in "${B_CASES[@]}"; do
  echo "$CONDA run --no-capture-output -n $GPU_ENV python exp2x.py run-b $(b_args "$ds") --output $RUNS/${ds}_B --case $c --device cuda >> $RUNS/logs/${ds}_B.log 2>&1" >> "$RUNS/b_cmds.txt"
done; done
execute "$RUNS/b_cmds.txt" "$GPU_JOBS"
B_WALL=$(( $(date +%s) - T0 )); log "   B line (GPU) took ${B_WALL} s"

log "9/9 summaries and provenance"
mkdir -p "$RES"
run "$BENCH_ENV" python exp2x.py summarise --runs "$RUNS" --output "$RES/a_line" --dataset "${A_SETS[@]}" | tee "$RUNS/logs/summary_a.log"
run "$BENCH_ENV" python exp2x.py summarise-b --runs "$RUNS" --output "$RES/b_line" --dataset "${B_SETS[@]}" | tee "$RUNS/logs/summary_b.log"
mkdir -p "$RES/cases" "$RES/gates" "$RES/b_line/steps"
for ds in "${A_SETS[@]}"; do cp "$RUNS/${ds}_A/cases.json" "$RES/cases/${ds}_A.json"; done
cp "$RUNS/regression.json" "$RUNS"/equivalence_*.json "$RES/gates/"
[[ -f "$RUNS/gpu_equivalence.json" ]] && cp "$RUNS/gpu_equivalence.json" "$RUNS/gpu_env.txt" "$RES/gates/"
for ds in "${B_SETS[@]}"; do
  cp "$RUNS/${ds}_B_paper_ref/paper_ref.json" "$RES/b_line/steps/${ds}_paper_ref.json"
  [[ -f "$RUNS/${ds}_B_eval/eval_b.json" ]] && cp "$RUNS/${ds}_B_eval/eval_b.json" "$RES/b_line/steps/${ds}_eval_b.json"
  [[ -f "$RUNS/${ds}_B_eval/rerun_vs_paper_signal.csv" ]] && cp "$RUNS/${ds}_B_eval/rerun_vs_paper_signal.csv" "$RES/b_line/steps/${ds}_rerun_vs_paper_signal.csv"
done
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN", "plan_commit": "$PLAN_COMMIT", "head_commit": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE, "resumed": $RESUME, "runs_dir": "$RUNS", "scores": "$SCORES", "exp2_scores": "$EXP2_SCORES",
  "a_datasets": "${A_SETS[*]}", "b_datasets": "${B_SETS[*]}", "jobs": $JOBS, "gpu_jobs": $GPU_JOBS, "batch": $BATCH,
  "a_wall_seconds": $A_WALL, "b_gpu_wall_seconds": $B_WALL,
  "run_started_utc": "$START", "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES"
