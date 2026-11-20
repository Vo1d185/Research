#!/usr/bin/env bash
# Experiment 2 (reports/stage2_ppi_02/plan.md): real slices, no S / raw V12 S / the advisor's S'.
#
#   cd /home/liu/Research/isoform_commot
#   bash reports/stage2_ppi_02/run_s.sh --smoke --runs analysis_s_smoke \
#       --post-input stage2_input_post --ant-input stage2_input_ant \
#       --post-genes stage2_genes_post_v2 --ant-genes stage2_genes_ant --db data/scseqcomm \
#       --scores-dir reports/stage2_ppi_01/scores --jobs 8 --gpu-env inject-gpu --gpu-jobs 4
#   # main (plan.md committed and unchanged): same arguments without --smoke, with
#   #   --runs analysis_s_main --results reports/stage2_ppi_02/summary
#   # interrupted: add --resume (finished cases are skipped)
#
# Environments: commot-research for tests, cases, the main setting (CPU) and the summary; --gpu-env (torch + CUDA,
# e.g. inject-gpu from stage2_inject_01 amendment_01) for the second setting with --device cuda.
# --runs is git-ignored by the analysis*/ rule; --results is committed (not results*/, which .gitignore ignores).
set -euo pipefail

PLAN=reports/stage2_ppi_02/plan.md
BENCH_ENV=commot-research
GPU_ENV=inject-gpu
RUNS= ; RES= ; POST= ; ANT= ; PGENES= ; AGENES= ; DB=data/scseqcomm ; SCORES=reports/stage2_ppi_01/scores
JOBS=8 ; GPU_JOBS=4 ; BATCH=5 ; SMOKE=0 ; RESUME=0 ; SKIP_SECOND=0 ; CONDA=${CONDA_EXE:-conda}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS=${2%/}; shift 2;;
    --results) RES=$2; shift 2;;
    --post-input) POST=${2%/}; shift 2;;
    --ant-input) ANT=${2%/}; shift 2;;
    --post-genes) PGENES=$2; shift 2;;
    --ant-genes) AGENES=$2; shift 2;;
    --db) DB=$2; shift 2;;
    --scores-dir) SCORES=$2; shift 2;;
    --jobs) JOBS=$2; shift 2;;
    --gpu-jobs) GPU_JOBS=$2; shift 2;;
    --batch) BATCH=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --gpu-env) GPU_ENV=$2; shift 2;;
    --smoke) SMOKE=1; shift;;
    --resume) RESUME=1; shift;;
    --skip-second) SKIP_SECOND=1; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
need() { [[ -n "$1" ]] || { echo "missing $2" >&2; exit 2; }; }
need "$RUNS" --runs; need "$POST" --post-input; need "$PGENES" --post-genes
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; SLICES=(post); SMOKE_FLAG=(--smoke)
else
  need "$RES" --results; need "$ANT" --ant-input; need "$AGENES" --ant-genes
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  git ls-files --error-unmatch "$PLAN" >/dev/null 2>&1 || { echo "$PLAN is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$PLAN" || { echo "$PLAN has uncommitted changes" >&2; exit 1; }
  SLICES=(post ant); SMOKE_FLAG=()
fi
if [[ $RESUME -eq 0 ]]; then [[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs or --resume" >&2; exit 1; }; fi
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
mkdir -p "$RUNS/logs" "$RUNS/cmds"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[exp2 $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
input_of() { if [[ $1 == post ]]; then echo "$POST"; else echo "$ANT"; fi; }
genes_of() { if [[ $1 == post ]]; then echo "$PGENES"; else echo "$AGENES"; fi; }
common_args() {   # slice setting
  local inp; inp=$(input_of "$1")
  echo "--input $inp --designs $SCORES/designs_$(basename "$inp").csv --ppi-scores $SCORES/ppi_scores.csv" \
       "--protein-map $SCORES/protein_map.csv --genes $(genes_of "$1") --db $DB" \
       "--cases $RUNS/$1_$2/cases.json --output $RUNS/$1_$2 --setting $2 --slice $1"
}
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "smoke=$SMOKE resume=$RESUME jobs=$JOBS gpu_jobs=$GPU_JOBS plan_commit=$PLAN_COMMIT"

log "1/7 unit tests"
run "$BENCH_ENV" python -m unittest test_stage2_s 2>&1 | tee "$RUNS/logs/unittest.log" | tail -3

log "2/7 case lists"
SETTINGS=(main); [[ $SKIP_SECOND -eq 0 ]] && SETTINGS+=(second)
for sl in "${SLICES[@]}"; do for st in "${SETTINGS[@]}"; do
  [[ -f "$RUNS/${sl}_${st}/cases.json" ]] || run "$BENCH_ENV" python stage2_s.py cases --setting "$st" --slice "$sl" \
      --output "$RUNS/${sl}_${st}/cases.json" ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"}
done; done

log "3/7 equivalence gate (posterior, main: this script vs evaluate_targets.py)"
if [[ ! -f "$RUNS/equivalence.json" ]]; then
  rm -rf "$RUNS/full" "$RUNS/eval_gate"
  run "$BENCH_ENV" python stage2_s.py run $(common_args post main) --case no_ppi raw --save-case "$RUNS/full" \
      > "$RUNS/logs/gate_run.log" 2>&1
  run "$BENCH_ENV" python evaluate_targets.py run --input "$POST" --genes "$PGENES" --db "$DB" \
      --results "$RUNS/full/no_ppi" "$RUNS/full/raw" --output "$RUNS/eval_gate" > "$RUNS/logs/gate_eval.log" 2>&1
  run "$BENCH_ENV" python stage2_s.py equivalence --runs "$RUNS/post_main" --eval-dir "$RUNS/eval_gate" \
      --case no_ppi raw --output "$RUNS/equivalence.json" | tee -a "$RUNS/run.log" \
    || { echo "equivalence gate failed; see $RUNS/equivalence.json" >&2; exit 1; }
  rm -rf "$RUNS/full"                                   # full transport files are large; the gate result stays
fi

batches() {   # slice setting batch env device -> command lines
  local names; names=$(run "$BENCH_ENV" python -c "import json,sys; print(' '.join(c['name'] for c in json.load(open(sys.argv[1]))['cases']))" "$RUNS/$1_$2/cases.json")
  local arr=($names) i
  for ((i = 0; i < ${#arr[@]}; i += $3)); do
    echo "$CONDA run --no-capture-output -n $4 python stage2_s.py run $(common_args "$1" "$2") --device $5 --case ${arr[*]:i:$3} >> $RUNS/logs/$1_$2.log 2>&1"
  done
}
execute() {   # command file, parallel
  local f=$1 k=0
  rm -f "$RUNS"/cmds/"$(basename "$f" .txt)"_*.sh
  while IFS= read -r line; do
    k=$((k + 1)); printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\n' "$line" > "$RUNS/cmds/$(basename "$f" .txt)_$(printf %04d $k).sh"
  done < "$f"
  ls "$RUNS"/cmds/"$(basename "$f" .txt)"_*.sh | xargs -P "$2" -n 1 bash || { echo "some runs failed; see $RUNS/logs" >&2; exit 1; }
}

log "4/7 main setting (CPU, $JOBS in parallel)"
T0=$(date +%s); : > "$RUNS/main_cmds.txt"
for sl in "${SLICES[@]}"; do batches "$sl" main "$BATCH" "$BENCH_ENV" cpu >> "$RUNS/main_cmds.txt"; done
execute "$RUNS/main_cmds.txt" "$JOBS"
MAIN_WALL=$(( $(date +%s) - T0 )); log "   main setting took ${MAIN_WALL} s"

SECOND_WALL=0
if [[ $SKIP_SECOND -eq 0 ]]; then
  log "5/7 GPU equivalence gate (posterior, eps 0.1, rho 10, 200 fixed iterations)"
  if [[ ! -f "$RUNS/gpu_equivalence.json" ]]; then
    run "$GPU_ENV" python -c "import torch; assert torch.cuda.is_available(); print(torch.__version__, torch.cuda.get_device_name(0))" \
        | tee "$RUNS/gpu_env.txt"
    run "$GPU_ENV" python gpu_unot.py equivalence --input "$POST" --variant none --eps 0.1 --rho 10 --iters 200 --diffusion-cutoff-um 363.0757 \
        --device cuda --output "$RUNS/gpu_equivalence.json" > "$RUNS/logs/gpu_gate.log" 2>&1 \
      || { echo "GPU equivalence gate failed; see $RUNS/logs/gpu_gate.log" >&2; exit 1; }
  fi
  log "6/7 second setting (GPU, $GPU_JOBS in parallel, one case per process)"
  T0=$(date +%s); : > "$RUNS/second_cmds.txt"
  for sl in "${SLICES[@]}"; do batches "$sl" second 1 "$GPU_ENV" cuda >> "$RUNS/second_cmds.txt"; done
  execute "$RUNS/second_cmds.txt" "$GPU_JOBS"
  SECOND_WALL=$(( $(date +%s) - T0 )); log "   second setting took ${SECOND_WALL} s"
else
  log "5-6/7 second setting skipped (--skip-second)"
fi

log "7/7 summary and provenance"
run "$BENCH_ENV" python stage2_s.py summarise --runs "$RUNS" --output "$RES" | tee "$RUNS/logs/summary.log"
mkdir -p "$RES/cases"
for d in "$RUNS"/*_main "$RUNS"/*_second; do [[ -f "$d/cases.json" ]] && cp "$d/cases.json" "$RES/cases/$(basename "$d").json"; done
cp "$RUNS/equivalence.json" "$RES/"; [[ -f "$RUNS/gpu_equivalence.json" ]] && cp "$RUNS/gpu_equivalence.json" "$RUNS/gpu_env.txt" "$RES/"
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "head_commit": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "resumed": $RESUME,
  "skip_second": $SKIP_SECOND,
  "runs_dir": "$RUNS",
  "post_input": "$POST", "ant_input": "$ANT", "post_genes": "$PGENES", "ant_genes": "$AGENES", "db": "$DB",
  "scores_dir": "$SCORES",
  "jobs": $JOBS, "gpu_jobs": $GPU_JOBS, "batch": $BATCH,
  "main_wall_seconds": $MAIN_WALL,
  "second_wall_seconds": $SECOND_WALL,
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
log "done -> $RES"
