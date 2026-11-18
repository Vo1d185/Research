#!/usr/bin/env bash
# Semi-synthetic injection test, reports/stage2_inject_01/plan.md (+ amendment_01.md for --device cuda).
#
#   cd /home/liu/Research/isoform_commot
#   # 0) region input (once; same pair list as stage2_input_post is checked below)
#   python stage2.py prepare --sample posterior --data data/visium_mouse_brain_sagittal_posterior \
#       --max-spots 400 --output stage2_input_post_roi400
#   # 1) smoke (2 seeds, oracle + permutation, raw, beta 0.1 and 1, 20000 iterations; not interpreted),
#   #    once per device; amendment_01.md section 5 picks the main-run device from ot_wall_seconds
#   bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 \
#       --check-lr stage2_input_post --runs analysis_inject_smoke --device cuda --gpu-env inject-gpu --jobs 4
#   bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 \
#       --check-lr stage2_input_post --runs analysis_inject_smoke_cpu --device cpu --jobs 24
#   # 2) main (plan.md and, for cuda, amendment_01.md must be committed and unchanged)
#   bash reports/stage2_inject_01/run_inject.sh --input stage2_input_post_roi400 --check-lr stage2_input_post \
#       --runs analysis_inject_01 --results reports/stage2_inject_01/summary --device cuda --gpu-env inject-gpu \
#       --jobs 4 --v12-report reports/stage2_ppi_01/scores/score_report.json
#
# Environments: commot-research for tests, truth, scores and summaries (CPU). OT runs use commot-research with
# --device cpu (official solver), or --gpu-env (a clone of commot-research plus a CUDA build of torch) with
# --device cuda (gpu_unot.py port; amendment_01.md). For cuda the script first checks the GPU env, runs
# test_gpu_unot on the GPU, and runs the equivalence gate (gpu_unot.py equivalence, max relative difference
# <= 1e-9 against the official solver on this input); any failure stops the script before the OT runs.
# --runs is git-ignored by the analysis*/ rule; --results is committed (do not name it results*/).
set -euo pipefail

PLAN=reports/stage2_inject_01/plan.md
AMEND=reports/stage2_inject_01/amendment_01.md
BENCH_ENV=commot-research
GPU_ENV=inject-gpu
DEVICE=cpu
EQ_ITERS=500
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
    --device) DEVICE=$2; shift 2;;
    --gpu-env) GPU_ENV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[[ -n "$INPUT" ]] || { echo "--input is required" >&2; exit 2; }
[[ -n "$RUNS" ]] || { echo "--runs is required" >&2; exit 2; }
[[ "$DEVICE" == cpu || "$DEVICE" == cuda ]] || { echo "--device must be cpu or cuda" >&2; exit 2; }
[[ -e "$RUNS" ]] && { echo "$RUNS exists; use a new --runs" >&2; exit 1; }
if [[ -n "$CHECK_LR" ]]; then
  cmp -s "$INPUT/lr.csv" "$CHECK_LR/lr.csv" || { echo "$INPUT/lr.csv differs from $CHECK_LR/lr.csv" >&2; exit 1; }
fi
committed() {   # file is tracked and has no uncommitted changes
  git ls-files --error-unmatch "$1" >/dev/null 2>&1 || { echo "$1 is not committed" >&2; exit 1; }
  git diff --quiet HEAD -- "$1" || { echo "$1 has uncommitted changes" >&2; exit 1; }
}
if [[ $SMOKE -eq 1 ]]; then
  RES="$RUNS/summary"; SEEDS=(1 2); VARIANTS=(q100 q100_perm); NORMS=(raw); BETAS="0.1 1"; NITER=20000
  SMOKE_FLAG=(--allow-unconverged)
else
  [[ -n "$RES" ]] || { echo "--results is required for the main run" >&2; exit 2; }
  [[ -e "$RES" ]] && { echo "$RES exists; use a new --results" >&2; exit 1; }
  committed "$PLAN"
  [[ "$DEVICE" == cuda ]] && committed "$AMEND"
  SEEDS=(1 2 3); VARIANTS=(q100 q100_perm q90 q90_perm q80 q80_perm q70 q70_perm q60 q60_perm)
  NORMS=(raw ligand); BETAS="0.05 0.1 0.2 0.5 1"; NITER=200000
  SMOKE_FLAG=()
fi
if [[ "$DEVICE" == cuda ]]; then RUN_ENV=$GPU_ENV; else RUN_ENV=$BENCH_ENV; fi
setting() { case "$1" in main) echo "--eps 0.1 --rho 10";; sens) echo "--eps 0.2 --rho 1";; esac; }
PLAN_COMMIT=$(git log -1 --format=%H -- "$PLAN" 2>/dev/null || true); PLAN_COMMIT=${PLAN_COMMIT:-none}
AMEND_COMMIT=$(git log -1 --format=%H -- "$AMEND" 2>/dev/null || true); AMEND_COMMIT=${AMEND_COMMIT:-none}
mkdir -p "$RUNS/logs"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[inject $(date -u +%H:%M:%S)] $*" | tee -a "$RUNS/run.log"; }
START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
log "device=$DEVICE run_env=$RUN_ENV jobs=$JOBS smoke=$SMOKE"

log "1/6 unit tests"
run "$BENCH_ENV" python -m unittest test_stage2_inject 2>&1 | tee "$RUNS/unittest.log" | tail -3
if [[ "$DEVICE" == cuda ]]; then
  run "$GPU_ENV" python -c "import torch, scipy, pandas, matplotlib, numpy; \
assert torch.cuda.is_available(), 'CUDA not available in $GPU_ENV'; \
print('torch', torch.__version__, 'cuda', torch.version.cuda, 'gpu', torch.cuda.get_device_name(0), \
'numpy', numpy.__version__, 'scipy', scipy.__version__, 'pandas', pandas.__version__)" | tee "$RUNS/gpu_env.txt"
  run "$GPU_ENV" python -m unittest test_stage2_inject 2>&1 | tee "$RUNS/unittest_gpu_env.log" | tail -3
  run "$GPU_ENV" python -m unittest -v test_gpu_unot 2>&1 | tee "$RUNS/unittest_gpu.log" | tail -3
  grep -qi "skip" "$RUNS/unittest_gpu.log" && { echo "test_gpu_unot was skipped in $GPU_ENV" >&2; exit 1; }
fi

log "2/6 truth and scores per seed"
for s in "${SEEDS[@]}"; do
  run "$BENCH_ENV" python stage2_inject.py truth --input "$INPUT" --output "$RUNS/truth_s$s" --seed "$s"
  run "$BENCH_ENV" python stage2_inject.py scores --truth "$RUNS/truth_s$s" --output "$RUNS/scores_s$s" --seed "$s" \
      > "$RUNS/logs/scores_s$s.log"
done

EQ="$RUNS/equivalence"
if [[ "$DEVICE" == cuda ]]; then
  log "3/6 GPU equivalence gate ($EQ_ITERS fixed iterations, seed 1, tolerance 1e-9)"
  mkdir -p "$EQ"
  for st in main sens; do
    run "$GPU_ENV" python gpu_unot.py equivalence --input "$INPUT" --variant none $(setting $st) \
        --iters $EQ_ITERS --device cuda --output "$EQ/${st}_none.json" > "$RUNS/logs/eq_${st}_none.log" 2>&1 \
      || { echo "GPU equivalence gate failed: $st none (see $RUNS/logs/eq_${st}_none.log)" >&2; exit 1; }
    run "$GPU_ENV" python gpu_unot.py equivalence --input "$INPUT" --scores "$RUNS/scores_s1" --variant q100 \
        --norm raw --beta 1 $(setting $st) --iters $EQ_ITERS --device cuda --output "$EQ/${st}_q100_raw_b1.json" \
        > "$RUNS/logs/eq_${st}_q100.log" 2>&1 \
      || { echo "GPU equivalence gate failed: $st q100 (see $RUNS/logs/eq_${st}_q100.log)" >&2; exit 1; }
  done
  run "$GPU_ENV" python -c "import json, glob; d = [json.load(open(f)) for f in sorted(glob.glob('$EQ/*.json'))]; \
print('equivalence', len(d), 'cases, all passed:', all(x['passed'] for x in d), \
'max_rel_diff:', max(x['max_rel_diff'] for x in d), 'entries:', d[0]['entries'], \
'speedup per process (CPU s / GPU s):', [x['speedup'] for x in d])" | tee -a "$RUNS/run.log"
else
  log "3/6 GPU equivalence gate: not needed (device cpu)"
fi

log "4/6 OT runs (${JOBS} in parallel, $DEVICE)"
CMDS="$RUNS/commands.txt"; : > "$CMDS"
for st in main sens; do
  for s in "${SEEDS[@]}"; do
    common="python stage2_inject.py run --input $INPUT --truth $RUNS/truth_s$s --scores $RUNS/scores_s$s --seed $s $(setting $st) --nitermax $NITER --device $DEVICE"
    echo "$CONDA run --no-capture-output -n $RUN_ENV $common --variant none --output $RUNS/$st/s${s}_none > $RUNS/logs/${st}_s${s}_none.log 2>&1" >> "$CMDS"
    echo "$CONDA run --no-capture-output -n $RUN_ENV $common --variant uniform --beta $BETAS --output $RUNS/$st/s${s}_uniform > $RUNS/logs/${st}_s${s}_uniform.log 2>&1" >> "$CMDS"
    for v in "${VARIANTS[@]}"; do
      for n in "${NORMS[@]}"; do
        echo "$CONDA run --no-capture-output -n $RUN_ENV $common --variant $v --norm $n --beta $BETAS --output $RUNS/$st/s${s}_${v}_${n} > $RUNS/logs/${st}_s${s}_${v}_${n}.log 2>&1" >> "$CMDS"
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
OT_START=$(date +%s)
ls "$RUNS"/cmds/*.sh | xargs -P "$JOBS" -n 1 bash || { echo "some runs failed; see $RUNS/logs" >&2; exit 1; }
OT_WALL=$(( $(date +%s) - OT_START ))
log "   OT runs took ${OT_WALL} s wall time ($DEVICE, ${JOBS} in parallel)"

log "5/6 summaries"
mkdir -p "$RES"
for st in main sens; do
  run "$BENCH_ENV" python stage2_inject.py summarise --runs "$RUNS/$st" --output "$RES/$st" \
      ${V12:+--v12-report "$V12"} ${SMOKE_FLAG[@]+"${SMOKE_FLAG[@]}"} | tee "$RUNS/logs/summary_$st.log"
done

log "6/6 provenance"
mkdir -p "$RES/seeds"
for s in "${SEEDS[@]}"; do
  cp "$RUNS/truth_s$s/truth.json" "$RES/seeds/truth_s$s.json"
  cp "$RUNS/truth_s$s/pairs.csv" "$RES/seeds/pairs_s$s.csv"
  cp "$RUNS/scores_s$s/scores.json" "$RES/seeds/scores_s$s.json"
done
if [[ "$DEVICE" == cuda ]]; then
  mkdir -p "$RES/equivalence"; cp "$EQ"/*.json "$RES/equivalence/"; cp "$RUNS/gpu_env.txt" "$RES/gpu_env.txt"
fi
cat > "$RES/run_info.json" <<EOF
{
  "plan": "$PLAN",
  "plan_commit": "$PLAN_COMMIT",
  "amendment": "$( [[ "$DEVICE" == cuda ]] && echo "$AMEND" || echo "not used (device cpu)" )",
  "amendment_commit": "$( [[ "$DEVICE" == cuda ]] && echo "$AMEND_COMMIT" || echo none )",
  "head_commit_at_start": "$(git rev-parse HEAD 2>/dev/null || echo none)",
  "smoke": $SMOKE,
  "device": "$DEVICE",
  "run_env": "$RUN_ENV",
  "bench_env": "$BENCH_ENV",
  "equivalence_iterations": $( [[ "$DEVICE" == cuda ]] && echo $EQ_ITERS || echo 0 ),
  "input": "$INPUT",
  "check_lr": "$CHECK_LR",
  "v12_report": "$V12",
  "seeds": "${SEEDS[*]}",
  "betas": "$BETAS",
  "nitermax": $NITER,
  "run_started_utc": "$START",
  "run_finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "runs_dir": "$RUNS",
  "parallel_jobs": $JOBS,
  "ot_wall_seconds": $OT_WALL
}
EOF
log "done -> $RES"
