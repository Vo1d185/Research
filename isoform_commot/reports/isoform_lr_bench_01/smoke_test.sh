#!/usr/bin/env bash
# Smoke test for the isoform ranking benchmark (plan.md section 9). Not interpreted.
#
#   cd /home/liu/Research/isoform_commot
#   bash reports/isoform_lr_bench_01/smoke_test.sh \
#       --out isoform_lr_bench_runs/smoke_01 \
#       --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
#       --base-model models/esm2_t33_650M_UR50D
#
# Environments: commot-research (unit tests, jobs, merge, evaluate);
#               plminteract (plm_score.py, GPU). Override with --bench-env / --plm-env.
# Size: 3 contrasts per tier + the contrast holding the longest pair (2790 tokens), 52 jobs.
set -euo pipefail

DATASET=reports/isoform_lr_bench_01/dataset/benchmark.csv.gz
BENCH_ENV=commot-research
PLM_ENV=plminteract
CKPT_NAME=danliu1226/PLM-interact-650M-humanV12
CKPT_REV=2d0e7fb3deef94efa3f217daccc79092a6a42ceb
OUT= ; CKPT= ; BASE= ; CONDA=${CONDA_EXE:-conda}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT=$2; shift 2;;
    --checkpoint) CKPT=$2; shift 2;;
    --base-model) BASE=$2; shift 2;;
    --dataset) DATASET=$2; shift 2;;
    --bench-env) BENCH_ENV=$2; shift 2;;
    --plm-env) PLM_ENV=$2; shift 2;;
    --checkpoint-name) CKPT_NAME=$2; shift 2;;
    --checkpoint-revision) CKPT_REV=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
: "${OUT:?--out is required}"; : "${CKPT:?--checkpoint is required}"; : "${BASE:?--base-model is required}"
[[ -e "$OUT" ]] && { echo "$OUT exists; use a new --out" >&2; exit 1; }
mkdir -p "$OUT"
run() { local env=$1; shift; "$CONDA" run --no-capture-output -n "$env" "$@"; }
log() { echo "[smoke $(date -u +%H:%M:%S)] $*"; }

log "1/6 unit tests ($BENCH_ENV)"
run "$BENCH_ENV" python -m unittest test_isoform_lr_bench 2>&1 | tee "$OUT/unittest.log" | tail -3

log "2/6 decision-rule selftest on the frozen dataset (CPU)"
run "$BENCH_ENV" python isoform_lr_bench.py selftest --dataset "$DATASET" --output-dir "$OUT/selftest"

log "3/6 jobs (3 contrasts per tier + longest pair)"
run "$BENCH_ENV" python isoform_lr_bench.py jobs --dataset "$DATASET" --output-dir "$OUT/jobs" \
    --limit-contrasts 3 --include-longest --seed 0 > "$OUT/jobs.log"

log "4/6 score with plm_score.py ($PLM_ENV, full length, no truncation)"
run "$PLM_ENV" python plm_score.py score --jobs "$OUT/jobs/plm_input.csv" \
    --checkpoint "$CKPT" --base-model "$BASE" \
    --checkpoint-name "$CKPT_NAME" --checkpoint-revision "$CKPT_REV" \
    --output "$OUT/scores.csv" --allow-long

log "5/6 merge"
run "$BENCH_ENV" python isoform_lr_bench.py merge --dataset "$DATASET" --jobs-dir "$OUT/jobs" \
    --scores "$OUT/scores.csv" --output "$OUT/merged.csv"

log "6/6 evaluate (exploratory: smoke subset, no training overlap)"
run "$BENCH_ENV" python isoform_lr_bench.py evaluate --scored "$OUT/merged.csv" \
    --output-dir "$OUT/eval" --exploratory --n-perm 1000 --n-boot 100 > "$OUT/eval.log"

CHECK='
import json, sys, pathlib
out = pathlib.Path(sys.argv[1])
m = json.loads((out/"merged.csv.meta.json").read_text())
s = json.loads((out/"selftest"/"selftest.json").read_text())
j = json.loads((out/"jobs"/"jobs.json").read_text())
ok = dict(selftest_passed=s["passed"], n_jobs=j["n_jobs"], max_tokens=j["max_tokens"],
          scored_rows=m["n_rows"], checkpoint_sha256=m["checkpoint_sha256"],
          checkpoint_is_registered=m["checkpoint_sha256"] ==
              "92565937d65a9b91a4e897c99653c972b99cfe37067eb42d9e783af688350547",
          device=m.get("device"))
(out/"smoke_summary.json").write_text(json.dumps(ok, indent=2))
print(json.dumps(ok, indent=2))
if not (ok["selftest_passed"] and ok["checkpoint_is_registered"] and ok["scored_rows"] * 2 == ok["n_jobs"]):
    sys.exit("SMOKE FAILED")
print("SMOKE PASSED (scores are not interpreted)")
'
run "$BENCH_ENV" python -c "$CHECK" "$OUT"
