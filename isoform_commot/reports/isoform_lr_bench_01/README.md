# isoform 排序基准 01

计划（预注册）：[plan.md](plan.md)。数据：[dataset/](dataset/)（冻结，sha256 见 plan.md 第 1 节）。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../isoform_lr_bench.py` | jobs / merge / evaluate / selftest | commot-research |
| `../../test_isoform_lr_bench.py` | 单元测试（`python -m unittest test_isoform_lr_bench`） | commot-research |
| `../../plm_score.py` | 打分（已有，未改动） | plminteract |
| `smoke_test.sh` | 冒烟测试：52 个作业，含最长一对 | 两个都用 |
| `run_main.sh` | 正式运行；检查 plan.md 已提交 | 两个都用 |

所有路径都是命令行参数；输出目录已存在时拒绝运行。环境名可用 `--bench-env`、`--plm-env` 覆盖；
`conda` 不在 PATH 时设置 `CONDA_EXE`。

## 在服务器上

```bash
cd /home/liu/Research/isoform_commot
# 1) 冒烟
bash reports/isoform_lr_bench_01/smoke_test.sh \
  --out isoform_lr_bench_runs/smoke_01 \
  --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
  --base-model models/esm2_t33_650M_UR50D
# 2) 正式（plan.md 必须已提交）
bash reports/isoform_lr_bench_01/run_main.sh \
  --runs isoform_lr_bench_runs/main_01 \
  --results reports/isoform_lr_bench_01/results_humanV12 \
  --checkpoint models/plm_650m_humanV12/pytorch_model.bin \
  --base-model models/esm2_t33_650M_UR50D \
  --training-root data/lr_benchmark/plm_training
```

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/isoform_lr_bench_01/plan.md，确认它已提交，记下提交号。
> 按第 7 节执行：先跑 smoke_test.sh，通过（第 9 节条件）后再跑 run_main.sh，参数照第 7 节。
> 不修改 plan.md 和代码；出错就停下并报告错误输出。
> 完成后提交 reports/isoform_lr_bench_01/results_humanV12/，汇报 verdict.json 内容和提交号。
