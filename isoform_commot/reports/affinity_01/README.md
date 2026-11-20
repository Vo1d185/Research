# 实验三:BindPred 能否给出可用作 S 的结合强度

计划(预注册):[plan.md](plan.md)。本实验只打分、不做 OT。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../affinity.py` | prepare / embed / predict / evaluate | prepare、evaluate 用 commot-research(需 biopython);embed 用 `--plm-python`;predict 用 `--predict-env`(需 catboost) |
| `../../test_affinity.py` | 单元测试(`python -m unittest test_affinity`,需 catboost) | commot-research |
| `run_affinity.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md 已提交;`--resume` 续跑 | 见上 |
| `inputs/` | 冻结的输入:GtoPdb 2026.3.1 的 interactions 与 peptides、UniProt 2026_04 序列;sha256 记在 `dataset.json` | — |

所有路径都是命令行参数;输出已存在时拒绝运行。运行目录用 `analysis_affinity_*`(被 `.gitignore` 忽略),
提交目录用 `summary`。

## 一次性准备

```bash
cd /home/liu/Research/isoform_commot
mkdir -p models/bindpred && cd models/bindpred
R=89a595d1c161a95783d904f5c5b3cce565e52848
for f in ESM2_BindPred.cbm embeddings_only.csv; do
  curl -L -o $f "https://huggingface.co/hbp5181/BindPred/resolve/$R/$f"
done
sha256sum ESM2_BindPred.cbm embeddings_only.csv
# 预期:f56ac4a1d63149d6e870f89621c2a63c46070d8b9725502230a14167afaeb47f  ESM2_BindPred.cbm
#       0896c9ee8eca42efba6340e1727c2a6c6cc96da1c86f232d76ffe905c75cd047  embeddings_only.csv
cd /home/liu/Research/isoform_commot
conda run -n commot-research pip install catboost biopython     # 纯 python 轮子,不动已有依赖
```

`embeddings_only.csv` 为 639 MB,只读其中两列序列用于查泄漏;放在 `models/`(已被 .gitignore 忽略)。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/affinity_01/plan.md,确认它已提交,记下提交号。
> 按 README 的"一次性准备"下载 BindPred 权重与训练表,核对两个 sha256;不一致就停下汇报。
> 先跑 run_affinity.sh 的冒烟(--smoke),通过后跑正式运行。
> 不修改 plan.md 和代码;出错就停下并报告错误输出。
> 完成后提交 reports/affinity_01/summary/,汇报:
> summary/decisions.json 中 primary_receptor_groups 的全部内容(rho_bar、ci95、paired_diff_ci95、A1/A2/A3、verdict);
> secondary_ligand_groups 的 rho_bar 与 verdict;pooled 全部内容;
> strata 中各层的 n_groups 与 rho_bar;group_summary;
> data/dataset.json 的 counts;embed 与 predict 的用时;各提交号。
