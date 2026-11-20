# 实验四:BindPred 的 S′ 进入 OT

计划(预注册):[plan.md](plan.md)。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../bindpred_ot.py` | jobs / import / cases / check / summarise | commot-research |
| `../../test_bindpred_ot.py` | 单元测试与端到端小例子(需要 `vendor/_unot.py`) | commot-research |
| `../../affinity.py` | embed(`--plm-python`)与 predict(commot-research,需 catboost),来自实验三,未修改 | 见左 |
| `../../stage2_s.py` | OT 求解与评价(`--setting main`),来自实验二,未修改 | commot-research |
| `run_bindpred_ot.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md 已提交;`--resume` 续跑 | 见上 |
| `inputs/mouse_panel_sequences.tsv` | 冻结的 384 条小鼠面板序列(UniProt 2026_04),sha256 与实验二一致 | — |

前提:实验三的一次性准备已完成(`models/bindpred/ESM2_BindPred.cbm`;commot-research 中已装 catboost、biopython)。
运行目录用 `analysis_bpot_*`(被 `.gitignore` 忽略),提交目录用 `summary`。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/bindpred_ot_01/plan.md,确认它已提交,记下提交号。
> 先跑冒烟,通过后跑正式运行(命令见 run_bindpred_ot.sh 开头的注释)。
> 不修改 plan.md 和代码;出错就停下并报告错误输出。复现门槛(第 6 步)不通过时停下汇报。
> 正式运行中断时,用同样的命令加 --resume 续跑。
> 完成后提交 reports/bindpred_ot_01/summary/,汇报:
> summary/decisions.json 的 verdict 全部内容;
> details.post 与 details.ant 中的 receptors_common、M、F1(p、perm_mean、perm_sd)、F2_vs_uniform、F3_vs_no_ppi、
> uniform_vs_no_ppi、sensitivity 全部、top_tf_5、top_tf_20;
> summary/scores/gate.json、import.json、jobs.json 的内容;
> run_info.json 的 embed_seconds、ot_wall_seconds;各提交号。
