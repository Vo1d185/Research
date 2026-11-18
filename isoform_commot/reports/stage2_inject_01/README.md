# 注入实验 01

计划（预注册）：[plan.md](plan.md)。不用 PDE、不用真实 PLM 分数、不看真实靶基因评价。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../stage2_inject.py` | truth / scores / run / summarise | commot-research |
| `../../test_stage2_inject.py` | 单元测试（`python -m unittest test_stage2_inject`；端到端小例子需先运行 `fetch_solver.py`） | commot-research |
| `run_inject.sh` | 冒烟（`--smoke`）与正式运行；正式运行会检查 plan.md 已提交 | commot-research |

所有路径都是命令行参数；输出已存在时拒绝运行。运行目录用 `analysis_inject_*`（被 `.gitignore` 忽略），提交目录用 `summary`（不要用 `results*` 命名）。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/stage2_inject_01/plan.md，确认它已提交，记下提交号。
> 按第 7 节执行：先准备区域输入（第 0 步），再跑冒烟测试，通过后跑正式运行，参数照第 7 节。
> 不修改 plan.md 和代码；出错就停下并报告错误输出。
> 完成后提交 reports/stage2_inject_01/summary/，汇报 summary/main/decisions.json 和 summary/sens/decisions.json 中的
> R1_oracle_usable、R2_best_beta_oracle_raw、R3_threshold_auroc、R4_v12、R5_receptor_metric_sees_oracle、runs_not_converged，
> 以及冒烟测试中每次求解的大致耗时和提交号。
