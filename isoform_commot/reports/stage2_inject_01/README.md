# 注入实验 01

计划(预注册):[plan.md](plan.md);修订 01(OT 求解可在 GPU 上运行,算法不变):[amendment_01.md](amendment_01.md);
修订 02(主设置的收敛标准改用 Stage 2 第 8.7 节,外加稳定性诊断):[amendment_02.md](amendment_02.md)。
不用 PDE、不用真实 PLM 分数、不看真实靶基因评价。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../stage2_inject.py` | truth / scores / run / summarise / stability / recheck-commands;`run --device cpu` 或 `cuda`;`summarise --convergence cap` 或 `checks` | commot-research(cuda 时 run 在 inject-gpu) |
| `../../gpu_unot.py` | COMMOT 求解器的 torch 实现;`equivalence` 子命令与官方求解器逐项对照 | inject-gpu |
| `../../test_stage2_inject.py` | 单元测试(`python -m unittest test_stage2_inject`;端到端小例子需先运行 `fetch_solver.py`) | commot-research |
| `../../test_gpu_unot.py` | torch 实现与官方求解器一致(GPU 可用时也在 cuda 上测) | inject-gpu |
| `run_amend02.sh` | 修订 02 的收尾:用 400,000 次上限重跑 8 条稳定性运行;汇总(主设置 `--convergence checks`,敏感性按原规则);稳定性比较;溯源 | commot-research(复查运行沿用原命令的环境) |
| `run_inject.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md(cuda 时还有 amendment_01.md)已提交;cuda 时先过等价门槛 | 见上 |

所有路径都是命令行参数;输出已存在时拒绝运行。运行目录用 `analysis_inject_*`(被 `.gitignore` 忽略),提交目录用 `summary`(不要用 `results*` 命名)。

## 给 VS Code 插件的指令

修订 02 的收尾(正式 OT 已在 11-19 跑完,不重跑):

> 读 isoform_commot/reports/stage2_inject_01/amendment_02.md,确认 plan.md、amendment_01.md、amendment_02.md 都已提交,记下提交号。
> 按 amendment_02.md 第 5 节执行:第 0 步移走残留目录时不要打开里面的文件;先跑冒烟,通过后再跑正式。
> 不修改 plan.md、两份修订和代码;出错就停下并报告错误输出。
> 完成后提交 reports/stage2_inject_01/summary/,汇报:
> summary/stability/stability_main.json 的 compared_solves、failed_checks_used、failed_checks_long、max_abs_diff、passed、verdict;
> summary/main/decisions.json 和 summary/sens/decisions.json 中的 convergence_rule、runs_at_cap、runs_not_converged、R1_oracle_usable、
> R2_best_beta_oracle_raw、R3_threshold_auroc、R4_v12、R5_receptor_metric_sees_oracle;
> summary/run_info.json 的 check_ot_wall_seconds;各提交号。

修订 01 时的原指令(已执行,留作记录):

> 读 plan.md 和 amendment_01.md,按 plan.md 第 7 节执行,第 2、3 步换成 amendment_01.md 第 8 节的命令;两次冒烟后按第 5 节选设备,再跑正式运行。
