# 注入实验 01

计划(预注册):[plan.md](plan.md);修订 01(OT 求解可在 GPU 上运行,算法不变):[amendment_01.md](amendment_01.md)。
不用 PDE、不用真实 PLM 分数、不看真实靶基因评价。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../stage2_inject.py` | truth / scores / run / summarise;`run --device cpu` 或 `cuda` | commot-research(cuda 时 run 在 inject-gpu) |
| `../../gpu_unot.py` | COMMOT 求解器的 torch 实现;`equivalence` 子命令与官方求解器逐项对照 | inject-gpu |
| `../../test_stage2_inject.py` | 单元测试(`python -m unittest test_stage2_inject`;端到端小例子需先运行 `fetch_solver.py`) | commot-research |
| `../../test_gpu_unot.py` | torch 实现与官方求解器一致(GPU 可用时也在 cuda 上测) | inject-gpu |
| `run_inject.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md(cuda 时还有 amendment_01.md)已提交;cuda 时先过等价门槛 | 见上 |

所有路径都是命令行参数;输出已存在时拒绝运行。运行目录用 `analysis_inject_*`(被 `.gitignore` 忽略),提交目录用 `summary`(不要用 `results*` 命名)。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/stage2_inject_01/plan.md 和 amendment_01.md,确认两者都已提交,记下提交号。
> 按 plan.md 第 7 节执行,其中第 2、3 步换成 amendment_01.md 第 8 节的命令。
> 先准备区域输入(第 0 步);再跑两次冒烟(GPU 一次、CPU 一次);按 amendment_01.md 第 5 节的耗时规则选设备;再跑正式运行。
> 不修改 plan.md、amendment_01.md 和代码;出错就停下并报告错误输出。GPU 冒烟的等价门槛不通过时停下汇报。
> 完成后提交 reports/stage2_inject_01/summary/,汇报:
> 两次冒烟 run_info.json 的 ot_wall_seconds,以及按第 5 节选了哪种设备;
> GPU 冒烟 run.log 中以 "equivalence" 开头的那一行(max_rel_diff、entries、speedup)和 gpu_env.txt 的内容;
> summary/main/decisions.json 和 summary/sens/decisions.json 中的 device、R1_oracle_usable、R2_best_beta_oracle_raw、
> R3_threshold_auroc、R4_v12、R5_receptor_metric_sees_oracle、runs_not_converged;
> 正式运行的 ot_wall_seconds;各提交号。
