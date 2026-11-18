# Stage 2 加入 PLM-interact 分数：第一步

计划：[plan.md](plan.md)。本步只打分和描述，不跑 OT。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../stage2_ppi.py` | genes / jobs / import / designs | commot-research |
| `../../test_stage2_ppi.py` | 单元测试（`python -m unittest test_stage2_ppi`） | commot-research |
| `../../proteins.py map` | 基因 → 小鼠蛋白（已有，未改动） | commot-research |
| `../../plm_score.py` | 打分（已有，未改动） | plminteract |
| `run_scoring.sh` | 冒烟（`--smoke`）与正式运行；正式运行检查 plan.md 已提交 | 两个都用 |

所有路径都是命令行参数；输出已存在时拒绝运行。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/stage2_ppi_01/plan.md，确认它已提交，记下提交号。
> 按第 5 节执行：先跑冒烟测试，通过后再跑正式运行，参数照第 5 节。
> 不修改 plan.md 和代码；出错就停下并报告错误输出。
> 完成后提交 reports/stage2_ppi_01/scores/，汇报 score_report.json 中的 n_pairs、n_scored、n_missing、score_distribution、panel_vs_unannotated，
> 以及 proteins.json 里的 unresolved 列表和提交号。
