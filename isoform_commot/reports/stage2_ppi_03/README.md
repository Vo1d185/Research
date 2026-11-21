# 实验二扩大版:COMMOT 论文的数据集

计划:[plan.md](plan.md)(先提交,再运行)。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../exp2x_data.py` | 下载(Zenodo、CellChatDB 人、scSeqComm 人)、A 线输入、B 线输入(论文流程) | commot-research |
| `../../exp2x.py` | A 线 cases / run / eval-commot / regression / summarise;B 线 paper-ref / run-b / eval-b / summarise-b | commot-research;run-b 用 inject-gpu |
| `../../test_exp2x.py` | 单元测试(`python -m unittest test_exp2x`) | commot-research |
| `run_prepare_x.sh` | 第 1 步:数据与 humanV12 打分(不跑 OT) | commot-research + plminteract |
| `run_x.sh` | 第 2 步:门槛、A 线、参照、B 线、汇总 | commot-research + inject-gpu |

已有文件只调用,没有改动:`stage2.py`、`stage2_s.py`、`stage2_ppi.py`、`proteins.py`、`plm_score.py`、`evaluate_targets.py`、
`commot_baseline.py`、`gpu_unot.py`。

所有路径都是命令行参数(默认值写在两个脚本开头);输出已存在时拒绝运行。

## 给 VS Code 插件的指令

**第 1 步**

> 读 isoform_commot/reports/stage2_ppi_03/plan.md,确认它已提交,记下提交号。
> 在 isoform_commot/ 下按第 6 节第 2 项执行 run_prepare_x.sh:先冒烟,检查通过条件,再正式运行。
> 不修改 plan.md 和代码;出错就停下,报告错误输出和对应的日志文件。
> 完成后提交 reports/stage2_ppi_03/scores/。汇报:
> - scores/data_report.json 中每个输入的 pairs_selected 或 pairs、receptors、paper_input_x_max_abs_diff、paper_eval.n_hvg_matches;
> - 每个打分组的 n_pairs、n_scored、n_missing、unresolved_genes;
> - 各步用时和提交号。

**第 2 步**(第 1 步提交之后)

> 按 plan.md 第 6 节第 3 项执行 run_x.sh:先冒烟,检查通过条件,再正式运行(中断后加 --resume)。
> 不修改 plan.md 和代码;门槛失败就停下并报告。
> 完成后提交 reports/stage2_ppi_03/summary/。汇报:
> - gates/ 下三个门槛的结果;
> - a_line/decisions.json 中 pooled.deciding 的 D1/D2/D3 标签与数值;
> - 每套数据的 receptors_common 和 M;
> - b_line/b_line.json 每套数据的 metric 表和 reproduction(step0a、step0b);
> - run_info.json 的用时和提交号。
