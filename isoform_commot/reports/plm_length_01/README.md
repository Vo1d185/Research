# 实验一:V12 在长、短蛋白配对上的区分能力

计划(预注册):[plan.md](plan.md)。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../plm_length.py` | fetch / prepare / work / merge / evaluate | commot-research;打分经 `--plm-python` 调用 PLM 环境 |
| `../../test_plm_length.py` | 单元测试(`python -m unittest test_plm_length`) | commot-research |
| `run_length.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md 已提交;`--resume` 续跑 | 见上 |

所有路径都是命令行参数;输出已存在时拒绝运行。运行目录用 `analysis_length_*`(被 `.gitignore` 忽略),
提交目录用 `summary`(不要用 `results*` 命名)。打分分成 74 块,每块独立完成;中断后加 `--resume`,已完成的块跳过。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/plm_length_01/plan.md,确认它已提交,记下提交号。
> 按第 7 节执行:先下载测试集(第 2 步),再跑冒烟(第 3 步),通过后跑正式运行(第 4 步)。
> 不修改 plan.md 和代码;出错就停下并报告错误输出。正式运行中断时,用同样的命令加 --resume 续跑。
> 完成后提交 reports/plm_length_01/summary/,汇报:
> summary/verdict.json 的全部内容;summary/dataset.json 中的 pairs、excluded_pair_seen、primary、e2、jobs;
> summary/run_info.json 的 scoring_wall_seconds;冒烟和正式运行的提交号。
