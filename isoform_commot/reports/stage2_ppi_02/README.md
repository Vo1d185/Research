# 实验二:真实切片上比较不加 S、原始 S 与老师的 S′

计划(预注册):[plan.md](plan.md)。

| 文件 | 作用 | 环境 |
|---|---|---|
| `../../stage2_s.py` | cases / run / equivalence / summarise | commot-research(主设置);inject-gpu + `--device cuda`(次设置) |
| `../../test_stage2_s.py` | 单元测试与端到端小例子(`python -m unittest test_stage2_s`;需要 `vendor/_unot.py`) | commot-research |
| `run_s.sh` | 冒烟(`--smoke`)与正式运行;正式运行检查 plan.md 已提交;`--resume` 续跑;两道等价门槛 | 见上 |

所有路径都是命令行参数;输出已存在时拒绝运行。运行目录用 `analysis_s_*`(被 `.gitignore` 忽略),提交目录用 `summary`。

## 给 VS Code 插件的指令

> 读 isoform_commot/reports/stage2_ppi_02/plan.md,确认它已提交,记下提交号。
> 按第 7 节执行:先跑冒烟(第 2 步),通过后跑正式运行(第 3 步)。
> 不修改 plan.md 和代码;出错就停下并报告错误输出。等价门槛不通过时停下汇报。正式运行中断时,用同样的命令加 --resume 续跑。
> 完成后提交 reports/stage2_ppi_02/summary/,汇报:
> summary/decisions.json 中 main 和 second 两部分的 verdict 全部内容;
> details.post 与 details.ant 中的 receptors_common、M、D1_raw、D1_snorm(p、perm_mean、perm_sd)、D2_raw、D2_snorm、D3_raw_vs_snorm;
> second 的 stability;summary/equivalence.json 与 gpu_equivalence.json 的结果;
> summary/run_info.json 的 main_wall_seconds、second_wall_seconds;冒烟时每个条件的大致用时;各提交号。
