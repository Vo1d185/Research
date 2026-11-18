# 注入实验 01 · 修订 01:OT 求解可以在 GPU 上运行

**状态:** 2026-11-18 起草,在任何冒烟或正式运行之前提交。`plan.md` 一字不改。
`run_inject.sh` 在 `--device cuda` 的正式运行中会检查本文件已提交且工作区无改动。

## 1. 改了什么

只改 OT 求解器运行在哪个设备上:

| 项目 | plan.md(提交 accf094) | 本修订 |
|---|---|---|
| 求解器 | `research.solve_sparse`(COMMOT `unot_sinkhorn_l1_sparse`,CPU) | `--device cpu`:不变;`--device cuda`:`gpu_unot.py::solve_sparse_gpu`(同一算法的 torch 实现) |
| OT 运行所在环境 | `commot-research` | `--device cuda` 时用新环境 `inject-gpu`(`commot-research` 的克隆加 CUDA 版 torch) |
| truth / scores / summarise | `commot-research` | 不变 |

**以下全部不变:** 真值规则、合成分数、10 个分数变体、两种归一化、β 网格、ε/ρ 两种设置、r₀ = 50 µm、截断 363.08 µm、
stopthr = 1e-11、nitermax = 200,000、种子 1–3、判定标准 R1–R5、汇总代码、冒烟与正式运行的内容。

## 2. 为什么是同一个算法

1. 迭代逐行照抄 `vendor/_unot.py::unot_sinkhorn_l1_sparse`(COMMOT 提交 d117445b):f、g 的更新顺序相同;
   停止规则相同(每 10 次迭代检查一次相对变化);迭代计数的含义相同(大于 nitermax 表示到达上限);全程 float64。
2. 外层 `solve_sparse_gpu` 与 `research.solve_sparse` 一样:同样的输入检查、只保留有预算的行列、同一个缩放因子。
   结果用**同一个** `sparse_mass_check` 和 `sparse_kkt_check` 判断是否收敛,所以"未收敛"的判定标准也相同。
3. 唯一的差别是浮点求和的顺序(GPU 并行累加)。这会带来约 1e-15 量级的舍入差;GPU 上同一设置重复运行也可能有这个量级的差别。
   极少数情况下,停止判定可能提前或推后一个检查点(10 次迭代);在 stopthr = 1e-11 下,这对指标的影响远小于汇总表的 5 位小数。
4. 不改用 float32。RTX 6000 Ada 的 float64 速度只有 float32 的约 1/64,所以 GPU 不一定更快(见第 5 节)。

## 3. 等价性门槛(运行前登记,不可放宽)

`--device cuda` 时,`run_inject.sh` 在任何正式 OT 运行之前执行下面三项,任何一项不通过就停止:

1. 环境预检:`inject-gpu` 中 `torch.cuda.is_available()` 为真,并记录 torch、CUDA、GPU 型号、numpy、scipy、pandas 版本(`gpu_env.txt`)。
2. `test_gpu_unot` 在 GPU 上运行,不允许被跳过:随机小问题,固定 300 次迭代,与官方求解器的相对差 < 1e-9;
   完整求解时迭代次数相同、是否到达上限相同、运输矩阵相对差 < 1e-9。`test_stage2_inject` 也在 `inject-gpu` 中再跑一遍。
3. 实际输入上的对照:种子 1,两种设置 × {`none`,`q100` raw β = 1},共 4 个问题。官方 CPU 求解器和 GPU 版本各跑固定 500 次迭代
   (stopthr = 0),要求条目顺序一致,运输矩阵最大相对差 ≤ 1e-9(`gpu_unot.py equivalence`)。

不通过就回到 `--device cpu`,按 plan.md 原样运行。门槛结果(`equivalence/*.json`、`gpu_env.txt`)随汇总一起提交。

## 4. 运行规则

1. 一次汇总只能来自一种设备。`summarise` 会检查,混用就报错。正式运行要么全部在 GPU,要么全部在 CPU。
2. 如果 GPU 正式运行中途失败(例如显存不足),换一个新的 `--runs` 目录,**全部**重新运行(GPU 或 CPU 都可以),不拼接两次运行的结果。
3. 设备只能根据冒烟测试的耗时选择,不能根据任何指标的数值选择。两种设备给出的数值在 1e-9 以内相同,所以选哪种不影响结论。
4. `run_info.json` 记录 device、运行环境、plan 与本修订的提交号、等价门槛的迭代次数。

## 5. 用哪种设备跑正式运行(按耗时预先规定)

GPU 不一定更快:float64 在这张卡上慢,而 plan.md 的 CPU 方案是 24 个进程并行。所以冒烟测试**两种设备各跑一次**,内容完全相同:

1. GPU 冒烟:`--device cuda --gpu-env inject-gpu --jobs 4`,目录 `analysis_inject_smoke`;
2. CPU 冒烟:`--device cpu --jobs 24`,目录 `analysis_inject_smoke_cpu`。

比较两次 `summary/run_info.json` 中的 `ot_wall_seconds`(第 4 步 OT 运行的总耗时)。
冒烟只有 16 条运行命令,CPU 在正式运行中可以用满 24 个并行,所以规则是:
**GPU 耗时 ≤ CPU 耗时 × 2/3 时正式运行用 GPU,否则用 CPU**(相等或接近时用 plan.md 原定的 CPU)。
GPU 冒烟没通过等价门槛时直接用 CPU。等价门槛 `equivalence/*.json` 中的 `cpu_seconds`、`port_seconds`、`speedup`
(单进程各跑 500 次迭代)只作参考,不参与这条规则。

## 6. 本地验证(已完成;不是实验结果)

在没有 GPU 的机器上,用同一份 torch 代码在 CPU 上运行(测试专用设备名 `torch-cpu`),与官方求解器对照:

- `test_gpu_unot` 3 项、`test_stage2_inject` 10 项、`test_stage2_ppi` 7 项全部通过(有 torch 与无 torch 的环境都测过);
- 玩具输入上的完整冒烟流程(GPU 分支):等价门槛 4 个问题全部通过,最大相对差 8.3e-15;
- 同一冒烟分别用官方求解器和 torch 版本运行,28 次求解逐行比较:迭代次数完全相同(包括哪 13 次到达 20,000 上限),
  指标最大差 6e-14,`decisions.json` 的所有判定字段相同。

GPU 上的实际数值和耗时由服务器上的冒烟测试给出。

## 7. 环境(新建,不改动现有环境)

```bash
conda create -n inject-gpu --clone commot-research -y
conda run -n inject-gpu pip install torch --index-url https://download.pytorch.org/whl/cu121
conda run -n inject-gpu python -c "import torch, scipy, pandas; print(torch.__version__, torch.cuda.is_available())"
```

最后一行必须打印 `True`。torch 的具体版本由 `gpu_env.txt` 记录。`commot-research` 与 `plminteract` 保持原样。

## 8. 执行(替换 plan.md 第 7 节中的第 2、3 步命令;第 0、1、4 步不变)

冒烟(两条都跑,结果不解读):
```bash
bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 --check-lr stage2_input_post \
    --runs analysis_inject_smoke --device cuda --gpu-env inject-gpu --jobs 4
bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 --check-lr stage2_input_post \
    --runs analysis_inject_smoke_cpu --device cpu --jobs 24
```
按第 5 节选设备。正式(GPU):
```bash
bash reports/stage2_inject_01/run_inject.sh --input stage2_input_post_roi400 --check-lr stage2_input_post \
    --runs analysis_inject_01 --results reports/stage2_inject_01/summary --device cuda --gpu-env inject-gpu \
    --jobs 4 --v12-report reports/stage2_ppi_01/scores/score_report.json
```
正式(CPU):把上面的 `--device cuda --gpu-env inject-gpu --jobs 4` 换成 `--device cpu --jobs 24`。
