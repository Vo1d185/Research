# 注入实验 01：修改 cost 公式能否让 PPI 分数发挥作用

**状态：预注册计划。** 2026-11-18 起草。本文件提交后才运行正式实验；`run_inject.sh` 会检查本文件已提交且工作区无改动。
不使用 PDE，不使用真实 PLM 分数（这里的 S 是按已知质量合成的），也不查看任何真实靶基因评价。

代码：`stage2_inject.py`（truth / scores / run / summarise），`test_stage2_inject.py`，`reports/stage2_inject_01/run_inject.sh`。
沿用 `stage2.py` 的输入读取、spot 连接和稀疏 cost，以及 `research.solve_sparse`（固定版本的 COMMOT 求解器），均不改动。

## 0. 问题

把 cost 写成一族：**C = −log F(d) − β·log S_eff**。
- β = 0：不用 S；
- β = ε：S 以一次方进入 Gibbs 核；
- β = 1：原公式，S 以 1/ε 次方进入核。

S_eff 取原始 S（`raw`），或按配体归一化的 S′ = S / Σ_{该配体的受体} S（`ligand`）。

改公式不会给 S 增加信息，只会改变 S 在多大程度上左右运输结果。所以要回答的问题是：

- **Q1：** S 准确时，是否存在 β > 0 让推断更接近真值？也就是这个框架能不能用上 S。
- **Q2：** S 要多准才有用？用 AUROC 表示的门槛。
- **Q3：** V12 在本面板上的 AUROC（`reports/stage2_ppi_01/plan.md` 第 4.1 节）是否越过门槛？
- **Q4：** 受体层的指标（真实靶基因评价的对应物）能否看出 S 的作用？

## 1. 输入

`stage2_input_post_roi400`：`stage2.py prepare --max-spots 400` 取后部切片中心的 400 个 spot。
它的 LR 对清单必须与 `stage2_input_post` 完全相同，`run_inject.sh --check-lr` 会检查。
使用真实坐标和真实的配体、受体表达（即 OT 预算 a、b）。

## 2. 真值（已知答案；解析公式，不解方程）

- **诱饵：** 每种类型随机取 1/3 的配对设为诱饵（K = 0），代表数据库中实际不发生的配对。
- **结合常数：** 真实配对的 log₁₀K 相对值 ~ U(−1, 1)，跨两个数量级。整体尺度 κ 的取法是：在有受体表达的 spot 上，Σ K·L 的中位数为 1，即中等程度的饱和。
- **配体到达：** L_i(y) = Σ_x G(d)·a_i(x)。
  - 分泌型：G = exp(−d/λ)，λ = 100 µm，截断 363.08 µm，包括 d = 0；
  - 接触型：G = 1，范围为自身加六个邻居。
  - G 刻意不同于 OT 用的 F = r₀/(r₀ + d)，避免“用什么生成就用什么还原”。
- **竞争结合：** C\*_ij(y) = b_j(y)·K_ij·L_i(y) / (1 + Σ_i′ K_i′j·L_i′(y))，不考虑配体被消耗。
- **种子：** 1、2、3。每个种子独立抽取诱饵、K 和 S 的噪声。

## 3. 合成的 S（已知质量）

- 潜变量 u：真实配对 u = 1 + 0.5·z(log K)（z 截断在 ±2），诱饵 u = −1。
  S = sigmoid(4·(u + σ·η))，η ~ N(0, 1)。系数 4 让分数呈两极化，与 V12 类似。
- σ 用 AUROC 期望值的闭式公式二分求解，使 AUROC（真实配对 vs 诱饵）分别等于：1.0（`q100`，σ = 0，排序与真值完全一致）、0.9、0.8、0.7、0.6。
- 每个 S 另做一份同类型内打乱的版本（`_perm`）：分布不变，对应关系被破坏。
- 均匀对照 `uniform`：S′ = 1/k，k 为该配体的受体数。

## 4. OT 与指标

- **参数设置（都预注册，各自判定，以主设置为准）：** 主设置 ε = 0.1、ρ = 10（Stage 2 主分析）；敏感性设置 ε = 0.2、ρ = 1。
- **其余参数：** r₀ = 50 µm（Stage 2 选出的 r₀\*），截断 363.08 µm，stopthr = 1e-11，nitermax = 200,000。
  未收敛的运行记为“未收敛”，不参与判定，单独报告。
- **网格：** β ∈ {0.05, 0.1, 0.2, 0.5, 1}，归一化 ∈ {raw, ligand}。
  每种设置、每个种子 22 次运行（none 1 次、uniform 1 次、10 个 S 各 2 种归一化），共 132 次运行、636 次求解。
- **主指标——合并 Spearman：** 在所有（配对, spot）上（只取受体有表达的 spot，包括诱饵，诱饵的真值为 0），计算推断信号 r̂_p(y) = Σ_x T_p(x, y) 与真值 C\*_p(y) 的 Spearman 相关。
- **次要指标——受体层 Spearman：** 对每个受体，计算 Σ_p r̂ 与 Σ_p C\* 跨 spot 的 Spearman，再对受体取平均。这是真实靶基因评价的对应物。
- **描述性指标：** 配体分配误差（总变差距离），诱饵拿到的运输量占比。
- **增益：** 指标减去同种子 `none` 的指标。**相对打乱的增益：** 指标减去同种子、同（归一化, β）的打乱版本的指标。
- **不用“单个配对内部跨 spot 的 Spearman”：** S 对每个配对是常数，基本不改变配对内部的空间排序。

## 5. 判定（预注册；每种参数设置分别计算，以主设置为准）

| 编号 | 内容 |
|---|---|
| R1 框架能用上 S（正对照） | `q100` 在某个（归一化, β > 0）格子里，3 个种子都收敛，合并 Spearman 都高于 `none`，也都高于自己的打乱版本 |
| R2 β 的预测（描述性） | `q100`、`raw` 时，三个种子平均的最优 β 是否落在 {ε/2, ε, 2ε}。依据：真值复合物与 K 成正比，对应核中 S 的一次方，即 β = ε |
| R3 质量门槛 | 对每个质量水平：用种子 1 选出（归一化, β），取相对打乱增益最大者，并列时取较小的 β；种子 2、3 在同一格子上相对打乱的增益都 > 0，则判为“可用”。门槛 = 可用水平中最低的 AUROC；同时报告是否单调 |
| R4 V12 的定位 | 取 `stage2_ppi_01` 中面板配对 vs 未注释组合的 AUROC：95% 区间下限 ≥ 门槛为“越过”，上限 < 门槛为“未越过”，其余为“不确定” |
| R5 受体层能否看出 | `q100` 在某个格子里，受体层指标相对 `none` 和相对打乱都 > 0（3 个种子） |

| R1 | R4 | R5 | 结论 |
|---|---|---|---|
| 否 | — | — | 在这组参数下，这种 cost 用不上 S，改公式没有用；需要改 S 进入模型的方式，或改 OT 参数 |
| 是 | 未越过 | — | V12 目前不够准，问题不在公式 |
| 是 | 越过 | 是 | 可以进入真实切片的置换检验（另写计划），使用 R3 选出的（归一化, β） |
| 是 | 越过 | 否 | 真实靶基因评价原则上看不出 S 的作用，真实切片实验需要换评价方法 |
| 是 | 不确定 | — | 需要更精确的 V12 质量估计 |

如果主设置的 R1 为“否”而敏感性设置为“是”，结论是：**ρ 是瓶颈**。边际约束强时，S 项只能重新分配运输量，不能把错误的配对关掉。

## 6. 本地预试（Claude Science 沙箱中的模拟输入，不是结果）

- 320 个 spot、145 对的模拟网格，主设置，`q100`、`raw`：β 从 0 到 1，合并 Spearman 从 0.416 升到 0.563。β ≥ 0.5 时 50,000 次迭代没有收敛，所以 nitermax 定为 200,000。β = 1 时诱饵仍拿到约 19% 的运输量。
- 84 个 spot、30 对的极小网格：敏感性设置约 100 次迭代收敛，`q100` 在 β = 1 时诱饵运输量为 0；主设置约需 44,000 次迭代。

以上只用来确定迭代上限、检查流程，不作为结论。

## 7. 执行

在 `isoform_commot/` 目录下：

0. 准备区域输入（只做一次）：
   `python stage2.py prepare --sample posterior --data data/visium_mouse_brain_sagittal_posterior --max-spots 400 --output stage2_input_post_roi400`
1. 提交本文件和代码，记下提交号。
2. 冒烟测试：
   `bash reports/stage2_inject_01/run_inject.sh --smoke --input stage2_input_post_roi400 --check-lr stage2_input_post --runs analysis_inject_smoke --jobs 8`
   内容：2 个种子，`q100` 和它的打乱版本，`raw`，β = 0.1 和 1，迭代上限 20,000，结果不解读。
   通过条件：退出码为 0；`analysis_inject_smoke/summary/main` 和 `.../sens` 下都有 `decisions.json`。从 `logs/` 里记下每次求解的耗时，用来估计正式运行的时间。
3. 正式运行：
   `bash reports/stage2_inject_01/run_inject.sh --input stage2_input_post_roi400 --check-lr stage2_input_post --runs analysis_inject_01 --results reports/stage2_inject_01/summary --jobs 24 --v12-report reports/stage2_ppi_01/scores/score_report.json`
   如果第一步打分（`stage2_ppi_01`）还没完成，可以先不加 `--v12-report`。之后再用 `python stage2_inject.py summarise --runs analysis_inject_01/main --output reports/stage2_inject_01/summary_r4/main --v12-report …` 计算 R4（敏感性设置同理）；`summarise` 只读已有的运行结果，不重跑 OT。
4. 提交 `reports/stage2_inject_01/summary/`。`analysis_inject_*` 由现有 `.gitignore` 的 `analysis*/` 规则忽略。

环境：所有步骤都在 `commot-research`（CPU，不需要 GPU）。
计算量：636 次求解。在本地 320 个 spot 的模拟输入上，主设置每次求解约 1.5–5 分钟；ROI400 的实际耗时由冒烟测试测出。并行数用 `--jobs` 设置。

## 8. 已知局限

1. 真值是人为规定的规则，结论只在这个规则下成立。
2. 合成 S 的噪声模型（潜变量加 sigmoid）是简化。V12 的误差结构可能不同，例如有长度偏好（`reports/isoform_lr_bench_01/results.md`）。
3. 只用了后部切片中心的 400 个 spot。
4. 诱饵比例 1/3、K 跨两个数量级、λ = 100 µm 都是固定的选择，没有扫描。
