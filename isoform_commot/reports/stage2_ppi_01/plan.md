# Stage 2 加入 PLM-interact 分数：第一步（打分与描述）

**状态：预注册计划。** 2026-11-18 起草。本文件提交后才运行正式打分；`run_scoring.sh` 会检查本文件已提交且工作区无改动。
**本步不运行任何 OT、不查看任何靶基因评价。** 第二步（带 S 的 OT）的设计见第 6 节，现在只是草案：第一步结果出来后定稿，另行提交，提交前不运行任何带 S 的 OT。

代码：`stage2_ppi.py`（genes / jobs / import / designs），`test_stage2_ppi.py`，`reports/stage2_ppi_01/run_scoring.sh`。
沿用已有代码：`proteins.py map`（基因 → 小鼠蛋白），`plm_score.py score`（打分），均不改动。

## 0. 背景

- Stage 2 主分析按预注册停用了 PPI 分数（S = 1，见 `HANDOFF.md`）。当时的依据是：humanV11 对面板中的已知配对给出接近 0 的分数，已知配对与未注释组合的 AUROC 只有 0.68；而且按 −log S 计算，低分会主导运输结果。这条记录保留，不覆盖。
- 本计划按导师的建议重新加入 S，做法是**基因层、按配体归一化**：S′(L,R) = S(L,R) / Σ_R′ S(L,R′)。模型换成 humanV12。
- 只用基因层：`reports/isoform_lr_bench_01/results.md` 显示，humanV12 区分不了同一基因的 isoform。
- 基因层参考：`reports/plm_lr_benchmark/report.md` 中，humanV12 在人类 H-M 上随机负样本 AUROC 0.812、保持度数负样本 0.783、正样本 > 0.5 的占 60.8%。它在小鼠面板上的表现未知，第 4 节只做描述。

## 1. 输入

`stage2_input_post`、`stage2_input_ant`（`stage2.py prepare` 的输出，使用其中的 `lr.csv`、`entities.csv`）。打分对象为两张切片 LR 对的并集。
已知规模（`reports/stage2_pairs_01`，不需要分数）：

| | 后部 | 前部 |
|---|---|---|
| LR 对 | 540 | 576 |
| 配体 | 201 | 222 |
| 只有 1 个受体的配体（S′ 恒为 1） | 95（对应 95 对） | 103（对应 103 对） |
| 每个配体的受体数：中位 / 最大 | 2 / 12 | 2 / 12 |
| 含复合物受体的对 | 161 | 168 |
| 同源配对（配体 = 受体） | 32 | 40 |

每个配体只属于一种类型（分泌型或接触型），所以按配体归一化不会把两种类型混在一起。

## 2. 打分规则（先定，不随结果修改）

1. **序列：** 小鼠 Swiss-Prot canonical 全长。`proteins.py map --genes … --species mouse --synonyms --allow-unresolved`。
   没有唯一已审阅条目的基因列为“未解决”，不猜，本步不加人工 override。例如有 α、β 两个条目的 Nrxn1、Nrxn2 很可能属于这种情况。
2. **模型：** `danliu1226/PLM-interact-650M-humanV12` @ `2d0e7fb3deef94efa3f217daccc79092a6a42ceb`，checkpoint sha256 `92565937…`。
   `stage2_ppi.py import` 会核对 sha256，不符时拒绝运行（冒烟测试除外）。
3. **打分：** `plm_score.py score`，逐对、不截断、`--allow-long`。两个输入顺序都打分，**取平均**。
4. **复合物：** 拆成亚基，每个“配体亚基 × 受体亚基”组合各打一个分，**配对分数取最大值**；同一个蛋白和自己的组合不计入。规则与 `plm_lr_benchmark_plan.md` 的 H-X 相同。
5. **同源配对：** 配体和受体是同一个蛋白时（拆亚基后只剩自身组合），按“蛋白与自身”打分，标记 `homophilic`，单独报告。
6. **缺失：** 任一亚基未解决的配对没有分数，**不插补**。

## 3. S′ 的定义（第二步使用；本步只计算和描述）

对每张切片分别用该切片的 `lr.csv` 计算（`stage2_ppi.py designs`）：

| 列 | 定义 |
|---|---|
| `s_raw` | S |
| `s_norm_ligand`（第二步主分析） | S(L,R) / Σ_{R′ ∈ R(L)} S(L,R′)，R(L) 为该切片中与 L 配对的受体；组内任一对缺分数，则整组缺失 |
| `s_uniform_ligand`（对照） | 1 / k(L)，k(L) = R(L) 的大小 |
| `s_norm_receptor`（敏感性） | S(L,R) / Σ_{L′} S(L′,R) |

**需要预先说明的偏差：** 若一个配体的受体分数都差不多，S′ ≈ 1/k，−log S′ ≈ log k。这会让受体多的配体整体 cost 更高，而这个效应来自受体数量，不来自 S。所以第二步必须用 `s_uniform_ligand` 作为对照（第 6 节）。

## 4. 本步输出与描述性统计（不作判定）

输出目录 `reports/stage2_ppi_01/scores/`（不用 `results*` 命名，因为 `.gitignore` 忽略 `results*/`）：
- 分数：`ppi_scores.csv`（每个配对一行）、`subunit_pair_scores.csv`；
- S′ 表：`designs_stage2_input_post.csv`、`designs_stage2_input_ant.csv`（各附 `.json` 统计）；
- 报告：`score_report.json`；
- 映射与作业记录：`protein_map.csv`、`proteins.json`、`unresolved_candidates.csv`（如有）、`entity_pairs.csv`、`jobs.json`、`scores.meta.json`、`run_info.json`。

报告以下内容，全部为描述，不作判定：

1. **覆盖率：** 按类型、复合物、同源配对分别报告有分数的配对数；列出未解决的基因。
2. **分数分布：** 分位数，以及 > 0.99、< 0.01、> 0.5 的比例（基准 01 中分数是两极化的）；两个输入顺序之差。
3. **每张切片：**
   - 组内任一对缺分数的配体数；
   - 多受体配体中，组内原始 S 的极差 < 0.01 的配体数，即 S′ 基本等于 1/k 的配体；
   - −log S′ 对 log k 的 R²，即 S′ 的变化有多少只是 1/k 造成的。

### 4.1 V12 在本面板上的区分能力（描述性，供注入实验使用）

打分时，额外从面板中随机抽取**未注释的组合**：
- 配体和受体都来自面板，且属于同一种类型；
- 排除面板中已有的配对、CellChatDB 中任何类别下有注释的配对（`data/cellchatdb`）、配体与受体是同一实体的组合，以及两边共享亚基的组合；
- 数量为有分数的面板配对数的 3 倍（`--cross-ratio 3`），随机种子 1。

它们**不是阴性**，只是没有注释。报告面板配对与这些组合的 AUROC（总体和按类型），附 bootstrap 95% 区间（2,000 次，种子 0）。
这个 AUROC 用来在 `reports/stage2_inject_01` 的质量梯度上给 V12 定位。
输出 `cross_scores.csv`，以及 `score_report.json` 的 `panel_vs_unannotated` 一项。

## 5. 执行

在 `isoform_commot/` 目录下：

1. 提交本文件和代码，记下提交号。
2. 冒烟测试（CellChatDB 目录默认 `data/cellchatdb`，可用 `--cellchat` 指定）：
   `bash reports/stage2_ppi_01/run_scoring.sh --smoke --input stage2_input_post --input stage2_input_ant --runs plm_jobs_stage2_ppi_smoke --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D`
   内容：单元测试、完整的基因映射、6 个配对加最长的一对，分数不解读。
   通过条件：全部步骤退出码为 0，`plm_jobs_stage2_ppi_smoke/results/` 下有 `ppi_scores.csv`、`score_report.json` 和两个 `designs_*.csv`。
3. 正式运行：
   `bash reports/stage2_ppi_01/run_scoring.sh --input stage2_input_post --input stage2_input_ant --runs plm_jobs_stage2_ppi_01 --results reports/stage2_ppi_01/scores --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D`
   如需复用冒烟测试的映射结果，加 `--proteins plm_jobs_stage2_ppi_smoke/proteins`。
4. 提交 `reports/stage2_ppi_01/scores/`。`plm_jobs_stage2_ppi_*` 由现有 `.gitignore` 的 `plm_jobs*/` 规则忽略。

环境：`commot-research`（单元测试、genes、`proteins.py map`、jobs、import、designs）；`plminteract`（`plm_score.py`，GPU）。
基因映射需要访问 UniProt（服务器可以）。GPU 时间未实测，按作业数估计在一小时以内。

## 6. 第二步草案（不在本次执行；第一步完成后定稿提交）

- **主分析 cost：** C = −log F(d) − log S′（用户选定：在原公式中把 S 换成按配体归一化的 S′）。
- **敏感性分析：** C = −log F(d) − ε·log S′（S′ 作为 KL 的参考测度，Gibbs 核 = F^(1/ε)·S′）；原始 S；按受体归一化。
- **对照：** `no_ppi`（S = 1，即 Stage 2 已有的 A 条件）；`uniform`（S′ = 1/k，分离受体数量效应）；同类型内打乱原始 S 后再归一化（若干种子）。
- **参数：** 与 Stage 2 主分析相同（ε = 0.1、ρ = 10、截断 363.08 µm、r₀\* = 50 µm）。
- **评价：** 与 Stage 2 相同的靶基因评价。只有当 S′ 同时优于 `uniform` 和打乱分布时，才能说“S′ 的逐对差异带来了信息”；具体阈值在第二步计划中定。
- **缺分数的配对：** 第二步计划决定。候选做法：(a) 所有条件（含 `no_ppi`）都只用有分数的配对重跑；(b) 为未解决的基因写有依据的 override 后补打分。
- **注意：** Stage 2 的评价对空间信息不太敏感（打乱坐标后指标仍有约 80%，`reports/stage2_section8/r0_validation.md`）。建议同时在 PDE 模拟中加入已知结合强度的场景（`HANDOFF.md` 下一步第 5 条），在有真值的情况下检验 S′ 的作用。
- **实现：** 在 `stage2.py run` 的稀疏 cost 上，对每个配对加上 −β·log S′；`score_mode` 写入 `summary.csv`。

## 7. 已知局限

1. humanV12 在人类数据上训练，这里用于小鼠序列（跨物种）。`HANDOFF.md` 记录过 humanV11 换成人源序列后结果几乎相同，但 humanV12 没有检查过。
2. 只有 1 个受体的配体约占一半，按配体归一化对它们不起作用。
3. 未解决的基因会让部分配对没有分数，可能集中在有多个已审阅条目的基因上。
4. 本步只给出分数和描述，不能说明 S 对通讯推断有没有帮助；那是第二步的问题。
