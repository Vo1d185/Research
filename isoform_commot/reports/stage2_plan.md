# 第二阶段计划

2026-11-14 确定，**在运行任何第二阶段 OT 之前提交**。第 1–5 节是已确定的决定；
为了能执行而补充的细节都列在第 7 节"待确认"，在确认之前按那里写的默认值实现，并在结果中标明。

**第二次决定（2026-11-14，在任何新的 OT 之前提交）：** 第 7 节的结论见该节；方法比较、ε/ρ、截断和 r0 选择规则见第 8 节。
第 3 节的 diffusion 截断随之修改（原文保留在注释中）。

## 1. 数据

| 用途 | 数据 | 来源与 SHA256 |
|---|---|---|
| 主数据集 | 10x Visium `V1_Mouse_Brain_Sagittal_Posterior`，Space Ranger 1.1.0，3355 个组织内 spot | `…/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/`；matrix `df962209143860488bd679ff7a1b3bb14985208d0929fb814934ba7571c488b6`；spatial `b878d7eaca487ffc17224a870f86599abf9d45a99033c1ef65661ff119cef0d5` |
| 独立重复 | 10x Visium `V1_Mouse_Brain_Sagittal_Anterior`，Space Ranger 1.1.0 | `…/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/`；matrix `56078d8d6fe6c13de248fdb1c518b691cdef78fb00021b659786b4a47c6656d5`；spatial `5f41a803e2bd69fa4dfca6abc8fa2d4e0d76aeb6c72d7038a5fdcf9cc50a36f8` |

SHA256 为本地首次下载时计算（2026-11-14），不是 10x 发布的校验值。
所有设置在后部切片上确定；前部切片只用于重复，**不在前部上调整任何参数**。

## 2. 预处理与 LR 筛选（沿用 COMMOT 的 Visium 教程）

依据：COMMOT `d117445` 的 `docs/notebooks/visium-mouse_brain.ipynb`。

- 表达：`var_names_make_unique` → `normalize_total`（默认目标值，即各 spot 总计数的中位数）→ `log1p`。
- 筛选：`filter_lr_database(min_cell_pct = 0.05)` 的规则，**每个亚基各自**在 ≥5% 的 spot 中检测到（count > 0），整张切片计算。
- 复合物表达：各亚基表达（归一化、log1p 之后）的**最小值**（COMMOT `heteromeric_rule = 'min'`）。
- LR 库：CellChatDB v2 小鼠（jinworks/CellChat `75253cd`，`dd59b196…`），Secreted Signaling 作 diffusion，
  Cell-Cell Contact 作 contact（见第 7 节第 1 条）。

## 3. 空间设置

- **坐标换算成 µm：** µm/像素 = 65 / `spot_diameter_fullres`（10x 定义该值为 65 µm 理论 spot 直径对应的全分辨率像素数）。
  后部切片 89.513 px → 0.7262 µm/px。用相邻 spot 中心距核对：约 138 px → 100.2 µm，与 Visium 的 100 µm 一致。
  准备输入时自动核对，偏差超过 2% 拒绝执行。
- **允许同一 spot 内通信**（d = 0）。
- **contact：** 本 spot 加上相邻的 6 个 spot。相邻关系由 `array_row`/`array_col` 的六边形网格确定
  （同一行列号 ±2；行号 ±1 且列号 ±1），只连组织内的 spot，所以边缘 spot 的邻居少于 6 个。
- **diffusion：** 截断内所有 spot 对（包括 d = 0）。**主分析截断 363.08 µm（= COMMOT 教程的 500 px），500 µm 为敏感性分析**
  （第二次决定；第一次决定为 500 µm）。换算见第 8.3 节。
- **r0 用 µm 表示**；diffusion 的核为 F = r0 / (r0 + d)，d = 0 时 F = 1。

**注意：** COMMOT 教程的 `dis_thr = 500` 直接作用在 `obsm['spatial']`（全分辨率像素坐标）上，
换算后约 **363 µm**，而它的文档写作"500 µm"。第一次决定使用 500 µm（约 688 px）；
第二次决定改为主分析用 363.08 µm（与教程相同的物理距离），500 µm 作敏感性分析。

## 4. 靶基因评价

- **只用 scSeqComm**（`95c9bb2`，KEGG + Reactome，小鼠符号）。
- **放弃 CCCvelo 库的理由**（盘点 02）：
  1. RecTFDB 是稠密表：762 个受体 × 1145 个 TF，每一行都有正分数，"受体有 TF 关联"只等于"受体在这 762 个之内"；
  2. TFTGDB 中前 10 个 TF 的靶基因并集，各组合的中位数为 777–3036 个，"≥5 个靶基因"从不起筛选作用，
     靶基因得分会变成几百到几千个基因的平均表达，失去受体特异性；
  3. 使用人类风格的大写符号，只能靠大写符号匹配小鼠基因，不是同源映射。
- **靶基因检测阈值 5%**（整张切片）。
- **复合物受体的 TF：** 在各亚基的关联中，每个 TF 取最高权重，再取前 10 个。
- 其余规则沿用 `evaluate_targets.py`：前 10 个 TF 为主分析，前 5、前 20 为敏感性分析；靶基因 ≥5 个；
  去掉所有参与 OT 的 LR 基因（包括亚基）；不报告 p 值。

## 5. 实现：稀疏版本

1. 每个 LR 对单独构建稀疏 cost（只存可连接的 spot 对，d = 0 和 contact 的零 cost 显式保存），再拼成整体问题。
2. 求解器：固定版本 COMMOT `unot_sinkhorn_l1_sparse`（`_unot.py`，commit `d117445`，SHA256 `4b50aca8…`，已在 `fetch_solver.py` 中校验）。
3. 质量守恒检查：该函数只返回 T，不返回对偶变量。利用"最后一步更新 g 之后列约束精确成立"，
   取 ν = b − Tᵀ1，由 T = μ·ν·exp((2ρ − C)/ε) 反推 μ，检查 |T1 + μ − a| < 1e-7（与稠密版本相同的检查），
   并检查每个存储元素都满足这一形式；禁止的连接不进入矩阵。
4. 在小例子上与稠密版本逐项对比（单元测试），通过标准：max |T_sparse − T_dense| ≤ 1e-6 · max T_dense。
5. case 结果以稀疏格式保存；`evaluate_targets.py` 和 `commot_baseline.py` 都支持稀疏格式。
6. 顺序：Visium 后部的一个小区域 → 与稠密版本逐项对比（同一标准）→ 整张切片。

## 6. 本阶段不做

- PPI 分数（已按预先登记停止，S = 1）。
- isoform 层面的表达分配。
- 在前部切片上选择或调整参数。

## 7. 待确认事项的结论（第二次决定）

| # | 问题 | 结论 |
|---|---|---|
| 1 | LR 库来源 | **已确认：** CellChatDB v2，Secreted + Contact |
| 2 | 参与 OT 的 LR 对 | **已确认：** 全部通过第 2 节筛选的 Secreted + Contact 对 |
| 3 | 我们模型的 ε、ρ | **由第 8 节取代：** 主分析 ε = 0.1、ρ = 10；ε = 0.2、ρ = 1 为敏感性分析 |
| 4 | r0 取值（µm） | **已确认：** 50、100、200、400。原"工作设置 200 µm"取消，改由第 8.5 节的规则选择 |
| 5 | 倒数距离对照 | 按默认执行（不包括；d = 0 时 d_ref/d 无定义），未单独确认 |
| 6 | COMMOT 基线 | **由第 8 节取代** |
| 7 | 复合物受体"有 TF 关联" | **已确认：** 任一亚基有关联即可 |
| 8 | 重复基因符号 | 按默认执行（`var_names_make_unique`），未单独确认 |
| 9 | 小区域对比的范围 | 已完成（128 个 spot，第 5 节的稠密对比全部通过） |

原第 7 节的默认值（第一次提交 `8348f7c`）：LR 库 CellChatDB v2 Secreted + Contact；全部筛选后配对；ε = 0.2、ρ = 1；
r0 = 50/100/200/400 µm，工作设置 200 µm；不含倒数距离对照；COMMOT 基线截断 500 µm 加 500 px；复合物"任一亚基"；
`var_names_make_unique`；小区域 128 个 spot。

## 8. 方法比较、参数与 r0 选择（第二次决定）

### 8.1 条件（后部切片）

主比较是我们的 cost 与 COMMOT 的 cost：两者 ε、ρ、预算归一化、运输方案（纯 collective）和求解设置都相同，**唯一的差别是 cost**。

| 编号 | 方法 | ε | ρ | 截断 | 运输方案 | 求解设置 | 用途 |
|---|---|---|---|---|---|---|---|
| A | 本研究 cost，r0 = 50、100、200、400 µm | 0.1 | 10 | diffusion 363.08 µm | collective | nitermax 50000，stopthr 1e-11 | 主比较 |
| B | COMMOT cost | 0.1 | 10 | 363.08 µm | 纯 collective（weights 1 0 0 0） | 同 A | 主比较 |
| C | COMMOT cost | 0.1 | 10 | 363.08 µm | 默认四方案混合（各 0.25） | COMMOT 默认（nitermax 1e4，stopthr 1e-8） | "按原文"基线 |
| A500、B500、C500 | 同 A、B、C | 0.1 | 10 | 500 µm | 同 A、B、C | 同 A、B、C | 敏感性：截断 |
| Aε、Bε | 同 A、B | 0.2 | 1 | 363.08 µm | 同 A、B | 同 A | 敏感性：ε、ρ |

敏感性分析每次只改变主设置中的一个因素。全部条件都做第 4 节的靶基因评价。

### 8.2 "唯一的差别是 cost"的具体做法

- **预算归一化相同：** 两者都在求解前把全部表达预算除以 max(Σa, Σb)（COMMOT `cot_sparse`：`max_amount = max(S.sum(), D.sum())`；
  本研究 `solve_sparse` 相同），只在预算非零的元素上求解，再乘回原单位。
- **求解器相同：** A 和 B 都用固定版本 `unot_sinkhorn_l1_sparse`，ε、ρ、nitermax、stopthr 相同。
  B 直接调用 COMMOT 的 `cot_sparse`，即 `cot_combine_sparse` 中 collective 方案所用的函数；weights 为 1 0 0 0 时
  另外三个方案的权重为 0，结果相同（单元测试核对），这样不必求解三个被丢弃的问题。
- **质量守恒检查相同：** B 用 COMMOT 内部归一化后的 cost 做与 A 相同的检查（|T1 + μ − a| < 1e-7 和元素的乘积形式），
  确认它是同一类问题的收敛解。C 按设计不满足单个实体的预算，只记录超出量。
- **cost（唯一不同之处）：**
  - 本研究：diffusion 为 log1p(d / r0)（d ≤ 截断）；contact 为 0（本 spot 加相邻 6 个 spot）；其余禁止。不做缩放。
  - COMMOT：d / cost_scale（d ≤ 截断），cost_scale 是截断内最大的 spot 间距离（COMMOT 内部的 cost 归一化，原样保留）；
    与 COMMOT 一样不区分 contact 和 diffusion，contact 配对也连接到截断距离。contact 的连接范围属于我们的
    cost（contact 的 F），因此这一差别算作 cost 的差别。

### 8.3 截断换算

- 363.08 µm = 500 px × 65 µm / 89.513 px（后部切片的 `spot_diameter_fullres`，0.72615 µm/px；第 3 节的邻距核对为 100.2 µm），
  即 COMMOT 教程把 `dis_thr = 500` 直接用在全分辨率像素坐标上所对应的物理距离。
- 后部切片上，363 µm 与 363.08 µm 给出完全相同的 spot 对：截断内最远的一对是 362.67 µm（六边形网格 √13 × 间距），
  截断外最近的一对是 399.38 µm（4 × 间距）。
- 前部切片使用同一个 µm 值。截断只比 √13 环高约 0.1%，前部切片运行时报告边数和这两个边界距离。

### 8.4 运输比例与问题规模

预算按总量归一化时，**ρ = 1 下的运输比例会随问题规模变化**：ε = 0.2、ρ = 1、截断 500 µm 时，后部切片 128 个 spot、12 对的小区域上
约 78% 的配体预算被运输，整张切片（540 对）上只有约 20–26%。原因是归一化后每个元素的预算随问题变大而变小，
而熵项和 ρ 不随之缩放。因此主分析两种方法都用 ε = 0.1、ρ = 10（COMMOT 默认），ε = 0.2、ρ = 1 只作敏感性分析。
每种条件都报告运输比例 ΣT / Σa 和 ΣT / Σb。

### 8.5 r0 选择规则

- 只在**后部切片**、条件 A 下选择。对每个 r0，取**全部可评价受体的平均偏 Spearman**：受体层面的接收信号（该受体所有配体之和），
  主靶基因集（前 10 个 TF，靶基因平均表达），控制受体表达（`evaluate_targets.py` 的 `partial_spearman`，level = receptor）。
  可评价受体 = 靶基因集状态为 evaluated，且四个 r0 下偏 Spearman 都有定义；报告受体数。
- 选平均值最大的 r0；并列时取较小的 r0。
- **前部切片只用来验证**：使用后部选出的 r0，不参与选择。
- r0 的选择和后部的方法比较用同一份数据，所以后部上 A（所选 r0）对 B 的比较偏向 A；前部切片是独立检验。
