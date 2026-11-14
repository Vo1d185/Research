# 第二阶段计划

2026-11-14 确定，**在运行任何第二阶段 OT 之前提交**。第 1–5 节是已确定的决定；
为了能执行而补充的细节都列在第 7 节"待确认"，在确认之前按那里写的默认值实现，并在结果中标明。

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
- **diffusion：** 所有 d ≤ 500 µm 的 spot 对（包括 d = 0）。
- **r0 用 µm 表示**；diffusion 的核为 F = r0 / (r0 + d)，d = 0 时 F = 1。

**注意：** COMMOT 教程的 `dis_thr = 500` 直接作用在 `obsm['spatial']`（全分辨率像素坐标）上，
换算后约 **363 µm**，而它的文档写作"500 µm"。本研究按决定使用 500 µm（约 688 px）。

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

## 7. 待确认（为能执行而补充的默认值）

| # | 问题 | 当前默认 | 说明 |
|---|---|---|---|
| 1 | LR 库来源 | CellChatDB v2（盘点所用），Secreted + Contact | COMMOT 教程用 COMMOT 内置的较早 CellChat 表，且只用 Secreted |
| 2 | 参与 OT 的 LR 对 | 全部通过第 2 节筛选的 Secreted + Contact 对 | 与 COMMOT 一样对筛选后的全部配对求解；靶基因评价只用其中可评价的受体 |
| 3 | 我们模型的 ε、ρ | 沿用第一阶段工作设置 ε = 0.2、ρ = 1 | |
| 4 | r0 取值（µm） | 50、100、200、400（第一阶段 r0/h = 0.5、1、2、4，Visium h = 100 µm）；工作设置 200 µm | |
| 5 | 倒数距离对照 | 不包括 | d = 0 时 d_ref/d 无定义 |
| 6 | COMMOT 基线的截断 | 500 µm（与本研究相同），另加 500 px（约 363 µm，教程原样）作对照 | ε = 0.1、ρ = 10，默认四方案混合和纯 collective（weights 1 0 0 0） |
| 7 | 复合物受体"有 TF 关联" | 任一亚基有关联即可 | 盘点 02：Contact 下两种规则结果相同，Secreted 略有差别 |
| 8 | 重复基因符号 | 按 COMMOT 用 `var_names_make_unique`（第一个副本保留原名） | 盘点 02 曾按同名相加；涉及 Ccl19、Ccl21a/b、Ccl27a、Il11ra2、Pcdha11 |
| 9 | 小区域对比的范围 | 后部切片中心附近 128 个 spot，Secreted 和 Contact 各取若干对 | 稠密版本的规模限制（≤ 4,000,000 个 cost 元素） |
