# isoform 部分方案草稿：按结合结构域判断兼容性

**状态：草稿，只是文档，没有运行任何东西。等与导师讨论后再执行。** 2026-11-15。
依据 [`reports/isoform_candidates_02/`](isoform_candidates_02/README.md)（后部面板 358 个基因，97 个候选；
02 在 [01](isoform_candidates_01/README.md) 的基础上把同一基因的多个已审阅条目当作 isoform，新增 Nrxn2）。
标"待讨论"的是需要决定的地方；标"待文献核对"的生物学说法在执行前要逐条核对原文。

## 0. 要回答的问题

Visium 3′ 计数只给基因层面的表达，分不出 isoform。若某个 isoform 丢失了与配对蛋白结合的结构域，它的表达不应进入这个 LR 对的预算。
本方案分两步：(1) 对每个"isoform × LR 对"判断是否保留结合能力（兼容、不兼容、无法判断）；
(2) 由于不知道各 isoform 的表达比例 π，在保持基因总量不变的前提下扫描 π，看 OT 结果对 isoform 组成有多敏感。
目标是敏感性分析，不是估计真实的 π。

## 1. isoform 的范围

- 每个基因 g 的 isoform 集合 I_g：UniProt 已审阅条目 ALTERNATIVE PRODUCTS 中的 isoform（canonical 加 alternative products）。
- **同一基因有多个已审阅条目时**（Nrxn1、Nrxn2、Nrxn3 的 α/β），把各条目及其 isoform 合并为同一基因的 isoform 集合。
  条目按主基因名和 MGI 交叉引用归入基因（盘点 02 的规则；本面板中没有只能通过 MGI 找到的条目）。
- **Nrxn2 在盘点 01 中被漏掉的原因：** 盘点 01 的候选规则只在单个条目内找带标记的 `Alternative sequence` 特征。
  Nrxn2 的 α（E9Q7X7）和 β（E9PUN2）是两个已审阅条目，各自的 ALTERNATIVE PRODUCTS 列出 1a、1b 两个 isoform，
  但另一个 isoform 的序列在另一个条目中，所以两个条目都没有 `Alternative sequence` 特征；01 列出了"有多个条目的基因"，
  但没有比较条目之间的序列。Nrxn1、Nrxn3 的 α 条目本身有 `Alternative sequence`，所以没有被漏掉。
  盘点 02 增加了条目间比较（与最长条目比较，相同的连续片段 ≥ 10 aa 算共有），Nrxn2 因 α 的 1–1135 在 β 中被替换而成为候选。
- 盘点中找不到条目的 11 个基因（Flt3l、Kitl、Pcdhb* 等）暂不处理（待讨论：是否用同义名查找）。

## 2. 兼容性判断（每个 isoform k × 每个含该基因的 LR 对）

### 2.1 证据

1. **isoform 序列与 canonical 的对应：** 由 UniProt isoform 序列与 canonical 全局比对（或直接用 `Alternative sequence` 特征及其"in isoform …"描述），
   得到每个 isoform 保留、删除、替换的 canonical 区段。
2. **拓扑：** `Signal`、`Transmembrane`、`Topological domain`（胞外 / 胞内）。
3. **结合决定区（BD）：** 依次取
   (a) UniProt `Region` 中描述提到该配对蛋白或其家族（如 "Interaction with NLGN1"）的区段，以及配体是该配对蛋白的 `Binding site`；
       配体是金属离子或小分子的 `Binding site`（如 neurexin LNS 结构域中的 Ca²⁺ 位点）不算结合决定区
       （盘点 01/02 的"结合"标记包含这些位点，只用于筛选候选）；
   (b) 没有 (a) 时，用 Pfam 结构域（经 InterPro 映射到 canonical 坐标）和 UniProt `Domain`，**但只有在人工整理的表中写明"该类相互作用由此结构域介导"并附文献时才用**；
   (c) 都没有：BD 未知。
4. 所有来源记录版本（UniProt release、InterPro / Pfam release）和查询时间，原始响应保存哈希。

### 2.2 分类

| 类别 | 条件 |
|---|---|
| **不兼容** | (i) 拓扑不满足：canonical 有信号肽而 isoform 没有（分泌或膜蛋白无法到达胞外），或 canonical 有跨膜区而 isoform 没有、且该基因在配对中作为需要在膜上发出信号的受体；(ii) BD 完全缺失，或保留部分 < 50%（待讨论阈值） |
| **兼容** | 拓扑保留，BD 完整保留且序列不变 |
| **无法判断** | BD 未知；或 BD 基本保留但其中有插入、替换（如短的剪接插入），且没有针对这个配对的文献证据；或 isoform 序列标注为片段 |
| **有文献依据的修饰** | BD 内的可变区对特定配对的影响有直接实验证据：按文献判为兼容或不兼容，逐条注明出处（只用于人工整理表中的条目） |

待讨论：
- 可溶形式（缺跨膜区、保留胞外结合区，如可溶受体）：结合但不发出信号，是记为不兼容、单独一类"诱饵"，还是在 OT 中作为不产生下游信号的竞争者。
- 复合物受体：判断落在哪个亚基上（只对含 isoform 的亚基判断；复合物表达仍取亚基最小值）。
- 50% 阈值和"短插入"的长度界限。

### 2.3 输出

每个"LR 对 × isoform"一行：类别、依据（BD 来源、保留比例、拓扑、文献）、各数据库版本。

## 3. 表达在 isoform 之间分配（π 扫描）

- 基因层面表达 x_g(s)（与 OT 预算同一尺度：normalize_total + log1p）。isoform k 的表达 x_{g,k}(s) = π_{g,k} · x_g(s)，
  π_{g,k} ≥ 0，Σ_k π_{g,k} = 1，**基因总量不变**。第一版 π 在空间上处处相同（没有空间分辨的 isoform 数据）。
- **按类别合并：** 对每个 LR 对，把基因 g 的 isoform 合并为兼容（C）、不兼容（I）、无法判断（U）三组，只扫描组的比例 π_C、π_I、π_U。
  建议网格：π_C ∈ {0, 0.25, 0.5, 0.75, 1}，其余按 I : U = 1 : 1 分配（待讨论）；U 组给出两个界：全部算兼容（上界）、全部算不兼容（下界）。
- **进入 OT 的方式：** 含 isoform 的基因拆成 isoform 组实体，各自带预算 π · x_g；LR 连接只连到该配对的兼容组（上界时也连 U 组）。
  这样不兼容的 isoform 仍占基因总量，但不参与这个配对，组之间在 collective OT 中照常竞争。其他基因不变；cost、ε、ρ、截断、r0 用第二阶段主分析的设置。
- **同一基因在多个配对中：** 一个 isoform 可能对配对 A 兼容、对配对 B 不兼容（如 neurexin 的 SS4 对不同配体），所以按 isoform 拆实体，而不是按配对拆。
- **报告：** 随 π 变化的受体层面接收信号、运输比例、靶基因评价指标；找出对 isoform 组成最敏感的 LR 对。结果解释为"如果 isoform 比例是 π，则……"。

## 4. 示例基因（先做 NRXN1、NRXN2、NRXN3）

| 基因 | UniProt | 盘点 01 中的可变区 | 后部面板中的配对 |
|---|---|---|---|
| Nrxn1 | α Q9CS84（1514 aa）、β P0DI97（468 aa） | α 上 6 个：379–393、387–393（Laminin G-like 2，即 SS2 位置，待核对）；1247–1276（Laminin G-like 6，30 aa，即 SS4 位置，待核对）；isoform 6 替换 1–1372（类似 β 的 N 端）。条目间：α 1–1129 被 β 1–83 替换，C 端 385 aa 相同 | Adgrl1、Clstn1–3、Dag1、Lrrtm1–4、Nlgn1–3（12 对） |
| Nrxn2 | α E9Q7X7（1710 aa）、β E9PUN2（660 aa） | 两个条目都没有 `Alternative sequence`（盘点 01 因此漏掉）。条目间：α 1–1135 被 β 1–85 替换，C 端 575 aa 相同 | 与 Nrxn1 相同的 12 个受体 |
| Nrxn3 | α Q6P9K9（1571 aa）、β Q8C985（567 aa） | α 上 3 个（1–373 替换，含 EGF-like 1、Laminin G-like 1–2；1047；1364–1470）。条目间：α 1–1088 被 β 1–81 替换；β 在 α 1364/1365 之间多 3 aa | Adgrl1、Lrrtm1–4、Nlgn1–3 |

**待文献核对的要点**（执行前逐条找原始文献）：
- α-neurexin 有 6 个 LNS（Laminin G-like）结构域和 3 个 EGF 样结构域，β-neurexin 只有与 α 的 LNS6 相同的一个 LNS 结构域和 β 特有的 N 端。
- 对 neuroligin、LRRTM、latrophilin（Adgrl1）的结合由 LNS6 介导，β 与 α 都能结合；SS4 插入调节结合：如 LRRTM2 只结合 SS4−，
  neuroligin-1 的结合受 SS4 和 neuroligin 自身剪接位点 B 共同影响。
- dystroglycan（Dag1）结合 α-neurexin 的 LNS 结构域，且受 SS2、SS4 插入影响；calsyntenin（Clstn）结合的部位待查。
- 若这些都成立，neurexin 的兼容性主要取决于 α/β（对需要 LNS1–5 的配对）和 SS4 状态（对 LNS6 配对），正好检验"BD 内短插入"的规则。

**第二批（结构域依据清楚、配对在面板中）：**
- Fgfr1、Fgfr2、Fgfr3：Ig 样结构域 D3 的 IIIb / IIIc 互斥外显子改变 FGF 配体特异性（待文献核对）。面板中的配对有 Fgf1、Fgf9、Fgf10（Fgfr1 还有 Ncam1）。
  这是"可变区落在结合结构域内、且对不同配体结论不同"的典型例子，可以作为规则的阳性对照。
- Ptprd、Ptprs：Ig 样结构域中的微外显子（meA、meB）调节与 Slitrk、Lrrc4b 等的结合（待文献核对）。面板中有 Slitrk1–5、Lrfn4/5、Lrrc4b 等配对。
- Cd44（13 个 isoform）：可变外显子在胞外近膜区，与 Spp1、Lgals9 结合的部位待查。可以作为"可变区在胞外但不在已知 BD 内"的例子。

## 5. 需要讨论后才能定的事

1. 第 2.2 节的阈值、可溶形式和复合物的处理。
2. π 网格与 U 组的处理（两个界还是单独扫描）。
3. 第二阶段 r0 选定前后：isoform 分析用哪个 r0（建议用后部选出的 r0\*，不再重新选择）。
4. 是否需要 isoform 层面的外部证据（例如小鼠脑长读长或剪接位点表达数据）来约束 π 的范围；没有的话只做扫描。
5. 先在哪张切片上做（建议只在后部做示例，前部用于检查结论是否一致）。

## 6. 执行顺序（讨论确认后）

1. 写成计划条目并提交（规则、阈值、π 网格、示例基因），在看任何 isoform 相关的 OT 结果之前。
2. 实现兼容性判断，人工核对 NRXN1–3 的结果和文献表，提交。
3. 实现 isoform 组实体与 π 扫描，先在小区域与不拆分的结果核对（π_C = 1 且 I、U 为空时应与第二阶段结果完全相同），提交。
4. 在后部切片上跑示例基因，报告。
