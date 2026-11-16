# PLM-interact 在 CellChatDB 配体–受体对上的评测：计划

2026-11-16 写定，**在给任何评测配对打分之前提交**。
第 5 节的判定标准是用户给的建议值，**待确认**；第 8 节列出所有待确认或待决定的事项。
确认之前只做两件事：构建数据集（第 1–2 节）和冒烟测试（第 7 节第 2 步，只用 2026-11-14 已经打过分、看过结果的配对）。

## 0. 问题和范围

- 问题：PLM-interact 照原样使用（不微调），能否给 CellChatDB 里真实的配体–受体对打出正确的互作分数。
- 这是一次评测，**不接入 OT**，也不改变 2026-11-14 的决定（PPI 作为 cost 的工作已按预先登记停止，见 `diagnostics/README.md`）。
  这次的结果不会自动让 PPI 回到 cost 中。
- GPU 打分，不影响正在运行的 PDE 队列（见第 3.6 节）。

## 1. 数据集

### 1.1 来源（全部固定版本，SHA256 在本地首次下载时计算并记录）

| 用途 | 来源 | 固定版本 |
|---|---|---|
| 主数据集 | CellChatDB v2 人类，`jinworks/CellChat` 仓库 `data/CellChatDB.human.rda` | commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`（与小鼠版相同） |
| 补充 | CellChatDB v2 小鼠，同一仓库 `data/CellChatDB.mouse.rda` | 同一 commit；已下载，`dd59b196…` |
| 负样本排除 | CellPhoneDB v5，`ventolab/cellphonedb-data`（服务器可以访问） | tag `v5.0.0` 对应的 commit，下载时记录 |
| 序列 | UniProt REST（reviewed 条目，canonical 序列） | 记录每个条目的 entry/sequence version 和 UniProt release |
| 训练集重叠 | HF `danliu1226/cross_species_benchmarking`（humanV11 的训练集）、`danliu1226/STRING_V12_TrainingSet`（humanV12）、`danliu1226/Bernett_benchmarking`（Leakage-Free） | 下载时记录 HF revision；`cross_species_benchmarking@19dbed2`、`Bernett_benchmarking@5d2ad03` 已查到 |

CellChatDB 是 GPL-3。按已有做法（`reports/cellchat_inventory_01/pairs.csv`），从中导出的配对表注明出处后提交；原始 `.rda` 不提交。

### 1.2 正样本规则

1. 去掉 `annotation = "Non-protein Signaling"`；剩下的类别为 Secreted Signaling / Cell-Cell Contact / ECM-Receptor。
2. **单基因配对：** `ligand` 和 `receptor` 都不是复合物名（不在 `complex` 表中），且 `ligand.symbol`、`receptor.symbol` 各只有一个符号。
3. **含复合物的配对**（次要分析）：复合物按 `complex` 表拆成亚基。每个"配体亚基 × 受体亚基"组合记录一个分数，
   配对分数取这些分数的最大值。两个亚基是同一个蛋白的组合不计入。
4. **去重：** 正样本的单位是无序的蛋白对。由于主分数是两个输入顺序的平均（第 3.2 节），(L, R) 和 (R, L) 是同一个样本。
   同一对出现在多行时合并，记录全部类别；分层时在它的每个类别中各计一次。
5. **同源配对（L = R，如 CDH1–CDH1）：** 保留在数据表中，但**不进入主分析和任何 AUROC 计算**，只单独报告分数分布。
   理由：负样本不能包含 L = R 的组合（未注释的同源二聚体无法当作阴性），保留会让正负样本不可比。（待确认，第 8 节）
6. 基因映射不到唯一 UniProt 条目（第 1.4 节）的配对去掉，并列出。

### 1.3 分析集

| 名称 | 物种 | 配对 | 长度 | 用途 |
|---|---|---|---|---|
| **H-M（主）** | 人类 | 单基因，非同源 | 两个蛋白的 canonical 全长都 ≤ 800 aa | 判定 |
| H-A | 人类 | 单基因，非同源 | 不限 | 长度分层（描述性） |
| H-X | 人类 | 单基因 + 含复合物，非同源 | 不限 | 复合物次要分析（描述性） |
| M-M | 小鼠 | 同 H-M | 同 H-M | 补充（描述性） |
| M-A | 小鼠 | 同 H-A | 不限 | 补充（描述性） |

每个分析集用自己的正样本池生成自己的负样本（第 2 节）。
H-M 的负样本不从 H-A 中截取，因为截取会破坏保持度数的性质。

### 1.4 基因 → UniProt

扩展 `proteins.py map`：加 `--species human|mouse` 和 `--genes`（一列基因符号的 CSV，替代 `--input`）。规则与小鼠完全相同：

- 该物种中**恰好一个** reviewed 条目的**主基因名**等于该符号（不区分大小写），取 canonical 全长序列；
- 零个或多个命中：写入 `<output>_candidates.csv`（候选条目名称和长度），**不猜**；
- 不同基因映射到同一个条目：两者都写入 candidates；
- 新增 `--allow-unresolved`：其余基因照常写出映射，未解决的基因记在 `proteins.json`，相关配对去掉。默认行为（报错退出）不变，小鼠主流程不受影响。
- override 只由用户提供（`lr_symbol,accession,reason`）。NRXN1 在 OT 面板里用 β 型是那个分析的决定，这里**不自动沿用**。

### 1.5 需要统计并汇报的数量

每个分析集：正样本数（按三个类别分开，以及合计）、去重后的蛋白数、两个蛋白都 ≤ 800 aa 的配对比例、
同源配对数、因映射失败去掉的配对数，以及 H-X 中每个配对的亚基组合数分布。

## 2. 负样本

### 2.1 排除集 E（无序，任何方向）

以下任一来源中出现过的蛋白对都不能作负样本：

1. CellChatDB 人类和小鼠的**全部**相互作用（四个类别都算，包括 Non-protein Signaling）；复合物拆成亚基，每个"配体亚基 × 受体亚基"都算。
   小鼠符号转为大写后与人类符号匹配（这是符号匹配，不是同源映射）。
2. CellChatDB `cofactor` 表：同一相互作用的配体 × 其 agonist、antagonist、co_A_receptor、co_I_receptor 基因。
   理由：拮抗物和共受体常与配体直接结合（如 FST–activin、NOG–BMP）。（待确认，第 8 节）
3. CellPhoneDB v5 的全部蛋白–蛋白相互作用（复合物拆成亚基），同时按 UniProt accession 和基因符号匹配。小分子伙伴忽略。

报告 E 的大小，以及每个来源去掉了多少候选组合。

### 2.2 候选池

对分析集 S：配体池 = S 的正样本中担任配体的基因，受体池 = 担任受体的基因。
候选负样本 = 配体池 × 受体池中满足 L ≠ R、{L, R} ∉ E、且与已选负样本不重复（无序）的组合。
H-X 以"实体"（单基因或复合物）为单位：实体对中任一亚基组合在 E 中，或实体对本身有注释，即排除。

### 2.3 随机负样本

数量 = 10 × 正样本数（与 PLM-interact 原文测试集的 1:10 一致），从全部候选组合中均匀、无放回地抽取。
这里的"每个正样本配 10 个"按总体比例实现，负样本不绑定到具体的正样本。随机种子 1。

### 2.4 保持度数的打乱负样本

在配体–受体二分图上做边交换（同一基因作配体和作受体时是两个不同的节点）：

1. 从正样本的边集出发。随机取两条边 (L1, R1)、(L2, R2)，换成 (L1, R2)、(L2, R1)。
   只有两条新边都满足 L ≠ R、不在 E 中、与当前边集不重复（无序），且两条新边互不相同时才接受。
2. 一直交换到：成功交换次数 ≥ 10 × 边数，**并且**没有任何边在 E 中（因此不再含原来的正样本）。尝试次数上限为 1000 × 边数。
3. 到上限仍有边在 E 中：去掉这些边，记录数量和受影响基因的度数偏差。
4. 结果每个配体、每个受体的出现次数与正样本相同，数量 1:1。种子 1、2、3 各生成一套。

### 2.5 提交

每个分析集一个目录 `reports/plm_lr_benchmark/dataset/<集>/`：`positives.csv`、`negatives_random_s1.csv`、
`negatives_degree_s{1,2,3}.csv`、`exclusion_summary.json`；生成代码与数据一起提交。

## 3. 打分

### 3.1 模型

| 版本 | HF 仓库 @ revision | 用途 |
|---|---|---|
| humanV11 | `danliu1226/PLM-interact-650M-humanV11` @ `e86e392dec13dd0c23252c94947b04a7a9821b0e` | **主** |
| humanV12 | `danliu1226/PLM-interact-650M-humanV12` @ `2d0e7fb3deef94efa3f217daccc79092a6a42ceb` | 描述 |
| Leakage-Free | `danliu1226/PLM-interact-650M-Leakage-Free-Dataset` @ `d995319a6e53510bfbb4eff6ab0292fedf2ec4db` | 描述 |

基础模型 `facebook/esm2_t33_650M_UR50D` @ `08e4846e537177426273712802403f7ba8261b6c`。三个 checkpoint 都是已经下载并做过自检的同一批文件。

### 3.2 主分数：`plm_score.py`

逐对打分、不截断（超过 1603 token 时用 `--allow-long`）、核对每对的 token 数。
(L, R) 和 (R, L) 两个顺序都打分，**主分数 = 两个概率的平均**。

### 3.3 官方脚本

`PLMinteract/inference/inference_PPI_singleGPU.py`（commit `ab2ae6ae1aa81accf4c1f7f6f341164baf97809a`），代码不改，
按 README 的单 GPU 命令运行：`torchrun --nproc_per_node=1 … --seed 2 --batch_size_val 16 --max_length 1603 --embedding_size 1280 --model_name esm2_t33_650M_UR50D`。

- 输入是官方要求的 `query,text` CSV。输出 `pred_scores.csv` 没有配对编号，按行序与输入对应，因此每个输入文件旁边保存同序的编号表。
- `len(query) + len(text) + 3 > 1603` 的配对会被官方脚本截断（`truncation='longest_first'`），逐对记录"是否截断"以及两条序列截断前后的长度。
- 环境：官方依赖中 `sentence_transformers`、`datasets`、`scikit-learn` 不在 `plminteract` 环境里。
  建一个继承 `plminteract` 的独立 venv（`--system-site-packages`），只装 `sentence-transformers==2.7.0`、`datasets==2.19.0`、
  `scikit-learn==1.4.1.post1`（官方 `environment.yml` 的版本；前两个官方没有固定版本），并用约束文件锁住已有的 torch、transformers 等版本。
  `plminteract` 环境本身不改动。
- 设置 `TORCH_FORCE_WEIGHTS_ONLY_LOAD=1`，让官方脚本里的 `torch.load` 以 weights-only 方式读 checkpoint
  （脚本代码不变；这三个文件都是纯 state dict，`plm_score.py` 已经用 weights-only 方式加载过）。
- 用途：一致性检查（3.4）；在长度分层中，作为长配对的"官方截断后分数"并列报告。**判定用 3.2 的分数。**

### 3.4 一致性检查

对两边都不截断的配对（≤ 1603 token），报告官方脚本与 `plm_score.py` 之间的最大差值，分别按单个顺序和两个顺序的平均计算。
容差 1e-3（待确认）。超出时报告并查原因（批内 padding 等）；判定所用分数不受影响。

### 3.5 打分范围

同一个有序蛋白对在同一模型、同一序列区段下只打一次分，各分析集共用。

| 分析集 | 配对 | 模型 | `plm_score.py` | 官方脚本 |
|---|---|---|---|---|
| H-M | 正样本、全部负样本、排名网格（配体池 × 受体池） | V11、V12、Leakage-Free | ✓ | ✓ |
| H-A | 正样本、全部负样本 | 三个 | ✓ | ✓ |
| H-X | 正样本、全部负样本（展开为亚基组合） | 三个 | ✓ | ✓ |
| M-M | 正样本、全部负样本、排名网格 | V11 | ✓ | ✓ |
| M-A | 正样本、全部负样本 | V11 | ✓ | ✓ |
| H-M 胞外段 | 两个蛋白都能取到胞外段的正样本和负样本 | V11 | ✓ | — |

胞外段用 `proteins.py` 的 `region_sequence(…, 'extracellular')` 规则（单一胞外拓扑域 → 单次跨膜 I 型的信号肽到跨膜段之间 → 最长成熟链 → 去信号肽的全长）。
多次跨膜蛋白按该规则拒绝，含有这类蛋白的配对不进入胞外段比较，列出数量。胞外段不做排名。

### 3.6 运行方式

- GPU（RTX 6000 Ada 48 GB，目前空闲）；PDE 队列只用 CPU，以 SCHED_IDLE 运行，不停、不调。
- 同时最多 3 个打分进程，每个 `OMP_NUM_THREADS=2`、`nice 10`；作业按每块 5000 对切分，`setsid nohup` 运行，每块一个 manifest（git HEAD、输入 SHA256、checkpoint SHA256）。
- 冒烟测试测出吞吐量后，在报告里给出总耗时估计。

## 4. 指标

每项都给 95% bootstrap 置信区间：2000 次重抽样，种子 1，百分位法。
AUROC/AUPR 把正样本和负样本分别有放回地重抽；排名指标按配体重抽；热门蛋白诊断按蛋白重抽。
配对之间共享蛋白、并不独立，所以这些区间偏窄。

1. **AUROC、AUPR。** AUROC 用 Mann–Whitney，平局计 1/2。AUPR 用平均精度（阶梯式），并给出随机基线（正样本比例：随机负样本为 1/11 ≈ 0.091，保持度数负样本为 0.5）。
   随机负样本和保持度数负样本分开报告；保持度数的三个种子分别报告，并给出三者的平均。
2. **分层（全部描述性）：**
   - 类别：第 k 类的正样本 vs 配体属于第 k 类配体池的负样本；
   - 长度：H-A 中两个蛋白都 ≤ 800 aa vs 其他（正负样本都按自身长度分组）；"其他"组同时列出官方截断后的分数；
   - 全长 vs 胞外段：H-M 中两个蛋白都能取到胞外段的同一批配对，V11；
   - 三个模型版本；
   - 复合物：H-X 全部，以及"至少一侧是复合物"的正样本 vs "至少一侧是复合物"的负样本；
   - 小鼠：M-M、M-A；
   - 训练集中见过 / 没见过（第 7 条）。
3. **排名**（H-M、M-M，每个模型）：对每个配体，候选受体 = 受体池，去掉配体本身，以及与该配体在 E 中有注释、但不是本集正样本的受体（这些既不算命中也不算竞争者）。
   真实受体 = 该配体在本集的正样本受体。最好排名 = 1 + 分数 ≥ 最高真实受体分数的非真实候选数（平局算输）。
   报告 top-1、top-5 命中率和最好排名的中位数，以及随机打分下的期望命中率（每个配体 1 − C(n−m, k)/C(n, k) 的平均，n 为候选数，m 为真实受体数）。
4. **分数分布：** 正样本、随机负样本、保持度数负样本的分数直方图；各组分数 > 0.5 的比例。同源配对单独画一组。
5. **热门蛋白诊断**（H-M，三个模型）：每个蛋白在 H-M 排名网格里全部未注释组合（不在 E 中、不是自身）上的平均分数，
   与该基因在人类 CellChatDB 中出现的相互作用行数（四个类别、作为配体或受体、包括作为复合物亚基）之间的 Spearman 相关；
   同时报告与它在 H-M 正样本中的度数的 Spearman 相关。
6. **参照：** 把 2026-11-14 已复现的自检结果放进同一张表。注意 V11 的 AUPR 0.913 / AUROC 0.980 是作者**小鼠跨物种测试集**
   （`mouse.ppi.qrels.seq.test.csv` 抽 2000 对，阳性 179，种子 2），不是人类测试集。同表列出 V12（0.838 / 0.955）、
   Leakage-Free 小鼠（0.238 / 0.636）和 Leakage-Free 自己的 Bernett 测试集（0.661 / 0.680）。
7. **训练集重叠：** 每个模型对照它自己的训练数据（V11：`cross_species_benchmarking` 的 train 和 val；V12：`STRING_V12_TrainingSet`；
   Leakage-Free：`Bernett_benchmarking` 的 train 和 val），按两条序列完全一致匹配无序蛋白对。
   - 每个配对标为：作为阳性见过 / 作为阴性见过 / 配对没见过但两个蛋白都见过 / 至少一个蛋白没见过。
   - 分"作为阳性见过"和"没见过"两组报告正样本的 AUROC（负样本相同，另报"只用没见过的负样本"）。
   - 负样本中作为阳性见过的数量单独列出。
   - 训练序列可能来自旧版 Ensembl/STRING，与当前 UniProt canonical 不一定完全一致，所以"见过"可能被低估。拿不到某个训练集时注明，不阻塞其余部分。

## 5. 判定标准（用户建议值，待确认）

- **主分析：** humanV11、H-M（单基因、非同源、两个蛋白都 ≤ 800 aa）、全长序列、`plm_score.py` 分数（两个顺序平均）。
- **"可用于区分 LR 配对"需同时满足三条：**
  1. 随机负样本 AUROC ≥ 0.85；
  2. 保持度数负样本 AUROC ≥ 0.75；
  3. 配体 top-5 命中率 ≥ 50%。
- 其他所有分层、模型版本和物种都是描述性的，不作判定依据，也不能弥补主分析的不通过。

为执行而补充的细节（我的建议，待确认）：
- 三条都按**点估计**判断，置信区间只报告；
- 第 2 条用三个种子 AUROC 的**平均值**，同时报告每个种子是否各自达到 0.75；
- 第 3 条的排名按 4.3 的定义（平局算输，去掉在其他来源中有注释的受体）。

## 6. 输出

- `reports/plm_lr_benchmark/dataset/`：第 1.5 节的数量表和第 2.5 节的文件。
- `reports/plm_lr_benchmark/smoke/`：冒烟测试。
- `reports/plm_lr_benchmark/`：`report.md`（数量表、主分析、分层表、热门蛋白诊断、一致性检查、训练集重叠、所有偏离计划之处），
  以及逐对分数（gzip CSV）、指标表、直方图。
- 原始下载放在 `data/lr_benchmark/`（不提交），打分中间文件放在 `plm_lr_benchmark_runs/`（不提交）。

## 7. 执行顺序

1. 提交本计划，用户推送。
2. **确认前可做：**
   - 构建数据集：下载、映射、正样本、排除集、负样本、数量统计，提交；
   - 冒烟测试：**只用 2026-11-14 已经打过分、看过结果的配对**（`diagnostics/runs/positive_controls_01`，15 × 15 组合 × 两个顺序 = 450 个作业，
     含超过 1603 token 的长配对）。三个模型各跑一次官方脚本，与已有的 `plm_score.py` 分数比较；测吞吐量；检查截断记录。
     **确认之前不给任何评测配对打分。**
3. 用户确认第 5 节和第 8 节 → 把确认结果和评测代码一起提交 → 正式打分 → 评价 → 写报告 → 提交，提醒推送。

## 8. 待确认或待决定的事项

1. 第 5 节的三个阈值，以及点估计、种子平均、排名定义这三个执行细节。
2. 同源配对（L = R）不进入主分析（1.2 第 5 条）。
3. 排除集包含 CellChatDB cofactor 关系（2.1 第 2 条）。
4. 类别分层按"配体所属类别"给负样本分组（4.2）。
5. 一致性容差 1e-3（3.4）。
6. 小鼠只跑 humanV11；胞外段只跑 H-M、只用 `plm_score.py`、不做排名（3.5）。
7. 映射后未解决的基因（包括 NRXN1 这类有多个 reviewed 条目的）：默认去掉相关配对；数据集构建后列出清单，由用户决定是否提供 override。
