# PLM-interact × CellChatDB 评测：确认前的准备工作

2026-11-16。计划见 `reports/plm_lr_benchmark_plan.md`（`88d1fe1`）。代码 `41bd254`：`lr_benchmark.py`、`plm_official.py`、`lr_eval.py`、`proteins.py`（扩展）。
这里只有数据集和冒烟测试，**没有给任何评测配对打分**。判定标准仍待确认。

## 1. 数据来源（SHA256 为本地首次下载时计算）

| 数据 | 版本 | SHA256 |
|---|---|---|
| CellChatDB v2 人类 `CellChatDB.human.rda` | jinworks/CellChat `75253cd` | `582e99db2bb1b1e05c0364553893d804952bd007fef870720a9cc659a7bd3ad8` |
| CellChatDB v2 小鼠 `CellChatDB.mouse.rda` | 同上（已有） | `dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413` |
| CellPhoneDB v5 `interaction_input.csv` 等 4 个文件 | ventolab/cellphonedb-data tag v5.0.0 = `71ffa8a` | 见 `dataset/sources/downloads.json` |
| 训练集（humanV11 / V12 / Leakage-Free） | HF `cross_species_benchmarking@19dbed2`、`STRING_V12_TrainingSet@819f938`、`Bernett_benchmarking@5d2ad03` | 同上 |
| UniProt | release 2026_04 | 每个条目的版本见 `dataset/sources/protein_map_*.csv` |

## 2. 正样本数量（`dataset/counts.csv`）

| 分析集 | 正样本 | Secreted | Contact | ECM | 含复合物 | 去重蛋白数 | 两个蛋白都 ≤ 800 aa |
|---|---|---|---|---|---|---|---|
| 人类全部蛋白配对（规则之前） | 2239 | 1280 | 535 | 424 | 983 | — | — |
| **H-M（主）** | **599** | 446 | 140 | 13 | 0 | 524 | 100%（按定义） |
| H-A（单基因，不限长度） | 1121 | 613 | 349 | 159 | 0 | 815 | 599（53.4%） |
| H-X（含复合物） | 2098 | 1263 | 411 | 424 | 977 | 1063 | 788（37.6%，全部亚基 ≤ 800） |
| 小鼠全部蛋白配对（规则之前） | 2377 | 1280 | 662 | 435 | 995 | — | — |
| M-M | 554 | 419 | 126 | 9 | 0 | 478 | 100% |
| M-A | 1053 | 568 | 336 | 149 | 0 | 749 | 554（52.6%） |

- 去掉的同源配对（L = R）：人类 88，小鼠 95（`homophilic_*.csv`）。
- 因基因映射不到唯一条目去掉的配对：人类 53（单基因 47），小鼠 354（单基因 284）。
- 同一对在多行出现时合并；本库中没有出现跨类别或正反两个方向重复的情况。
- 主分析的 ECM 只有 13 对（ECM 配体多为长链胶原/层粘连蛋白，受体多为整合素复合物）。
- CellChatDB 有两行符号列与主列不一致：`CXCL8_CXCR1` 的 `ligand.symbol` 写成 CXCL1，`CD99_PILRB` 是小鼠大小写。正样本用主列；排除集两者都收。

## 3. 基因 → UniProt（规则与小鼠主流程相同，没有猜）

- 人类 1148 个基因，映射 1138 个；**未解决 10 个**：CALCA（2 个条目）、CCL3L3（条目的第一主基因名是 CCL3L1）、GPR1（现名 CMKLR2）、
  MT-RNR2（7 个小肽）、NCL（UniProt 主基因名为 "NUCLEOLIN"）、NPY6R（人类无 reviewed 条目）、NRXN1/2/3（α、β 两个条目）、PRNP（2 个条目）。
- 小鼠 1199 个基因，映射 1046 个；**未解决 153 个**（142 个无命中，11 个多命中）。大多是 CellChat 用的符号在 UniProt 只是同义名
  （如 Il4ra → Il4r、Kitl → Kitlg、Fasl → Faslg、Gh → Gh1），或只有 unreviewed 条目（如 Col6a3、Ifngr2），以及 H2 基因。
- 明细：`dataset/sources/unresolved_{human,mouse}.csv`。

## 4. 负样本（`dataset/<集>/`）

- 排除集 E：符号层面 4565 对，accession 层面 1943 对。来源：CellChatDB 人类和小鼠（四类，复合物拆亚基，含 cofactor 关系）、
  CellPhoneDB v5 的 1878 条蛋白–蛋白相互作用（`is_ppi = True`；另 1033 条小分子相互作用按计划忽略）。
- 各集配体池 × 受体池中有注释、但不是本集正样本的组合（排名时去掉）：H-M 51 对（其中 43 对只来自 CellPhoneDB）。

| 集 | 配体池 × 受体池 | 排除后候选 | 随机负样本（1:10，种子 1） | 保持度数（1:1，种子 1/2/3） |
|---|---|---|---|---|
| H-M | 270 × 271（排名网格 73,153 对） | 72,367 | 5,990 | 599 × 3 |
| H-A | 455 × 386 | 174,077 | 11,210 | 1,121 × 3 |
| H-X | 608 × 545（实体） | 326,392 | 20,980 | 2,098 × 3 |
| M-M | 231 × 262（网格 60,507 对） | 59,818 | 5,540 | 554 × 3 |
| M-A | 405 × 371 | 148,763 | 10,530 | 1,053 × 3 |

保持度数的打乱：每套都做满 10 × 边数次成功交换，没有残留在 E 中的边，度数与正样本完全相同（`build.json`）。
在提交 `41bd254` 上重新生成，结果与第一次生成逐字节相同。

## 5. 冒烟测试（`smoke/`，只用 2026-11-14 已经看过的 15 × 15 阳性对照，450 个作业）

| 模型 | 官方脚本截断 | 比较的作业 | 官方 vs `plm_score.py` 最大差值（单个顺序） | 两顺序平均的最大差值 |
|---|---|---|---|---|
| humanV11 | 40 | 410 | 1.0e-5 | 7e-6 |
| humanV12 | 40 | 410 | 2.5e-5 | 1.3e-5 |
| Leakage-Free | 40 | 410 | 3e-6 | 2e-6 |

- 三个模型都在 1e-3 容差内。被截断的 40 个作业都截到 1603 token。
- `plm_score.py` 今天重跑 V11，与 11-14 的分数完全相同（最大差值 0）。
- 官方脚本在单独的 venv（`/home/liu/venvs/plminteract-official`，继承 `plminteract`，加 sentence-transformers 2.7.0、datasets 2.19.0、
  scikit-learn 1.4.1.post1）中运行，脚本未改，`TORCH_FORCE_WEIGHTS_ONLY_LOAD=1`。官方仓库跟踪了 `__pycache__/*.pyc`，第一次运行时被 Python 重写，
  已恢复为提交版本，之后的运行不再写字节码。

## 6. 工作量（`plm_lr_benchmark_runs/jobs/jobs.json`）

| 作业集 | 无序蛋白对 | 作业（两个顺序） | 超过 1603 token | 模型 × 打分方式 |
|---|---|---|---|---|
| 人类全长（H-M 含网格、H-A、H-X） | 109,207 | 218,414 | 23,268（10.7%） | 3 × 2 |
| 小鼠全长（M-M 含网格、M-A） | 68,952 | 137,904 | 6,968 | 1 × 2 |
| 人类胞外段（H-M 正负样本） | 3,628 | 7,256 | 0 | 1 × 1 |

- 胞外段：H-M 正负样本的 8,194 对中只有 3,628 对（44%）两个蛋白都能取到胞外段；174 个蛋白被规则拒绝（173 个多次跨膜，大多是 GPCR）。
- 按冒烟测试的速度粗估：每个模型在人类全长上 `plm_score.py` 约 7 小时、官方脚本约 6 小时；总计约 50 GPU 小时，3 个进程并行估计 **1.5–2.5 天**。
