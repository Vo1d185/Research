# PLM-interact × CellChatDB 评测：数据集、冒烟测试和参照值

2026-11-16。计划见 `reports/plm_lr_benchmark_plan.md`（`88d1fe1`，用户确认见第 9 节，`55e9027`）。
代码：`lr_benchmark.py`、`plm_official.py`、`lr_eval.py`、`proteins.py`（扩展）。
本文件的数字是按第 9.3 节同义名规则重新映射后的数据集（在 `4a0bea9` 上生成）。写这份文件时，**还没有给任何评测配对打分**。
第一版数据集（只用主基因名规则，`b8ca0e8`）的数字在 git 历史里。

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
| **H-M（主）** | **605** | 452 | 140 | 13 | 0 | 527 | 100%（按定义） |
| H-A（单基因，不限长度） | 1127 | 619 | 349 | 159 | 0 | 818 | 605（53.7%） |
| H-X（含复合物） | 2104 | 1269 | 411 | 424 | 977 | 1066 | 794（37.7%，全部亚基 ≤ 800） |
| 小鼠全部蛋白配对（规则之前） | 2377 | 1280 | 662 | 435 | 995 | — | — |
| M-M | 575 | 428 | 138 | 9 | 0 | 494 | 100% |
| M-A | 1080 | 582 | 349 | 149 | 0 | 774 | 575（53.2%） |

- 同源配对（L = R）不进 AUROC，单独描述性报告：人类 88 对，小鼠 95 对，其中映射成功、会打分的人类 88 对、小鼠 45 对（`homophilic_*.csv`）。
- 因基因映射不到唯一条目去掉的配对：人类 47（单基因 41），小鼠 317（单基因 257）。
  同义名规则救回人类 6 对、小鼠 37 对（第一版分别去掉 53 和 354）。
- 同一对在多行出现时合并；本库中没有出现跨类别或正反两个方向重复的情况。
- 主分析的 ECM 只有 13 对（ECM 配体多为长链胶原/层粘连蛋白，受体多为整合素复合物）。
- CellChatDB 有两行符号列与主列不一致：`CXCL8_CXCR1` 的 `ligand.symbol` 写成 CXCL1，`CD99_PILRB` 是小鼠大小写。正样本用主列；排除集两者都收。

## 3. 基因 → UniProt（主基因名规则 + 第 9.3 节同义名规则，没有手工指定）

- 人类 1148 个基因，映射 1141 个（主基因名 1138，via synonym 3：GPR1 → CMKLR2、NCL → 主基因名 "NUCLEOLIN" 的条目、NPY6R → NPY6RP）。
  **仍未解决 7 个**：CALCA、NRXN1、NRXN2、NRXN3、PRNP（主基因名多个条目）、MT-RNR2（7 个小肽，按确认直接去掉）、
  CCL3L3（它是 CCL3L1 条目中的第二个主基因名，不是同义名，严格按规则不采用）。
- 小鼠 1199 个基因，映射 1062 个（主基因名 1046，via synonym 16）；**仍未解决 137 个**：122 个无命中（多为只有 unreviewed 条目，以及 H2 基因），
  11 个主基因名多命中，1 个同义名多命中，3 个同义名指向的条目已被另一个基因按主基因名占用。
- 小鼠 via synonym 的 16 个：Cd244a → Cd244、Cd8b1 → Cd8b、Fasl → Faslg、Fcer2a → Fcer2、Flt3l → Flt3lg、Fpr3 → Fpr-s1、Gh → Gh1、
  Hc → C5、Icosl → Icoslg、Igfl3 → Igfl、Il4ra → Il4r、Kitl → Kitlg、Klk5 → Klk1b5、Ncl → 主基因名 "Nucleolin" 的条目、Ntf5 → Ntf4、Pilrb1 → Pilrb。
  **其中 Klk5 → Klk1b5 值得怀疑**：Klk1b5 是组织激肽释放酶家族的另一个成员，只是把 "Klk5" 列为同义名；按规则保留，只影响小鼠补充分析。
- 映射表中 `resolution` 列标注 `via synonym`；明细：`dataset/sources/protein_map_*.csv`、`dataset/sources/unresolved_{human,mouse}.csv`。

## 4. 负样本（`dataset/<集>/`）

- 排除集 E：符号层面 4565 对，accession 层面 1943 对。来源：CellChatDB 人类和小鼠（四类，复合物拆亚基，含 cofactor 关系）、
  CellPhoneDB v5 的 1878 条蛋白–蛋白相互作用（`is_ppi = True`；另 1033 条小分子相互作用按计划忽略）。
- 各集配体池 × 受体池中有注释、但不是本集正样本的组合（排名时去掉）：H-M 51 对（其中 43 对只来自 CellPhoneDB）。

| 集 | 配体池 × 受体池 | 排除后候选 | 随机负样本（1:10，种子 1） | 保持度数（1:1，种子 1/2/3） |
|---|---|---|---|---|
| H-M | 270 × 274（排名网格 73,963 对） | 73,171 | 6,050 | 605 × 3 |
| H-A | 455 × 389 | 175,436 | 11,270 | 1,127 × 3 |
| H-X | 608 × 548（实体） | 328,210 | 21,040 | 2,104 × 3 |
| M-M | 238 × 271（网格 64,483 对） | 63,773 | 5,750 | 575 × 3 |
| M-A | 417 × 384 | 158,606 | 10,800 | 1,080 × 3 |

保持度数的打乱：每套都做满 10 × 边数次成功交换，没有残留在 E 中的边，度数与正样本完全相同（`build.json`）。
同一提交上重新生成，结果逐字节相同。`build.json` 中的 `tracked_changes_paths` 只有本目录的旧文件（生成时已移走）。

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

## 6. 工作量（`plm_lr_benchmark_runs/jobs/jobs.json`，按第 9.4 节：官方脚本不跑排名网格）

| 作业集 | 无序蛋白对 | 作业 | 超过 1603 token | 打分 |
|---|---|---|---|---|
| `human_full`：H-M（含网格）、H-A、H-X、同源配对 | 110,232 | 220,376 | 23,220 | `plm_score.py`，三个模型 |
| `human_pn`：同上，不含网格 | 55,101 | 110,114 | 23,220 | 官方脚本，三个模型 |
| `mouse_full`：M-M（含网格）、M-A、同源配对 | 73,224 | 146,403 | 6,991 | `plm_score.py`，humanV11 |
| `mouse_pn`：同上，不含网格 | 21,524 | 43,003 | 6,991 | 官方脚本，humanV11 |
| `human_ecd`：H-M 正负样本的胞外段 | 3,664 | 7,328 | 0 | `plm_score.py`，humanV11 |

- 同源配对每对一个作业（两个顺序相同）。
- 胞外段：H-M 正负样本的 8,276 对中只有 3,664 对（44%）两个蛋白都能取到胞外段；176 个蛋白被规则拒绝（多次跨膜，大多是 GPCR）。

## 7. 参照值（`reference/`）

| 数据 | 模型 | 抽样 | AUPR | AUROC |
|---|---|---|---|---|
| 作者小鼠跨物种测试集（`mouse.ppi.qrels.seq.test.csv`，2000 对抽样，阳性 179，种子 2） | humanV11 | 2026-11-14 | 0.913 | 0.980 |
| 同上 | humanV12 | 2026-11-14 | 0.838 | 0.955 |
| 同上 | Leakage-Free | 2026-11-14 | 0.238 | 0.636 |
| Bernett 金标准测试集（Leakage-Free 自己的基准，2000 对，阳性 988） | Leakage-Free | 2026-11-14 | 0.661 | 0.680 |
| **作者人类验证集**（HF `cross_species_benchmarking` 的 `val/human.ppi.qrels.seq.test.csv`，2000 对，阳性 196，种子 2） | humanV11 | 2026-11-16 | 0.948 | 0.991 |

作者的人类数据只有 `train/` 和 `val/`，没有单独的人类测试集；上表最后一行位于 `val/` 目录，可能参与过模型选择，不能当作独立测试集。
