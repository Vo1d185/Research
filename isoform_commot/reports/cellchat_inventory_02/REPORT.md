# CellChatDB 扩展盘点：数据集 × 阈值 × 靶基因库 × 复合物规则

2026-11-14。只计数，**不选面板，不跑 OT**。代码 `cellchat_inventory_grid.py`；明细见本目录 CSV。
表中每格为"配对数 / 不同受体数"（复合物受体算一个受体）。

## 0. 数据与来源

| 项目 | 内容 |
|---|---|
| Stereo-seq | `data/mousebrain_bin60.h5ad`，`layer:count`，整张组织 7765 spot、25691 个基因（与盘点 01 相同） |
| Visium | 10x `V1_Mouse_Brain_Sagittal_Posterior`，Space Ranger 1.1.0，Spatial 3' v1，mm10；3355 个组织内 spot，32285 个特征 |
| Visium 来源 | `https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/` |
| Visium SHA256 | `…_filtered_feature_bc_matrix.h5`：`df962209143860488bd679ff7a1b3bb14985208d0929fb814934ba7571c488b6`；`…_spatial.tar.gz`：`b878d7eaca487ffc17224a870f86599abf9d45a99033c1ef65661ff119cef0d5`（本地首次下载时计算，非 10x 发布） |
| 为什么是这份 | COMMOT `d117445` 的 `docs/notebooks/Basic_usage.ipynb` 和 `visium-mouse_brain.ipynb` 都用 `sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')` 读取；scanpy 1.10.3 对 `V1_` 样本使用 Space Ranger 1.1.0 的地址 |
| CellChatDB | 与盘点 01 相同（CellChat v2 `75253cd`，`dd59b196…`） |
| scSeqComm | 与盘点 01 相同（`95c9bb2`，KEGG + Reactome PPR，小鼠符号） |
| CCCvelo 库 | `Databases/Database.rar`，commit `fac213c600e7929b5d14e7d95dbdaabd722f9bc9`，SHA256 `0afe64a0f0c7c5db14fba2a8c50e031d534233556a03254c46a2c1e9a3c016df`（RAR5，用 bsdtar 3.7.7 / libarchive 解压） |

**Database.rar 的内容**

| 文件 | 大小 | SHA256 | 内容 |
|---|---|---|---|
| `Database/LigRecDB.csv` | 84,013 B | `2f10e3bda1e5…` | 3992 个 LR 对；与已固定的 `diff_LigRecDB.csv` + `cont_LigRecDB.csv` **完全相同** |
| `Database/RecTFDB.csv` | 32,202,572 B | `f145b3370b42…` | 872,490 行 = 762 个受体 × 1145 个 TF，每行一个正分数（稠密表） |
| `Database/TFTGDB.csv` | 10,531,843 B | `584be2b7f107…` | 577,019 条 TF→靶基因边，1145 个 TF，23,648 个靶基因；分数 1–4（未使用） |

完整哈希在 `data/cccvelo_db/db.json`。CCCvelo 库用大写符号，按小鼠符号转大写匹配（符号匹配，不是同源映射）。

**Visium 的重复基因符号：** 38 个符号对应多个特征，按"同符号相加"处理（任一副本检测到即算检测到）。
其中 CellChat 用到的有 Ccl19、Ccl21a、Ccl21b、Ccl27a、Il11ra2、Pcdha11。COMMOT 的 notebook 用
`var_names_make_unique`，只保留第一个副本的原名，因此这几个基因的检测比例可能与 COMMOT 略有不同。

## 1. 规则

| 维度 | 取值 |
|---|---|
| 阈值 | 1%（CCCvelo 皮层脚本 `select_LRTG(..., pct.ct = 0.01)`）和 5%。**同一阈值同时用于配体、受体和靶基因**（CCCvelo 的 pct.ct 也同时筛 L、R 和 TG） |
| 靶基因库 | scSeqComm；CCCvelo RecTFDB + TFTGDB |
| 复合物规则 | `single_gene_only`：盘点 01 的规则，含复合物的配对在 (a) 去掉 |
| | `complex_any_subunit_tf`：复合物表达取各亚基最小值；(c1) 只要**任一亚基**有 TF 关联 |
| | `complex_all_subunits_tf`：同上，但 (c1) 要求**全部亚基**都有 TF 关联 |

步骤：(0) 全部 →(a) 单基因（只对 `single_gene_only` 起作用，复合物规则下 (a) = 全部）→(b) 配体和受体检测比例 ≥ 阈值
→(c1) 受体有 TF 关联 →(c) 前 10 个 TF 的靶基因并集，去掉所有通过 (b) 的配对（两类合并、含亚基）中的基因，
检测比例 ≥ 阈值的靶基因 ≥5 个。复合物受体的 TF 排序：在所有亚基的关联中，每个 TF 取最大权重，再取前 10 个。

**复合物的检测比例（主表）：** 按你给的定义，先对每个 spot 取各亚基计数的最小值，再算最小值 > 0 的 spot 比例，
也就是"全部亚基在同一个 spot 中都检测到"。

**与 COMMOT 实际代码的差别：** COMMOT（`d117445`）运输时确实用亚基最小值作为复合物表达量
（`commot/tools/_spatial_communication.py` 第 49 行 `heteromeric_rule='min'`；配体第 101 行、受体第 108 行取最小值），
但它的配对筛选函数 `filter_lr_database`（`commot/preprocessing/_ligand_receptor_database.py` 第 125–128 行）
不是按最小值算比例，而是要求**每个亚基各自**的检测比例 ≥ 阈值。后者更宽松，结果放在第 3 节对照。

**一致性检查：** Stereo-seq / 5% / scSeqComm / `single_gene_only` 一行与盘点 01 完全一致
（Secreted 23 / 17 → 12 / 9 → 8 / 6；Contact 110 / 59 → 41 / 25 → 14 / 10）。

## 2. 主表

(0) 全部：Secreted 1280 / 385，Contact 662 / 275（与其他维度无关）。

### Secreted Signaling（复合物检测：同一 spot 共检出）

| 数据集 | 阈值 | 靶基因库 | 复合物规则 | a | b | c1 | c |
|---|---|---|---|---|---|---|---|
| Stereo-seq bin60 | 1% | CCCvelo | single_gene_only | 628 / 250 | 84 / 47 | 81 / 45 | 80 / 44 |
| Stereo-seq bin60 | 1% | CCCvelo | complex_any_subunit_tf | 1280 / 385 | 122 / 62 | 119 / 60 | 118 / 59 |
| Stereo-seq bin60 | 1% | CCCvelo | complex_all_subunits_tf | 1280 / 385 | 122 / 62 | 119 / 60 | 118 / 59 |
| Stereo-seq bin60 | 1% | scSeqComm | single_gene_only | 628 / 250 | 84 / 47 | 61 / 35 | 53 / 29 |
| Stereo-seq bin60 | 1% | scSeqComm | complex_any_subunit_tf | 1280 / 385 | 122 / 62 | 99 / 50 | 91 / 44 |
| Stereo-seq bin60 | 1% | scSeqComm | complex_all_subunits_tf | 1280 / 385 | 122 / 62 | 90 / 47 | 82 / 41 |
| Stereo-seq bin60 | 5% | CCCvelo | single_gene_only | 628 / 250 | 23 / 17 | 22 / 16 | 22 / 16 |
| Stereo-seq bin60 | 5% | CCCvelo | complex_any_subunit_tf | 1280 / 385 | 24 / 18 | 23 / 17 | 23 / 17 |
| Stereo-seq bin60 | 5% | CCCvelo | complex_all_subunits_tf | 1280 / 385 | 24 / 18 | 23 / 17 | 23 / 17 |
| Stereo-seq bin60 | 5% | scSeqComm | single_gene_only | 628 / 250 | 23 / 17 | 12 / 9 | 8 / 6 |
| Stereo-seq bin60 | 5% | scSeqComm | complex_any_subunit_tf | 1280 / 385 | 24 / 18 | 13 / 10 | 8 / 6 |
| Stereo-seq bin60 | 5% | scSeqComm | complex_all_subunits_tf | 1280 / 385 | 24 / 18 | 13 / 10 | 8 / 6 |
| Visium posterior | 1% | CCCvelo | single_gene_only | 628 / 250 | 210 / 107 | 201 / 99 | 200 / 98 |
| Visium posterior | 1% | CCCvelo | complex_any_subunit_tf | 1280 / 385 | 428 / 162 | 419 / 154 | 418 / 153 |
| Visium posterior | 1% | CCCvelo | complex_all_subunits_tf | 1280 / 385 | 428 / 162 | 418 / 153 | 417 / 152 |
| Visium posterior | 1% | scSeqComm | single_gene_only | 628 / 250 | 210 / 107 | 158 / 73 | 137 / 60 |
| Visium posterior | 1% | scSeqComm | complex_any_subunit_tf | 1280 / 385 | 428 / 162 | 376 / 128 | 355 / 115 |
| Visium posterior | 1% | scSeqComm | complex_all_subunits_tf | 1280 / 385 | 428 / 162 | 354 / 122 | 333 / 109 |
| Visium posterior | 5% | CCCvelo | single_gene_only | 628 / 250 | 129 / 75 | 122 / 69 | 121 / 68 |
| Visium posterior | 5% | CCCvelo | complex_any_subunit_tf | 1280 / 385 | 199 / 99 | 192 / 93 | 191 / 92 |
| Visium posterior | 5% | CCCvelo | complex_all_subunits_tf | 1280 / 385 | 199 / 99 | 192 / 93 | 191 / 92 |
| Visium posterior | 5% | scSeqComm | single_gene_only | 628 / 250 | 129 / 75 | 95 / 52 | 82 / 43 |
| Visium posterior | 5% | scSeqComm | complex_any_subunit_tf | 1280 / 385 | 199 / 99 | 165 / 76 | 152 / 67 |
| Visium posterior | 5% | scSeqComm | complex_all_subunits_tf | 1280 / 385 | 199 / 99 | 148 / 71 | 135 / 62 |

### Cell-Cell Contact（复合物检测：同一 spot 共检出）

| 数据集 | 阈值 | 靶基因库 | 复合物规则 | a | b | c1 | c |
|---|---|---|---|---|---|---|---|
| Stereo-seq bin60 | 1% | CCCvelo | single_gene_only | 591 / 244 | 202 / 95 | 131 / 63 | 131 / 63 |
| Stereo-seq bin60 | 1% | CCCvelo | complex_any_subunit_tf | 662 / 275 | 213 / 100 | 136 / 68 | 136 / 68 |
| Stereo-seq bin60 | 1% | CCCvelo | complex_all_subunits_tf | 662 / 275 | 213 / 100 | 136 / 68 | 136 / 68 |
| Stereo-seq bin60 | 1% | scSeqComm | single_gene_only | 591 / 244 | 202 / 95 | 96 / 43 | 96 / 43 |
| Stereo-seq bin60 | 1% | scSeqComm | complex_any_subunit_tf | 662 / 275 | 213 / 100 | 99 / 46 | 99 / 46 |
| Stereo-seq bin60 | 1% | scSeqComm | complex_all_subunits_tf | 662 / 275 | 213 / 100 | 99 / 46 | 99 / 46 |
| Stereo-seq bin60 | 5% | CCCvelo | single_gene_only | 591 / 244 | 110 / 59 | 71 / 43 | 71 / 43 |
| Stereo-seq bin60 | 5% | CCCvelo | complex_any_subunit_tf | 662 / 275 | 112 / 61 | 73 / 45 | 73 / 45 |
| Stereo-seq bin60 | 5% | CCCvelo | complex_all_subunits_tf | 662 / 275 | 112 / 61 | 73 / 45 | 73 / 45 |
| Stereo-seq bin60 | 5% | scSeqComm | single_gene_only | 591 / 244 | 110 / 59 | 41 / 25 | 14 / 10 |
| Stereo-seq bin60 | 5% | scSeqComm | complex_any_subunit_tf | 662 / 275 | 112 / 61 | 41 / 25 | 14 / 10 |
| Stereo-seq bin60 | 5% | scSeqComm | complex_all_subunits_tf | 662 / 275 | 112 / 61 | 41 / 25 | 14 / 10 |
| Visium posterior | 1% | CCCvelo | single_gene_only | 591 / 244 | 303 / 130 | 190 / 80 | 190 / 80 |
| Visium posterior | 1% | CCCvelo | complex_any_subunit_tf | 662 / 275 | 339 / 146 | 217 / 96 | 217 / 96 |
| Visium posterior | 1% | CCCvelo | complex_all_subunits_tf | 662 / 275 | 339 / 146 | 217 / 96 | 217 / 96 |
| Visium posterior | 1% | scSeqComm | single_gene_only | 591 / 244 | 303 / 130 | 147 / 54 | 147 / 54 |
| Visium posterior | 1% | scSeqComm | complex_any_subunit_tf | 662 / 275 | 339 / 146 | 170 / 68 | 170 / 68 |
| Visium posterior | 1% | scSeqComm | complex_all_subunits_tf | 662 / 275 | 339 / 146 | 170 / 68 | 170 / 68 |
| Visium posterior | 5% | CCCvelo | single_gene_only | 591 / 244 | 240 / 103 | 169 / 72 | 169 / 72 |
| Visium posterior | 5% | CCCvelo | complex_any_subunit_tf | 662 / 275 | 260 / 113 | 184 / 82 | 184 / 82 |
| Visium posterior | 5% | CCCvelo | complex_all_subunits_tf | 662 / 275 | 260 / 113 | 184 / 82 | 184 / 82 |
| Visium posterior | 5% | scSeqComm | single_gene_only | 591 / 244 | 240 / 103 | 129 / 49 | 129 / 49 |
| Visium posterior | 5% | scSeqComm | complex_any_subunit_tf | 662 / 275 | 260 / 113 | 141 / 57 | 141 / 57 |
| Visium posterior | 5% | scSeqComm | complex_all_subunits_tf | 662 / 275 | 260 / 113 | 141 / 57 | 141 / 57 |

## 3. 复合物检测比例：同一 spot 共检出 vs COMMOT 的逐亚基规则

`complex_any_subunit_tf` 规则下：

| 类别 | 数据集 | 阈值 | 靶基因库 | b 同一 spot 共检出 | b 逐亚基（COMMOT） | c 同一 spot 共检出 | c 逐亚基（COMMOT） |
|---|---|---|---|---|---|---|---|
| Contact | Stereo-seq bin60 | 1% | CCCvelo | 213 / 100 | 225 / 107 | 136 / 68 | 146 / 75 |
| Contact | Stereo-seq bin60 | 1% | scSeqComm | 213 / 100 | 225 / 107 | 99 / 46 | 108 / 53 |
| Contact | Stereo-seq bin60 | 5% | CCCvelo | 112 / 61 | 112 / 61 | 73 / 45 | 73 / 45 |
| Contact | Stereo-seq bin60 | 5% | scSeqComm | 112 / 61 | 112 / 61 | 14 / 10 | 14 / 10 |
| Contact | Visium posterior | 1% | CCCvelo | 339 / 146 | 347 / 151 | 217 / 96 | 225 / 101 |
| Contact | Visium posterior | 1% | scSeqComm | 339 / 146 | 347 / 151 | 170 / 68 | 178 / 73 |
| Contact | Visium posterior | 5% | CCCvelo | 260 / 113 | 271 / 118 | 184 / 82 | 192 / 87 |
| Contact | Visium posterior | 5% | scSeqComm | 260 / 113 | 271 / 118 | 141 / 57 | 149 / 62 |
| Secreted | Stereo-seq bin60 | 1% | CCCvelo | 122 / 62 | 169 / 80 | 118 / 59 | 165 / 77 |
| Secreted | Stereo-seq bin60 | 1% | scSeqComm | 122 / 62 | 169 / 80 | 91 / 44 | 138 / 62 |
| Secreted | Stereo-seq bin60 | 5% | CCCvelo | 24 / 18 | 26 / 20 | 23 / 17 | 25 / 19 |
| Secreted | Stereo-seq bin60 | 5% | scSeqComm | 24 / 18 | 26 / 20 | 8 / 6 | 8 / 6 |
| Secreted | Visium posterior | 1% | CCCvelo | 428 / 162 | 542 / 182 | 418 / 153 | 532 / 173 |
| Secreted | Visium posterior | 1% | scSeqComm | 428 / 162 | 542 / 182 | 355 / 115 | 469 / 135 |
| Secreted | Visium posterior | 5% | CCCvelo | 199 / 99 | 269 / 115 | 191 / 92 | 261 / 108 |
| Secreted | Visium posterior | 5% | scSeqComm | 199 / 99 | 269 / 115 | 152 / 67 | 222 / 83 |

通过 (b) 的配对中含复合物的数量（共检出 → 逐亚基）：Stereo-seq 1% 时 Secreted 38 → 85、Contact 11 → 23；
5% 时 Secreted 1 → 3、Contact 2 → 2。Visium 1% 时 Secreted 218 → 332、Contact 36 → 44；5% 时 Secreted 70 → 140、Contact 20 → 31。

## 4. 从表中能直接读出的情况（不涉及选择）

- **数据集：** Visium 通过 (b) 的配对远多于 Stereo-seq bin60。例如 Secreted、5%、单基因：Visium 129 对，Stereo-seq 23 对。
- **阈值：** Stereo-seq 上从 5% 降到 1%，Secreted 单基因 (b) 从 23 对增加到 84 对，Contact 从 110 对增加到 202 对。
- **复合物：** 在 Stereo-seq 5% 下，复合物几乎不增加配对（Secreted +1、Contact +2），因为多个亚基很少在同一个 bin60 spot 中同时检出；
  在 Visium 1% 下，Secreted 增加 218 对。
- **任一亚基 vs 全部亚基：** Contact 在所有组合下完全相同；Secreted 只在部分组合有差别，主要是 scSeqComm
  （如 Visium 1% 的 (c)：355 / 115 vs 333 / 109；Stereo-seq 1%：91 / 44 vs 82 / 41）。
- **靶基因库：** 用 CCCvelo 库时 (c) 几乎等于 (c1)。CCCvelo 的 RecTFDB 是稠密表，(c1) 实际上只检查受体是否在它的 762 个受体中；
  而 TFTGDB 很大，前 10 个 TF 的靶基因并集在各组合中位数为 777–3036 个，"≥5 个"这一条从不起作用。
  scSeqComm 的中位数为 4–98 个。用 CCCvelo 库时，靶基因得分会是几百到几千个基因的平均表达。
- **scSeqComm 在 Stereo-seq 5% 下 (c1)→(c) 下降明显**（Contact 41 → 14），因为靶基因也要满足 5% 检测；1% 时 Contact 没有下降（96 → 96）。

## 5. 限制

- 检测比例按整张组织 / 整张切片计算；只看检测比例，没有用平均表达。
- Visium spot 直径 55 µm、中心间距 100 µm，包含多个细胞；Stereo-seq bin60 的物理尺寸还未核实。两者的检测比例不能直接等同。
- 靶基因的检测阈值与 L/R 阈值相同，这是本次的设定；TFTGDB 的分数（1–4）没有使用。
- CCCvelo 库与小鼠数据之间是大写符号匹配。
- 复合物受体的 TF 排序（各亚基取最大权重）是本次的设定，COMMOT 和 CCCvelo 都没有对应规则。

## 6. 文件

| 文件 | 内容 |
|---|---|
| `steps.csv` | 全部组合（含 `detection = per_subunit`）每一步的配对数、受体数、配体数和共享数量 |
| `pairs_passing_b.csv` | 每个组合中通过 (b) 的配对：检测比例、是否通过 (c1)/(c)、靶基因数 |
| `grid_tables.md`、`commot_filter_table.md` | 第 2、3 节的表（由 CSV 自动生成） |
| `inventory.json` | 数据集规模、重复符号、来源哈希、规则 |

```bash
python cellchat_inventory_grid.py fetch-cccvelo-db
python cellchat_inventory_grid.py fetch-visium
python cellchat_inventory_grid.py report --output reports/cellchat_inventory_02
```
