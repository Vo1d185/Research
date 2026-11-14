# CellChatDB 小鼠库在整张组织上的盘点

2026-11-14。只报告数量，**不选面板**；面板规则等确认后再定。
代码 `cellchat_inventory.py`；本目录的 CSV 是全部明细。

## 0. 数据与来源

| 项目 | 内容 |
|---|---|
| 表达数据 | `data/mousebrain_bin60.h5ad`，`layer:count`，**整张组织 7765 个 spot**（不是 128 spot 的 ROI） |
| CellChatDB | CellChat v2，jinworks/CellChat commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`；文件最后修改于 `592ba49cf4`（2025-04-16） |
| 文件 | `data/CellChatDB.mouse.rda`，SHA256 `dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413`（本地首次下载时计算，非作者发布）；GPL-3，不随仓库分发 |
| 规模 | 3379 条：Secreted Signaling 1280、Cell-Cell Contact 662、ECM-Receptor 435、Non-protein Signaling 1002 |
| 靶基因库 | scSeqComm `95c9bb2`（与 `evaluate_targets.py` 相同） |
| CCCvelo 表 | 已固定的 `diff_LigRecDB.csv` / `cont_LigRecDB.csv`（来源字段：stMLnet 3985 条，CellSig 7 条） |

**h5ad 的 `area` 列：** 只有一个取值 3600，7765 个 spot 全部相同。3600 = 60²，与 bin60（60×60 个 DNB 合并为一个 spot）
和 ROI 中位最近邻间距 h = 60 一致，说明坐标很可能以 DNB 为单位。Stereo-seq 的 DNB 中心间距通常是 500 nm
（部分芯片为 715 nm），据此 h 约为 30 µm（或约 43 µm）。**这是推断，需要对照数据来源的芯片说明确认**，
在确认前坐标单位仍记为 `source_unit_unverified`。

## 1. 规则

| 步骤 | 条件 |
|---|---|
| 0 | CellChatDB 中该类别的全部配对 |
| a | 配体和受体都是单基因；含复合物的配对单独计数（表 1） |
| b | 配体和受体在**整张组织** ≥5% 的 spot 中检测到（count > 0）。h5ad 中没有该基因也算不通过 |
| c1 | 受体在 scSeqComm（KEGG + Reactome）中有 receptor→TF 关联 |
| c | 按 `evaluate_targets.py` 的规则：PPR 最高的前 10 个 TF，取其靶基因并集，保留在 ≥5% spot 中检测到的，**≥5 个** |

(c) 中从靶基因里去掉的基因：所有通过 (b) 的配对（两类合并）中出现的配体和受体基因。以后从这个范围里选出的
任何面板，要去掉的基因都是它的子集，所以这里通过 (c) 的受体在任何这样的面板中都仍可评价。
作为敏感性分析，只去掉该配对自己的配体和受体：**结果完全相同**（两种去除方式下每一步数量一致）。

"共享"按同一类别内计算：共享配体 = 该配对的配体在剩余配对中出现 ≥2 次；受体同理；"任一"是两者之一。

## 2. 结果

**表 1：步骤 (a) 前后的复合物情况**

| 类别 | 单基因 | 仅配体是复合物 | 仅受体是复合物 | 两者都是 | 合计 |
|---|---|---|---|---|---|
| Secreted Signaling | 628 | 3 | 641 | 8 | 1280 |
| Cell-Cell Contact | 591 | 14 | 57 | 0 | 662 |

Secreted Signaling 中一半配对的受体是复合物（如多亚基细胞因子受体），步骤 (a) 去掉了其中约 51%。

**表 2：每一步剩下的配对**

| 类别 | 步骤 | 配对 | 不同受体 | 不同配体 | 共享配体 | 共享受体 | 共享任一 |
|---|---|---|---|---|---|---|---|
| Secreted | 0 全部 | 1280 | 385 | 386 | 1124 | 1125 | 1214 |
| | a 单基因 | 628 | 250 | 260 | 523 | 524 | 586 |
| | b 表达 ≥5% | **23** | 17 | 16 | 11 | 11 | 18 |
| | c1 有 TF 关联 | 12 | 9 | 10 | 4 | 5 | 8 |
| | c 靶基因 ≥5 | **8** | 6 | 7 | 2 | 3 | 5 |
| Contact | 0 全部 | 662 | 275 | 295 | 506 | 505 | 588 |
| | a 单基因 | 591 | 244 | 266 | 450 | 454 | 518 |
| | b 表达 ≥5% | **110** | 59 | 45 | 88 | 76 | 95 |
| | c1 有 TF 关联 | 41 | 25 | 24 | 25 | 25 | 37 |
| | c 靶基因 ≥5 | **14** | 10 | 12 | 3 | 6 | 9 |
| 两类合并 | b | 133 | 74 | 61 | 99 | 87 | 113 |
| | c | 22 | 16 | 19 | 5 | 9 | 14 |

两类合并后，(b) 有 2 个受体同时出现在两类中；(c) 的共享数量等于两类之和，没有跨类共享的实体。

**步骤 (b) 不通过的原因**

| 类别 | 通过 (a) | h5ad 中缺少其中一个基因 | 检测比例 <5% | 通过 (b) |
|---|---|---|---|---|
| Secreted | 628 | 116 | 489 | 23 |
| Contact | 591 | 89 | 392 | 110 |

缺少的基因（例如 Fgf21、Ins1、Ccl8）是没有测到/被上游过滤，还是符号不一致，这里没有区分。

供定规则时参考（**不是建议的阈值**）：只改变检测比例阈值时，通过 (b) 的配对数 / 不同受体数

| 阈值 | 1% | 2% | 5% | 10% |
|---|---|---|---|---|
| Secreted | 84 / 47 | 56 / 37 | 23 / 17 | 12 / 9 |
| Contact | 202 / 95 | 160 / 77 | 110 / 59 | 70 / 37 |

**步骤 (c) 不通过的原因**（在通过 (b) 的配对中）：Secreted 有 11 对受体没有 TF 关联、4 对靶基因不足；
Contact 有 69 对没有 TF 关联、27 对靶基因不足。很多黏附分子类受体不在 KEGG/Reactome 信号通路里，
所以 (c) 天然偏向经典信号受体。

**通过 (c) 的 22 对**（只是列出，不是面板）

| 类别 | 通路 | 配对（配体检测比例 / 受体检测比例，靶基因数） |
|---|---|---|
| Secreted | NRG | Nrg1–Erbb4（0.78 / 0.58，68）、Nrg2–Erbb4（0.23 / 0.58，68）、Nrg3–Erbb4（0.91 / 0.58，68） |
| | PTN | Ptn–Alk（0.27 / 0.15，31）、Ptn–Sdc2（0.27 / 0.08，34） |
| | PTPR | Ptprs–Ntrk3（0.38 / 0.49，65） |
| | EGF | Tgfa–Egfr（0.05 / 0.06，41） |
| | VEGF | Vegfa–Flt1（0.06 / 0.10，15） |
| Contact | NCAM | Ncam1–Ncam1（0.57 / 0.57，55）、Ncam1–Ncam2（0.57 / 0.54，55）、Ncam1–Fgfr1（0.57 / 0.06，42） |
| | PTPR / NGL | Lrfn5–Ptprf（0.52 / 0.11，18）、Lrrc4c–Ptprf（0.74 / 0.11，18）、Lrrc4b–Ptprf（0.12 / 0.11，18） |
| | CADM / NECTIN | Cadm1–Nectin3（0.59 / 0.11，14）、Cadm3–Nectin3（0.19 / 0.11，14）、Nectin1–Nectin3（0.08 / 0.11，14） |
| | FLRT | Flrt2–Flrt2（0.26 / 0.26，47） |
| | PTPRM | Ptprm–Ptprm（0.36 / 0.36，14） |
| | JAM | Jam3–Jam3（0.09 / 0.09，15） |
| | CNTN | Cntn1–Notch1（0.51 / 0.05，8） |
| | CLDN | Cldn11–Cldn11（0.07 / 0.07，5） |

## 3. CellChat 与 CCCvelo 的类型标注

匹配方式：CellChat 小鼠符号转大写后与 CCCvelo 匹配（只比较单基因配对）。这是符号匹配，不是同源映射。
同时检查两种方向：CCCvelo 可能把同一对写成"受体→配体"（例如 CCCvelo 的 NLGN1→NRXN1 对应 CellChat 的
Nrxn1→Nlgn1）。同基因配对（如 Negr1–Negr1）只计一次。

| 关系 | 条数 |
|---|---|
| 一致（Secreted↔diffusion、Contact↔contact） | 593（同方向 578，反方向 15） |
| **不一致：CellChat Contact，CCCvelo diffusion** | 180（同方向 134，反方向 46） |
| **不一致：CellChat Secreted，CCCvelo contact** | 2（Tnfsf8–Tnfrsf8、Tnfsf4–Tnfrsf4；都不通过 (b)） |
| CellChat 是 ECM-Receptor / Non-protein，CCCvelo 是 diffusion | 161 |

不一致共 182 条，对应 **148 个不同的 CellChat 配对**，其中 52 个通过 (b)，11 个通过 (c)。
CCCvelo 的 contact 表只有 109 条，大部分膜蛋白黏附配对在 CCCvelo 中被放进了 diffusion 表。

**CCCvelo 自身的方向不一致：** 有 17 对基因，CCCvelo 两个方向都收录，但类型不同（例如 CD80→CD28 是 diffusion，
CD28→CD80 是 contact）：Cd46–Jag1、Cd80–Cd28、Cd80–Ctla4、Cd86–Cd28、Cd86–Ctla4、Efnb1–Epha4、Icam1–Spn、
Nectin2–Tigit、Nectin3–Nectin2、Ptprc–Cd22、Ptprc–Mrc1、Sele–Cd44、Selplg–Sele、Selplg–Sell、Selplg–Selp、
Sema4d–Cd72、Tnfsf4–Tnfrsf4。

**通过 (b) 的 52 个不一致配对**（全部是 CellChat Contact / CCCvelo diffusion；\* 表示也通过 (c)）

| 通路 | 数量 | 配对（CellChat 方向） |
|---|---|---|
| ADGRB | 1 | Rtn4r–Adgrb1 |
| CADM | 5 | Cadm1–Cadm1、Cadm1–Nectin3\*、Cadm3–Cadm1、Cadm3–Cadm3、Cadm3–Nectin3\* |
| CDH | 2 | Cdh2–Cdh2、Cdh4–Cdh4 |
| CLDN | 1 | Cldn11–Cldn11\* |
| CNTN | 5 | Cntn1–Notch1\*、Cntn1–Notch2、Cntn1–Nrcam、Cntn2–Cntn2、Cntn2–L1cam |
| EPHB | 6 | Efnb2–Epha4、Efnb2–Ephb1、Efnb2–Ephb2、Efnb3–Epha4、Efnb3–Ephb1、Efnb3–Ephb2 |
| JAM | 1 | Jam3–Jam3\* |
| L1CAM | 1 | L1cam–L1cam |
| MAG | 1 | Mag–Mag |
| NCAM | 5 | Ncam1–Fgfr1\*、Ncam1–L1cam、Ncam1–Ncam1\*、Ncam1–Ncam2\*、Ncam2–L1cam |
| NECTIN | 2 | Nectin1–Cadm3、Nectin1–Nectin3\* |
| NEGR | 1 | Negr1–Negr1 |
| NGL | 3 | Lrrc4–Ntng2、Lrrc4b–Ptprf\*、Lrrc4c–Ntng1 |
| NRXN | 9 | Nrxn1/2/3 – Nlgn1/2/3（全部 9 种组合） |
| PTPR | 2 | Lrrc4b–Ptprd、Lrrc4b–Ptprs |
| PTPRM | 1 | Ptprm–Ptprm\* |
| SEMA4 | 1 | Sema4d–Plxnb1 |
| SEMA5 | 1 | Sema5b–Plxna1 |
| SEMA6 | 4 | Sema6a–Plxna2、Sema6a–Plxna4、Sema6b–Plxna2、Sema6b–Plxna4 |

全部 148 个配对（包括未通过 (b) 的）见 `cellchat_vs_cccvelo.csv`（`relation = conflict`）。

**当前第一阶段面板在 CellChat 中的情况**

| 面板配对 | CCCvelo | CellChat |
|---|---|---|
| IL1RAPL1→PTPRD | diffusion | 不在库中 |
| FGF14→FGFR2 | diffusion | 不在库中 |
| NEGR1→NEGR1 | diffusion | Cell-Cell Contact（NEGR） |
| NLGN1→NRXN1 | diffusion | 方向相反（Nrxn1→Nlgn1），Cell-Cell Contact（NRXN） |
| EFNA5→EPHA4 / EPHA7 / EPHA5 | contact | Cell-Cell Contact（EPHA） |
| SEMA5A→PLXNB3 | contact | Cell-Cell Contact（SEMA5） |

当前 4 个 "diffusion" 配对，在 CellChat 中要么不存在，要么是 Cell-Cell Contact。

## 4. 限制

- 检测比例按整张组织计算，局部区域里的比例会不同。
- (b) 只看检测比例，没有用平均表达量。
- 多亚基复合物（表 1）暂未处理；ECM-Receptor 和 Non-protein Signaling 不在两类之内。
- 与 CCCvelo 的比较是大写符号匹配，不是经过验证的人鼠同源映射。

## 5. 文件

| 文件 | 内容 |
|---|---|
| `steps.csv` | 表 2 的完整数据（含 c 的敏感性行和"≥2 个配对的配体/受体数"） |
| `complexes.csv` | 表 1 |
| `pairs.csv` | 每个 Secreted/Contact 配对：各步是否通过、检测比例、是否在 h5ad 中、TF 关联、靶基因数 |
| `receptor_targets.csv` | 通过 (b) 的每个受体：使用的 TF、候选靶基因、保留的靶基因 |
| `cellchat_vs_cccvelo.csv` | 所有能匹配到 CCCvelo 的单基因配对：方向、关系、是否通过 (b)/(c) |
| `current_panel_in_cellchat.csv` | 当前面板 8 对在 CellChat 中的方向与类别 |
| `inventory.json` | 来源、哈希、规则 |

```bash
python cellchat_inventory.py fetch
python cellchat_inventory.py report --panel stage1_input_01/lr.csv --output reports/cellchat_inventory_01
```
