# isoform 归一化的可行性：H-A isoform 统计与空间 isoform 公开数据调研

2026-11-17。只用公开数据库和文献，没有运行新的打分，没有下载大文件。

## 0. 背景和范围

导师建议的归一化：同一对 LR 基因内，把每个 isoform 组合 (i, j) 的 PLM 分数除以该基因对全部 isoform 组合分数之和：

  S'(i, j) = S(i, j) / Σ_{i', j'} S(i', j')

- 配体和受体都只有 1 个 isoform 时，只有一个组合，S' ≡ 1。所以只有"至少一方 ≥ 2 个 isoform"的配对，归一化后的值才不等于 1。
- 同一基因对内的 S' 之和为 1，只反映各组合之间的相对份额。不同基因对之间的绝对分数高低，不再保留在 S' 中。

本报告只做两件事：
1. 第 1 节：统计 H-A 中满足"至少一方 ≥ 2 个 isoform"的配对。
2. 第 2 节：调研能提供空间 isoform 表达量的公开数据集。

## 1. H-A 配对的 UniProt isoform 统计

### 1.1 方法

- **配对：** PLM LR 评测的 H-A 集（人类、单基因、非同源、不限长度，1,127 对，818 个蛋白；`reports/plm_lr_benchmark/dataset/H-A/positives.csv`）。
- **每个基因的条目：** 评测用的映射（`dataset/sources/protein_map_human.csv`），每个基因恰好一个 reviewed 条目。
- **isoform 来源：** 条目的 ALTERNATIVE PRODUCTS 注释，UniProt release 2026_04，REST 查询于 2026-11-17 03:41 UTC。没有这项注释的条目计为 1 个 isoform。
- **三种计数**（主计数是第一种，另两种作敏感性对照）：

  | 计数 | 含义 |
  |---|---|
  | `n_with_sequence`（主） | UniProt 中有序列的 isoform：状态为 Displayed（canonical）或 Described。PLM 能打分的就是这些 |
  | `n_listed` | 注释中列出的全部 isoform，也包括 External（序列在另一个条目里）和 Not described（没有序列） |
  | `n_splicing` | 只有注释的事件包含 Alternative splicing 时，才取 `n_with_sequence`，否则计为 1 |

- **代码与输出：** 代码 `isoform_counts.py`（`e6dedcb`），输出 `reports/isoform_counts/`：
  - `genes.csv`：每个基因一行；
  - `pairs_H-A.csv`：每个配对一行；
  - `summary.csv`：汇总；
  - `manifest.json`：运行记录。

  先提交代码，再在该提交上重新生成输出。结果与提交前的试运行逐字节相同。

### 1.2 蛋白层面

H-A 的 818 个蛋白中：
- 424 个有 ALTERNATIVE PRODUCTS 注释；
- **420 个（51.3%）有 ≥ 2 个带序列的 isoform**。

各状态的 isoform 数：Displayed 424，Described 1,045，External 5，Not described 2。

带序列的 isoform 数分布：

| 个数 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 11 | 12 | 16 | 17 | 18 | 19 | 21 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 蛋白数 | 398 | 202 | 91 | 46 | 28 | 19 | 13 | 6 | 2 | 4 | 2 | 1 | 3 | 1 | 1 | 1 |

- isoform 最多的蛋白：FGFR1 21，CD44 19，OPRM1 18，FGFR2 / FN1 / VEGFA 各 17，CD46 16。
- 注释中的产生方式：
  - 411 个条目只有 Alternative splicing；
  - 8 个条目同时含 alternative promoter usage 或 alternative initiation；
  - 5 个条目只有 alternative promoter usage 或 alternative initiation。

### 1.3 配对层面（主计数 `n_with_sequence`）

| 子集 | 配对 | **至少一方 ≥ 2 个 isoform** | 比例 | 配体 ≥ 2 | 受体 ≥ 2 | 两方都 ≥ 2 | 组合数中位数 | 组合数最大 | 组合数合计 |
|---|---|---|---|---|---|---|---|---|---|
| **H-A 全部** | 1,127 | **862** | **76.5%** | 515 | 685 | 338 | 3 | 323 | 8,194 |
| 两个蛋白都 ≤ 800 aa | 605 | 390 | 64.5% | 235 | 284 | 129 | 2 | 132 | 2,846 |
| 至少一个 > 800 aa | 522 | 472 | 90.4% | 280 | 401 | 209 | 4 | 323 | 5,348 |
| Secreted | 619 | 425 | 68.7% | 227 | 336 | 138 | 2 | 136 | 3,909 |
| Contact | 349 | 303 | 86.8% | 204 | 241 | 142 | 4 | 132 | 2,405 |
| ECM | 159 | 134 | 84.3% | 84 | 108 | 58 | 4 | 323 | 1,880 |

敏感性对照（H-A 全部）：
- `n_listed`：864 对（76.7%），组合合计 8,226；
- `n_splicing`：859 对（76.2%），组合合计 8,007。

三种计数下的配对数在 859–864 之间。

其他数量：
- 组合数恰好为 2 的配对有 266 个，≥ 10 的有 195 个。
- 组合最多的五对：FN1–CD44 323（17 × 19，ECM）、VEGFA–FLT1 136、CD209–CEACAM1 132、NCAM1–FGFR1 126、LGALS9–CD44 114。
- 8,194 个组合中，65.3% 来自"至少一个 > 800 aa"的配对。
  在 humanV11 评测中，这一组的 AUROC 明显较低（随机负样本 0.660，保持度数负样本 0.611；见 `reports/plm_lr_benchmark/report.md` 第 1.1 节）。
- **规模（只是估计，没有打分）：** 如果给 H-A 全部组合打分，两个输入顺序共约 16,400 个作业，约是 humanV11 本次 human_full 作业量（220,376）的 7%。
  isoform 序列目前没有下载。

### 1.4 注意事项

1. UniProt 只收录经过整理的 isoform，比 GENCODE/Ensembl 的转录本少得多。这里的数是"有蛋白序列、PLM 能直接打分"的 isoform 数，不是转录本多样性的全貌。
2. 有多个 reviewed 条目的基因按评测的映射规则已经去掉，不在 H-A 中：NRXN1/2/3、CALCA、PRNP，相关配对见 `dataset/sources/unresolved_human.csv`。
   这些基因的 isoform 分在不同条目里（例如 NRXN 的 α/β），本统计没有计入。
3. 计数包括 canonical isoform。"≥ 2"指 canonical 之外至少还有一个带序列的 isoform。
4. 小鼠情况（只作对照，来自已有的[清单 02](isoform_candidates_02/README.md)）：后部切片 OT 面板 358 个基因中，347 个有 reviewed 条目，其中 120 个（34.6%）有 ≥ 2 个 isoform 标识，97 个是清单的候选基因。
   清单 02 的计数是全部列出的 isoform 标识，与本节的 `n_with_sequence` 不完全相同。

## 2. 能提供空间 isoform 表达量的公开数据集（只调研）

我们现有的数据是 10x Visium 小鼠脑矢状切片（`V1_Mouse_Brain_Sagittal_Posterior/Anterior`，3′ 短读长）。它只有基因层面的计数（见 `reports/isoform_plan_draft.md`），不能区分编码区不同的 isoform。

要在空间上区分 isoform，需要对空间条形码 cDNA 做全长（长读长）测序。下表只列有长读长空间数据、且能确认获取方式的数据集。

**在本次调研中逐项核对的内容：**
- 文章的 Data availability 段落：Europe PMC 全文 XML，存在 `data/literature/isoform_datasets_epmc/`，不提交。
- GEO 条目的公开状态、样本名和文件大小：GEO 页面。
- 数据页面本身。

没有下载任何数据文件。

### 2.1 数据集

| # | 数据集 | 物种 / 组织 | 平台和空间分辨率 | 能拿到的 isoform 数据 | 获取方式 |
|---|---|---|---|---|---|
| 1 | **SiT**，Lebrigand et al. 2023，*Nucleic Acids Res* 51(8):e47，doi:10.1093/nar/gkad169 | 小鼠，C57BL/6J，8–12 周：左半脑冠状切片 2 张（CBS1、CBS2，相距 50 µm），嗅球 1 张（MOB） | Visium（55 µm spot）+ ONT PromethION，有 UMI | 每个样本一个 Seurat 对象，含 `ISO` assay（isoform 层面计数）、长读长基因层面 assay 和短读长 assay | GEO **GSE153859**，公开（2020-08-26）。长读长处理文件 `GSM4656182/3/4_*-SiT.tar.gz` 分别为 54 / 51 / 16 MB；短读长 GSM4656179–81；原始数据在 SRA。ISS 数据和代码在 Zenodo 10.5281/zenodo.7664367；在线浏览 www.isomics.eu |
| 2 | Joglekar et al. 2021，*Nat Commun* 12:463，doi:10.1038/s41467-020-20343-5 | 小鼠 P8，矢状切片 1 张（n = 1）；同时有 P7 海马和前额叶的单细胞长读长数据 | Visium（55 µm）+ ONT，以及 Visium + PacBio（文中称 sliso-Seq） | 只有原始读长：GEO 中 Visium_PacBio（GSM4800809）和 Visium_ONT（GSM4800810）没有处理后文件。文中用外显子层面的 Ψ 分析，SiT 文章评价其覆盖低、没有 UMI | GEO **GSE158450**，公开（2020-11-17），原始数据在 SRA；短读长 Visium 矩阵为 GSM4800808。STOmicsDB STDS0000004 有该空间数据的镜像（3,024 spots）；单细胞 isoform 图谱见 isoformAtlas.com |
| 3 | ONT 公开数据 "Visium HD 3′ mouse brain"（EPI2ME，2025-06） | 小鼠，8 周雄性脑，新鲜冰冻切片 1 张 | Visium HD 3′（2 µm 条形码，常用 8 µm bin）+ ONT PromethION（4 个 flow cell） | 有 EPI2ME `wf-single-cell` 的输出（含全长转录本信息），是否含转录本计数矩阵未核实；另有 Space Ranger 输出 | AWS S3 `s3://ont-open-data/visium_hd_2025.06`，匿名下载，不需要登录（`aws s3 sync --no-sign-request …`）。**体量很大**：POD5 5.0 TB、BAM 707 GB、预处理 1.5 TB、Space Ranger 219 GB，`wf-single-cell` 大小未列出 |
| 4 | PacBio 公开数据 "Kinnex Visium HD mouse brain" | 小鼠 P56 雄性脑（按 Spl-ISO-Seq2 文章的描述） | Visium HD 3′ + PacBio Revio（Kinnex） | 目录名为 `1-Sreads`（分段后的读长）；是否有处理后的 isoform 矩阵未核实 | `https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-RevioSPRQ-Kinnex-VisiumHD-mouseBrain/`（路径取自 Spl-ISO-Seq2 文章 Methods；目录内容和大小未核实） |
| 5 | **Spl-ISO-Seq2**，Tilgner 实验室，bioRxiv 10.1101/2025.06.25.661563；正式版在 *Nature Methods*（2026） | 小鼠 P56 雄性脑冠状切片：Stereo-seq 2 张（S1、S2，相邻），Visium HD 1 张（V1）；另外重新分析了上面的第 3、4 个数据集（文中 V2、V3） | Stereo-seq（约 220 nm spot，500 nm 分辨率，按细胞分割）+ ONT / PacBio；Visium HD 3′ + ONT / PacBio（按细胞分割或 8 µm bin） | 7 个长读长空间数据集（5 张切片）的处理后结果，可以在线浏览 | bioRxiv 版的 Data availability 只写了在线浏览 spatialisoforms.weill.cornell.edu，**没有找到原始或处理后数据的登录号**；Nature Methods 版未核实 |
| 6 | **scNaST**，Boileau et al. 2022，*Front Genet* 13:912572，doi:10.3389/fgene.2022.912572 | 小鼠心脏，心梗后，短轴切片 4 张 | Visium（55 µm）+ ONT | 7,616 个 spot、19,794 个转录本 isoform（de novo 转录组） | SRA BioProject **PRJNA843953**（8 个 run），补充数据在 Zenodo 10.5281/zenodo.6546361，代码 github.com/dieterich-lab/ScNaST |
| 7 | 人类胶质瘤，Ren et al. 2023，*Nat Commun* 14，doi:10.1038/s41467-023-36707-6 | 人类，11 个样本：GEO 描述为 5 DMG、5 GBM、1 个瘤周；SPLISOSM 的脚本写作 6 GBM + 5 DMG | Visium（55 µm）+ ONT | 处理后数据在 GEO 和 Figshare，是否含 isoform 矩阵未核实；SPLISOSM 使用的是该研究提供的处理后 isoform 定量 | 处理后数据：GEO **GSE194329**（公开）和 Figshare 10.6084/m9.figshare.20653908。**长读长原始数据受控**：GSA-Human HRA001960，需要申请，文章称 3–6 周 |
| 8 | **SPLISOSM**，Su et al. 2026，*Nat Biotechnol*，doi:10.1038/s41587-025-02965-6（方法文章，重新分析已有数据） | 小鼠脑（SiT、Visium HD + ONT）；人类胶质瘤（第 7 项）；另有短读长 Visium 和 Xenium Prime 5K | 同原始数据 | 整理好的逐样本 isoform 层面 AnnData，以及空间变异检验结果 | Zenodo 10.5281/zenodo.16905935（处理后数据，检验结果 zip 约 60 MB）和 10.5281/zenodo.15556390（`data/` 目录）；文档中的 SiT-MOB 示例约 100 MB；代码 github.com/JiayuSuPKU/SPLISOSM_paper |

### 2.2 与我们数据的对应关系（只列事实）

- **物种和组织：**
  - 与我们的小鼠成年脑 Visium 最接近的是 SiT（成年小鼠脑，Visium，55 µm spot，有 isoform 层面计数，处理后文件共约 120 MB）。但 SiT 是冠状切片，我们的是矢状切片。
  - Joglekar 是矢状切片，但是 P8 幼鼠，而且只有原始读长。
  - 两套 Visium HD（ONT、PacBio）和 Spl-ISO-Seq2 都是成年小鼠冠状切片，分辨率更高（2 µm / 8 µm bin，Stereo-seq 500 nm），数据量大，或只有在线浏览。
- **没有任何公开数据集与我们的两张切片（`V1_Mouse_Brain_Sagittal_Posterior/Anterior`）是同一张组织。**
  要使用空间 isoform 比例，只有两种做法：一是在不同切片之间按脑区对应；二是换用这些数据集本身做通讯推断。
- **人类：** 唯一找到的人类长读长空间数据是胶质瘤（第 7 项），原始数据需要申请。PLM 评测用的是人类序列，而 OT 分析用的是小鼠切片，这一点不变。
- **短读长 3′ 数据**（包括我们自己的切片）最多只能给出 3′ 端多样性，例如 SPLISOSM 用 Sierra 做的 TREND 定量。它区分不了大多数编码区不同的 isoform，所以没有列入上表。

### 2.3 没有核实的事项

- 第 3、4 项公开目录中是否有现成的转录本/isoform 计数矩阵，以及第 4 项的大小。
- 第 5 项 Nature Methods 正式版的 Data availability。
- 第 7 项 GEO/Figshare 处理后文件的具体内容。
- 各数据集里 isoform 的注释体系（GENCODE/Ensembl 转录本、de novo 转录本）与 UniProt isoform 的对应关系。
  这决定了第 1 节的 UniProt isoform 能否直接对上空间数据中的转录本，需要拿到文件后才能核对。

## 3. 小结（只陈述观察到的结果）

1. H-A 的 1,127 对中，**862 对（76.5%）**至少一方有 ≥ 2 个带序列的 UniProt isoform，归一化后的值会不等于 1；另外 265 对（23.5%）恒等于 1。
   - 这个比例在含 > 800 aa 蛋白的配对中是 90.4%，在两个都 ≤ 800 aa 的配对中是 64.5%。
   - 换用另外两种计数，比例为 76.2%–76.7%。
2. 能提供空间 isoform 表达量的公开数据集主要是**小鼠脑**：
   - SiT（Visium + ONT）：成年、冠状，处理后文件小、公开；
   - Joglekar：P8、矢状，只有原始读长；
   - ONT 和 PacBio 的 Visium HD 3′ 数据：成年，数据量大；
   - Spl-ISO-Seq2：成年，Stereo-seq / Visium HD，只找到在线浏览。

   此外还有小鼠心脏（scNaST）和人类胶质瘤（原始数据受控），没有和我们切片同一张组织的数据。

是否下载其中某个数据集（例如 SiT 的三个处理后文件，约 120 MB），以及下一步做什么，需要你决定。本报告没有做任何下载或分析。
