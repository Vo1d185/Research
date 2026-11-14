# COMMOT 复现代码核对（Zenodo 10.5281/zenodo.7272562）

2026-11-14。只读代码和数据文件，没有运行其中任何代码。

## 来源

记录：`https://zenodo.org/api/records/7272562`，版本 0.0.2，发布日期 2022-11-02，标题
"COMMOT: Screening cell-cell communication in spatial transcriptomics via collective optimal transport"。
下载了 15 个文件中的 3 个，MD5 与 Zenodo 记录一致：

| 文件 | 大小（字节） | MD5（Zenodo） | SHA256（本地计算） |
|---|---|---|---|
| `downstream_gene_correlation.zip` | 101,570,759 | `2945fc722a492ac048d80ac032ae399a` | `5d507c5e1cc7a9237955fc6ddaa443f5f0a3574f8d53bcc572c67d916b07765c` |
| `pde_model_comparison.zip` | 741,268,637 | `945636f2cf17f753ca61d1bcc333088a` | `ec8eb74c878f299b0b7c017a1ca758777c188c2ea1c672120e62b2de42ccc50f` |
| `visium-mouse-brain.zip` | 320,249,561 | `17b161cfa13f3e5af189eb32d05adaec` | `52cb41cb889ce7c72a617394505b20b7b86a6e02922517b783bb388ef488d363` |

存放在 `data/zenodo_7272562/`（不进版本库）。只解压了代码和文本文件（`extracted_code/`）；
PDE 参数文件 `meta_data.pkl` 用只允许 numpy 数组重建的受限反序列化读取，没有执行其他对象。

**图号：** 代码和笔记本中没有任何图号。下面 (a) 对应的是 `downstream_gene_correlation/` 中计算
"接收信号与已知下游基因相关性"的笔记本，(b) 对应 `pde_model_comparison/`；它们与论文补充图 21–24、
Extended Data Fig. 1 的对应关系是按内容判断的，**没有在存档中得到确认**。

## (a) LR 对的"已知靶基因"

代码：`cell-level_downstream_gene_correlation.ipynb`（spot 层面）和
`cluster-level_downstream_gene_correlation_comparison.ipynb`（聚类层面，另有与 CellChat、Giotto、CellPhoneDB 的比较）。
两者用同一个函数 `get_downstream_genes`，调用参数 `tf_ppr_thresh = 1e-5`（函数默认值 1e-3）。

**规则（小鼠）：**

1. **受体 → TF：** 只用 scSeqComm 的 `TF_PPR_REACTOME_mouse`（不用 KEGG）。保留 `tf_PPR ≥ 1e-5` 的全部行，
   受体的 TF = 这些行中该受体的全部 TF（**不取前 k 个**；同一 TF 在多个通路中出现会重复列出，但靶基因按集合使用，重复没有影响）。
2. **复合物受体：** 把 COMMOT 的受体名按 `_` 拆开、用 `,` 连接，与 scSeqComm 的受体名**整体精确匹配**
   （scSeqComm 有 113 个这样的复合物条目，如 `Acvr1,Acvr2a`）；不按亚基分别查找。匹配不上的配对被跳过。
3. **TF → 靶基因：** `TF_TG_TRRUSTv2_RegNetwork_High_mouse`；只用有靶基因表的 TF；靶基因 = 这些 TF 的靶基因并集。
   受体没有任何可用 TF 时跳过该配对。
4. **没有**：靶基因检测比例阈值、最少靶基因数、去掉 LR 基因、受体表达控制。
5. **靶基因表达：** 全转录组表达先 `filter_genes(min_cells = 3)`，在 log1p(原始计数) 上用
   `highly_variable_genes(min_mean = 0.0125, max_mean = 3, min_disp = 0.5)` 取高变基因，靶基因只在高变基因中计。
   - spot 层面：靶基因得分 = 高变基因中靶基因表达的平均值。Visium 鼠脑一节**没有**对高变基因矩阵做 log1p
     （乳腺癌一节做了），即用原始计数的平均。
   - 聚类层面：每个 leiden 聚类中，靶基因里 Wilcoxon 检验 p ≤ 0.05 且得分 > 0 的比例。
6. **指标：** 每个 LR 对一个 Spearman 相关（spot 层面：所有 spot 上接收信号对靶基因得分；聚类层面：各聚类的平均接收信号对上述比例），
   NaN 记为 0；不做偏相关。
7. **接收信号来源：** Visium 鼠脑用 `dis_thr = 1000`（像素坐标，约 726 µm）、CellChat 小鼠库（`ligand_receptor_database` 默认）、
   `filter_lr_database(min_cell_pct = 0.05)`、`heteromeric = True`、默认 COT 参数（四方案混合）的结果
   （`visium-mouse-brain/1-lr_signaling.ipynb` → `adata_filteredcellchat_dis1000.h5ad`）；聚类层面另外只保留每个基因在 ≥10% spot 中表达的配对。

**数据一致性：** 存档中的 `TF_PPR_REACTOME_mouse.csv`（5,324 行）与本研究从 scSeqComm `95c9bb2` 转换的 REACTOME 行完全相同
（PPR 最大差 5e-16）；TF 列表（857 个）与本研究 `tf_tg.csv` 的 TF 集合相同。

**与本研究规则（计划第 4 节、`evaluate_targets.py`）的差别：**

| | COMMOT | 本研究 |
|---|---|---|
| PPR 来源 | 只用 REACTOME | KEGG + REACTOME，TF 权重取各通路最大值 |
| TF 选择 | PPR ≥ 1e-5 的全部 TF | 权重最高的前 10 个（5、20 为敏感性分析） |
| 复合物受体 | 与 scSeqComm 复合物条目整体匹配 | 按亚基查找，TF 取各亚基中的最高权重（"任一亚基"） |
| 靶基因筛选 | 只限高变基因；无检测阈值、无最少个数；不去 LR 基因 | 检测比例 ≥5%；≥5 个；去掉全部 LR 基因（含亚基） |
| 表达尺度 | 鼠脑：高变基因原始计数 | normalize_total + log1p |
| 指标 | 每个 LR 对的 Spearman | Spearman、控制受体表达的偏 Spearman；LR 对和受体两个层面 |

结论：(a) 与本研究的规则不同。按要求不改主分析，把 COMMOT 的规则作为额外的敏感性分析写进计划（第 8.8 节），在靶基因评价之前提交。

另外记录：论文 Visium 鼠脑分析用的是 `dis_thr = 1000` 像素（约 726 µm），不是教程的 500 像素。本研究的截断（计划第 8.3 节）没有改动。

## (b) Extended Data Fig. 1 的 PDE 模拟

**代码：** 存档中**没有** PDE 求解代码。`benchmark_cot.py` 等从模块 `synthetic_ligand_receptor_interaction`
导入 `plot_fields`，这个模块不在存档中。存档里有每次模拟的结果和输入：`result.h5`（py-pde 输出）、`grid.npy`、
`ligand_prod_*.npy`（配体产生场）、`receptor_init_*.npy`（受体初值）、`meta_data.pkl`（参数）、`result_labels.txt`。

**PDE 设置（`result.h5` 属性）：** py-pde 0.13.3，PDE 类 `LigandReceptorsPDE`，`ScipySolver`，初始步长 0.025，
t = 0 → 50，保存 21 个时间点；`CartesianGrid` 50 × 50，范围 [−200, 200]²，非周期边界（网格中心 −196, −188, …, 196）。

**参数（`meta_data.pkl`，10 种情形 × 10 次模拟，同一情形内 10 次模拟参数相同，只有随机的产生场和受体初值不同）：**

| 参数 | 取值（全部情形） |
|---|---|
| 配体扩散系数 | 200 |
| 结合速率 k_on | 1（有相互作用的对），否则 0 |
| 解离速率 k_off | 0.001（有相互作用的对） |
| 配体降解 | 0.2 |
| 受体降解、复合物降解 | 0 |

| 情形 | 配体 × 受体 | 相互作用对数 |
|---|---|---|
| 1 | 2 × 1 | 2 |
| 2 | 2 × 3 | 4 |
| 3 | 3 × 3 | 9（全连接） |
| 4 | 5 × 2 | 10（全连接） |
| 5 | 3 × 5 | 15（全连接） |
| 6 | 5 × 5 | 10（环状，每个配体 2 个受体） |
| 7 | 5 × 5 | 15（每个配体 3 个受体） |
| 8 | 5 × 5 | 20（每个配体 4 个受体） |
| 9 | 10 × 10 | 20（环状） |
| 10 | 10 × 10 | 100（全连接） |

**OT 一侧（`benchmark_cot.py`）：** `spatial_communication`，`dis_thr = 150`，`cot_eps_p = 0.25`，`cot_rho = 1.0`，
`cost_type = 'euc_square'`（平方欧氏距离），collective 用 weights (1, 0, 0, 0)，"pairwise OT"用 (0, 0, 0, 1)；
另有 UOT、POT 变体（需要存档中的 `commot_v1_ot_compare` 包）。

**评价（`evaluate_result_*.py`）：** 每个相互作用对，在受体初值 > 0 的网格点上，比较 PDE 最终时刻的复合物场 `Complex_i_j`
与 OT 的接收信号场：Spearman、Pearson（NaN 记 0）、均方误差、R²、解释方差。

**说明：** 网格坐标单位在存档中没有写明。参数可以从存档完整读出，但没有求解器代码，不能直接重跑 PDE；
如需重跑，要按 py-pde 0.13.3 和上面的参数重写反应-扩散方程，这部分属于重新实现，不是复现。
