# 下一阶段：结果解释 → 参数稳健性 → 真实数据

## 1. 先解释现有模拟结果

运行 README 中的 analyze.py，保留三个图和配对比较表。应回答：

- r0 变化是否改变 diffusion 的距离分布？contact 的 T 是否因竞争而变化？
- 固定 family、r0 后，PPI 主要改变总量，还是改变连接上的分配？
- 哪些 LR 对获得/失去通信量，是否符合输入的模拟 score？

这一步用于理解模型机制，不能把模拟 score 的结果称为 PLM-interact 的真实预测效果。
不要用最大 total_T 或最大平均距离选择“最好”的 r0。

接下来固定少量代表性 r0 条件，改变 eps/rho，检查结论是否只在一个正则化/惩罚设定下出现。
可先探索 eps={0.1,0.2,0.5}、rho={0.5,1,2}，每次使用独立输入/输出目录并保留其他条件。
这些只是建议的敏感性范围，不是原论文推荐值。模拟中的小 score 在较小 eps 下可能产生很强的分配偏好。

## 2. 已核对的 CCCvelo 文献依据

论文：Decoding cell state transitions driven by dynamic cell-cell communication in spatial transcriptomics.
DOI：https://doi.org/10.1038/s43588-025-00934-2
Nature Computational Science 6, 193–207 (2026).

本次读取用户上传的扫描 PDF；定位为 PDF 第 10、11、13 页，印刷页码 202、203、205。
公式页已进行图像核对。以下是文献摘记，不是已经实现的适配器。

### Methods：预处理和先验网络（印刷页 202）

- 描述 Seurat 4.2.0 的 NormalizeData、RunPCA、FindNeighbors；缺少细胞注释时使用 FindClusters。
- 在发送/接收细胞群中筛选候选 ligand/receptor；举出的默认平均表达阈值为 0.05，表达细胞比例阈值为 10%。
- 使用自己的 LigRecDB、RecTFDB、TFTGDB。LigRecDB 含 3659 个非冗余 LR 交互。
- 人源先验数据库分析小鼠时进行同源基因映射。
- 本项目只研究 LR cost 与 OT；无需为了跑这一阶段直接复现 TF/TG、PINN 或 latent time 模块。

### Methods：空间 LR 信号（印刷页 203，公式 1–3）

- diffusion：对发送者求和的 (1/d_sr)*L_s^j*R_r^j。
- contact：对接收者的直接邻居求和 L_s^j*R_r^j。
- 论文用 Delaunay triangulation 构建直接邻居。
- 原文把 d_sr 称为 relative distance。具体代码里的坐标缩放、零距离处理还需核对，不能自行假定为微米。
- 原论文这里计算的是 LRS，不是 COMMOT 的 OT cost；本项目的 -log(F*S) 是另行提出的模型。

当前脚本的半径邻接、d_ref 归一化、diffusion 截断半径、S 融合和 OT 都应与论文原方法区分。
类型来源应核对 LigRecDB 的实际字段；当前接受 CellChat 标签只是通用输入便利功能，不证明与 CCCvelo 的分组相同。
目前尚未检验该 prior database 文件的具体类型编码或下载完成真实数据。

### Data availability（印刷页 205）

| 数据 | 论文列出的来源 | 对本研究的初步用途 |
|---|---|---|
| 小鼠成年皮层 Stereo-seq | 文中提供 mousebrain_bin60.h5ad 的 Dropbox 链接 | 优先评估的 gene-level 空间基线 |
| 小鼠胚胎 Slide-seqV2 | GEO GSE197353 | 后续跨组织/发育情景验证 |
| 人前列腺癌 MERFISH | Vizgen MERSCOPE FFPE Human Immuno-oncology datasets | 人源数据对照，需检查靶向基因面板的 LR 覆盖 |

先验数据库及分析归档：https://doi.org/10.5281/zenodo.17384457
作者软件：https://github.com/SunXQlab/CCCvelo
作者分析：https://github.com/SunXQlab/CCCvelo-analysis

作者分析仓库已读版本：7c3c95657cd89c17a548c2bf9659ef1b47c1d754。
已读取该版本 README 和
[apply_in_cortex/1_select_LRTG.R](https://github.com/SunXQlab/CCCvelo-analysis/blob/7c3c95657cd89c17a548c2bf9659ef1b47c1d754/apply_in_cortex/1_select_LRTG.R)。

重要差异：这个 cortex 脚本出现 SCTransform、run_Imputation 调用，
以及 select_LRTG(..., pct.ct=0.01, expr.ct=0.01)。
不能把正文中举出的默认值当作所有数据集实际使用的参数。
脚本中计算了 imputed matrix，并不单凭这一行就能断言它被用作最终 L/R 表达；
还需检查 select_LRTG、表达层与后续 Python 输入链路。
该版本 README 所列 0_proprecess_data.R 的直接 URL 返回 404，需定位实际文件，不能当作已复现。

## 3. 推荐的第一轮真实数据范围

优先核实小鼠皮层 Stereo-seq 文件及其表达层、空间坐标和注释。
先选一个连续组织区域（例如约 100–200 个位置）和少量已确认类型的 LR 对；
保留邻域边界信息并报告裁剪影响，不以随机打散空间点代替区域选择。
当前稠密实现必须满足 n_ligand*n_receptor*n_spots^2 <= 4000000。
该小区域试验不等于原文使用整个发送细胞群的完整复现。

实现真实适配器之前逐项确认：

1. 原始 counts、归一化层和坐标单位；细胞注释及 sender/receiver 范围。
2. LigRecDB 实体、物种映射、类型标签、多亚基复合物策略。
3. 单独建立 Delaunay 邻接对照，与当前半径方案区分；不引入分割流程。
4. 先运行 gene-level 的距离实验，再加入有来源记录的真实 PLM-interact 分数。
5. 若要提出 isoform 特异结论，需要 isoform 表达或明确的比例假设及其敏感性实验。
   这些 CCCvelo 空间数据本身不能被当成已测得的 isoform 丰度。

验收输出应包括输入来源、转换日志、LR 类型/缺失/复合物处理报告、固定预算下的消融实验，
以及参数稳定性；在没有外部生物学验证时不宣称某个 r0 最优或估计了真实结合频率。
