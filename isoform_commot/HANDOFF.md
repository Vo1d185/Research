# 交接说明（给在 VS Code 中继续工作的 Claude 和本人）

本文件总结此前网页版会话中的决定，在新的会话中先读它。

## 研究目标

在 COMMOT 的 collective OT 框架中，把通信类型、空间距离和 PLM-interact 预测的蛋白相互作用分数结合起来，设计新的 cost：
`C = -log F(d) - log max(S, 1e-8)`。研究这些因素如何改变细胞通信的推断结果，再扩展到蛋白 isoform。
T 是模型推断的通信量，不是测得的结合频率；S 是模型预测分数，不是亲和力。

## 已确定的事项

| 事项 | 决定 |
|---|---|
| 工作设置 | r0 = 2h = 120（坐标单位尚未核实；若为 Stereo-seq DNB 单位，约 60 µm），ε=0.2，ρ=1。这是比较用的基线，不是生物学最优值 |
| 固定输入 | `stage1_input_01`：128 个 spot、8 个 LR 对（每类 4 对），PPI 实验期间不改动 |
| 主要评价 | 靶基因相关性（COMMOT 的方法）：Spearman、受体表达对照、控制受体表达后的偏 Spearman |
| 靶基因库 | scSeqComm 小鼠库（COMMOT 用的同一个库），受体→TF 取 PPR 最高的前 10 个 TF（敏感性分析用前 5 和前 20） |
| COMMOT 基线 | 原版 `cot_combine_sparse`，ε=0.1、ρ=10，dis_thr 取 5h（另有 3h、10h）；纯 collective 方案作为另一个对照 |
| PLM 序列 | 小鼠 UniProt canonical 全长序列，两个方向取平均，禁止截断；主模型 650M-humanV11 |
| 对照 | S=1、同类型均值常数、同类型内全部排列 |

## 下一步

1. 按 `PPI_EVAL.md` 在服务器上跑完第 1–7 步（第一次在真实数据上运行）。
2. 下载 COMMOT 的 Zenodo 复现代码（10.5281/zenodo.7272562），核对它从 LR 对到靶基因的规则。
3. 解压 CCCvelo 的 `Database.rar`，把 RecTFDB 和 TFTGDB 转成 `receptor_tf.csv` 和 `tf_tg.csv`，作为第二套靶基因库（没有权重时会自动使用全部 TF）。
4. 扩大规模：用 COMMOT 的 `unot_sinkhorn_l1_sparse` 做稀疏求解，扩大 ROI，增加 LR 对，在独立 ROI 上重复。
5. PDE 模拟（参考 COMMOT Extended Data Fig. 1）：扫描 λ=√(D/c) 来确定 r0；设定不同的结合速率 a_ij，验证 S 的作用机制。
6. isoform：先用现有数据做情景分析（比例 π 从 0 扫到 1，保持总表达量不变），再考虑空间长读长数据（Spl-ISO-Seq2、Visium HD 长读长、SiT）。

## 文件对照

| 文件 | 作用 |
|---|---|
| `stage1.py` | 下载 LR 表和数据、准备输入、运行 OT；新增 `attach-ppi` 和 `export-genes` |
| `research.py` | cost、求解器、各种对照条件（`permutation_mode`、`constant_control`、`include_inverse`） |
| `proteins.py` | UniProt 映射、生成打分任务、导入分数 |
| `plm_score.py` | 官方 PLM-interact 推理（单独环境） |
| `commot_baseline.py` | COMMOT 原版基线 |
| `evaluate_targets.py` | 下载靶基因库、运行靶基因评价 |
| `analyze.py` | 运输方案分析；新增 `control_matched.csv` 和 `shared_entity_split.csv` |
