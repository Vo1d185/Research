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

## 真实数据上的发现（2026-11-14，蛋白映射）

- **NRXN1**：有两个已审核的小鼠条目，α 型 Q9CS84（1514 aa）和 β 型 P0DI97（468 aa）。面板中的配对是 `NLGN1 → NRXN1`，neuroligin 结合两者共有的 LNS6 结构域，所以主分析用 β 型（`override.csv`）。α 型保留下来，作为以后 isoform 情景分析的第一个真实例子。gene-level 的 Nrxn1 表达同时包含 α 和 β。
- **NLGN1 → NRXN1 在 CCCvelo 表中被标为 diffusion。** 两者都是膜蛋白，靠跨突触接触结合，生物学上更接近 contact。第一阶段固定沿用 CCCvelo 的分类，不改；第二阶段和 CellChatDB 的注释逐对比较。
- **FGF14**：属于 FGF11–14（FHF），细胞内蛋白，一般认为不分泌、不激活 FGFR。如果面板中有 `FGF14 → FGFR2`，其数据库注释值得怀疑，可以用来检验 PLM 分数会不会给它打低分。
- **超长蛋白**：PTPRD 1912 aa、PLXNB3 1902 aa、SEMA5A 1077 aa 等。主分析用全长序列加 `--allow-long`；敏感性分析用 `--region extracellular`。
- **面板**：IL1RAPL1→PTPRD、FGF14→FGFR2、NEGR1→NEGR1、NLGN1→NRXN1（CCCvelo 标为 diffusion）；EFNA5→EPHA4/EPHA5/EPHA7、SEMA5A→PLXNB3（contact）。四个 "diffusion" 对中，三个是膜蛋白之间的粘附（更接近 contact），FGF14 是细胞内蛋白。所以第一阶段 r0 的结果只能当作方法本身的行为，不能解读为分泌信号的扩散距离。第二阶段改用 CellChatDB 的 Secreted Signaling 配对。
- **竞争结构**：只有 EFNA5 被三个 contact 对共用，所以 PPI 带来的预算重新分配预计主要出现在这一组；四个 diffusion 对互不共享实体，PPI 对它们基本只是整体缩放。
- **长度**：全长时，IL1RAPL1–PTPRD 和 SEMA5A–PLXNB3 超过 1603 token；用胞外段时，SEMA5A–PLXNB3 仍然超过。所以两次打分都要加 `--allow-long`，并在结果里标明。

## PLM-interact 打分结果（2026-11-14，650M-humanV11）

- 环境：plminteract（Python 3.10，torch 2.4.1+cu121，transformers 4.40.1，RTX 6000 Ada）。checkpoint 提交 e86e392dec13dd0c23252c94947b04a7a9821b0e。
- 自检：官方小鼠测试集 n=2000（179 个正样本），AUPR 0.913，AUROC 0.980，安装和权重正常。
- **8 个已知配对几乎全部接近 0**（全长：EFNA5–EPHA 约 1e-4，NLGN1–NRXN1β 2e-4，IL1RAPL1–PTPRD 0.002，NEGR1–NEGR1 0.024）。只有 SEMA5A–PLXNB3 在胞外段时为 0.34，全长时为 4e-5，对截取区域极其敏感。换成人源序列，结果几乎一样，所以不是物种问题。
- 已知配对和未注释组合的区分度：AUROC 0.68（全长）、0.62（胞外段），区分能力很弱。
- **结论：这批分数不进 OT。** 按 −log S 计算，S=1e-4 时每条连接要额外加约 9 个 cost 单位，而放弃运输只需要约 2（2ρ），几乎所有连接都会被关掉，结果完全由分数决定，和生物学无关。不允许事后调整 S 的变换方式让结果"看起来合理"，那是循环论证。
- 下一步诊断：(1) 换 650M-humanV12 checkpoint 再打一次分；(2) 阳性对照：经典的胞外配体–受体对，按受体长度分成 ≤800 aa 和 >800 aa 两组，用来区分"长度超出训练范围"和"这一类相互作用本身学不到"这两种原因。

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
