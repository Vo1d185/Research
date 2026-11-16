# 交接说明（给在 VS Code 中继续工作的 Claude 和本人）

本文件总结此前网页版会话中的决定，在新的会话中先读它。

## 当前状态（2026-11-16 08:05 UTC）：先读这一节

分支 `research/isoform-commot-v1`，代码在 `/home/liu/Research/isoform_commot`，主环境 `/home/liu/.conda/envs/commot-research/bin/python`。
PPI 作为 cost 的方向已按预先登记的检查停止（S = 1）；现在的 cost 是 C = −log F(d)。下面的"第一阶段"各节是历史记录。

### 必须遵守的规则

- **先登记，后计算：** 任何新的选择先写进计划并提交，再运行；已经看到相关结果时，停下来问用户。
- **第二阶段两张切片的靶基因评价都已经做完并看过结果**（后部 11-16 03:46，前部 05:20，各一次）。
  之后任何新的分析（换评价方法、换参数、加条件）都属于"看过结果之后"，必须先由用户决定并登记。
- **提交与推送：** 用 `git -c user.name=Claude -c user.email=noreply@anthropic.com commit`；不推送，每次提交后提醒用户推送。
- **运行方式：** 长运行用 `setsid nohup` 与会话分离，manifest 记录启动时的 git HEAD 和工作区状态；不删除已有输出（需要时改名移走）。
- **PDE 队列不能停：** 它在 SCHED_IDLE 下运行，第二阶段的任务优先。

### 计划在哪里

| 内容 | 位置 |
|---|---|
| 第二阶段方法比较、r₀ 选择（8.5）、收敛（8.7）、前部确认（8.9）、配对比较（8.10）、稳健性（8.11）、C726（8.12） | `reports/stage2_plan.md` 第 8 节 |
| **r₀ 验证**（条件 F、I、P，判断标准 (a)–(d)，阴性对照的解读） | `reports/stage2_plan.md` 第 9 节（`45ee06a`，运行记录 9.8） |
| PDE 模拟中 r₀ 与 λ 的关系 | `reports/pde_r0_plan.md`（运行记录第 14 节） |
| isoform 方案草稿（等与导师讨论，不执行） | `reports/isoform_plan_draft.md`，盘点 `reports/isoform_candidates_02/` |
| 第 8 节报告框架（结果部分只有标题） | `reports/stage2_section8/report.md` |

### 已经完成并提交的

- 收敛判断（`7ff2bc7`，计划第 9.9 节）：两张切片全部条件的硬性要求通过；稳定性诊断 A（r₀ = 50、400）、B、F 都满足，50000 次的结果全部接受。
- 后部评价和 r₀ 选择（`d1aac56`，计划第 9.10 节）：**r₀\* = 50 µm**（网格最小值，四个 r₀ 之间差别很小）。
- 前部确认性检验（`fd75068`，计划第 9.11 节）：**第 9.5 节判为"部分满足"**，(a)、(b) 不满足，(c)、(d) 只按字面规则满足；
  第 8.9 节的严格确认和容差确认都不通过。
- 报告：`reports/stage2_section8/r0_validation.md`（r₀ 验证，完整）和 `reports/stage2_section8/report.md`（第 8 节，除第 11 节外已填）。
- 结果目录：`reports/stage2_{post,ant}_01/`（第 8 节）、`reports/stage2_r0val_{post,ant}/`（第 9 节）、`reports/stage2_stability_01/`、
  `reports/stage2_convergence_{post,ant}_01/`；评价输出在 `target_eval_stage2_{post,ant}/`（不入库）。

### 还在运行的

| 运行 | 状态 |
|---|---|
| PDE 队列（`pde_r0_results/`，`simulation.jobqueue`） | 11-16 07:58 起并发上限改为 28（`pde_r0_results/limit`）；4404 个作业完成约 2360 个，剩下的都是 ε = 0.1、ρ = 10 |

第二阶段没有还在运行的任务。稳健性（8.11）已完成：`reports/stage2_robustness_01/`，报告第 11 节。

### 接下来的顺序（不用等用户）

1. **PDE：** 全部作业结束、稳定性诊断（D1–D3）和第 5.1 节的检查通过后，`simulation.benchmark evaluate` → `analyze`，写 `reports/pde_r0_01/`。
2. **等用户决定第二阶段的下一步。** r₀ 验证的结论是前部没有重复出后部的结果；报告列了规则之外的观察
   （坐标打乱后指标仍保留约 80%；ε = 0.2、ρ = 1 下 r₀ 的作用更明显且两张切片一致）。是否据此改评价方法或加条件，由用户决定并先登记。
3. isoform 方案草稿仍在等用户与导师讨论，不执行。

## 研究目标（第一阶段起的历史记录，以下各节同）

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
