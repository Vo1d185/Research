# PLM-interact 分数诊断（不进入 OT）

目的：在决定是否把 PLM-interact 分数放进 cost 之前，检查它对已知配体–受体对是否给出高分。
这里的结果只用于判断分数质量，**不生成 ppi.csv，不参与 OT**。

## 内容

| 文件 | 作用 |
|---|---|
| `human_pairs.py` | 配对表 → 人源 UniProt（reviewed、主基因名完全一致，否则要求 override）→ 两个输入顺序的 `plm_input.csv`；按受体长度 ≤800 / >800 aa 分组 |
| `summarize.py` | 核对分数文件来自本次输入，输出 `pair_scores.csv`、`group_summary.csv`、`ligand_ranks.csv` |
| `positive_controls.csv` | 15 对经典人源配体–受体阳性对照 |
| `lr_panel_human.csv` | 小鼠面板 8 对，改用人源同源序列（物种检查） |
| `human_overrides.csv` | NRXN1 用 β 型 P58400，与小鼠主分析一致（见 HANDOFF.md） |
| `known_cross_reactivity.csv` | 已知的交叉结合（GH1–PRLR），报告但不计入 AUROC |
| `runs/` | 全部结果，纳入 git |

## 复现（2026-11-14）

```bash
# commot-research 环境
python diagnostics/human_pairs.py --pairs diagnostics/positive_controls.csv \
  --known-cross diagnostics/known_cross_reactivity.csv --output diagnostics/runs/positive_controls_01
python diagnostics/human_pairs.py --pairs diagnostics/lr_panel_human.csv \
  --overrides diagnostics/human_overrides.csv --output diagnostics/runs/lr_panel_human_01
# plminteract 环境，V11 / V12 各一次
python plm_score.py score --jobs diagnostics/runs/<run>/plm_input.csv \
  --checkpoint models/plm_650m_humanV1X/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D \
  --checkpoint-name danliu1226/PLM-interact-650M-humanV1X --checkpoint-revision <sha> \
  --allow-long --output diagnostics/runs/<run>/scores_humanV1X.csv
# commot-research 环境
python diagnostics/summarize.py --run diagnostics/runs/<run>
```

Checkpoint：humanV11@`e86e392dec13dd0c23252c94947b04a7a9821b0e`，
humanV12@`2d0e7fb3deef94efa3f217daccc79092a6a42ceb`；基础模型 esm2_t33_650M_UR50D@`08e4846e537177426273712802403f7ba8261b6c`。
官方源码 commit `ab2ae6ae1aa81accf4c1f7f6f341164baf97809a`。

`runs/mouse_panel_01/` 是主流程小鼠面板分数（`ppi_01*` 目录）的存档副本：
V11/V12 × 全长/胞外段，以及两个 checkpoint 的自检结果。

## 结果摘要

两个 checkpoint 的安装自检都正常（小鼠官方测试集抽 2000 对，阳性 179）：
V11 AUPR 0.913 / AUROC 0.980；V12 AUPR 0.838 / AUROC 0.955。

**小鼠面板（两个顺序取平均）**

| 配对 | V11 全长 | V11 胞外 | V12 全长 | V12 胞外 |
|---|---|---|---|---|
| IL1RAPL1–PTPRD | 0.0023 | 0.0014 | 0.0015 | 0.025 |
| FGF14–FGFR2 | 0.0048 | 0.00003 | 0.99 | 0.0058 |
| NEGR1–NEGR1 | 0.024 | 0.0032 | 1.0 | 1.0 |
| NLGN1–NRXN1 | 0.0002 | 0.00006 | 0.0045 | 0.002 |
| EFNA5–EPHA4 | 0.0001 | 0.0001 | 1.0 | 0.17 |
| EFNA5–EPHA7 | 0.00008 | 0.00005 | 1.0 | 0.079 |
| EFNA5–EPHA5 | 0.00006 | 0.00003 | 0.19 | 0.43 |
| SEMA5A–PLXNB3 | 0.00004 | 0.34 | 0.0005 | 0.68 |

已知对 vs 面板未注释交叉对 AUROC：V11 0.68 / 0.62，V12 0.86 / 0.76（全长 / 胞外）。
人源同源序列的结果与小鼠几乎相同，物种不是主要原因。

**人源阳性对照（全长）**

| 组 | checkpoint | 已知对中位数 | 已知对 ≥0.5 | 交叉对中位数 | 交叉对 ≥0.5 | AUROC |
|---|---|---|---|---|---|---|
| 受体 ≤800 aa（10 对） | V11 | 0.97 | 6/10 | 0.002 | 24% | 0.82 |
| | V12 | 0.999 | 10/10 | 0.42 | 45% | 0.89 |
| 受体 >800 aa（5 对） | V11 | 0.013 | 2/5 | 0.0002 | 19% | 0.77 |
| | V12 | 0.997 | 3/5 | 0.78 | 66% | 0.64 |

- V11 的失败集中在长受体，也包括 IL2–IL2RA（0.0005）、IL6–IL6R（0.014）、EPO–EPOR（0.048）。
- V12 对已知对给高分，但对未注释组合也常给高分：IL6R、EGFR、IFNGR1、INSR 与几乎任意配体的
  中位分数约 0.9；V11 对 PRLR 也有同样现象（0.99）。这类分数更像是"该蛋白容易被判为有相互作用"，
  而不是配对特异的。
- V12 给 FGF14–FGFR2 0.99。FGF14 是胞内 FHF，一般认为不激活 FGFR，这个高分不能当作结合证据。
- 已知对很可能出现在 STRING 训练数据中；高分不一定是独立证据。

结论（暂定）：在当前面板上，两个 checkpoint 的分数对 checkpoint、序列区段都很敏感，
且存在蛋白层面的偏置。分数不能直接作为 `-log S` 进入 cost。

## 限制

- UniProt canonical 全长包含前体（pro-EGF 1207 aa、preproinsulin、pro-TGF-β1），不是成熟配体。
  EGF–EGFR 中配体本身也超过 800 aa。
- 多亚基受体只用单链（如 IL2 高亲和受体需要 IL2RA/IL2RB/IL2RG）。
- 交叉对只是未注释，不是确定的阴性；阳性对照每组只有 5–10 对。

## 预先登记的最后一次尝试（2026-11-14 写定，计算前提交）

本节在计算 Leakage-Free 分数和任何双中心化分数**之前**写定并提交，代码 `preregistered.py` 一同提交。
只做这一次；看到结果后不修改阈值、变换、候选方法或输入。

### 候选方法（4 个）

| 编号 | checkpoint | 变换 |
|---|---|---|
| M1 | 650M-Leakage-Free-Dataset @ `d995319a6e53510bfbb4eff6ab0292fedf2ec4db` | 原始分数 |
| M2 | humanV11 @ `e86e392dec13dd0c23252c94947b04a7a9821b0e` | 双中心化 |
| M3 | humanV12 @ `2d0e7fb3deef94efa3f217daccc79092a6a42ceb` | 双中心化 |
| M4 | Leakage-Free（同上） | 双中心化 |

V11、V12 的原始分数一并输出作参照，不参与判定。

### 定义

- **原始分数**：两个输入顺序的概率取平均（与 `ppi.csv` 一致）。
- **双中心化**：每个顺序的概率先截断到 [1e-6, 1−1e-6] 再取 logit，两个顺序的 logit 取平均得到 x；
  在该集合**完整的配体 × 受体网格**上（人源对照 15×15 = 225 对，小鼠面板 6×8 = 48 对，
  已知对和已知交叉结合都包括在内）计算 z = x − 配体均值 − 受体均值 + 总均值。
- **固定输入**：人源对照 `runs/positive_controls_01`（plm_input.csv sha256 `7ebe478c…`）；
  小鼠面板全长 `plm_jobs_01`（`7eac12f1…`）、胞外段 `plm_jobs_01_ecd`（`402cf988…`），NRXN1 为 β 型。
  V11、V12 使用已有分数，不重新打分；Leakage-Free 对这三组输入各打分一次，设置与之前相同
  （不截断、`--allow-long`、两个输入顺序）。
- **安装检查**：Leakage-Free 先跑同样的小鼠自检（2000 对）。AUPR < 0.5 视为安装问题，修好后重跑，
  不算作一次方法结果。

### 通过标准（三条同时满足）

1. **AUROC ≥ 0.85**：人源阳性对照，15 个已知对 vs 全部未注释交叉对（209 对，不含 GH1–PRLR），
   两个受体长度组合并计算，平局计 1/2。
2. **正确受体排第一 ≥ 10/15**：每个配体在 15 个受体中，已知受体的分数**严格**高于其他所有受体
   （平局不算第一；GH1 的竞争者不含 PRLR）。
3. **小鼠面板方向一致**：8 个已知对中**每一对**在全长和胞外段的方向都相同（8/8；方向为 0 算不一致）。
   方向的定义：双中心化用 z 的符号；原始分数用"分数 − 同一区段面板全部 48 对分数的中位数"的符号。

### 判定和后续

- 4 个候选中有至少一个三条全部满足：报告该方法，**不自动进入 OT**，讨论后再决定。
  因为一共比较了 4 个方法，单个方法通过时证据要打折扣。
- 4 个全部不通过：在本文件写成负面结果，PPI 部分停止，不再换 checkpoint、变换、序列区段或阈值。
- 结果目录 `runs/preregistered_01/`：`criteria.csv`、`ligand_ranks.csv`、`panel_directions.csv`、`verdict.json`。

### 结果（2026-11-14，计算于预先登记提交 `b262538` 之后）

**判定：4 个候选全部不通过，而且没有任何一个候选满足三条中的任意一条。按预先登记，PPI 部分到此停止。**

| 方法 | 候选 | AUROC（≥0.85） | 排第一（≥10/15） | 方向一致（8/8） |
|---|---|---|---|---|
| V11 原始 | 参照 | 0.796 | 4 | 4 |
| V11 双中心化 | M2 | 0.841 | 6 | 2 |
| V12 原始 | 参照 | 0.807 | 5 | 5 |
| V12 双中心化 | M3 | 0.819 | 4 | 6 |
| Leakage-Free 原始 | M1 | 0.525 | 2 | 6 |
| Leakage-Free 双中心化 | M4 | 0.539 | 0 | 4 |

- 双中心化让 V11 的 AUROC 和排名略有提高（0.796 → 0.841，4 → 6），V12 没有改善；两者的方向一致性都不达标。
  也就是说，扣掉单个蛋白的平均效应之后，剩下的"配对特异"部分在全长和胞外段之间仍然不稳定。
- Leakage-Free 对几乎所有组合都给 0.6–0.7 左右，已知对和交叉对分不开（两组 AUROC 均约 0.53）。
- V11 双中心化的 AUROC 接近 0.85，但排名和方向两条都差得远，不存在"差一点通过"的情况。

**与预先登记的偏离（安装检查，不涉及三条通过标准）：** Leakage-Free 的小鼠自检 AUPR 0.238 / AUROC 0.636，
低于登记的 0.5。规则假定这只可能是安装问题。为了判断，在该 checkpoint 自己的训练基准
（Bernett 金标准测试集，HF `danliu1226/Bernett_benchmarking`@`5d2ad03…`，1:1，抽 2000 对，不截断）上自检：
AUPR 0.661 / AUROC 0.680（随机约 0.49）。权重以 strict 方式加载到与 V11/V12 相同的基础模型和代码，
V11/V12 自检正常。据此判断安装没有问题，小鼠上的低分是模型本身的性质，没有可修复的安装错误，于是按原计划完成评价。
这一偏离只可能让 Leakage-Free 更难通过，不会偏向通过。

**负面结果（写入论文时的表述）：** 在 15 对经典人源配体–受体阳性对照和 8 对小鼠面板上，
PLM-interact 的三个官方 650M checkpoint（humanV11、humanV12、Leakage-Free），无论原始分数还是
logit 尺度双中心化，都不能同时做到：区分已知配对与未注释组合、为每个配体把正确受体排在第一、
以及在全长与胞外段之间保持一致的方向。因此这些分数不能作为细胞间通信 cost 中的配对相容性项 S。
依据 `runs/preregistered_01/`（`criteria.csv`、`ligand_ranks.csv`、`panel_directions.csv`、`verdict.json`，
输入副本在 `inputs/`）。

不再进行的工作：更换 checkpoint、更换或组合变换、更换序列区段、调整阈值。
