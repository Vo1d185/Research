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
