# 第二阶段：前部切片的配对清单及与后部的交集

由 `python stage2.py pairs --posterior stage2_input_post --anterior stage2_input_ant --output reports/stage2_pairs_01`
生成（代码 `fa848c3`，工作区无修改；见 `pairs.json`）。两张切片按计划第 2 节的同一规则**各自独立**预处理和筛选
（COMMOT Visium 教程：`normalize_total` + `log1p`，每个亚基在本切片 ≥5% 的 spot 中检测到，复合物取亚基最小值），
每张切片的 OT 只用它自己的配对清单。

## 数据

| | 后部（主数据集） | 前部（独立重复） |
|---|---|---|
| 样本 | `V1_Mouse_Brain_Sagittal_Posterior` | `V1_Mouse_Brain_Sagittal_Anterior` |
| 原始文件 SHA256 | matrix `df962209…`，spatial `b878d7ea…` | matrix `56078d8d…`，spatial `5f41a803…`（与计划第 1 节一致） |
| 组织内 spot | 3355 | 2695 |
| µm/像素 | 0.72615 | 0.72648 |
| 相邻 spot 中心距 | 100.21 µm | 100.25 µm |

来源：10x Space Ranger 1.1.0 公开数据（`https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/<样本>/`），2026-11-14 下载。

## 配对数

| 类型 | 后部 | 前部 | 交集 | 并集 | Jaccard | 受体：后部 / 前部 / 交集 |
|---|---|---|---|---|---|---|
| Secreted（diffusion） | 269 | 283 | 235 | 317 | 0.74 | 115 / 129 / 114 |
| Contact | 271 | 293 | 266 | 298 | 0.89 | 118 / 127 / 117 |
| 合计 | 540 | 576 | 501 | 615 | 0.81 | 224 / 247 / 222 |

只在一张切片上通过的配对，在另一张切片上都是因为某个亚基低于 5% 被排除（`pairs_union.csv` 的 `*_filter` 列）：
只在后部的 39 对（Secreted 34、Contact 5），只在前部的 75 对（Secreted 48、Contact 27）。

## 截断处的 spot 对（计划第 8.3 节）

| 切片 | 截断 | 0 < d ≤ 截断的 spot 对 | 截断内最远 | 截断外最近 |
|---|---|---|---|---|
| 后部 | 363.08 µm | 84,757 | 362.67 µm | 399.38 µm |
| 前部 | 363.08 µm | 67,231 | 362.84 µm | 399.57 µm |
| 后部 | 500 µm | 131,078 | 499.59 µm | 500.32 µm |
| 前部 | 500 µm | 104,419 | 499.82 µm | 500.55 µm |

两张切片上 363 µm 与 363.08 µm 给出完全相同的 spot 对。

## 文件

`pairs_posterior.csv`、`pairs_anterior.csv`（各自的配对清单）；`pairs_union.csv`（并集，含是否在每张切片上通过及原因）；
`counts.csv`；`spot_pairs.csv`；`pairs.json`（输入哈希与启动时的 git 状态）。
