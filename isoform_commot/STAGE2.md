# 第二阶段：Visium、稀疏求解

决定和待确认事项见 [reports/stage2_plan.md](reports/stage2_plan.md)（在任何第二阶段 OT 之前提交）。
本文件只说明怎样运行，以及稀疏实现是怎样验证的。

## 运行

```bash
conda activate commot-research
cd /home/liu/Research/isoform_commot

# 数据（固定 SHA256）
python stage2.py fetch --sample posterior
python stage2.py fetch --sample anterior

# 小区域：中心附近 128 个 spot，每类 6 对（只用于与稠密版本对比）
python stage2.py prepare --data data/visium_mouse_brain_sagittal_posterior --output stage2_input_post_roi128 \
    --max-spots 128 --pairs-per-type 6
python stage2.py run --input stage2_input_post_roi128 --output stage2_results_post_roi128
python stage2.py compare-dense --input stage2_input_post_roi128 --results stage2_results_post_roi128

# 整张切片：输入与全转录组表达
IN=stage2_input_post
python stage2.py prepare --data data/visium_mouse_brain_sagittal_posterior --output $IN
python stage2.py export-genes --input $IN --output stage2_genes_post

# 计划第 8 节的条件（--eps、--rho、--diffusion-cutoff-um 覆盖 config.json，实际取值写入 summary.csv 和 manifest.json）
python stage2.py run --input $IN --output stage2_results_post_A    --eps 0.1 --rho 10 --diffusion-cutoff-um 363.0757
python stage2.py run --input $IN --output stage2_results_post_A500 --eps 0.1 --rho 10 --diffusion-cutoff-um 500
python stage2.py run --input $IN --output stage2_results_post_Aeps --eps 0.2 --rho 1  --diffusion-cutoff-um 363.0757
SAME="--weights 1 0 0 0 --nitermax 50000 --stopthr 1e-11"    # 纯 collective，求解设置与 A 相同
python commot_baseline.py --input $IN --output commot_stage2_post_B    --dis-thr-um 363.0757 $SAME
python commot_baseline.py --input $IN --output commot_stage2_post_B500 --dis-thr-um 500 $SAME
python commot_baseline.py --input $IN --output commot_stage2_post_Beps --dis-thr-um 363.0757 $SAME --eps 0.2 --rho 1
python commot_baseline.py --input $IN --output commot_stage2_post_C    --dis-thr-um 363.0757   # 按原文：默认混合
python commot_baseline.py --input $IN --output commot_stage2_post_C500 --dis-thr-um 500

# 靶基因评价（scSeqComm；复合物受体"任一亚基有关联"）
python evaluate_targets.py run --input $IN --genes stage2_genes_post --output target_eval_stage2_post \
    --results stage2_results_post_A stage2_results_post_A500 stage2_results_post_Aeps \
              commot_stage2_post_B commot_stage2_post_B500 commot_stage2_post_Beps commot_stage2_post_C commot_stage2_post_C500

# 计划第 8 节的汇总：条件表（运输比例、检查、靶基因指标）、r0 选择（第 8.5 节规则）、受体配对比较、图
python stage2.py report --eval target_eval_stage2_post --output reports/stage2_post_01
```

## 输入目录（`prepare` 输出）

| 文件 | 内容 |
|---|---|
| `expression.csv` | spot × 实体；复合物为各亚基（normalize_total + log1p 后）的最小值 |
| `spatial.csv` | x、y（µm），`array_row`、`array_col`（六边形邻接） |
| `lr.csv` | 通过 COMMOT 筛选的配对：ligand、receptor、type（diffusion / contact）、通路 |
| `entities.csv` | 实体及其亚基 |
| `filter_audit.csv` | 每个 Secreted / Contact 配对是否通过筛选及原因 |
| `config.json`、`preparation.json` | 参数；来源、归一化目标值、µm/px 及邻距核对、文件哈希 |

## 稀疏结果格式（`plans.py`）

`case_*.npz`：`format='coo'`、`shape`、`row`、`col`（只存可连接的元素）、`C`、`T`、`a`、`b`、`spots`、`ligands`、
`receptors`；第二阶段的结果另有 `D`（距离，µm）和 `pair`（每个元素属于哪个 LR 对）。
行号 = 配体序号 × n + 发送 spot，列号 = 受体序号 × n + 接收 spot。没有存储的元素即禁止连接。
发送或接收预算为 0 的元素不存储（它们不能承载运输）。`plans.load` 同时读取第一阶段的稠密文件和稀疏文件。

## 验证

- **求解器：** `research.solve_sparse` 调用固定版本 `unot_sinkhorn_l1_sparse`，与稠密版本 `solve` 使用相同的预算缩放、
  零预算行列剔除和收敛参数（`nitermax = 50000`、`stopthr = 1e-11`）。
- **质量守恒检查：** 稀疏函数不返回对偶变量。最后一步更新后列约束精确成立，因此取 ν = b − Tᵀ1，由
  T = μ·ν·exp((2ρ − C)/ε) 反推 μ，检查 |T1 + μ − a| < 1e-7 和每个元素的乘积形式。ν 低于舍入误差的列
  （ρ/ε 很大时，如 COMMOT 的 ρ = 10、ε = 0.1）无法检查，结果中报告实际检查的行列数（`check_rows_checked`）。
- **覆盖全部行列的最优性检查（`research.sparse_kkt_check`）：** 在 ρ = 10、ε = 0.1 下，未运输质量约为 1e-46，
  上面的检查在 128 spot 小区域上只能覆盖 828 行中的 217 行。因此另加一项不需要分辨 μ、ν 的检查：
  求解器返回 T = exp((f + g − C)/ε)，于是 log T + C/ε 在每个连通分量内可加（φ_i + ψ_j），沿生成树求出 φ、ψ
  （每个分量差一个常数 γ），其余元素检验可加形式（`kkt_shape_error`，对数尺度 < 1e-8）；
  μ = exp(φ + γ − ρ/ε)、ν = exp(ψ − γ − ρ/ε)。每一行 |T1 + μ − a| ≤ tol、每一列 |Tᵀ1 + ν − b| ≤ tol 都给 γ 一个上界或下界，
  每个分量的上下界相交才通过（tol = 1e-7，与原检查相同）。`kkt_tol` 报告 1e-7、1e-9、…、1e-15 中能通过的最小值。
  COMMOT 纯 collective（条件 B）用 COMMOT 内部归一化后的 cost（d / 截断内最大距离）做同样两项检查。
- **单元测试（`test_stage2.py`）：** 三组 ε/ρ 下与稠密版本逐项一致（≤ 1e-10 相对误差）；两项检查都能发现偏离最优的行和错误的 cost，
  全行列检查能发现未收敛的解（3 次迭代），且在 ρ = 10、ε = 0.1 下覆盖原检查跳过的行；
  COMMOT `cot_sparse` 与 `cot_combine_sparse(weights = 1 0 0 0)` 结果逐项相同；
  `run` 的 ε、ρ、截断覆盖值被 compare-dense 正确使用；
  每对稀疏 cost 与从完整距离矩阵独立构建的稠密 cost 完全一致；端到端 run + compare-dense；稀疏和稠密计划读出的接收信号一致。
- **Visium 小区域（128 spot、12 对）：** 见下方结果。

## 第一次决定下的探索性结果（2026-11-14，后部切片，S = 1，ε = 0.2，ρ = 1，截断 500 µm）

这些运行在第二次决定（计划第 8 节）之前完成，只用于验证稀疏实现；它们暴露了第 8.4 节的规模问题。

**小区域与稠密版本对比**（128 spot，6 个 Secreted + 6 个 Contact 配对，7 个配体 × 12 个受体，35,767 个元素）：
四个 r0（50、100、200、400 µm）全部通过——可连接元素集合完全相同，cost 差 0，
T 的最大相对差 1.7e-15 – 2.2e-15（标准 ≤ 1e-6），禁止连接上的运输为 0；828 行全部通过质量守恒检查。

**整张切片**（3355 spot，筛选后 540 对：269 Secreted、271 Contact；201 个配体 × 224 个受体；
6,449,494 个元素，其中 diffusion 5,378,091、contact 1,071,403；106 个配体、126 个受体被多个配对共享）：

| r0 (µm) | 总运输量 | diffusion | contact | diffusion 平均距离 (µm) | 求解时间 |
|---|---|---|---|---|---|
| 50 | 42,653 | 5,967 | 36,686 | 4.0 | 8 s |
| 100 | 43,724 | 7,055 | 36,669 | 33.1 | 8 s |
| 200 | 48,059 | 11,462 | 36,596 | 111.9 | 8 s |
| 400 | 56,047 | 19,608 | 36,439 | 192.4 | 8 s |

187,622 行、215,705 列全部可检查并通过（行误差 1.6e-16，乘积形式误差 1.1e-20）。每个 case 文件约 83 MB。

**需要讨论（计划第 7 节第 3 条 ε、ρ）：** 运输比例随问题规模变化很大。小区域上约 78%（配体预算），
整张切片上约 20–26%。原因是预算先除以总量再求解（COMMOT 也这样做），问题越大，每个元素的质量越小，
而熵项和 ρ 不随之缩放。因此在第一阶段小区域上选定的 ε = 0.2、ρ = 1，在整张切片上的含义不同；
在比较方法或解读整张切片的结果之前，应先确定 ε、ρ 在整张切片上的取法。
