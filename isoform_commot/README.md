# Isoform-aware COMMOT：第一版实验代码

这是可复现的数值原型。**尚未完成 CCCvelo 论文数据集/QC/预处理复现，也没有运行真实 PLM-interact 模型。**
模拟 PPI 分数只用于检查接入。T 是模型通信量，不是实际结合次数或频率。

现在进入 [第一阶段：固定 CCCvelo LR 分类的参数实验](STAGE1.md)，包含下载、h5ad 检查、ROI 输入适配及 OT 参数敏感性运行。
这是为本研究 cost 设计的适配器，不是完整 CCCvelo SCT/ALRA 预处理复现。

## 在 VS Code 的 SSH 远程终端运行

代码和 Python 都运行在你的 SSH 服务器上，不需要连接 Codex。
使用 Python 3.11 或 3.12（依赖在当前环境中验证，旧 Python 请另建环境）。

```bash
conda activate commot-research
cd Research/isoform_commot
python -m pip install -r requirements-stage1.txt
python fetch_solver.py
python -m unittest -v
python research.py --demo --input demo_input --output results_demo
```

输入和输出目录已存在时程序拒绝覆盖，请给重复实验使用新目录名。
COMMOT 源文件按固定提交下载并验证 SHA256；上游源文件不随本仓库重新分发。
来源：https://github.com/zcang/COMMOT/blob/d117445bc07eaa19109c7609e97d9e35b26e99ca/commot/_optimal_transport/_unot.py

## 三组实验

1. 倒数距离对照 `F=d_ref/d`。`d_ref` 是最小正空间距离，用于使 F 无量纲。
2. `F=r0/(r0+d)`，r0 为中位最近邻间距 h 的 0.5、1、2、4 倍。
3. 形状对照 `F=(r0+d_ref)/(r0+d)`：把不同核在同一距离处归一化。

非平衡 OT 对 cost 的常数偏移敏感，所以不能把 `d_ref/d` 与未经单位处理的 `1/d` 称为完全相同的实验。
r0 同时影响衰减形状和整体 cost；第三组帮助区分这些影响。
扫描范围是研究假设，不是已确定的最佳生物学长度尺度。

每组比较 S=1、输入 PPI 分数、同类型内打乱分数。少量 permutation 只用于流程检查，不能支持显著性结论。
存在多个共享 ligand/receptor 的候选对时才能充分检查分配竞争；一个独立 LR 对的常数 score 可能主要改变总量。

## Cost 和求解器

合法连接采用 `C=-log(F)-log(max(S,1e-8))`；非法连接为 +inf。
contact-type 先通过邻接限制，再令 F=1；零 cost 的合法连接保留。
默认 contact 半径 1.1h、diffusion 半径 5h 都是可调整的初始假设，**不是已核实的 CCCvelo 参数**。
同 spot 连接默认排除；spot 邻接并不证明两个细胞的膜实际接触。

调用官方 `unot_sinkhorn_l1_dense`，目标是：

`min <C,T> + rho*(sum(mu)+sum(nu)) + eps*(h(T)+h(mu)+h(nu))`

约束 `T 1 + mu = a`、`T^T 1 + nu = b`，且非负；`h(x)=sum(x*(log(x)-1))`。
三项熵系数相同。表达量统一除以 max(sum(a),sum(b)) 后求解，T 再还原尺度；此尺度在全部条件中相同。
每次实验不重新按 cost 最大值缩放，以免改变 eps/rho 的实际含义。
这实现了 collective OT 核心，不等价于 COMMOT 高层 API 默认的多种运输方案混合与全部下游步骤。
程序直接检查 dual potentials 对应的质量平衡；测试还使用独立的 SLSQP 原始目标求解作对照。

## 真实输入

目录中放置以下文件，可参考 `--demo` 生成的格式：

| 文件 | 格式 |
|---|---|
| expression.csv | 第一列 spot ID，列名 gene/isoform ID，数值非负 |
| spatial.csv | 第一列 spot ID，以及 x、y；与表达矩阵 ID 完全一致 |
| lr.csv | ligand,receptor,type；明确指定 diffusion/contact |
| ppi.csv | 可选；ligand,receptor,score；分数在 [0,1] 且覆盖全部使用的 LR 对 |
| config.json | 见下方配置 |

```json
{
  "coordinate_unit": "micrometer",
  "expression_level": "gene",
  "expression_state": "raw_counts",
  "full_transcriptome": true,
  "score_provenance": "填写实际模型版本、序列或结果来源",
  "r0_h": [0.5, 1, 2, 4],
  "contact_radius_h": 1.1,
  "diffusion_radius_h": 5.0,
  "eps": 0.2,
  "rho": 1.0,
  "seed": 20261112,
  "n_permutations": 3
}
```

raw_counts 使用全转录组的 library-size 归一化到 10000 后 log1p；这是通用初始选项，尚未与 CCCvelo 核对。
已经处理过的非负表达使用 `processed`，避免重复 log。不要输入 z-score 后的负表达。
类型可接受 CellChat 的 `Secreted Signaling` / `Cell-Cell Contact` 标签；ECM 等其他类别需显式
`type_mapping` 与 `type_mapping_provenance`，不自动猜测。
未表达/缺少实体的 LR 对输出到 excluded_lr.csv；多亚基复合物未实现自动处理，不能把缺少的复合物当作单蛋白。
isoform 数据使用 `expression_level=isoform` 并提供匹配的实体 ID；不要把一个 gene 的完整表达复制到多个 isoform。
PPI score 是先验相容性信息，不是蛋白质亲和力或量化的结合频率。序列评分模型的校准和定位信息仍需研究。

```bash
python research.py --input my_data --output results_real_01
```

## 输出和限制

- summary.csv：各 cost/score 条件的总通信量、未匹配量、平均距离、相对对照的 L1 变化。
- case_*.npz：实际使用的 C、T、表达量预算、spot 和实体顺序；行/列索引为 entity_index*n_spots+spot_index。
- case_*_spots.csv：按实体求和后的 spot 通信网络。
- manifest.json：完成标记、输入哈希、源代码哈希与参数；没有 complete manifest 的目录不能当作完整实验。
- 没有 ppi.csv 时只运行距离实验，并明确记录 PPI 尚未测试。

采用稠密矩阵，限制 400 万 cost 元素。先用小型 ROI 和少量 LR 对，完整组织需要后续稀疏实现。
目前没有最佳 r0、生物学有效性、真实 isoform 通信发现或 CCCvelo 完整复现的结论。

## 运行完模拟实验之后

使用原来的输入目录与结果目录；分析脚本不会重新计算 OT，也不需要重新下载求解器。

    python analyze.py --results results_demo --input demo_input --output analysis_demo

如果输入/结果目录叫 demo_input_02、results_demo_02，请相应修改参数。输出目录必须是新目录。
脚本核对原输入哈希、spot/实体顺序和预算；输入后来被编辑过时，应恢复原输入或重新运行实验。

依次打开以下图片：

1. distance_by_type.png：没有 PPI 时，分别看 diffusion/contact 总量和平均距离。
   contact 的核不随 r0 变化，但其 T 可以通过共享 receptor/ligand 的竞争而改变。
   这个等间距直线 demo 的 contact 平均距离为 1，不是对真实组织的普遍结论。
2. ppi_matched.png：固定同一个 family、r0，比较 PPI 和 no_ppi，区分总量变化与分配变化。
3. lr_pair_changes.png：展示 r0_shape 中最接近 r0/h=1 的条件下各 LR 对的量。
   该条件仅用于展示，不表示最优 r0；表格包含所有条件。

表格：

- by_type.csv：每个条件、每种类型的总 T、平均距离和总量占比。
- by_lr.csv：每个条件、每个 LR 对的总 T 和平均距离。
- ppi_matched.csv：PPI 与同一距离设置下 no_ppi 的配对比较。
- lr_pair_changes.csv：同一距离设置下，各 LR 对获得/失去的通信量。
- identical_plans.csv：同一距离设置下完全相同的运输矩阵；相同 T 本身不能证明输入 score 相同。
- analysis_manifest.json：完成标记与分析所读取结果文件的哈希。

paired_relative_L1 = sum(abs(T_ppi-T_no_ppi))/sum(T_no_ppi)，
分母是**同一个距离设置**的 no_ppi，与旧 summary.csv 的全局 case_000 基准不同。
allocation_TV = 0.5*sum(abs(T_ppi/sum(T_ppi)-T_no_ppi/sum(T_no_ppi)))，
取值 0–1：只把 T 整体缩放时为 0；它不是准确率或改变的边数比例。
任一矩阵总量为 0 时，TV 留空，因为归一化分配没有定义。
一个实体可能参与两类通信，因此不把未匹配预算武断地拆分到两种类型。
