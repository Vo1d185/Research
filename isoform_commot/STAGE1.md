# 第一阶段：固定 CCCvelo LR 分类，测试距离参数

本阶段使用 CCCvelo 已分好的 diffusion/contact LR 数据库，不新增 CellChat LR 对。
还需要独立下载论文的空间表达数据；两张 LR 表不是表达矩阵或空间坐标。
第二阶段再引入 CellChatDB 注释和更多样本；CellChat 不会生成真实测量的表达或坐标。

目前实现的是**第一阶段输入适配与参数实验程序**。真实组织的参数结论须等下载、检查并运行数据后获得。
先得到稳定的候选参数范围，不根据最大的 total_T 宣称 r0 最优。

## 1. 更新 GitHub 代码，沿用你的 conda 环境

在 VS Code SSH 终端逐行运行，不要把两条命令粘成一行：

```bash
cd /home/liu/Research
git fetch origin
git switch research/isoform-commot-v1
git pull --ff-only
conda activate commot-research
cd isoform_commot
python -m pip install -r requirements-stage1.txt
python fetch_solver.py
```

本项目测试使用 Python 3.12。不要在同一个终端叠加激活旧 `.venv`。
云端开发使用项目虚拟环境，不影响你的 conda 使用方式。

## 2. 下载两类 LR 表和论文的小鼠脑数据

```bash
python stage1.py fetch-lr
python stage1.py download-data
```

下载位置：

- `data/cccvelo_lr/diff_LigRecDB.csv`
- `data/cccvelo_lr/cont_LigRecDB.csv`
- `data/mousebrain_bin60.h5ad`

LR 表固定为 CCCvelo 提交 `fac213c600e7929b5d14e7d95dbdaabd722f9bc9`，验证 SHA256。
这版有 3883 个 diffusion 对、109 个 contact 对，两表没有重复 LR 对。
这是**当前固定的公开表**，不要与论文正文所述 3659 条先验 LR 或某次分析的表达筛选结果混为一谈。

h5ad 链接来自论文 Data availability（印刷页 205），并已按原始提供方
[Spateo 的下载代码](https://github.com/aristoteleo/spateo-release/blob/1bd8a35e6e7b0dfc712d76b1c32987dc06652774/spateo/sample_data.py#L157)
逐字符核对：
<https://www.dropbox.com/s/c5tu4drxda01m0u/mousebrain_bin60.h5ad?dl=1>。
修正了初版扫描件录入错误：首字符是 `c`，`xda0` 后面是数字 `1`，不是字母 `l`。
该链接在当前云环境未完成下载验证；需要在能访问 Dropbox 的 SSH 服务器或本机执行。
下载器验证 HDF5 文件签名、响应长度（若提供）、HDF5 可读取性和 obs/var 元数据，
防止把 HTML 页面或不完整文件当数据，并记录本地 SHA256。下载中每 64 MiB 打印进度。
本地哈希只用于复现，不是作者发布的真实性校验值。
若链接失效，不能用模拟矩阵代替论文数据；从论文/原作者核实新地址后可传 `--url '新地址'`。
也可以浏览器下载后自行放到上述路径，跳过 `download-data`。

使用旧代码时，也可以通过已支持的 `--url` 立即指定核实后的正确地址：

```bash
python stage1.py download-data --url 'https://www.dropbox.com/s/c5tu4drxda01m0u/mousebrain_bin60.h5ad?dl=1'
```

## 3. 先检查文件，不猜测它的表达层

```bash
python stage1.py inspect --h5ad data/mousebrain_bin60.h5ad --output stage1_inspect.json
```

检查报告列出 X、raw、layers、空间坐标键、基因 ID、细胞注释和小样本数值范围。
可以先把报告发来确认。整数样本不一定证明是 raw counts，`raw` 这个名称也不证明存的是原始计数。
已有日志归一化/插补表达用 `processed`；确认原始计数后才用 `counts`。
坐标来源没有核实时用 `source_unit_unverified`，不能直接写成 micrometer。
本阶段 r0/h 是无量纲比例，因此可以先做计算敏感性实验，但不能将结果报告为真实传播微米数。

## 4. 生成固定的小区域实验输入

以下示例**以确认 X 是原始 counts、空间键为 spatial 为前提**：

```bash
python stage1.py prepare --h5ad data/mousebrain_bin60.h5ad --matrix X --expression-state counts --spatial-key spatial --coordinate-unit source_unit_unverified --symbol-policy author-uppercase --output stage1_input_01
```

若计数在 `layers['counts']`，改为 `--matrix layer:counts`。
若选的是已处理的非负表达，改为 `--expression-state processed`，不再归一化。
若坐标键是 `X_spatial`，改为 `--spatial-key X_spatial`。
若 var_names 是 Ensembl ID 而基因符号在 var 的列中，使用 `--gene-key 实际列名`。

`author-uppercase` 对照作者 Demo 的 `toupper` 做符号匹配，**不等于经过验证的人鼠同源映射**。
第一版仅用于计算流程和参数敏感性。若已有可靠的一对一同源基因表，请替换为
`--gene-map orthologs.csv`，列为 `input_gene,lr_gene`；歧义映射会拒绝执行。
多亚基复合物目前排除并记录，不能把缺少的复合物当成单个蛋白。

默认选择靠近组织坐标中位位置的真实 spot 为中心，再取最近的 128 个位置，保留所有细胞类型。
这不是随机抽样，也不保证中心就是皮层。**先看生成的 roi.png**，必要时用
`--center X Y` 选择皮层中的连续区域，在新的输出目录生成输入。
外部发送者被裁掉，不能当作整张组织或 CCCvelo 整个发送群体的复现。

在该 ROI 中，先按全转录组库大小归一化，再选 LR；默认表达比例 >=5%、平均处理后表达 >=0.05。
每类最多取 4 对，按表达支持比例、基因名确定，完全不使用 OT 结果或 PPI 分数选对。
这些筛选值是本实验设置，不是宣称照搬论文统一参数。固定面板后比较全部 r0。
如果某一类型没有可用对，会输出 `status=blocked` 和 `lr_audit.csv`，不猜测类型或编造表达。
稠密实现限制 `n_ligand*n_receptor*n_spots^2 <= 4000000`。

输入目录中的关键文件：

| 文件 | 用途 |
|---|---|
| expression.csv、spatial.csv、lr.csv | OT 所需的固定表达、坐标和类型 |
| config.json | 实验参数 |
| preparation.json | 数据哈希、ROI、映射、归一化、候选数量和限制 |
| lr_audit.csv | 每对 LR 入选或被排除的原因 |
| cell_meta.csv、roi.png | 检查区域和注释 |

空间约束采用 Delaunay 邻接并去掉长度 >2h 的 contact 边；h 为 ROI 中位最近邻间距。
2h 是可修改的 pilot 设置（`--contact-radius-h`），不是作者代码 200/0.73 的物理阈值。
排除 self edges，不使用细胞分割或 shared cell boundary。
diffusion 暂不加额外截断，在整个选定 ROI 中考虑非自身发送者。
邻接只代表空间代理条件，不证明实际膜接触。

## 5. 距离实验与 OT 参数检查

第一轮固定 eps=0.2、rho=1：

```bash
python stage1.py run --input stage1_input_01 --output stage1_results_01
```

比较五个条件：

- 对照：`F=d_ref/d`，d_ref 为最小正距离，避免对有量纲的 1/d 直接取 log。
- `F=r0/(r0+d)`：r0/h = 0.5、1、2、4。

cost 为 `-log(F*S)`，contact 合法边 F=1，非法边禁止运输。没有 PPI 文件时 S=1。
此模型和无量纲化是我们的研究设计，不是 CCCvelo 的原始 LRS 或 COMMOT 原始距离 cost。
非平衡 OT 对 cost 的常数偏移敏感，r0 会同时改变整体强度和曲线形状。
已有 `r0_shape` 作为可选的匹配幅度对照保留在通用引擎中，第一轮默认不额外展开。

检查通过后，再在同一个输入上运行九组 eps/rho：

```bash
python stage1.py run --input stage1_input_01 --output stage1_results_sensitivity_01 --sensitivity
```

eps={0.1,0.2,0.5}，rho={0.5,1,2}，共 45 个距离条件（无 PPI 时）。
查看：

- `all_conditions.csv`：总 T、未匹配预算、平均距离和距离对照变化。
- `all_types.csv`：diffusion/contact 分开后的结果。
- `ot_parameter_stability.csv`：相同距离和分数条件，改变 eps/rho 后相对于 eps=0.2、rho=1 的分配 TV 和 L1。
- `ot_00/analysis/distance_by_type.png`：第一组 OT 设置下的距离曲线。
- 每个 `ot_*/results/manifest.json`：完成状态、固定输入哈希、求解器版本。

评估：质量守恒/禁止边是否满足；diffusion 距离分布怎样随 r0 变化；contact 是否因共享预算竞争改变；
结论是否只在某组 eps/rho 下出现；换独立 ROI 后方向是否一致。
先把第一轮视为 pilot，另选验证 ROI 再缩小候选范围，不把最大总 T、最低 cost 或最平滑曲线当成准确率。

## 6. 接入真实 PPI 分数

拿到真实 PLM-interact 结果后，用新输入目录执行同样的 prepare，额外添加：

```text
--ppi real_scores.csv --score-provenance "实际模型版本和序列来源"
```

文件列为 `ligand,receptor,score`，必须覆盖所选面板，不能用 1 填缺失真实分数。
会自动比较无 PPI、真实分数、同类型内置换分数；固定距离设置后解释差异。
若要测试竞争，检查 preparation.json 中 shared_ligands/shared_receptors；只有互不共享实体的 LR 时，
分数变化可能主要影响各 LR 的通信总量。
当前没有运行 PLM-interact，也没有把模拟分数写进真实数据输入。
用 gene-level 输入完成上述步骤，仍不能证明 isoform 特异作用或真实结合频率。

## 方法依据

- CCCvelo：https://doi.org/10.1038/s43588-025-00934-2
- 固定软件版本：https://github.com/SunXQlab/CCCvelo/tree/fac213c600e7929b5d14e7d95dbdaabd722f9bc9
- 公开 Demo 使用 SCT data -> ALRA，LRS 读取 Imputate 层；本项目第一阶段使用明确记录的非负表达适配。
- 本项目保留官方 COMMOT L1-UNOT 核心，仅对允许连接上的 cost 做实验；并非整个 CCCvelo 或 COMMOT 高层默认流程复现。
