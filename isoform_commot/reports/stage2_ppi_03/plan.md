# 实验二扩大版:COMMOT 论文的数据集(stage2_ppi_03)

**状态:预注册计划。** 2026-11-21 定稿。本文件提交后才运行第 1 步(数据与 V12 打分)和第 2 步(OT 与评价);
两个脚本的正式运行都检查本文件已提交且没有改动。

**依据与用户决定(2026-11-21):**
- 实验二(`reports/stage2_ppi_02`,计划 4e3874f)只用了两张鼠脑切片;用户要求用更多数据检验结论,并且使用 COMMOT 论文
  (Cang et al., *Nat Methods* 20, 218–228, 2023)里的数据,便于和论文比较。
- 设计由用户选定:**第 0 步复现 + A 线判定 + B 线对照**;A 线加入**长度对照 S**(描述性)。
- β 已删除(固定为 1,用户 11-20 决定);继续使用 PLM-interact humanV12。

## 0. 问题

与实验二相同的三个问题,在新数据上回答:
- **D1 信息:** 原始 S、S′ 各自是否带有配对信息(与同类型内打乱的分数比较)?
- **D2 改进:** 原始 S、S′ 各自是否比不加 S 更好?
- **D3:** 原始 S 与 S′ 哪个更好?

另外描述:1/k、长度对照、COMMOT 原版在同一评价下的位置;以及在论文自己的流程和指标下,各条件与论文报告值的关系(B 线)。

## 1. 数据

论文用下游靶基因做定量验证的只有 3 套数据(正文 "Robust identification of CCC direction and downstream target",
补充图 21;代码 `cell-level_downstream_gene_correlation.ipynb`)。全部来自论文的 Zenodo 存档 10.5281/zenodo.7272562(v0.0.2)。

| 代号 | 数据 | 物种 | 单位 | 数量 | 基因 | 本实验中的角色 |
|---|---|---|---|---|---|---|
| `breast` | Visium 人乳腺癌 V1_Breast_Cancer_Block_A_Section_1 | 人 | spot | 3,798 | 36,601 | **A 线判定**;B 线 |
| `seqfish` | seqFISH+ 小鼠体感皮层(Eng et al. 2019),论文用的 523 个皮层细胞 | 小鼠 | 单细胞 | 523 | 10,000 | **A 线判定**;B 线 |
| `post` | Visium 鼠脑后部(实验二的判定切片) | 小鼠 | spot | 3,355 | 32,285 | A 线参照;B 线 |
| `ant` | Visium 鼠脑前部(实验二的确认切片;论文没用它做靶基因验证) | 小鼠 | spot | 2,695 | 32,285 | A 线参照 |

论文其余数据不用:Slide-seqV2 海马(论文没用于靶基因验证;53,173 个 bead)、STARmap(903 / 1,020 个基因)、
MERFISH(155 个基因)、果蝇胚胎与人皮肤(模拟数据)。

**来源与校验:**
- Zenodo 文件按 API 公布的 MD5 校验,并记录 SHA256(`exp2x_data.py fetch-zenodo`)。
- 乳腺癌的 10x 文件取自 `downstream_gene_correlation.zip` 内的副本:`filtered_feature_bc_matrix.h5` 的 SHA256 为 `b3cbb6ea…`,
  `spatial.tar.gz` 为 `415530a9…`(全文见脚本常量)。同一存档内的后部切片文件与本仓库已固定的 10x 文件**逐字节相同**,两个 SHA256 都一致。
- seqFISH+ 取自 `seqfishplus-mouse-sscortex.zip` 的 `original_data`:`expression.npz`、`genes.txt`、`positions.csv`(SHA256 见脚本)。
  坐标是像素,**103 nm/像素**(Eng et al. 2019 图注:"X-Y coordinates are in pixels (103 nm per pixel)")。
  中位最近邻距离为 103.1 像素,即 10.62 µm。
- CellChatDB v2:与小鼠相同的提交 `75253cd`;人库 `CellChatDB.human.rda` 的 SHA256 为 `582e99db…`。
- scSeqComm 人类表:与小鼠相同的提交 `95c9bb2`,包括 `TF_PPR_KEGG_human.rda`、`TF_PPR_REACTOME_human.rda`、`TF_TG_TRRUSTv2.RData`。
  人用 TRRUSTv2 是 COMMOT 的选择;小鼠用 TRRUSTv2 + RegNetwork High,与实验二相同。
  - 核对结果:`TF_TG_TRRUSTv2` 与论文存档导出的版本完全相同(795 个 TF,靶基因集合一致)。
  - 人 REACTOME 表比论文存档多:24,306 行,对论文的 21,565 行。论文的行全部包含在内,多出 2,741 行。
  - 小鼠两张表完全相同。
  - 处理:A 线用 `95c9bb2`(与小鼠来源一致);B 线用论文存档导出的版本(与论文一致)。

## 2. 第 0 步:复现论文的数值(B 线的前提)

| 步骤 | 做法 | 通过标准(预先固定) |
|---|---|---|
| 0a 指标 | 用论文的评价函数(照搬 notebook 的规则:REACTOME PPR ≥ 1e-5、复合物受体按 `A_B → A,B` 整体精确匹配、只在高变基因中取靶基因、NaN 记 0),作用于**论文自己输出文件中的接收信号** | 评价对数等于论文(192 / 164 / 519);中位数与论文正文相差 ≤ 0.001(0.237 / 0.180 / 0.230);均值与 notebook 打印值相差 ≤ 0.001(0.23194 / 0.19682 / 0.23122) |
| 0b 求解 | 用固定版本的 COMMOT 求解器(`commot_baseline.py`,四方案混合,ε = 0.1,ρ = 10,COMMOT 默认迭代设置),在论文的 LR 列表、论文的预处理和论文的截断下重跑 | 逐 LR 对比较接收信号(spot 上的 Spearman),中位数 ≥ 0.99 |

- 论文的截断是坐标单位 1000,换成物理距离并不统一:后部约 726 µm(0.726 µm/像素),乳腺癌约 366 µm(0.366 µm/像素),seqFISH+ 为 103 µm。
- 预处理按论文写出 dis1000 文件的那个 notebook:
  - Visium:normalize_total + log1p;
  - seqFISH+:**只做 log1p,不做 normalize_total**(`1-lr_signaling.ipynb` 第 12 格)。
  - `prepare-b` 把我们算出的 LR 基因表达与论文输出文件里的 X 逐值比较,最大差 > 1e-4 就停止。
- 评价用的表达尺度也按论文:鼠脑为高变基因原始计数,另两套为 log1p。
- 0a 或 0b 不通过时,B 线照常报告,但标记"未复现",不和论文数值直接对照。A 线不受影响。

## 3. A 线(判定)

### 3.1 cost 与参数

- C = log(1 + d/r₀) − β·log S_eff,r₀ = 50 µm,**β = 1**(不加 S 时 β = 0)。
- ε = 0.2,ρ = 1,扩散截断 363.0757 µm(物理距离,四套数据相同),nitermax 50,000,stopthr 1e-11,CPU 固定版本求解器。
  其余与实验二主设置相同。

### 3.2 新数据的输入规则(与 Stage 2 / 实验二相同,只补单细胞需要的部分)

- CellChatDB v2 的分泌型与接触型配对,按 COMMOT 的过滤规则(每个亚基在 ≥ 5% 的 spot 或细胞中检出)。
  - 表达:normalize_total + log1p;复合物取亚基最小值。
  - 无分数的配对:取同类型有分数配对的中位数(与实验二相同)。
- 乳腺癌:接触型配对连自身和六个相邻 spot(与鼠脑切片相同)。
- seqFISH+:接触型配对连自身和**距离 ≤ 2 倍中位最近邻距离**的细胞,即 2 × 10.62 = 21.24 µm。
  - 本地准备时,每个细胞的中位邻居数为 4,7 个细胞除自身外没有邻居。
  - 坐标 = 像素 × 0.103 µm。
- 靶基因评价与实验二相同:前 10 个 TF、平均表达、检出 ≥ 5%、至少 5 个靶基因、去掉 LR 基因、复合物"任一亚基"、控制受体表达的偏 Spearman。
  人用 §1 的人类表。

### 3.3 条件(每套数据相同,共 205 个)

| 代号 | S | 用途 |
|---|---|---|
| `no_ppi` | S = 1 | 不加 S |
| `raw` | V12 原始分数 | 原始 S |
| `snorm` | S′ = S / Σ_R S(同一配体的所有受体) | 老师的 S′ |
| `uniform` | 1/k | 描述性对照,不带配对信息 |
| `length_raw` | **长度对照**:同一配体内,把该配体有分数配对的 V12 分数按受体长度重新分配(最长的受体得到最大值;受体长度 = 最长亚基的 UniProt 规范序列长度;同长按受体名);分布不变 | 描述性:V12 带来的变化是否只是长度 |
| `raw_perm_k`、`snorm_perm_k`,k = 1…100 | 分数在同类型、有分数的配对之间打乱(与实验二相同),S′ 在打乱后再归一化 | 零分布 |
| `ref_commot_B` | COMMOT 纯 collective,ε = 0.2,ρ = 1,363.0757 µm(实验二的 Bε) | 参照(不分接触 / 扩散,在它自己有定义的受体上比较) |

- 新的两套数据的 V12 分数在第 1 步重新打分。
- 鼠脑两张切片沿用实验二的分数(`reports/stage2_ppi_01/scores`),这样参照部分与实验二完全可比。

### 3.4 指标

- **M**:共同受体集合上受体层偏 Spearman 的平均值。共同受体集合 = 该数据集中核心、打乱和长度对照条件下偏 Spearman 都有定义的受体;报告个数。
- 区间:对受体做 bootstrap,1000 次,`numpy.random.default_rng(1)`。

### 3.5 判定

**每套数据**(与实验二规则相同,β = 1):
- D1:p = (1 + #{M_perm ≥ M}) / (1 + N)。
- D2:Δ =(S − `no_ppi`)的受体平均差,bootstrap 95% 区间。
- D3:`raw` − `snorm`。

**判定只用两套新数据(`breast`、`seqfish`),合并:**
- **D1 合并:** 每套数据算 z = (M − 打乱均值) / 打乱标准差,两套相加得 Z_obs;每个种子 k 用同样的方式得到 Z_k(两套数据用同一组种子)。
  p = (1 + #{Z_k ≥ Z_obs}) / (1 + N)。p ≤ 0.05 为"有信息"。
- **D2 合并:** 两套数据所有共同受体的平均差(按受体数加权)。区间用分层 bootstrap:每次在每套数据内部对受体有放回抽样,1000 次,种子 1。
  - "改进":区间下限 > 0,**且**两套数据各自的平均差都 > 0;
  - "变差":区间上限 < 0,且两套各自都 < 0;
  - 其余为"不确定"。
- **D3 合并:** 同 D2 的规则,标签为"原始 S 更好 / S′ 更好 / 无法区分"。
- 结论表与实验二相同:D1 有信息且 D2 改进 → 改进了推断;D1 有信息、D2 不确定或变差 → 带有配对信息但不足以改进;D1 无信息 → 看不出配对信息。
- 鼠脑两张切片(实验二看过结果)只作参照。四套数据的合并结果单列,**标为描述性**。
- **失败规则:** 新数据集的核心条件(`no_ppi`、`raw`、`snorm`、`uniform`)求解检查不通过 → 合并判定记为"无法判定"。
  打乱条件失败则剔除并报告,N 相应减少。

### 3.6 门槛(在全部求解之前)

1. **回归门槛:** 在实验二的后部输入上重跑 `no_ppi`、`raw`,逐受体与实验二提交的
   `reports/stage2_ppi_02/summary/post_main_receptor_partial_spearman.csv` 比较,差 ≤ 1e-12,否则停止。
   (实验二的 `raw` 就是 β = 1,所以应当完全相同。)
2. **等价门槛:** 乳腺癌和 seqFISH+ 的 `no_ppi`、`raw` 另存完整结果,用 `evaluate_targets.py run` 原样评价,
   与本脚本逐受体比较,差 ≤ 1e-12(同实验二)。
3. **GPU 门槛(B 线):** `gpu_unot.py equivalence` 在后部 B 线输入上(ε = 0.1,ρ = 10,200 次固定迭代)与固定版本求解器比较,
   最大相对差 ≤ 1e-9。

### 3.7 预先记录的预期

- 实验二:D1 p = 0.79 / 0.99(原始 S)、0.98 / 0.96(S′);打乱 M 的标准差为 0.002–0.003。所以新数据上 D1 很可能仍为"无信息"。
- 乳腺癌是人类数据,V12 是用人类 STRING 训练的,这是 V12 最有利的情况;如果 D1 在这里也不成立,就不能归因于跨物种。
- seqFISH+ 只有 523 个细胞,噪声预计大于 Visium。

## 4. B 线(描述性,不判定)

### 4.1 协议

- LR 列表 = 论文输出文件中的 `df_ligrec`:COMMOT 自带的 CellChatDB v1,只有分泌型;过滤后 300 / 250 / 598 对。
- 截断 = 坐标单位 1000;预处理和评价尺度见 §2;指标 = 每个 LR 对一个 spot 层 Spearman,报告中位数和均值。

### 4.2 条件(每套数据)

| 代号 | 内容 |
|---|---|
| `paper_output` | 论文自己的接收信号(0a) |
| `commot_rerun` | COMMOT 原版重跑(0b) |
| `receptor_only` | 用受体表达代替接收信号(看指标中有多少来自受体表达本身) |
| `ours_no_ppi`、`ours_raw`、`ours_snorm`、`ours_uniform` | 我们的 cost(r₀ = 50 µm,β = 1),纯 collective,ε = 0.1,ρ = 10,论文的截断;GPU,nitermax 50,000 |

- V12 分数在第 1 步对 v1 列表单独打分(人、鼠各一组)。
- 报告每个条件的 n、均值、中位数,以及逐 LR 对的配对差(相对 `ours_no_ppi` 和 `commot_rerun`)及 bootstrap 区间。都只作描述。

## 5. 描述性分析(不判定)

1. 前 5、前 20 个 TF 时的 D2(每套数据)。
2. `uniform`、`length_raw` 相对 `no_ppi`;`raw` 相对 `length_raw`(每套数据、两套新数据合并、四套合并)。
3. `ref_commot_B` 相对 `no_ppi`(在它自己有定义的受体上)。
4. 四套数据合并的 D1 / D2 / D3,标为描述性。

## 6. 执行(在 `isoform_commot/` 目录下)

1. 提交本文件和代码,记下提交号。
2. **第 1 步,数据与打分**(`reports/stage2_ppi_03/run_prepare_x.sh`,环境 commot-research + plminteract):
   - 冒烟:`bash reports/stage2_ppi_03/run_prepare_x.sh --smoke --runs analysis_x_prep_smoke --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D`
     - 通过条件:退出码 0;五个输入目录都有 `preparation.json`;`analysis_x_prep_smoke/scores/data_report.json` 已生成;
       `stage2_input_xb_*/preparation.json` 中 `paper_input_x_max_abs_diff` ≤ 1e-4。
   - 正式:去掉 `--smoke`,`--runs analysis_x_prep --results reports/stage2_ppi_03/scores`。提交 `reports/stage2_ppi_03/scores/`。
3. **第 2 步,实验**(`reports/stage2_ppi_03/run_x.sh`):
   - 冒烟:`bash reports/stage2_ppi_03/run_x.sh --smoke --runs analysis_x_smoke --jobs 8 --gpu-env inject-gpu --gpu-jobs 4`。
     冒烟内容:回归门槛、两个等价门槛、seqFISH+ 的 9 个冒烟条件与参照、seqFISH+ 的 0a / 0b、GPU 门槛、seqFISH+ 的 `ours_no_ppi`。结果不解读。
     - 通过条件:退出码 0;`analysis_x_smoke/summary/a_line/decisions.json`、`b_line/b_line.json`、`gates/` 下三个门槛文件都已生成且 passed。
     - 记下各步用时(`run.log`)。
   - 正式:`bash reports/stage2_ppi_03/run_x.sh --runs analysis_x_main --results reports/stage2_ppi_03/summary --jobs 8 --gpu-env inject-gpu --gpu-jobs 4`。
     中断后加 `--resume`。提交 `reports/stage2_ppi_03/summary/`。

**计算量:**
- A 线:4 套数据 × 205 次求解,CPU。本地每次约 6 s(seqFISH+)到 11 s(乳腺癌),单核。
- 参照:4 次 COMMOT 纯 collective,CPU。
- B 线:3 次 COMMOT 混合重跑(CPU;本地 seqFISH+ 约 4 分钟),12 次 GPU 求解(后部截断 726 µm,问题最大)。
- 第 1 步:4 组 V12 打分。

## 7. 已知局限

1. 评价仍是靶基因表达的间接证据,分辨率约 0.01(实验二)。
2. 论文的指标(B 线)不控制受体表达。在本地 seqFISH+ 上,`receptor_only` 的中位数(0.2305)与 COMMOT 原版(0.2301)几乎相同(见 §8),
   所以 B 线的数值主要反映受体表达,不能用来区分不同 S。
3. seqFISH+ 的"接触"邻居是按距离定义的(§3.2),不是细胞膜接触的测量;FOV 之间按论文的拼接坐标计算距离。
4. V12 对输入顺序敏感,分数取两个顺序的平均(与实验二相同)。
5. 人类 REACTOME 表与论文版本不同(§1),A 线乳腺癌的靶基因集合因此与论文不同;B 线用论文版本。

## 8. 写计划时已经看过的东西(透明记录)

开发和本地测试在提交前完成,用的是**随机替代分数**,不是 V12 分数。没有任何真实 S 的求解或评价。看过的内容:

- **A 线:** 输入规模。
  - seqFISH+:791 对(扩散 607、接触 184)、190 个配体、224 个受体,每次求解约 528 万个元素。
  - 乳腺癌:523 对(扩散 309、接触 214)、216 个配体、245 个受体。
  - 不加 S 时可评价的受体数:乳腺癌 164,seqFISH+ 174;`no_ppi` 的 M 为 0.105 和 0.043。这两个数与 S 无关。
  - 等价门槛在乳腺癌上通过(差 0.0)。
- **B 线(seqFISH+):**
  - 0a 通过:519 对,均值 0.23117,中位数 0.23007。
  - 0b 通过:接收信号与论文逐值一致,相对 L1 最大 1.4e-8。
  - `receptor_only`:均值 0.2317,中位数 0.2305。
  - `ours_no_ppi`(CPU,50,000 次迭代):均值 0.2306,中位数 0.2279。求解到达迭代上限,检查通过,与实验二次设置的情况相同(到上限只记录)。
    每次约 13–15 分钟(CPU,单核);正式运行在 GPU 上。
- **单元测试:** `test_exp2x` 11 项通过,其中包括在 `test_stage2_s` 的 toy 输入上与 `stage2_s.py run` 结果逐值相同。
