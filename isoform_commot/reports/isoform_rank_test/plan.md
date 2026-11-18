# Isoform 排序测试：计划（打分前提交）

2026-11-18。本计划在对测试集打分之前提交；提交后不再修改方法和指标。
代码：`isoform_rank_test.py`（`fetch`、`build`、`evaluate`）。打分：`plm_score.py`，自 `f132e4f` 起未改动。

## 1. 范围

- **输入：** 用户构建的 `reports/isoform_rank_test/testset.csv`，与本计划在同一个提交中（SHA256 `c35e1167…`，完整值记录在 `build.json`）。
  - 共 26 行、13 组，每组是同一配体对同一受体基因的 IIIb 和 IIIc 两行：
    - directional 10 组：expected 为 bind + nobind；
    - control 2 组：FGF1 × FGFR2、FGF1 × FGFR3；
    - weak 1 组：FGF2 × FGFR2，即 IIIc 为 bind、IIIb 为 weak。
  - FGFR1 的行在本计划之前已由用户删除。
  - `expected`、`confidence`、`rule`、`citation` 都按用户给出的使用，这里不核对、不修改。
- **性质：** 本测试只作描述。用户没有给出通过阈值，所以不作"可用"判定。如果需要判定，必须在打分前另行登记。
- 不涉及 OT 运行。

## 2. 序列（已在打分前构建）

- **来源：** UniProt 2026_04，2026-11-18 05:41 UTC 获取。
  - 原始文件放在 `data/isoform_rank_test/uniprot/`（不提交），每个文件的 SHA256 记在 `fetch.json`，同时复制进 `build.json`。
  - 配体：条目的 canonical 序列（条目 JSON）。
  - 受体：用该 isoform 自己的 ID 取 FASTA，例如 `P21802-3.fasta`，而不是条目的 canonical。
- **自动检查**（`build` 中实现，全部通过）：
  - isoform ID 属于该行的 `receptor_acc`；
  - 该 ID 列在受体条目的 ALTERNATIVE PRODUCTS 中；
  - 状态为 Displayed 的 isoform 序列等于 canonical 序列，状态为 Described 的不等于；
  - 同一受体的 IIIb 和 IIIc 序列不相同。
- 逐条序列的条目版本、序列版本、长度和 SHA256 见 `sequences.csv`。

受体（`receptor_check.csv`）。UniProt 名称和同义名照原样列出，与 testset 的 IIIb/IIIc 标签一致：

| 受体 | IIIb：ID，长度，UniProt 名称（同义名） | IIIc：ID，长度，UniProt 名称（同义名） | 长度差 | 两序列不同的区段（IIIb / IIIc） |
|---|---|---|---|---|
| FGFR2 | P21802-3，822 aa，isoform 3（BFR-1; FGFR2IIIb; KGFR），Described | P21802-1，821 aa，isoform 1（BEK; FGFR2IIIc），Displayed，即 canonical | 1 | 314–362 / 314–361 |
| FGFR3 | P22607-2，808 aa，isoform 2（IIIb），Described | P22607-1，806 aa，isoform 1（IIIc），Displayed，即 canonical | 2 | 311–360 / 311–358 |

两条序列在区段之外完全相同：FGFR2 前 313 位和后 460 位相同，FGFR3 前 310 位和后 448 位相同。

配体（canonical，长度单位 aa）：

| 配体 | 长度 | 配体 | 长度 |
|---|---|---|---|
| FGF1（P05230） | 155 | FGF9（P31371） | 208 |
| FGF2（P09038） | 288 | FGF10（O15520） | 208 |
| FGF3（P11487） | 239 | FGF17（O60258） | 216 |
| FGF7（P21781） | 194 | FGF18（O76093） | 207 |
| FGF8（P55075） | 233 | FGF22（Q9HCT0） | 170 |

- FGF2 的 UniProt canonical 是 288 aa 的序列，按要求照样使用。
- 四条受体序列都是 806–822 aa，比 humanV11 训练数据中最长的蛋白（800 aa）还长；humanV12 的训练数据含最长 1,996 aa 的蛋白（LR 评测报告 11.2 节）。
- 最长的一对是 1,113 token，没有超过官方的 1,603 上限，所以不存在截断问题。仍然照评测的做法传 `--allow-long`，对结果没有影响。

## 3. 打分

- `plm_score.py`，在 `plminteract` 环境中运行，基础模型是 ESM-2 650M。
- checkpoint 与 LR 评测相同：humanV11 `e86e392…`，humanV12 `2d0e7fb…`。
- `jobs.csv`：26 个不同的（配体，isoform）配对，每个配对两个输入顺序，共 52 个作业。
- 不截断（逐对核对 token 数）。每行分数 S 取两个顺序的平均。
- 输出：`reports/isoform_rank_test/scores/{humanV11,humanV12}.csv`，每个都带 `.meta.json`，记录输入和 checkpoint 的 SHA256。
- **冒烟测试（已做，在本计划提交之前）：** 用 `plm_score.py` 重打 LR 评测中的两个作业（`P0DMC3__Q8IVK1` 的两个顺序）。
  - humanV11 和 humanV12 的分数都与评测中的原值完全相同，最大差值为 0。
  - 结果在 `smoke/`。
  - 本计划提交前，没有对测试集中的任何配对打分。

## 4. 归一化与指标（固定）

1. **归一化：** S'(行) = S(行) / [S(IIIb 行) + S(IIIc 行)]，在组内计算。
   - 若组内两行分数之和为 0，S' 无定义。这样的 directional 组记为未命中，并单独注明。
   - 注意：每组只有两行，所以 S'(bind) > S'(nobind) 与 S(bind) > S(nobind) 等价。归一化不改变命中与否，只改变差值的尺度：S'(bind) − S'(nobind) = (S_bind − S_nobind) / (S_bind + S_nobind)，范围是 [−1, 1]。
2. **directional 组：**
   - S'(bind) > S'(nobind) 严格成立才算命中，相等算未命中；
   - 报告命中数 / 组数、命中率，以及每组的 S'(bind) − S'(nobind)；
   - 随机基线为 50%。
3. **control 组：**
   - 报告 |S'(IIIb) − S'(IIIc)|，预期较小；
   - 没有给定阈值，所以不作判定，也不计入命中率。
4. **weak 组：** 只报告两行的 S 和 S'，不作判定。
5. **高置信子集：** `confidence = high` 的 directional 组单独再报一次命中率。共 4 组：FGF7、FGF10、FGF8 各 × FGFR2，以及 FGF8 × FGFR3。
6. 以上全部对 humanV11 和 humanV12 各报一份。
7. **补充的描述量**（我加的，现在固定，不作判定）：
   - **精确二项尾概率** P(X ≥ k)，X ~ Binomial(n, 0.5)。n 只有 10 和 4，而且各组并不独立（见第 5 节），只作参考。
   - **固定偏向的参照命中率。** 如果模型不论配体是谁，总把同一个 isoform 打得更高，命中率就等于预期该 isoform 较高的组所占的比例：
     - 全部 10 组中，预期 IIIc 较高的有 6 组，所以"总是 IIIc 较高"得 60%，"总是 IIIb 较高"得 40%；
     - 高置信 4 组中两者各 2 组，所以都是 50%。
   - 同时按受体统计 S'(IIIb) > S'(IIIc) 的组数（`receptor_sign.csv`）。如果同一受体的所有组方向都相同，命中率反映的是受体 isoform 本身的偏向，而不是配体特异性。

## 5. 预先写明的局限

- **组数少，而且组与组不独立：**
  - 8 个 FGFR2 组用的是同一对受体序列；
  - FGF8、FGF9、FGF1 各出现在两个组中。
- **差异区段占比小：** IIIb 和 IIIc 只在约 50 个残基上不同，约占全长 806–822 aa 的 6%。
- **长度：** 受体序列长于 humanV11 的全部训练蛋白。LR 评测中，humanV11 在含 > 800 aa 蛋白的配对上区分能力明显较弱（评测报告 1.1 节）。
- **没有检查训练集重叠**（不在本次要求之内）。
- **标签来源：** 预期标签来自用户的 testset，本测试不核对。

## 6. 输出与流程

1. 提交本计划，以及 `isoform_rank_test.py`、`testset.csv`、`sequences.csv`、`receptor_check.csv`、`jobs.csv`、`build.json`、`smoke/`。
2. 正式打分：humanV11 和 humanV12 各 52 个作业。
3. `isoform_rank_test.py evaluate` 输出 `results/`：`rows.csv`、`groups.csv`、`summary.csv`、`receptor_sign.csv`、`tables.md`、`manifest.json`。
4. 写 `results.md`，内容包括：
   - 每组的分数和 S'；
   - 命中率（全部组 / 高置信组，humanV11 和 humanV12 各一份），连同 directional 组数和随机基线 50%；
   - 第 4.7 条的补充描述量。
5. 提交。推送由用户执行。
