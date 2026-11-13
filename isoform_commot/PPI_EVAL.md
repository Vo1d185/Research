# 第二阶段：接入 PLM-interact 分数，并用靶基因评价

本阶段在**完全不变的输入**上（`stage1_input_01`：128 个 spot、8 个 LR 对、r0=2h、ε=0.2、ρ=1）完成四件事：

1. 小鼠基因 → 小鼠蛋白序列 → 官方 PLM-interact 打分 → `ppi.csv`
2. 在 S=1、真实分数、同类型均值常数、同类型内全部排列之间比较 OT 结果
3. 在同一输入上运行 COMMOT 原版（作为基线）
4. 靶基因评价（COMMOT 论文用的评价方式）：接收信号是否和下游靶基因表达一致

所有命令都在 SSH 服务器上，在 `Research/isoform_commot` 目录下执行。每个输出目录都必须是新的，程序不会覆盖已有目录。

## 0. 更新代码

新代码在分支 `claude/exciting-faraday-oxza82`，它是 `research/isoform-commot-v1` 加上新的提交。

```bash
cd /home/liu/Research
```

```bash
git fetch origin
```

```bash
git switch research/isoform-commot-v1
```

```bash
git merge --ff-only origin/claude/exciting-faraday-oxza82
```

```bash
conda activate commot-research
```

```bash
cd isoform_commot
```

```bash
python -m pip install -r requirements-eval.txt
```

```bash
python fetch_solver.py
```

```bash
python -m unittest
```

应显示 25 个测试全部通过。`fetch_solver.py` 现在还会下载 COMMOT 原版的 `_cot.py`（固定提交，并校验 SHA256）。

## 1. 小鼠基因 → 小鼠蛋白（commot 环境）

```bash
python proteins.py map --input stage1_input_01 --output proteins_01
```

- 程序从 h5ad 找回原始小鼠符号（例如 `APOE` → `Apoe`），再到 UniProt 查询**小鼠、已审核（Swiss-Prot）、主基因名完全一致**的条目。
- 每个基因必须恰好对应 1 个条目。程序会先检查完所有基因，把有歧义或找不到条目的基因及其候选（accession、蛋白名、长度）一次性写到 `proteins_01_candidates.csv`，然后停止，不生成输出目录。你逐个确认后写一个 `override.csv`（列为 `lr_symbol,accession,reason`，`lr_symbol` 用 `lr.csv` 里的大写名称，`reason` 必填），再加上 `--override override.csv` 重新运行。
- 一个基因对应多个已审核条目，通常是同一基因的不同蛋白形式，例如 `Nrxn1` 的 α 型（Q9CS84）和 β 型（P0DI97）。gene-level 数据分不出它们，选哪个必须写明理由。没选的那个保留下来，以后做 isoform 分析时用。
- **请逐行检查 `proteins_01/protein_map.csv`**（最多 16 个蛋白）：accession、蛋白名、长度是否合理。`longer_than_training_max=True` 表示该蛋白超过 800 aa，超出了官方训练和测试数据的长度范围。
- 用的是 UniProt canonical 全长序列（包含信号肽），和训练数据的形式一致。原始 JSON 保存在 `proteins_01/uniprot/`，以后做胞外段敏感性分析时可以直接用。

## 2. 生成打分任务

```bash
python proteins.py pairs --input stage1_input_01 --proteins proteins_01 --lr-dir data/cccvelo_lr --output plm_jobs_01
```

- 已知配对：8 对 × 2 个方向（L→R 和 R→L）。
- 面板内的非配对组合（ligand × receptor 中不在 `lr.csv` 里的），同样 2 个方向，**只用来检查分数，不进入 OT**。`--lr-dir` 会标出其中哪些其实在完整的 LR 库里有记录。
- 查看 `plm_jobs_01/jobs.json` 中的 `max_tokens` 和 `jobs_over_official_1603_tokens`。

## 3. PLM-interact 打分（单独的 Python 3.10 环境）

建立环境（只需做一次，**不要**装进 commot-research 环境）：

```bash
conda create -n plm-interact python=3.10 -y
```

```bash
conda activate plm-interact
```

先按服务器的 CUDA 版本安装 torch，例如 CUDA 11.8：

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

```bash
pip install -r requirements-plm.txt
```

查询 checkpoint 当前的提交号，**记下这个值**：

```bash
python -c "from huggingface_hub import HfApi; print(HfApi().model_info('danliu1226/PLM-interact-650M-humanV11').sha)"
```

用这个提交号下载权重和对应的 ESM-2 基础模型（把 `<sha>` 换成上一步的输出）：

```bash
huggingface-cli download danliu1226/PLM-interact-650M-humanV11 --revision <sha> --local-dir models/plm_650m_humanV11
```

```bash
huggingface-cli download facebook/esm2_t33_650M_UR50D --local-dir models/esm2_t33_650M_UR50D
```

如果提示 `huggingface-cli` 已弃用或找不到这个命令，把它换成 `hf download`，其余参数不变。

**先做安装自检。** 用官方的小鼠跨物种测试集，确认环境和权重加载都没问题：

```bash
git clone https://github.com/liudan111/PLM-interact.git models/PLM-interact-src
```

```bash
git -C models/PLM-interact-src checkout ab2ae6aae1aa81accf4c1f7f6f341164baf97809
```

```bash
python plm_score.py selftest --test-file models/PLM-interact-src/data/test/mouse.ppi.qrels.seq.test.csv --checkpoint models/plm_650m_humanV11/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D --checkpoint-name danliu1226/PLM-interact-650M-humanV11 --checkpoint-revision <sha> --output plm_selftest_01.json
```

把得到的 AUPR 和论文里小鼠跨物种的结果对比。差距很大，通常说明基础模型或 checkpoint 不匹配。默认抽 2000 对，没有 GPU 时可以加 `--n 300`。

**正式打分：**

```bash
python plm_score.py score --jobs plm_jobs_01/plm_input.csv --checkpoint models/plm_650m_humanV11/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D --checkpoint-name danliu1226/PLM-interact-650M-humanV11 --checkpoint-revision <sha> --output plm_scores_01.csv
```

- 每次只算一对，`max_length` 按这一对的实际长度设置，**禁止截断**，并逐对核对 token 数（官方批量脚本会静默截断）。
- 如果有蛋白对超过 1603 个 token，程序会拒绝执行并列出这些对。这时需要先决定规则，比如受体只用胞外段，或者明确记录"按全长打分"（加 `--allow-long`）。**不要看了 OT 结果再决定。**
- 稳健性对照：35M 模型用 `facebook/esm2_t12_35M_UR50D` 作为基础模型，另外也可以用 `PLM-interact-650M-humanV12`。只报告结果是否一致，不从中挑最好的。

## 4. 导入分数（回到 commot 环境）

```bash
conda activate commot-research
```

```bash
python proteins.py import --input stage1_input_01 --jobs plm_jobs_01 --scores plm_scores_01.csv --output ppi_01
```

**进 OT 之前，先看 `ppi_01/score_report.json`：**

| 字段 | 看什么 |
|---|---|
| `cognate` | 8 个分数的范围。如果全部接近 1，PPI 项约等于 0，OT 结果应该和 S=1 几乎一样 |
| `order_asymmetry` | 两个方向的差异。差异大说明分数不稳定 |
| `cognate_vs_unannotated_cross_auroc` | 已知配对是否比未注释的组合得分更高。接近 0.5 说明模型对这些小鼠蛋白区分不出来 |

`ppi.csv` 中的分数是两个方向的平均值（这是事先定好的规则）。

## 5. 加入 PPI 并运行 OT

```bash
python stage1.py attach-ppi --input stage1_input_01 --ppi-dir ppi_01 --output stage1_input_01_ppi
```

```bash
python stage1.py run --input stage1_input_01_ppi --output stage1_results_ppi_01
```

- `attach-ppi` 原样复制 `stage1_input_01` 并加上 `ppi.csv`，同时**核对表达、坐标和 LR 文件的哈希都没变**。
- 默认只跑 r0=2h（`--r0-h 2`）。比较的条件包括：`no_ppi`（S=1）、`ppi`（真实分数）、`ppi_constant`（每类都用该类分数的均值）、`perm_diffusion_*` 和 `perm_contact_*`（一类内部的全部排列，另一类保持真实分数）。
- 在 `stage1_results_ppi_01/ot_00/analysis/` 里看：
  - `control_matched.csv`：每种条件相对 S=1 的变化。`between_pair_TV` 是各 LR 对之间份额的变化，`within_pair_TV` 是各对内部空间分配的变化。只做整体缩放时，后者为 0。
  - `shared_entity_split.csv`：被多个 LR 对共享的配体或受体，它的通信量在这些对之间怎么分。如果这个表是空的，说明面板里没有共享实体，竞争效应在这个面板上测不到。
  - `lr_pair_changes.png`：现在画的是 r0_raw、2h 这个条件。

## 6. COMMOT 原版基线

```bash
python commot_baseline.py --input stage1_input_01 --output commot_results_01
```

```bash
python commot_baseline.py --input stage1_input_01 --output commot_results_01_pure --weights 1 0 0 0
```

- 调用原版 `cot_combine_sparse`，默认参数为 ε=0.1、ρ=10，`dis_thr` 取 3h、5h、10h（**5h 是主要设置**）。
- **已发现的情况：** COMMOT 默认的四方案混合中，有三个方案是分开求解的，每个都会重复使用预算，所以单个配体或受体的总运输量**可能超过它的表达预算**。超出量记录在 `summary.csv` 的 `max_*_budget_excess` 列。`--weights 1 0 0 0` 只用 collective 方案，严格满足预算约束。两者都保留，写论文时要说明。
- COMMOT 不区分 contact 和 diffusion，并且允许同一个 spot 内部的运输。这些都是原版的行为，不做修改。

## 7. 靶基因评价

导出同一批 spot 的全部基因表达。程序会重新计算 LR 基因的表达并和 `expression.csv` 逐值核对：

```bash
python stage1.py export-genes --input stage1_input_01 --output genes_01
```

下载并转换 scSeqComm 的小鼠靶基因库（固定提交，并校验 SHA256）：

```bash
python evaluate_targets.py fetch-db
```

运行评价。`--covariate` 填 `cell_meta.csv` 里的细胞类型列名，先打开那个文件确认：

```bash
python evaluate_targets.py run --input stage1_input_01 --genes genes_01 --results stage1_results_01/ot_00/results stage1_results_ppi_01/ot_00/results commot_results_01 commot_results_01_pure --covariate scc_anno --output target_eval_01
```

事先定好的规则（详见 `evaluate_targets.py` 开头的说明）：

- 受体 → TF：在 scSeqComm 小鼠 KEGG 和 Reactome 中取 PPR 分数最高的前 10 个 TF（主分析），前 5 个和前 20 个作为敏感性分析。
- 靶基因：这些 TF 的靶基因并集，去掉面板里的全部 LR 基因，只保留在 ≥5% 的 ROI spot 中有表达的基因，每个受体至少需要 5 个靶基因。
- 靶基因得分：平均表达（和 COMMOT 一致）；z-score 平均作为敏感性分析。
- 指标：
  - Spearman 相关；
  - 只用受体表达时的 Spearman（对照）；
  - 控制受体表达（以及细胞类型）后的偏 Spearman。
  - 不报告 p 值，因为 spot 之间存在空间自相关。

输出文件：

| 文件 | 内容 |
|---|---|
| `target_sets.csv` | 每个受体用了哪些 TF、多少个靶基因，或者为什么没有被评价 |
| `method_summary.csv` | S=1 时，各距离设计、inverse 和 COMMOT 原版在靶基因指标上的比较（这一项也能为 r0 的选择提供证据） |
| `ppi_comparison.csv` | 真实分数和各对照相对 S=1 的变化。`rank_corr_vs_no_ppi = 1` 表示该单元只是被整体缩放，Spearman 不可能改变 |
| `permutation_summary.csv` | 真实分数在全部排列中的位置，`p_ge` 越小越好。只有 4 对时，最小只能到 1/24 左右 |
| `target_eval.png` | 左图：各方法的比较；右图：真实分数和各对照的比较 |

## 结果应该怎样解读

- **真实分数只有在明显好于打乱分数和均值常数时**，才能说"具体是哪一对拿到哪个分数"起了作用。
- 偏 Spearman 不高于"只用受体表达"的对照时，说明推断出的通信没有提供受体表达以外的信息。
- 128 个 spot、8 对 LR 只够检查流程是否正确，不够得出统计结论。正式评价需要扩大区域（稀疏求解）、增加 LR 对，并在独立 ROI 上重复。
- COMMOT 论文从 LR 对到靶基因的具体规则写在它的 Zenodo 复现代码里（10.5281/zenodo.7272562），我们还没有核对。下载下来以后，要和"前 10 个 TF"这条规则对齐。

## 在云端验证过和没有验证过的部分

- 已验证：25 个单元测试和集成测试（合成数据，用假的 UniProt 和假的 PLM 模型）；scSeqComm 真实文件的下载、校验和转换（857 个 TF、731 个受体）；COMMOT 原版求解器在示例数据上的运行和预算行为。
- 没有验证：UniProt、HuggingFace、Zenodo 和 h5ad 都在云端无法访问，所以真实蛋白映射、真实模型推理、真实数据上的评价都要在服务器上第一次运行。
