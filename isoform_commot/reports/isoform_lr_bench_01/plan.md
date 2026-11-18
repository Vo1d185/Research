# isoform 排序基准 01：PLM-interact humanV12 能否区分同一基因的不同 isoform

**状态：预注册计划。** 2026-11-18 起草。本文件提交（commit）之后才运行正式评估；
`run_main.sh` 会检查本文件已提交且工作区无改动，并把提交号写入结果。
正式运行之后对本文件的任何修改、以及第 6 节之外的分析，一律标注为“探索性”。

代码：`isoform_lr_bench.py`（jobs / merge / evaluate / selftest），`test_isoform_lr_bench.py`，
`reports/isoform_lr_bench_01/smoke_test.sh`，`reports/isoform_lr_bench_01/run_main.sh`。
打分沿用 `plm_score.py`，不改动。

## 0. 要回答的问题

固定一个配对蛋白（partner），把同一组候选序列（variant）分别和它配对打分，
看 PLM 分数能否把“有互作记录的那个 variant”排在“没有记录的 variant”前面。两层：

| 层 | 候选序列 | 与参照的序列一致性（中位数） | 作用 |
|---|---|---|---|
| **Tier-H（主）** | 同一基因的不同剪接 isoform（UniProt `-1/-2/...`） | 80.2% | 要回答的问题本身，难 |
| **Tier-E（对照）** | 同一 UniProt 家族的不同基因（旁系同源） | 19.6% | 检验这种排序设定在容易的情况下是否有效 |

与已有工作的关系：
- `reports/isoform_feasibility.md` 第 0 节导师建议的归一化 S'(i,j) = S(i,j) / Σ S(i',j')，前提是 PLM 分数在同一基因对的 isoform 之间有可靠差别。本基准直接检验这个前提。
- `reports/isoform_plan_draft.md` 用结构域规则判断 isoform 兼容性。若本基准判定 isoform 层不可用，isoform 部分继续只走结构域规则路线。
- 本基准是人类序列（与 humanV12 的训练物种一致）；OT 分析用小鼠切片，这一点不变。

## 1. 数据（已冻结）

文件 `dataset/benchmark.csv.gz`，**sha256 `d3738f564f48fa5792e7e0ee0d2e7f2e5246407a83e754a975cdca3f84bedb3b`**，
1,817 行。来源、查询式、比对参数和每一步的计数见 `dataset/provenance.json`。代码中写死了这个 sha256，文件不符时拒绝运行
（除非加 `--exploratory`）。

- **范围：** UniProt release 2026_04（2026-11-04），人类、已审阅条目。配体 = GO:0048018 且定位为分泌或膜；受体 = GO:0004888。
- **证据（“有记录”）：** UniProt INTERACTION 注释中带 isoform 后缀的互作记录（如 `P21802-3` 与某配对蛋白）。
- **Tier-H 对照：** 每个（基因 × 配对蛋白）一个对照。参照 = 有记录的 isoform（有多个时取实验数最多的）；
  竞争者 = 同一基因中没有记录的其他 isoform。若某个“没有记录”的 isoform 在 IntAct 中有针对同一配对蛋白的
  isoform 层记录，则移出竞争者（`excluded_other_evidence`）。
- **Tier-E 对照：** 每个（家族 × 方向 × 配对蛋白）一个对照。参照 = 基因层有记录的成员；竞争者 = 同家族其他成员。
  若某成员与该配对蛋白在 CellChatDB v2（`reports/plm_lr_benchmark/dataset/positives_all_human.csv`，含复合物亚基）中是已知配对，
  则移出竞争者（`excluded_known_lr`）。
- **长度：** 每个蛋白 24–1996 aa（humanV12 训练范围）。
- 每个对照恰有 1 个参照、≥ 1 个竞争者。

| | Tier-H | Tier-E |
|---|---|---|
| 对照数 | 146 | 154 |
| 参照基因数 | 73 | 103 |
| 参与排序的行（参照 + 竞争者） | 146 + 666 = 812 | 154 + 731 = 885 |
| 不参与排序的行 | 35（另有证据） | 51（CellChatDB 已知）+ 34（另有证据） |
| 每个对照的候选数 k | 中位 4（2–19） | 中位 5（2–12） |
| 随机排序的 hit@1 期望 = mean(1/k) | 0.265 | 0.249 |

Tier-H 参照中，canonical isoform 49 个、非 canonical 97 个。竞争者与参照的一致性分布：
[0,50) 132，[50,70) 105，[70,90) 207，[90,99) 194，[99,100] 28。

## 2. 标签的含义和局限（写在前面）

“没有记录”**不等于**“不结合”。没有可用的 isoform 层阴性数据（本次查询的 IntAct 记录中没有标为阴性的）。因此：
- 本基准测的是“模型是否把有记录的 isoform 排在前面”，不是二分类准确率；不报告 AUROC。
- 若竞争者其实也结合，模型把它排在前面会被算错。这会让结果偏向“区分不出”。所以判定为 FAIL 时，结论是
  “在现有证据下没有显示出 isoform 层的排序能力”，不是“模型一定区分不了”。
- 记录本身偏向研究得多的 isoform，可能和“canonical”“最长”这类与互作无关的属性相关。所以设置第 4.4 节的平凡基线，
  并要求模型显著超过它们（C3）。

## 3. 打分

- 模型：`danliu1226/PLM-interact-650M-humanV12` @ `2d0e7fb3deef94efa3f217daccc79092a6a42ceb`，
  checkpoint sha256 `92565937d65a9b91a4e897c99653c972b99cfe37067eb42d9e783af688350547`
  （与 `reports/plm_lr_benchmark/reference/selftest_humanV12_mouse_test.json` 相同的文件）；
  底座 `models/esm2_t33_650M_UR50D`。代码核对 checkpoint sha256，不符时拒绝给出预注册结果。
- 工具：`plm_score.py score`，逐对、不截断、核对 token 数，与 `plm_lr_benchmark_plan.md` 第 3.2 节相同。
- 两个顺序（variant, partner）和（partner, variant）都打分，**分数 = 两个概率的平均**。共 3,634 个作业。
- 全长打分（`--allow-long`）：312 行超过 1603 token（其中参与排序的 286 行），最长 2790 token。按 token 分层只做描述（第 6 节）。

## 4. 指标

对每个对照：把 k 个候选按分数从高到低排序（并列取平均名次），参照的名次记为 rank。

1. **hit@1：** rank = 1 的对照比例（并列第一不算命中）。
2. **平均归一化名次：** (rank − 1)/(k − 1) 的平均；0 = 总排第一，随机期望 0.5。
3. **置换检验 p 值（单侧）：** 零假设下参照在 k 个名次中等概率出现、各对照独立；10,000 次模拟，p = (1 + #{零分布 ≥ 观测}) / (1 + 10,000)，
   随机种子 20261118。
4. **平凡基线（只 Tier-H）：** (a) 总选 canonical isoform；(b) 总选最长的 isoform。用同一套 hit@1 计算。
   在冻结数据上：canonical 0.336，最长 0.384。
5. **按基因的整群 bootstrap 95% 区间：** 以参照基因为单位重抽（同一基因的多个对照不独立，如 BAG6 有 10 个），2,000 次。只用于报告，不参与判定。
6. **训练集重叠：** 与 `lr_eval.py` 相同的规则（`lr_eval.training_pairs`，精确序列 SHA1，无序对，label = 1），
   训练文件 `data/lr_benchmark/plm_training/STRING_V12_TrainingSet/protein.pairs_human_V12_train.tsv`。
   对每个对照标记“参照对是否是训练正样本”，并在“参照对不在训练正样本中”的 Tier-H 子集上重算指标（C5）。

## 5. 判定标准（预注册）

| 编号 | 条件 | 数据 |
|---|---|---|
| C1 | hit@1 的 p < 0.05 **且** 平均归一化名次的 p < 0.05 | Tier-H 全部 |
| C2 | hit@1 ≥ 0.45 **且** 平均归一化名次 ≤ 0.40 | Tier-H 全部 |
| C3 | hit@1 − max(canonical 基线, 最长基线) ≥ 0.10（按冻结数据即 hit@1 ≥ 0.484） | Tier-H 全部 |
| C5 | 平均归一化名次的 p < 0.05；该子集 < 30 个对照时记为“无法评估”，不参与判定 | Tier-H，参照对不在训练正样本中 |
| C4 | hit@1 ≥ 0.50 **且** hit@1 的 p < 0.05 | Tier-E 全部 |

- **isoform 层判定** = PASS 当且仅当 C1、C2、C3 都成立，且 C5 成立或无法评估；否则 FAIL。没有运行训练集重叠时输出 “incomplete”，不算最终判定。
- **旁系层判定** = PASS 当且仅当 C4 成立。

| 旁系层 | isoform 层 | 解读 | 对后续的含义 |
|---|---|---|---|
| PASS | PASS | V12 在固定配对蛋白时能把有记录的 isoform 排在前面，且超过 canonical / 最长基线 | 可以在 isoform 层使用 S；归一化 S' 的前提成立；进入 isoform 层 OT 实验的设计 |
| PASS | FAIL | 只有基因层的排序能力 | isoform 层不用 PLM 分数，继续 `isoform_plan_draft.md` 的结构域规则；S' 归一化主要在重新分配噪声 |
| FAIL | 任意 | 这种排序设定或标签在容易的情况下都不起作用 | isoform 层结果不可解释；先改进证据来源或标签，再重做 |

阈值的理由：随机 hit@1 约 0.25–0.27，C2、C4 要求高出约 0.2；C3 要求比不看配对蛋白的规则至少高 0.10，
排除“模型只是偏好 canonical 或长序列”。**如需修改阈值，请在提交本文件之前修改**（同时修改 `isoform_lr_bench.py` 中的 `CRITERIA`）。

## 6. 探索性分析（不参与判定）

- Tier-H 竞争者按与参照的一致性分箱的成对胜率（[0,50)、[50,70)、[70,90)、[90,99)、[99,100]，随机 = 0.5）。
- 按 token 数（≤ 1603 / > 1603）分层。
- 单一顺序的分数；不参与排序的行（`excluded_*`）的分数分布。
- 同一流程用 humanV11、Leakage-Free 打分（`--model` 只改变训练集重叠所用的训练文件；需加 `--exploratory`，因为预注册只登记了 humanV12）。
- 按家族、按配体 / 受体方向分层。

## 7. 执行

在 `isoform_commot/` 目录下：

1. 提交本文件和代码，记下提交号。
2. 冒烟测试（第 9 节）：
   `bash reports/isoform_lr_bench_01/smoke_test.sh --out isoform_lr_bench_runs/smoke_01 --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D`
3. 正式运行：
   `bash reports/isoform_lr_bench_01/run_main.sh --runs isoform_lr_bench_runs/main_01 --results reports/isoform_lr_bench_01/results_humanV12 --checkpoint models/plm_650m_humanV12/pytorch_model.bin --base-model models/esm2_t33_650M_UR50D --training-root data/lr_benchmark/plm_training`
4. 提交 `reports/isoform_lr_bench_01/results_humanV12/`（metrics.json、verdict.json、per_contrast.csv、pairwise.csv、row_scores.csv、
   jobs.json、scores.meta.json、merged.meta.json、run_info.json）。`isoform_lr_bench_runs/` 不提交（已加入 .gitignore）。
5. 结果解读写在 `reports/isoform_lr_bench_01/results.md`，先报告第 5 节判定，再报告第 6 节探索性结果。

环境：`commot-research`（单元测试、jobs、merge、evaluate、训练集重叠）；`plminteract`（`plm_score.py`，GPU）。
预计 GPU 时间为十分钟量级（3,634 个作业，未实测）。

## 8. 已知局限

1. 阳性证据只来自 UniProt / IntAct 的 isoform 层记录，覆盖有限，并偏向研究得多的基因（BAG6、CLEC2D、APP 等占多个对照）。
2. 没有 isoform 层阴性数据（第 2 节）。
3. 数据集是交互式构建后冻结的，没有附构建脚本；冻结文件 + sha256 是预注册的输入。
4. 人类序列；对小鼠 OT 的适用性需要另外论证。
5. 旁系层的竞争者只用 CellChatDB v2 排除了已知配对，没有用 CellPhoneDB v5；可能有少量真实配对被当作竞争者，会让 Tier-E 偏难。

## 9. 冒烟测试的通过条件（不解读分数）

`smoke_test.sh` 依次运行：单元测试；在冻结数据上用合成分数做的判定规则自检（oracle → PASS，随机 → FAIL，
只认 canonical / 只认最长 → C3 不成立）；3 个对照 / 层 + 含最长一对（2790 token）的对照共 52 个作业的真实打分；merge；evaluate（探索性）。
通过条件：全部步骤退出码为 0，`smoke_summary.json` 中 `selftest_passed` 与 `checkpoint_is_registered` 为 true，
且打分行数 × 2 = 作业数。冒烟测试的分数不看、不报告。
