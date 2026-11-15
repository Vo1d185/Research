# isoform 候选清单 02（后部切片面板，只盘点）

由下面的命令生成（代码 `14307cf`，工作区无修改，2026-11-15 05:41:13 UTC 启动；UniProt release 2026_04，与清单 01 相同；见 `inventory.json`）：

```
python isoform_inventory.py --input stage2_input_post --output reports/isoform_candidates_02 \
    --previous reports/isoform_candidates_01
```

只做计数和标记，不做任何分析。与[清单 01](../isoform_candidates_01/README.md) 相比，只补充了"同一基因有多个已审阅条目"的规则，其余规则不变。

## 规则

- **基因：** 后部切片 OT 面板（540 对）中的全部配体和受体基因，复合物拆成亚基，共 358 个（同 01）。
- **基因的 UniProt 条目（新）：**
  1. 小鼠、已审阅、主基因名等于该基因的条目（同 01）；
  2. 再取这些条目的 MGI 交叉引用，查询交叉引用到同一 MGI 基因的全部小鼠已审阅条目（`xref:"mgi-MGI:…"`），一并归入该基因。
     这样主基因名不同、但属于同一个 MGI 基因的条目也会被找到。

  同一基因的多个条目**视为该基因的不同 isoform**。基因的 isoform 数 = 各条目 ALTERNATIVE PRODUCTS 中 isoform 标识的并集
  （`merged_isoforms`；没有这项注释的条目计为一个）。
- **可变区，两个来源：**
  - `alternative_sequence`：条目内的 `Alternative sequence` 特征（同 01）；
  - `between_entries`（新）：多条目基因中，以最长的条目为参照（长度相同按登录号），把其余每个条目与参照比较。
    两条序列中按顺序排列、完全相同且长度 ≥ 10 aa 的片段算作共有（`difflib.SequenceMatcher`，关闭 autojunk）；
    共有片段之间的参照区段记为可变区，类型为缺失（另一条目没有对应残基）或替换（另一条目在此处有自己的残基）。
    另一条目在两个共有片段之间多出、而参照没有对应残基的部分，记为插在参照该位置的插入。坐标用参照条目的 canonical 坐标。
- **标记**（同 01，在描述可变区的那个条目上判断）：胞外区、`Domain`、结合区（描述含 bind / interact 的 `Region`，以及 `Binding site`）。
- **候选：** (a) 同 01：某个条目有 ≥ 2 个 isoform，且至少一个 `alternative_sequence` 可变区带有标记；
  或 (b)（新）：基因有 ≥ 2 个已审阅条目，且至少一个 `between_entries` 可变区带有标记。
  NRXN1 α/β 作为已知例子仍放在 `candidates.csv` 第一行。

## 计数

| 项目 | 01 | 02 |
|---|---|---|
| 面板基因 / 有已审阅条目 | 358 / 347 | 358 / 347 |
| 没有符合规则的条目 | 11 | 11（同 01：Flt3l、Kitl、Klk6、Ncl、Pcdhb7/9/12/17/20/22、Pvr；没有条目就没有 MGI 号，新规则对它们不起作用） |
| 只通过 MGI 交叉引用找到的条目 | — | 0 |
| 有多个已审阅条目的基因 | 3 | 3：Nrxn1（Q9CS84、P0DI97）、Nrxn2（E9Q7X7、E9PUN2）、Nrxn3（Q6P9K9、Q8C985） |
| 按规则 (a) 的候选 | 96 | 96（与 01 完全相同） |
| 按规则 (b) 的候选 | — | 3：Nrxn1、Nrxn2、Nrxn3 |
| **候选基因** | **96** | **97**（新增 Nrxn2；受体 58、配体 30、两者都是 9；配对类型：Contact 50、Secreted 38、两类都有 9） |

与 01 的全部差别见 `changes_vs_previous.csv`：只有 Nrxn2 一行（条目相同，从非候选变为候选）。

## 多条目基因的条目间可变区

| 基因 | 参照（α） | 另一条目（β） | 参照中的可变区 | 共有 | 标记 |
|---|---|---|---|---|---|
| Nrxn1 | Q9CS84，1514 aa | P0DI97，468 aa | 1–1129 被 β 的 1–83 替换 | C 端 385 aa | 胞外、结构域（EGF-like 1–3、Laminin G-like 1–5）、结合 |
| Nrxn2 | E9Q7X7，1710 aa | E9PUN2，660 aa | 1–1135 被 β 的 1–85 替换 | C 端 575 aa | 同上 |
| Nrxn3 | Q6P9K9，1571 aa | Q8C985，567 aa | 1–1088 被 β 的 1–81 替换；β 在参照 1364/1365 之间多 3 aa | 483 aa | 替换区同上；插入点在胞外区，不在已注释结构域内 |

三个基因的 Laminin G-like 6（LNS6）都在共有部分，α 与 β 的 canonical 序列在 LNS6 及之后完全相同（Nrxn3 的 3 aa 插入除外）。

**注意：** 三个基因的"结合"标记都来自 `Binding site` 特征，这些特征的配体是 Ca²⁺（Q9CS84 有 10 个，E9Q7X7 有 10 个，分布在 LNS 结构域中），
不是与配对蛋白的结合区；`Binding site` 特征没有文字描述，所以 `variable_regions.csv` 的 `binding` 列为空。
这是清单 01 沿用下来的标记规则，本次没有改；isoform 方案中判断结合决定区时要把金属离子结合位点排除（见方案草稿第 2.1 节）。

## Nrxn2 在清单 01 中是怎样漏掉的

Nrxn2 的 α（E9Q7X7）和 β（E9PUN2）由不同启动子产生，在 UniProt 中是两个已审阅条目。
每个条目的 ALTERNATIVE PRODUCTS 都列出 2 个 isoform（1a、1b），但其中一个的序列在另一个条目里（外部序列），
所以两个条目都没有 `Alternative sequence` 特征。清单 01 的候选规则只在单个条目内寻找带标记的 `Alternative sequence`，
虽然列出了"有多个已审阅条目的基因"，但没有比较条目之间的序列，因此 Nrxn2 不是候选。
Nrxn1 和 Nrxn3 的 α 条目本身有 `Alternative sequence`（剪接位点等），所以按规则 (a) 已经是候选；Nrxn1 另外作为已知例子列在第一行。

## 文件

| 文件 | 内容 |
|---|---|
| `candidates.csv` | 候选基因（Nrxn1 第一行）：条目、参照条目、合并后的 isoform 数、两类可变区的数目和标记、满足的规则（`candidate_rule`）、角色、配对类型、参与的 LR 对 |
| `genes.csv` | 全部 358 个基因的同样信息和 `candidate`、`candidate_alternative_sequence`、`candidate_several_entries` 列 |
| `entries.csv` | 每个 UniProt 条目：找到的途径（`found_by`）、MGI 号、长度、isoform 名称、胞外区的依据和区段、条目版本、序列 SHA256 |
| `variable_regions.csv` | 每个可变区：来源（`source`）、所在条目和比较的条目、位置、描述、三项标记、重叠的结构域和结合区名称 |
| `changes_vs_previous.csv` | 与清单 01 相比候选状态或条目有变化的基因 |
| `inventory.json` | 输入哈希、UniProt release、查询规则、最小共有片段长度、每个条目 JSON 的 SHA256、启动时的 git 状态 |
