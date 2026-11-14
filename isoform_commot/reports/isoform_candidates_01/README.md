# isoform 候选清单 01（后部切片面板，只盘点）

由 `python isoform_inventory.py --input stage2_input_post --output reports/isoform_candidates_01` 生成
（代码 `4fd2422`，工作区无修改；UniProt release 2026_04；见 `inventory.json`）。只做计数和标记，不做任何分析。

## 规则

- **基因：** 后部切片 OT 面板（540 对）中的全部配体和受体基因，复合物拆成亚基，共 358 个。
- **UniProt 条目：** 小鼠（10090）、已审阅（Swiss-Prot）、主基因名等于该基因的条目（与 `proteins.py` 相同的查询；不用同义名）。
- **isoform：** 条目 ALTERNATIVE PRODUCTS 中列出的全部 isoform（canonical 加 alternative products）；没有这项注释的计为 1。
- **可变区：** `Alternative sequence` 特征（canonical 坐标）。
- **标记**（可变区与以下区域有重叠即标记）：
  - 胞外区：`Topological domain` 中描述为 Extracellular 的区段；没有拓扑注释、亚细胞定位含 Secreted 的条目，取信号肽之后的整条链（没有信号肽注释则取整条链）；
  - 已注释的结构域：`Domain` 特征；
  - 已注释的结合区：描述含 bind / interact 的 `Region` 特征，以及 `Binding site` 特征。
- **候选：** 某个条目有 ≥2 个 isoform，且至少一个可变区带有上述任一标记。NRXN1 α/β 作为已知例子放在 `candidates.csv` 第一行。

## 计数

| 项目 | 数量 |
|---|---|
| 面板基因 | 358（配体 151、受体 162、两者都是 45） |
| 有已审阅小鼠条目 | 347 |
| 没有符合规则的条目 | 11：Flt3l、Kitl、Klk6、Ncl、Pcdhb7/9/12/17/20/22、Pvr（主基因名不一致或只有未审阅条目，未用同义名查找） |
| 有多个已审阅条目的基因 | 3：Nrxn1（Q9CS84 α、P0DI97 β）、Nrxn2（E9PUN2、E9Q7X7）、Nrxn3（Q8C985、Q6P9K9） |
| ≥2 个 isoform 的基因 | 120 |
| 其中可变区落在胞外区 / 结构域 / 结合区 | 88 / 54 / 23 |
| **候选基因** | **96**（受体 58、配体 29、两者都是 9；参与的配对：Contact 49 个基因、Secreted 38、两类都有 9） |

**NRXN1：** α 条目 Q9CS84（1514 aa）有 6 个可变区，落在胞外区、结构域和结合区；β 条目 P0DI97（468 aa）没有拓扑注释
（胞外区标记为 `none_annotated`）。两个条目的 ALTERNATIVE PRODUCTS 都列出 7 个 isoform。α/β 由不同启动子产生，
在 UniProt 中是两个条目，不表现为 `Alternative sequence`，所以按上面的规则不会被自动识别，单独列出。
Nrxn2、Nrxn3 也是 α/β 分为两个条目。

## 文件

| 文件 | 内容 |
|---|---|
| `candidates.csv` | 候选基因（Nrxn1 第一行）：条目、isoform 数、可变区数、三项标记、角色、配对类型、参与的 LR 对 |
| `genes.csv` | 全部 358 个基因的同样信息和 `candidate` 列 |
| `entries.csv` | 每个 UniProt 条目：长度、isoform 名称、胞外区的依据和区段、条目版本、序列 SHA256 |
| `variable_regions.csv` | 每个可变区：位置、描述、三项标记、重叠的结构域和结合区名称 |
| `inventory.json` | 输入哈希、UniProt release、查询规则、每个条目 JSON 的 SHA256、启动时的 git 状态 |
