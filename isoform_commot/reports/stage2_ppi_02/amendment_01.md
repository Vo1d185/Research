# 实验二 修订 01:扩散截断固定为计划值 363.0757 µm(2026-11-20)

**写于任何带 S 的 OT 求解之前。** 原计划提交为 4e3874f。两次冒烟都在求解前停下(见下),没有产生任何评价结果。分析设计、判定规则和 β 均不变。

## 1. 发现的问题

- 计划第 2 节(plan.md 第 35 行)规定扩散截断为 363.08 µm,与 Stage 2 条件 A 相同。
- `stage2_s.py run` 直接读输入目录 `config.json` 里的 `diffusion_cutoff_um`,并检查它等于 363.0757。
- `stage2_input_post/config.json` 存的是 500.0(11-14 按当时的默认值准备)。Stage 2 跑后部时,是用 `--diffusion-cutoff-um 363.0757` 在命令行覆盖的(STAGE2.md 第 33 行)。
- `stage2_input_ant` 存的是 363.0757,不受影响。
- 冒烟第 3 步(等价门槛)因此报 `Unexpected diffusion cutoff`,停在求解之前。这个检查起了作用;错的是代码假定输入目录存的就是计划值。

## 2. 改动(只涉及截断从哪里来)

| 文件 | 改动 |
|---|---|
| `stage2_s.py` | 新增常量 `CUTOFF_UM = 363.0757` 和函数 `load_input`:读入输入后,若 config 存的截断与计划值不同,就按 363.0757 重建扩散连接(`stage2.spatial_edges`,与 Stage 2 相同);删掉原来的检查。每个条件的 `case.json` 记录 `cutoff_um`(实际使用)和 `config_cutoff_um`(config 里存的) |
| `gpu_unot.py` | `equivalence` 新增可选参数 `--diffusion-cutoff-um`;不给时行为不变(注入实验的用法不受影响)。结果 json 记录实际截断 |
| `run_s.sh` | GPU 等价门槛传入 `--diffusion-cutoff-um 363.0757`,让门槛检查的就是正式运行要解的问题 |
| `test_stage2_s.py` | 新增 `Cutoff` 测试:config 存 363.0757 和 500 时,扩散连接都与按 363.0757 构建的完全相同,`no_ppi`/`raw` 的评价结果逐值相同。把修复去掉后,这个测试会失败(已验证) |

不改输入目录:`stage2_input_post` 被 Stage 2 早先的结果引用,原地改 config 会破坏溯源;另建输入目录则要换计划里写的输入路径。

## 3. 本地验证

- 单元测试:test_stage2_s 7 项、test_gpu_unot 3 项、test_stage2_inject 全部通过。
- toy 端到端冒烟:后部 config 存 500 时,`run_s.sh --smoke` 跑通;case.json 记录 `cutoff_um = 363.0757`、`config_cutoff_um = 500.0`;等价门槛差 0.0。
- GPU 门槛路径(torch CPU 模拟):不给参数时元素数 576(按 500),给 363.0757 时 572,与 `stage2_s.py` 构建的问题一致。

## 4. 运行环境修复(不涉及代码)

服务器的 git 只在 `/opt/conda/bin`,`conda run` 时不在 PATH 上,第一次冒烟因此在单元测试停下。已把 git 软链接到 `~/.local/bin/git`。

## 5. 两次失败的冒烟目录

改名保留,不删除:`analysis_s_smoke_failed_nogit/`、`analysis_s_smoke_failed_cutoff/`。都没有求解过带 S 的 OT。
