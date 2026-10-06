---
title: '补洞'
sidebar:
  label: '补洞'
  order: 6
---


## 中文
**作用**：用两条独立路线填补草图中的剩余 gap：TGS-GapCloser2 用 ONT UL reads，quarTeT GapFiller 用 contig 集作填充序列。
**入口**：`gapcloser.sh` 与 `quarTeT.sh`（互不调用，各自独立提交）
**输入**：`chr2.fa`（本地草图，两脚本共用）；ONT UL reads ``{PROJ_LOTUS_ZC}`/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa`；填充序列 ``{PROJ_LOTUS_ZC}`/Gifu_hifionly/04.tel_check/contigs.fa`。
**输出**：TGS-GapCloser2 的 `ont_fill_out.*`（填充结果与日志）；quarTeT 的 `Gifu_gapfill*`（填充后的基因组与绘图文件）。
**运行**：`bash gapcloser.sh`；`bash quarTeT.sh`
**工具**：TGS-GapCloser2 `--scaff chr2.fa --output ont_fill_out --ne --thread 64`；quarTeT（`quartet.py GapFiller`）`-d chr2.fa -g <contigs.fa> -p Gifu_gapfill -f 20000 -i 90 -t 88`。

<details>
<summary>English</summary>
**Purpose**: Close the remaining gaps in the draft by two independent routes: TGS-GapCloser2 with ONT UL reads, and quarTeT GapFiller using the contig set as filling sequence.
**Entry point**: `gapcloser.sh` and `quarTeT.sh` (neither invokes the other; each is submitted on its own)
**Inputs**: `chr2.fa` (local draft, shared by both scripts); ONT UL reads ``{PROJ_LOTUS_ZC}`/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa`; filling sequence ``{PROJ_LOTUS_ZC}`/Gifu_hifionly/04.tel_check/contigs.fa`.
**Outputs**: TGS-GapCloser2 `ont_fill_out.*` (filled sequence and logs); quarTeT `Gifu_gapfill*` (the filled genome and drawing files).
**Run**: `bash gapcloser.sh`; `bash quarTeT.sh`
**Tools**: TGS-GapCloser2 `--scaff chr2.fa --output ont_fill_out --ne --thread 64`; quarTeT (`quartet.py GapFiller`) `-d chr2.fa -g <contigs.fa> -p Gifu_gapfill -f 20000 -i 90 -t 88`.

</details>


