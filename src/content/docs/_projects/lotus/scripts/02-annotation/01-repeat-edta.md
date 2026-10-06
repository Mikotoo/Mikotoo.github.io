---
title: '重复序列屏蔽（EDTA）'
sidebar:
  label: '重复序列屏蔽（EDTA）'
  order: 1
---


## 中文

**作用**：用 EDTA 从头注释并屏蔽 Gifu 基因组的转座元件，再用 LTR_retriever 计算 LAI。
**入口**：`process.sh`（屏蔽）；`LAI.sh`（LAI 打分）。
**输入**：`Gifu_v1.0.fasta`；`LAI.sh` 另需 `<genome>.mod.harvest.combine.scn` 与 `<genome>.mod.finder.combine.scn`。
**输出**：`Gifu_v1.0.fasta.mod.*`，其中 `Gifu_v1.0.fasta.mod.MAKER.masked` 供 `04_braker/Gifu_braker.sh` 使用；LAI 结果与作业日志 `LAI.out`。
**运行**：在存放基因组与 EDTA 输出的目录中执行 `EDTA.pl --genome Gifu_v1.0.fasta --overwrite 1 --sensitive 1 --anno 1 --threads 264 --force 1`，随后 `bash LAI.sh`（须先把脚本中的大豆基因组名 `ChiHei_v1.0.fasta` 及其 `-inharvest`/`-infinder` 文件名换成 Gifu 的对应文件）。
**工具**：EDTA `EDTA.pl --sensitive 1 --anno 1 --threads 264 --force 1`；LTR_retriever `-genome <genome> -inharvest <scn> -infinder <scn> -threads 88`。

<details>
<summary>English</summary>
**Purpose**: De novo TE annotation and masking of the Gifu assembly with EDTA, then LAI scoring of the EDTA output with LTR_retriever.
**Entry point**: `process.sh` (masking); `LAI.sh` (LAI scoring).
**Inputs**: `Gifu_v1.0.fasta`; `LAI.sh` additionally needs `<genome>.mod.harvest.combine.scn` and `<genome>.mod.finder.combine.scn`.
**Outputs**: `Gifu_v1.0.fasta.mod.*`, including `Gifu_v1.0.fasta.mod.MAKER.masked`, consumed by `04_braker/Gifu_braker.sh`; the LAI result with job stdout in `LAI.out`.
**Run**: from the directory holding the genome and the EDTA output, `EDTA.pl --genome Gifu_v1.0.fasta --overwrite 1 --sensitive 1 --anno 1 --threads 264 --force 1`, then `bash LAI.sh` (first replace the soybean genome name `ChiHei_v1.0.fasta` and its `-inharvest`/`-infinder` files with the Gifu ones).
**Tools**: EDTA `EDTA.pl --sensitive 1 --anno 1 --threads 264 --force 1`; LTR_retriever `-genome <genome> -inharvest <scn> -infinder <scn> -threads 88`.

</details>


