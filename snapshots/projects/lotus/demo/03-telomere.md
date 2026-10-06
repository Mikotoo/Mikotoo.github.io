---
title: '示例：端粒识别'
sidebar:
  label: '示例：端粒识别'
  order: 3
---


## 中文
**作用**：在一条短的合成染色体上运行自研端粒判定代码 `telomere.py find-tel`，无需真实组装即可验证算法。
**入口**：`run_demo.sh`
**输入**：`input/demo_chr.fa`（2,034 bp：0–280 为 40×`TTTAGGG`，280–1720 为填充序列，1720–2000 为 40×`CCCTAAA`；重复单元 7 bp，每条臂 40 个单元共 280 bp）与 `input/demo_chr.fa.fai`（脚本要求的 FASTA 索引）。
**输出**：`work/` 下的坐标表 `demo.telomere_coords.tsv`、去重单元表 `demo_tel_unit.bed` 与两个分箱密度文件 `demo_tel.100.result`、`demo_tel.100.detail.tsv`；参考结果在 `expected_output/`。
**运行**：`bash run_demo.sh`（重建 `work/`、运行该步骤并打印坐标表）。
**工具**：Python 3 标准库，无第三方依赖；运行时间远低于 1 秒。

<details>
<summary>English</summary>
**Purpose**: Run the custom telomere-calling code `telomere.py find-tel` on a short synthetic chromosome so the algorithm can be verified without a genome assembly.
**Entry point**: `run_demo.sh`
**Inputs**: `input/demo_chr.fa` (2,034 bp: 0–280 is 40×`TTTAGGG`, 280–1720 filler, 1720–2000 is 40×`CCCTAAA`; the repeat unit is 7 bp and each arm holds 40 units spanning 280 bp) and `input/demo_chr.fa.fai`, the FASTA index the script requires.
**Outputs**: in `work/`, the coordinate table `demo.telomere_coords.tsv`, the deduplicated unit list `demo_tel_unit.bed` and the two bin-density files `demo_tel.100.result` and `demo_tel.100.detail.tsv`; reference copies are in `expected_output/`.
**Run**: `bash run_demo.sh` (rebuilds `work/`, runs the step and prints the coordinate table).
**Tools**: Python 3 standard library, no third-party dependency; runtime is well under a second.

</details>


