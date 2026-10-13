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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>expected_output/demo.telomere_coords.tsv</code> · text · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/03_telomere/expected_output/demo.telomere_coords.tsv)

```text
chromosome	end	bin_region_start	bin_region_end	telomere_start	telomere_end	length_bp	dominant_motif	n_units	mean_bin_density	max_bin_density
demo_chr1	5p	0	300	0	280	280	TTTAGGG	40	0.9333	1.0000
demo_chr1	3p	1700	2000	1720	2000	280	CCCTAAA	40	0.9333	1.0000
```

</details>

<details open>
<summary><code>expected_output/demo_tel.100.detail.tsv</code> · text · 22 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/03_telomere/expected_output/demo_tel.100.detail.tsv)

```text
chrom	bin_start	bin_end	count	density
demo_chr1	0	100	15	1
demo_chr1	100	200	15	1
demo_chr1	200	300	12	0.8
demo_chr1	300	400	0	0
demo_chr1	400	500	0	0
demo_chr1	500	600	0	0
demo_chr1	600	700	0	0
demo_chr1	700	800	0	0
demo_chr1	800	900	0	0
demo_chr1	900	1000	0	0
demo_chr1	1000	1100	0	0
demo_chr1	1100	1200	0	0
demo_chr1	1200	1300	0	0
demo_chr1	1300	1400	0	0
demo_chr1	1400	1500	0	0
demo_chr1	1500	1600	0	0
demo_chr1	1600	1700	0	0
demo_chr1	1700	1800	12	0.8
demo_chr1	1800	1900	15	1
demo_chr1	1900	2000	15	1
demo_chr1	2000	2034	0	0
```

</details>

<details open>
<summary><code>run_demo.sh</code> · Shell · 30 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/03_telomere/run_demo.sh)

```bash
#!/bin/bash
# Demonstration of the telomere pipeline on a short synthetic chromosome.
#
# input/demo_chr.fa        2 kb sequence: 40 x TTTAGGG, filler, 40 x CCCTAAA
# input/demo_chr.fa.fai    FASTA index required by telomere.py
# expected_output/         the four files the run produces
#
# Only the find-tel subcommand is exercised; coverage, alignments and plot need BAM
# files and a full assembly, which are not part of a small demo.
#
# Usage:  bash run_demo.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SCRIPT="$REPO/scripts/03_telomere/telomere.py"

rm -rf "$HERE/work"
mkdir -p "$HERE/work"

python3 "$SCRIPT" find-tel \
  --fasta-asm "$HERE/input/demo_chr.fa" \
  --workdir "$HERE/work" \
  --prefix demo

echo
echo "--- demo.telomere_coords.tsv ---"
cat "$HERE/work/demo.telomere_coords.tsv"
echo
echo "Compare with $HERE/expected_output/"
```

</details>
