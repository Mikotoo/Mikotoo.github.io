---
title: '03 端粒'
sidebar:
  label: '03 端粒'
  order: 3
---


## 中文

**作用**：用自研代码表征 Gifu T2T 组装的端粒：扫描 7 bp 端粒重复单元、换算为 100 bp 窗口密度、在每个染色体末端调用端粒区间、统计其上的 ONT/HiFi 深度并提取跨端粒读段。
**入口**：`run_telomere.sh`（内部调用 `telomere.py` 的 `find-tel`/`coverage`/`alignments`/`plot` 子命令）。
**输入**：`00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta`（无 `.fai` 时先 `samtools faidx`）；`02_assembly/Gifu/08.quality/01.mapping/01.ont/ont_q30.bam` 与 `.../02.hifi/ccs_q30.bam`（无 `.bai` 时先 `samtools index`）；环境变量 `TELOMERE_*`（工作目录、前缀、步骤、Python/samtools/bedtools 路径、输入 FASTA/BAM、bin 大小、绘图侧翼与线程，默认值均为 Gifu 生产值）。
**输出**：工作目录中以前缀 `${PREFIX}`（默认 `GifuT2T`）命名的表：`*.tel_unit.bed`、`*_tel.100.result`+`*_tel.100.detail.tsv`、`*.telomere_coords.tsv`、`*.telomere_coverage.100bp.tsv`+`.summary.tsv`、`*.telomere_alignments.tsv`；`plots/` 下的端粒长度、末端、覆盖度汇总与读段分布 PDF。
**运行**：`csub run_telomere.sh`（全流程）；`TELOMERE_STEPS=plot bash run_telomere.sh`（仅重绘）；顺序为 `find-tel` → `coverage` → `alignments` → `plot`；脚本内 `SCRIPT_DIR` 写死为 `${PROJECT_ROOT}/03_annotation/Gifu/01.tel`（`csub` 会复制脚本），须按实际路径修改。
**工具**：`telomere.py`（自研）：bin 100 bp、密度 0.5、单元最大间隔 14 bp、最多 5 个凹陷箱、绘图侧翼 10,000 bp，跨端读段 `--min-overlap 50`、`--max-reads-per-arm 2500`、`--min-mapq 0`；samtools 1.18 `faidx`、`index -@16`；Python ``{CONDA_PY3}`/bin/python3`；bedtools 经 `--bedtools` 传入；`alignments` 需 pysam，`plot` 需 matplotlib。

<details>
<summary>English</summary>
**Purpose**: Characterise the telomeres of the Gifu T2T assembly with in-house code: scan 7-bp telomere repeat units, convert them to 100-bp window density, call the telomere tract at each chromosome end, measure ONT/HiFi depth over the calls and extract spanning reads.
**Entry point**: `run_telomere.sh` (which calls the `find-tel`/`coverage`/`alignments`/`plot` subcommands of `telomere.py`).
**Inputs**: `00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta` (`samtools faidx` runs when the `.fai` is missing); `02_assembly/Gifu/08.quality/01.mapping/01.ont/ont_q30.bam` and `.../02.hifi/ccs_q30.bam` (`samtools index` runs when the `.bai` is missing); the `TELOMERE_*` environment variables (workdir, prefix, steps, Python/samtools/bedtools paths, input FASTA/BAMs, bin size, plot flank and threads), each defaulted to the production Gifu value.
**Outputs**: tables in the workdir prefixed with `${PREFIX}` (default `GifuT2T`): `*.tel_unit.bed`, `*_tel.100.result`+`*_tel.100.detail.tsv`, `*.telomere_coords.tsv`, `*.telomere_coverage.100bp.tsv`+`.summary.tsv`, `*.telomere_alignments.tsv`; PDFs under `plots/` for telomere lengths, chromosome ends, coverage summary and read distribution.
**Run**: `csub run_telomere.sh` (full pipeline); `TELOMERE_STEPS=plot bash run_telomere.sh` (re-plot only); the order is `find-tel` → `coverage` → `alignments` → `plot`; `SCRIPT_DIR` is hardcoded to `${PROJECT_ROOT}/03_annotation/Gifu/01.tel` (because `csub` copies the script), so it must be edited to the real path.
**Tools**: `telomere.py` (in-house): bin 100 bp, density 0.5, max unit gap 14 bp, max dip bins 5, plot flank 10,000 bp, spanning reads `--min-overlap 50`, `--max-reads-per-arm 2500`, `--min-mapq 0`; samtools 1.18 `faidx`, `index -@16`; Python ``{CONDA_PY3}`/bin/python3`; bedtools passed through `--bedtools`; `alignments` requires pysam and `plot` requires matplotlib.

</details>


