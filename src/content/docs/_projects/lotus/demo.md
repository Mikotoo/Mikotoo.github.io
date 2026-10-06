---
title: '可运行示例'
sidebar:
  label: '可运行示例'
  order: -1
---


## 中文
**作用**：本仓库自研代码的最小可运行演示：每个 demo 自包含、只用手工或合成的小输入，并随附一次正确运行应产出的结果。
**入口**：`03_telomere/run_demo.sh`、`05_centromere/run_demo.sh`；每个 demo 的目录含 `input/`、`expected_output/` 与 `run_demo.sh`。
**输入**：`03_telomere/input/` 的 2 kb 合成染色体 FASTA 及其 `.fai` 索引；`05_centromere/input/` 的 `clusters.list`（8 个重复簇）与 `pairs.jaccard_ge0.7.tsv`（6 条 Jaccard 对）。
**输出**：各 demo 的 `work/`，即本次运行的结果，用于与 `expected_output/` 比较。
**运行**：`bash demo/03_telomere/run_demo.sh`；`bash demo/05_centromere/run_demo.sh`（脚本按自身位置定位仓库，可从任意工作目录调用）。
**工具**：只需 Python 3 解释器与标准库，不用第三方包、参考基因组或测序数据；两个 demo 都在 1 秒内跑完。

<details>
<summary>English</summary>
**Purpose**: Small, runnable demonstrations of this repository's custom code; each demo is self-contained, uses a tiny hand-made or synthetic input, and ships the output a correct run produces.
**Entry point**: `03_telomere/run_demo.sh`, `05_centromere/run_demo.sh`; each demo directory holds `input/`, `expected_output/` and `run_demo.sh`.
**Inputs**: the 2 kb synthetic chromosome FASTA and its `.fai` index in `03_telomere/input/`; `clusters.list` (8 repeat clusters) and `pairs.jaccard_ge0.7.tsv` (6 Jaccard pairs) in `05_centromere/input/`.
**Outputs**: each demo's `work/` directory, holding the result of the current run for comparison with `expected_output/`.
**Run**: `bash demo/03_telomere/run_demo.sh`; `bash demo/05_centromere/run_demo.sh` (the scripts resolve the repository layout relative to their own location, so they can be invoked from any working directory).
**Tools**: only a Python 3 interpreter and the standard library — no third-party packages, no reference genome and no sequencing data; both demos run in under a second.

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [示例：端粒识别](/_projects/lotus/demo/03-telomere/) | — |
| [示例：着丝粒串联重复](/_projects/lotus/demo/05-centromere/) | — |


