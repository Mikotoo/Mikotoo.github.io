---
title: '示例：着丝粒串联重复'
sidebar:
  label: '示例：着丝粒串联重复'
  order: 5
---


## 中文
**作用**：在一张极小的手工簇图上运行着丝粒流程的并查集合并步骤 `build_supercluster_unionfind.py`。
**入口**：`run_demo.sh`
**输入**：`input/clusters.list`（8 个重复簇，每行一个）与 `input/pairs.jaccard_ge0.7.tsv`（6 条不低于 0.7 截断的 Jaccard 对，制表符分隔 `A B Jaccard`：c1-c2 0.91、c1-c3 0.88、c2-c3 0.95、c4-c5 0.77、c6-c7 0.83、c7-c8 0.81）。
**输出**：`work/cluster2super.tsv` 与 `work/super_summary.tsv`；参考结果在 `expected_output/`，内容为 `super_0001 3 c1,c2,c3`、`super_0002 2 c4,c5`、`super_0003 3 c6,c7,c8`。
**运行**：`bash run_demo.sh`（把输入复制进 `work/` 后在其中运行该步骤，并打印两个输出文件）。
**工具**：Python 3 标准库，无第三方依赖；运行时间远低于 1 秒。

<details>
<summary>English</summary>

**Purpose**: Run the union-find merging step of the centromere pipeline (`build_supercluster_unionfind.py`) on a tiny hand-made cluster graph.
**Entry point**: `run_demo.sh`
**Inputs**: `input/clusters.list` (8 repeat clusters, one per line) and `input/pairs.jaccard_ge0.7.tsv` (6 pairwise Jaccard indices at or above the 0.7 cutoff, tab-separated `A B Jaccard`: c1-c2 0.91, c1-c3 0.88, c2-c3 0.95, c4-c5 0.77, c6-c7 0.83, c7-c8 0.81).
**Outputs**: `work/cluster2super.tsv` and `work/super_summary.tsv`; reference copies are in `expected_output/` and read `super_0001 3 c1,c2,c3`, `super_0002 2 c4,c5`, `super_0003 3 c6,c7,c8`.
**Run**: `bash run_demo.sh` (copies the inputs into `work/`, runs the step there and prints both output files).
**Tools**: Python 3 standard library, no third-party dependency; runtime is well under a second.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>expected_output/cluster2super.tsv</code> · text · 8 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/05_centromere/expected_output/cluster2super.tsv)

```text
c1	super_0001
c2	super_0001
c3	super_0001
c4	super_0002
c5	super_0002
c6	super_0003
c7	super_0003
c8	super_0003
```

</details>

<details open>
<summary><code>expected_output/super_summary.tsv</code> · text · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/05_centromere/expected_output/super_summary.tsv)

```text
super_0001	3	c1,c2,c3
super_0002	2	c4,c5
super_0003	3	c6,c7,c8
```

</details>

<details open>
<summary><code>input/clusters.list</code> · text · 8 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/05_centromere/input/clusters.list)

```text
c1
c2
c3
c4
c5
c6
c7
c8
```

</details>

<details open>
<summary><code>input/pairs.jaccard_ge0.7.tsv</code> · text · 6 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/05_centromere/input/pairs.jaccard_ge0.7.tsv)

```text
c1	c2	0.91
c1	c3	0.88
c2	c3	0.95
c4	c5	0.77
c6	c7	0.83
c7	c8	0.81
```

</details>

<details open>
<summary><code>run_demo.sh</code> · Shell · 28 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/demo/05_centromere/run_demo.sh)

```bash
#!/bin/bash
# Demonstration of build_supercluster_unionfind.py on a three-supercluster example.
#
# input/    clusters.list and pairs.jaccard_ge0.7.tsv (hand-made, tiny)
# work/     scratch directory the script writes into
# expected_output/  the two files the script produces, kept for comparison
#
# Usage:  bash run_demo.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SCRIPT="$REPO/scripts/05_centromere/01_supercluster/blast/build_supercluster_unionfind.py"

rm -rf "$HERE/work"
mkdir -p "$HERE/work"
cp "$HERE/input/clusters.list" "$HERE/input/pairs.jaccard_ge0.7.tsv" "$HERE/work/"

python3 "$SCRIPT" --outdir "$HERE/work" --jac 0.7

echo
echo "--- cluster2super.tsv ---"
cat "$HERE/work/cluster2super.tsv"
echo
echo "--- super_summary.tsv ---"
cat "$HERE/work/super_summary.tsv"
echo
echo "Compare with $HERE/expected_output/"
```

</details>
