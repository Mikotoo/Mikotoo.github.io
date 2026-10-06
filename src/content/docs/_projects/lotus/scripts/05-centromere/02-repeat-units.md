---
title: '重复单元'
sidebar:
  label: '重复单元'
  order: 2
---


## 中文

**作用**：把 TRF 记录变成每条染色体的候选着丝粒单体序列并跨染色体连接：对每个重复单元做 canonical 化（旋转与反向互补）、选出最能解释该家族 k×L 长度结构的单体长度 L，再自比对单体得到染色体对连接数与连通分量。
**入口**：`run.sh`
**输入**：`gifu_cent.bed`、`mg20_cent.bed`（TRF 记录，单元序列在指定列）；`gifu_chr_top_monomer.tsv`、`mg20_chr_top_monomer.tsv`（chr、单元长度、单元序列、span）；`cluster.py` 用 `self.sim90.tsv`。
**输出**：`run.sh` 写 `lotus_chr_top_monomer.fa`、`lotus_monomer.*` BLAST 库、`self.allvsall.tsv`、`self.noSelf.tsv`、`self.sim80.tsv`、`self.sim90.tsv`、`top_hit_per_query.tsv`、`chr_pair_links.tsv`；`canon_unit_and_length.py` 写扩展 TSV 与 `length_summary.tsv`；`pick_monomer_from_clusters.py` 写 `cluster_monomer.tsv`；`cluster.py` 把簇打印到 stdout。
**运行**：`bash run.sh`；`python3 canon_unit_and_length.py -i mg20_cent.bed -o mg20.units.canonical.tsv --summary mg20.length_summary.tsv --seq-col 9`；`python3 pick_monomer_from_clusters.py -i clusters.tsv -o cluster_monomer.tsv --tol 6 --max-k 6 --top 5`；`python3 cluster.py`（`run.sh` 第 62 行以上为注释记录需单独运行；`cluster.py` 硬编码 `self.sim90.tsv`，须在本目录运行）。
**工具**：BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn-short -dust no -soft_masking false -word_size 7 -evalue 1e-10 -max_target_seqs 10000`，`outfmt 6`；awk/sort 过滤 `pident>=80 && len>=100`（`self.sim80.tsv`）与 `pident>=90 && len>=140`（`self.sim90.tsv`）；三个 Python 脚本仅用标准库；`#CSUB` 队列 c01、88 槽。

<details>
<summary>English</summary>
**Purpose**: Turn TRF records into a per-chromosome set of candidate centromeric monomer sequences and link them across chromosomes: canonicalise each repeat unit (rotation and reverse complement), pick the monomer length L that best explains the family's k×L length structure, and self-align the monomers to build chromosome-pair link counts and connected components.
**Entry point**: `run.sh`
**Inputs**: `gifu_cent.bed`, `mg20_cent.bed` (TRF records, unit sequence in the configured column); `gifu_chr_top_monomer.tsv`, `mg20_chr_top_monomer.tsv` (chr, unit length, unit sequence, span); `self.sim90.tsv` for `cluster.py`.
**Outputs**: `run.sh` writes `lotus_chr_top_monomer.fa`, the `lotus_monomer.*` BLAST database, `self.allvsall.tsv`, `self.noSelf.tsv`, `self.sim80.tsv`, `self.sim90.tsv`, `top_hit_per_query.tsv` and `chr_pair_links.tsv`; `canon_unit_and_length.py` writes the extended TSV and `length_summary.tsv`; `pick_monomer_from_clusters.py` writes `cluster_monomer.tsv`; `cluster.py` prints clusters to stdout.
**Run**: `bash run.sh`; `python3 canon_unit_and_length.py -i mg20_cent.bed -o mg20.units.canonical.tsv --summary mg20.length_summary.tsv --seq-col 9`; `python3 pick_monomer_from_clusters.py -i clusters.tsv -o cluster_monomer.tsv --tol 6 --max-k 6 --top 5`; `python3 cluster.py` (everything above line 62 of `run.sh` is a commented record to be run separately; `cluster.py` hardcodes `self.sim90.tsv` and must run in this directory).
**Tools**: BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn-short -dust no -soft_masking false -word_size 7 -evalue 1e-10 -max_target_seqs 10000` with `outfmt 6`; awk/sort filters `pident>=80 && len>=100` (`self.sim80.tsv`) and `pident>=90 && len>=140` (`self.sim90.tsv`); all three Python scripts are standard library only; `#CSUB` queue c01, 88 slots.

</details>


