---
title: '一致性序列'
sidebar:
  label: '一致性序列'
  order: 3
---


## 中文

**作用**：为四个着丝粒重复单元 `unit_161`、`unit_172`、`unit_186`、`unit_330` 构建共识序列：从每个基因组提取该单元的命中、用 MAFFT 比对、折叠成共识，再把 Gifu 与 MG20 的共识相互比对得到跨材料共识。
**入口**：`run_unit_mafft_consensus.sh`
**输入**：`Lotus_GifuT2T_v1.0.fasta`、`Lotus_MG20T2T_v1.0.fasta`；每个单元的 BED `gifu_unit_<unit>.bed` 与 `mg20_unit_<unit>.bed`（BED6）。
**输出**：`separate_then_merge_consensus/` 下 `01_seqs/<tag>_<unit>.fa`、`02_msa/<tag>_<unit>.aln.fa`、`03_consensus/<tag>_<unit>.cons.fa`、`03_dist/<tag>_<unit>.dist.tsv`、`04_merge2cons_msa/<unit>.Gifu_MG20.cons.{fa,aln.fa}`、`05_final/<unit>.final.cons.fa` 与 `<unit>.final.dist.tsv`。
**运行**：`bash run_unit_mafft_consensus.sh`（须在含两个基因组 FASTA 与 `gifu_`/`mg20_` 前缀 BED 的目录中运行；`blast.sh` 除末行外均为注释）；单次调用：`python get_consensus.py -i aln.fa -o cons.fa -d dist.tsv -t 0.5 --include-gap --name <name>`。
**工具**：MAFFT `--thread 32 --auto`；bedtools `getfasta -name -s`、`jaccard`、`intersect`、`sort`、`merge`；BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn-short -word_size 7 -evalue 1e-10`（`outfmt 6`）；`get_consensus.py` 用 Biopython `AlignIO`，`Levenshtein` 可选（内置纯 Python 回退）；`THRESH=0.5`。

<details>
<summary>English</summary>
**Purpose**: Build consensus sequences for the four centromeric repeat units `unit_161`, `unit_172`, `unit_186` and `unit_330`: hits of each unit are extracted from each genome, aligned with MAFFT and collapsed into a consensus, and the Gifu and MG20 consensuses are then aligned to each other to give a cross-accession consensus.
**Entry point**: `run_unit_mafft_consensus.sh`
**Inputs**: `Lotus_GifuT2T_v1.0.fasta`, `Lotus_MG20T2T_v1.0.fasta`; the per-unit BEDs `gifu_unit_<unit>.bed` and `mg20_unit_<unit>.bed` (BED6).
**Outputs**: under `separate_then_merge_consensus/`: `01_seqs/<tag>_<unit>.fa`, `02_msa/<tag>_<unit>.aln.fa`, `03_consensus/<tag>_<unit>.cons.fa`, `03_dist/<tag>_<unit>.dist.tsv`, `04_merge2cons_msa/<unit>.Gifu_MG20.cons.{fa,aln.fa}`, `05_final/<unit>.final.cons.fa` and `<unit>.final.dist.tsv`.
**Run**: `bash run_unit_mafft_consensus.sh` (must run from a directory holding the two genome FASTAs and the `gifu_`/`mg20_`-prefixed BEDs; `blast.sh` is commented out except its last line); a single call: `python get_consensus.py -i aln.fa -o cons.fa -d dist.tsv -t 0.5 --include-gap --name <name>`.
**Tools**: MAFFT `--thread 32 --auto`; bedtools `getfasta -name -s`, `jaccard`, `intersect`, `sort`, `merge`; BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn-short -word_size 7 -evalue 1e-10` (`outfmt 6`); `get_consensus.py` uses Biopython `AlignIO` with `Levenshtein` optional (a pure-Python fallback is built in); `THRESH=0.5`.

</details>


