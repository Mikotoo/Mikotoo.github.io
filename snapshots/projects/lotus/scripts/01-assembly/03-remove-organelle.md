---
title: '去除细胞器序列'
sidebar:
  label: '去除细胞器序列'
  order: 3
---


## 中文
**作用**：将 contig 比对 NCBI `nt`，经 GI → taxid → 学名汇总每个物种的命中数，据此识别并去除质体/线粒体来源的 contig。
**入口**：`03.1.blast.sh`（内部调用 `03.2.get_gi.py`、`03.3.get_name.py`），随后运行 `03.4.get_ratio.py`
**输入**：`contigs.fa`；`nt` 库 ``{SOFTWARE_ZC}`/nt_lib/nt`；`cutted_nucl_gb.accession2taxid` 与 `names.dmp`（须放在工作目录）。
**输出**：`nt.blast.xml`、`tiqu_gi.txt`、`scientific_name.txt`、`final_result.txt`（`Name/Hit_reads/percentage1/percentage2`，按命中数排序）。
**运行**：`bash 03.1.blast.sh`，然后 `python 03.4.get_ratio.py`
**工具**：`blastn -outfmt 5 -evalue 1e-5 -num_threads 88`；Python 仅用标准库（`re`、`collections`），以 `python` 调用。

<details>
<summary>English</summary>
**Purpose**: BLAST the contigs against NCBI `nt`, resolve GI → taxid → scientific name and count hits per species to identify and remove plastid/mitochondrial (organelle) contigs.
**Entry point**: `03.1.blast.sh` (calls `03.2.get_gi.py` and `03.3.get_name.py`), after which `03.4.get_ratio.py` runs separately
**Inputs**: `contigs.fa`; `nt` database ``{SOFTWARE_ZC}`/nt_lib/nt`; `cutted_nucl_gb.accession2taxid` and `names.dmp` (must be in the working directory).
**Outputs**: `nt.blast.xml`, `tiqu_gi.txt`, `scientific_name.txt`, `final_result.txt` (`Name/Hit_reads/percentage1/percentage2`, sorted by hit count).
**Run**: `bash 03.1.blast.sh`, then `python 03.4.get_ratio.py`
**Tools**: `blastn -outfmt 5 -evalue 1e-5 -num_threads 88`; Python with the standard library only (`re`, `collections`), invoked as `python`.

</details>


