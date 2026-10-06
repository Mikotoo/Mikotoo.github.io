---
title: '突变体与公共数据'
sidebar:
  label: '突变体与公共数据'
  order: 2
---


## 中文
**作用**：定量 *Lotus japonicus* 结瘤突变体以及公共大豆、*Medicago truncatula* 数据集中的 LjSYN1/LjSYN2 同源基因。
**入口**：`index/index_GifuT2T.sh`；各数据集目录下的 `align.sh`
**输入**：每个数据集目录的 `samples.conf`（`sample`、`fq1`、可选 `fq2`，`#` 开头行跳过）及其指向的 FASTQ；`align.sh` 头部变量指定的索引与注释 GTF。
**输出**：每样本 `<sample>.bam(.bai)`、`align/<sample>/<sample>.gtf|.tsv`、`<sample>_hisat.log`；`prepDE.list`、`count/prepDE_counts.csv`、`count/prepDE_transcript.csv`；`index/index_GifuT2T.sh` 就地写 `GifuT2T.*.ht2`。
**运行**：`bash index/index_GifuT2T.sh`；在每个数据集目录内 `bash align.sh`（输出写当前目录，且须先把脚本内 ``{HOME_ZC}``、``{CONDA_PY3}`` 的硬编码绝对路径改成自己的）。
**工具**：HISAT2（`hisat2-build -p 32`；比对 `-p 32`）、samtools（`view`、`sort`、`index`）、StringTie（`-p -G -e -B -A`）、`prepDE.py`（Python 3）；集群投递头 `#CSUB -q c01 -n 32`。

<details>
<summary>English</summary>
**Purpose**: Quantify the LjSYN1/LjSYN2 homologues in *Lotus japonicus* nodulation mutants and in public soybean and *Medicago truncatula* datasets.
**Entry point**: `index/index_GifuT2T.sh`; the `align.sh` in each dataset directory
**Inputs**: each dataset directory's `samples.conf` (`sample`, `fq1`, optional `fq2`; `#`-prefixed lines are skipped) and the FASTQ files it points to; the index and annotation GTF set in the header variables of `align.sh`.
**Outputs**: per sample `<sample>.bam(.bai)`, `align/<sample>/<sample>.gtf|.tsv`, `<sample>_hisat.log`; `prepDE.list`, `count/prepDE_counts.csv`, `count/prepDE_transcript.csv`; `index/index_GifuT2T.sh` writes the `GifuT2T.*.ht2` index next to itself.
**Run**: `bash index/index_GifuT2T.sh`; `bash align.sh` inside each dataset directory (output goes to the current directory, and the hardcoded ``{HOME_ZC}`` / ``{CONDA_PY3}`` absolute paths in the script must be replaced first).
**Tools**: HISAT2 (`hisat2-build -p 32`; alignment `-p 32`), samtools (`view`, `sort`, `index`), StringTie (`-p -G -e -B -A`), `prepDE.py` (Python 3); cluster job header `#CSUB -q c01 -n 32`.

</details>


