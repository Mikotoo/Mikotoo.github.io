---
title: '端粒检查'
sidebar:
  label: '端粒检查'
  order: 4
---


## 中文
**作用**：扫描端粒重复单元并列出富集它们的 1 kb 窗口，另独立运行 BUSCO 评估基因区完整度。
**入口**：`find_tel.sh`（调用 `get_tel.py`）；`busco.sh` 独立运行
**输入**：`contigs.fa`（端粒扫描与 `busco.sh` 共用）。
**输出**：`Gifu.tel_unit.bed`、`contigs.fa.fai`、`Gifu.length`、`Gifu_1k.bed`、`Gifu_tel.result`（每窗口 >50 命中）；BUSCO 输出目录 `Gifu_lotus`。
**运行**：`bash find_tel.sh`；`bash busco.sh`
**工具**：Python + Biopython `SeqIO`（基序 `CCCTAAA`、`TTTAGGG` 及各位置的单碱基替换）；samtools `faidx`；bedtools `makewindows -w 1000`、`intersect -wa`；BUSCO `-l embryophyta_odb10 -m genome --cpu 64`。

<details>
<summary>English</summary>

**Purpose**: Scan the contig set for telomeric repeat units and list the 1 kb windows dense in them, and run BUSCO separately for gene-space completeness.
**Entry point**: `find_tel.sh` (calls `get_tel.py`); `busco.sh` runs independently
**Inputs**: `contigs.fa` (shared by the telomere scan and `busco.sh`).
**Outputs**: `Gifu.tel_unit.bed`, `contigs.fa.fai`, `Gifu.length`, `Gifu_1k.bed`, `Gifu_tel.result` (>50 hits per window); BUSCO output directory `Gifu_lotus`.
**Run**: `bash find_tel.sh`; `bash busco.sh`
**Tools**: Python with Biopython `SeqIO` (motifs `CCCTAAA` and `TTTAGGG` plus all single-base substitutions); samtools `faidx`; bedtools `makewindows -w 1000`, `intersect -wa`; BUSCO `-l embryophyta_odb10 -m genome --cpu 64`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>busco.sh</code> · Shell · 7 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/04_telomere_check/busco.sh)

```bash
#4.1 busco初步检查组装完整性

busco -i contigs.fa \
      -l {SOFTWARE}/busco_downloads/lineages/embryophyta_odb10 \
      -o Gifu_lotus \
      -m genome \
      --cpu 64
```

</details>

<details open>
<summary><code>find_tel.sh</code> · Shell · 8 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/04_telomere_check/find_tel.sh)

```bash
#4.2 初步检查每条contig是否组装到端粒，辅助后续挂载排序

python get_tel.py
sort -u Gifu.tel_unit.bed |sort -k 1,1 -k 1n > tmp && mv tmp Gifu.tel_unit.bed
samtools faidx contigs.fa
cat contigs.fa.fai |awk '{print $1"\t"$2}' > Gifu.length
bedtools makewindows -g Gifu.length -w 1000 > Gifu_1k.bed
bedtools intersect -a Gifu_1k.bed -b Gifu.tel_unit.bed -wa |uniq -c|awk '$1>50{print $2"\t"$3"\t"$4"\t"$1}' > Gifu_tel.result
```

</details>

<details open>
<summary><code>get_tel.py</code> · Python · 18 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/04_telomere_check/get_tel.py)

```python
from Bio import SeqIO

fasta_file = "contigs.fa"
bed_file = "Gifu.tel_unit.bed"
strings = ['CCCTAAA', 'TTTAGGG']
print("Total strings:", len(strings))
x = ['A', 'T', 'G', 'C']

with open(bed_file, 'w') as f:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        for string in strings:
            for i in range(len(seq) - 6):
                for j in range(7):
                    for m in range(len(x)):
                        modified_string = string[:j] + x[m] + string[j+1:]
                        if seq[i:i + 7] == modified_string:
                            f.write(record.id + "\t" + str(i) + "\t" + str(i + 6) + "\t" + modified_string + "\n")
```

</details>
