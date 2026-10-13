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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>01.Lj_nodulation_mutants/align.sh</code> · Shell · 67 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/01.Lj_nodulation_mutants/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_01
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH

INDEX="{PROJ_04_LOTUS_GENOME}/05_RNA/03.mutant/index/GifuT2T"
REF="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>02.Lj_roothair_inoculation/align.sh</code> · Shell · 67 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/02.Lj_roothair_inoculation/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_02
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH

INDEX="{PROJ_04_LOTUS_GENOME}/05_RNA/03.mutant/index/GifuT2T"
REF="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>03.Lj_nodule_cpop1/align.sh</code> · Shell · 67 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/03.Lj_nodule_cpop1/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_03
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH

INDEX="{PROJ_04_LOTUS_GENOME}/05_RNA/03.mutant/index/GifuT2T"
REF="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>05.Gm_root_zones/align.sh</code> · Shell · 68 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/05.Gm_root_zones/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_05
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH


INDEX="{HOME_ZC}/Project/01_ChiDou_genome/08_RNAseq/01.index/CH"
REF="{HOME_ZC}/Project/01_ChiDou_genome/00_assembly_result/ChiHei.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>06.Gm_roothair/align.sh</code> · Shell · 68 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/06.Gm_roothair/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_06
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH


INDEX="{HOME_ZC}/Project/01_ChiDou_genome/08_RNAseq/01.index/CH"
REF="{HOME_ZC}/Project/01_ChiDou_genome/00_assembly_result/ChiHei.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>07.Mt_nodfactor/align.sh</code> · Shell · 67 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/07.Mt_nodfactor/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat_07
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH

INDEX="{HOME_ZC}/Project/Ref/Medicago_v4/Medic_v4"
REF="{HOME_ZC}/Project/Ref/Medicago_v4/Medic_truncatula.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
```

</details>

<details open>
<summary><code>index/index_GifuT2T.sh</code> · Shell · 19 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/02_mutant_and_public_data/index/index_GifuT2T.sh)

```bash
#!/bin/bash
#CSUB -J hisat_build
#CSUB -q c01
#CSUB -o hisat_build.out
#CSUB -e hisat_build.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# 构建 Lotus japonicus GifuT2T 的 hisat2 索引

set -euo pipefail
export PATH={CONDA_PY3}/bin:$PATH

FASTA="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta"
INDEX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${INDEX_DIR}/GifuT2T"

hisat2-build -p 32 "${FASTA}" "${OUT}"
echo "[DONE] index built at ${OUT}"
```

</details>
