---
title: 'rDNA 拷贝数'
sidebar:
  label: 'rDNA 拷贝数'
  order: 3
---


## 中文

**作用**：从测序深度估计 Gifu 与 MG20 的 45S 与 5S rDNA 拷贝数；把 T2T 装配的 rDNA 阵列屏蔽后追加一个 45S 与一个 5S 代表重复构成折叠参考，WGS/HiFi reads 比对后用 mosdepth 量深度，组件深度除以单拷贝 BUSCO 位点深度中位数即得拷贝数，并与装配来源计数并列。
**入口**：`run_rDNA_CN.sh`
**输入**：装配 `{PROJ_04_LOTUS_GENOME}/00_assembly_result/{GifuT2T/Lotus_GifuT2T_v1.0.fasta,MG20T2T/Lotus_MG20T2T_v1.0.fasta}`；WGS `01_data/<acc>/WGS/<acc>_Clean_R1.fq` 与 `_R2.fq`；HiFi `01_data/<acc>/hifi/<acc>_hifi.fastq.gz`；BUSCO `02_assembly/<acc>/08.quality/02.busco/<acc>_lotus/run_embryophyta_odb10/full_table.tsv`；rDNA 阵列 BED `03_annotation/<acc>/rDNA/rRNA.bed`。
**输出**：`${OUTDIR}` = `{PROJ_04_LOTUS_GENOME}/03_annotation/rDNA_CN`，每材料一个子目录：`02_rdna/`（`rDNA_arrays.bed`、代表重复 `.fa`/`.gff`）、`03_busco/`（`*_BUSCO_singlecopy.bed`）、`04_ref/`（`*_CN_reference.fa` 及索引）、`05_map/`（`*_{WGS,HiFi}_CN.bam`）、`06_depth/`（mosdepth `*.regions.bed.gz`）、`07_cn/`（`*_{WGS,HiFi}_CN_summary.tsv`），另有合并表 `final_rDNA_CN.tsv`。
**运行**：`bash run_rDNA_CN.sh`（`#CSUB -J rDNA_CN`、队列 c01、`-n 64`、`span[hosts=1]`；一次提交内先 Gifu 后 MG20）。
**工具**：bwa `mem -t 50`；minimap2 `-ax map-hifi -t 50`；samtools `faidx`、`sort -@ 16`、`index`、`flagstat`；bedtools `maskfasta`、`intersect -v`；mosdepth `-t 8 --by <bed>`；Barrnap `--kingdom euk` 校验代表重复；Python 3（pandas、numpy）。

<details>
<summary>English</summary>

**Purpose**: Estimate the 45S and 5S rDNA copy number of Gifu and MG20 from read depth; the rDNA arrays of the T2T assembly are masked and one representative 45S and one representative 5S repeat are appended to form a collapsed reference, WGS and HiFi reads are mapped and depth is measured with mosdepth, and each rDNA component's depth divided by the median depth over single-copy BUSCO loci gives the copy number, set beside the assembly-derived counts.
**Entry point**: `run_rDNA_CN.sh`
**Inputs**: the assemblies `{PROJ_04_LOTUS_GENOME}/00_assembly_result/{GifuT2T/Lotus_GifuT2T_v1.0.fasta,MG20T2T/Lotus_MG20T2T_v1.0.fasta}`; WGS `01_data/<acc>/WGS/<acc>_Clean_R1.fq` and `_R2.fq`; HiFi `01_data/<acc>/hifi/<acc>_hifi.fastq.gz`; BUSCO `02_assembly/<acc>/08.quality/02.busco/<acc>_lotus/run_embryophyta_odb10/full_table.tsv`; the rDNA array BED `03_annotation/<acc>/rDNA/rRNA.bed`.
**Outputs**: `${OUTDIR}` = `{PROJ_04_LOTUS_GENOME}/03_annotation/rDNA_CN`, one subdirectory per accession: `02_rdna/` (`rDNA_arrays.bed`, the representative repeats `.fa`/`.gff`), `03_busco/` (`*_BUSCO_singlecopy.bed`), `04_ref/` (`*_CN_reference.fa` and indexes), `05_map/` (`*_{WGS,HiFi}_CN.bam`), `06_depth/` (mosdepth `*.regions.bed.gz`), `07_cn/` (`*_{WGS,HiFi}_CN_summary.tsv`), plus the merged table `final_rDNA_CN.tsv`.
**Run**: `bash run_rDNA_CN.sh` (`#CSUB -J rDNA_CN`, queue c01, `-n 64`, `span[hosts=1]`; Gifu then MG20 in one submission).
**Tools**: bwa `mem -t 50`; minimap2 `-ax map-hifi -t 50`; samtools `faidx`, `sort -@ 16`, `index`, `flagstat`; bedtools `maskfasta`, `intersect -v`; mosdepth `-t 8 --by <bed>`; Barrnap `--kingdom euk` to verify the representative repeats; Python 3 (pandas, numpy).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details>
<summary><code>run_rDNA_CN.sh</code> · Shell · 354 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/04_rDNA/03_copy_number/run_rDNA_CN.sh)

```bash
#!/bin/bash
#CSUB -J rDNA_CN
#CSUB -q c01
#CSUB -o rDNA_CN.out
#CSUB -e rDNA_CN.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# ==========================================================================
# rDNA copy-number estimation pipeline
# Coverage-based CNV (Illumina WGS) + HiFi independent validation
# For Gifu and MG20 T2T assemblies
# ==========================================================================

set -eo pipefail

THREADS=50                                             # max threads for mapping
MOSDEPTH_THREADS=8
PROJ={PROJ_04_LOTUS_GENOME}
OUTDIR=${PROJ}/03_annotation/rDNA_CN
mkdir -p ${OUTDIR}

# ── Paths ──────────────────────────────────────────────────────────────
GIFU_FA=${PROJ}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta
MG20_FA=${PROJ}/00_assembly_result/MG20T2T/Lotus_MG20T2T_v1.0.fasta

GIFU_WGS_R1=${PROJ}/01_data/Gifu/WGS/Gifu_Clean_R1.fq
GIFU_WGS_R2=${PROJ}/01_data/Gifu/WGS/Gifu_Clean_R2.fq
MG20_WGS_R1=${PROJ}/01_data/MG20/WGS/MG20_Clean_R1.fq
MG20_WGS_R2=${PROJ}/01_data/MG20/WGS/MG20_Clean_R2.fq

GIFU_HIFI=${PROJ}/01_data/Gifu/hifi/Gifu_hifi.fastq.gz
MG20_HIFI=${PROJ}/01_data/MG20/hifi/MG20_hifi.fastq.gz

GIFU_BUSCO=${PROJ}/02_assembly/Gifu/08.quality/02.busco/Gifu_lotus/run_embryophyta_odb10/full_table.tsv
MG20_BUSCO=${PROJ}/02_assembly/MG20/08.quality/02.busco/MG20_lotus/run_embryophyta_odb10/full_table.tsv

GIFU_RDNA_BED=${PROJ}/03_annotation/Gifu/rDNA/rRNA.bed
MG20_RDNA_BED=${PROJ}/03_annotation/MG20/rDNA/rRNA.bed

# ── Representative repeat coordinates (1-based, extracted from T2T) ────
# Gifu Chr2 45S: first complete unit (~10 kb, from first 5.8S to before second 5.8S)
GIFU_45S_CHR=Chr2
GIFU_45S_START=537436
GIFU_45S_END=547437
# Gifu Chr2 5S: first 5S cluster (~300 bp)
GIFU_5S_CHR=Chr2
GIFU_5S_START=24476202
GIFU_5S_END=24476511

# MG20 Chr2 45S: first complete unit (~13 kb)
MG20_45S_CHR=Chr2
MG20_45S_START=1528120
MG20_45S_END=1541119
# MG20 Chr2 5S: first 5S cluster (~300 bp)
MG20_5S_CHR=Chr2
MG20_5S_START=26673363
MG20_5S_END=26673672

# ── rDNA array regions for masking (BED) ───────────────────────────────
make_rdna_bed() {
    local out=$1 species=$2
    case $species in
        Gifu)
            cat > ${out} <<'BEDEOF'
Chr2	537436	949316	45S_NOR_I
Chr2	24476202	25027540	5S_array
Chr5	14943604	15344676	45S_Chr5
Chr6	8844950	8943185	45S_Chr6_I
Chr6	26125598	27223832	45S_Chr6_II
BEDEOF
            ;;
        MG20)
            cat > ${out} <<'BEDEOF'
Chr2	1528119	2857257	45S_NOR_I
Chr2	26673363	27897859	5S_array
Chr5	15013664	15659889	45S_Chr5
Chr6	15057396	15451961	45S_Chr6_I
Chr6	26255619	26622326	45S_Chr6_II
BEDEOF
            ;;
    esac
    echo "  rDNA arrays BED -> ${out}"
}


# ==========================================================================
#  process_one_species  <name>  <fa>  <wgs_R1>  <wgs_R2>  <hifi>  <busco_tsv>
# ==========================================================================
process_one_species() {
    SP_NAME=$1
    SP_FA=$2
    SP_WGS_R1=$3
    SP_WGS_R2=$4
    SP_HIFI=$5
    SP_BUSCO=$6
    SP_C45=$7
    SP_C45s=$8
    SP_C45e=$9
    SP_C5=${10}
    SP_C5s=${11}
    SP_C5e=${12}

    echo ""
    echo "================================================================"
    echo "  Processing ${SP_NAME}"
    echo "================================================================"

    WD=${OUTDIR}/${SP_NAME}
    mkdir -p ${WD}/02_rdna ${WD}/03_busco ${WD}/04_ref ${WD}/05_map ${WD}/06_depth ${WD}/07_cn

    # ── Step 1: rDNA arrays BED ─────────────────────────────────────
    echo "[${SP_NAME}] Step 1: rDNA arrays BED"
    make_rdna_bed ${WD}/02_rdna/rDNA_arrays.bed ${SP_NAME}

    # ── Step 2: Extract representative 45S & 5S repeats ──────────────
    echo "[${SP_NAME}] Step 2: Extract representative repeats"
    echo "  Extracting 45S: ${SP_C45}:${SP_C45s}-${SP_C45e}"
    echo "  Extracting 5S:  ${SP_C5}:${SP_C5s}-${SP_C5e}"
    samtools faidx ${SP_FA} "${SP_C45}:${SP_C45s}-${SP_C45e}" \
        | sed "s/^>.*/>${SP_NAME}_45S_repeat/" \
        > ${WD}/02_rdna/${SP_NAME}_45S_repeat.fa
    samtools faidx ${SP_FA} "${SP_C5}:${SP_C5s}-${SP_C5e}" \
        | sed "s/^>.*/>${SP_NAME}_5S_repeat/" \
        > ${WD}/02_rdna/${SP_NAME}_5S_repeat.fa

    # Verify with barrnap
    barrnap --kingdom euk ${WD}/02_rdna/${SP_NAME}_45S_repeat.fa \
        > ${WD}/02_rdna/${SP_NAME}_45S_repeat.gff 2>/dev/null
    barrnap --kingdom euk ${WD}/02_rdna/${SP_NAME}_5S_repeat.fa \
        > ${WD}/02_rdna/${SP_NAME}_5S_repeat.gff 2>/dev/null
    echo "  45S repeat: $(grep -c 'rRNA' ${WD}/02_rdna/${SP_NAME}_45S_repeat.gff || echo 0) rRNA hits"
    echo "  5S repeat:  $(grep -c 'rRNA' ${WD}/02_rdna/${SP_NAME}_5S_repeat.gff || echo 0) rRNA hits"

    # ── Step 3: Build collapsed rDNA reference ───────────────────────
    echo "[${SP_NAME}] Step 3: Build collapsed reference"
    bedtools maskfasta \
        -fi ${SP_FA} \
        -bed ${WD}/02_rdna/rDNA_arrays.bed \
        -fo ${WD}/04_ref/${SP_NAME}_T2T.rDNA_masked.fa
    cat ${WD}/04_ref/${SP_NAME}_T2T.rDNA_masked.fa \
        ${WD}/02_rdna/${SP_NAME}_45S_repeat.fa \
        ${WD}/02_rdna/${SP_NAME}_5S_repeat.fa \
        > ${WD}/04_ref/${SP_NAME}_CN_reference.fa
    samtools faidx ${WD}/04_ref/${SP_NAME}_CN_reference.fa
    bwa index ${WD}/04_ref/${SP_NAME}_CN_reference.fa 2> ${WD}/04_ref/bwa_index.log
    echo "  Reference: $(wc -c < ${WD}/04_ref/${SP_NAME}_CN_reference.fa) bp"

    # ── Step 4: BUSCO single-copy baseline ───────────────────────────
    echo "[${SP_NAME}] Step 4: BUSCO single-copy baseline"
    awk -F'\t' '
    BEGIN{OFS="\t"}
    !/^#/ && $2=="Complete" {
        s=$4; e=$5;
        if(s>e){tmp=s; s=e; e=tmp}
        print $3, s-1, e, $1
    }' ${SP_BUSCO} > ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.raw.bed

    bedtools intersect -v \
        -a ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.raw.bed \
        -b ${WD}/02_rdna/rDNA_arrays.bed \
        > ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.bed

    BUSCO_N=$(wc -l < ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.bed)
    echo "  Single-copy BUSCO loci: ${BUSCO_N}"

    # ── Step 5: Illumina WGS mapping ─────────────────────────────────
    echo "[${SP_NAME}] Step 5: Illumina WGS mapping (bwa mem, ${THREADS} threads)"
    bwa mem -t ${THREADS} -R "@RG\tID:${SP_NAME}\tSM:${SP_NAME}" \
        ${WD}/04_ref/${SP_NAME}_CN_reference.fa \
        ${SP_WGS_R1} ${SP_WGS_R2} \
        | samtools sort -@ 16 -o ${WD}/05_map/${SP_NAME}_WGS_CN.bam -
    samtools index ${WD}/05_map/${SP_NAME}_WGS_CN.bam
    samtools flagstat ${WD}/05_map/${SP_NAME}_WGS_CN.bam \
        > ${WD}/05_map/${SP_NAME}_WGS_CN.flagstat.txt
    echo "  WGS BAM done"

    # ── Step 6: rDNA component BED from representative repeats ───────
    echo "[${SP_NAME}] Step 6: rDNA component BED"
    {
        grep 'rRNA' ${WD}/02_rdna/${SP_NAME}_45S_repeat.gff | \
        awk -F'\t' '{
            type=""
            if($9 ~ /5_8S|5\.8S/) type="5.8S"
            else if($9 ~ /18S/) type="18S"
            else if($9 ~ /28S/) type="28S"
            if(type) print $1"\t"$4-1"\t"$5"\t"type
        }'
        grep 'rRNA' ${WD}/02_rdna/${SP_NAME}_5S_repeat.gff | \
        awk -F'\t' '{
            type=""
            if($9 ~ /5S/ && $9 !~ /5_8S|5\.8S/) type="5S"
            if(type) print $1"\t"$4-1"\t"$5"\t"type
        }'
    } > ${WD}/06_depth/${SP_NAME}_rDNA_components.bed

    N_comp=$(wc -l < ${WD}/06_depth/${SP_NAME}_rDNA_components.bed)
    echo "  rDNA components: ${N_comp} regions"

    # ── Step 7: mosdepth for rDNA + BUSCO (WGS) ─────────────────────
    echo "[${SP_NAME}] Step 7: mosdepth (WGS)"
    mosdepth -t ${MOSDEPTH_THREADS} --by ${WD}/06_depth/${SP_NAME}_rDNA_components.bed \
        ${WD}/06_depth/${SP_NAME}_WGS_rDNA \
        ${WD}/05_map/${SP_NAME}_WGS_CN.bam 2>/dev/null
    mosdepth -t ${MOSDEPTH_THREADS} --by ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.bed \
        ${WD}/06_depth/${SP_NAME}_WGS_BUSCO \
        ${WD}/05_map/${SP_NAME}_WGS_CN.bam 2>/dev/null
    echo "  mosdepth done"

    # ── Step 8: Calculate CN (WGS) ───────────────────────────────────
    echo "[${SP_NAME}] Step 8: Calculate CN (WGS)"
    python3 - "${SP_NAME}" "${WD}" "WGS" <<'PYEOF'
import sys, gzip, numpy as np, pandas as pd

def read_mosdepth(path):
    opener = gzip.open if path.endswith('.gz') else open
    rows = []
    with opener(path, 'rt') as f:
        for line in f:
            x = line.rstrip().split('\t')
            rows.append({'chrom':x[0],'start':int(x[1]),'end':int(x[2]),
                         'name':x[3] if len(x)>=5 else '.','depth':float(x[-1])})
    return pd.DataFrame(rows)

name, wd, tag = sys.argv[1], sys.argv[2], sys.argv[3]
busco = read_mosdepth(f'{wd}/06_depth/{name}_{tag}_BUSCO.regions.bed.gz')
rdna  = read_mosdepth(f'{wd}/06_depth/{name}_{tag}_rDNA.regions.bed.gz')
D_single = np.median(busco['depth'])

print(f"  D_single (BUSCO median) = {D_single:.3f}x")
for _,r in rdna.iterrows():
    cn = r['depth'] / D_single if D_single>0 else 0
    print(f"  {r['name']:6s}  depth={r['depth']:9.3f}  CN={cn:.2f}")

sub = rdna[rdna['name'].isin(['18S','5.8S','28S'])]
cn_45s = np.median(sub['depth']) / D_single if D_single>0 and len(sub)>0 else np.nan
sub5 = rdna[rdna['name']=='5S']
cn_5s = np.median(sub5['depth'])/D_single if D_single>0 and len(sub5)>0 else np.nan

with open(f'{wd}/07_cn/{name}_{tag}_CN_summary.tsv','w') as f:
    f.write("accession\ttype\testimated_CN\tsingle_copy_depth\n")
    f.write(f"{name}\t45S\t{cn_45s:.3f}\t{D_single:.3f}\n")
    f.write(f"{name}\t5S\t{cn_5s:.3f}\t{D_single:.3f}\n")
print(f"  => {name} {tag}: 45S={cn_45s:.2f}  5S={cn_5s:.2f}")
PYEOF

    # ── Step 9: HiFi mapping (independent validation) ────────────────
    echo "[${SP_NAME}] Step 9: HiFi mapping (minimap2, ${THREADS} threads)"
    minimap2 -ax map-hifi -t ${THREADS} -R "@RG\tID:${SP_NAME}_hifi\tSM:${SP_NAME}" \
        ${WD}/04_ref/${SP_NAME}_CN_reference.fa ${SP_HIFI} \
        | samtools sort -@ 16 -o ${WD}/05_map/${SP_NAME}_HiFi_CN.bam -
    samtools index ${WD}/05_map/${SP_NAME}_HiFi_CN.bam
    echo "  HiFi BAM done"

    # ── Step 10: mosdepth + CN for HiFi ─────────────────────────────
    echo "[${SP_NAME}] Step 10: CN (HiFi)"
    mosdepth -t ${MOSDEPTH_THREADS} --by ${WD}/06_depth/${SP_NAME}_rDNA_components.bed \
        ${WD}/06_depth/${SP_NAME}_HiFi_rDNA \
        ${WD}/05_map/${SP_NAME}_HiFi_CN.bam 2>/dev/null
    mosdepth -t ${MOSDEPTH_THREADS} --by ${WD}/03_busco/${SP_NAME}_BUSCO_singlecopy.bed \
        ${WD}/06_depth/${SP_NAME}_HiFi_BUSCO \
        ${WD}/05_map/${SP_NAME}_HiFi_CN.bam 2>/dev/null

    python3 - "${SP_NAME}" "${WD}" "HiFi" <<'PYEOF'
import sys, gzip, numpy as np, pandas as pd
def read_mosdepth(path):
    opener = gzip.open if path.endswith('.gz') else open
    rows = []
    with opener(path, 'rt') as f:
        for line in f:
            x = line.rstrip().split('\t')
            rows.append({'chrom':x[0],'start':int(x[1]),'end':int(x[2]),
                         'name':x[3] if len(x)>=5 else '.','depth':float(x[-1])})
    return pd.DataFrame(rows)
name, wd, tag = sys.argv[1], sys.argv[2], sys.argv[3]
busco = read_mosdepth(f'{wd}/06_depth/{name}_{tag}_BUSCO.regions.bed.gz')
rdna  = read_mosdepth(f'{wd}/06_depth/{name}_{tag}_rDNA.regions.bed.gz')
D_single = np.median(busco['depth'])
sub = rdna[rdna['name'].isin(['18S','5.8S','28S'])]
cn_45s = np.median(sub['depth'])/D_single if D_single>0 and len(sub)>0 else np.nan
sub5 = rdna[rdna['name']=='5S']
cn_5s = np.median(sub5['depth'])/D_single if D_single>0 and len(sub5)>0 else np.nan
with open(f'{wd}/07_cn/{name}_{tag}_CN_summary.tsv','w') as f:
    f.write("accession\ttype\testimated_CN\tsingle_copy_depth\n")
    f.write(f"{name}\t45S\t{cn_45s:.3f}\t{D_single:.3f}\n")
    f.write(f"{name}\t5S\t{cn_5s:.3f}\t{D_single:.3f}\n")
print(f"  => {name} {tag}: 45S={cn_45s:.2f}  5S={cn_5s:.2f}")
PYEOF

    echo "  ${SP_NAME} done."
}


# ==========================================================================
# Main
# ==========================================================================
echo "========================================"
echo " rDNA copy-number estimation pipeline"
echo " started at $(date)"
echo "========================================"

# ── Gifu ──────────────────────────────────────────────────────────────
process_one_species "Gifu" \
    "${GIFU_FA}" "${GIFU_WGS_R1}" "${GIFU_WGS_R2}" "${GIFU_HIFI}" "${GIFU_BUSCO}" \
    "${GIFU_45S_CHR}" ${GIFU_45S_START} ${GIFU_45S_END} \
    "${GIFU_5S_CHR}"  ${GIFU_5S_START}  ${GIFU_5S_END}

# ── MG20 ──────────────────────────────────────────────────────────────
process_one_species "MG20" \
    "${MG20_FA}" "${MG20_WGS_R1}" "${MG20_WGS_R2}" "${MG20_HIFI}" "${MG20_BUSCO}" \
    "${MG20_45S_CHR}" ${MG20_45S_START} ${MG20_45S_END} \
    "${MG20_5S_CHR}"  ${MG20_5S_START}  ${MG20_5S_END}


# ==========================================================================
# Final summary table
# ==========================================================================
echo ""
echo "========================================"
echo " Final summary"
echo "========================================"

SUMMARY=${OUTDIR}/final_rDNA_CN.tsv
cat > ${SUMMARY} <<'HEADER'
accession	type	method	assembly_CN	coverage_CN	single_copy_depth
HEADER

for sp in Gifu MG20; do
    for typ in 45S 5S; do
        # Assembly CN from previous run (rDNA_array_summary.tsv)
        case "${sp}_${typ}" in
            Gifu_45S)  asm_cn="40+16+10+76=142" ;;
            Gifu_5S)   asm_cn="1480" ;;
            MG20_45S)  asm_cn="129+19+38+20=206" ;;
            MG20_5S)   asm_cn="3491" ;;
        esac
        for tag in WGS HiFi; do
            f=${OUTDIR}/${sp}/07_cn/${sp}_${tag}_CN_summary.tsv
            if [ -f "${f}" ]; then
                val=$(awk -F'\t' -v t="${typ}" '$2==t{print $3}' ${f})
                d=$(awk -F'\t' -v t="${typ}" '$2==t{print $4}' ${f})
                echo "${sp}	${typ}	${tag}	${asm_cn}	${val}	${d}" >> ${SUMMARY}
            fi
        done
    done
done

echo ""
echo "Final table:  ${SUMMARY}"
cat ${SUMMARY}
echo ""
echo "========================================"
echo " Pipeline finished at $(date)"
echo "========================================"
```

</details>
