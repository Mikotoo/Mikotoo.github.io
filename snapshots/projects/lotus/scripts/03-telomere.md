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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>run_telomere.sh</code> · Shell · 104 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/03_telomere/run_telomere.sh)

```bash
#!/bin/bash
#CSUB -J telomere
#CSUB -q c01
#CSUB -o run_telomere.out
#CSUB -e run_telomere.error
#CSUB -n 32
#CSUB -R span[hosts=1]
#
# Gifu telomere analysis — single entry point.
#
#   csub run_telomere.sh                    # full pipeline (compute + one plot pass)
#   TELOMERE_STEPS=plot bash run_telomere.sh   # re-plot only
#   TELOMERE_STEPS=compute bash run_telomere.sh  # find-tel + coverage + alignments
#
# Python: telomere.py  (subcommands: find-tel | coverage | alignments | plot | all)

set -euo pipefail

PROJECT_ROOT="{PROJ_04_LOTUS_GENOME}"
# csub copies the script to ~/.cbsched; BASH_SOURCE then points there, not at telomere.py.
SCRIPT_DIR="${PROJECT_ROOT}/03_annotation/Gifu/01.tel"
WORKDIR="${TELOMERE_WORKDIR:-${SCRIPT_DIR}}"
PREFIX="${TELOMERE_PREFIX:-GifuT2T}"
STEPS="${TELOMERE_STEPS:-all}"

PYTHON="${TELOMERE_PYTHON:-{CONDA_PY3}/bin/python3}"
SAMTOOLS="${TELOMERE_SAMTOOLS:-{SOFTWARE_ZC}/TeloComp-1.0.0/Dependencies/samtools-1.18/samtools}"

FASTA_ASM="${TELOMERE_FASTA_ASM:-${PROJECT_ROOT}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta}"
FASTA_PLOT="${TELOMERE_FASTA_PLOT:-${PROJECT_ROOT}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta}"
ONT_BAM="${TELOMERE_ONT_BAM:-${PROJECT_ROOT}/02_assembly/Gifu/08.quality/01.mapping/01.ont/ont_q30.bam}"
HIFI_BAM="${TELOMERE_HIFI_BAM:-${PROJECT_ROOT}/02_assembly/Gifu/08.quality/01.mapping/02.hifi/ccs_q30.bam}"
BEDTOOLS="${TELOMERE_BEDTOOLS:-{CONDA_PY3}/bin/bedtools}"

BIN_SIZE="${TELOMERE_BIN_SIZE:-100}"
PLOT_FLANK="${TELOMERE_PLOT_FLANK:-10000}"
THREADS="${TELOMERE_THREADS:-16}"

BASE_ARGS=(
  --workdir "${WORKDIR}"
  --prefix "${PREFIX}"
  --bin-size "${BIN_SIZE}"
  --plot-flank "${PLOT_FLANK}"
  --fasta-asm "${FASTA_ASM}"
  --fasta-plot "${FASTA_PLOT}"
  --ont-bam "${ONT_BAM}"
  --hifi-bam "${HIFI_BAM}"
  --bedtools "${BEDTOOLS}"
)

index_bam() {
  local bam="$1"
  if [[ ! -f "${bam}.bai" ]]; then
    echo "Indexing ${bam} ..."
    "${SAMTOOLS}" index -@"${THREADS}" "${bam}"
  fi
}

cd "${WORKDIR}"
echo "WORKDIR : ${WORKDIR}"
echo "PREFIX  : ${PREFIX}"
echo "STEPS   : ${STEPS}"

case "${STEPS,,}" in
  all)
    if [[ ! -f "${FASTA_ASM}.fai" ]]; then
      echo "Building FAI: ${FASTA_ASM}.fai"
      "${SAMTOOLS}" faidx "${FASTA_ASM}"
    fi
    index_bam "${ONT_BAM}"
    index_bam "${HIFI_BAM}"
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" all "${BASE_ARGS[@]}" --reuse-units
    ;;
  compute)
    index_bam "${ONT_BAM}"
    index_bam "${HIFI_BAM}"
    for step in find-tel coverage alignments; do
      "${PYTHON}" "${SCRIPT_DIR}/telomere.py" "${step}" "${BASE_ARGS[@]}" --reuse-units
    done
    ;;
  plot)
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" plot "${BASE_ARGS[@]}"
    ;;
  find-tel|find_tel|coverage|alignments|align|reads)
    step="${STEPS,,}"
    step="${step//_/-}"
    [[ "${step}" == "find_tel" ]] && step="find-tel"
    [[ "${step}" == "align" || "${step}" == "reads" ]] && step="alignments"
    if [[ "${step}" != "find-tel" ]]; then
      index_bam "${ONT_BAM}"
      index_bam "${HIFI_BAM}"
    fi
    extra=()
    [[ "${step}" == "find-tel" ]] && extra=(--reuse-units)
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" "${step}" "${BASE_ARGS[@]}" "${extra[@]}"
    ;;
  *)
    echo "ERROR: unknown TELOMERE_STEPS=${STEPS}" >&2
    echo "  Use: all | compute | plot | find-tel | coverage | alignments" >&2
    exit 1
    ;;
esac

echo "Done. Tables: ${WORKDIR}/${PREFIX}.*  Figures: ${WORKDIR}/plots/"
```

</details>

<details>
<summary><code>telomere.py</code> · Python · 2714 行（内嵌前 120 行）</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/03_telomere/telomere.py)

```python
#!/usr/bin/env python3
"""
Gifu telomere analysis (motif calling, coverage, spanning reads, plots).

Use subcommands: find-tel, coverage, alignments, plot, all

Steps:
  1. Scan 7-bp telomere units (TTTAGGG / CCCTAAA, <=1 mismatch) on the plus strand.
  2. Compute 100-bp window density as sequence coverage (merged unit spans / bin width).
  3. Call telomere bins from each chromosome end inward; bridge up to 5 consecutive
     low-density bins when high-density bins continue on both sides.
  4. Refine coordinates from the outermost and innermost telomere bins: walk
     consecutive units from each bin edge toward the chromosome end; final span
     is the union of both walks (gaps in the middle are ignored).
  5. Write coordinates TSV; optional PDF figures if matplotlib is available.

Units are deduplicated by (chrom, start) before density and coordinate calling.
"""

from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
import argparse
import bisect
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

MOTIFS = ("CCCTAAA", "TTTAGGG")
BASES = "ATGC"
UNIT_LEN = 7
DEFAULT_BIN = 100
DEFAULT_DENSITY = 0.5
DEFAULT_MAX_UNIT_GAP = 14
DEFAULT_MAX_DIP_BINS = 5
DEFAULT_PLOT_FLANK = 10_000
MAX_DENSITY = 1.0
MOTIF_COLORS = {"CCCTAAA": "#2E86AB", "TTTAGGG": "#C73E1D"}


@dataclass
class TelUnit:
    chrom: str
    start: int
    end: int
    motif: str


@dataclass
class BinRow:
    chrom: str
    bin_start: int
    bin_end: int
    count: int
    density: float


@dataclass
class TelomereCall:
    chrom: str
    end: str
    bin_start: int
    bin_end: int
    tel_start: int
    tel_end: int
    length_bp: int
    dominant_motif: str
    n_units: int
    mean_bin_density: float
    max_bin_density: float


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    chrom: Optional[str] = None
    chunks: List[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if chrom is not None:
            yield chrom, "".join(chunks).upper()


def scan_telomere_units(fasta_path: Path) -> List[TelUnit]:
    units: List[TelUnit] = []
    for chrom, seq in parse_fasta(fasta_path):
        n = len(seq)
        if n < UNIT_LEN:
            continue
        for motif in MOTIFS:
            for i in range(n - UNIT_LEN + 1):
                for j in range(UNIT_LEN):
                    for base in BASES:
                        variant = motif[:j] + base + motif[j + 1 :]
                        if seq[i : i + UNIT_LEN] == variant:
                            units.append(TelUnit(chrom, i, i + UNIT_LEN, variant))
    return dedupe_units(units)


def dedupe_units(units: Sequence[TelUnit]) -> List[TelUnit]:
    """One record per (chrom, start); prefer exact TTTAGGG / CCCTAAA motifs."""
    best: Dict[Tuple[str, int], TelUnit] = {}
    for u in units:
        key = (u.chrom, u.start)
        cur = best.get(key)
```

_完整文件共 2714 行，其余部分请在 [GitHub](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/03_telomere/telomere.py) 上查看。_

</details>
