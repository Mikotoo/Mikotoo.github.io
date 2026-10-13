---
title: 'gap 区域分析'
sidebar:
  label: 'gap 区域分析'
  order: 999
---


## 中文

**作用**：比较已发表的 Gifu (v1.2) 与 MG20 (gnm3) 组装的 gap 数目与长度，并写出统计表与对比图。
**入口**：`gap_analysis.py`
**输入**：与 `gap_analysis/` 同级的 `Gifu/Gifu_gap.gff`、`Gifu/Lotusjaponicus_Gifu_v1.2_genome.fa.fai`、`MG20/Lotus_MG20_gap.gff`、`MG20/Lotus_japonicus.fasta.fai`（GFF3 中 type 为 `gap`，1-based 闭区间）；染色体长度取自 `.fai`。
**输出**：`results/`（`{Gifu,MG20}_gap_summary_per_chromosome.csv`、`{key}_gap_positions.csv`、`{key}_gap_size_distribution.csv`、`gap_genome_comparison.csv`）与 `figures/`（`{key}_gap_per_chromosome.png`、`{key}_gap_size_distribution.png`、`{key}_gap_map.png`、`comparison_genome_level.png`、`comparison_per_chromosome.png`、`comparison_size_ecdf.png`）。
**运行**：`python gap_analysis.py`（无绝对路径：数据、`results/`、`figures/` 均由脚本位置推导，`Gifu/`、`MG20/` 必须与 `gap_analysis/` 同级）。
**工具**：Python 3 + pandas、numpy、matplotlib（`Agg`）；字体优先 `WenQuanYi Micro Hei`，回退 `DejaVu Sans`；染色体排序正则 `(?:chr|Lj)(\d+|[A-Za-z]+)$`。

<details>
<summary>English</summary>

**Purpose**: Compare gap number and length between the published Gifu (v1.2) and MG20 (gnm3) assemblies and write the statistics tables and comparison plots.
**Entry point**: `gap_analysis.py`
**Inputs**: in the directory next to `gap_analysis/`: `Gifu/Gifu_gap.gff`, `Gifu/Lotusjaponicus_Gifu_v1.2_genome.fa.fai`, `MG20/Lotus_MG20_gap.gff`, `MG20/Lotus_japonicus.fasta.fai` (the GFF3 records with type `gap`, 1-based closed intervals); chromosome lengths come from the `.fai`.
**Outputs**: `results/` (`{Gifu,MG20}_gap_summary_per_chromosome.csv`, `{key}_gap_positions.csv`, `{key}_gap_size_distribution.csv`, `gap_genome_comparison.csv`) and `figures/` (`{key}_gap_per_chromosome.png`, `{key}_gap_size_distribution.png`, `{key}_gap_map.png`, `comparison_genome_level.png`, `comparison_per_chromosome.png`, `comparison_size_ecdf.png`).
**Run**: `python gap_analysis.py` (no absolute paths: the data, `results/` and `figures/` are derived from the script location, so `Gifu/` and `MG20/` must sit next to `gap_analysis/`).
**Tools**: Python 3 with pandas, numpy and matplotlib (`Agg`); font preference `WenQuanYi Micro Hei`, falling back to `DejaVu Sans`; chromosome order regex `(?:chr|Lj)(\d+|[A-Za-z]+)$`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details>
<summary><code>gap_analysis.py</code> · Python · 431 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/01_published_assemblies/gap_analysis/gap_analysis.py)

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gap_analysis.py — 分析两个百脉根(Lotus japonicus)基因组组装中的 gap 数目与长度，并绘制可视化图片。

输入（GFF3 格式的 gap 注释，坐标 1-based 闭合区间，长度 = end - start + 1）：
    Gifu: 04_compGenome/00.published_lotus/Gifu/Gifu_gap.gff
    MG20: 04_compGenome/00.published_lotus/MG20/Lotus_MG20_gap.gff
染色体长度取自对应 FASTA 的 .fai 索引文件。

输出：
    results/  : CSV 汇总表（每个基因组、每条染色体的 gap 统计；大小分布；gap 位置表）
    figures/  : 每个基因组单独的可视化 + 两个基因组对比图
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator

# ---------------------------------------------------------------------------
# 全局设置
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(BASE)          # 00.published_lotus
OUT_R = os.path.join(BASE, 'results')
OUT_F = os.path.join(BASE, 'figures')
os.makedirs(OUT_R, exist_ok=True)
os.makedirs(OUT_F, exist_ok=True)

plt.rcParams['font.family'] = ['WenQuanYi Micro Hei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 100

COLOR_GIFU = '#4C72B0'   # 蓝
COLOR_MG20 = '#DD8452'   # 橙

GENOMES = [
    dict(key='Gifu', label='Gifu (v1.2)',
         gff=os.path.join(DATA, 'Gifu', 'Gifu_gap.gff'),
         fai=os.path.join(DATA, 'Gifu', 'Lotusjaponicus_Gifu_v1.2_genome.fa.fai')),
    dict(key='MG20', label='MG20 (gnm3)',
         gff=os.path.join(DATA, 'MG20', 'Lotus_MG20_gap.gff'),
         fai=os.path.join(DATA, 'MG20', 'Lotus_japonicus.fasta.fai')),
]


# ---------------------------------------------------------------------------
# 解析
# ---------------------------------------------------------------------------
def parse_gap_gff(path):
    """解析 gap GFF，返回 DataFrame(seqid, start, end, length)。"""
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gap':
                continue
            seqid, start, end = p[0], int(p[3]), int(p[4])
            rows.append((seqid, start, end, end - start + 1))
    return pd.DataFrame(rows, columns=['seqid', 'start', 'end', 'length'])


def parse_fai(path):
    """解析 .fai 索引 -> {seqid: length}。"""
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            d[p[0]] = int(p[1])
    return d


def chrom_key(name):
    """染色体排序键：主要染色体 1..6 在前，然后是 0(unplaced)，再是细胞器。"""
    m = re_chr.search(name)
    if not m:
        return (9, 0, name)
    num = m.group(1)
    if num.isdigit():
        n = int(num)
        return (0, n, name) if n >= 1 else (1, 0, name)
    if num in ('C', 'M', 'chloro', 'mito'):
        return (2, 0, name)
    return (9, 0, name)


re_chr = re.compile(r'(?:chr|Lj)(\d+|[A-Za-z]+)$')


def summarize(gff_path, fai_path):
    """返回 (df_gaps, df_chr, genome_stats)。"""
    df = parse_gap_gff(gff_path)
    lens = parse_fai(fai_path)
    order = sorted([c for c in lens if c in set(df.seqid)], key=chrom_key)

    rows = []
    for c in order:
        sub = df[df.seqid == c]
        L = lens[c]
        n = len(sub)
        tot = int(sub.length.sum()) if n else 0
        rows.append(dict(
            chromosome=c, chr_length=L, n_gaps=n, total_gap_bp=tot,
            mean_gap_bp=round(sub.length.mean(), 2) if n else np.nan,
            median_gap_bp=float(sub.length.median()) if n else np.nan,
            max_gap_bp=int(sub.length.max()) if n else 0,
            min_gap_bp=int(sub.length.min()) if n else 0,
            gap_density_per_Mb=round(n / (L / 1e6), 3) if L else np.nan,
            gap_fraction_pct=round(tot / L * 100, 4) if L else np.nan,
        ))
    df_chr = pd.DataFrame(rows)

    n_all = int(df.shape[0])
    tot_all = int(df.length.sum())
    size_all = df.length.astype(float)
    genome_size = sum(lens.values())   # 完整组装大小（含无 gap 的序列）
    stats = dict(
        genome_size_bp=genome_size,
        n_gaps=n_all,
        total_gap_bp=tot_all,
        mean_gap_bp=round(size_all.mean(), 2) if n_all else np.nan,
        median_gap_bp=float(size_all.median()) if n_all else np.nan,
        max_gap_bp=int(size_all.max()) if n_all else 0,
        min_gap_bp=int(size_all.min()) if n_all else 0,
        gap_fraction_pct=round(tot_all / genome_size * 100, 4),
        gap_density_per_Mb=round(n_all / (genome_size / 1e6), 3),
        n_chromosomes_with_gap=len(order),
    )
    return df, df_chr, stats


# ---------------------------------------------------------------------------
# 绘图工具
# ---------------------------------------------------------------------------
def fmt_k(x, pos=None):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    if x >= 1e3:
        return f'{x/1e3:.0f}K'
    return f'{x:.0f}'


def fmt_bp(x, pos=None):
    if x >= 1e6:
        return f'{x/1e6:.2f} Mb'
    if x >= 1e3:
        return f'{x/1e3:.1f} kb'
    return f'{x:.0f} bp'


def bar_chr(ax, df_chr, col, ylabel, title, color, err=None, log=False):
    """按染色体顺序画柱状图，并在柱顶标注数值。"""
    x = np.arange(len(df_chr))
    labels = [short_name(c) for c in df_chr['chromosome']]
    vals = df_chr[col].astype(float)
    bars = ax.bar(x, vals, color=color, width=0.62, edgecolor='white')
    if log:
        ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.margins(y=0.12)
    for xi, v in zip(x, vals):
        if pd.isna(v):
            continue
        txt = fmt_k(v)
        ax.text(xi, v, txt, ha='center', va='bottom', fontsize=7.5, color='#333333')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return bars


def short_name(chrom):
    """把染色体名缩短用于图例：LjG1.1_chr1 -> chr1；lotja.MG20.gnm3.Lj1 -> Lj1"""
    m = re.search(r'(chr\d+|Lj\d+|Lj[A-Za-z]+)$', chrom)
    return m.group(1) if m else chrom


def gap_map(ax, df, chrlens, color, title, max_mark=8000):
    """染色体 ideogram：每条染色体画一条按长度缩放的线，gap 位置打点。"""
    order = sorted([c for c in chrlens if c in set(df.seqid)], key=chrom_key)
    y = 0
    yticks, ylabels = [], []
    for c in order:
        L = chrlens[c]
        sub = df[df.seqid == c]
        ax.plot([0, L], [y, y], color='#BBBBBB', lw=8, solid_capstyle='butt', zorder=1)
        if len(sub):
            sizes = sub['length'].values
            pos = sub['start'].values
            # 小 gap 用细点，>=1000 的用粗点/深色
            small = sizes < 1000
            if small.any():
                ax.scatter(pos[small], np.full(small.sum(), y), s=2.5,
                           color=color, alpha=0.75, zorder=3, linewidths=0)
            big = ~small
            if big.any():
                ax.scatter(pos[big], np.full(big.sum(), y), s=14,
                           color='#B22222', alpha=0.9, zorder=4,
                           marker='|', linewidths=1.2)
        ax.text(-L * 0.015, y, short_name(c), ha='right', va='center', fontsize=9)
        y += 1
        yticks.append(y - 0.5)
        ylabels.append(f'{L/1e6:.1f} Mb')
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel('染色体位置 (bp)')
    ax.set_ylabel('染色体 (长度)')
    ax.set_title(title, fontsize=11)
    ax.set_ylim(-0.6, y - 0.4)
    ax.margins(x=0.02)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # 图例
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=4, label='gap < 1 kb'),
        Line2D([0], [0], marker='|', color='none', markeredgecolor='#B22222', markersize=10, markeredgewidth=1.5, label='gap ≥ 1 kb'),
    ]
    ax.legend(handles=handles, loc='lower right', fontsize=8, frameon=False)


def size_distribution(ax_hist, ax_ecdf, lengths, color, label, n, tot_bp):
    """gap 长度分布：log 直方图 + ECDF。"""
    l = np.asarray(lengths, dtype=float)
    # 直方图（log10 分箱）
    bins = np.logspace(np.log10(max(l.min(), 1)), np.log10(l.max()), 40)
    ax_hist.hist(l, bins=bins, color=color, alpha=0.85, edgecolor='white', linewidth=0.3)
    ax_hist.set_xscale('log')
    ax_hist.set_xlabel('gap 长度 (bp, log)')
    ax_hist.set_ylabel('gap 数目')
    ax_hist.set_title(f'{label}：gap 长度分布', fontsize=11)
    ax_hist.text(0.97, 0.95, f'n = {n:,}\n总长 = {tot_bp:,} bp\n平均 = {l.mean():.1f} bp\n中位数 = {np.median(l):.0f} bp',
                 transform=ax_hist.transAxes, ha='right', va='top', fontsize=9,
                 bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='0.7'))
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)
    # ECDF
    srt = np.sort(l)
    ecdf = np.arange(1, len(srt) + 1) / len(srt)
    ax_ecdf.plot(srt, ecdf, color=color, lw=1.8)
    ax_ecdf.set_xscale('log')
    ax_ecdf.set_xlabel('gap 长度 (bp, log)')
    ax_ecdf.set_ylabel('累积比例 (ECDF)')
    ax_ecdf.set_title(f'{label}：gap 长度累积分布', fontsize=11)
    ax_ecdf.set_ylim(0, 1.02)
    ax_ecdf.spines['top'].set_visible(False)
    ax_ecdf.spines['right'].set_visible(False)


def per_genome_figures(key, label, df, df_chr, stats):
    color = COLOR_GIFU if key == 'Gifu' else COLOR_MG20
    chrlens = dict(zip(df_chr['chromosome'], df_chr['chr_length']))

    # ---- 每条染色体的 gap 统计 ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    bar_chr(axes[0, 0], df_chr, 'n_gaps', 'gap 数目', f'{label}：各染色体 gap 数目',
            color=color)
    bar_chr(axes[0, 1], df_chr, 'total_gap_bp', '总长度 (bp)', f'{label}：各染色体 gap 总长度',
            color=color)
    bar_chr(axes[1, 0], df_chr, 'gap_density_per_Mb', 'gap 数 / Mb', f'{label}：各染色体 gap 密度',
            color=color)
    bar_chr(axes[1, 1], df_chr, 'gap_fraction_pct', 'gap 占比 (%)', f'{label}：各染色体 gap 占比',
            color=color)
    fig.suptitle(f'{label} — gap 统计（全基因组共 {stats["n_gaps"]:,} 个 gap，总长 {stats["total_gap_bp"]:,} bp，'
                 f'占基因组 {stats["gap_fraction_pct"]:.3f}%）', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(OUT_F, f'{key}_gap_per_chromosome.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ---- gap 长度分布 ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    size_distribution(axes[0], axes[1], df['length'], color, label,
                      stats['n_gaps'], stats['total_gap_bp'])
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_F, f'{key}_gap_size_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ---- 染色体 gap 位置图 ----
    fig, ax = plt.subplots(figsize=(12, max(4, 0.9 * len(df_chr))))
    gap_map(ax, df, chrlens, color, f'{label}：gap 在染色体上的分布')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_F, f'{key}_gap_map.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def comparison_figures(all_dfs, all_chrs, all_stats):
    keys = [g['key'] for g in GENOMES]
    colors = [COLOR_GIFU, COLOR_MG20]
    labels = [g['label'] for g in GENOMES]

    # ---- 基因组水平对比 ----
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    # 总数
    ax = axes[0]
    x = np.arange(2)
    counts = [all_stats[k]['n_gaps'] for k in keys]
    bars = ax.bar(x, counts, color=colors, width=0.55, edgecolor='white')
    for xi, v in zip(x, counts):
        ax.text(xi, v, f'{v:,}', ha='center', va='bottom', fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel('gap 数目')
    ax.set_title('gap 总数', fontsize=11)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    # 总长度
    ax = axes[1]
    tots = [all_stats[k]['total_gap_bp'] for k in keys]
    bars = ax.bar(x, tots, color=colors, width=0.55, edgecolor='white')
    for xi, v in zip(x, tots):
        ax.text(xi, v, fmt_bp(v), ha='center', va='bottom', fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel('总长度 (bp)')
    ax.set_title('gap 总长度', fontsize=11)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    # 占比 + 密度
    ax = axes[2]
    fracs = [all_stats[k]['gap_fraction_pct'] for k in keys]
    dens = [all_stats[k]['gap_density_per_Mb'] for k in keys]
    bars = ax.bar(x - 0.18, fracs, width=0.36, color=colors, edgecolor='white', label='gap 占比 (%)')
    for xi, v in zip(x - 0.18, fracs):
        ax.text(xi, v, f'{v:.3f}%', ha='center', va='bottom', fontsize=9)
    bars = ax.bar(x + 0.18, dens, width=0.36, color=[c + '99' for c in colors],
                  edgecolor='white', label='gap 密度 (/Mb)')
    for xi, v in zip(x + 0.18, dens):
        ax.text(xi, v, f'{v:.1f}', ha='center', va='bottom', fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel('数值')
    ax.set_title('gap 占基因组比例 / 密度', fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.suptitle('Gifu vs MG20：基因组 gap 总体对比', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(OUT_F, 'comparison_genome_level.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ---- 按染色体 rank 对比 gap 密度 ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, col, ylab, ttl in [
        (axes[0], 'n_gaps', 'gap 数目', '各染色体 gap 数目对比'),
        (axes[1], 'gap_fraction_pct', 'gap 占比 (%)', '各染色体 gap 占比对比'),
    ]:
        for k, c, lab in zip(keys, colors, labels):
            dc = all_chrs[k]
            ax.bar(np.arange(len(dc)) + (0.0 if k == keys[0] else 0.45),
                   dc[col].astype(float), width=0.45, color=c, label=lab, edgecolor='white')
        ax.set_xticks(np.arange(len(all_chrs[keys[0]])) + 0.225)
        ax.set_xticklabels([short_name(c) for c in all_chrs[keys[0]]['chromosome']], rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(ylab); ax.set_title(ttl, fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_F, 'comparison_per_chromosome.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ---- 长度分布对比（ECDF）----
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for k, c, lab in zip(keys, colors, labels):
        l = np.sort(all_dfs[k]['length'].astype(float).values)
        ecdf = np.arange(1, len(l) + 1) / len(l)
        ax.plot(l, ecdf, color=c, lw=2, label=f"{lab} (n={len(l):,}, 总长 {all_stats[k]['total_gap_bp']:,} bp)")
    ax.set_xscale('log')
    ax.set_xlabel('gap 长度 (bp, log)')
    ax.set_ylabel('累积比例 (ECDF)')
    ax.set_title('Gifu vs MG20：gap 长度累积分布对比', fontsize=12)
    ax.legend(fontsize=9, frameon=False)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25, which='both')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_F, 'comparison_size_ecdf.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    all_dfs, all_chrs, all_stats = {}, {}, {}

    for g in GENOMES:
        key, label = g['key'], g['label']
        print(f'[1/2] 分析 {label} ...')
        df, df_chr, stats = summarize(g['gff'], g['fai'])
        all_dfs[key], all_chrs[key], all_stats[key] = df, df_chr, stats

        # 保存结果表
        df_chr.to_csv(os.path.join(OUT_R, f'{key}_gap_summary_per_chromosome.csv'),
                      index=False, encoding='utf-8-sig')
        df.to_csv(os.path.join(OUT_R, f'{key}_gap_positions.csv'),
                  index=False, encoding='utf-8-sig')

        # 大小分布（分箱）表
        bins = [0, 10, 100, 1000, 10000, np.inf]
        binlab = ['1-10', '11-100', '101-1000', '1001-10000', '>10000']
        s = pd.cut(df['length'], bins=bins, labels=binlab)
        dist = (s.value_counts().reindex(binlab).rename_axis('size_class')
                .to_frame('n_gaps'))
        dist['total_bp'] = dist.index.map(
            lambda b: int(df.loc[s == b, 'length'].sum()))
        dist.to_csv(os.path.join(OUT_R, f'{key}_gap_size_distribution.csv'),
                    encoding='utf-8-sig')

        print(f'    gap 数 = {stats["n_gaps"]:,}，总长 = {stats["total_gap_bp"]:,} bp，'
              f'平均 {stats["mean_gap_bp"]} bp，中位数 {stats["median_gap_bp"]} bp，'
              f'最大 {stats["max_gap_bp"]:,} bp')
        print(f'    基因组 {stats["genome_size_bp"]:,} bp，gap 占比 {stats["gap_fraction_pct"]:.4f}%，'
              f'密度 {stats["gap_density_per_Mb"]} 个/Mb')

        print(f'[2/2] 绘制 {label} 图片 ...')
        per_genome_figures(key, label, df, df_chr, stats)

    # 基因组水平对比表
    comp = pd.DataFrame(all_stats).T
    comp.insert(0, 'genome', [g['label'] for g in GENOMES])
    comp.to_csv(os.path.join(OUT_R, 'gap_genome_comparison.csv'), encoding='utf-8-sig')

    print('绘制对比图 ...')
    comparison_figures(all_dfs, all_chrs, all_stats)
    print('全部完成。')


if __name__ == '__main__':
    main()
```

</details>
