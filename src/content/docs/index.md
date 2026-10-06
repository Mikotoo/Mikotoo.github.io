---
title: 学习，学习，学习
description: 植物基因组与单细胞分析的技术笔记、分析流程与论文复现
template: splash
hero:
  tagline: 植物基因组 · 单细胞与空间转录组 · 分析流程与复现笔记
  actions:
    - text: 站点说明
      link: /intro/
      icon: right-arrow
    - text: Lotus T2T 分析流程
      link: /lotus/
      variant: minimal
---

<div class="hero-cards">
  <a href="/lotus/">
    <strong>Lotus T2T 分析流程</strong>
    <span>Gifu / MG20 端粒到端粒基因组：组装、注释、端粒、rDNA、着丝粒、比较基因组、转录组、单核</span>
  </a>
  <a href="/rscript/">
    <strong>论文分析脚本（Rscript）</strong>
    <span>结构变异与 PAV、同源基因、新基因分类、bulk RNA-seq、单核 RNA-seq 与 scTenifoldKnk</span>
  </a>
  <a href="/notes/soybean-snrna/">
    <strong>笔记与复现</strong>
    <span>基因组获取、数据下载、sRNA / RNA-seq / BS-seq 流程、大豆根瘤单细胞 + 空转复现</span>
  </a>
  <a href="/index-all/">
    <strong>内容索引</strong>
    <span>本站全部页面清单，按来源分组</span>
  </a>
</div>

## 这个站点是什么

两类内容，分开维护：

- **分析流程** —— 直接来自 `Lotus_genome` 与 `Rscript` 两个代码仓库的 `README.md`，构建时自动抓取。改流程说明时只改代码仓库，不存在「网站上一份、仓库里一份」的问题。
- **笔记与复现** —— 手写在 `_posts/` 下的 Markdown，构建时自动转成站点页面，配图自动本地化。

技术栈是 [Astro](https://astro.build/) + [Starlight](https://starlight.astro.build/)，静态构建、自带全文搜索，部署在 GitHub Pages。细节见[站点说明](/intro/)。
