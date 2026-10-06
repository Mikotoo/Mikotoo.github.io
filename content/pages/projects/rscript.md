---
title: 多组学分析与论文图表
description: Lotus T2T 项目中的结构变异、基因对应、转录组与单核分析代码。
template: splash
prev: false
next: false
hero:
  tagline: 把基因组间的差异，与基因和细胞中的表达联系起来。
  actions:
    - text: 查看分析细节
      link: /_projects/rscript/
      icon: right-arrow
    - text: 单核分析
      link: /_projects/rscript/singlecell-rnaseq-analysis/
      variant: minimal
---

<div class="about-grid">
  <div><p class="eyebrow">PROJECT CONTEXT</p><h2>连接不同类型的数据</h2><p>这组分析围绕 Lotus T2T 项目展开，覆盖 Gifu、MG20 新旧基因组间的结构变异和基因对应，并结合 bulk 与单核 RNA-seq 探索表达特征。</p></div>
  <div><p class="eyebrow">SCOPE</p><h2>聚焦图表背后的计算</h2><p>收录与论文 Fig. 2–4 及相关补充数据对应的主要分析脚本，涵盖 Shell、Python 和 R。基因组组装部分可在 <a href="/projects/lotus/">Lotus T2T 项目</a> 中继续阅读。</p></div>
</div>

<div class="section-heading"><p class="eyebrow">CONNECTED ANALYSES</p><h2>四个相互连接的分析方向</h2></div>
<div class="code-grid">
  <article class="code-card"><p class="eyebrow">GENOME STRUCTURE</p><h3>结构变异与 PAV</h3><p>通过全基因组比对、SyRI 与未比对区段分析，比较材料和组装版本之间的差异。</p><a class="text-link" href="/_projects/rscript/">阅读方法 →</a></article>
  <article class="code-card"><p class="eyebrow">GENE CORRESPONDENCE</p><h3>同源对应与新增基因</h3><p>建立不同版本之间的基因对应，区分补洞、注释修正与新注释的基因。</p><a class="text-link" href="/_projects/rscript/">阅读方法 →</a></article>
  <article class="code-card"><p class="eyebrow">EXPRESSION NETWORKS</p><h3>表达模式与共表达</h3><p>结合 WGCNA、Mfuzz 与表达相关性分析，观察基因表达的组织特征。</p><a class="text-link" href="/_projects/lotus/scripts/07-transcriptome/">查看表达分析 →</a></article>
  <article class="code-card"><p class="eyebrow">SINGLE NUCLEUS</p><h3>细胞类型与基因网络</h3><p>单核数据质控、聚类、注释，以及基于 scTenifoldKnk 的下游分析。</p><a class="text-link" href="/_projects/rscript/singlecell-rnaseq-analysis/">查看单核分析 →</a></article>
</div>
