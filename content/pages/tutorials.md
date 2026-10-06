---
title: 生物信息学基础教程
description: 从 Linux、参考基因组文件，到 Python 与 R 数据可视化，再到 RNA-seq 序列比对与定量，循序学习生物信息学的基础操作。
template: splash
prev: false
next: false
hero:
  tagline: 从读懂一份基因组文件，到把测序读长放回基因组。跟着实例，动手理解每一步。
  actions:
    - text: 从第一节开始
      link: /tutorials/01-linux-genome/
      icon: right-arrow
    - text: RNA-seq 与序列比对
      link: /tutorials/03-rnaseq-alignment/
      variant: minimal
---

<div class="section-heading"><p class="eyebrow">LEARN BY DOING</p><h2>把基础操作，放进真实的数据里</h2><p>适合刚接触生物信息学、希望熟悉命令行、绘图与测序数据分析的学习者。以大豆参考基因组为起点，把文件查看、结果统计、可视化与 RNA-seq 定量连成一条学习路径。</p></div>

<div class="topic-grid">
  <a class="topic-card" href="/tutorials/01-linux-genome/"><strong>读懂数据</strong><span>认识 FASTA、GFF3 与基因结构，用 Linux 命令查看和统计文件。</span></a>
  <a class="topic-card" href="/tutorials/02-python-r-visualization/"><strong>表达结果</strong><span>用 Python 和 R 把基因组统计表转化为清晰的图形。</span></a>
  <a class="topic-card" href="/tutorials/03-rnaseq-alignment/"><strong>比对与定量</strong><span>先跑通 FASTQ → 比对 → 定量，再读 SAM/BAM 字段与比对原理。</span></a>
  <a class="topic-card" href="/tutorials/bioinformatics/lotus-circos-data.zip"><strong>动手练习</strong><span>下载百脉根圈图数据与 RNA-seq 教学数据，跟着正文完成练习。</span></a>
</div>

<div class="section-heading"><p class="eyebrow">COURSE CHAPTERS</p><h2>课程目录</h2><p>三节内容前后衔接：第三节的绘图沿用第二节的工具，比对则用第一节下载的参考基因组。已熟悉前面内容的读者可以按需跳读。</p></div>

<!-- AUTO:LESSON_CARDS -->

<div class="about-grid">
  <div><p class="eyebrow">BEFORE YOU START</p><h2>准备一个练习环境</h2><p>第一节与第三节使用 Linux 命令行，可在 Linux、服务器或 WSL 中练习。第三节需要安装 <code>hisat2</code>、<code>samtools</code>、<code>stringtie</code> 等工具，正文给出 conda 安装命令；大豆 FASTA 与 GFF3 的获取方式在第一节正文中介绍。</p></div>
  <div><p class="eyebrow">FROM TABLES TO FIGURES</p><h2>让三节内容衔接起来</h2><p>第二节包含 Python 和 R 的环境安装与绘图代码，基础图形使用第一节生成的统计表，进阶圈图使用单独的 Lotus 数据包。第三节提供一套教学模拟数据，可在完全离线的情况下跑通「比对 → BAM → 定量」。</p><a class="text-link" href="/tutorials/bioinformatics/rnaseq-demo-data.zip" download>下载 RNA-seq 教学数据 · ZIP ↓</a></div>
</div>

<div class="profile-strip"><div><h2>学完基础，再看一个完整分析</h2><p>继续阅读单细胞、转录组与基因组分析的项目和笔记。</p></div><a class="text-link" href="/notes/">浏览分析笔记 →</a></div>
