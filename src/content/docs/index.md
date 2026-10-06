---
title: 学习，学习，学习
description: 个人文档与代码站：分析流程、论文复现、生信工具与踩坑笔记
template: splash
hero:
  tagline: 个人文档与代码站 · 分析流程 · 论文复现 · 工具笔记
  actions:
    - text: 项目文档
      link: /projects/
      icon: right-arrow
    - text: 笔记与复现
      link: /notes/soybean-snrna/
      variant: minimal
    - text: 站点说明
      link: /intro/
      variant: minimal
---

<div class="hero-cards">
  <a href="/projects/">
    <strong>项目文档</strong>
    <span>已收录的代码仓库：分析流程、参数与步骤说明，构建时直接从仓库 README 同步</span>
  </a>
  <a href="/index-all/">
    <strong>内容索引</strong>
    <span>本站全部页面清单，按项目与笔记分组</span>
  </a>
  <a href="/code/">
    <strong>代码与资源</strong>
    <span>随文脚本、分析图件与可下载资料</span>
  </a>
</div>

## 这里放什么

一个持续更新的个人文档 + 代码站，目前的内容分三类：

- **项目文档** —— 每个项目的分析流程与步骤说明。这些内容**不在本站维护**，而是构建时从对应代码仓库的 `README.md` 自动同步：仓库里改一句，站点跟着变，不会出现两份说法。
- **笔记与复现** —— 手写的教程、流程记录、论文复现与踩坑。新增一篇就是新建一个 Markdown 文件，不需要动任何配置。
- **代码与资源** —— 随文脚本、图件，以及可下载的资料。

技术栈是 [Astro](https://astro.build/) + [Starlight](https://starlight.astro.build/)，静态构建、自带全文搜索（<kbd>Ctrl</kbd> / <kbd>⌘</kbd> + <kbd>K</kbd>），部署在 GitHub Pages。想了解内容怎么组织、怎么加新内容，看[站点说明](/intro/)。
