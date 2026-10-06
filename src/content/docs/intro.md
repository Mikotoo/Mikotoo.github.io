---
title: 站点说明
description: 这个站点是什么、内容怎么组织、为什么这样分工
---

import { CardGrid, LinkCard } from '@astrojs/starlight/components';

这里放两类东西：**能跑的分析流程**（跟着代码走，和 git 仓库同步）和**读过的笔记**（自己写的心得与论文复现）。两者刻意分开，因为它们的维护方式完全不同。

<CardGrid>
  <LinkCard
    title="Lotus T2T 分析流程"
    href="/lotus/"
    description="Gifu / MG20 端粒到端粒基因组：组装、注释、端粒、rDNA、着丝粒、比较基因组、转录组、单核"
  />
  <LinkCard
    title="论文分析脚本（Rscript）"
    href="/rscript/"
    description="结构变异与 PAV、同源基因、新基因分类、bulk RNA-seq、单核 RNA-seq 与 scTenifoldKnk"
  />
  <LinkCard
    title="笔记与复现"
    href="/notes/genome/"
    description="基因组获取、数据下载、sRNA / RNA-seq / BS-seq 流程、大豆根瘤单细胞+空转复现"
  />
  <LinkCard
    title="代码与资源"
    href="/code/"
    description="随文脚本片段、图件，以及可供下载的资料"
  />
</CardGrid>

## 为什么这样组织

| 内容 | 放在哪 | 谁说了算 |
|---|---|---|
| 分析流程与参数 | `Lotus_genome/`、`Rscript/` 两个代码仓库 | **代码仓库**。每个步骤目录的 `README.md` 是唯一事实来源，站点在构建时自动抓取，不手工复制 |
| 笔记、教程、论文复现 | 本仓库 `_posts/` 目录 | **这里的 Markdown 文件**。构建时自动转成站点页面 |
| 站点外观与导航 | 本仓库 `astro.config.mjs`、`src/` | 本仓库 |

这样做的好处是：改流程说明时只需要改代码仓库里的 README，不会出现「网站上一份、仓库里一份」的不一致；写笔记时也只需要新建一个 Markdown 文件，不碰任何模板代码。

## 我是怎么更新的

```bash
# 1. 流程说明有更新时，先拉取代码仓库
git -C ../Lotus_genome pull
git -C ../Rscript pull

# 2. 本地预览（构建前会自动抓取 README 和迁移笔记）
pnpm dev            # http://localhost:4321

# 3. 提交博客内容
git add -A && git commit -m "update notes" && git push
```

`pnpm dev` 与 `pnpm build` 都会先执行两个生成脚本：

- `scripts/sync-docs.mjs` —— 抓取 `../Lotus_genome` 与 `../Rscript` 里所有 `README.md`，生成「分析流程」文档区（相对链接改写为 GitHub 链接，图片指向 raw 地址）。
- `scripts/migrate-posts.mjs` —— 把 `_posts/*.md` 转成站点笔记页，并把正文里引用本站图片的远程地址换成本地资源，交给 Astro 压缩优化。

生成结果写在 `src/content/docs/` 下，已在 `.gitignore` 中忽略——**它们每次构建重新生成**，不需要提交。

## 写作约定

## 写作约定

### 写一篇笔记

在 `_posts/` 新建 `YYYY-MM-DD-标题.md`，格式与旧站一致：

```markdown
---
title: 文章标题
date: 2026-01-01
category: 笔记
---

正文……
```

构建时它会自动出现在「笔记与复现」里。图片放在 `downloads/image/` 下，用本仓库的 GitHub 原始地址引用即可，构建时会自动本地化并压缩。

### 更新流程说明

直接改对应代码仓库里的 `README.md`，然后重新构建。站点导航按目录结构自动生成：

- `Lotus_genome/scripts/01_assembly/02_purge_haplotigs/README.md` → `/lotus/scripts/01-assembly/02-purge-haplotigs/`
- `Lotus_genome/scripts/01_assembly/README.md` → `/lotus/scripts/01-assembly/`
- `Rscript/README.md` → `/rscript/`

### 调整导航或外观

导航在 `astro.config.mjs` 的 `sidebar` 里，配色与排版在 `src/styles/custom.css` 里。新增笔记后需要同时在 `sidebar` 的「笔记与复现」里补一行。

## 一些说明

- **双语 README**：上游 README 是中文在前、英文在后。站点默认展开中文，英文收进可折叠区，需要时点开。
- **大文件不进站点**：仓库历史里提交过单细胞分析缓存（`.h5ad`，合计约 220 MB）。构建后的 `scripts/prune-dist.mjs` 会把它们从站点产物中剔除；这些文件仍在 GitHub 上，可按原始链接访问。
- **全文搜索**：Starlight 自带，快捷键 <kbd>Ctrl</kbd> / <kbd>⌘</kbd> + <kbd>K</kbd>。
