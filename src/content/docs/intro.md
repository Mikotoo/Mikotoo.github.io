---
title: 站点说明
description: 内容怎么组织、新内容怎么加、构建时发生了什么
---

import { CardGrid, LinkCard } from '@astrojs/starlight/components';

这是一个**持续更新的个人文档 + 代码站**。三类内容，维护方式完全不同：

<CardGrid>
  <LinkCard
    title="项目文档"
    href="/projects/"
    description="每个项目的分析流程与步骤说明。内容来自代码仓库的 README，构建时自动同步"
  />
  <LinkCard
    title="笔记与复现"
    href="/index-all/"
    description="手写的教程、流程记录、论文复现与踩坑。新增就是新建一个 Markdown 文件"
  />
  <LinkCard
    title="代码与资源"
    href="/code/"
    description="随文脚本、分析图件与可下载资料"
  />
</CardGrid>

## 谁说了算

| 内容 | 放在哪 | 谁说了算 |
| --- | --- | --- |
| 项目流程与参数 | 各项目自己的代码仓库 | **代码仓库**：根目录与各步骤目录的 `README.md` 是唯一事实来源，站点构建时抓取 |
| 笔记、教程、论文复现 | 本仓库 `_posts/` | **这里的 Markdown 文件**：构建时自动转成站点页面 |
| 站点外观、导航、首页 | 本仓库 `astro.config.mjs`、`src/` | 本仓库 |

这样分工的好处：改流程说明时只改代码仓库，站点不会和仓库说法不一致；写笔记时不碰模板和配置。

## 加一个新项目

项目文档**不在站点里手写**，而是构建时从项目仓库抓取。两种收录方式：

| 方式 | 做法 | 适用场景 |
| --- | --- | --- |
| 登记（推荐） | 在 `src/data/projects.json` 里加一条，用 `dir`/`repo` 指向本地仓库 | 自己的项目，需要定制标题、排序、关键词 |
| 自动发现 | 在仓库根目录放一个 `content.config.ts`（内容随意，只是标记） | 外部工具或临时克隆的仓库 |

无论哪种方式，都要求：

1. **仓库根目录有 `README.md`**（每个步骤/子目录也建议一份，讲清做法、输入、参数）；
2. **仓库放在本站同级目录**或 `vendor/` 下（后者适合作为 submodule 收进站点仓库）；
3. 重新构建：`pnpm gen && pnpm dev`。

站内路径规则：`<仓库>/scripts/01_assembly/02_purge_haplotigs/README.md` → `/_projects/<key>/scripts/01-assembly/02-purge-haplotigs/`。侧边栏的层级与「阶段 → 步骤」顺序都按目录自动生成。

> 想让某个目录显示成中文标签，在 `scripts/sync-docs.mjs` 的 `LABELS` 里加一条映射即可；不加也能跑，只是标签会是目录名。

## 加一篇笔记

在 `_posts/` 新建 `YYYY-MM-DD-标题.md`：

```markdown
---
title: 文章标题
date: 2026-01-01
category: 笔记
---

正文……
```

侧边栏与内容索引会自动带上它，不需要改配置。

图片放在 `downloads/image/<目录>/` 下，正文里用本仓库的 GitHub 原始地址引用：

```markdown
![说明][1]

[1]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/xxx/yyy.png
```

构建时会自动复制到 `public/media/` 并改写为站点绝对路径，**线上不依赖 GitHub raw 是否可访问**。

## 构建时发生了什么

```bash
pnpm install          # 首次
pnpm gen              # 三个生成脚本（dev / build 会自动执行，无需手跑）
pnpm dev              # http://localhost:4321
pnpm build            # 产出 dist/ 并做体积检查
```

| 脚本 | 作用 |
| --- | --- |
| `scripts/sync-docs.mjs` | 发现带 `content.config.ts` 的仓库，抓取全部 `README.md` 生成项目文档区；相对链接改写为 GitHub 链接，图片指向 raw |
| `scripts/migrate-posts.mjs` | 把 `_posts/*.md` 转成站点页面，并把正文里本仓库的图片地址本地化到 `public/media/` |
| `scripts/gen-index.mjs` | 生成项目总览页、内容索引、代码页，以及**完整侧边栏**（`src/sidebar.generated.mjs`） |

生成结果都写在 `src/content/docs/` 与 `src/sidebar.generated.mjs`，已在 `.gitignore` 中忽略——**每次构建重新生成，不需要提交**。

## 一些说明

- **双语 README**：项目 README 若是中文在前、英文在后，站点默认展开中文，英文收进可折叠区。
- **大文件不进站点**：仓库历史里提交过单细胞分析缓存（`.h5ad`，合计约 220 MB），构建后由 `scripts/prune-dist.mjs` 从产物中剔除，仍可在 GitHub 上访问。
- **全文搜索**：Starlight 自带，快捷键 <kbd>Ctrl</kbd> / <kbd>⌘</kbd> + <kbd>K</kbd>。
- **换站名/副标题**：改 `src/site.config.mjs`。首页文案在 `src/content/docs/index.md`。
