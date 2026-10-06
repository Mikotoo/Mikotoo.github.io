# Mikotoo.github.io

个人文档与代码站，部署在 GitHub Pages：<https://mikotoo.github.io>

不是在介绍某一个项目，而是一个**持续更新的集合**：每个项目一份文档区，另加手写笔记。技术栈是 **Astro + Starlight**（静态站点，构建于 GitHub Actions，本地无需 Ruby）。

## 内容来源

| 区域 | 来源 | 维护方式 |
| --- | --- | --- |
| 项目文档 | 各个项目仓库的 `README.md` | 构建时自动抓取，**不在本站维护** |
| 笔记与复现 | 本仓库 `_posts/*.md` | 沿用旧站写法，构建时自动转换 |
| 项目清单与卡片文案 | `src/data/projects.json` | 新增项目时加一条 |
| 站点外观、导航、首页 | `astro.config.mjs`、`src/site.config.mjs`、`src/styles/custom.css` | 手工维护 |

两条原则：

1. **流程说明只在代码仓库里写一遍**，站点构建时抓取，不会出现「网站一份、仓库一份」。
2. **写笔记只新建一个 Markdown 文件**，不需要改模板或配置。

## 本地开发

需要 Node.js ≥ 22.12（`node -v` 确认）与 pnpm。

```bash
pnpm install          # 首次
pnpm gen              # 抓取 README + 迁移笔记 + 生成索引与侧边栏（dev/build 会自动执行）
pnpm dev              # http://localhost:4321
pnpm build            # 产出 dist/，并做体积检查
pnpm preview          # 预览构建结果
```

> 受限沙箱中构建时需允许 esbuild 启动子进程；可用 `ASTRO_TELEMETRY_DISABLED=1` 关闭遥测。

## 收录一个新项目

项目文档**不在站点里手写**，构建时从项目仓库抓取。两种方式：

| 方式 | 做法 | 适用场景 |
| --- | --- | --- |
| 登记（推荐） | 在 `src/data/projects.json` 加一条，用 `dir`/`repo` 指向本地仓库 | 自己的项目，需定制标题、排序、关键词 |
| 自动发现 | 在仓库根目录放一个 `content.config.ts`（内容随意，只是标记） | 外部工具或临时克隆的仓库 |

两种方式都要求：仓库根目录有 `README.md`（每个步骤目录也建议一份），并且仓库位于本仓库的**同级目录**或 `vendor/` 下。

站内路径：`<仓库>/scripts/01_assembly/02_purge_haplotigs/README.md` → `/_projects/<key>/scripts/01-assembly/02-purge-haplotigs/`。

## 目录结构

```text
astro.config.mjs            站点配置与侧边栏骨架
src/site.config.mjs         站名、副标题、描述
src/data/projects.json      项目清单（展示顺序、卡片文案）
src/content/docs/           手写页面（index / intro / about / resources）
src/content/docs/_projects/ 项目文档，由仓库 README 生成（已 gitignore）
src/content/docs/notes/     由 _posts 生成（已 gitignore）
src/sidebar.generated.mjs   侧边栏（自动生成）
src/styles/custom.css       自定义样式
scripts/sync-docs.mjs       发现项目仓库并聚合 README
scripts/migrate-posts.mjs   迁移 _posts 并本地化配图
scripts/gen-index.mjs       生成项目总览 / 内容索引 / 代码页 / 侧边栏
scripts/check-size.mjs      构建体积检查
_posts/                     笔记源文件（手写）
downloads/                  图片与可下载资料
code/                       随文脚本与图件
```

## 写一篇新笔记

在 `_posts/` 新建 `YYYY-MM-DD-标题.md`：

```markdown
---
title: 文章标题
date: 2026-01-01
category: 笔记
---

正文……
```

图片放到 `downloads/image/<目录>/` 下，正文用本仓库的 GitHub 原始地址引用：

```markdown
![说明][1]

[1]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/xxx/yyy.png
```

构建时会自动复制到 `public/media/` 并改写为站点绝对路径，因此线上不依赖 GitHub raw 的可用性。新增文章后侧边栏由 `pnpm gen` 自动更新，无需手改配置。

## 部署

`.github/workflows/deploy.yml`：推送到 `main` 后由 GitHub Actions 构建并发布到 Pages（Settings → Pages → Source 需选择 **GitHub Actions**）。仓库名是 `Mikotoo.github.io`，属于用户站点，站点根路径为 `/`。

`.github/workflows/preview.yml`：`preview` 分支只验证构建、产物作 artifact 上传，不发布。

CI 上没有同级项目仓库，`sync-docs.mjs` 会跳过它们并给出提示，因此线上只发布笔记与手写页面；本地的完整项目文档在构建时自然包含。

## 已知的体积问题

仓库历史中提交过单细胞分析缓存（`code/single_cell/02_QC/cache/*.h5ad`，约 217 MB）与 `web_summary.html`。它们不会进入站点产物，但会让 `git clone` 变慢。要清理可把这些数据迁到 Zenodo 等外部存储，再用 `git filter-repo` 重写历史。

## 授权

本站文章与自研脚本采用 MIT 许可；引用的第三方资源保留原作者条款。
