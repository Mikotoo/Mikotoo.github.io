# Mikotoo.github.io

Mikotoo 的个人研究与技术展示站，基于 Astro + Starlight，发布于 <https://mikotoo.github.io>。

访客入口展示研究方向、项目、笔记与代码；本文及 `docs/` 面向维护者，不作为访客页面发布。内容与布局设计见 `docs/presentation-design.md`。

## 日常编辑入口

| 要修改的内容 | 文件或目录 |
| --- | --- |
| 首页大标题、副标题、按钮、卡片与正文 | `content/pages/index.md` |
| 项目总览介绍、代码页介绍、关于与资源说明 | `content/pages/` 对应页面 |
| 笔记与复现 | `content/notes/*.md` |
| 基础教程总览与课程卡片 | `content/pages/tutorials.md`、`config/tutorials.json` |
| 基础教程正文 | 在原始 `生物信息学教程/` 目录编辑，再运行 `pnpm import:tutorials`；详见 `docs/tutorials.md` |
| 图片与附件 | `public/media/` |
| 顶栏站名、作者、站点与仓库地址 | `config/site.mjs` |
| 文章侧栏固定分组 | `config/navigation.mjs` |
| 顶部主导航 | `src/components/Header.astro` |
| 项目名称、简介、来源路径、排序、首页精选 | `config/projects.json` |
| 笔记卡片标题、摘要、分类、首页精选 | `config/readings.json` |
| 项目步骤中文标题 | `config/project-labels.json` |
| 代码页列出的随文脚本 | `config/code-files.json` |
| 颜色、行距、卡片样式 | `src/styles/custom.css` |

`content/pages/` 中的 `<!-- AUTO:... -->` 是自动区域占位符。可自由修改周围正文；每个标记保留一次，构建会填入项目或笔记卡片、代码清单及内容索引。首页两个精选区的文案分别来自 `config/projects.json`、`config/readings.json`，由 `featured: true` 选入。

## 本地运行

使用 Node.js 22.12+（建议 Node 24）与项目指定版本的 pnpm。

```bash
pnpm install --frozen-lockfile
pnpm dev       # 本地预览，监听 content/、config/、snapshots/ 的保存
pnpm build     # 装配、构建、体积与严格站内文件链接检查
pnpm preview   # 查看 dist/ 构建结果
pnpm test      # 隔离目录中的同步和链接检查回归测试
```

`pnpm gen` 只装配内容，不更新外部文档、不读取同级项目目录，也不修改内容源或快照。

## 写笔记

在 `content/notes/` 新建 `my-note.md`：

```markdown
---
title: 我的笔记
date: 2026-01-01
---

正文……

![示例](/media/my-note/example.png)
```

图片放到 `public/media/my-note/example.png`。笔记网址为 `/notes/my-note/`；文件名应保持稳定。普通 Markdown 页面不要使用组件 import；当前站点未配置 MDX。

## 更新项目快照

项目正文仍在各项目仓库的 README 中维护；本站发布已提交的快照。

```bash
pnpm sync      # 显式读取 config/projects.json 中登记的源目录
pnpm build     # 验证快照能构建且链接正常
```

审阅并提交 `snapshots/` 的修改。只改项目 README 不会自动更新线上网站，普通 `pnpm build` 也不会偷偷更新快照。

- `key` 决定 `/_projects/<key>/`，发布后不要随意修改。
- `dir` 相对本站父目录解析，也支持绝对路径。例如同级项目填 `Lotus_genome`；站点内 vendor 目录可填 `Mikotoo.github.io/vendor/项目名`。
- `href` 指明 GitHub 仓库，`branch` 可指定链接分支，默认 `main`。
- 根目录必须有 README；同步前检查所有项目源。源缺失或读取失败会报错，保留现有快照。
- 同步先写 `.generated/sync-staging/`，所有页面生成成功后才发布快照。
- `snapshots/manifest.json` 记录项目页数及 README 路径与原始字节的 SHA256；不记录每次运行变化的时间戳。该哈希表示 README 来源，不是整个 Git 提交或网站产物的哈希。
- 添加项目需要登记配置并运行同步；删除登记后，普通构建即不再发布该项目，下次同步清除过期快照。

## 目录职责

```text
content/pages/           手写页面与带自动占位符的页面
content/notes/           手写笔记（已从旧 _posts 迁移）
config/                  手工配置
snapshots/projects/      项目文档快照，必须提交，不手改
snapshots/manifest.json  快照来源摘要，必须提交
public/                  正式静态资源，直接复制到网站
src/styles/              自定义样式
src/content.config.ts    Starlight 内容加载配置
src/content/docs/        自动装配结果，已忽略，不编辑
.generated/              导航与同步暂存等，已忽略，不编辑
scripts/                 同步、装配、开发监听与检查工具
code/                    随文分析代码（只发布登记的链接）
downloads/               历史下载资源，保留既有 GitHub 链接
docs/                    维护记录与预览资料
```

新配图统一放 `public/media/`。现有笔记引用的 53 个资源已固定在此目录，构建不再从 `downloads/` 复制。历史 `downloads/` 保留以避免破坏外部引用；大型分析数据和 Git 历史清理另行处理。

## 发布与检查范围

- `preview` 分支：Actions 构建与测试，上传预览产物。
- `main` 分支：Actions 构建与测试，成功后发布 GitHub Pages。
- CI 只需本仓库内容、配置、静态资源与快照，不需要 Lotus_genome/Rscript 同级目录。
- 站内链接检查默认严格，缺文件返回非零退出码，阻止部署。检查 HTML 的 href/src 文件目标，包括相对路径；不检查外部网站可用性、锚点 ID 或 JavaScript 动态请求。
- Starlight 的 `404.html` canonical `/404/` 元数据单独豁免，实际指向 `/404/` 的导航链接仍会报错。

详细迁移映射与验收记录见 `docs/maintenance.md`。

## 授权说明

本次整理不新增或变更授权。引用的第三方资料保留其原作者条款；本站内容的统一许可需由作者明确指定。
