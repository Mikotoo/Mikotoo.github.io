# 内容结构整理与维护记录

本次整理保留 Astro + Starlight 和已有页面网址，调整内容来源、配置位置与构建职责。

## 编辑入口迁移

| 旧位置 | 新位置 | 说明 |
| --- | --- | --- |
| `src/content/docs/index.md` | `content/pages/index.md` | 首页所有文案、按钮与卡片 |
| `src/content/docs/intro.md` | `content/pages/intro.md` | 更新写作说明，移除普通 Markdown 中无效的组件导入 |
| `src/content/docs/about.md` / `resources.md` | `content/pages/` 同名文件 | 手写页面 |
| `scripts/gen-index.mjs` 内嵌的大段页面文案 | `content/pages/projects.md`、`code.md`、`index-all.md` | 自动区域用 `<!-- AUTO:... -->` 标记，文字可直接编辑 |
| `_posts/*.md` | `content/notes/*.md` | 使用已转换的笔记正文，停止每次构建重复迁移 |
| `src/site.config.mjs` | `config/site.mjs` | 移除未使用的首页副标题副本 |
| `astro.config.mjs` 内的导航骨架 | `config/navigation.mjs` | Astro 配置只负责集成 |
| `src/data/projects.json` | `config/projects.json` | 展示标题以手写配置为准 |
| `sync-docs.mjs` 内 LABELS | `config/project-labels.json` | 中文步骤标题独立维护 |
| `src/content/docs/_projects/` | `snapshots/projects/` | 项目快照，共 56 页 |
| `src/data/projects.generated.json` | `snapshots/manifest.json` | 稳定来源哈希代替运行时间戳和环境相关标志 |
| `src/sidebar.generated.mjs` | `.generated/sidebar.mjs` | 自动导航，不入库 |
| `public/media/`（以前自动复制） | `public/media/`（现在正式源文件） | 53 个已有媒体文件随站点提交 |

旧笔记与新文件的对应关系：

| 旧文件 | 新文件（位于 content/notes/） |
| --- | --- |
| `2022-04-06-genome.md` | `genome.md` |
| `2022-04-11-download.md` | `download.md` |
| `2022-04-12-sirna.md` | `sirna.md` |
| `2022-04-14-repeatmask.md` | `repeatmask.md` |
| `2022-04-19-RNAseq.md` | `rnaseq.md` |
| `2022-05-13-bsseq.md` | `bsseq.md` |
| `2024-4-7-steroSeq.md` | `soybean-snrna.md` |

旧源码保留在 Git 历史中。原来的 `downloads/` 和 `code/` 保留，不在本次整理中清理数据或重写历史。新配图使用 `public/media/`。

## 构建职责

```text
项目仓库 README
    │ pnpm sync（手动执行、所有源可读才更新）
    ▼
snapshots/（审阅并提交）
    │
    ├── content/（页面、笔记）
    ├── config/（站点、项目、标签、代码清单）
    │       pnpm gen
    ▼
src/content/docs/ + .generated/sidebar.mjs（自动装配，不编辑）
    │       Astro / Starlight
    ▼
dist/ + Pagefind 搜索索引
    │       体积检查与严格站内文件链接检查
    ▼
GitHub Pages
```

- 正常构建不读取外部项目目录，也不改变内容源或已提交快照。
- 同步只处理登记的项目，不再通过同级目录的 `content.config.ts` 自动发现。
- 项目源缺失时报错，不默默发布空文档；普通构建只要求对应快照存在。
- 快照先在 `.generated/sync-staging/` 全部生成，再替换发布目录。此过程防止读取或转换失败破坏旧快照，但不是跨文件系统事务；提交历史仍是恢复依据。
- 删除配置中的项目后，构建不再装配它；下次成功同步移除快照中的旧项目。
- `config/code-files.json` 显式列出本站代码页的文件，避免未跟踪的本地分析文件影响构建。
- 修正了 `00-prepare` 的数字排序，使第 0 步位于第 1 步之前。

## 页面与资源约定

`content/pages/` 的文件名映射到同名 URL；`index.md` 映射 `/`。`content/notes/` 文件名映射 `/notes/<文件名>/`。项目继续使用 `/_projects/<key>/...`，目录整理不改变这些 URL。

普通页面使用 `.md`，不直接导入 JSX/MDX 组件。需要交互组件时，再显式引入 MDX 支持。

页面自动标记必须各保留一次：

- `projects.md`：`<!-- AUTO:PROJECT_CARDS -->`
- `index-all.md`：`<!-- AUTO:CONTENT_INDEX -->`
- `code.md`：`<!-- AUTO:CODE_FILES -->`、`<!-- AUTO:PROJECT_REPOS -->`

## 验收记录

- 本地完整构建成功，71 个 HTML 页面、约 29.1 MB。
- 所有原有 HTML 页面路径与重构前逐一对比一致。
- 无同级项目目录的隔离副本构建通过；依赖复用本机 node_modules，内容目录、配置和快照均为独立复制。
- 隔离副本连续两次内容装配逐文件哈希一致，导航文件一致。
- 构建前后 content/config/snapshots/public 的逐文件哈希不变。
- 站内 href/src 文件目标检查通过，包含相对路径；默认严格失败。
- `pnpm test` 在临时目录验证同步可重复、README 改动改变来源哈希、源缺失不修改旧快照，以及真实死链能阻止检查通过。
- 本地 dev 服务实测：首页源文件保存后自动装配，HTTP 响应包含新内容；测试标记已恢复。

链接检查不覆盖外部 URL 可用性、页内锚点 ID、JavaScript 动态加载请求。Starlight 内置 404 canonical 元数据单独豁免；实际指向 `/404/` 的导航链接仍失败。

当前只是本地整理与验证，尚未提交或推送本次改动。
