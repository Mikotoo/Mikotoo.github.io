# 生物信息学基础教程维护

## 当前课程

- 课程总览：`/tutorials/`
- 第 1 节：`/tutorials/01-linux-genome/`
- 第 2 节：`/tutorials/02-python-r-visualization/`

正文完整保留原稿的教学内容、代码块、综合练习与思考题。导入时去掉原稿的顶层标题（网页标题由 frontmatter 显示），规范正文标题层级，添加课程目录、相邻章节与下载入口。代码块中的注释不作为 Markdown 标题处理。长讲义的页内目录只显示主要章节。

## 编辑与更新

`config/tutorials.json` 登记原始目录、两节讲义的文件名、站点标题、简介、标签和顺序。`sourceDir` 相对本站根目录解析，目前指向同级的 `../生物信息学教程`。

原始讲义仍在用户的教程目录中维护：

- `01_linux_genome_intro_practice.md`
- `02_visualization_simple sample/02_visualization_python_R.md`
- `02_visualization_simple sample/lotus_circos_data/` 下的 9 个数据文件

修改原稿后执行：

```bash
pnpm import:tutorials
pnpm build
```

导入产物：

- `snapshots/tutorials/`：转换后的完整网页正文。
- `snapshots/tutorials.manifest.json`：原稿、正文、数据和 ZIP 的 SHA-256。
- `public/tutorials/bioinformatics/`：原始 Markdown 讲义、9 个练习数据文件、圈图 ZIP。

以上产物与配置和代码一并提交。普通 `pnpm gen`、`pnpm build` 和 CI 只读取本站已保存的产物，不需要外部教程目录；不会在构建时自动更新原稿。不要直接修改快照，下一次导入会覆盖它们。

新增章节时，将原始讲义登记到配置，再显式导入。课程卡片、侧栏、章节之间的前后链接由顺序生成；总览页的导语和学习路径可在 `content/pages/tutorials.md` 调整。

## 下载资料

- 两节原始 Markdown 讲义均可下载。
- `lotus-circos-data.zip` 解压后得到 `lotus_circos_data/` 及其中 9 个原始文件，文件名不变。
- ZIP 使用固定日期和 STORE 格式，重复导入字节一致，约 449 KiB。
- 基础绘图使用第一节大豆 FASTA/GFF3 生成的统计表；ZIP 中的 Lotus 数据用于 circlize 进阶练习，两者不混用。
- 已填写的 Word 作业、两本教材 PDF 和含教材的原始大 ZIP 未纳入网页公开附件。正文中原有的练习题完整保留。

维护说明仅保存在仓库 docs/，不进入教程的公开正文。
