# 生物信息学基础教程维护

## 当前课程

- 课程总览：`/tutorials/`
- 第 1 节：`/tutorials/01-linux-genome/`
- 第 2 节：`/tutorials/02-python-r-visualization/`
- 第 3 节：`/tutorials/03-rnaseq-alignment/`

正文完整保留原稿的教学内容、代码块、综合练习与思考题。导入时去掉原稿的顶层标题（网页标题由 frontmatter 显示），规范正文标题层级，添加课程目录、相邻章节与下载入口。代码块中的注释不作为 Markdown 标题处理。长讲义的页内目录只显示主要章节。

### 第 3 节的章节顺序

第 3 节按「**先做 → 再看懂 → 最后懂原理**」组织，这是刻意的教学设计，调整内容时请保持：

| 章 | 内容 | 作用 |
| --- | --- | --- |
| 三 | RNA-seq 数据是怎么产生的 | 背景：FASTQ、质量值，动手前必须知道 |
| 四 | 动手实操：从 FASTQ 到表达矩阵 | 先把流程跑通并保存输出，**不解释字段** |
| 五 | 读懂比对结果：SAM/BAM 文件 | 对照第四章自己生成的输出逐字段解释；5.11 节回头核对全部输出 |
| 六 | 比对算法原理 | 解释第五章见到的现象（`N`、`MAPQ 0`、`NH:i:2`） |

第四、五、六章之间存在互相引用（如「见 6.9 节」「第 4.10 节的实操」）。重排章节时必须同步更新这些编号引用；练习章第 7 章的题目顺序也按这条线排列。

## 编辑与更新

`config/tutorials.json` 登记原始目录、各节讲义的文件名、站点标题、简介、标签、顺序与练习数据。`sourceDir` 相对本站根目录解析，目前指向同级的 `../生物信息学教程`。

原始讲义仍在用户的教程目录中维护：

- `01_linux_genome_intro_practice.md`
- `02_visualization_simple sample/02_visualization_python_R.md`
- `02_visualization_simple sample/lotus_circos_data/` 下的 9 个数据文件
- `03_rnaseq_alignment/03_rnaseq_alignment_practice.md`
- `03_rnaseq_alignment/rnaseq_demo_data/`：`demo_genome.fa`、`demo_annotation.gtf`、`README.md`

修改原稿后执行：

```bash
pnpm import:tutorials
pnpm build
```

导入产物：

- `snapshots/tutorials/`：转换后的完整网页正文。
- `snapshots/tutorials.manifest.json`：按节记录原稿、正文、练习数据与 ZIP 的 SHA-256，并标记哪些文件是生成的。
- `public/tutorials/bioinformatics/`：各节原始 Markdown 讲义、练习数据目录与 ZIP。

以上产物与配置和代码一并提交。普通 `pnpm gen`、`pnpm build` 和 CI 只读取本站已保存的产物，不需要外部教程目录；不会在构建时自动更新原稿。不要直接修改快照或 `public/tutorials/`，下一次导入会覆盖它们。

新增章节时，将原始讲义登记到配置，再显式导入。课程卡片、侧栏、章节之间的前后链接由顺序生成；总览页的导语和学习路径可在 `content/pages/tutorials.md` 调整。

## 练习数据

配置中的 `practice` 字段声明某节的练习数据：

| 字段 | 作用 |
| --- | --- |
| `dir` | 原始数据目录，必须**扁平**（不含子目录），单文件上限 5 MiB |
| `assetName` | 发布目录名：`public/tutorials/bioinformatics/<assetName>/` |
| `archive`、`archiveRoot` | 生成的 ZIP 文件名与解压后的根目录名 |
| `label`、`summary` | 页面顶部下载提示的文案 |
| `generateReads`、`reads` | 由原稿的 FASTA + GTF **现场生成**教学 reads（不提交 FASTQ） |

第 2 节的 `lotus_circos_data/` 为原始数据直接发布；第 3 节的 `rnaseq_demo_data/` 中，`demo_reads_R1/R2.fastq.gz` 由导入器根据 `demo_genome.fa` 与 `demo_annotation.gtf` 用固定随机种子生成，因此**不提交 FASTQ**，重复导入字节一致。生成时刻意放入三类读长：跨外显子接头的、来自重复区的、基因组中不存在的，分别用于观察 CIGAR 的 `N`、`NH:i:2`/`MAPQ 0` 和未比对记录。

## 下载资料

- 三节原始 Markdown 讲义均可下载。
- `lotus-circos-data.zip`：解压得到 `lotus_circos_data/`，9 个原始文件，文件名不变。
- `rnaseq-demo-data.zip`：解压得到 `rnaseq_demo_data/`，含模拟参考序列、注释、280 对教学读长与说明。
- ZIP 使用固定日期与 STORE 格式，重复导入字节一致。
- 第 2 节基础绘图使用第 1 节大豆 FASTA/GFF3 生成的统计表；ZIP 中的 Lotus 数据用于 circlize 进阶练习，两者不混用。第 3 节演示数据是**教学模拟数据**，不可当作真实生物学结果。
- 已填写的 Word 作业、两本教材 PDF 和含教材的原始大 ZIP 未纳入网页公开附件。正文中原有的练习题完整保留。

## 校验

`.generated/` 下的临时脚本（不入库）可复现以下检查：三节正文与代码块逐字节保留、原稿下载与源文件一致、练习数据与清单哈希一致、重复导入字节一致，以及演示数据的教学属性（每个基因都有跨接头 read、重复区 read 恰好两个位点、噪声 read 在基因组中不存在）。

维护说明仅保存在仓库 docs/，不进入教程的公开正文。
