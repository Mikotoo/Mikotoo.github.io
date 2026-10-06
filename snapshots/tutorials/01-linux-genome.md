---
title: "第 1 节：Linux 与参考基因组"
description: "从 Linux 基础命令出发，认识 FASTA、GFF3 与基因结构，并完成大豆参考基因组的查看、提取和统计练习。"
sidebar:
  label: "第 1 节：Linux 与参考基因组"
  order: 1
prev: {"label":"课程目录","link":"/tutorials/"}
next: {"label":"第 2 节：Python 与 R 数据可视化","link":"/tutorials/02-python-r-visualization/"}
tableOfContents:
  minHeadingLevel: 2
  maxHeadingLevel: 2
---

[← 课程目录](/tutorials/) · 第 1 节 / 共 2 节 · [下载本节原始讲义（Markdown）](/tutorials/bioinformatics/01-linux-genome.md)

## 一、章节导入

生物信息学分析通常从“数据文件”开始。对于基因组学研究来说，最基础、最常用的两类文件是：

1. **参考基因组序列文件**，通常是 FASTA 格式；
2. **基因结构注释文件**，通常是 GFF3 或 GTF 格式。

例如我们研究大豆基因时，首先需要知道：

- 大豆参考基因组有哪些染色体？
- 每条染色体的序列是什么？
- 基因位于哪条染色体上？
- 基因的起始位置和终止位置在哪里？
- 一个基因由哪些外显子、CDS、UTR 等结构组成？
- 大豆基因组中一共有多少个基因？

这些问题都可以通过 Linux 基础命令对 FASTA 和 GFF 文件进行查看、提取和统计来完成。

本节课将以大豆参考基因组为例，学习如何从 Phytozome 下载基因组 FASTA 文件和 GFF 注释文件，并使用 Linux 命令进行初步分析。

> 课堂提醒：Phytozome、SoyBase、Ensembl Plants、NCBI 等数据库中，同一物种可能存在多个参考基因组版本。下载数据时，应确保 FASTA 文件和 GFF 文件来自同一个物种、同一个品种、同一个版本。

---

## 二、本章学习目标

完成本章学习后，学生应能够：

1. 理解 Linux 命令行在生物信息学中的基本作用；
2. 掌握常用 Linux 基础命令，如 `pwd`、`ls`、`cd`、`mkdir`、`cp`、`mv`、`rm`、`cat`、`less`、`head`、`tail`、`grep`、`wc`、`awk`；
3. 了解参考基因组数据库中常见文件类型；
4. 理解 FASTA 文件和 GFF3 文件的基本结构；
5. 能够查看大豆参考基因组文件；
6. 能够从 GFF 文件中提取特定基因信息；
7. 能够使用 `awk` 统计基因、mRNA、exon、CDS 等数量；
8. 初步理解“基因由哪些组成部分构成”。

---

## 三、Linux 基础命令介绍

### 3.1 目录和文件操作命令

| 命令 | 作用 | 示例 |
|---|---|---|
| `pwd` | 显示当前所在目录 | `pwd` |
| `ls` | 查看当前目录下的文件 | `ls` |
| `ls -lh` | 以易读方式显示文件大小 | `ls -lh` |
| `cd` | 进入某个目录 | `cd genome` |
| `cd ..` | 返回上一级目录 | `cd ..` |
| `mkdir` | 新建目录 | `mkdir soybean_genome` |
| `mkdir -p` | 递归创建目录，如果上级目录不存在也会自动创建 | `mkdir -p ~/bioinfo_course/01_linux_genome` |
| `cp` | 复制文件 | `cp file1 file2` |
| `mv` | 移动或重命名文件 | `mv old.txt new.txt` |
| `rm` | 删除文件 | `rm test.txt` |
| `rm -r` | 删除目录及其内容 | `rm -r test_dir` |

示例：

```bash
pwd
mkdir soybean_genome
cd soybean_genome
pwd
```

---

### 3.2 查看文件内容的命令

| 命令 | 作用 | 示例 |
|---|---|---|
| `cat` | 一次性显示整个文件内容 | `cat test.txt` |
| `less` | 分页查看大文件 | `less genome.fa` |
| `head` | 查看文件前几行 | `head file.txt` |
| `tail` | 查看文件后几行 | `tail file.txt` |
| `head -n 20` | 查看前 20 行 | `head -n 20 file.txt` |
| `tail -n 20` | 查看后 20 行 | `tail -n 20 file.txt` |

对于基因组文件，通常不建议直接用 `cat` 查看，因为文件很大，容易刷屏。推荐使用：

```bash
less 文件名
head 文件名
tail 文件名
```

如果文件是压缩文件，例如 `.gz` 格式，可以使用：

```bash
zcat 文件名.gz | head
zcat 文件名.gz | less
```

---

### 3.3 搜索和统计命令

| 命令 | 作用 | 示例 |
|---|---|---|
| `grep` | 按关键词搜索 | `grep "gene" annotation.gff3` |
| `grep -v` | 排除某些行 | `grep -v "^#" annotation.gff3` |
| `wc -l` | 统计行数 | `wc -l file.txt` |
| `cut` | 按列提取 | `cut -f 1 file.txt` |
| `sort` | 排序 | `sort file.txt` |
| `uniq` | 去重复 | `uniq file.txt` |
| `awk` | 按列处理和统计 | `awk '$3=="gene"' annotation.gff3` |

`awk` 是生物信息学中非常重要的文本处理工具，尤其适合处理以制表符分隔的文件，例如 GFF、BED、VCF 等。

---

## 四、基因组数据库中常见文件类型

在 Phytozome、Ensembl Plants、NCBI、SoyBase 等数据库中，一个物种的基因组页面通常会包含多种文件。

### 4.1 Genome FASTA 文件

常见后缀：

```text
.fa
.fasta
.fa.gz
.fasta.gz
```

作用：

FASTA 文件保存的是参考基因组的 DNA 序列。例如大豆有 20 条染色体，FASTA 文件中通常包含 20 条主要染色体序列，以及一些未定位的 scaffold 或 contig。

FASTA 文件格式示例：

```text
>Gm01
ATGCTAGCTAGCTACGATCGATCG...
>Gm02
ATCGATCGATCGTACGTAGCTAGC...
```

特点：

- 以 `>` 开头的行是序列名称；
- 后面的多行是 DNA 序列；
- DNA 序列通常由 A、T、C、G、N 组成；
- `N` 表示未知碱基。

---

### 4.2 GFF3 注释文件

常见后缀：

```text
.gff
.gff3
.gff3.gz
```

作用：

GFF3 文件记录基因在基因组上的位置和结构信息，例如：

- gene；
- mRNA；
- exon；
- CDS；
- five_prime_UTR；
- three_prime_UTR。

GFF3 文件一般有 9 列：

| 列号 | 名称 | 含义 |
|---|---|---|
| 第 1 列 | seqid | 染色体或序列名称 |
| 第 2 列 | source | 注释来源 |
| 第 3 列 | type | 注释类型，如 gene、mRNA、exon、CDS |
| 第 4 列 | start | 起始位置 |
| 第 5 列 | end | 终止位置 |
| 第 6 列 | score | 得分，若无则为 `.` |
| 第 7 列 | strand | 正负链，`+` 或 `-` |
| 第 8 列 | phase | CDS 阅读框信息 |
| 第 9 列 | attributes | 注释属性，如 ID、Parent、Name |

示例：

```text
Gm01  phytozome  gene  1000  5000  .  +  .  ID=Glyma.01G000100
Gm01  phytozome  mRNA  1000  5000  .  +  .  ID=Glyma.01G000100.1;Parent=Glyma.01G000100
Gm01  phytozome  exon  1000  1200  .  +  .  Parent=Glyma.01G000100.1
Gm01  phytozome  CDS   1050  1200  .  +  0  Parent=Glyma.01G000100.1
```

---

### 4.3 CDS FASTA 文件

常见后缀：

```text
.cds.fa
.cds.fasta
```

作用：

CDS 文件保存基因的编码序列，也就是可以翻译成蛋白质的 DNA 序列。

---

### 4.4 Protein FASTA 文件

常见后缀：

```text
.pep.fa
.protein.fa
.aa.fa
```

作用：

Protein FASTA 文件保存蛋白质序列。后续做 BLAST、结构域分析、同源基因搜索、系统发育树分析时经常使用。

---

### 4.5 Transcript FASTA 文件

常见后缀：

```text
.transcript.fa
.mrna.fa
```

作用：

Transcript FASTA 文件保存转录本序列，通常包括 UTR 和 CDS 区域。

---

### 4.6 README 或 metadata 文件

这些文件通常说明：

- 数据来源；
- 版本号；
- 文件命名规则；
- 物种信息；
- 下载日期；
- 使用说明。

课堂中要特别强调：

> **FASTA 文件和 GFF 文件必须来自同一个基因组版本。**

例如不能把 Wm82.a2 的 FASTA 和 Wm82.a4 的 GFF 混合使用，否则基因坐标可能对不上。

---

## 五、基因的基本组成

从基因组注释角度看，一个典型真核基因通常包括：

```text
gene
 └── mRNA / transcript
      ├── 5' UTR
      ├── exon
      ├── CDS
      ├── intron
      ├── exon
      ├── CDS
      └── 3' UTR
```

### 5.1 gene

表示一个基因在基因组上的完整区域。

### 5.2 mRNA / transcript

表示该基因产生的转录本。一个基因可能只有一个转录本，也可能有多个转录本。

### 5.3 exon

外显子，是成熟 mRNA 中保留下来的片段。

### 5.4 intron

内含子，是转录后会被剪接去除的片段。

注意：在很多 GFF 文件中，intron 不一定直接写出来，可以根据相邻 exon 的间隔推断。

### 5.5 CDS

CDS 是 coding sequence 的缩写，表示编码区，是最终可以翻译成蛋白质的序列部分。

### 5.6 UTR

UTR 是 untranslated region 的缩写，表示非翻译区，包括：

- 5' UTR；
- 3' UTR。

UTR 不翻译成蛋白质，但可能影响 mRNA 稳定性、翻译效率和调控。

---

## 六、课堂操作练习

### 练习 1：建立工作目录

#### 任务

在自己的 Linux 环境中建立一个本章练习目录。

#### 命令

```bash
mkdir -p ~/bioinfo_course/01_linux_genome
cd ~/bioinfo_course/01_linux_genome
pwd
```

#### 思考题

1. `mkdir -p` 中的 `-p` 有什么作用？
2. `pwd` 输出的是什么？
3. `~` 代表什么目录？

---

### 练习 2：从 Phytozome 下载大豆参考基因组文件

#### 任务

进入 Phytozome 网站，下载大豆参考基因组的 FASTA 文件和 GFF3 注释文件。

Phytozome 网站：

```text
https://phytozome-next.jgi.doe.gov/
```

#### 操作步骤

1. 打开 Phytozome 网站；
2. 在搜索框中搜索：

```text
Glycine max
```

或：

```text
soybean
```

3. 选择大豆 Williams 82 参考基因组版本；
4. 课堂建议统一选择同一个版本，例如：

```text
Glycine max Wm82.a4.v1
```

5. 下载以下两类文件：

```text
Genome FASTA 文件
Annotation GFF3 文件
```

下载后文件名可能类似：

```text
Gmax_508_v4.0.fa.gz
Gmax_508_Wm82.a4.v1.gene.gff3.gz
```

实际文件名以网站下载页面为准。

---

### 练习 3：查看下载后的文件

#### 任务

查看当前目录下有哪些文件，并观察文件大小。

#### 命令

```bash
ls
ls -lh
```

#### 问题

1. FASTA 文件的后缀是什么？
2. GFF 文件的后缀是什么？
3. 哪个文件更大？
4. 为什么基因组 FASTA 文件通常比较大？

---

### 练习 4：查看 FASTA 文件前几行

假设下载的基因组文件名为：

```text
Gmax_genome.fa.gz
```

如果你的文件名不同，请替换成自己的文件名。

#### 命令

```bash
zcat Gmax_genome.fa.gz | head
```

如果文件未压缩：

```bash
head Gmax_genome.fa
```

#### 观察内容

你应该能看到类似结构：

```text
>Gm01
NNNNNNNNNNNNNNNNNNNN
ATGCTAGCTAGCTAGCTAGC
```

#### 问题

1. FASTA 文件中，哪一类行表示序列名称？
2. DNA 序列由哪些字母组成？
3. `N` 代表什么？

---

### 练习 5：统计 FASTA 文件中有多少条序列

#### 命令

```bash
zgrep "^>" Gmax_genome.fa.gz | head
```

统计序列数量：

```bash
zgrep "^>" Gmax_genome.fa.gz | wc -l
```

#### 解释

在 FASTA 文件中，每一条序列都以 `>` 开头。因此统计 `>` 的数量，就可以知道 FASTA 文件中有多少条序列。

#### 思考题

1. 如果结果是 20，说明什么？
2. 如果结果大于 20，可能说明什么？
3. 染色体、scaffold、contig 有什么区别？

---

### 练习 6：提取 FASTA 文件中的序列名称

#### 命令

```bash
zgrep "^>" Gmax_genome.fa.gz | head -n 20
```

只保留序列名称的第一部分：

```bash
zgrep "^>" Gmax_genome.fa.gz | sed 's/>//' | awk '{print $1}' | head -n 20
```

#### 问题

1. 大豆主要染色体是否是 20 条？
2. 序列名称是 `Gm01`、`Chr01`，还是其他格式？
3. 为什么后续分析中要注意染色体命名格式？

---

### 练习 7：查看 GFF3 文件前几行

假设注释文件名为：

```text
Gmax_annotation.gff3.gz
```

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | head
```

跳过以 `#` 开头的注释行：

```bash
zcat Gmax_annotation.gff3.gz | grep -v "^#" | head
```

#### 问题

1. GFF3 文件有几列？
2. 每列之间用什么分隔？
3. 第 3 列表示什么？
4. 第 4 列和第 5 列表示什么？

---

### 练习 8：查看 GFF 文件中有哪些注释类型

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | grep -v "^#" | cut -f 3 | sort | uniq
```

统计每种类型的数量：

```bash
zcat Gmax_annotation.gff3.gz | grep -v "^#" | cut -f 3 | sort | uniq -c
```

#### 可能看到的类型

```text
gene
mRNA
exon
CDS
five_prime_UTR
three_prime_UTR
```

#### 问题

1. `gene` 和 `mRNA` 是同一个概念吗？
2. `exon` 和 `CDS` 有什么区别？
3. 为什么一个基因可能对应多个 mRNA？

---

### 练习 9：统计大豆基因数量

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | awk '$3=="gene"' | wc -l
```

更严谨地跳过注释行：

```bash
zcat Gmax_annotation.gff3.gz | awk '$0 !~ /^#/ && $3=="gene"' | wc -l
```

#### 解释

在 GFF3 文件中，第 3 列表示 feature 类型。当第 3 列等于 `gene` 时，该行代表一个基因。

#### 问题

1. 大豆基因组中一共有多少个 gene？
2. 这个数量是否等于 mRNA 数量？
3. 如果 mRNA 数量多于 gene 数量，说明什么？

---

### 练习 10：分别统计 gene、mRNA、exon、CDS 数量

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | awk '$0 !~ /^#/ {print $3}' | sort | uniq -c
```

只统计指定类型：

```bash
zcat Gmax_annotation.gff3.gz | awk '$3=="gene"' | wc -l
zcat Gmax_annotation.gff3.gz | awk '$3=="mRNA"' | wc -l
zcat Gmax_annotation.gff3.gz | awk '$3=="exon"' | wc -l
zcat Gmax_annotation.gff3.gz | awk '$3=="CDS"' | wc -l
```

#### 练习记录表

| 类型 | 数量 |
|---|---|
| gene |  |
| mRNA |  |
| exon |  |
| CDS |  |

#### 思考题

1. 为什么 exon 数量远大于 gene 数量？
2. 为什么 CDS 数量可能和 exon 数量不完全相同？
3. UTR 是否属于 CDS？

---

### 练习 11：提取前 10 个基因的信息

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | awk '$3=="gene"' | head
```

提取前 10 个基因的染色体、起始位置、终止位置和属性信息：

```bash
zcat Gmax_annotation.gff3.gz | awk '$3=="gene" {print $1,$4,$5,$7,$9}' | head -n 10
```

如果想让输出以制表符分隔：

```bash
zcat Gmax_annotation.gff3.gz | awk 'BEGIN{OFS="\t"} $3=="gene" {print $1,$4,$5,$7,$9}' | head -n 10
```

#### 问题

1. 第 1 列是什么？
2. 第 4、5 列是什么？
3. 第 7 列的 `+` 和 `-` 代表什么？
4. 第 9 列中是否包含基因 ID？

---

### 练习 12：自动提取第一个基因 ID

不同版本的 GFF 文件中，基因 ID 格式可能略有差异。可以先自动提取第一个 gene 的 ID。

#### 命令

```bash
GENE_ID=$(zcat Gmax_annotation.gff3.gz | awk '$3=="gene"{print $9; exit}' | sed -E 's/.*ID=([^;]+).*/\1/')

echo $GENE_ID
```

#### 解释

这条命令完成了三件事：

1. 找到第一个 `gene` 行；
2. 提取第 9 列属性信息；
3. 从 `ID=xxx` 中提取基因 ID。

---

### 练习 13：提取某一个特定基因的全部注释信息

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | grep "$GENE_ID"
```

如果你已经知道某个基因 ID，例如：

```text
Glyma.01G000100
```

可以使用：

```bash
zcat Gmax_annotation.gff3.gz | grep "Glyma.01G000100"
```

#### 观察内容

你可能会看到：

```text
gene
mRNA
exon
CDS
UTR
```

等多种注释行。

#### 问题

1. 这个基因位于哪条染色体？
2. 起始位置是多少？
3. 终止位置是多少？
4. 位于正链还是负链？
5. 这个基因有几个 exon？
6. 这个基因有几个 CDS？

---

### 练习 14：统计某个基因有多少个 exon 和 CDS

以 `$GENE_ID` 为例：

```bash
zcat Gmax_annotation.gff3.gz | grep "$GENE_ID" | awk '$3=="exon"' | wc -l
```

```bash
zcat Gmax_annotation.gff3.gz | grep "$GENE_ID" | awk '$3=="CDS"' | wc -l
```

#### 问题

1. exon 数量是多少？
2. CDS 数量是多少？
3. exon 数量和 CDS 数量是否相同？
4. 如果不同，可能是什么原因？

---

### 练习 15：计算每个基因的长度

基因长度可以用：

```text
end - start + 1
```

计算。

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | awk 'BEGIN{OFS="\t"} $3=="gene" {print $1,$4,$5,$5-$4+1,$9}' | head
```

输出列含义：

```text
染色体    起始位置    终止位置    基因长度    属性信息
```

#### 问题

1. 为什么要加 1？
2. 基因长度是否等于 CDS 长度？
3. 一个基因很长，是否一定编码很长的蛋白？

---

### 练习 16：统计每条染色体上的基因数量

#### 命令

```bash
zcat Gmax_annotation.gff3.gz | awk '$3=="gene"{print $1}' | sort | uniq -c
```

整理成更清楚的格式：

```bash
zcat Gmax_annotation.gff3.gz | awk 'BEGIN{OFS="\t"} $3=="gene"{count[$1]++} END{for(chr in count) print chr,count[chr]}' | sort -k1,1
```

#### 记录表

| 染色体 | 基因数量 |
|---|---|
| Gm01 |  |
| Gm02 |  |
| Gm03 |  |
| ... |  |

#### 思考题

1. 每条染色体上的基因数量是否相同？
2. 基因数量多的染色体一定更长吗？
3. 染色体长度和基因数量之间可能有什么关系？

---

## 七、综合练习

### 综合任务 1：完成大豆参考基因组文件初步检查

请完成以下任务，并将命令和结果记录下来。

#### 任务要求

1. 创建工作目录；
2. 下载大豆参考基因组 FASTA 和 GFF3 文件；
3. 查看两个文件的大小；
4. 查看 FASTA 文件前 10 行；
5. 统计 FASTA 文件中有多少条序列；
6. 查看 GFF3 文件前 10 行非注释内容；
7. 统计 GFF3 中 gene、mRNA、exon、CDS 的数量；
8. 提取一个基因的完整注释信息；
9. 统计每条染色体上的基因数量。

---

### 综合任务 2：整理结果表格

请学生整理如下表格：

| 项目 | 结果 |
|---|---|
| 使用的物种 | Glycine max |
| 使用的品种/参考基因组 | Williams 82 |
| 使用的版本 | 例如 Wm82.a4.v1 |
| FASTA 文件名 |  |
| GFF3 文件名 |  |
| FASTA 中序列数量 |  |
| gene 数量 |  |
| mRNA 数量 |  |
| exon 数量 |  |
| CDS 数量 |  |
| 示例基因 ID |  |
| 示例基因所在染色体 |  |
| 示例基因起始位置 |  |
| 示例基因终止位置 |  |
| 示例基因正负链 |  |
| 示例基因 exon 数量 |  |
| 示例基因 CDS 数量 |  |

---

## 八、课后思考题

### 8.1 概念题

1. 什么是参考基因组？
2. FASTA 文件主要保存什么信息？
3. GFF3 文件主要保存什么信息？
4. gene、mRNA、exon、CDS、UTR 分别是什么意思？
5. 为什么基因组 FASTA 和 GFF 文件必须来自同一个版本？

---

### 8.2 命令题

解释以下命令的含义：

```bash
zgrep "^>" genome.fa.gz | wc -l
```

```bash
zcat annotation.gff3.gz | awk '$3=="gene"' | wc -l
```

```bash
zcat annotation.gff3.gz | grep -v "^#" | cut -f 3 | sort | uniq -c
```

```bash
zcat annotation.gff3.gz | awk 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4,$5,$5-$4+1,$9}' | head
```

---

### 8.3 拓展题

1. 如果一个基因有多个 mRNA，说明什么？
2. 为什么 exon 数量通常远多于 gene 数量？
3. 为什么 CDS 不等于 exon？
4. 为什么有些基因在负链上？
5. 如果一个基因位于负链，它的 start 和 end 是否会反过来写？


## 十、参考网址

- Phytozome: https://phytozome-next.jgi.doe.gov/
- SoyBase: https://www.soybase.org/
- NCBI Genome: https://www.ncbi.nlm.nih.gov/genome/
- Ensembl Plants: https://plants.ensembl.org/
