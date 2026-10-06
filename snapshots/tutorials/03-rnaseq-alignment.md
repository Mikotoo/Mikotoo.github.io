---
title: "第 3 节：RNA-seq 与序列比对"
description: "先动手跑通「FASTQ → 比对 → BAM → 定量」，再对照自己生成的 SAM/BAM 逐个字段读懂输出，最后理解序列比对与剪接比对的原理。"
sidebar:
  label: "第 3 节：RNA-seq 与序列比对"
  order: 3
prev: {"label":"第 2 节：Python 与 R 数据可视化","link":"/tutorials/02-python-r-visualization/"}
next: false
tableOfContents:
  minHeadingLevel: 2
  maxHeadingLevel: 2
---

[← 课程目录](/tutorials/) · 第 3 节 / 共 3 节 · [下载本节原始讲义（Markdown）](/tutorials/bioinformatics/03-rnaseq-alignment.md)

> **练习数据**：[下载RNA-seq 教学演示数据（ZIP）](/tutorials/bioinformatics/rnaseq-demo-data.zip)。解压后得到 `rnaseq_demo_data/`。其中 `demo_reads_R1/R2.fastq.gz` 由脚本从 `demo_genome.fa` 与 `demo_annotation.gtf` 离线生成，**是教学模拟数据**，只用于跑通命令与解读 SAM/BAM 字段。

<details>
<summary>单独查看或下载数据文件</summary>

- [demo_annotation.gtf](/tutorials/bioinformatics/rnaseq_demo_data/demo_annotation.gtf)
- [demo_genome.fa](/tutorials/bioinformatics/rnaseq_demo_data/demo_genome.fa)
- [demo_reads_R1.fastq](/tutorials/bioinformatics/rnaseq_demo_data/demo_reads_R1.fastq)
- [demo_reads_R2.fastq](/tutorials/bioinformatics/rnaseq_demo_data/demo_reads_R2.fastq)
- [README.md](/tutorials/bioinformatics/rnaseq_demo_data/README.md)

</details>

## 一、章节导入

第一课我们学会了看懂参考基因组：`FASTA` 是序列，`GFF3` 是基因的位置和结构。第二课我们学会了把统计结果画成图。

但这两课看的都是「地图」。地图上没有回答一个关键问题：**在某个组织、某个时期，哪些基因正在被表达？**

RNA-seq 就是回答这个问题的常用手段。它的产物不是一条完整的染色体，而是几千万条短序列——**读长（read）**。分析的第一步，就是把这些读长放回参考基因组上，也就是**序列比对（alignment）**。

本课按「先做、再看懂、最后懂原理」的顺序推进：

| 顺序 | 你要做的事 | 对应章节 |
| --- | --- | --- |
| 1 | **先跑一遍**：从 FASTQ 出发，完成比对、排序、索引，一路做到基因表达量，跑完再看自己生成了哪些文件 | 第四章 |
| 2 | **再看懂输出**：对照你刚刚产生的 SAM/BAM，把每一列、每个数字的含义弄清楚 | 第五章 |
| 3 | **最后懂原理**：回到最初的问题——软件凭什么能在几十分钟内把几千万条读长各就各位，为什么 RNA 的比对还要特殊处理 | 第六章 |

这样安排的原因很直接：`FLAG=99`、`CIGAR=21M300N29M` 这类字段，如果先讲定义再去看文件，很容易变成死记硬背；先在自己的结果里见到它们，再回头解释，印象会深得多。算法部分同理——先用过 `hisat2`，再理解它内部做了什么，比凭空想象索引结构要容易。

**所以第四章请先把命令跑通，不要停下来纠结每一列的意思**；遇到看不懂的输出，拍照或复制下来，第五章会逐个解释。

> 课堂提醒：RNA-seq 的读长来自**成熟 mRNA**，而参考基因组是**DNA**。mRNA 已经剪掉了内含子，所以一条读长可能横跨两个外显子。这是 RNA-seq 比对区别于 DNA 重测序比对的核心难点，也是本课的重点。

---

## 二、本章学习目标

完成本章学习后，学生应能够：

1. 说清 RNA-seq 从建库到 FASTQ 的基本流程，理解 read、fragment、insert size 的区别；
2. 读懂 FASTQ 的四行结构，会计算 Phred 质量值，理解 Q20 与 Q30 的实际含义；
3. 掌握 `hisat2-build`、`hisat2`、`samtools`、`stringtie`、`featureCounts` 的基本用法；
4. 能够独立完成「FASTQ → SAM → BAM → 排序 → 索引 → 定量」的完整流程；
5. 会用 `samtools flagstat`、`idxstats`、`view` 检查比对质量并排查常见报错；
6. 逐列读懂 SAM 格式，会用位运算解释 FLAG，会拆解 CIGAR 字符串；
7. 说明 MAPQ 与多重比对的含义，会判断一条读长的比对是否可信；
8. 理解 count、FPKM、TPM 三种表达量的区别与适用场景；
9. 说出序列比对要解决的问题，理解打分、动态规划、种子-扩展、索引四种思路的关系；
10. 解释为什么 DNA 比对器不能直接用于 RNA-seq，理解剪接比对（splice-aware alignment）的原理。

---

## 三、RNA-seq 数据是怎么产生的

### 3.1 从 RNA 到测序读长

一次典型的真核 RNA-seq 建库流程：

```text
提取总 RNA
   ↓  （poly-A 富集 mRNA，或 rRNA 去除）
mRNA
   ↓  （片段化 fragmentation，把长转录本打断成 200–500 bp 小片段）
短片段
   ↓  （反转录 reverse transcription → cDNA）
cDNA
   ↓  （末端修复、加测序接头 adapter、PCR 扩增）
测序文库
   ↓  （上机测序，边合成边测序）
FASTQ 文件（几千万条 read）
```

这里有两个概念要分清：

| 名词 | 含义 |
| --- | --- |
| **fragment / insert** | 打断后的一段 cDNA，也就是「插入片段」，长度约 200–500 bp |
| **read** | 测序仪一次读出的一段碱基，长度 50–150 bp（双端就是两端各读一条） |

所以**双端测序（paired-end）**并不是把一条 fragment 从头读到尾，而是从它的两端各读一小段。两个 read 之间的那段没有被读到的序列，需要靠比对回参考序列来补全。

### 3.2 FASTQ 文件：四行一条记录

FASTQ 是测序数据的标准交付格式，**一条 read 占四行**：

```text
@demo_00001
TGTGTAGCGCGGGGTCGTTCTCTTTGGTGATACTCAAATTCGAGTCCCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

| 行 | 内容 | 说明 |
| --- | --- | --- |
| 第 1 行 | `@` 开头 | read 名称，同一 fragment 的两条 read 名字相同（可能带 `/1`、`/2` 后缀） |
| 第 2 行 | 碱基序列 | 这条 read 读到的 A/T/C/G（偶尔有 N） |
| 第 3 行 | `+` 开头 | 分隔行，可以省略名称 |
| 第 4 行 | 质量字符串 | **与第 2 行逐位对应**，每个字符代表该碱基的测序质量 |

**第 2 行和第 4 行的长度必须完全相等**，这是检查 FASTQ 是否损坏最直接的方法。

双端数据通常交付两个文件：

```text
sample_R1.fastq.gz   # read1，第 1、3、5……条记录
sample_R2.fastq.gz   # read2，与 R1 一一对应
```

两个文件的行数必须一致，第 N 条记录必须来自同一个 fragment。

### 3.3 碱基质量值：Phred 分数与 ASCII 编码

测序仪不只给出碱基，还给出「这个碱基我有多大把握」。把握用**错误率**表示，再换算成 Phred 质量值：

```text
Q = -10 × log10(P)

P 是该碱基测错的概率
```

| Phred 质量 Q | 错误率 P | 通俗说法 | 1 万条 read（150 bp）中的错误碱基数 |
| ---: | ---: | --- | ---: |
| 10 | 1/10 | 很差 | 150000 |
| 20 | 1/100 | 可用 | 15000 |
| 30 | 1/1000 | 常用门槛 | 1500 |
| 40 | 1/10000 | 很好 | 150 |

质量值在文件里不是数字，而是**字符**，规则是：

```text
字符的 ASCII 码 = Q + 33
```

| 字符 | ASCII | Phred Q |
| --- | ---: | ---: |
| `!` | 33 | 0 |
| `5` | 53 | 20 |
| `?` | 63 | 30 |
| `I` | 73 | 40 |

所以上面示例中一排 `I` 表示每个碱基的错误率约为万分之一。**Q30 的比例**（常写作 `%Q30`）是评价一批数据质量最常用的指标，`fastp`、`FastQC` 报告里都会给出。

质量值在后续分析中的作用：比对软件用它区分「真实的错配」和「测序错误」，变异检测软件更是直接依赖它。

### 3.4 双端测序为什么更常用

| 优点 | 说明 |
| --- | --- |
| 更容易唯一比对 | 两端序列同时对上同一个位置，随机巧合的概率大幅降低 |
| 能估计插入片段长度 | 由两个 read 的比对位置推算，可用于检查建库质量 |
| 能发现结构变异 | 两端比对到相距很远的位置，提示缺失、倒位、融合 |
| 定量更可靠 | 片段两端信息互相约束，减少多重比对 |

比对结果里的 `TLEN`（观察到的片段长度）和 `RNEXT`（另一条 read 所在的染色体）就是为这些用途准备的，我们在第五章详细讲。

### 3.5 RNA-seq 数据的四个特点

理解这四点，才能理解后面所有的算法和参数选择：

1. **读长来自外显子**：内含子已被剪接，读长在基因组上的比对是**不连续**的。
2. **存在剪接位点**：一条 50 bp 的 read 可能前 21 bp 在外显子 1、后 29 bp 在外显子 2。
3. **丰度跨度极大**：高表达基因可占全部读长的百分之几，低表达基因可能只有几条 read。
4. **可能有链特异性**：dUTP 等建库方法会保留 RNA 的方向信息，计数时必须选择正确的链参数，否则结果会明显偏低甚至接近 0。

---

## 四、动手实操：从 FASTQ 到表达矩阵

这一章的任务只有一个：**把流程跑通**。

从 FASTQ 开始，依次完成建索引、比对、排序、索引、统计、定量，最后得到一张基因表达表。中途你会看到一堆当时还看不懂的输出——`FLAG`、`CIGAR`、`NH:i:2`、`MAPQ`——**先不要停下来研究它们**，把典型的几行复制到笔记里，第五章会拿着它们逐个字段解释。

跑完之后你应该拥有：

```text
index/demo_index.*.ht2      比对索引（一组文件）
align/demo.sam              比对结果（文本）
align/demo.bam + .bai       排序并建索引后的比对结果（二进制）
logs/demo_hisat.log         比对摘要
logs/demo.flagstat          比对统计
quant/demo.tsv              每个基因的 FPKM / TPM
quant/featureCounts.txt     每个基因的 read 计数
```

### 4.1 两条练习路线

| 路线 | 数据 | 目的 | 是否需要联网 |
| --- | --- | --- | --- |
| **A. 离线快速演示** | `rnaseq_demo_data/`（教学模拟数据，2 条模拟染色体、4 个模拟基因、280 对 read） | 完整跑通命令链，观察 `N`、`NH`、`MAPQ`、FLAG | 不需要 |
| **B. 真实大豆数据** | 第一课下载的大豆 FASTA/GFF3 + 教师提供的 RNA-seq reads | 完成一次真实的转录组定量 | 需要（下载数据） |

> 路线 A 的数据是**教学模拟数据**，序列为脚本随机生成，只用于练习命令和解读字段，**不代表任何真实的生物学结果**。

### 4.2 环境准备

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
  hisat2 samtools stringtie subread fastp gffread -y
conda activate rnaseq
```

逐个确认版本：

```bash
hisat2 --version | head -1
samtools --version | head -1
stringtie --version
featureCounts -v 2>&1 | head -2
```

**为什么要记录版本**：不同版本的索引格式、默认参数可能不同，结果复现必须能说清用的是哪个版本。

### 4.3 建立工作目录

```bash
mkdir -p ~/bioinfo_course/03_rnaseq/{raw,index,align,quant,logs}
cd ~/bioinfo_course/03_rnaseq
```

把数据放进 `raw/`：

```bash
# 路线 A：把演示数据解压到 raw/
unzip rnaseq-demo-data.zip
cp rnaseq_demo_data/demo_genome.fa       raw/
cp rnaseq_demo_data/demo_annotation.gtf  raw/
cp rnaseq_demo_data/demo_reads_R1.fastq  raw/
cp rnaseq_demo_data/demo_reads_R2.fastq  raw/
ls -lh raw/
```

检查 FASTQ 是否完整（最实用的两条命令）：

```bash
wc -l raw/demo_reads_R1.fastq          # 行数应能被 4 整除
head -4 raw/demo_reads_R1.fastq        # 看第一条记录的四行
```

> 演示数据为了便于查看，FASTQ **按明文存放**。真实测序数据几乎都是 `.gz` 压缩格式，把上面两条换成
> `zcat raw/sample_R1.fastq.gz | wc -l` 和 `zcat raw/sample_R1.fastq.gz | head -4` 即可；
> `hisat2`、`fastp` 都能直接读取 `.gz`，不需要先解压。

### 4.4 建立比对索引

```bash
hisat2-build -p 4 raw/demo_genome.fa index/demo_index
```

**参数说明**：

| 参数 | 含义 |
| --- | --- |
| `-p 4` | 使用 4 个线程（按机器调整，不要写满全部核心） |
| 第 1 个位置参数 | 参考基因组 FASTA |
| 第 2 个位置参数 | **索引前缀**（不是文件名）：会生成 `demo_index.1.ht2` … `demo_index.8.ht2` |

生成的索引是**一组文件**，后续比对时用 `-x index/demo_index` 指向这个**前缀**。这是最常见的报错来源 —— 写成了 `-x demo_index.1.ht2` 或 `-x demo_genome.fa`。

真实数据（路线 B）：接第一课下载的 FASTA。

```bash
hisat2-build -p 8 raw/Gmax_508_v4.0.fa index/Wm82
```

### 4.5 质控（可选，但真实数据建议做）

```bash
fastp -i raw/sample_R1.fastq.gz -I raw/sample_R2.fastq.gz \
      -o raw/clean_R1.fastq.gz -O raw/clean_R2.fastq.gz \
      -w 8 -j logs/fastp.json -h logs/fastp.html
```

`fastp` 会去掉接头、低质量碱基和过短 read，并输出 JSON/HTML 报告。打开 HTML 报告关注三件事：

1. `Q30` 比例（一般要求 > 85%）；
2. 是否检测到接头污染；
3. `duplication` 比例（过高提示起始量不足或 PCR 循环过多）。

### 4.6 比对：生成 SAM

```bash
hisat2 -p 4 -x index/demo_index \
       -1 raw/demo_reads_R1.fastq -2 raw/demo_reads_R2.fastq \
       -S align/demo.sam --new-summary 2> logs/demo_hisat.log
```

**参数说明**：

| 参数 | 含义 |
| --- | --- |
| `-x` | 索引**前缀** |
| `-1` / `-2` | 双端 read 文件（顺序必须与建库一致） |
| `-S` | 输出的 SAM 文件 |
| `--new-summary` | 以更易读的格式输出比对摘要 |
| `2> logs/...` | **把摘要和日志重定向出来**。HISAT2 的统计信息走标准错误，不重定向会在屏幕上刷过 |

先看摘要：

```bash
cat logs/demo_hisat.log
```

重点关注：

| 指标 | 含义 | 偏低时的可能原因 |
| --- | --- | --- |
| `Overall alignment rate` | 总比对率 | 参考序列不匹配、污染、接头未去除 |
| `Aligned concordantly` | 双端一致比对率 | RNA 数据正常应占主要部分 |

### 4.7 SAM → BAM → 排序 → 索引

```bash
# 1. SAM 转 BAM（-b 输出 BAM）
samtools view -b -o align/demo.raw.bam align/demo.sam

# 2. 按坐标排序（-@ 指定线程）
samtools sort -@ 4 -o align/demo.bam align/demo.raw.bam

# 3. 建立索引（生成 align/demo.bam.bai）
samtools index align/demo.bam

# 4. 检查文件完整性
samtools quickcheck align/demo.bam && echo "文件完整"
```

也可以一步完成（`samtools sort` 能直接读 SAM）：

```bash
samtools sort -@ 4 -o align/demo.bam align/demo.sam
samtools index align/demo.bam
```

### 4.8 保存比对统计

```bash
samtools flagstat align/demo.bam | tee logs/demo.flagstat
samtools idxstats align/demo.bam > logs/demo.idxstats
```

你会看到 `mapped` 大约是 **96%**，不是 100%——因为演示数据里刻意放了 10 对（20 条）基因组中不存在的读长，正是为了让这一行有内容可看。

`flagstat` 的每一行分别统计什么、`idxstats` 的四列是什么，**第五章会逐一解释**。这里先把输出存进 `logs/`，学完第五章再回来对照一遍，你会发现同样一份输出读起来完全是另一回事。

### 4.9 把三类典型记录留下来

比对的输出成千上万行，但你只需要关注三类**特殊情况**。请把每一类的典型行复制到自己的笔记里，第五章和第六章都会用到它们。

**（1）跨外显子接合点的 read（CIGAR 含 `N`）**

```bash
samtools view align/demo.bam | awk '$6 ~ /N/' | head -5
```

应该能看到类似 `21M300N29M` 的写法：前后两段匹配之间，跳过了 300 bp。
**先记下这一行的 `POS`（第 4 列）和完整的 `CIGAR`（第 6 列）**，为什么中间要跳过去，第六章 6.9 节解释。

**（2）多重比对**

```bash
samtools view align/demo.bam | grep "NH:i:2" | head -6
```

这些行通常 `MAPQ`（第 5 列）为 `0`，并带 `NH:i:2` 标签。

**注意 HISAT2 的默认行为**：默认只报告**一个主比对**，其他候选位置不写进文件，只用 `NH` 标签告知「还有别的位置」。所以同一个 read 名在这里通常**只出现一行**。

想在文件里同时看到两个候选位置，需要显式要求报告多条比对：

```bash
hisat2 -p 4 -x index/demo_index -k 2 \
       -1 raw/demo_reads_R1.fastq -2 raw/demo_reads_R2.fastq \
       -S align/demo.k2.sam --new-summary 2> logs/demo_k2.log
samtools sort -@ 4 -o align/demo.k2.bam align/demo.k2.sam
samtools index align/demo.k2.bam
samtools view align/demo.k2.bam | grep "NH:i:2" | head -6
```

这时同一个 read 名会出现两行，分别落在两条染色体上，其中一条的 `FLAG` 含 `0x100`（256，次要比对），`HI` 标签分别为 1 和 2。**为什么一条 read 会有两个同样好的位置？** 看看数据说明里那段重复序列，第六章 6.11 节解释。

**（3）未比对上的 read**

```bash
samtools view -f 4 align/demo.bam | head -3
```

`-f 4` 表示「FLAG 中含未比对位」，因此无论这条 read 还带了哪些其他状态都能命中。

> 不要写成 `awk '$2==4'`：双端数据中未比对的 read 通常同时带有 `0x1`、`0x40`（或 `0x80`）、`0x8` 等位，FLAG 往往是 77、141 这类数字，**恰好等于 4 的情况很少**。判断是否含某一位，要用位运算（`samtools` 的 `-f`/`-F`）而不是数值相等。

**（4）再随手抽一条普通记录**

```bash
samtools view align/demo.bam | head -1 | awk '{for(i=1;i<=11;i++) printf "%d\t%s\n", i, $i}'
```

上面这些字段现在都还是「一堆数字和字母」。第五章会用你自己跑出来的这几行，把 11 列、`FLAG`、`CIGAR` 逐个拆开。

### 4.10 只保留唯一比对

三种常见写法，按需选择：

```bash
# 写法一：精确筛「只有一条比对」的 read（笔记中的做法）
samtools view -h align/demo.bam | grep -E "^@|NH:i:1" > align/demo.unique.sam
samtools sort -@ 4 -o align/demo.unique.bam align/demo.unique.sam
samtools index align/demo.unique.bam

# 写法二：按 MAPQ 阈值过滤（更快，但阈值含义依赖比对器）
samtools view -b -q 1 -o align/demo.q1.bam align/demo.bam

# 写法三：排除次要比对与补充比对
samtools view -b -F 0x900 -o align/demo.primary.bam align/demo.bam
```

> 注意 `-h`：`samtools view` 默认不输出头部，而 `samtools sort` 需要头部中的 `@SQ` 才能知道参考序列长度，因此**保留头部是必须的**。
>
> 写法一与写法二并不等价：`-q 1` 只是「MAPQ 不小于 1」，而 `NH:i:1` 才是「确实只有一条比对」。两者取舍取决于你对灵敏度与严格度的要求，报告中应写清楚用了哪一种。

### 4.11 定量：StringTie

```bash
stringtie -p 4 -G raw/demo_annotation.gtf -e -B \
          -o quant/demo.gtf -A quant/demo.tsv \
          align/demo.unique.bam
```

**参数说明**：

| 参数 | 含义 |
| --- | --- |
| `-G` | 参考注释（GTF），`-G` 与 `-e` 同时使用时只对已知基因定量，不做新转录本组装 |
| `-e` | 只估计已知转录本的表达量 |
| `-B` | 为每个样本额外输出 Ballgown 格式的表（供后续差异分析使用） |
| `-o` | 输出 GTF |
| `-A` | 输出**基因水平的丰度表**（TPM/FPKM），这是最常用的结果 |

查看结果：

```bash
head -3 quant/demo.tsv
column -t quant/demo.tsv | head -8
```

`-A` 输出的关键列：

| 列名 | 含义 |
| --- | --- |
| `Gene ID` | 基因 ID（与 GTF 的 `gene_id` 对应） |
| `Gene Name` | 基因名 |
| `Reference` | 所在参考序列（染色体） |
| `Strand` | 链 |
| `Start` / `End` | 坐标 |
| `Coverage` | 覆盖该基因的 read 数（可粗略当作 count 使用） |
| `FPKM` | FPKM 值 |
| `TPM` | TPM 值 |

> `Coverage`、`count`、`FPKM`、`TPM` 到底有什么区别、做差异分析该用哪一个，第五章 5.10 节专门讲。这里先记住：**这张表就是本章的最终产出**。

### 4.12 从 GTF 取 count：prepDE.py

差异表达分析需要各样本的 count 矩阵。若已有每个样本的 `-B` 输出，可用 StringTie 附带的脚本汇总：

```bash
# 准备清单：每行「样本名 <Tab> 该样本的 gtf 路径」
printf 'demo\t%s\n' "$PWD/quant/demo.gtf" > quant/prepDE.list

prepDE.py -i quant/prepDE.list -g quant/counts.csv -t quant/transcripts.csv
sed 's/,/\t/g' quant/counts.csv > quant/counts.txt
head -3 quant/counts.txt
```

多样本时把每行都写进 `prepDE.list`，一次生成全部样本的矩阵。

**另一种更直接的做法**（推荐用于已知注释的定量）：

```bash
featureCounts -T 4 -p -a raw/demo_annotation.gtf \
              -o quant/featureCounts.txt align/demo.unique.bam
```

- `-p`：数据是双端，**按 fragment 计数**。忘记加会大致把计数翻倍；
- `-a`：注释文件（GTF/GFF）；
- `-s`：链特异性，`0` 非链特异、`1` 正向、`2` 反向（dUTP 常用 `2`）。**不确定时先用 `0` 试跑，再对照已知基因的表达排查**；
- 输出 `featureCounts.txt.summary` 里给出「已分配 / 未分配」的统计，**未分配比例过高说明注释、链参数或比对有问题**。

### 4.13 合并多样本表达矩阵

拿到每个样本的 TPM 表后，按第一列基因 ID 合并：

```bash
cd quant
for f in *.tsv; do
  s=$(basename "$f" .tsv)
  awk -v s="$s" 'BEGIN{OFS="\t"} FNR==1{next} {print $1, $NF}' "$f" | sort -k1,1 > "$s.tpm"
done
# 逐步合并为矩阵：gene_id + 各样本 TPM
join -t $'\t' -a1 -e NA -o 0,2.2 demo.tpm sample2.tpm > merge.tmp
```

第一课学过的 `cut`、`sort`、`awk`、`join` 在这里都用得上；样本数量多时，改用 Python 的 `pandas`（第二课）更稳妥。

### 4.14 常见报错与排查

这张表是**查阅用**的：跑的过程中报错就来查一行。其中涉及字段含义的原因（如链特异性），会在第五章解释。

| 报错 / 现象 | 原因 | 处理 |
| --- | --- | --- |
| `Could not locate a HISAT2 index with prefix "..."` | `-x` 没指向索引前缀 | 确认 `index/demo_index.*.ht2` 存在，`-x` 只写到 `demo_index` |
| `samtools index: "x.bam" is not sorted` | 未排序或排序中断 | 先 `samtools sort`，检查 `-o` 是否写全 |
| `[E::idx_find_and_load] Could not retrieve index file` | 没有 `.bai` | 运行 `samtools index`（原理见 5.8 节） |
| `fail to read the header from ...` | BAM 缺少头部（常见于 `grep` 过滤时漏了 `-h`） | 过滤 SAM 时加 `-h` 或改用 `samtools view -b` |
| 比对率极低（< 30%） | 参考基因组不匹配 / 污染 / 接头未去 / read 方向错误 | 检查参考版本、跑 `fastp`、确认 `-1`/`-2` 未颠倒 |
| `featureCounts` 计数约为预期一半 | 双端数据漏加 `-p` | 加 `-p` 重新计数 |
| 计数普遍接近 0 | 链特异性参数错误（见 5.10 节） | 改用 `-s 1` 或 `-s 2` 重试，并对照已知高表达基因 |
| `stringtie` 报 `no such file` | `-G` 给的是 GFF3 而非 GTF | 用 `gffread x.gff3 -T -o x.gtf` 转换 |
| `samtools quickcheck` 失败 | 文件被截断（下载/写入中断） | 重新生成或重新下载 |

---

## 五、读懂比对结果：SAM/BAM 文件

上一章你已经把这些文件生成出来了。现在把 `align/demo.sam`、`align/demo.bam` 打开，对着真实内容逐个字段看——**这一章解释的每一个字段，都能在你自己的结果里找到对应的例子**。

### 5.1 从 SAM 到 BAM：为什么需要二进制格式

比对结果的标准格式是 **SAM（Sequence Alignment/Map）**，它是**纯文本**，方便阅读和调试。但一次实验的 SAM 文件动辄几十 GB，因此实际存储和交换使用它的二进制压缩版本 **BAM**。

| 格式 | 特点 |
| --- | --- |
| **SAM** | 纯文本，可以直接 `less`、`grep` 查看 |
| **BAM** | SAM 内容的 **BGZF 压缩**（分块 gzip），体积约为 SAM 的 1/5–1/10，可被程序随机访问 |
| **CRAM** | 参考序列压缩，体积更小，适合长期归档，需要参考基因组才能还原 |
| **`bai` / `.csi`** | BAM 的索引，使程序可以只读取某个区间而不必扫描整个文件 |

**BGZF 的关键意义**：普通 gzip 必须从头解压，而 BGZF 把数据切成独立的块并建立虚拟偏移，因此可以「跳到」任意位置读取。这就是**必须先排序、再建索引**的原因 —— 只有当 read 按坐标有序排列时，索引才能告诉你「染色体 X 第 1000–2000 位的数据在第几块」。

### 5.2 SAM 的两部分

SAM 文件由两部分组成：

```text
头部（header）：每行以 @ 开头，描述参考序列、比对软件、样本等信息
比对记录（alignment）：每行一条 read 的比对结果，制表符分隔，至少 11 列
```

用 `samtools view -H file.bam` 只看头部，`samtools view -h file.bam` 同时输出头部与记录。

### 5.3 头部行

| 标签 | 含义 | 典型内容 |
| --- | --- | --- |
| `@HD` | 文件级信息 | `VN:1.6`（格式版本）、`SO:coordinate`（已按坐标排序） |
| `@SQ` | 参考序列信息，**每个序列一行** | `SN:demo_chr1`（名称）、`LN:4000`（长度） |
| `@RG` | 测序文库与样本信息 | `ID:demo`、`SM:demo_sample`、`PL:ILLUMINA` |
| `@PG` | 生成该文件的程序及命令行 | `PN:hisat2`、`VN:2.2.1`、`CL:"hisat2 -x ..."` |
| `@CO` | 任意注释文本 | 无格式要求 |

**为什么读 BAM 必须能读到头部**：`@SQ` 定义了参考序列，比对记录中的 `RNAME` 和 `POS` 才有意义。如果 BAM 缺少头部或头部与参考不一致，很多工具会直接报错。

`@PG` 记录的命令行非常有用：几个月后回看结果，能立刻知道当时用的是哪个软件、哪些参数。

### 5.4 比对记录的 11 个必选列

这是本课的核心表格，务必逐列理解。

| # | 字段 | 名称 | 含义 |
| ---: | --- | --- | --- |
| 1 | `QNAME` | read 名称 | 与 FASTQ 中的名称一致；双端的两条 read 名字相同 |
| 2 | `FLAG` | 标记 | 一个整数，用二进制位表示多种状态（见 5.5） |
| 3 | `RNAME` | 参考序列名 | 比对到哪条染色体；未比对时为 `*` |
| 4 | `POS` | 位置 | 比对上的**最左**坐标，**1-based**；未比对时为 0 |
| 5 | `MAPQ` | 比对质量 | 见 6.11；255 表示不可用 |
| 6 | `CIGAR` | 比对描述 | 用一串「数字+字母」描述 read 如何贴合参考（见 5.6） |
| 7 | `RNEXT` | 另一条 read 的染色体 | 双端时 mate 所在的染色体；`=` 表示与本条相同；`*` 表示无信息 |
| 8 | `PNEXT` | 另一条 read 的位置 | mate 的最左坐标 |
| 9 | `TLEN` | 插入片段长度 | 两条 read 覆盖的总跨度；**最左的那条为正数，另一条为负数** |
| 10 | `SEQ` | 序列 | 见下方说明 |
| 11 | `QUAL` | 质量 | 与 `SEQ` 逐位对应的 Phred+33 字符 |

关于 `SEQ` 有一个**极易出错**的细节：

> `SEQ` 存的是**与参考序列方向一致**的序列。当 FLAG 含 `0x10`（比对到负链）时，这一列存的是原始 read 的**反向互补**序列。`QUAL` 也按同样的方向**倒序**存放。

也就是说，**不能直接把 SAM 里的 `SEQ` 当成 FASTQ 里读到的原始序列**。需要原始序列时应结合 FLAG 做反向互补还原。

**示例**（`demo_00001` 的一对 read，来自演示数据）：

```text
demo_00001  99   demo_chr1  101  60  50M  =  281   230  TGTGTAGCGCGGGGTCGTTCTCTTTGGTGATACTCAAATTCGAGTCCCAT  IIII...（50 个 I）
demo_00001  147  demo_chr1  281  60  50M  =  101  -230  AAAATAGAGTTCATTTGGCGTGGACGTCCAAAGACCCCAACTTCATTCAA  IIII...（50 个 I）
```

逐项对照：

| 字段 | 第 1 行 | 第 2 行 | 说明 |
| --- | --- | --- | --- |
| `QNAME` | `demo_00001` | `demo_00001` | 同一个 fragment，名字相同 |
| `FLAG` | 99 | 147 | 99 = 成对+正常配对+mate 在负链+第一条；147 = 成对+正常配对+本条负链+第二条 |
| `POS` | 101 | 281 | read1 从 101 开始，read2 从 281 开始 |
| `MAPQ` | 60 | 60 | 唯一比对 |
| `CIGAR` | 50M | 50M | 两条都是 50 bp 连续匹配，没有剪接 |
| `RNEXT` | `=` | `=` | 两条在同一染色体 |
| `PNEXT` | 281 | 101 | 互为对方的 mate 位置 |
| `TLEN` | 230 | −230 | 片段从 101 到 330，长度 230；最左者为正 |

（质量字符串在示例中简写为 50 个 `I`。）

### 5.5 FLAG：一个数字表示一串状态

`FLAG` 是一个位掩码，每一位（bit）代表一种状态。**解读方法**：把数字转成二进制，看哪些位是 1。

| 位（十六进制） | 十进制 | 含义 |
| --- | ---: | --- |
| `0x1` | 1 | 该 read 是成对测序的一条 |
| `0x2` | 2 | 与 mate 形成「正常配对」（方向、距离合理） |
| `0x4` | 4 | **本条 read 未比对** |
| `0x8` | 8 | mate 未比对 |
| `0x10` | 16 | 本条比对到**负链** |
| `0x20` | 32 | mate 比对到负链 |
| `0x40` | 64 | 本条是 read1（第一条） |
| `0x80` | 128 | 本条是 read2（第二条） |
| `0x100` | 256 | **次要比对**（secondary，同一 read 的第 2 条及以后的比对） |
| `0x200` | 512 | 未通过质控（QC fail） |
| `0x400` | 1024 | **PCR 重复**（duplicate） |
| `0x800` | 2048 | 补充比对（supplementary，如跨接头的分段比对） |

常见取值速查：

| FLAG | 拆解 | 含义 |
| ---: | --- | --- |
| 0 | — | 单端、正链、已比对 |
| 4 | `0x4` | 未比对 |
| 16 | `0x10` | 单端、负链、已比对 |
| 99 | 1+2+32+64 | 双端正常配对，read1，mate 在负链 |
| 147 | 1+2+16+128 | 双端正常配对，read2，本条在负链 |
| 256 | `0x100` | 次要比对（常与 MAPQ 0 同时出现） |
| 1024 | `0x400` | 标记为 PCR 重复 |
| 2048 | `0x800` | 补充比对 |

**验证 99**：`99 = 64 + 32 + 2 + 1 = 0x40 + 0x20 + 0x02 + 0x01`，即 read1、正常配对、mate 在负链、成对测序 —— 与示例完全一致。

**验证 147**：`147 = 128 + 16 + 2 + 1 = 0x80 + 0x10 + 0x02 + 0x01`，即 read2、本条在负链、正常配对、成对测序，是上面 read1 的 mate。

**常用过滤写法**（`samtools` 支持十进制和十六进制）：

```bash
samtools view -F 4     demo.bam    # 排除未比对的（-F 表示「去掉含该位的」）
samtools view -F 0x900 demo.bam    # 排除次要比对与补充比对（256 + 2048）
samtools view -f 2     demo.bam    # 只看正常配对的（-f 表示「必须含该位」）
samtools view -c -F 4  demo.bam    # 统计已比对的 read 数
```

### 5.6 CIGAR：read 是怎么贴合参考序列的

`CIGAR` 由若干「数字 + 操作符」组成，例如 `50M`、`21M300N29M`、`10S40M`。

| 操作符 | 名称 | 是否消耗参考 | 是否消耗 read | 含义 |
| --- | --- | :---: | :---: | --- |
| `M` | alignment match | 是 | 是 | 比对匹配。**注意：可以包含错配** |
| `=` | sequence match | 是 | 是 | 完全匹配（碱基相同） |
| `X` | sequence mismatch | 是 | 是 | 错配（碱基不同） |
| `I` | insertion | 否 | 是 | read 中有、参考中没有的碱基 |
| `D` | deletion | 是 | 否 | 参考中有、read 中没有的碱基 |
| `N` | skipped region | 是 | 否 | **跳过参考上的这一段（如内含子）** |
| `S` | soft clip | 否 | 是 | read 两端未参与比对的碱基，**序列仍保留在 SEQ 中** |
| `H` | hard clip | 否 | 否 | 被截掉的碱基，**SEQ 中不含这些碱基** |
| `P` | padding | 否 | 否 | 占位，极少使用 |

两个必须记住的规则：

1. **`SEQ` 的长度 = `M + I + S + = + X` 之和**（`D`、`N`、`H` 不消耗 read）。
   例如 `21M300N29M` → `SEQ` 长度为 50。
2. **`I`/`D` 与 `N` 的区别**：`I`/`D` 通常很小（几 bp），多来自测序错误或真实小变异；`N` 可以长达几千甚至几十万 bp，在 RNA-seq 中代表**内含子**。

**实例对照**（演示数据中的一条跨接合点 read）：

```text
demo_junction01  0  demo_chr1  380  60  21M300N29M  *  0  0  CGTGGAACTGGCCTGCCAACTTGTACTTACAAAGTTGTCGTACATGTGTC  IIII...  NH:i:1  NM:i:0
```

- `POS = 380`：比对起点；
- `21M`：380–400 共 21 bp 与外显子 1 匹配；
- `300N`：跳过 401–700，正是 300 bp 的内含子；
- `29M`：701–729 与外显子 2 匹配。

`M` 中的错配数量可以从可选标签 `NM`（编辑距离）和 `MD`（错配描述）读到，见下一节。

### 5.7 可选标签

`11` 列之后是可选标签，格式为 `TAG:TYPE:VALUE`。HISAT2 常用以下标签：

| 标签 | 类型 | 含义 |
| --- | --- | --- |
| `NH:i:` | 整数 | 报告的比对条数（1 = 唯一） |
| `HI:i:` | 整数 | 当前是第几个候选 |
| `NM:i:` | 整数 | 编辑距离（错配 + 插入 + 缺失碱基数） |
| `MD:Z:` | 字符串 | 错配与匹配的压缩描述，变异检测常用 |
| `AS:i:` | 整数 | 比对得分 |
| `XS:i:` | 整数 | 次优比对得分（与 `AS` 接近说明该位置不够唯一） |
| `RG:Z:` | 字符串 | 所属的 read group，与头部 `@RG` 的 `ID` 对应 |
| `YT:Z:` | 字符串 | 链特异性状态（`UU`/`CP`/`UC` 等），HISAT2 特有 |

**判断多重比对**最直接的方式就是看 `NH`：

```bash
samtools view demo.bam | grep -c "NH:i:1"      # 唯一比对条数
samtools view demo.bam | awk '$5==0' | wc -l   # MAPQ 为 0 的行数
```

### 5.8 BAM 索引与「必须排序」

```text
sample.sam  ──samtools sort──▶  sample.bam（按坐标排序）──samtools index──▶  sample.bam.bai
```

- `samtools sort` 把比对记录按 `RNAME` + `POS` 排序。**排序要求正是由索引决定的**；
- `samtools index` 生成 `.bai` 索引，之后 `samtools view sample.bam demo_chr1:101-400` 才能只读取该区间；
- 未排序的 BAM 仍可 `view`、`flagstat`，但 `index` 会报错，区间查询也会退化为全表扫描；
- `@HD` 中的 `SO:coordinate` 就是「已排序」的标记，可用 `samtools view -H` 检查。

### 5.9 samtools 常用命令速查

| 命令 | 用途 |
| --- | --- |
| `samtools view -b -o out.bam in.sam` | SAM 转 BAM |
| `samtools view -h in.bam \| head` | 查看头部与前几条记录 |
| `samtools view -H in.bam` | 只查看头部 |
| `samtools view -c in.bam` | 统计记录数 |
| `samtools view -F 4 -c in.bam` | 统计已比对条数 |
| `samtools view in.bam demo_chr1:101-400` | 提取某区间的比对 |
| `samtools sort -@ 4 -o sorted.bam in.bam` | 按坐标排序 |
| `samtools index sorted.bam` | 建立索引（生成 `.bai`） |
| `samtools flagstat in.bam` | 汇总比对统计（总数、比对率、配对率、重复率） |
| `samtools idxstats in.bam` | 每条参考序列上的比对条数 |
| `samtools stats in.bam` | 生成详细统计报告 |
| `samtools quickcheck in.bam` | 快速检查文件是否完整（截断文件会失败） |
| `samtools faidx genome.fa` | 为 FASTA 建索引 |
| `samtools depth in.bam` | 逐位点覆盖深度 |

**`flagstat` 报告怎么读**（各项含义）：

| 行 | 含义 |
| --- | --- |
| `total` | 文件中的总记录数（**注意：包含次要比对，可能大于 read 总数**） |
| `primary` | 主比对条数 |
| `secondary` | 次要比对条数 |
| `supplementary` | 补充比对条数 |
| `duplicates` | 被标记为 PCR 重复的条数 |
| `mapped (%)` | 已比对比例 |
| `properly paired (%)` | 正常配对比例（双端数据的核心质量指标） |
| `singletons (%)` | 只有一条 read 比对上的比例 |

### 5.10 从 BAM 到表达量

BAM 告诉我们「每条 read 比对到哪里」，但生物学问题问的是「每个基因有多少表达量」。中间还需要一步：**计数（counting / quantification）**。

| 输出 | 含义 | 用途 |
| --- | --- | --- |
| **count（原始计数）** | 落在该基因上的 read（或 fragment）条数 | 差异表达分析（DESeq2、edgeR）的**唯一推荐输入** |
| **FPKM** | 每千碱基每百万比对 read 的 fragments | 同一样本内、不同基因间比较尚可；**跨样本比较不可靠** |
| **TPM** | transcripts per million，先按长度归一化、再按总量归一化 | 同一样本内的基因间比较，跨样本比较优于 FPKM |

$$\text{FPKM} = \frac{\text{比对到该基因的 fragment 数}}{\text{基因长度 (kb)} \times \text{总比对 fragment 数 (百万)}}$$

$$\text{TPM} = \frac{\text{reads per kilobase}}{\sum(\text{reads per kilobase})} \times 10^6$$

**为什么差异分析要用 count 而不是 TPM**：DESeq2 等工具需要原始计数来建模「测序深度」和「生物学变异」的离散分布；TPM/FPKM 已经做过除法，会破坏统计模型所需的方差结构。

常用计数工具：

| 工具 | 特点 |
| --- | --- |
| `featureCounts`（subread） | 快、显存小，直接吃 GTF/GFF 与 BAM；`-p` 用于双端，`-s` 指定链特异性 |
| `htseq-count` | 经典、直观，速度较慢 |
| `StringTie -e -B -A` | 与转录本组装同一套流程，可同时给出 TPM/FPKM |
| `prepDE.py`（StringTie 附带） | 从各样本 GTF 汇总成 count 矩阵，供 DESeq2 使用 |

**链特异性参数是新手最容易出错的地方**：dUTP 建库的 read1 通常与转录本反向，计数时 `-s 2`（reverse）；若误设 `-s 0`，计数会明显偏低甚至接近 0。遇到「结果异常少」时，应优先怀疑链参数。

### 5.11 回到你在第四章跑出的结果

现在把第四章存下来的输出重新读一遍，应该都能对上号了：

| 你保存的输出 | 现在怎么读它 |
| --- | --- |
| `logs/demo.flagstat` | 逐行对照 5.9 节的 `flagstat` 表：`mapped 96.43%` 对应那 20 条刻意放进来的读长；`0 duplicates` 是因为 HISAT2 不标记重复 |
| `logs/demo.idxstats` | 四列分别为参考序列名、长度、该序列上的比对条数、未比对条数。**用处**：大量 read 堆在少数序列上提示污染或重复序列；某条染色体几乎没有 read，提示该染色体组装或注释有问题 |
| 带 `N` 的 `CIGAR` 行 | 5.6 节的 `N` 操作符：跳过的那段就是内含子（原理见第六章 6.9 节） |
| `NH:i:2` 且 `MAPQ 0` 的行 | 5.7 节的 `NH` 标签：这条 read 在基因组上有两个同样好的位置，也就是 `demo_chr1:3001-3350` 与 `demo_chr2:1650-2000` 那段重复序列（原理见第六章 6.11 节） |
| `-f 4` 找到的行 | 5.4 节的 `*`／`0` 取值：未比对时 `RNAME`、`CIGAR` 无法给出有意义的值 |
| `quant/demo.tsv` | 5.10 节的 TPM/FPKM 两列，就是你在第四章最后生成的表达量 |
| `quant/featureCounts.txt.summary` | `Unassigned_*` 各行分别对应未比对、多重比对、落在基因间区等情况，都是前几节讲过的概念 |

如果哪一行还对不上，回到对应小节查一遍——**能用自己的数据解释每一行，比记住定义有用得多**。

---

## 六、比对算法原理：软件内部发生了什么

到这里，你已经用过 `hisat2`，也读懂了它输出的 SAM/BAM。现在回到最初的问题：**`hisat2` 到底做了什么，能在几十分钟内把几千万条 read 放到 10 亿 bp 的基因组上？**

先回想上一章你亲眼见到的两件事，这一章就是解释它们：

- 为什么有些 read 的 `CIGAR` 中间会出现一长段 `N`？（见 6.9 剪接比对）
- 为什么有些 read 的 `MAPQ` 是 0、还带着 `NH:i:2`？（见 6.11 多重比对）

### 6.1 比对要解决什么问题

给定：

- **参考序列** `R`，长度 `n`（大豆基因组约 10 亿 bp，人类约 31 亿 bp）；
- **查询序列**（read）`Q`，长度 `m`（通常 50–150 bp）；

要求：找出 `Q` 在 `R` 中**最可能出现的位置**，并说明它是怎么对上的（哪些碱基匹配、哪些错配、哪里插入了空位）。

这是一个「在超长文本中做近似字符串搜索」的问题，难点有三个：

1. 允许**错配**和**空位**（测序错误、真实变异、剪接都会造成不完全匹配）；
2. 参考序列极长，不可能对每个位置都做一次完整比较；
3. 数据量极大，一次实验几千万到上亿条 read。

### 6.2 基本术语

| 术语 | 含义 |
| --- | --- |
| **match / mismatch** | 匹配 / 错配：两个位置碱基相同 / 不同 |
| **gap（空位）** | 一方有碱基、另一方没有。在 read 中出现叫插入（I），在参考中出现叫缺失（D） |
| **identity** | 一致率，通常指比对区间内匹配碱基数 ÷ 比对长度 |
| **coverage（覆盖度）** | 参考序列某个位置被多少条 read 覆盖 |
| **seed（种子）** | 用于快速定位的短精确匹配片段 |
| **reference bias** | 参考序列本身缺失某段序列时，来自该段的 read 无法比对上的现象 |

### 6.3 打分：怎么判断「哪个比对更好」

比对软件需要一个统一标准来比较不同方案，这就是**打分函数（scoring scheme）**。最常用的一组：

```text
匹配   match    +1
错配   mismatch -1
空位   gap      -2（或：开空位 -5，延伸空位 -2）
```

第三种写法叫**仿射空位罚分（affine gap penalty）**，因为它更符合生物学：一次造成 3 bp 的缺失，通常是一个事件，不应该按 3 个独立空位同等惩罚。所以「开一个空位」罚得重，「把已有空位延长」罚得轻。

### 6.4 动态规划：Needleman–Wunsch 全局比对

**动态规划（dynamic programming, DP）** 是序列比对的经典解法。以全局比对 Needleman–Wunsch（1970）为例。

设 `F(i, j)` 为「`Q` 的前 `i` 个碱基」与「`R` 的前 `j` 个碱基」的最优得分，递推式为：

```text
F(i, j) = max(
    F(i-1, j-1) + s(Qi, Rj)      # 两个碱基对齐（匹配或错配）
    F(i-1, j)   + g              # Qi 对应一个空位
    F(i, j-1)   + g              # Rj 对应一个空位
)

边界：F(i, 0) = i × g      F(0, j) = j × g      F(0, 0) = 0
```

其中 `s()` 是匹配/错配得分，`g` 是空位罚分。

**小例子**：`Q = GATT`，`R = GAT`，匹配 +1、错配 −1、空位 −1。填表：

|  | j=0 | j=1 (G) | j=2 (A) | j=3 (T) |
| --- | ---: | ---: | ---: | ---: |
| **i=0** | 0 | −1 | −2 | −3 |
| **i=1 (G)** | −1 | **1** | 0 | −1 |
| **i=2 (A)** | −2 | 0 | **2** | 1 |
| **i=3 (T)** | −3 | 1 | 1 | **3** |
| **i=4 (T)** | −4 | 0 | 0 | **2** |

从右下角 `F(4,3) = 2` 回溯，得到的最优比对是：

```text
G A T T
G A - T
```

3 个匹配（+3）加 1 个空位（−1）= 2 分。

**要点**：

- 表格的每个格子只依赖左、上、左上三个格子，因此可以逐行填充；
- 时间复杂度 `O(m × n)`，空间复杂度 `O(m × n)`；
- 回溯（traceback）决定最终的比对写法，这正是后面 `CIGAR` 字符串的来源。

### 6.5 Smith–Waterman 局部比对

Needleman–Wunsch 是**全局比对**，要求两条序列从头到尾都参与比对，适合长度相近的序列。

而 read 比对是**局部**的：read 只对应参考序列里的某一段，其余部分完全无关。Smith–Waterman（1981）的改动很简洁：

```text
F(i, j) = max( 0,                       # 关键：不允许出现负分
               F(i-1, j-1) + s(Qi, Rj),
               F(i-1, j)   + g,
               F(i, j-1)   + g )
```

- 引入 `0` 作为下限，任何负分部分都被「截断」；
- 回溯从**全表最大值**开始，遇到 0 停止，得到的就是局部最优比对。

### 6.6 为什么不能直接对全基因组做动态规划

假设用 Smith–Waterman 把一条 150 bp 的 read 比对到 10 亿 bp 的大豆基因组：

```text
150 × 1,000,000,000 = 1.5 × 10^11 次格子计算   ← 仅一条 read
```

一次实验有 3000 万条 read，这是天文数字。所以实际比对软件采用了三类**加速策略**：

| 策略 | 思路 | 代表 |
| --- | --- | --- |
| 索引（indexing） | 预先把参考序列加工成可快速查找的结构，避免扫描全部位置 | BWA、Bowtie2、HISAT2、STAR |
| 种子-扩展（seed-and-extend） | 先用短精确匹配快速定位候选位置，只对候选位置做精确比对 | BLAST、BWA-MEM |
| 启发式（heuristic） | 提前放弃明显不可能更优的方向，牺牲理论最优换取速度 | BLAST 的 X-drop、带状 DP |

### 6.7 种子与扩展：从 BLAST 到 BWA

**种子-扩展**是最直观的加速思路，分两步：

```text
第一步（种子）：把 read 切成若干短的 k-mer（例如 15 bp），
                在参考序列索引里快速找到这些 k-mer 出现的位置。
第二步（扩展）：只在这些位置附近做带罚分的精确比对（DP），
                得到最终得分和 CIGAR。
```

BLAST 是最早把这一思路工程化的工具，用「两次命中」规则过滤随机匹配，并用 **E-value**（在随机序列中期望出现该得分的次数）衡量显著性。E-value 越小越可信，这也是第一课里 `blastn`、`blastp` 常用 `-evalue 1e-5` 的原因。

后续的 BWA、Bowtie2 把「种子」做得更精细（如 BWA-MEM 的**超级最大精确匹配 SMEM**），并配合高效的索引结构，才达到今天的速度。

### 6.8 索引：把「搜索」变成「查表」

问题：如何在一部 10 亿字的书里，瞬间找到某个 15 字短语的所有出现位置？

答案：**提前把所有后缀排好序**。

| 结构 | 含义 | 作用 |
| --- | --- | --- |
| **后缀数组（suffix array, SA）** | 把参考序列所有后缀按字典序排序后的起始位置数组 | 查找一个短串只需二分查找，`O(m log n)` |
| **BWT（Burrows–Wheeler 变换）** | 把序列做循环移位排序后取最后一列 | 让相同上下文聚在一起，便于压缩和检索 |
| **FM-index** | BWT + 秩表（rank/checkpoint）+ 抽样后缀数组 | 支持**后向搜索**，`O(m)` 完成查找，且体积远小于原序列 |

关键点：**FM-index 让「在基因组中查找」变成了一次字符接一个字符的查表操作**，不需要保存整个基因组的多份拷贝，因此索引体积可控（人类基因组约几个 GB），查询速度极快。BWA、Bowtie2、HISAT2 都属于这一家族。

`hisat2-build` 做的事情，就是把 `genome.fa` 编译成这种索引 —— 这也解释了为什么它耗时长、占用空间大，而之后的每次比对都很快。

### 6.9 RNA-seq 的剪接比对：为什么 DNA 比对器不够用

这是本课最重要的一个知识点，也是你在第四章 4.9 节看到 `21M300N29M` 这种 `CIGAR` 的原因。

回顾 3.5：read 来自成熟 mRNA，内含子已被剪掉。例如演示数据中 `demo_gene1` 的结构：

```text
参考基因组（demo_chr1）:
  外显子1              内含子（300 bp）          外显子2
  101 ── 400                                    701 ── 1000
        └──────────── 剪接 ────────────────────────┘

成熟 mRNA:
  外显子1（300 bp） + 外显子2（300 bp）
```

如果有一条 read 正好跨越接合点（junction），它的两半在基因组上相隔 300 bp：

```text
read（50 bp）:  CGTGGAACTGGCCTGCCAACT | TGTACTTACAAAGTTGTCGTACATGTGTC
                       ↑ 21 bp                    ↑ 29 bp
基因组位置:        380─400                 701─729
                                    中间跳过 401─700（内含子）
```

- **DNA 比对器**（如 BWA-MEM）默认只允许很短的 gap，遇到这种情况只能把这条 read 局部比对上（soft-clip），或者干脆判为未比对；
- **剪接感知比对器**（splice-aware aligner，如 HISAT2、STAR）会把 read 切成几段分别比对，再把它们「缝」在一起，用一个**长 gap** 表示被跳过的内含子。

这个「长 gap」在 SAM 文件里写成 CIGAR 的 **`N`** 操作符。上面这条 read 的 CIGAR 就是：

```text
21M300N29M
```

读作：前 21 bp 与参考匹配（`M`），**跳过参考序列上的 300 bp（`N`）**，最后 29 bp 继续匹配。

因此有两条实用结论：

1. **做 RNA-seq 必须用剪接感知比对器**，用 BWA 直接比对会大量损失跨接头的 read；
2. **看到 `N` 不要当成错误**，它往往正是内含子，反而是「这条 read 来自成熟 mRNA」的证据。

真核生物的内含子通常以 `GT` 开头、`AG` 结尾（GT–AG 规则），高质量比对器会优先选择符合该规则的剪接位点；如果提供 GTF 注释，`HISAT2` 还可以直接从注释中提取已知剪接位点，进一步提高准确率。

### 6.10 HISAT2 与 STAR 的思路差异

两者都是主流剪接比对器，但权衡不同：

| | **HISAT2** | **STAR** |
| --- | --- | --- |
| 索引结构 | 层级式图 FM-index（GFM index）：全基因组索引 + 大量局部索引，把剪接位点、SNP 作为「图」上的分支 | 全基因组**后缀数组** + 最大可映射前缀（MMP）搜索 |
| 内存占用 | 较小（人类基因组约 4–5 GB） | 很大（人类基因组常需 30 GB 量级） |
| 速度 | 快 | 更快 |
| 剪接处理 | 在索引层面表达已知/候选剪接位点 | 先做种子比对，再通过 MMP 把种子「缝合」成长比对 |
| 常见场景 | 内存有限、需要稳定可控 | 服务器内存充足、追求速度 |

> 参数与内存需求随版本、基因组和参数变化，实际部署时应以所用版本的官方文档和本机测试为准。

### 6.11 比对质量 MAPQ 与多重比对

一条 read 可能在基因组上有多个同样好的位置，例如它来自重复序列、旁系同源基因或保守结构域。这种情况叫**多重比对（multi-mapping）**——你在第四章 4.9 节筛出的那些 `MAPQ 0`、带 `NH:i:2` 的记录就是这一类，演示数据里的重复区正是刻意为此准备的。

| 字段 | 含义 |
| --- | --- |
| **MAPQ** | 比对质量，`Q = -10 × log10(P)`，`P` 是「该比对位置是错的」的概率。取值 0 表示无法区分；255 表示该值不可用 |
| **NH** | 该 read 在文件中报告的比对条数（`NH:i:2` 表示有 2 条） |
| **HI** | 当前这条是第几个候选（从 1 开始） |
| **secondary / supplementary** | FLAG 中的标记位，表示「同一 read 的第 2 条比对」/「补充比对（如跨接头的分段比对）」 |

实践中常见的处理方式：

| 策略 | 做法 | 代价 |
| --- | --- | --- |
| 只保留唯一比对 | 过滤掉 `NH > 1` 的 read | 简单可靠，但会丢弃来自多拷贝基因的真实信息 |
| 按 MAPQ 阈值过滤 | 只保留 `MAPQ ≥ 1` 或 `≥ 10` | 更快，但阈值含义依赖比对器，不如 `NH` 直接 |
| 保留并分配权重 | 按 1/NH 分配计数 | 保留信息，但需要额外脚本，且影响可比性 |

第 4.10 节的实操里会看到，第一篇笔记里的做法正是用 `grep "NH:i:1"` 精确筛出唯一比对。这个做法简单、含义明确，代价是要多扫描一遍 SAM 文件。

### 6.12 常用比对工具怎么选

| 工具 | 主要用途 | 是否剪接感知 |
| --- | --- | --- |
| **BWA-MEM** | DNA 重测序、变异检测的默认选择 | 否 |
| **Bowtie2** | 短读长、ChIP-seq、小基因组 | 否 |
| **HISAT2** | RNA-seq 剪接比对，内存友好 | **是** |
| **STAR** | RNA-seq 剪接比对，追求速度与灵敏度 | **是** |
| **minimap2** | 长读长（PacBio / Nanopore）；也支持 `-ax sr` 短读长和 `-ax splice` 剪接模式 | 部分 |

选择原则：**先看数据类型，再看资源**。DNA 用 BWA-MEM，RNA 用 HISAT2 或 STAR；基因组小、内存紧就用 HISAT2。

---

## 七、课堂练习

练习按「先跑命令、再读字段、最后定量」的顺序排列，与前四、五、六章一致：练习 1–5 走流程，练习 6–10 拆解 SAM/BAM 字段，练习 11–14 做过滤、定量与交叉验证。

### 练习 1：检查 FASTQ 的完整性

```bash
wc -l raw/demo_reads_R1.fastq
wc -l raw/demo_reads_R2.fastq
```

**问题**：两个文件的行数是否相同？是否都能被 4 整除？一共多少对 read？

---

### 练习 2：解读质量值

```bash
head -4 raw/demo_reads_R1.fastq
```

**问题**：

1. 第 4 行第 1 个字符对应的 Phred 质量是多少？
2. 换算成错误率是多少？
3. 用 `awk` 统计整个文件中质量字符的 ASCII 范围：

```bash
awk 'NR%4==0' raw/demo_reads_R1.fastq | fold -w1 | sort -u | tr '\n' ' '
```

---

### 练习 3：建立索引并理解「前缀」

```bash
hisat2-build -p 4 raw/demo_genome.fa index/demo_index
ls -lh index/
```

**问题**：

1. 一共生成了几个索引文件？总大小是多少？
2. 为什么比对时 `-x` 要写 `index/demo_index` 而不是某个具体文件？
3. 索引文件比原始 FASTA 大还是小？

---

### 练习 4：完成比对并读懂摘要

```bash
hisat2 -p 4 -x index/demo_index \
       -1 raw/demo_reads_R1.fastq -2 raw/demo_reads_R2.fastq \
       -S align/demo.sam --new-summary 2> logs/demo_hisat.log
cat logs/demo_hisat.log
```

**问题**：

1. 总比对率是多少？
2. 为什么不会是 100%？（提示：看看演示数据里放了什么）
3. 如果这是真实数据，比对率 60% 时你会先检查什么？

---

### 练习 5：SAM 转 BAM 并统计

```bash
samtools sort -@ 4 -o align/demo.bam align/demo.sam
samtools index align/demo.bam
samtools flagstat align/demo.bam
```

**问题**：

1. `total` 是多少？与「read 总数」是否相同？
2. `properly paired` 的比例是多少？
3. `duplicates` 为什么是 0？（提示：HISAT2 是否标记重复？）

---

### 练习 6：逐列拆解一条记录

```bash
samtools view align/demo.bam | head -1 | awk '{for(i=1;i<=11;i++) printf "%d\t%s\n", i, $i}'
```

**问题**：把 11 列的名称逐一写出，并解释这条记录的 FLAG 是怎样加出来的。

---

### 练习 7：用位运算验证 FLAG

```bash
samtools view align/demo.bam | awk '{print $2}' | sort -n | uniq -c | sort -rn | head
```

**问题**：

1. 出现最多的 FLAG 是哪几个？
2. 把它们拆成二进制位，说明每种状态；
3. 下面两条记录分别代表什么？

```bash
samtools view align/demo.bam | awk '$2==99' | head -1
samtools view align/demo.bam | awk '$2==147' | head -1
```

---

### 练习 8：找出跨接合点的 read

```bash
samtools view align/demo.bam | awk '$6 ~ /N/' | head -5
```

**问题**：

1. 这些 CIGAR 中的 `N` 前面和后面的数字分别是多少？
2. 用 `POS + 前面的数字` 和 `POS + 前面的数字 + N + 后面的数字` 算出两段所在的坐标，与 `demo_annotation.gtf` 中的外显子坐标对照；
3. 验证 `SEQ` 长度是否等于 `M + I + S` 之和：

```bash
samtools view align/demo.bam | awk '$6 ~ /N/ {print length($10), $6}' | head -3
```

---

### 练习 9：观察多重比对

先看默认输出中的多重比对：

```bash
samtools view align/demo.bam | grep "NH:i:2" | head -6
```

再要求 HISAT2 报告全部候选位置，重新比对一次：

```bash
hisat2 -p 4 -x index/demo_index -k 2 \
       -1 raw/demo_reads_R1.fastq -2 raw/demo_reads_R2.fastq \
       -S align/demo.k2.sam --new-summary 2> logs/demo_k2.log
samtools sort -@ 4 -o align/demo.k2.bam align/demo.k2.sam
samtools index align/demo.k2.bam
samtools view align/demo.k2.bam | grep "NH:i:2" | head -6
```

**问题**：

1. 默认输出中，这些行的 `MAPQ` 是多少？为什么？同一个 read 名出现了几行？
2. 加上 `-k 2` 后，同一个 read 名出现了几行？分别落在哪条参考序列、哪个坐标？
3. 其中一行的 `FLAG` 为什么含 `0x100`？`HI` 标签的两个值分别是什么？
4. 用 `samtools idxstats` 的输出，说明为什么重复区会产生这种情况：

```bash
samtools idxstats align/demo.bam
```

---

### 练习 10：找到未比对上的 read

```bash
samtools view -f 4 align/demo.bam | head -3
samtools view -c -f 4 align/demo.bam
```

**问题**：

1. 未比对记录的 `RNAME`、`POS`、`CIGAR` 分别是什么？与第 5.4 节的说明是否一致？
2. 这些记录的 `FLAG` 具体是哪些数值？把它们拆成二进制位，说明除了「未比对」还标了什么；
3. 为什么不能用 `awk '$2==4'` 来找它们？

---

### 练习 11：只保留唯一比对并比较

```bash
samtools view -h align/demo.bam | grep -E "^@|NH:i:1" > align/demo.unique.sam
samtools sort -@ 4 -o align/demo.unique.bam align/demo.unique.sam
samtools index align/demo.unique.bam
samtools view -c align/demo.bam
samtools view -c align/demo.unique.bam
```

**问题**：

1. 过滤后少了多少条？这些正是多重比对与未比对的 read；
2. 如果改成 `samtools view -b -q 1`，结果是否相同？为什么？（提示：注意两者判定标准不同）

---

### 练习 12：定量并解读 TPM

```bash
stringtie -p 4 -G raw/demo_annotation.gtf -e -B \
          -o quant/demo.gtf -A quant/demo.tsv align/demo.unique.bam
column -t quant/demo.tsv | head
```

**问题**：

1. 表中有哪几列？哪一列最接近「表达量」？
2. 4 个模拟基因的 TPM 分别是多少？它们的总和应该是多少？
3. 为什么 count 与 TPM 的排序可能不同？

---

### 练习 13：用 featureCounts 交叉验证

```bash
featureCounts -T 4 -p -a raw/demo_annotation.gtf \
              -o quant/featureCounts.txt align/demo.unique.bam
head -3 quant/featureCounts.txt
cat quant/featureCounts.txt.summary
```

**问题**：

1. `Assigned` 与 `Unassigned_Unmapped` 分别是多少？
2. 两个工具的计数是否一致？不一致时可能是什么原因？
3. 去掉 `-p` 再跑一次，计数如何变化？

---

### 练习 14：读懂 `flagstat` 的每一项

```bash
samtools flagstat align/demo.bam > logs/demo.flagstat
samtools stats align/demo.bam > logs/demo.stats
grep -E "^SN" logs/demo.stats | head -20
```

**问题**：从 `stats` 中找出「平均读长」「插入片段平均长度」「比对率」，并说明它们与 `flagstat` 的哪些行对应。

---

## 八、综合练习

### 综合任务 1：完成一次完整的 RNA-seq 定量

按以下顺序，从原始 FASTQ 走到基因表达表，并把每一步的命令和关键输出记录下来：

1. 检查 FASTQ 完整性；
2. 建立（或复用）参考基因组索引；
3. 质控（真实数据必做，模拟数据可跳过）；
4. 比对并保存日志；
5. SAM → BAM → 排序 → 索引；
6. 用 `flagstat`、`idxstats` 检查比对质量；
7. 过滤唯一比对；
8. 定量（StringTie 与 featureCounts 各做一次）；
9. 整理成表达表。

### 综合任务 2：填写比对质量报告表

| 项目 | 结果 |
| --- | --- |
| 使用的参考基因组文件 |  |
| 索引文件个数与总大小 |  |
| 使用的比对软件与版本 |  |
| read 总对数 |  |
| 总比对率 |  |
| 双端一致比对率 |  |
| 多重比对条数 |  |
| 未比对条数 |  |
| 带 `N` 的 CIGAR 条数 |  |
| 定量工具与参数（含链参数） |  |
| 表达量最高的基因 |  |
| 该基因的 count / TPM |  |

### 综合任务 3：排查一个「坏结果」

下面是某同学的报告，请指出至少三处可疑之处并给出排查步骤：

```text
比对率：31%
properly paired：8%
featureCounts Assigned：4%
表达量最高基因 TPM：0.7
```

> 提示：考虑参考基因组版本、read 文件顺序、接头污染、链特异性参数、是否漏加 `-p`。

---

## 九、课后思考题

### 9.1 概念题

1. 为什么 RNA-seq 的 read 可以在基因组上「不连续」地比对？
2. CIGAR 中的 `N`、`D`、`I` 分别在什么情况下出现？为什么 `N` 可以很长，而 `D` 通常很短？
3. `M` 与 `=` 有什么区别？为什么大多数比对器默认输出 `M` 而不是 `=`/`X`？
4. FLAG 为 `83` 的 read 处于什么状态？（提示：`83 = 64 + 16 + 2 + 1`）
5. MAPQ 0 与 MAPQ 60 的区别是什么？为什么 `NH:i:1` 比 `MAPQ > 0` 更适合判断唯一比对？
6. 为什么差异表达分析要用 count，而不能直接用 TPM？
7. FPKM 与 TPM 的区别在哪里？做跨样本比较时哪个更合适？
8. 为什么 `samtools index` 要求 BAM 先按坐标排序？
9. `SEQ` 列在 FLAG 含 `0x10` 时与 FASTQ 中的原始序列是什么关系？
10. 为什么说 `fastp` 报告的 `%Q30` 比「平均质量值」更能反映数据质量？

### 9.2 命令题

解释以下命令的含义，并说明每条命令的输出里应该重点看什么：

```bash
hisat2-build -p 8 genome.fa index/Wm82
```

```bash
hisat2 -p 8 -x index/Wm82 -1 clean_R1.fq.gz -2 clean_R2.fq.gz -S sample.sam --new-summary 2> sample.log
```

```bash
samtools sort -@ 8 -o sample.bam sample.sam && samtools index sample.bam
```

```bash
samtools flagstat sample.bam
```

```bash
samtools view -h sample.bam | grep -E "^@|NH:i:1" | samtools view -b -o unique.bam -
```

```bash
featureCounts -T 8 -p -s 2 -a annotation.gtf -o counts.txt unique.bam
```

```bash
stringtie -p 8 -G annotation.gtf -e -B -o sample.gtf -A sample.tsv unique.bam
```

### 9.3 拓展题

1. 一条 150 bp 的 read 覆盖了三个外显子，它的 CIGAR 可能长什么样？请写出一个合理的例子并解释每一段。
2. 如果参考基因组缺少某个基因的第二个外显子（组装缺失），来自该外显子的 read 会得到什么样的比对结果？这会怎样影响该基因的表达量估计？
3. 同样是多重比对，来自「重复序列」和来自「旁系同源基因」的 read 在生物学含义上有什么不同？处理策略是否应该一样？
4. 为什么用 STAR 需要 30 GB 内存，而 HISAT2 只要几个 GB？用第 6.8 节的索引概念解释这个差距。
5. 如果一个基因的 TPM 很高但 count 很低，可能是什么原因？（提示：基因长度）
6. 长读长测序（PacBio/Nanopore）的 RNA-seq 为什么可以绕开剪接比对的部分难题？它带来了什么新问题？

---

## 十、参考网址

- HISAT2 官方文档：<https://daehwankimlab.github.io/hisat2/>
- STAR 官方仓库：<https://github.com/alexdobin/STAR>
- SAM 格式规范（SAMv1）：<https://samtools.github.io/hts-specs/SAMv1.pdf>
- SAM 格式说明（htslib）：<https://www.htslib.org/doc/sam.html>
- samtools 文档：<https://www.htslib.org/doc/samtools.html>
- StringTie：<https://ccb.jhu.edu/software/stringtie/>
- featureCounts / subread：<https://subread.sourceforge.net/>
- Bowtie2：<https://bowtie-bio.sourceforge.net/bowtie2/>
- BWA：<https://bio-bwa.sourceforge.net/>
- fastp：<https://github.com/OpenGene/fastp>
- FastQC：<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
- ENA（欧洲核苷酸档案，下载公共测序数据）：<https://www.ebi.ac.uk/ena/browser/home>
- NCBI SRA：<https://www.ncbi.nlm.nih.gov/sra>
- Phytozome（植物参考基因组）：<https://phytozome-next.jgi.doe.gov/>
