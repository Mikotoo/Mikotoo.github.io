# RNA-seq 教学演示数据

**这不是真实测序数据。** 参考序列由脚本随机生成，reads 由脚本从这些模拟外显子切取，只用于练习
`hisat2` / `samtools` / `stringtie` 命令，以及观察 SAM/BAM 中的 FLAG、CIGAR、MAPQ、NH 等字段。
请不要把这里的坐标、比对率或表达量当作任何真实的生物学结论。

## 文件

| 文件 | 说明 |
| --- | --- |
| `demo_genome.fa` | 2 条模拟染色体：`demo_chr1`（4000 bp）、`demo_chr2`（2000 bp） |
| `demo_annotation.gtf` | 4 个模拟基因的注释，覆盖正负链与 2–3 个外显子 |
| `demo_reads_R1.fastq`、`demo_reads_R2.fastq` | 双端读长，各 280 条，读长 50 bp |
| `README.md` | 本说明 |

## 基因结构

| 基因 | 位置 | 链 | 外显子（1-based） | 转录本长度 |
| --- | --- | --- | --- | --- |
| `demo_gene1` | `demo_chr1:101-1000` | + | 101-400、701-1000 | 600 bp |
| `demo_gene2` | `demo_chr1:1501-2400` | − | 1501-1800、2101-2400 | 600 bp |
| `demo_gene3` | `demo_chr2:101-1150` | + | 101-350、501-750、901-1150 | 750 bp |
| `demo_gene4` | `demo_chr2:1201-1600` | − | 1201-1375、1451-1600 | 325 bp |

`demo_chr1:3001-3350` 与 `demo_chr2:1650-2000` 是同一段 350 bp 序列的复制，用来制造多重比对。

## reads 里刻意放了什么

| 类型 | 数量 | 比对后用来观察 |
| --- | --- | --- |
| 基因来源，其中每个外显子接头上固定有 2 对跨接头片段 | 每个基因 60 对 | CIGAR 中的 `N`：内含子被跳过 |
| 来自重复区 | 30 对 | 多重比对：`NH:i:2`、`MAPQ` 0 |
| 随机序列，基因组中不存在 | 10 对 | 未比对：FLAG 含 0x4，RNAME/CIGAR 为 `*` |
| 合计 | 280 对 | — |

读长 50 bp，片段长度 150 bp 起，因此双端读长之间通常不相接。质量值随机落在 Q30–Q40。
同一个片段的两条 read 使用相同的 read 名（`demo_00001` 等），符合 HISAT2 对双端文件的要求。

这里的 FASTQ **按明文存放**，方便直接查看。真实测序数据几乎都是 `sample_R1.fastq.gz` 这样的压缩格式，
查看时把 `head`、`wc -l` 换成 `zcat 文件名 | head`、`zcat 文件名 | wc -l` 即可；`hisat2`、`fastp` 都能直接读取 `.gz`。

## 建议的练习顺序

```bash
hisat2-build -p 4 demo_genome.fa demo_index
hisat2 -p 4 -x demo_index -1 demo_reads_R1.fastq -2 demo_reads_R2.fastq -S demo.sam 2> demo_hisat.log
samtools sort -@ 4 -o demo.bam demo.sam
samtools index demo.bam
samtools flagstat demo.bam
samtools view demo.bam | head
```

然后找出带 `N` 的 CIGAR、带 `NH:i:2` 的行，以及用 `samtools view -f 4` 找到的未比对行，
对照讲义第四章（先跑通流程）与第五章（逐字段解释）逐个拆开。
