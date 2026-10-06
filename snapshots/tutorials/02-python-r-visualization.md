---
title: "第 2 节：Python 与 R 数据可视化"
description: "把上一节的统计结果转化为柱状图、箱线图、热图与散点图，并使用 Lotus 示例数据练习基因组圈图。"
sidebar:
  label: "第 2 节：Python 与 R 数据可视化"
  order: 2
prev: {"label":"第 1 节：Linux 与参考基因组","link":"/tutorials/01-linux-genome/"}
next: false
tableOfContents:
  minHeadingLevel: 2
  maxHeadingLevel: 2
---

[← 课程目录](/tutorials/) · 第 2 节 / 共 2 节 · [下载本节原始讲义（Markdown）](/tutorials/bioinformatics/02-python-r-visualization.md)

> **练习数据**：[下载 Lotus 圈图数据包（ZIP）](/tutorials/bioinformatics/lotus-circos-data.zip)。解压后得到 `lotus_circos_data/` 文件夹，放在本节工作目录中即可。基础绘图仍使用第一节大豆数据生成的统计表；这个数据包用于后面的 `circlize` 进阶练习。

<details>
<summary>单独查看或下载圈图数据文件</summary>

- [Cent.txt](/tutorials/bioinformatics/lotus_circos_data/Cent.txt)
- [Chr_length.txt](/tutorials/bioinformatics/lotus_circos_data/Chr_length.txt)
- [DEL_density.txt](/tutorials/bioinformatics/lotus_circos_data/DEL_density.txt)
- [GC.0.5M.bed](/tutorials/bioinformatics/lotus_circos_data/GC.0.5M.bed)
- [gene_density.txt](/tutorials/bioinformatics/lotus_circos_data/gene_density.txt)
- [GifuT2T_MG20T2T_collinear.txt](/tutorials/bioinformatics/lotus_circos_data/GifuT2T_MG20T2T_collinear.txt)
- [INS_density.txt](/tutorials/bioinformatics/lotus_circos_data/INS_density.txt)
- [SNP_density.txt](/tutorials/bioinformatics/lotus_circos_data/SNP_density.txt)
- [TE_density.txt](/tutorials/bioinformatics/lotus_circos_data/TE_density.txt)

</details>

## 一、章节导入

在上一节课中，我们学习了如何从 Phytozome 下载大豆参考基因组文件，并使用 Linux 基础命令查看 FASTA 和 GFF3 注释文件。同时，我们还学习了使用 `grep`、`awk`、`cut`、`sort`、`uniq`、`wc` 等命令对基因组文件进行初步统计，例如：

- FASTA 文件中有多少条序列；
- GFF3 文件中有多少种注释类型；
- GFF3 文件中有多少个 gene、mRNA、exon、CDS；
- 每条染色体上有多少个基因；
- 每个基因的长度是多少；
- 每条染色体上不同 feature 类型的数量是多少。

这些统计结果如果只用表格展示，虽然准确，但不够直观。在生物信息学研究、课程汇报和论文写作中，我们通常会将统计结果转化为图形。图形可以帮助我们更快地发现数据规律，例如：

- 哪条染色体基因最多；
- 基因长度主要集中在哪个范围；
- 不同染色体上的基因长度分布是否相似；
- 染色体长度和基因数量是否相关；
- GFF3 注释文件中哪类 feature 数量最多。

本节课将基于上一节课得到的大豆基因组统计结果，学习如何分别使用 **Python** 和 **R** 绘制常见的生物信息学图形，包括柱状图、箱线图、热图、饼图、直方图、密度图、散点图，以及进阶的基因组圈图。

---

## 二、本章学习目标

完成本章学习后，学生应能够：

1. 理解生物信息学可视化的基本目的；
2. 知道不同图形适合展示什么类型的数据；
3. 掌握 Python 中 `pandas`、`numpy`、`matplotlib`、`scipy` 的基础用法；
4. 掌握 R 中 `tidyverse`、`pheatmap` 和 `circlize` 的基础用法；
5. 能够根据基因组统计表绘制柱状图、箱线图、饼图、直方图、密度图、热图和散点图；
6. 能够将图形保存为 PDF 或 PNG 文件；
7. 能够对绘图结果进行简单解释。

---

## 三、可视化前的数据准备

本节课继续使用上一节课下载的大豆参考基因组文件。

假设我们有两个文件：

```bash
genome.fa.gz
annotation.gff3.gz
```

其中：

```text
genome.fa.gz        大豆参考基因组 FASTA 文件
annotation.gff3.gz  大豆基因组 GFF3 注释文件
```

如果你的实际文件名不同，可以使用软链接统一名称。例如：

```bash
ln -s 你的基因组文件.fa.gz genome.fa.gz
ln -s 你的注释文件.gff3.gz annotation.gff3.gz
```

这样后续命令就可以统一使用 `genome.fa.gz` 和 `annotation.gff3.gz`。

---

## 四、建立第二课工作目录

```bash
mkdir -p ~/bioinfo_course/02_visualization
cd ~/bioinfo_course/02_visualization
```

如果第一课的数据保存在上一节课目录中，可以复制或链接过来：

```bash
ln -s ../01_linux_genome/genome.fa.gz .
ln -s ../01_linux_genome/annotation.gff3.gz .
```

查看当前目录：

```bash
ls -lh
```

---

## 五、生成绘图所需的统计表

在正式使用 Python 和 R 画图之前，我们先使用 Linux 命令从 FASTA 和 GFF3 文件中生成几个简单的统计表。

这些表格是本节课所有图形的输入文件。

---

### 1. 统计 GFF3 文件中不同 feature 类型的数量

GFF3 文件的第 3 列表示 feature 类型，例如：

- gene；
- mRNA；
- exon；
- CDS；
- five_prime_UTR；
- three_prime_UTR。

生成统计表：

```bash
echo -e "feature\tcount" > feature_count.tsv

zcat annotation.gff3.gz | \
awk '$0 !~ /^#/ {count[$3]++} END{for (f in count) print f"\t"count[f]}' | \
sort -k1,1 >> feature_count.tsv
```

查看结果：

```bash
cat feature_count.tsv
```

结果示例：

```text
feature count
CDS     300000
exon    320000
gene    55000
mRNA    60000
```

这个文件可以用于绘制：

- 柱状图；
- 饼图。

---

### 2. 统计每条染色体上的基因数量

```bash
echo -e "chr\tgene_count" > chr_gene_count.tsv

zcat annotation.gff3.gz | \
awk '$0 !~ /^#/ && $3=="gene"{count[$1]++} END{for (chr in count) print chr"\t"count[chr]}' | \
sort -k1,1V >> chr_gene_count.tsv
```

查看结果：

```bash
head chr_gene_count.tsv
```

结果示例：

```text
chr     gene_count
Gm01    3000
Gm02    2800
Gm03    2600
```

这个文件可以用于绘制：

- 柱状图；
- 散点图。

---

### 3. 计算每个基因的长度

基因长度计算公式为：

```text
gene_length = end - start + 1
```

生成基因长度表：

```bash
zcat annotation.gff3.gz | \
awk 'BEGIN{OFS="\t"; print "gene_id","chr","start","end","strand","length"} 
$0 !~ /^#/ && $3=="gene" {
    id=$9;
    sub(/.*ID=/,"",id);
    sub(/;.*/,"",id);
    print id,$1,$4,$5,$7,$5-$4+1
}' > gene_length.tsv
```

查看结果：

```bash
head gene_length.tsv
```

结果示例：

```text
gene_id           chr     start   end     strand  length
Glyma.01G000100   Gm01    1000    5000    +       4001
Glyma.01G000200   Gm01    8000    9500    -       1501
```

这个文件可以用于绘制：

- 直方图；
- 密度图；
- 箱线图。

---

### 4. 统计每条染色体上不同 feature 的数量

```bash
echo -e "chr\tfeature\tcount" > chr_feature_count.tsv

zcat annotation.gff3.gz | \
awk '$0 !~ /^#/ {count[$1"\t"$3]++} END{for (k in count) print k"\t"count[k]}' | \
sort -k1,1V -k2,2 >> chr_feature_count.tsv
```

查看结果：

```bash
head chr_feature_count.tsv
```

结果示例：

```text
chr     feature count
Gm01    CDS     15000
Gm01    exon    16000
Gm01    gene    3000
Gm01    mRNA    3300
Gm02    CDS     14000
```

这个文件可以用于绘制：

- 热图。

---

### 5. 统计每条染色体长度

从 FASTA 文件中统计每条序列长度：

```bash
zcat genome.fa.gz | \
awk 'BEGIN{OFS="\t"; print "chr","length"} 
/^>/ {
    if(NR>1) print chr,len;
    chr=substr($1,2);
    len=0;
    next
}
{
    len+=length($0)
}
END{
    print chr,len
}' > chr_length.tsv
```

查看结果：

```bash
head chr_length.tsv
```

这个文件可以与 `chr_gene_count.tsv` 结合，用于绘制：

- 染色体长度与基因数量关系散点图。

---

### 6. 本节课需要的输入文件汇总

完成以上命令后，当前目录中应包含以下统计表：

```text
feature_count.tsv
chr_gene_count.tsv
gene_length.tsv
chr_feature_count.tsv
chr_length.tsv
```

它们分别表示：

| 文件名 | 内容 | 可用于绘制的图形 |
|---|---|---|
| `feature_count.tsv` | 不同 feature 类型数量 | 柱状图、饼图 |
| `chr_gene_count.tsv` | 每条染色体上的基因数量 | 柱状图、散点图 |
| `gene_length.tsv` | 每个基因的长度 | 直方图、密度图、箱线图 |
| `chr_feature_count.tsv` | 每条染色体上不同 feature 数量 | 热图 |
| `chr_length.tsv` | 每条染色体长度 | 散点图 |

---

## 六、不同图形适合展示什么数据

不同类型的数据适合用不同的图形表示。

| 图形 | 适合展示的数据 | 本节课示例 |
|---|---|---|
| 柱状图 | 不同类别之间的数量比较 | 每条染色体基因数量 |
| 饼图 | 各部分占整体的比例 | GFF3 中不同 feature 的比例 |
| 直方图 | 连续变量的分布 | 基因长度分布 |
| 密度图 | 连续变量的平滑分布趋势 | 基因长度密度分布 |
| 箱线图 | 不同组之间连续变量分布比较 | 不同染色体上的基因长度 |
| 热图 | 矩阵型数据 | 染色体 × feature 数量矩阵 |
| 散点图 | 两个连续变量之间的关系 | 染色体长度与基因数量关系 |

需要注意：

- 柱状图用于类别变量之间的比较；
- 直方图用于连续数值变量的分布；
- 箱线图适合比较不同分组的数据分布；
- 热图适合展示二维矩阵；
- 散点图适合观察两个数值变量之间是否相关。

---

## 七、Python 绘图模块

### 1. 安装 Python 环境

如果使用 conda，可以建立一个专门的绘图环境：

```bash
conda create -n bio_plot python=3.10 pandas numpy matplotlib scipy -y
conda activate bio_plot
```

如果使用 mamba，可以使用：

```bash
mamba create -n bio_plot python=3.10 pandas numpy matplotlib scipy -y
mamba activate bio_plot
```

本节课 Python 部分主要使用以下包：

| Python 包 | 作用 |
|---|---|
| `pandas` | 读取和整理表格 |
| `numpy` | 数值计算和 log 转换 |
| `matplotlib` | 基础绘图 |
| `scipy` | 用于绘制密度曲线 |

---

### 2. Python 脚本运行方式

本节课的 Python 代码不建议直接一行一行粘贴到交互式 Python 中运行，而是统一采用下面的方式：

```text
第一步：编辑一个 .py 脚本文件
第二步：把完整代码写入脚本
第三步：在命令行中使用 python 脚本名.py 执行
第四步：检查输出的 PDF 或 PNG 图片
```

例如：

```bash
nano plot_python_check_data.py
```

写好代码后保存退出，然后运行：

```bash
python plot_python_check_data.py
```

如果使用 VS Code、Notepad++ 或其他编辑器，也可以在本地编辑好 `.py` 文件，再上传到服务器运行。

---

### 3. Python 读取表格

新建脚本：

```bash
nano plot_python_check_data.py
```

写入以下内容：

```python
import pandas as pd

feature = pd.read_csv("feature_count.tsv", sep="\t")
chr_gene = pd.read_csv("chr_gene_count.tsv", sep="\t")
gene_len = pd.read_csv("gene_length.tsv", sep="\t")
chr_feature = pd.read_csv("chr_feature_count.tsv", sep="\t")
chr_length = pd.read_csv("chr_length.tsv", sep="\t")

print("feature_count.tsv:")
print(feature.head())

print("\nchr_gene_count.tsv:")
print(chr_gene.head())

print("\ngene_length.tsv:")
print(gene_len.head())
```

运行脚本：

```bash
python plot_python_check_data.py
```

这一步的目的是确认 Python 能够正常读取前面生成的 `.tsv` 表格。

---

### 4. Python 绘制柱状图：不同 feature 类型数量

新建脚本：

```bash
nano plot_python_feature_barplot.py
```

写入以下内容：

```python
import pandas as pd
import matplotlib.pyplot as plt

feature = pd.read_csv("feature_count.tsv", sep="\t")

plt.figure(figsize=(6, 4))
plt.bar(feature["feature"], feature["count"])
plt.xlabel("Feature type")
plt.ylabel("Count")
plt.title("Number of different feature types in GFF3")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_feature_barplot.pdf")
plt.savefig("python_feature_barplot.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_feature_barplot.py
```

运行后会生成：

```text
python_feature_barplot.pdf
python_feature_barplot.png
```

---

### 5. Python 绘制柱状图：每条染色体上的基因数量

新建脚本：

```bash
nano plot_python_chr_gene_barplot.py
```

写入以下内容：

```python
import pandas as pd
import matplotlib.pyplot as plt

chr_gene = pd.read_csv("chr_gene_count.tsv", sep="\t")

plt.figure(figsize=(8, 4))
plt.bar(chr_gene["chr"], chr_gene["gene_count"])
plt.xlabel("Chromosome")
plt.ylabel("Gene number")
plt.title("Gene number on each chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_chr_gene_barplot.pdf")
plt.savefig("python_chr_gene_barplot.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_chr_gene_barplot.py
```

这个图可以回答：

```text
不同染色体上的基因数量是否相同？
哪条染色体上的基因最多？
哪条染色体上的基因最少？
```

---

### 6. Python 绘制直方图：基因长度分布

新建脚本：

```bash
nano plot_python_gene_length_histogram.py
```

写入以下内容：

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

gene_len = pd.read_csv("gene_length.tsv", sep="\t")

# 普通基因长度直方图
plt.figure(figsize=(6, 4))
plt.hist(gene_len["length"], bins=50)
plt.xlabel("Gene length")
plt.ylabel("Gene count")
plt.title("Distribution of gene length")
plt.tight_layout()
plt.savefig("python_gene_length_histogram.pdf")
plt.savefig("python_gene_length_histogram.png", dpi=300)
plt.close()

# log10 转换后的基因长度直方图
gene_len["log10_length"] = np.log10(gene_len["length"])

plt.figure(figsize=(6, 4))
plt.hist(gene_len["log10_length"], bins=50)
plt.xlabel("log10(Gene length)")
plt.ylabel("Gene count")
plt.title("Distribution of log10 gene length")
plt.tight_layout()
plt.savefig("python_gene_length_log_histogram.pdf")
plt.savefig("python_gene_length_log_histogram.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_gene_length_histogram.py
```

由于有些基因特别长，普通坐标下可能不容易看清主体分布，因此常对基因长度进行 `log10` 转换。

---

### 7. Python 绘制密度图：基因长度分布趋势

新建脚本：

```bash
nano plot_python_gene_length_density.py
```

写入以下内容：

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

gene_len = pd.read_csv("gene_length.tsv", sep="\t")
x = np.log10(gene_len["length"])

density = gaussian_kde(x)
x_grid = np.linspace(x.min(), x.max(), 200)

plt.figure(figsize=(6, 4))
plt.plot(x_grid, density(x_grid))
plt.xlabel("log10(Gene length)")
plt.ylabel("Density")
plt.title("Density plot of gene length")
plt.tight_layout()
plt.savefig("python_gene_length_density.pdf")
plt.savefig("python_gene_length_density.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_gene_length_density.py
```

密度图比直方图更平滑，适合展示整体分布趋势。

---

### 8. Python 绘制箱线图：不同染色体基因长度比较

新建脚本：

```bash
nano plot_python_gene_length_boxplot.py
```

写入以下内容：

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

gene_len = pd.read_csv("gene_length.tsv", sep="\t")
gene_len["log10_length"] = np.log10(gene_len["length"])

chr_list = sorted(gene_len["chr"].unique())

data = [
    gene_len.loc[gene_len["chr"] == chr_name, "log10_length"]
    for chr_name in chr_list
]

plt.figure(figsize=(10, 4))
plt.boxplot(data, labels=chr_list, showfliers=False)
plt.xlabel("Chromosome")
plt.ylabel("log10(Gene length)")
plt.title("Gene length distribution on each chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_gene_length_boxplot.pdf")
plt.savefig("python_gene_length_boxplot.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_gene_length_boxplot.py
```

其中：

```python
showfliers=False
```

表示隐藏极端离群点，使图中主体分布更清楚。

---

### 9. Python 绘制饼图：不同 feature 类型比例

新建脚本：

```bash
nano plot_python_feature_piechart.py
```

写入以下内容：

```python
import pandas as pd
import matplotlib.pyplot as plt

feature = pd.read_csv("feature_count.tsv", sep="\t")

plt.figure(figsize=(6, 6))
plt.pie(
    feature["count"],
    labels=feature["feature"],
    autopct="%1.1f%%",
    startangle=90
)
plt.title("Proportion of feature types in GFF3")
plt.tight_layout()
plt.savefig("python_feature_piechart.pdf")
plt.savefig("python_feature_piechart.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_feature_piechart.py
```

饼图适合展示整体组成比例，但如果类别太多，饼图会变得不清楚。

---

### 10. Python 绘制热图：染色体 × feature 数量矩阵

新建脚本：

```bash
nano plot_python_chr_feature_heatmap.py
```

写入以下内容：

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

chr_feature = pd.read_csv("chr_feature_count.tsv", sep="\t")

matrix = chr_feature.pivot(
    index="chr",
    columns="feature",
    values="count"
).fillna(0)

matrix_log = np.log10(matrix + 1)

plt.figure(figsize=(8, 6))
plt.imshow(matrix_log, aspect="auto")
plt.colorbar(label="log10(count + 1)")
plt.xticks(
    ticks=np.arange(len(matrix_log.columns)),
    labels=matrix_log.columns,
    rotation=45,
    ha="right"
)
plt.yticks(
    ticks=np.arange(len(matrix_log.index)),
    labels=matrix_log.index
)
plt.xlabel("Feature type")
plt.ylabel("Chromosome")
plt.title("Chromosome-feature count heatmap")
plt.tight_layout()
plt.savefig("python_chr_feature_heatmap.pdf")
plt.savefig("python_chr_feature_heatmap.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_chr_feature_heatmap.py
```

这里对数据进行了：

```python
log10(count + 1)
```

转换，是因为不同 feature 数量差异可能很大，直接画图时颜色可能被少数大值主导。

---

### 11. Python 绘制散点图：染色体长度与基因数量关系

新建脚本：

```bash
nano plot_python_chr_length_gene_scatter.py
```

写入以下内容：

```python
import pandas as pd
import matplotlib.pyplot as plt

chr_gene = pd.read_csv("chr_gene_count.tsv", sep="\t")
chr_length = pd.read_csv("chr_length.tsv", sep="\t")

df = pd.merge(chr_length, chr_gene, on="chr", how="inner")
df["length_Mb"] = df["length"] / 1_000_000

plt.figure(figsize=(6, 4))
plt.scatter(df["length_Mb"], df["gene_count"])

for _, row in df.iterrows():
    plt.text(row["length_Mb"], row["gene_count"], row["chr"], fontsize=8)

plt.xlabel("Chromosome length (Mb)")
plt.ylabel("Gene number")
plt.title("Relationship between chromosome length and gene number")
plt.tight_layout()
plt.savefig("python_chr_length_gene_scatter.pdf")
plt.savefig("python_chr_length_gene_scatter.png", dpi=300)
plt.close()
```

运行脚本：

```bash
python plot_python_chr_length_gene_scatter.py
```

这个图可以用来观察：

```text
染色体越长，基因数量是否越多？
是否存在染色体很长但基因数量并不多的情况？
```

---

## 八、R 绘图模块

### 1. 安装 R 包

R 语言中常用的生物信息学绘图包包括 `ggplot2`、`dplyr`、`readr`、`tidyr` 等。实际上，这些包都包含在 `tidyverse` 这个集合包中。

因此，本节课不需要分别安装 `ggplot2`、`dplyr`、`readr` 和 `tidyr`，只需要安装 `tidyverse` 即可。普通热图使用 `pheatmap`，进阶圈图使用 `circlize`，这两个包需要额外安装。

进入 R：

```bash
R
```

安装包：

```r
install.packages(c("tidyverse", "pheatmap", "circlize"))
```


说明：

| R 包 | 作用 |
|---|---|
| `tidyverse` | 包含 `ggplot2`、`dplyr`、`readr`、`tidyr` 等常用数据整理和绘图包 |
| `pheatmap` | 用于绘制热图 |
| `circlize` | 用于绘制基因组圈图、环形热图、共线性连线图等 |

`tidyverse` 中常用组件包括：

| 组件 | 主要作用 |
|---|---|
| `ggplot2` | 绘图 |
| `dplyr` | 数据筛选、分组、汇总 |
| `readr` | 读取表格文件 |
| `tidyr` | 数据长宽格式转换 |
| `tibble` | 更友好的数据框格式 |

---

### 2. R 脚本运行方式

本节课的 R 代码不建议直接一行一行粘贴到 R 交互界面中运行，而是统一采用下面的方式：

```text
第一步：编辑一个 .R 脚本文件
第二步：把完整代码写入脚本
第三步：在命令行中使用 Rscript 脚本名.R 执行
第四步：检查输出的 PDF 或 PNG 图片
```

例如：

```bash
nano plot_R_check_data.R
```

写好代码后保存退出，然后运行：

```bash
Rscript plot_R_check_data.R
```

---

### 3. R 读取表格

新建脚本：

```bash
nano plot_R_check_data.R
```

写入以下内容：

```r
library(tidyverse)

feature <- read_tsv("feature_count.tsv", show_col_types = FALSE)
chr_gene <- read_tsv("chr_gene_count.tsv", show_col_types = FALSE)
gene_len <- read_tsv("gene_length.tsv", show_col_types = FALSE)
chr_feature <- read_tsv("chr_feature_count.tsv", show_col_types = FALSE)
chr_length <- read_tsv("chr_length.tsv", show_col_types = FALSE)

cat("feature_count.tsv:\n")
print(head(feature))

cat("\nchr_gene_count.tsv:\n")
print(head(chr_gene))

cat("\ngene_length.tsv:\n")
print(head(gene_len))
```

运行脚本：

```bash
Rscript plot_R_check_data.R
```

这一步用于确认 R 能正常读取 `.tsv` 表格。

---

### 4. R 绘制柱状图：不同 feature 类型数量

新建脚本：

```bash
nano plot_R_feature_barplot.R
```

写入以下内容：

```r
library(tidyverse)

feature <- read_tsv("feature_count.tsv", show_col_types = FALSE)

p <- ggplot(feature, aes(x = feature, y = count)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Feature type",
    y = "Count",
    title = "Number of different feature types in GFF3"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("R_feature_barplot.pdf", p, width = 6, height = 4)
ggsave("R_feature_barplot.png", p, width = 6, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_feature_barplot.R
```

---

### 5. R 绘制柱状图：每条染色体上的基因数量

新建脚本：

```bash
nano plot_R_chr_gene_barplot.R
```

写入以下内容：

```r
library(tidyverse)

chr_gene <- read_tsv("chr_gene_count.tsv", show_col_types = FALSE)

p <- ggplot(chr_gene, aes(x = chr, y = gene_count)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = "Gene number",
    title = "Gene number on each chromosome"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("R_chr_gene_barplot.pdf", p, width = 8, height = 4)
ggsave("R_chr_gene_barplot.png", p, width = 8, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_chr_gene_barplot.R
```

---

### 6. R 绘制直方图：基因长度分布

新建脚本：

```bash
nano plot_R_gene_length_histogram.R
```

写入以下内容：

```r
library(tidyverse)

gene_len <- read_tsv("gene_length.tsv", show_col_types = FALSE)

# 普通基因长度直方图
p <- ggplot(gene_len, aes(x = length)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(
    x = "Gene length",
    y = "Gene count",
    title = "Distribution of gene length"
  )

ggsave("R_gene_length_histogram.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_histogram.png", p, width = 6, height = 4, dpi = 300)

# log10 转换后的基因长度直方图
gene_len <- gene_len %>%
  mutate(log10_length = log10(length))

p <- ggplot(gene_len, aes(x = log10_length)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(
    x = "log10(Gene length)",
    y = "Gene count",
    title = "Distribution of log10 gene length"
  )

ggsave("R_gene_length_log_histogram.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_log_histogram.png", p, width = 6, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_gene_length_histogram.R
```

---

### 7. R 绘制密度图：基因长度分布趋势

新建脚本：

```bash
nano plot_R_gene_length_density.R
```

写入以下内容：

```r
library(tidyverse)

gene_len <- read_tsv("gene_length.tsv", show_col_types = FALSE) %>%
  mutate(log10_length = log10(length))

p <- ggplot(gene_len, aes(x = log10_length)) +
  geom_density() +
  theme_bw() +
  labs(
    x = "log10(Gene length)",
    y = "Density",
    title = "Density plot of gene length"
  )

ggsave("R_gene_length_density.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_density.png", p, width = 6, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_gene_length_density.R
```

---

### 8. R 绘制箱线图：不同染色体基因长度比较

新建脚本：

```bash
nano plot_R_gene_length_boxplot.R
```

写入以下内容：

```r
library(tidyverse)

gene_len <- read_tsv("gene_length.tsv", show_col_types = FALSE) %>%
  mutate(log10_length = log10(length))

p <- ggplot(gene_len, aes(x = chr, y = log10_length)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = "log10(Gene length)",
    title = "Gene length distribution on each chromosome"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("R_gene_length_boxplot.pdf", p, width = 10, height = 4)
ggsave("R_gene_length_boxplot.png", p, width = 10, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_gene_length_boxplot.R
```

其中：

```r
outlier.shape = NA
```

表示不显示极端离群点，使箱线图主体分布更清楚。

---

### 9. R 绘制饼图：不同 feature 类型比例

新建脚本：

```bash
nano plot_R_feature_piechart.R
```

写入以下内容：

```r
library(tidyverse)

feature <- read_tsv("feature_count.tsv", show_col_types = FALSE)

p <- ggplot(feature, aes(x = "", y = count, fill = feature)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(
    title = "Proportion of feature types in GFF3",
    fill = "Feature"
  )

ggsave("R_feature_piechart.pdf", p, width = 6, height = 6)
ggsave("R_feature_piechart.png", p, width = 6, height = 6, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_feature_piechart.R
```

说明：R 中的饼图通常是通过柱状图加极坐标转换得到的：

```r
geom_col() + coord_polar(theta = "y")
```

---

### 10. R 绘制热图：染色体 × feature 数量矩阵

新建脚本：

```bash
nano plot_R_chr_feature_heatmap.R
```

写入以下内容：

```r
library(tidyverse)
library(pheatmap)

chr_feature <- read_tsv("chr_feature_count.tsv", show_col_types = FALSE)

matrix_data <- chr_feature %>%
  pivot_wider(
    names_from = feature,
    values_from = count,
    values_fill = 0
  )

mat <- as.data.frame(matrix_data)
rownames(mat) <- mat$chr
mat$chr <- NULL

mat_log <- log10(as.matrix(mat) + 1)

pdf("R_chr_feature_heatmap.pdf", width = 8, height = 6)
pheatmap(
  mat_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Chromosome-feature count heatmap"
)
dev.off()

png("R_chr_feature_heatmap.png", width = 2400, height = 1800, res = 300)
pheatmap(
  mat_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Chromosome-feature count heatmap"
)
dev.off()
```

运行脚本：

```bash
Rscript plot_R_chr_feature_heatmap.R
```

这里使用 `pivot_wider()` 将长格式数据转换成矩阵格式。

原始表格是长格式：

```text
chr     feature count
Gm01    gene    3000
Gm01    exon    16000
Gm02    gene    2800
```

热图需要的是矩阵格式：

```text
        gene    exon    CDS
Gm01    3000    16000   15000
Gm02    2800    15000   14000
```

---

### 11. R 绘制散点图：染色体长度与基因数量关系

新建脚本：

```bash
nano plot_R_chr_length_gene_scatter.R
```

写入以下内容：

```r
library(tidyverse)

chr_gene <- read_tsv("chr_gene_count.tsv", show_col_types = FALSE)
chr_length <- read_tsv("chr_length.tsv", show_col_types = FALSE)

df <- inner_join(chr_length, chr_gene, by = "chr") %>%
  mutate(length_Mb = length / 1000000)

p <- ggplot(df, aes(x = length_Mb, y = gene_count, label = chr)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw() +
  labs(
    x = "Chromosome length (Mb)",
    y = "Gene number",
    title = "Relationship between chromosome length and gene number"
  )

ggsave("R_chr_length_gene_scatter.pdf", p, width = 6, height = 4)
ggsave("R_chr_length_gene_scatter.png", p, width = 6, height = 4, dpi = 300)
```

运行脚本：

```bash
Rscript plot_R_chr_length_gene_scatter.R
```

---

## 九、Python 综合练习

### 练习 1：读取数据

使用 Python 读取以下文件，并输出前 5 行：

```text
feature_count.tsv
chr_gene_count.tsv
gene_length.tsv
chr_feature_count.tsv
chr_length.tsv
```

要求使用：

```python
pd.read_csv()
```

---

### 练习 2：绘制 feature 数量柱状图

使用 `feature_count.tsv` 绘制柱状图。

要求：

1. x 轴为 feature 类型；
2. y 轴为数量；
3. 添加标题；
4. x 轴文字旋转 45 度；
5. 保存为 PDF 和 PNG。

---

### 练习 3：绘制每条染色体基因数量柱状图

使用 `chr_gene_count.tsv` 绘制柱状图。

思考：

1. 哪条染色体基因最多？
2. 哪条染色体基因最少？
3. 是否所有染色体基因数量相近？

---

### 练习 4：绘制基因长度分布图

使用 `gene_length.tsv` 绘制：

1. 普通基因长度直方图；
2. log10 转换后的基因长度直方图；
3. log10 转换后的基因长度密度图。

思考：

1. 为什么需要 log10 转换？
2. 基因长度分布是否均匀？
3. 大部分基因集中在哪个长度范围？

---

### 练习 5：绘制不同染色体的基因长度箱线图

使用 `gene_length.tsv` 绘制箱线图。

要求：

1. x 轴为染色体；
2. y 轴为 `log10(gene length)`；
3. 隐藏极端离群点；
4. 保存为 PDF。

---

### 练习 6：绘制染色体 × feature 热图

使用 `chr_feature_count.tsv` 绘制热图。

要求：

1. 行为染色体；
2. 列为 feature；
3. 数值为 count；
4. 使用 `log10(count + 1)` 转换。

---

### 练习 7：绘制染色体长度与基因数量散点图

使用 `chr_length.tsv` 和 `chr_gene_count.tsv` 绘制散点图。

要求：

1. x 轴为染色体长度，单位为 Mb；
2. y 轴为基因数量；
3. 每个点标注染色体名称；
4. 保存为 PDF 和 PNG。

---

## 十、R 综合练习

### 练习 1：读取数据

使用 R 读取以下文件：

```text
feature_count.tsv
chr_gene_count.tsv
gene_length.tsv
chr_feature_count.tsv
chr_length.tsv
```

要求使用：

```r
library(tidyverse)
read_tsv()
```

---

### 练习 2：使用 ggplot2 绘制柱状图

绘制：

1. feature 类型数量柱状图；
2. 每条染色体基因数量柱状图。

要求：

1. 使用 `geom_col()`；
2. 设置 x 轴标签倾斜；
3. 使用 `theme_bw()`；
4. 保存为 PDF 和 PNG。

---

### 练习 3：使用 ggplot2 绘制基因长度分布

绘制：

1. 基因长度直方图；
2. log10 基因长度直方图；
3. log10 基因长度密度图。

要求使用：

```r
geom_histogram()
geom_density()
```

---

### 练习 4：使用 ggplot2 绘制箱线图

使用 `gene_length.tsv` 绘制不同染色体上的基因长度箱线图。

要求：

1. 使用 `geom_boxplot()`；
2. x 轴为染色体；
3. y 轴为 `log10(length)`；
4. 不显示极端离群点。

---

### 练习 5：使用 ggplot2 绘制饼图

使用 `feature_count.tsv` 绘制不同 feature 类型比例图。

要求：

1. 使用 `geom_col()`；
2. 使用 `coord_polar(theta = "y")` 转换为饼图；
3. 使用 `fill = feature` 显示不同类别。

---

### 练习 6：使用 pheatmap 绘制热图

使用 `chr_feature_count.tsv` 生成染色体和 feature 的矩阵，然后绘制热图。

要求：

1. 使用 `pivot_wider()` 转换矩阵；
2. 使用 `log10(count + 1)` 转换；
3. 使用 `pheatmap()` 绘图。

---

### 练习 7：绘制染色体长度与基因数量散点图

使用 `chr_length.tsv` 和 `chr_gene_count.tsv` 绘制散点图。

要求：

1. 使用 `inner_join()` 合并两个表格；
2. 将染色体长度转换成 Mb；
3. 使用 `geom_point()` 绘制散点；
4. 使用 `geom_text()` 标注染色体名称。

---

## 十一、课后思考题

### 1. 图形选择题

请判断下面的数据适合用什么图表示。

| 数据类型 | 推荐图形 |
|---|---|
| 每条染色体基因数量 |  |
| 不同 feature 类型的比例 |  |
| 所有基因长度分布 |  |
| 不同染色体基因长度比较 |  |
| 染色体和 feature 的数量矩阵 |  |
| 染色体长度和基因数量的关系 |  |

---

### 2. 概念题

1. 柱状图和直方图有什么区别？
2. 箱线图中的中位数表示什么？
3. 箱线图中的上下四分位数表示什么？
4. 为什么基因长度分布常常需要 log 转换？
5. 热图适合展示什么类型的数据？
6. 饼图在类别很多时为什么不适合使用？
7. PDF 和 PNG 图像格式有什么区别？
8. 为什么论文作图通常推荐保存 PDF？
9. Python 和 R 都能画图，它们各自有什么优点？
10. 为什么绘图前通常要先整理成规范的表格？

---

### 3. 结果解释题

根据你画出的图回答：

1. 大豆不同染色体上的基因数量是否相同？
2. 哪条染色体基因数量最多？
3. 哪条染色体基因数量最少？
4. 大部分基因长度集中在哪个范围？
5. 染色体长度和基因数量之间是否有明显正相关？
6. GFF3 文件中哪一种 feature 数量最多？
7. exon 数量为什么通常多于 gene 数量？
8. CDS 数量和 exon 数量是否相同？为什么？
9. 如果某条染色体很长但基因数量不多，可能说明什么？
10. 为什么热图中常使用 `log10(count + 1)` 转换？

---

## 十三、学生最终需要提交的内容


```text
02_visualization/
├── feature_count.tsv
├── chr_gene_count.tsv
├── gene_length.tsv
├── chr_feature_count.tsv
├── chr_length.tsv
├── plot_python.py
├── plot_R.R
├── python_feature_barplot.pdf
├── python_chr_gene_barplot.pdf
├── python_gene_length_histogram.pdf
├── python_gene_length_log_histogram.pdf
├── python_gene_length_density.pdf
├── python_gene_length_boxplot.pdf
├── python_feature_piechart.pdf
├── python_chr_feature_heatmap.pdf
├── python_chr_length_gene_scatter.pdf
├── R_feature_barplot.pdf
├── R_chr_gene_barplot.pdf
├── R_gene_length_histogram.pdf
├── R_gene_length_log_histogram.pdf
├── R_gene_length_density.pdf
├── R_gene_length_boxplot.pdf
├── R_feature_piechart.pdf
├── R_chr_feature_heatmap.pdf
└── R_chr_length_gene_scatter.pdf
```

同时提交一份简单报告，回答：

```text
1. 本次使用了哪些输入文件？
2. 这些输入文件分别来自 FASTA 还是 GFF3？
3. 绘制了哪些图？
4. 每种图展示了什么信息？
5. 你认为哪一种图最适合展示基因组统计信息？
6. Python 和 R 作图的区别是什么？
7. 本次绘图过程中遇到了什么问题？如何解决？
```

---

## 十四、本章总结

本节课从上一节课得到的基因组统计结果出发，学习了如何使用 Python 和 R 进行基础可视化。

核心流程是：

```text
FASTA / GFF3 文件
        ↓
Linux 命令统计
        ↓
生成 TSV 表格
        ↓
Python / R 读取表格
        ↓
绘制图形
        ↓
解释生物学意义
```

对于生物信息学初学者来说，绘图不仅是“让结果好看”，更重要的是帮助我们发现数据中的规律。例如：

- 不同染色体上的基因数量是否均匀；
- 基因长度是否集中在某个范围；
- 染色体长度是否影响基因数量；
- GFF3 注释中不同 feature 的组成比例如何；
- 不同 feature 在不同染色体上的分布是否一致。

这一节课建立的是从“数据统计”到“结果展示”的基本能力，是后续学习 RNA-seq 表达热图、差异表达火山图、GO 富集气泡图、基因组共线性图、群体遗传结构图等高级生物信息学图形的基础。

---

## 十五、附录：一个完整的 Python 绘图脚本

如果希望一次性生成所有 Python 图，可以新建一个完整脚本：

```bash
nano plot_all_python.py
```

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

feature = pd.read_csv("feature_count.tsv", sep="\t")
chr_gene = pd.read_csv("chr_gene_count.tsv", sep="\t")
gene_len = pd.read_csv("gene_length.tsv", sep="\t")
chr_feature = pd.read_csv("chr_feature_count.tsv", sep="\t")
chr_length = pd.read_csv("chr_length.tsv", sep="\t")

# 1. feature 柱状图
plt.figure(figsize=(6, 4))
plt.bar(feature["feature"], feature["count"])
plt.xlabel("Feature type")
plt.ylabel("Count")
plt.title("Number of different feature types in GFF3")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_feature_barplot.pdf")
plt.savefig("python_feature_barplot.png", dpi=300)
plt.close()

# 2. 染色体基因数量柱状图
plt.figure(figsize=(8, 4))
plt.bar(chr_gene["chr"], chr_gene["gene_count"])
plt.xlabel("Chromosome")
plt.ylabel("Gene number")
plt.title("Gene number on each chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_chr_gene_barplot.pdf")
plt.savefig("python_chr_gene_barplot.png", dpi=300)
plt.close()

# 3. 基因长度直方图
plt.figure(figsize=(6, 4))
plt.hist(gene_len["length"], bins=50)
plt.xlabel("Gene length")
plt.ylabel("Gene count")
plt.title("Distribution of gene length")
plt.tight_layout()
plt.savefig("python_gene_length_histogram.pdf")
plt.savefig("python_gene_length_histogram.png", dpi=300)
plt.close()

# 4. log10 基因长度直方图
gene_len["log10_length"] = np.log10(gene_len["length"])
plt.figure(figsize=(6, 4))
plt.hist(gene_len["log10_length"], bins=50)
plt.xlabel("log10(Gene length)")
plt.ylabel("Gene count")
plt.title("Distribution of log10 gene length")
plt.tight_layout()
plt.savefig("python_gene_length_log_histogram.pdf")
plt.savefig("python_gene_length_log_histogram.png", dpi=300)
plt.close()

# 5. 密度图
x = gene_len["log10_length"]
density = gaussian_kde(x)
x_grid = np.linspace(x.min(), x.max(), 200)
plt.figure(figsize=(6, 4))
plt.plot(x_grid, density(x_grid))
plt.xlabel("log10(Gene length)")
plt.ylabel("Density")
plt.title("Density plot of gene length")
plt.tight_layout()
plt.savefig("python_gene_length_density.pdf")
plt.savefig("python_gene_length_density.png", dpi=300)
plt.close()

# 6. 箱线图
chr_list = sorted(gene_len["chr"].unique())
data = [gene_len.loc[gene_len["chr"] == chr_name, "log10_length"] for chr_name in chr_list]
plt.figure(figsize=(10, 4))
plt.boxplot(data, labels=chr_list, showfliers=False)
plt.xlabel("Chromosome")
plt.ylabel("log10(Gene length)")
plt.title("Gene length distribution on each chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("python_gene_length_boxplot.pdf")
plt.savefig("python_gene_length_boxplot.png", dpi=300)
plt.close()

# 7. 饼图
plt.figure(figsize=(6, 6))
plt.pie(feature["count"], labels=feature["feature"], autopct="%1.1f%%", startangle=90)
plt.title("Proportion of feature types in GFF3")
plt.tight_layout()
plt.savefig("python_feature_piechart.pdf")
plt.savefig("python_feature_piechart.png", dpi=300)
plt.close()

# 8. 热图
matrix = chr_feature.pivot(index="chr", columns="feature", values="count").fillna(0)
matrix_log = np.log10(matrix + 1)
plt.figure(figsize=(8, 6))
plt.imshow(matrix_log, aspect="auto")
plt.colorbar(label="log10(count + 1)")
plt.xticks(ticks=np.arange(len(matrix_log.columns)), labels=matrix_log.columns, rotation=45, ha="right")
plt.yticks(ticks=np.arange(len(matrix_log.index)), labels=matrix_log.index)
plt.xlabel("Feature type")
plt.ylabel("Chromosome")
plt.title("Chromosome-feature count heatmap")
plt.tight_layout()
plt.savefig("python_chr_feature_heatmap.pdf")
plt.savefig("python_chr_feature_heatmap.png", dpi=300)
plt.close()

# 9. 散点图
df = pd.merge(chr_length, chr_gene, on="chr", how="inner")
df["length_Mb"] = df["length"] / 1_000_000
plt.figure(figsize=(6, 4))
plt.scatter(df["length_Mb"], df["gene_count"])
for _, row in df.iterrows():
    plt.text(row["length_Mb"], row["gene_count"], row["chr"], fontsize=8)
plt.xlabel("Chromosome length (Mb)")
plt.ylabel("Gene number")
plt.title("Relationship between chromosome length and gene number")
plt.tight_layout()
plt.savefig("python_chr_length_gene_scatter.pdf")
plt.savefig("python_chr_length_gene_scatter.png", dpi=300)
plt.close()
```

运行：

```bash
python plot_all_python.py
```

---

## 十六、附录：一个完整的 R 绘图脚本

如果希望一次性生成所有 R 图，可以新建一个完整脚本：

```bash
nano plot_all_R.R
```

```r
library(tidyverse)
library(pheatmap)

feature <- read_tsv("feature_count.tsv")
chr_gene <- read_tsv("chr_gene_count.tsv")
gene_len <- read_tsv("gene_length.tsv")
chr_feature <- read_tsv("chr_feature_count.tsv")
chr_length <- read_tsv("chr_length.tsv")

# 1. feature 柱状图
p <- ggplot(feature, aes(x = feature, y = count)) +
  geom_col() +
  theme_bw() +
  labs(x = "Feature type", y = "Count", title = "Number of different feature types in GFF3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("R_feature_barplot.pdf", p, width = 6, height = 4)
ggsave("R_feature_barplot.png", p, width = 6, height = 4, dpi = 300)

# 2. 染色体基因数量柱状图
p <- ggplot(chr_gene, aes(x = chr, y = gene_count)) +
  geom_col() +
  theme_bw() +
  labs(x = "Chromosome", y = "Gene number", title = "Gene number on each chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("R_chr_gene_barplot.pdf", p, width = 8, height = 4)
ggsave("R_chr_gene_barplot.png", p, width = 8, height = 4, dpi = 300)

# 3. 基因长度直方图
p <- ggplot(gene_len, aes(x = length)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(x = "Gene length", y = "Gene count", title = "Distribution of gene length")
ggsave("R_gene_length_histogram.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_histogram.png", p, width = 6, height = 4, dpi = 300)

# 4. log10 基因长度直方图
gene_len <- gene_len %>% mutate(log10_length = log10(length))
p <- ggplot(gene_len, aes(x = log10_length)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(x = "log10(Gene length)", y = "Gene count", title = "Distribution of log10 gene length")
ggsave("R_gene_length_log_histogram.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_log_histogram.png", p, width = 6, height = 4, dpi = 300)

# 5. 密度图
p <- ggplot(gene_len, aes(x = log10_length)) +
  geom_density() +
  theme_bw() +
  labs(x = "log10(Gene length)", y = "Density", title = "Density plot of gene length")
ggsave("R_gene_length_density.pdf", p, width = 6, height = 4)
ggsave("R_gene_length_density.png", p, width = 6, height = 4, dpi = 300)

# 6. 箱线图
p <- ggplot(gene_len, aes(x = chr, y = log10_length)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(x = "Chromosome", y = "log10(Gene length)", title = "Gene length distribution on each chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("R_gene_length_boxplot.pdf", p, width = 10, height = 4)
ggsave("R_gene_length_boxplot.png", p, width = 10, height = 4, dpi = 300)

# 7. 饼图
p <- ggplot(feature, aes(x = "", y = count, fill = feature)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Proportion of feature types in GFF3", fill = "Feature")
ggsave("R_feature_piechart.pdf", p, width = 6, height = 6)
ggsave("R_feature_piechart.png", p, width = 6, height = 6, dpi = 300)

# 8. 热图
matrix_data <- chr_feature %>%
  pivot_wider(names_from = feature, values_from = count, values_fill = 0)
mat <- as.data.frame(matrix_data)
rownames(mat) <- mat$chr
mat$chr <- NULL
mat_log <- log10(as.matrix(mat) + 1)

pdf("R_chr_feature_heatmap.pdf", width = 8, height = 6)
pheatmap(mat_log, cluster_rows = FALSE, cluster_cols = FALSE, main = "Chromosome-feature count heatmap")
dev.off()

png("R_chr_feature_heatmap.png", width = 2400, height = 1800, res = 300)
pheatmap(mat_log, cluster_rows = FALSE, cluster_cols = FALSE, main = "Chromosome-feature count heatmap")
dev.off()

# 9. 散点图
df <- inner_join(chr_length, chr_gene, by = "chr") %>%
  mutate(length_Mb = length / 1000000)

p <- ggplot(df, aes(x = length_Mb, y = gene_count, label = chr)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = "Chromosome length (Mb)", y = "Gene number", title = "Relationship between chromosome length and gene number")
ggsave("R_chr_length_gene_scatter.pdf", p, width = 6, height = 4)
ggsave("R_chr_length_gene_scatter.png", p, width = 6, height = 4, dpi = 300)
```

运行：

```bash
Rscript plot_all_R.R
```

---

## 十七、进阶：使用 circlize 复现百脉根基因组圈图

### 1. 本节学习目标

前面的内容主要使用大豆参考基因组统计结果学习基础可视化，包括柱状图、箱线图、热图、饼图、密度图和散点图。作为进阶内容，本节不再使用简化示例，而是使用一个真实的百脉根基因组圈图案例，让学生复现一张多轨道基因组圈图。

本节示例数据来自两个百脉根基因组：

```text
GifuT2T
MG20T2T
```

本节使用 R 语言中的 `circlize` 包绘制圈图，展示内容包括：

1. 两个基因组的 12 条染色体；
2. 着丝粒区域；
3. GC 含量；
4. SNP 密度；
5. InDel 密度；
6. 基因密度；
7. TE 密度；
8. GifuT2T 与 MG20T2T 之间的共线性连线。

完成本节学习后，学生应能够理解：

- 什么类型的数据适合用圈图展示；
- 圈图中的 sector、track 和 link 分别代表什么；
- 如何准备基因组圈图所需的输入文件；
- 如何运行一个完整的 `circlize` 绘图脚本；
- 如何根据输出图解释不同轨道展示的信息。

---

### 2. 圈图适合展示什么信息

在基因组学中，很多数据都具有明确的染色体坐标，例如：

- 染色体长度；
- GC 含量；
- 基因密度；
- TE 或重复序列密度；
- SNP / InDel / SV 密度；
- 着丝粒、端粒等特殊区域；
- 不同基因组之间的共线性关系。

这类数据如果都放在普通二维图中，往往需要很多张图才能展示完整。圈图的优势是可以把多种信息放在同一个基因组坐标系统中，让读者从整体上观察基因组结构、变异分布和共线性关系。

| 圈图元素 | 含义 | 本节案例 |
|---|---|---|
| sector | 圆环上的一个扇区，通常代表一条染色体 | G01–G06、M06–M01 |
| track | 圆环上的一层数据轨道 | GC、SNP、InDel、gene、TE 等 |
| link | 圆心区域的连线 | GifuT2T 与 MG20T2T 的共线性关系 |

---

### 3. 本节需要提供给学生的配套文件

本节圈图不需要自己从 GFF 或 FASTA 中生成数据，而是提前提供了一个完整的数据文件夹。文件夹命名为：

```text
lotus_circos_data
```


```text
lotus_circos_data/
    ├── Chr_length.txt
    ├── GC.0.5M.bed
    ├── gene_density.txt
    ├── TE_density.txt
    ├── Cent.txt
    ├── INS_density.txt
    ├── DEL_density.txt
    ├── SNP_density.txt
    └── GifuT2T_MG20T2T_collinear.txt
```

进入数据目录后运行脚本：

```bash
cd ~/bioinfo_course/02_visualization/lotus_circos_data
Rscript ../plot_lotus_circlize.R
```

也可以把脚本复制到数据目录中运行：

```bash
cd ~/bioinfo_course/02_visualization/lotus_circos_data
cp ../plot_lotus_circlize.R .
Rscript plot_lotus_circlize.R
```

运行结束后会生成：

```text
circos_plot_custom_v3_centromere.pdf
```

---

### 4. 每个输入文件的作用和格式

#### 4.1 `Chr_length.txt`

作用：定义染色体名称和坐标范围，是整个圈图的坐标基础。

脚本中读取方式：

```r
chr_length <- read_tsv("Chr_length.txt")
```

因此该文件需要有表头，至少包含以下三列：

```text
Name    chromStart    chromEnd
G01     0             67100000
G02     0             52000000
G03     0             48000000
```

字段说明：

| 列名 | 含义 |
|---|---|
| `Name` | 染色体名称，必须与其他文件第一列保持一致 |
| `chromStart` | 染色体起始坐标，通常为 0 |
| `chromEnd` | 染色体终止坐标，也就是染色体长度 |

本脚本设定的染色体顺序为：

```r
chr_order <- c(
  "G01","G02","G03","G04","G05","G06",
  "M06","M05","M04","M03","M02","M01"
)
```

其中 `G01–G06` 表示 GifuT2T 的 6 条染色体，`M01–M06` 表示 MG20T2T 的 6 条染色体。这里将 MG20T2T 按 `M06` 到 `M01` 的顺序排列，是为了让两个基因组在圈图中形成更清楚的对应关系。

---

#### 4.2 `GC.0.5M.bed`

作用：提供每个窗口的 GC 含量，用折线轨道展示。

脚本中读取方式：

```r
GC_density <- read_tsv("GC.0.5M.bed", col_names = FALSE)
```

该文件没有表头，包含 4 列：

```text
G01     0        500000    0.36
G01     500000   1000000   0.38
G01     1000000  1500000   0.35
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 窗口起始位置 |
| 第 3 列 | `X3` | 窗口终止位置 |
| 第 4 列 | `X4` | GC 含量 |

文件名中的 `0.5M` 表示这个文件通常是按照 0.5 Mb 窗口统计的。

---

#### 4.3 `gene_density.txt`

作用：提供基因密度，用热图轨道展示。

脚本中读取方式：

```r
gene_density <- read_tsv("gene_density.txt", col_names = FALSE)
```

该文件没有表头，包含 4 列：

```text
G01     0        500000    12
G01     500000   1000000   18
G01     1000000  1500000   9
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 窗口起始位置 |
| 第 3 列 | `X3` | 窗口终止位置 |
| 第 4 列 | `X4` | 该窗口中的基因数量或基因密度 |

脚本会将第 4 列归一化到 0–1：

```r
gene_density <- gene_density %>% mutate(gene_norm = normalize01(X4))
```

因此学生不需要提前手动归一化。

---

#### 4.4 `TE_density.txt`

作用：提供 TE 或重复序列密度，用热图轨道展示。

脚本中读取方式：

```r
TE_density <- read_tsv("TE_density.txt", col_names = FALSE)
```

该文件没有表头，格式与 `gene_density.txt` 一致：

```text
G01     0        500000    0.42
G01     500000   1000000   0.55
G01     1000000  1500000   0.48
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 窗口起始位置 |
| 第 3 列 | `X3` | 窗口终止位置 |
| 第 4 列 | `X4` | 该窗口中的 TE 密度或 TE 覆盖比例 |

脚本同样会将第 4 列归一化到 0–1：

```r
TE_density <- TE_density %>% mutate(te_norm = normalize01(X4))
```

---

#### 4.5 `Cent.txt`

作用：提供着丝粒区域，用最外层红色矩形块展示。

脚本中读取方式：

```r
cent <- read_tsv("Cent.txt", col_names = FALSE)
```

该文件没有表头，包含 4 列：

```text
G01     25000000    31000000    1
G02     22000000    28000000    1
M01     24000000    30000000    1
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 着丝粒起始位置 |
| 第 3 列 | `X3` | 着丝粒终止位置 |
| 第 4 列 | `X4` | 标记值，通常填 1 |

脚本中着丝粒颜色由下面参数控制：

```r
CENT_COL <- "#B30000"
```

---

#### 4.6 `SNP_density.txt`

作用：提供 SNP 密度，用柱状轨道展示。

脚本中读取方式：

```r
snp_density <- read_tsv("SNP_density.txt", col_names = FALSE)
```

该文件没有表头，包含 4 列：

```text
G01     0        500000    0.00024
G01     500000   1000000   0.00031
G01     1000000  1500000   0.00028
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 窗口起始位置 |
| 第 3 列 | `X3` | 窗口终止位置 |
| 第 4 列 | `X4` | SNP 密度 |

脚本会先用窗口长度把密度转换为数量：

```r
snp_count = X4 * win_len
```

然后对极端高值进行裁剪，并归一化到 0–1：

```r
snp_clip = clip_upper(snp_count, q = CLIP_Q)
snp_norm = normalize01(snp_clip)
```

这里的 `CLIP_Q` 是极端值裁剪分位数，默认值为：

```r
CLIP_Q <- 0.99
```

---

#### 4.7 `INS_density.txt` 和 `DEL_density.txt`

作用：分别提供插入和缺失变异密度。脚本会将两者合并为 InDel 密度轨道。

脚本中读取方式：

```r
ins_density <- read_tsv("INS_density.txt", col_names = FALSE)
del_density <- read_tsv("DEL_density.txt", col_names = FALSE)
```

这两个文件都没有表头，格式相同：

```text
G01     0        500000    0.00005
G01     500000   1000000   0.00008
G01     1000000  1500000   0.00004
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | 染色体名称 |
| 第 2 列 | `X2` | 窗口起始位置 |
| 第 3 列 | `X3` | 窗口终止位置 |
| 第 4 列 | `X4` | INS 或 DEL 密度 |

脚本会将 INS 和 DEL 转换为数量后相加：

```r
indel_count = ins_count + del_count
```

再进行极端值裁剪和 0–1 归一化。

---

#### 4.8 `GifuT2T_MG20T2T_collinear.txt`

作用：提供 GifuT2T 和 MG20T2T 之间的共线性区段，用圈图中心的连线展示。

脚本中读取方式：

```r
collinearity <- read_tsv("GifuT2T_MG20T2T_collinear.txt", col_names = FALSE)
```

该文件没有表头，至少包含 6 列：

```text
G01     1000000    1500000    M01     1200000    1700000
G02     5000000    5600000    M02     5100000    5700000
G03     8000000    8500000    M03     7800000    8300000
```

字段说明：

| 列号 | 脚本列名 | 含义 |
|---|---|---|
| 第 1 列 | `X1` | GifuT2T 染色体 |
| 第 2 列 | `X2` | GifuT2T 区段起始位置 |
| 第 3 列 | `X3` | GifuT2T 区段终止位置 |
| 第 4 列 | `X4` | MG20T2T 染色体 |
| 第 5 列 | `X5` | MG20T2T 区段起始位置 |
| 第 6 列 | `X6` | MG20T2T 区段终止位置 |

脚本中连线颜色根据 GifuT2T 染色体编号设置：

```r
g_chr <- substr(collinearity$X1[i], 2, 3)
col = link_cols[g_chr]
```

因此共线性文件第一列建议使用 `G01`、`G02`、`G03` 这类格式，而不要写成 `Chr01` 或其他格式。

---

### 5. 安装 R 包

前面 R 绘图部分已经介绍过 `tidyverse`。本节圈图还需要额外安装 `circlize`。

进入 R 后安装：

```r
install.packages(c("tidyverse", "pheatmap", "circlize"))
```

如果已经安装过，就不需要重复安装。

脚本中会加载：

```r
library(tidyverse)
library(circlize)
library(dplyr)
```

其中 `dplyr` 实际上已经包含在 `tidyverse` 中，这里保留 `library(dplyr)` 是为了让初学者清楚脚本中使用了 `mutate()`、`select()`、`arrange()`、`full_join()` 等数据整理函数。

---

### 6. 运行方式

本课程统一采用“编辑脚本文件，然后用命令行运行脚本”的方式。

#### 6.1 新建 R 脚本

```bash
cd ~/bioinfo_course/02_visualization
nano plot_lotus_circlize.R
```

将下面完整代码复制进去，保存并退出。

#### 6.2 准备数据目录

```

确认目录中包含以下文件：

```bash
ls lotus_circos_data
```

应能看到：

```text
Chr_length.txt
GC.0.5M.bed
gene_density.txt
TE_density.txt
Cent.txt
INS_density.txt
DEL_density.txt
SNP_density.txt
GifuT2T_MG20T2T_collinear.txt
```

#### 6.3 运行脚本

因为脚本中直接读取当前目录下的文件，所以要进入数据目录运行：

```bash
cd ~/bioinfo_course/02_visualization/lotus_circos_data
Rscript ../plot_lotus_circlize.R
```

运行成功后，会看到类似提示：

```text
Done: circos_plot_custom_v3_centromere.pdf  (CLIP_Q=0.99)
```

输出文件为：

```text
circos_plot_custom_v3_centromere.pdf
```

---

### 7. 完整 R 脚本：`plot_lotus_circlize.R`

```r
# ============================================================
# Circos plot (customized v3) - centromere version
# SNP & INDEL: NO log2, normalize (0-1) with outlier clipping
# Outer track: centromere blocks (from cent.txt), NOT Old_region
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(dplyr)
})

# ------------------------------
# 0. Params
# ------------------------------
CLIP_Q <- 0.99   # 0.99 / 0.995 / 0.999
CENT_COL <- "#B30000"  # 着丝粒块颜色，可改

# ------------------------------
# 1. Read inputs
# ------------------------------
chr_length   <- read_tsv("Chr_length.txt")

GC_density   <- read_tsv("GC.0.5M.bed", col_names = FALSE)
gene_density <- read_tsv("gene_density.txt", col_names = FALSE)
TE_density   <- read_tsv("TE_density.txt", col_names = FALSE)

# 用 cent.txt 替代 Old.bed
cent <- read_tsv("Cent.txt", col_names = FALSE)   # X1 chr, X2 start, X3 end, X4=1

ins_density  <- read_tsv("INS_density.txt", col_names = FALSE)
del_density  <- read_tsv("DEL_density.txt", col_names = FALSE)
snp_density  <- read_tsv("SNP_density.txt", col_names = FALSE)

collinearity <- read_tsv("GifuT2T_MG20T2T_collinear.txt", col_names = FALSE)

# ------------------------------
# 2. Chromosome order
# ------------------------------
chr_order <- c(
  "G01","G02","G03","G04","G05","G06",
  "M06","M05","M04","M03","M02","M01"
)

chr_length$Name <- factor(chr_length$Name, levels = chr_order)
chr_length <- chr_length %>% arrange(Name)

# ------------------------------
# 3. Helpers
# ------------------------------
to_numeric_bed <- function(df) {
  df %>%
    mutate(
      X1 = as.character(X1),
      X2 = as.numeric(X2),
      X3 = as.numeric(X3),
      X4 = as.numeric(X4)
    )
}

normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (is.infinite(rng[1]) || is.infinite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

# 只裁剪极端高值（winsorize upper tail）
clip_upper <- function(x, q = 0.995) {
  cap <- as.numeric(quantile(x, probs = q, na.rm = TRUE, names = FALSE, type = 7))
  pmin(x, cap)
}

# numeric convert
GC_density   <- to_numeric_bed(GC_density)
gene_density <- to_numeric_bed(gene_density)
TE_density   <- to_numeric_bed(TE_density)

cent        <- to_numeric_bed(cent)
ins_density <- to_numeric_bed(ins_density)
del_density <- to_numeric_bed(del_density)
snp_density <- to_numeric_bed(snp_density)

# factor + order
GC_density$X1   <- factor(GC_density$X1,   levels = chr_order); GC_density   <- GC_density   %>% arrange(X1)
gene_density$X1 <- factor(gene_density$X1, levels = chr_order); gene_density <- gene_density %>% arrange(X1)
TE_density$X1   <- factor(TE_density$X1,   levels = chr_order); TE_density   <- TE_density   %>% arrange(X1)

cent$X1        <- factor(cent$X1,        levels = chr_order); cent        <- cent        %>% arrange(X1)
ins_density$X1 <- factor(ins_density$X1, levels = chr_order); ins_density <- ins_density %>% arrange(X1)
del_density$X1 <- factor(del_density$X1, levels = chr_order); del_density <- del_density %>% arrange(X1)
snp_density$X1 <- factor(snp_density$X1, levels = chr_order); snp_density <- snp_density %>% arrange(X1)

# ------------------------------
# 4. SNP: count -> clip -> normalize (NO log2)
# ------------------------------
snp_density <- snp_density %>%
  mutate(
    win_len   = (X3 - X2),
    snp_count = X4 * win_len
  ) %>%
  mutate(
    snp_clip = clip_upper(snp_count, q = CLIP_Q),
    snp_norm = normalize01(snp_clip)
  )

# ------------------------------
# 5. INDEL: (INS+DEL) -> clip -> normalize (NO log2)
# ------------------------------
ins2 <- ins_density %>%
  mutate(win_len = (X3 - X2),
         ins_count = X4 * win_len) %>%
  select(X1, X2, X3, ins_count)

del2 <- del_density %>%
  mutate(win_len = (X3 - X2),
         del_count = X4 * win_len) %>%
  select(X1, X2, X3, del_count)

indel_density <- full_join(ins2, del2, by = c("X1","X2","X3")) %>%
  mutate(
    ins_count   = replace_na(ins_count, 0),
    del_count   = replace_na(del_count, 0),
    indel_count = ins_count + del_count
  ) %>%
  mutate(
    indel_clip = clip_upper(indel_count, q = CLIP_Q),
    indel_norm = normalize01(indel_clip)   # <<< 修正：用裁剪后的值归一化
  )

# ------------------------------
# 6. Gene / TE normalize (keep heatmap)
# ------------------------------
gene_density <- gene_density %>% mutate(gene_norm = normalize01(X4))
TE_density   <- TE_density   %>% mutate(te_norm   = normalize01(X4))

# ------------------------------
# 7. Colors
# ------------------------------
col_gc     <- "#2A9D8F"
col_snp    <- "#E76F51"
col_indel  <- "#457B9D"
col_gene   <- colorRamp2(c(0, 1), c("white", "#1B9E77"))
col_te     <- colorRamp2(c(0, 1), c("white", "#6A3D9A"))

link_cols <- c(
  "#4E79A7",
  "#59A14F",
  "#76B7B2",
  "#9C755F",
  "#B07AA1",
  "#F28E2B"
)
names(link_cols) <- paste0("0", 1:6)

# ------------------------------
# 8. Plot
# ------------------------------
pdf("circos_plot_custom_v3_centromere.pdf", width = 10, height = 10)

circos.clear()

# gap: group separation larger
gap_degree <- c(rep(4, 5), 30, rep(4, 5), 30)

circos.par(
  start.degree = 82,
  clock.wise = TRUE,
  gap.degree = gap_degree,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(
  factors = chr_order,
  xlim = cbind(chr_length$chromStart, chr_length$chromEnd)
)

# ---- Track 1: Centromere blocks ----
cent_gr <- cent %>% select(chr = X1, start = X2, end = X3, value = X4)

circos.genomicTrack(
  cent_gr, ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    # 用整条track高度的块标记着丝粒
    circos.genomicRect(
      region, value,
      ybottom = 0, ytop = 1,
      col = CENT_COL,
      border = NA, ...
    )
  },
  track.height = 0.085, bg.border = NA
)

# chr labels + axis on the outer track
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$cell.ylim[2] + mm_y(6),
              CELL_META$sector.index,
              facing = "outside", niceFacing = TRUE, cex = 1)
  circos.axis(h = "top", labels.cex = 0.8,
              major.tick.length = mm_y(0.8))
}, bg.border = NA)

# ---- Track 2: GC line (thin) ----
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    df = GC_density %>% filter(X1 == CELL_META$sector.index)
    if (nrow(df) > 0)
      circos.lines(df$X2, df$X4, col = col_gc, lwd = 1)
  }
)

# ---- Track 3: SNP bars (clipped + normalized) ----
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.10,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    df = snp_density %>% filter(X1 == CELL_META$sector.index)
    if (nrow(df) > 0)
      circos.rect(df$X2, 0, df$X3, df$snp_norm,
                  col = col_snp, border = NA)
  }
)

# ---- Track 4: INDEL bars (clipped + normalized) ----
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.10,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    df = indel_density %>% filter(X1 == CELL_META$sector.index)
    if (nrow(df) > 0)
      circos.rect(df$X2, 0, df$X3, df$indel_norm,
                  col = col_indel, border = NA)
  }
)

# ---- Track 5: Gene density heatmap ----
gene_gr <- gene_density %>% select(chr = X1, start = X2, end = X3, value = gene_norm)

circos.genomicTrack(
  gene_gr, ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
                       col = col_gene(value[[1]]),
                       border = NA, ...)
  },
  track.height = 0.08, bg.border = NA
)

# ---- Track 6: TE density heatmap ----
te_gr <- TE_density %>% select(chr = X1, start = X2, end = X3, value = te_norm)

circos.genomicTrack(
  te_gr, ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
                       col = col_te(value[[1]]),
                       border = NA, ...)
  },
  track.height = 0.08, bg.border = NA
)

# ---- Links ----
for (i in 1:nrow(collinearity)) {
  g_chr <- substr(collinearity$X1[i], 2, 3)
  circos.link(
    collinearity$X1[i], c(collinearity$X2[i], collinearity$X3[i]),
    collinearity$X4[i], c(collinearity$X5[i], collinearity$X6[i]),
    col = link_cols[g_chr], border = NA
  )
}

dev.off()
message("Done: circos_plot_custom_v3_centromere.pdf  (CLIP_Q=", CLIP_Q, ")")

```

---

### 8. 圈图轨道说明

这张图从外到内大致包括以下部分：

| 层级 | 数据来源 | 图形形式 | 生物学含义 |
|---|---|---|---|
| 最外层 | `Chr_length.txt` | 染色体标签和坐标轴 | 展示两个百脉根基因组的染色体组成 |
| Track 1 | `Cent.txt` | 红色矩形块 | 标记着丝粒区域 |
| Track 2 | `GC.0.5M.bed` | 折线 | 展示沿染色体分布的 GC 含量 |
| Track 3 | `SNP_density.txt` | 柱状图 | 展示 SNP 变异密度 |
| Track 4 | `INS_density.txt` + `DEL_density.txt` | 柱状图 | 展示 InDel 变异密度 |
| Track 5 | `gene_density.txt` | 热图 | 展示基因密度 |
| Track 6 | `TE_density.txt` | 热图 | 展示 TE 或重复序列密度 |
| 中心连线 | `GifuT2T_MG20T2T_collinear.txt` | 连线 | 展示 GifuT2T 与 MG20T2T 的共线性关系 |

---

### 9. 代码重点解释

#### 9.1 为什么设置 `chr_order`？

```r
chr_order <- c(
  "G01","G02","G03","G04","G05","G06",
  "M06","M05","M04","M03","M02","M01"
)
```

圈图中染色体的排列顺序不会自动符合我们的展示需求，因此需要手动设置。这个顺序先放 GifuT2T 的 6 条染色体，再反向放 MG20T2T 的 6 条染色体，有利于观察两个基因组之间的共线性连线。

---

#### 9.2 为什么要对 SNP 和 InDel 做裁剪和归一化？

```r
CLIP_Q <- 0.99
```

SNP 和 InDel 在某些区域可能特别集中，如果直接绘图，少数极端高值会让大部分区域看起来都很低，不利于观察整体分布。因此脚本使用上分位数裁剪极端高值，再把结果归一化到 0–1。

相关函数为：

```r
clip_upper <- function(x, q = 0.995) {
  cap <- as.numeric(quantile(x, probs = q, na.rm = TRUE, names = FALSE, type = 7))
  pmin(x, cap)
}

normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (is.infinite(rng[1]) || is.infinite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}
```

---

#### 9.3 为什么要用 `circos.clear()`？

`circlize` 会保存当前图形参数。如果连续绘制多张圈图，不先清空可能会造成参数冲突。因此脚本在绘图前运行：

```r
circos.clear()
```

---

#### 9.4 `circos.initialize()` 的作用

```r
circos.initialize(
  factors = chr_order,
  xlim = cbind(chr_length$chromStart, chr_length$chromEnd)
)
```

这一句用于初始化圈图的染色体坐标系统：

- `factors` 指定有哪些染色体；
- `xlim` 指定每条染色体的起始和终止坐标；
- 后续所有 track 都会基于这个坐标系统绘制。

---

#### 9.5 `circos.link()` 的作用

```r
circos.link(
  collinearity$X1[i], c(collinearity$X2[i], collinearity$X3[i]),
  collinearity$X4[i], c(collinearity$X5[i], collinearity$X6[i]),
  col = link_cols[g_chr], border = NA
)
```

这一部分用于绘制 GifuT2T 与 MG20T2T 之间的共线性连线。每一行共线性文件都会生成一条连接线。

---

### 10. 常见报错和检查方法

#### 10.1 找不到文件

报错示例：

```text
Error: 'Chr_length.txt' does not exist in current working directory
```

原因：没有在数据目录中运行脚本，或者文件名写错。

检查方法：

```bash
pwd
ls
```

确认当前目录中有：

```text
Chr_length.txt
GC.0.5M.bed
gene_density.txt
TE_density.txt
Cent.txt
INS_density.txt
DEL_density.txt
SNP_density.txt
GifuT2T_MG20T2T_collinear.txt
```

---

#### 10.2 染色体名称不一致

如果某个数据文件中写的是 `G1`，而 `chr_order` 中写的是 `G01`，就可能导致该染色体的数据无法正确显示。

检查第一列染色体名称：

```bash
cut -f 1 GC.0.5M.bed | sort | uniq
cut -f 1 gene_density.txt | sort | uniq
cut -f 1 TE_density.txt | sort | uniq
cut -f 1 SNP_density.txt | sort | uniq
```

应与脚本中的名称一致：

```text
G01 G02 G03 G04 G05 G06 M01 M02 M03 M04 M05 M06
```

---

#### 10.3 表格列数不对

检查文件列数：

```bash
awk '{print NF}' GC.0.5M.bed | sort | uniq -c
awk '{print NF}' GifuT2T_MG20T2T_collinear.txt | sort | uniq -c
```

其中密度文件通常应为 4 列，共线性文件至少应为 6 列。

---

#### 10.4 R 包没有安装

报错示例：

```text
there is no package called 'circlize'
```

解决方法：

```r
install.packages("circlize")
```

如果服务器不能联网，可以让教师提前在教学环境中安装，或使用 conda 环境安装 R 包。

---

### 11. 学生练习任务

#### 练习 1：复现圈图

使用教师提供的百脉根数据和 `plot_lotus_circlize.R` 脚本，运行：

```bash
cd ~/bioinfo_course/02_visualization/lotus_circos_data
Rscript ../plot_lotus_circlize.R
```

提交输出文件：

```text
circos_plot_custom_v3_centromere.pdf
```

---

#### 练习 2：查看输入文件

分别查看每个输入文件的前 5 行：

```bash
head Chr_length.txt
head GC.0.5M.bed
head gene_density.txt
head TE_density.txt
head Cent.txt
head SNP_density.txt
head INS_density.txt
head DEL_density.txt
head GifuT2T_MG20T2T_collinear.txt
```

回答：

1. 哪个文件有表头？
2. 哪些文件是 4 列？
3. 哪个文件用于绘制中心连线？
4. 哪个文件用于标记着丝粒？

---

#### 练习 3：修改参数观察图形变化

将脚本中的：

```r
CLIP_Q <- 0.99
```

分别改为：

```r
CLIP_Q <- 0.995
```

或：

```r
CLIP_Q <- 0.999
```

重新运行脚本，比较 SNP 和 InDel 轨道是否发生变化。

思考：

1. `CLIP_Q` 越大，保留的极端值越多还是越少？
2. 为什么变异密度轨道需要裁剪极端值？
3. 如果不裁剪，图形可能出现什么问题？

---

#### 练习 4：修改着丝粒颜色

将脚本中的：

```r
CENT_COL <- "#B30000"
```

改成其他颜色，例如：

```r
CENT_COL <- "#4D4D4D"
```

重新运行脚本，观察最外层着丝粒轨道颜色变化。

---

#### 练习 5：解释图形

根据生成的圈图回答：

1. 图中一共有多少个 sector？
2. 哪些 sector 属于 GifuT2T？
3. 哪些 sector 属于 MG20T2T？
4. 红色块表示什么区域？
5. 哪一层表示 GC 含量？
6. 哪一层表示 SNP 密度？
7. 哪一层表示 InDel 密度？
8. 哪两层是热图轨道？
9. 中心连线表示什么？
10. GifuT2T 和 MG20T2T 之间是否整体保持较强共线性？

---

### 12. 本节进阶小结

本节使用真实的百脉根 GifuT2T 与 MG20T2T 数据复现了一个多轨道基因组圈图。与前面的大豆基础统计图不同，圈图更适合展示沿染色体分布的多层基因组信息。

本节需要掌握三个核心概念：

```text
sector = 染色体或分组
track  = 圆环上的一层数据轨道
link   = 圆心区域的连线关系
```

同时需要理解，圈图本质上不是单纯为了“好看”，而是为了在统一的基因组坐标系统中整合多类数据。对于比较基因组研究而言，圈图可以同时展示基因组结构、变异分布、重复序列分布和共线性关系，是非常常见的综合性可视化方法。
