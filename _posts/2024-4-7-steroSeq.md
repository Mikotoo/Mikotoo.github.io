---
title: 【文章复现】 2023-Nature plants-大豆根瘤单细胞与空转
author:
date: 2024-04-07
category: Jekyll
layout: post
---

> 单细胞、空转分析流程学习，整理备查，顺便记录 <br>
> 本文作者是南方科技大学翟继先老师，文章相关的数据和[代码][1]已经公开。


## 1、文献概述
#### 1.1 研究目的
<p>大豆的共生固氮器官根瘤作为一种高度异质的组织，其各种类型的细胞都具有不同的生理特点和功能。</p>
<p>本文的研究目的是了解根瘤中不同类型的细胞在根瘤成熟过程中的具体贡献，以及细胞之间的关系。</p>

#### 1.2 建库测序
1）为了揭示根瘤成熟过程中的cell-type-specific动态基因表达，选取3个样本进行了单细胞核测序(10X Genomics Chromium): <br>

1. 接种12天的根瘤 (sn_12dpi)
2. 接种21天的根瘤 (sn_21dpi)
3. 作为对照，接种21天根瘤附近的根 (sn_root)

2）为了更好的进行细胞注释，又对上述几个时期的材料进行了空间转录组测序(stereo-seq) <br>

<i><b>两种技术结合，共同完成了细胞分类 </b></i><br>

![pic2][2]


#### 1.3 聚类与注释

###### 1）单细胞核RNA测序分析
1. 对3样本数据分别比对，统计细胞数、基因数等信息
2. 将3组dataset整合，得到了 15个cell clusters（cluster0-cluster14），并得到了每个cluster中上调的基因

	![pic3][3]

3. 利用拟南芥和其它豆科植物的标记基因对细胞进行注释，注释出了15个cluster中的5个，包括
	1. 根表皮细胞 cluster 5
	2. 根维管束 cluster 3
	3. 根瘤维管束 cluster 9
	4. 根瘤皮层 cluster 1
	5. 中心感染区	cluster 12

###### 2）Stereo-seq
<p>由于大豆根瘤中标记基因的缺乏，很多cluster没有注释成功。</p>
<p>为了克服这个问题，对12dpi，21dpi的根瘤做空转，根据空间表达信息对细胞进行注释，将其分为6个区域。</p>

![pic4][4]

> 感染区，内皮层，外皮层(cluster 2 and 4)，表皮，微管


<p> 基于去卷积(deconvolution)的方法，验证了之前注释的cluster，并对其他未注释的cluster完成了注释 </p>

1. 中心感染区 cluster 0; cluster7; cluster 11
2. 根瘤外皮层 cluster 2; cluster 4


![pic5][5]

<p>
	为了验证注释的准确性，对细胞特异性基因进行了GUS染色和RNA原位杂交，结果符合预测
</p>

![pic6][6]


#### 1.4 感染区亚型


<p><b>有4个cluster被定位到中心感染区，分别为0，7，11，12。</b></p>

<p>在2022年一篇百脉根<i>Lotus japonocus</i>的单细胞工作，手动分离了被细菌感染的细胞和没有被感染的细胞</p>

<p>因此，本文根据百脉根不同细胞的基因表达模式，将cluster 0，7，11确定为未感染的细胞(UC)，cluster 12为根瘤菌感染的细胞(IC)</p>

![pic7][7]

###### 1）UC

<p> 在UC中，cluster 0 在dpi 12，dpi 21两个发育时期都存在，而cluster 7 和 11 基本只存在于dpi 21的根瘤中。</p>

<p>
	为了揭示这几个不同UC cluster的分化轨迹，对三个cluster做拟时序分析，推断出分化方向是从cluster 0到 7和11，说明这两个cluster是在根瘤成熟过程中从cluster 0发展而来的。
</p>

<p>在大豆中，根瘤固定的氮主要用于合成脲类物质，合成过程主要发生在UC中。 
</p>

<p>本文发现，和脲合成有关的尿酸酶和天冬氨酸转移酶基因在3个UC cluster中都表达，特别是在cluster 7中上调。 而脲类物质转运相关的基因主要在cluster 0中表达。
</p>

![pic8][8]

**以上结果揭示了根瘤细胞内脲类化合物生产和运输的区室化**

<p>另外，β-淀粉酶在 cluster 11 中显著上调，说明cluster 11 参与共生固氮的能量供应。
</p>

<p>
总之，这些结果表明，UC可以分为不同的功能特化亚细胞类型，其中两种在根瘤发育的后期出现，这可以促进共生所需的营养和能量来源的交换。
</p>

###### 2）IC

对于根瘤菌侵染的细胞。<br>
首先，豆血红蛋白和nodulin基因在cluster 12中显著上调。<br>

<p>
	另外，编码糖转运蛋白和多个苹果酸合成酶家族的基因也在cluster 12中显著上调，说明大豆根瘤细胞核根瘤菌之间活跃的碳氮交换。
</p>

<p>
	对cluster 12进一步聚类，可以将其分为12-0和12-1两种亚细胞类型，其中12-0包含492个细胞，12-1包含38个细胞。<br>
	12-0在12dpi和21dpi两个发育时期中存在，而12-1几乎完全被12dpi的未成熟根瘤占据。RNA原位杂交可以证实这是两种不同的细胞亚型。
</p>

![pic9][9]

![pic10][10]


<p>编码共生体膜蛋白的基因在12-0中远高于12-1以及其他的cluster，说明溶质在12-0类型细胞的共生体之间更活跃的移动。
</p>

![pic11][11]

<p>
	检查了12-1 cluster中的特异性基因，发现50个基因中的9个是先前报道的SNF基因，这一比例显著高于其他cluster的细胞。<br>
	(6个20年PC综述的，3个近期发现的)
</p>

![pic12][12]

<p>
	进一步发现，这9个基因全都参与了侵染线的形成。说明cluster12-1的细胞可能参与了侵染线延申和根瘤菌释放等功能。
</p>

<p>
	最后，本文探索了cluster12-1中一个特异表达基因GLYMA_02G004800，敲除该基因后，50%植株根瘤数目增加，且根瘤的中心感染区呈白色，说明了cluster12-1的细胞在根瘤发育过程中的重要作用。
</p>

![pic13][13]



## 2、测序原理
### 2.1 单细胞转录组测序
### 2.2 空间转录组测序

## 3、软件安装

###### 0. snakemake

> 用于执行snakefile任务

###### 1. Cell Ranger
> 用于处理单细胞原始数据

官网[下载][14]后直接解压

```
pip install loguru pyranges
```

###### 2. velocyto

> 细胞速率分析
conda create -n velocyto python=3.8
conda activate velocyto
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install velocyto


###### 3. scanpy

> 用于单细胞数据的下游分析

```
conda create -n singleCell python=3.8
conda activate singleCell
conda install -c conda-forge scanpy python-igraph leidenalg
or
pip install scanpy
```

###### 4. Seurat

> R包，用于单细胞数据的下游分析

```
#R
BiocManager::install("Seurat") 
```

###### 5. scDblFinder

> R包，用于去除双峰

```
BiocManager::install("scDblFinder")
``` 

###### 6. scvi-tools; scanorama; harmonypy

> 用于集成单细胞数据

```
pip install scvi-tools annoy==1.16.0 scanorama harmonypy
```

###### 7. SAW
> 用于处理stereo-seq原始数据

使用Apptainer容器下载

```
#shell
conda install apptainer
apptainer build SAW_7.1.sif docker://stomics/saw:07.1.0
```

###### 8. stereopy

> 用于空转数据的下游分析

```
# 在新建的python3.8环境中安装

conda create -n st python=3.8
conda activate st
pip install stereopy IPython
```

###### 9. sctransform

> R包，用于空转数据的标准化

```
install.packages("sctransform")
```

###### 10. muon

> 用于空转数据的标准化

```
pip install muon
```

###### 11. scVelo

> 用于构建细胞轨迹

```
pip install scVelo
```

###### 12. monocle3

> R包，用于构建细胞轨迹

首先，手动安装[rtools][15]

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
options(download.file.method = "wininet")
devtools::install_github('cole-trapnell-lab/monocle3')
```
上述方法本地安装成功，但是服务器安装失败。使用conda安装成功。`conda install -c bioconda -y r-monocle3`

###### 13. cellrank

> 用于构建细胞轨迹

```
conda create -n cellrank python=3.8
conda activate cellrank
pip install cellrank
```

###### 14. diffxpy

> 用于鉴定差异表达基因

先下载[batchglm][16]源码目录，解压打开目录 `pip install -e .`

然后下载[diffxpy][17]源码目录，解压打开后 `pip install -e .`


###### 15. AUCell

> R包，用于计算基因集表达分数

```
BiocManager::install("AUCell")
```

###### 16. pyscenic 

> 用于计算基因集表达分数

```
pip install pyscenic
```

###### 17. rpy2

> 用于在python环境中调用R包

```
conda install -c r rpy2
```

###### 18. OrthoFinder 

###### 19. clusterprofiler

## 4、复现流程

### 4.1 参考基因组

**本文使用的基因组为 Soybean [Wm82 a2.v1][19], [Arabidopsis v11][18]**

<p>准备流程需要的gtf和bed文件：</p>

```
# gff2gtf
gffread Wm82v2.gff -T -o Wm82v2.gtf

# gff2bed
paftools.js gff2bed Wm82v2.gtf -j > Wm82v2_junc.bed
```

> paftools.js是minimap2的一部分，如果报错，尝试从github下载最新版的paftools.js文件,替换bin目录`which paftools.js`中现有的paftools.js。

### 4.2 为cellRanger构建index

```
cellranger mkref --genome=Wm82v2 --nthreads=48 --fasta=Wm82v2_genome.fa --genes=Wm82v2.gtf
```

成功之后的index文件：<br>
![pic20][20]

**注意：构建索引后，必须把索引目录中的genes/genes.gtf.gz解压，否则后续流程无法正常运行**

### 4.3 cellRanger获取表达矩阵

> 项目路径 `/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/`

**[作者已经把程序打包][21]** <br>

用snakemake封装了3个step:

1. cellranger count 比对
2. 自定义[python脚本][25]，从比对结果的bam文件中获取UMI信息
3. velocyto, 从比对结果中获取loom文件

由于数据目录、软件版本、conda环境等的不同，我对snakefile文件进行了[修改][22]

#### 4.3.1 项目目录的组织

<p>
	Snakemake是基于python3的数据分析流程构建工具，可以通过它将一系列任务组织起来，并通过config文件构建一个可重复、可扩展的pipeline<br>
	它可以结合conda，将流程扩展到服务器、集群环境中使用<br>
	还可以根据任务之间的依赖关系，智能的并行可以并行执行的任务
</p>

<p>
	Snakemake将以在提交任务时激活的conda环境为默认环境，如果流程中的不同步骤依赖不同的conda环境，可以在rule中加入conda yaml作为参数，snakemake将据此在.snakemake/conda目录中创建所需的环境。<br>
	本项目的第三步velocyto需要一个单独的conda环境，因此需要导出它的yaml文件备用
</p>

```
conda env export -n velocyto -f velocyto.yaml
```

**程序被编写到snakefile文件中；<br>
config.yaml文件作为索引，存储了工作目录、参考基因组路径、原始数据路径等信息;<br>**

![pic23][23]

每个样本的原始数据分别存储到不同的目录，并以特定的格式命名： **sample_S1_L001_R1(R2)_001.fastq.gz** 

![pic23][24]

**准备完成后，在snakefile所在目录下输入`snakemake`命令，snakemake会自动运行**

#### 4.3.2 snakefile流程注释

##### 1 config.yaml

1. cellRangerPath: cellranger的路径，似乎没有用到
2. cellRangerIndex: cellranger注释所在目录，见4.2
3. resultDir: 结果文件保存路径，需要先makedir建好
4. pipelineDir: 流程中间文件保存路径，这里设置和resultDir一致
5. genomeBed: 基因组bed12格式的注释文件，第二步提取UMI所用，见4.1
6. Samples: 单细胞测序原始数据信息

##### 2 构建df

snakemake的开头部分(rule之前),主要通过**df.pipe**和**df.assign**方法，将**config.yaml**中的信息输入，为后续步骤构建了dataframe索引

##### 3 cellranger比对

使用**cellranger**进行细胞定量，输出文件保存在 `/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step1_cellRanger/nodule_large/nodule_large/outs`

##### 4 提取UMI

输出文件保存在`/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step2_parseUmiDr`

##### 5 velocyto获取loom文件

输出文件保存在`/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step1_cellRanger/nodule_large/nodule_large/velocyto`



[1]: https://github.com/ZhaiLab-SUSTech/soybean_sn_st
[2]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/Schematic_diagram.png
[3]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/sn_fig3.jpg
[4]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/st_fig4.png
[5]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/cluster_fig5.jpg
[6]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/GUS_fig6.jpg
[7]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/CIZ_fig7.png
[8]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/CIZ_fig8.png
[9]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/IC_fig9.png
[10]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/IC_fig10.png
[11]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/IC_fig11.png
[12]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/IC_fig12.png
[13]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/IC_fig13.png
[14]: https://www.10xgenomics.com/support/software/cell-ranger/downloads
[15]: https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/Rtools/history.html
[16]: https://github.com/theislab/batchglm
[17]: https://github.com/theislab/diffxpy
[18]: https://data.jgi.doe.gov/refine-download/phytozome?genome_id=447&expanded=Phytozome-447
[19]: https://data.jgi.doe.gov/refine-download/phytozome?genome_id=275&expanded=Phytozome-275
[20]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/indexResult_fig20.png
[21]: https://github.com/ZhaiLab-SUSTech/soybean_sn_st/blob/main/main/snakemake_cellranger/snakefile
[22]: https://github.com/Mikotoo/Mikotoo.github.io/blob/main/code/cellRanger/
[23]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/dir_fig23.png
[24]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/radata_fig24.png
[25]: https://github.com/liuzj039/jpy_tools/blob/master/tools/singleCell/parseUmiDirectionFromCellrangerBam.py