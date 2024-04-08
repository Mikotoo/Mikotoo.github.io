---
title: 【文章复现】 2023-Nature plants-大豆根瘤单细胞与空转
author:
date: 2024-04-07
category: Jekyll
layout: post
---

> 单细胞、空转分析流程学习，整理备查，顺便记录 <br>
> 本文作者是南方科技大学翟继先老师，文章相关的数据和[代码][1]已经公开。


### 1、文献概述
##### 1.1 研究目的
<p>大豆的共生固氮器官根瘤作为一种高度异质的组织，其各种类型的细胞都具有不同的生理特点和功能。</p>
<p>本文的研究目的是了解根瘤中不同类型的细胞在根瘤成熟过程中的具体贡献，以及细胞之间的关系。</p>

##### 1.2 建库测序
1）为了揭示根瘤成熟过程中的cell-type-specific动态基因表达，选取3个样本进行了单细胞核测序(10X Genomics Chromium): <br>

1. 接种12天的根瘤 (sn_12dpi)
2. 接种21天的根瘤 (sn_21dpi)
3. 作为对照，接种21天根瘤附近的根 (sn_root)

2）为了更好的进行细胞注释，又对上述几个时期的材料进行了空间转录组测序(stereo-seq) <br>

<i><b>两种技术结合，共同完成了细胞分类 </b></i><br>

![pic2][2]


##### 1.3 聚类与注释

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


##### 1.4 感染区亚型


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



### 2、测序原理
##### 2.1 单细胞转录组测序
##### 2.2 空间转录组测序

### 3、软件安装
##### 3.1 Cell Ranger
> 





















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