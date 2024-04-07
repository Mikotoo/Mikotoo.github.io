---
title: 【文章复现】 2023-Nature plants-大豆根瘤单细胞与空转
author:
date: 2024-04-07
category: Jekyll
layout: post
---

> 单细胞、空转分析流程学习，整理备查，顺便记录 <br>
> 本文作者是南方科技大学翟继先老师，文章相关的数据和[代码][1]已经公开。


#### 1、文献概述
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


##### 1.3 结果

###### 1）单细胞核RNA测序分析
1. 对3样本数据分别比对，统计细胞数、基因数等信息
2. 将3组dataset整合，得到了 15个cell clusters（cluster0-cluster14），并得到了每个cluster中上调的基因

	![pic3][3]

3. 利用拟南芥和其它豆科植物的标记基因对细胞进行注释，注释出了15个cluster中的5个，包括
	1. 根表皮细胞 cluster 5
	2. 根维管束 cluster 3
	3. 根瘤维管束 cluster 9
	4. 根瘤皮层 cluster 1
	5. 感染细胞	cluster 12

###### 2）Stereo-seq
<p>由于大豆根瘤中标记基因的缺乏，很多cluster没有注释成功。</p>
<p>为了克服这个问题，对12dpi，21dpi的根瘤做空转，根据空间表达信息对细胞进行注释，将其分为6个区域。</p>

![pic4][4]

> 感染区，内皮层，外皮层(cluster 2 and 4)，表皮，微管


<p> 基于去卷积(deconvolution)的方法，验证了之前注释的cluster，并对其他未注释的cluster完成了注释 </p>

1. 未感染细胞 cluster 0; cluster7; cluster 11
2. 根瘤外皮层 cluster 2; cluster 4


![pic5][5]

<p>
	为了验证注释的准确性，对细胞特异性基因进行了GUS染色和RNA原位杂交，结果符合预测
</p>

![pic6][6]



[1]: https://github.com/ZhaiLab-SUSTech/soybean_sn_st
[2]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/Schematic_diagram.png
[3]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/sn_fig3.jpg
[4]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/st_fig4.png
[5]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/cluster_fig5.jpg
[6]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog7_soybean_snRNA/GUS_fig6.jpg