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
大豆的共生固氮器官根瘤作为一种高度异质的组织，其各种类型的细胞都具有不同的生理特点和功能。<br>
本文的研究目的是了解根瘤中不同类型的细胞在根瘤成熟过程中的具体贡献，以及细胞之间的关系。<br>

##### 1.2 建库测序
为了揭示根瘤成熟过程中的cell-type-specific动态基因表达，选取3个样本进行了单细胞核测序(10X Genomics Chromium): <br>

1. 接种12天的根瘤 (sn_12dpi)
2. 接种21天的根瘤 (sn_21dpi)
3. 作为对照的根 (sn_root) 
> 根的生长时期和21天根瘤相同

和空间转录组测序(stereo-seq)


[1]: https://github.com/ZhaiLab-SUSTech/soybean_sn_st