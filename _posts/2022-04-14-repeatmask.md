---
title: RepeatMasker的安装和使用
author: 
date: 2022-04-14
category: Jekyll
layout: post
---
> 识别并屏蔽基因组中的重复序列<br>
需求:<br>
**python3**环境<br>
TRF 4.09或更高  `conda install trf=4.09`<br>
h5py    `pip install h5py`<br>
[rmblast][2]    下载解压


##### 1、下载RepeatMasker安装包并解压
[官网链接][1]

##### 2、RepBase library
[baidu网盘][3] 提取码 b4to
下载，解压并复制到 `~/software/RepeatMasker/Libraries/` 目录中

##### 3、安装
```
cd ~/softwareRepeatMasker
perl ./configure
```
配置选2，输入rmblast路径，之后一路回车后键入5，完成配置


#### 5、使用
`RepeatMasker -pa 4 -species soybean -poly -html -gff -dir rm genome.fa 1>pro.log 2>err.log`

示例
> 大豆基因组

![pic4][4]

> Repeatmask之后的大豆基因组

![pic5][5]

> 识别到的重复序列

![pic6][6]


[1]: https://www.repeatmasker.org/RepeatMasker/
[2]: http://www.repeatmasker.org/RMBlast.html
[3]: https://pan.baidu.com/s/1I5D1K5S4UeLXNrJFrEed2A
[4]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog3_repeatmask/rm1.png
[5]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog3_repeatmask/rm2.png
[6]: https://github.com/Mikotoo/Mikotoo.github.io/raw/main/downloads/image/blog3_repeatmask/rm3.png