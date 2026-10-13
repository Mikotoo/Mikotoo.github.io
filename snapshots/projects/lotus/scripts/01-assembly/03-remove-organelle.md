---
title: '去除细胞器序列'
sidebar:
  label: '去除细胞器序列'
  order: 3
---


## 中文
**作用**：将 contig 比对 NCBI `nt`，经 GI → taxid → 学名汇总每个物种的命中数，据此识别并去除质体/线粒体来源的 contig。
**入口**：`03.1.blast.sh`（内部调用 `03.2.get_gi.py`、`03.3.get_name.py`），随后运行 `03.4.get_ratio.py`
**输入**：`contigs.fa`；`nt` 库 ``{SOFTWARE_ZC}`/nt_lib/nt`；`cutted_nucl_gb.accession2taxid` 与 `names.dmp`（须放在工作目录）。
**输出**：`nt.blast.xml`、`tiqu_gi.txt`、`scientific_name.txt`、`final_result.txt`（`Name/Hit_reads/percentage1/percentage2`，按命中数排序）。
**运行**：`bash 03.1.blast.sh`，然后 `python 03.4.get_ratio.py`
**工具**：`blastn -outfmt 5 -evalue 1e-5 -num_threads 88`；Python 仅用标准库（`re`、`collections`），以 `python` 调用。

<details>
<summary>English</summary>

**Purpose**: BLAST the contigs against NCBI `nt`, resolve GI → taxid → scientific name and count hits per species to identify and remove plastid/mitochondrial (organelle) contigs.
**Entry point**: `03.1.blast.sh` (calls `03.2.get_gi.py` and `03.3.get_name.py`), after which `03.4.get_ratio.py` runs separately
**Inputs**: `contigs.fa`; `nt` database ``{SOFTWARE_ZC}`/nt_lib/nt`; `cutted_nucl_gb.accession2taxid` and `names.dmp` (must be in the working directory).
**Outputs**: `nt.blast.xml`, `tiqu_gi.txt`, `scientific_name.txt`, `final_result.txt` (`Name/Hit_reads/percentage1/percentage2`, sorted by hit count).
**Run**: `bash 03.1.blast.sh`, then `python 03.4.get_ratio.py`
**Tools**: `blastn -outfmt 5 -evalue 1e-5 -num_threads 88`; Python with the standard library only (`re`, `collections`), invoked as `python`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>03.1.blast.sh</code> · Shell · 6 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/03_remove_organelle/03.1.blast.sh)

```bash
#3.1 排除长度较长（10Mb）的染色体序列，提取其余序列（可能是质体），比对到nt库（需要从ncbi下载到本地），输出fmt5格式

blastn -query contigs.fa -out nt.blast.xml -db {SOFTWARE_ZC}/nt_lib/nt -outfmt 5 -evalue 1e-5 -num_threads 88

python 03.2.get_gi.py
python 03.3.get_name.py
```

</details>

<details open>
<summary><code>03.2.get_gi.py</code> · Python · 31 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/03_remove_organelle/03.2.get_gi.py)

```python
# -*- coding:utf-8 -*-
import re
from collections import defaultdict

xmlfile=open("nt.blast.xml","r")
outfile=open("tiqu_gi.txt","w")

dict1=defaultdict(list)
for lines in xmlfile:
	line=lines.strip()
	read_id = re.match('<Iteration_query-def>.*</Iteration_query-def>',line)
	Hit_id = re.match('<Hit_id>.*</Hit_id>',line)
	Hit_def = re.match('<Hit_def>.*</Hit_def>',line)
	
	if read_id !=None:
		read_id=read_id.group()
		read_id = read_id.split("<")[1].split(">")[1]
		key=read_id

	elif Hit_id !=None:
		Hit_id = Hit_id.group()
		Hit_id = Hit_id.split("<")[1].split(">")[1]
		dict1[key].append(Hit_id)

	elif Hit_def !=None:
		Hit_def = Hit_def.group()
		Hit_def = Hit_def.split("<")[1].split(">")[1]
		dict1[key].append(Hit_def)

for key in dict1:
	outfile.write(key + "\t" + "\t".join(dict1[key])+"\n")
```

</details>

<details open>
<summary><code>03.3.get_name.py</code> · Python · 42 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/03_remove_organelle/03.3.get_name.py)

```python
# -*- coding:utf-8 -*-

tiqu_gi = open("tiqu_gi.txt","r")
gi2taxid = open("cutted_nucl_gb.accession2taxid","r")
taxid2name = open("names.dmp","r")
get_name = open("scientific_name.txt","w")
from collections import defaultdict

taxid_name_dict={}
for lines in taxid2name:
	if "scientific name" in lines:
		line = lines.strip().split("|")
		taxid = line[0].strip()
		name = line[1].strip()
		taxid_name_dict[taxid]=name
		
tiqu_dict=defaultdict(list)
for lines in tiqu_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	tiqu_dict[gi].append("\t".join(line))


gi_taxid_dict={}
for lines in gi2taxid:
	line = lines.strip().split("\t")
	GI = line[1]
	taxid = line[0]
	gi_taxid_dict[GI]=taxid

jiaoji=set(tiqu_dict.keys())&set(gi_taxid_dict.keys())

tax_list=taxid_name_dict.keys()

tiqu_gi = open("tiqu_gi.txt","r")
for lines in tiqu_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	if gi in jiaoji:
		taxid=gi_taxid_dict[gi]
		if taxid in tax_list:
			get_name.write("\t".join(line)+"\t"+taxid_name_dict[taxid]+"\n")
```

</details>

<details open>
<summary><code>03.4.get_ratio.py</code> · Python · 24 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/03_remove_organelle/03.4.get_ratio.py)

```python
# -*- coding:utf-8 -*-
from collections import Counter

scientific_name=open("scientific_name.txt","r")
final_result =open("final_result.txt","w")

name_list_all=[]
for lines in scientific_name:
	line = lines.strip().split("\t")
	name = line[-1]
	name_list_all.append(name)

count_result = Counter(name_list_all)
count_list = list(count_result.items())
count_list.sort(key=lambda x:x[1],reverse=True)

final_result.write("Name\tHit_reads\tpercentage1\tpercentage2\n")
for i in count_list:
	name = i[0]
	number = i[1]
	reads_num = 100000
	percentage1 = "%.2f%%"%(100*float(number)/float(reads_num))
	percentage2 ="%.2f%%"%(100*float(number)/float(len(name_list_all)))
	final_result.write(name+"\t"+str(number)+"\t"+str(percentage1)+"\t"+str(percentage2)+"\n")
```

</details>
