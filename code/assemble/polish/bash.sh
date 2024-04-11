#nohup ./polish.sh 10 2 soybean_chr_v0.8.fa /pub5/ZC/soybean/ONT/ont.fasta /pub5/ZC/soybean/WGS/WGS.merylDB ontPolish 1>ontPolish.log 2>&1 &
nohup ./polish.sh 10 2 ontPolish.iter_2.consensus.fasta /pub4/zhangchao/soybean/hifi_reads/soybean_ccs.fasta /pub4/zhangchao/T2T/assemble/data/BGI.meryl hifiPolish 1>hifiPolish.log 2>&1 &
