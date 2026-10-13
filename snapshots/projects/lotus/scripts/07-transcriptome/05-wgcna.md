---
title: 'WGCNA 共表达模块'
sidebar:
  label: 'WGCNA 共表达模块'
  order: 5
---


## 中文
**作用**：由 Gifu 与 MG20 组织图谱 TPM 矩阵构建加权基因共表达网络，检测并合并模块，把模块特征基因与组织性状关联，并提取 hub 基因。
**入口**：`WGCNA.R`
**输入**：上一级目录的 `Gifu_all_samples.tpm.tsv`、`MG20_all_samples.tpm.tsv`（`gene_id` 为行名、每样本一列）；样本分组从列名解析。
**输出**：`result/` 下的样本聚类、软阈值、基因聚类、模块检测与合并 PDF；`06.modules_expr/module_<color>.csv`（每模块的基因 × 样本矩阵）；`07.module_trait_correlation.csv`、`08.module_trait_pvalue.csv`、`07.trait_*.csv`、`09.module_trait_heatmap.pdf`；`10.hub_gene_results/`、`11.cytoscape/`，以及每模块一张 `12.gene_trait_combined_network_<module>.pdf`。
**运行**：`Rscript WGCNA.R`（无参数；在打算生成 `result/` 的目录内运行，两张 TPM 矩阵放在该目录的上一级）。
**工具**：R 4.4.3 + WGCNA（`adjacency`、`TOMsimilarity`、`cutreeDynamic`、`mergeCloseModules`、`exportNetworkToCytoscape`）、tidyverse/dplyr、tidygraph/igraph/ggraph、pheatmap；参数：软阈值 power 14、`minModuleSize` 30、`deepSplit` 2、合并 `cutHeight` 0.3、hub 阈值 MM ≥ 0.7 且 GS ≥ 0.5、基因过滤为至少一个条件组均值 TPM > 0.5 后取 MAD 前 75%。

<details>
<summary>English</summary>

**Purpose**: Build a weighted gene co-expression network from the Gifu and MG20 tissue-atlas TPM matrices, detect and merge modules, relate module eigengenes to tissue traits and extract hub genes.
**Entry point**: `WGCNA.R`
**Inputs**: `Gifu_all_samples.tpm.tsv` and `MG20_all_samples.tpm.tsv` one directory above the working directory (`gene_id` as row names, one column per sample); sample groups are parsed from the column names.
**Outputs**: under `result/`, the sample-clustering, soft-threshold, gene-clustering, module-detection and module-merge PDFs; `06.modules_expr/module_<color>.csv` (per-module gene × sample matrix); `07.module_trait_correlation.csv`, `08.module_trait_pvalue.csv`, `07.trait_*.csv`, `09.module_trait_heatmap.pdf`; `10.hub_gene_results/`, `11.cytoscape/`, and one `12.gene_trait_combined_network_<module>.pdf` per module.
**Run**: `Rscript WGCNA.R` (no arguments; run it from the directory intended to hold `result/`, with the two TPM matrices in that directory's parent).
**Tools**: R 4.4.3 with WGCNA (`adjacency`, `TOMsimilarity`, `cutreeDynamic`, `mergeCloseModules`, `exportNetworkToCytoscape`), tidyverse/dplyr, tidygraph/igraph/ggraph, pheatmap; settings as coded: soft power 14, `minModuleSize` 30, `deepSplit` 2, merge `cutHeight` 0.3, hub thresholds MM ≥ 0.7 and GS ≥ 0.5, gene filter = group mean TPM > 0.5 in at least one condition, then the top 75% by MAD.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details>
<summary><code>WGCNA.R</code> · R · 464 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/05_wgcna/WGCNA.R)

```r
#=================== step0：初始化 ===================
rm(list = ls())
library(WGCNA)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(pheatmap)
library(igraph)
library(scales) 
library(ggraph)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  # 启用多线程

# 设置工作目录
dir.create('result', showWarnings = FALSE)
setwd("result")

#=================== step1：导入表达矩阵 ===================
# TPM文件，第一列是gene_id，其余列是样本
Gifu <- read.table("../Gifu_all_samples.tpm.tsv",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# root_dpi02和root_hpi48_rep1是一样，保留hpi48和MG20保持一致
Gifu <- Gifu[, !(colnames(Gifu) %in% c("Gifu_root_dpi02_rep1",
                                       "Gifu_root_dpi02_rep2",
                                       "Gifu_root_dpi02_rep3"))]

MG20 <- read.table("../MG20_all_samples.tpm.tsv",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

datExpr0 <- cbind(Gifu, MG20)

# 进行过滤
# 保留至少在一种样本中，3次重复均值大于0.5的基因
# 获取分组信息
sample_groups <- sub("_rep[0-9]+$", "", colnames(datExpr0))

# 按 sample group 计算均值
group_means <- sapply(
  unique(sample_groups),
  function(g) {
    rowMeans(datExpr0[, sample_groups == g, drop = FALSE])
  }
)

# 确保矩阵结构正确
group_means <- as.matrix(group_means)
# 至少有一个样本条件，3 重复均值 > 0.5
expr_filter <- apply(group_means, 1, function(x) any(x > 0.5))
datExpr_filt <- datExpr0[expr_filter, ]
# 计算基因方差（跨所有样本）
gene_var <- apply(datExpr_filt, 1, mad)
#选取方差最大的前 75% 个基因
top_n <- min(length(gene_var)*0.75, length(gene_var))

top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]

datExpr1 <- datExpr_filt[top_genes, ]

# 检查是否有缺失值
gsg <- goodSamplesGenes(datExpr1, verbose = 3)
if (!gsg$allOK) {
  datExpr1 <- datExpr1[gsg$goodSamples, gsg$goodGenes]
}
datExpr <- apply(datExpr1,1,function(x)log2(x+1))
#=================== step2：层次聚类检测异常样本 ===================
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("01.sample_clustering.pdf", width = 10, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()

# 删除离群样本 MG20_flower_rep3
outlierSample <- "MG20_flower_rep3"
if (outlierSample %in% rownames(datExpr)) {
  datExpr <- datExpr[!rownames(datExpr) %in% outlierSample, ]
  cat(paste("Removed outlier sample:", outlierSample, "\n"))
}

# 重新绘制聚类图
sampleTree = hclust(dist(datExpr), method="average")
pdf("01_SampleClustering_no_outlier.pdf", width=12, height=8)
plot(sampleTree, main="Sample clustering (outlier removed)", sub="", xlab="")
dev.off()

#=================== step3：读取样本分组信息 ===================
# 从样本名中提取组织信息，例如 Gifu_leaf_rep1 → Gifu_leaf
sampleNames <- rownames(datExpr)
sampleType <- sub("_rep[0-9]+", "", sampleNames)
sampleType <- sub("_dpi[0-9]+", "", sampleType)
sampleType <- sub("_hpi[0-9]+", "", sampleType)
sampleType <- sub("Gifu_", "", sampleType)
sampleType <- sub("MG20_", "", sampleType)
traitData <- data.frame(row.names = sampleNames, tissue = sampleType)

# 将类别变量转为数值矩阵（适用于WGCNA相关性计算）
traitColors <- data.frame(model.matrix(~ 0 + tissue, data = traitData))
colnames(traitColors) <- sub("tissue", "", colnames(traitColors))

#=================== step4：选择软阈值参数 β ===================
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("02.pickSoftThreshold.pdf", width = 9, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.9, col="red")  # 通常选择R^2>0.8或0.9的power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# 设定软阈值（根据上图手动选择）
# softPower <- sft$powerEstimate
softPower <- 14

#=================== step5：构建加权网络并检测模块 ===================
adjacency <- adjacency(datExpr, power = softPower)

# 计算TOM矩阵
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf("03.gene_clustering.pdf", width = 12, height = 8)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# 动态剪切聚类
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,  #deepSplit取值范围[0,4]，越大切割越精细，模块越多
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

# 绘图
pdf("04.module_detection.pdf", width = 12, height = 8)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=================== step6：模块合并 ===================
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

pdf("05.module_merge.pdf", width = 10, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.3, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.3, verbose = 3) # cutHeight：合并阈值（越低合并越保守）
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
moduleColors <- mergedColors
MEs <- mergedMEs

names(mergedColors) <- colnames(datExpr)
output_dir <- "06.modules_expr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
moduleColors <- unique(mergedColors)

for (color in moduleColors) {
  geneIDs <- mergedColors == color
  geneNames <- names(mergedColors)[geneIDs]  # 对应的基因名

  # 从 datExpr 中按列名提取对应基因表达数据（样本 × 基因）
  moduleExpr <- datExpr[, geneNames]

  # 转置为 基因 × 样本，输出为 CSV
  moduleExpr <- t(moduleExpr)
  write.csv(moduleExpr,
            file = file.path(output_dir, paste0("module_", color, ".csv")),
            quote = FALSE)
}

#=================== step7：模块与样本性状的相关性 ===================
moduleTraitCor <- cor(MEs, traitColors, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# 输出结果
write.csv(moduleTraitCor, file = "07.module_trait_correlation.csv")
write.csv(moduleTraitPvalue, file = "08.module_trait_pvalue.csv")

# 提取每种性状相关性最高的module
moduleTraitCor <- as.data.frame(moduleTraitCor)

# 对每一列取最大的行名和值
traitModuleBest <- apply(moduleTraitCor, 2, function(x) {
  idx <- which.max(x)
  data.frame(
    Trait = names(x)[1],              # 列名（后面会被覆盖）
    module = rownames(moduleTraitCor)[idx],
    R = x[idx]
  )
})

# 合并结果
traitModuleBest <- do.call(rbind, traitModuleBest)

# 修正 Trait 名（列名）
traitModuleBest$Trait <- colnames(moduleTraitCor)

# 查看结果
write.csv(traitModuleBest, file = "07.trait_most_relevant_module.csv", row.names = F)

# 提取每个性状相关性大于0.6的module
# 对每一列取最大的行名和值
traitModuleOptional <- do.call(
  rbind,
  lapply(colnames(moduleTraitCor), function(trait) {
    x <- moduleTraitCor[, trait]
    idx <- which(x > 0.6)

    if (length(idx) == 0) return(NULL)

    data.frame(
      Trait  = trait,
      module = rownames(moduleTraitCor)[idx],
      R      = x[idx],
      row.names = NULL
    )
  })
)

write.csv(traitModuleOptional, file = "07.trait_Optional_high_correlation_modules.csv", row.names = F)

# 可视化热图
pdf("09.module_trait_heatmap.pdf", width = 10, height = 9)
moduleTraitCor <- as.matrix(moduleTraitCor)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #需要显示数值再放
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),   # 用实际 trait 名称
  yLabels = rownames(moduleTraitCor),   # 用模块名称
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  mar = c(16, 18, 3, 3),     # ⬅️ 增大下边距和左边距
  cex.text = 1.0,              
  cex.lab.x = 1.0,           # X轴标签大小
  cex.lab.y = 1.2,           # Y轴标签稍大一点
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)
dev.off()
#=================== step8：Hub gene 分析、绘图与输出 ===================
# 设置阈值
mm_threshold <- 0.7
gs_threshold <- 0.5

# 设置输出目录
output_dir <- "10.hub_gene_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
dir.create(plot_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

# ===== 样本和基因数量 =====
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

# 模块与性状名称
modNames <- substring(names(MEs), 3)
traitNames <- names(traitColors)

# 计算模块归属度（MM）和 p 值
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM.", modNames)
names(MMPvalue) <- paste0("p.MM.", modNames)

# 计算性状相关性（GS）和 p 值
geneTraitSignificance <- as.data.frame(cor(datExpr, traitColors, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste0("GS.", traitNames)
names(GSPvalue) <- paste0("p.GS.", traitNames)

# 汇总信息
geneInfo <- data.frame(
  Gene = colnames(datExpr),
  ModuleColor = mergedColors
) %>%
  bind_cols(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)

# 遍历每个模块 × 性状，画图 + 筛选 hub gene
for (trait in traitNames) {
  for (module in modNames) {
    inModule <- mergedColors == module
    moduleGenes <- geneInfo[inModule, ]

    mm_col <- paste0("MM.", module)
    gs_col <- paste0("GS.", trait)

    if (nrow(moduleGenes) > 1) {
      # 画图
      pdf(file = file.path(plot_dir, paste0("MM_vs_GS_", trait, "_", module, ".pdf")),
          width = 7, height = 7)

      verboseScatterplot(
        x = abs(moduleGenes[[mm_col]]),
        y = abs(moduleGenes[[gs_col]]),
        xlab = paste("Module Membership in", module, "module"),
        ylab = paste("Gene significance for", trait),
        main = paste("MM vs. GS\nModule:", module, "| Trait:", trait),
        col = module,
        abline = FALSE
      )
      abline(h = gs_threshold, col = "red", lty = 2)
      abline(v = mm_threshold, col = "red", lty = 2)
      dev.off()

      # 筛选 hub genes
      hub_genes <- moduleGenes %>%
        filter(abs(!!sym(mm_col)) >= mm_threshold & abs(!!sym(gs_col)) >= gs_threshold)
      cat("Trait:", trait, "Module:", module, "Hub genes found:", nrow(hub_genes), "\n")
      if (nrow(hub_genes) > 0) {
        write.csv(hub_genes,
                  file = file.path(table_dir, paste0("hubgenes_", trait, "_", module, ".csv")),
                  row.names = FALSE)
      }
    }
  }
}

### 每个模块的hub基因热图
# 设置输入输出路径
hubgene_dir <- "10.hub_gene_results/tables"
heatmap_dir <- "10.hub_gene_results/heatmaps_1"
if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir)

# 加载表达矩阵（样本 × 基因）
expr_mat <- datExpr  # log2(TPM+1) 表达量矩阵，列为基因，行为样本

# 遍历每个 hub gene 文件
hubgene_files <- list.files(hubgene_dir, pattern = "^hubgenes_.*\\.csv$", full.names = TRUE)

for (file in hubgene_files) {
  hub_df <- read.csv(file)
  if (nrow(hub_df) == 0) next  # 跳过空文件

  gene_list <- hub_df$Gene
  gene_list <- intersect(gene_list, colnames(expr_mat))  # 确保存在于表达矩阵中

  if (length(gene_list) >= 2) {
    mat <- expr_mat[, gene_list, drop = FALSE]

    pdf(file = file.path(heatmap_dir, paste0(tools::file_path_sans_ext(basename(file)), ".pdf")),
        width = 20, height = max(4, ncol(mat) * 0.2))

    pheatmap(t(mat),  # 行是基因，列是样本
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "none",  # 按行标准化
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 6)

    dev.off()
  }
}


##
# 创建输出目录
dir.create('11.cytoscape', recursive = TRUE)

moduleColorsWW <- mergedColors
TOM_matrix <- TOM  # TOM 是前面计算得到的相似性矩阵

# 获取所有模块名称
modules_all <- names(table(moduleColorsWW))

# 遍历每个模块
for (mod in modules_all) {
  # 获取该模块中基因
  probes <- colnames(datExpr)  # 原始表达矩阵的基因名（列名）
  inModule <- (moduleColorsWW == mod)
  modProbes <- probes[inModule]
  modGenes <- modProbes  # 可替换为其它ID，比如注释后的基因名

  # 获取该模块对应的 TOM 相似性矩阵
  modTOM <- TOM_matrix[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)

  # 输出文件路径
  outEdge <- file.path("11.cytoscape", paste0(mod, ".edge_list.txt"))
  outNode <- file.path("11.cytoscape", paste0(mod, ".node_list.txt"))

  # 导出到 Cytoscape：权重阈值 threshold 可调节边数量
  exportNetworkToCytoscape(
    modTOM,
    edgeFile = outEdge,
    nodeFile = outNode,
    weighted = TRUE,
    threshold = 0.3,
    nodeNames = modProbes,
    altNodeNames = modGenes,
    nodeAttr = moduleColorsWW[inModule]
  )
}

for (mod in modules_all) {
  mod_genes <- names(mergedColors[mergedColors == mod])
  modExpr <- datExpr[, mod_genes]

  cor_trait_gene <- cor(traitColors, modExpr, method = "spearman", use = "p")

  trait_gene_edges <- as.data.frame(as.table(cor_trait_gene)) %>%
    rename(from = Var1, to = Var2, weight = Freq) %>%
    filter(abs(weight) > 0.4)  # 降低阈值

  mod_idx <- colnames(datExpr) %in% mod_genes
  tom_sub <- TOM[mod_idx, mod_idx]
  gene_gene_edges <- as.data.frame(as.table(tom_sub)) %>%
    rename(from = Var1, to = Var2, weight = Freq) %>%
    filter(from != to, weight > 0.05)  # 降低阈值

  edge_all <- bind_rows(
    mutate(trait_gene_edges, type = "trait-gene"),
    mutate(gene_gene_edges, type = "gene-gene")
  )

  if(nrow(edge_all) == 0) {
    message(mod, " skipped due to empty network")
    next
  }

  all_nodes <- unique(c(edge_all$from, edge_all$to))
  nodes <- tibble(name = all_nodes) %>%
    mutate(type = ifelse(name %in% colnames(traitColors), "Trait", "Gene"))

  if(nrow(nodes) == 0){
    message(mod, " skipped due to empty nodes")
    next
  }

  graph_all <- tbl_graph(nodes = nodes, edges = edge_all, directed = FALSE)

  pdf(file.path(getwd(), paste0("12.gene_trait_combined_network_", mod, ".pdf")), width = 10, height = 8)
  ggraph(graph_all, layout = "fr") +
    geom_edge_link(aes(edge_alpha = abs(weight), color = type), show.legend = TRUE) +
    geom_node_point(aes(color = type), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c("trait-gene" = "red", "gene-gene" = "steelblue")) +
    theme_void()
  dev.off()
  message(mod, " done")
}
```

</details>
