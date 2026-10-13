---
title: 'Mfuzz 时序聚类'
sidebar:
  label: 'Mfuzz 时序聚类'
  order: 6
---


## 中文
**作用**：把 WGCNA brown 与 lightcyan1 模块的基因沿结瘤发育序列（root_uni → root_hpi48 → nodule_dpi10 → dpi21 → dpi45）做软聚类，再把同一批簇投影到根序列上。
**入口**：`Mfuzz.R`
**输入**：工作目录中的 `Gifu_all_samples.tpm.tsv`、`MG20_all_samples.tpm.tsv`，以及 `./result/06.modules_expr/module_brown.csv` 与 `module_lightcyan1.csv`。
**输出**：`result/12.Mfuzz_out/` 下的 `Dmin_plot.pdf`、`Mfuzz_stability.pdf`、`Mfuzz_test.pdf`、`cluster1.genes.list`…`cluster4.genes.list`、`Mfuzz_membership.tsv`、`Mfuzz_plot.pdf`、`Mfuzz_plot2.pdf`、`Mfuzz_root_projection_plot.pdf`、`Mfuzz_root_projection_plot2.pdf`。
**运行**：`Rscript Mfuzz.R`（无参数；工作目录须同时持有两张 TPM 矩阵与 `WGCNA.R` 生成的 `result/` 树）。
**工具**：R 4.4.3 + Mfuzz（`mestimate`、`Dmin`、`mfuzz`、`mfuzz.plot`、`mfuzz.plot2`）、Biobase、tidyverse；参数：阶段均值 TPM 取 `log1p`、`filter.NA(thres = 0.25)`、`fill.NA(mode = "mean")`、`filter.std(min.std = 0)`、`Dmin(crange = seq(4, 20), repeats = 3)`、最终 `c = 4`（`set.seed(1234)`）。

<details>
<summary>English</summary>

**Purpose**: Soft-cluster the genes of the WGCNA brown and lightcyan1 modules along the nodule developmental series (root_uni → root_hpi48 → nodule_dpi10 → dpi21 → dpi45) and project the same clusters onto the root series.
**Entry point**: `Mfuzz.R`
**Inputs**: `Gifu_all_samples.tpm.tsv` and `MG20_all_samples.tpm.tsv` in the working directory, plus `./result/06.modules_expr/module_brown.csv` and `module_lightcyan1.csv`.
**Outputs**: under `result/12.Mfuzz_out/` — `Dmin_plot.pdf`, `Mfuzz_stability.pdf`, `Mfuzz_test.pdf`, `cluster1.genes.list` … `cluster4.genes.list`, `Mfuzz_membership.tsv`, `Mfuzz_plot.pdf`, `Mfuzz_plot2.pdf`, `Mfuzz_root_projection_plot.pdf` and `Mfuzz_root_projection_plot2.pdf`.
**Run**: `Rscript Mfuzz.R` (no arguments; the working directory must hold both TPM matrices and the `result/` tree produced by `WGCNA.R`).
**Tools**: R 4.4.3 with Mfuzz (`mestimate`, `Dmin`, `mfuzz`, `mfuzz.plot`, `mfuzz.plot2`), Biobase and tidyverse; settings as coded: `log1p` on stage-mean TPM, `filter.NA(thres = 0.25)`, `fill.NA(mode = "mean")`, `filter.std(min.std = 0)`, `Dmin(crange = seq(4, 20), repeats = 3)`, and the final clustering `c = 4` with `set.seed(1234)`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details>
<summary><code>Mfuzz.R</code> · R · 374 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/06_mfuzz/Mfuzz.R)

```r
library(Mfuzz)
library(tidyverse)

## 1. 导入表达矩阵 ===================
# TPM文件，第一列是gene_id，其余列是样本
Gifu <- read.table("Gifu_all_samples.tpm.tsv",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# root_dpi02和root_hpi48_rep1是一样，保留hpi48和MG20保持一致
Gifu <- Gifu[, !(colnames(Gifu) %in% c("Gifu_root_dpi02_rep1",
                                       "Gifu_root_dpi02_rep2",
                                       "Gifu_root_dpi02_rep3"))]

MG20 <- read.table("MG20_all_samples.tpm.tsv",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

datExpr0 <- cbind(Gifu, MG20)

# 读取目标模块的基因
module1 <- read.csv("./result/06.modules_expr/module_brown.csv", header = T)
genelist1 <- module1[,1]
module2 <- read.csv("./result/06.modules_expr/module_lightcyan1.csv", header = T)
genelist2 <- module2[,1]

# 目标模块基因的表达矩阵
datExpr1 <- datExpr0[rownames(datExpr0) %in% genelist1 | rownames(datExpr0) %in% genelist2 , ]

# 只保留根瘤、接菌根和根的表达数据进行Mfuzz
noduleTraitList <- c("root_hpi48","root_uni","nodule") 

# 筛选列名中含有 noduleTraitList 中任一关键词的列
keep_cols <- grep(paste(noduleTraitList, collapse = "|"), colnames(datExpr0), value = TRUE)

# 提取子集
datExpr2 <- datExpr1[, keep_cols]

# 确保是 matrix
datExpr2 <- as.matrix(datExpr2)

# -------------------------
# 1) 样本信息解析
# -------------------------
stage_order <- c(
  "root_uni",
  "root_hpi48",
  "nodule_dpi10",
  "nodule_dpi21",
  "nodule_dpi45"
)

meta <- tibble(sample = colnames(datExpr2)) %>%
  mutate(
    stage = str_extract(sample, paste0(stage_order, collapse = "|"))
  ) %>%
  filter(!is.na(stage))

# 保留目标列
datExpr2 <- datExpr2[, meta$sample, drop = FALSE]

# stage 设定顺序
meta$stage <- factor(meta$stage, levels = stage_order)

# -------------------------
# 2) 按 stage 求均值（跨品种 + 重复）
# -------------------------
expr_mean <- sapply(stage_order, function(st) {
  cols <- meta$sample[meta$stage == st]
  if (length(cols) == 0) {
    stop(paste("ERROR: no samples found for stage:", st))
  }
  rowMeans(datExpr2[, cols, drop = FALSE], na.rm = TRUE)
})

expr_mean <- as.matrix(expr_mean)
rownames(expr_mean) <- rownames(datExpr2)

# -------------------------
# 3) 表达矩阵过滤 & 转换
# -------------------------

# 去掉全 0 基因
expr_mean <- expr_mean[rowSums(expr_mean) > 0, , drop = FALSE]

# 如果是 TPM / FPKM，建议 log1p
expr_mean <- log1p(expr_mean)

# -------------------------
# 4) 构建 ExpressionSet
# -------------------------
eset <- new("ExpressionSet", exprs = expr_mean)

# NA 处理
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = "mean")

# 去掉标准差为 0 的基因
eset <- filter.std(eset, min.std = 0)

# Mfuzz 必须标准化
eset <- standardise(eset)

# -------------------------
# 5) Mfuzz 最佳簇选取
# -------------------------
dir.create("result/12.Mfuzz_out", recursive = TRUE, showWarnings = FALSE)

set.seed(1234)
m <- mestimate(eset)

# Step1: Dmin
pdf("result/12.Mfuzz_out/Dmin_plot.pdf", width = 6, height = 5)

Dmin <- Dmin(eset, m = m, crange = seq(4, 20), repeats = 3, visu = FALSE)

plot(4:20, Dmin, type = "b",
     xlab = "Cluster number (c)",
     ylab = "Min. centroid distance")

dev.off()

# Step2: stability
pdf("result/12.Mfuzz_out/Mfuzz_stability.pdf", width = 6, height = 5)

stability <- sapply(4:20, function(k) {
  cl1 <- mfuzz(eset, c = k, m = m)
  cl2 <- mfuzz(eset, c = k, m = m)
  mean(cl1$cluster == cl2$cluster)
})

plot(4:20, stability, type = "b",
     xlab = "Cluster number (c)",
     ylab = "Clustering stability",
     main = "Mfuzz stability analysis")

dev.off()

# Step3: 查看最优簇聚类
c <- 8

cl <- mfuzz(eset, c = c, m = m)
print(cl$size)

pdf(file.path("result/12.Mfuzz_out", "Mfuzz_test.pdf"), width = 24, height = 8)

mfuzz.plot2(
  eset, cl,
  mfrow = c(2, 4),
  time.labels = colnames(exprs(eset)),
  ylim.set = c(0, 0),
  xlab = "Stage",
  ylab = "Expression changes",
  x11 = FALSE,
  ax.col = "black", bg = "white",
  colo = "fancy",          # ✅ 关键修改
  col.axis = "black",
  col.lab = "black",
  col.main = "black",
  col.sub = "black",
  col = "black",
  centre = FALSE,
  centre.col = "black",
  centre.lwd = 2,
  cex.main = 3,
  cex.lab = 2.5,
  cex.axis = 2,
  Xwidth = 5,
  Xheight = 5
)

dev.off()

# -------------------------
# 6) 根据最佳簇结果，合并为4簇
# -------------------------

c <- 4
m <- mestimate(eset)

set.seed(1234)
cl <- mfuzz(eset, c = c, m = m)

# 查看每个 cluster 的基因数
print(cl$size)

# cluster 基因列表
for (k in 1:c) {
  genes_k <- names(cl$cluster[cl$cluster == k])
  write.table(
    genes_k,
    file = file.path("result/12.Mfuzz_out", paste0("cluster", k, ".genes.list")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# membership 表
write.table(
  data.frame(
    GeneID = rownames(cl$membership),
    cl$membership,
    Cluster = cl$cluster
  ),
  file = "result/12.Mfuzz_out/Mfuzz_membership.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# 7) 可视化
# -------------------------
pdf("result/12.Mfuzz_out/Mfuzz_plot.pdf", width = 12, height = 8)

mfuzz.plot(
  eset,
  cl,
  mfrow = c(2, 2),
  time.labels = colnames(exprs(eset)),
  new.window = FALSE
)

dev.off()

pdf(file.path("result/12.Mfuzz_out", "Mfuzz_plot2.pdf"), width = 12, height = 8)

mfuzz.plot2(
  eset, cl,
  mfrow = c(2, 2),
  time.labels = colnames(exprs(eset)),
  ylim.set = c(0, 0),
  xlab = "Stage",
  ylab = "Expression changes",
  x11 = FALSE,
  ax.col = "black", bg = "white",
  colo = "fancy",          # ✅ 关键修改
  col.axis = "black",
  col.lab = "black",
  col.main = "black",
  col.sub = "black",
  col = "black",
  centre = FALSE,
  centre.col = "black",
  centre.lwd = 2,
  cex.main = 1.6,
  cex.lab = 1.4,
  cex.axis = 1.1,
  Xwidth = 5,
  Xheight = 5
)

dev.off()

# -------------------------
# 8) 染菌根分析
# -------------------------

# 基于上述聚类的基因簇，观察各基因在染菌根的表达趋势
# 只保留根瘤、接菌根和根的表达数据进行Mfuzz
rootTraitList <- c("root_hpi48","root") 

# 筛选列名中含有 rootTraitList 中任一关键词的列
keep_cols1 <- grep(paste(rootTraitList, collapse = "|"), colnames(datExpr0), value = TRUE)

# 提取子集
datExpr3 <- datExpr1[, keep_cols1]

# 确保是 matrix
datExpr3 <- as.matrix(datExpr3)

# root 时间顺序
stage_order1 <- c(
  "root_uni",
  "root_hpi48",
  "root_dpi10",
  "root_dpi21",
  "root_dpi45"
)

meta_root <- tibble(sample = colnames(datExpr3)) %>%
  mutate(
    stage = str_extract(sample, paste0(stage_order1, collapse = "|"))
  ) %>%
  filter(!is.na(stage))

# 设置 stage 顺序
meta_root$stage <- factor(meta_root$stage, levels = stage_order1)

# 按 root stage 求均值（跨品种 + 重复）
expr_root_mean <- sapply(stage_order1, function(st) {
  cols <- meta_root$sample[meta_root$stage == st]
  if (length(cols) == 0) {
    stop(paste("ERROR: no samples found for stage:", st))
  }
  rowMeans(datExpr3[, cols, drop = FALSE], na.rm = TRUE)
})

expr_root_mean <- as.matrix(expr_root_mean)
rownames(expr_root_mean) <- rownames(datExpr3)

# 去掉全 0
expr_root_mean <- expr_root_mean[rowSums(expr_root_mean) > 0, , drop = FALSE]

# log 转换（TPM）
expr_root_mean <- log1p(expr_root_mean)

# nodule Mfuzz 的基因集合
mfuzz_genes <- names(cl$cluster)

# 在 root 表达矩阵中取交集
expr_root_mfuzz <- expr_root_mean[
  rownames(expr_root_mean) %in% mfuzz_genes,
  ,
  drop = FALSE
]

# 保证顺序与 cl 一致（非常重要）
expr_root_mfuzz <- expr_root_mfuzz[mfuzz_genes, ]

# 不重新聚类，构建eset_root
eset_root <- new("ExpressionSet", exprs = expr_root_mfuzz)

# NA 处理
eset_root <- filter.NA(eset_root, thres = 0.25)
eset_root <- fill.NA(eset_root, mode = "mean")

# 去掉标准差为 0
eset_root <- filter.std(eset_root, min.std = 0)

# 标准化（仅用于画趋势）
eset_root <- standardise(eset_root)

# 可视化
pdf("result/12.Mfuzz_out/Mfuzz_root_projection_plot.pdf",
    width = 12, height = 8)

mfuzz.plot(
  eset_root,
  cl,
  mfrow = c(2, 2),
  time.labels = colnames(exprs(eset_root)),
  new.window = FALSE
)

dev.off()

pdf("result/12.Mfuzz_out/Mfuzz_root_projection_plot2.pdf",
    width = 12, height = 8)

mfuzz.plot2(
  eset_root, cl,
  mfrow = c(2, 2),
  time.labels = colnames(exprs(eset)),
  ylim.set = c(0, 0),
  xlab = "Stage",
  ylab = "Expression changes",
  x11 = FALSE,
  ax.col = "black", bg = "white",
  colo = "fancy",          # ✅ 关键修改
  col.axis = "black",
  col.lab = "black",
  col.main = "black",
  col.sub = "black",
  col = "black",
  centre = FALSE,
  centre.col = "black",
  centre.lwd = 2,
  cex.main = 1.6,
  cex.lab = 1.4,
  cex.axis = 1.1,
  Xwidth = 5,
  Xheight = 5
)

dev.off()
```

</details>
