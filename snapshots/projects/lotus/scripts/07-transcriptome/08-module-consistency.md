---
title: '模块一致性'
sidebar:
  label: '模块一致性'
  order: 8
---


## 中文
**作用**：在 12 个配对条件上比较 WGCNA brown 与 lightcyan1 模块基因在 Gifu 与 MG20 之间的逐基因表达模式，以 Spearman 为主、Pearson 为辅，写出逐基因表与相关分布图。
**入口**：`module_gene_expression_correlation.R`
**输入**：`Gifu_all_samples_mean.tsv`、`MG20_all_samples_mean.tsv`（`gene_id` 列加 12 个条件列，脚本剥离 `Gifu_`/`MG20_` 前缀）；`result/06.modules_expr/module_brown.csv` 与 `module_lightcyan1.csv`。
**输出**：`result/13.module_gene_expression_correlation/` 下的 `brown_lightcyan1_gene_expression_correlations.tsv` 与 `.filtered.tsv`、`correlation_threshold_summary.tsv`、`spearman_correlation_distribution.pdf|.png`、`pearson_correlation_distribution.pdf|.png` 与 `README.txt`。
**运行**：`Rscript module_gene_expression_correlation.R`（无参数；脚本按自身目录解析路径，须与两张均值矩阵及 `result/` 树并列并从该目录运行）。
**工具**：R 4.4.3 + ggplot2/scales（`ggsave(device = cairo_pdf)`）与 base R（`cor.test(method = "spearman"/"pearson")`）；参数：相关基于 12 个条件下的 log2(均值 TPM + 1)，高相关阈值 R ≥ 0.8，基因需在至少 2 个条件下均值 TPM ≥ 0.5 且两向量 SD > 0。

<details>
<summary>English</summary>

**Purpose**: Compare the per-gene expression pattern of the WGCNA brown and lightcyan1 module genes between Gifu and MG20 across the 12 matched conditions, using Spearman as the primary statistic and Pearson as an auxiliary one, and write the per-gene tables and correlation distribution plots.
**Entry point**: `module_gene_expression_correlation.R`
**Inputs**: `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv` (a `gene_id` column plus the 12 condition columns; the `Gifu_`/`MG20_` prefix is stripped in the script); `result/06.modules_expr/module_brown.csv` and `module_lightcyan1.csv`.
**Outputs**: under `result/13.module_gene_expression_correlation/` — `brown_lightcyan1_gene_expression_correlations.tsv` and `.filtered.tsv`, `correlation_threshold_summary.tsv`, `spearman_correlation_distribution.pdf|.png`, `pearson_correlation_distribution.pdf|.png` and `README.txt`.
**Run**: `Rscript module_gene_expression_correlation.R` (no arguments; all paths are resolved relative to the script's own directory, so the script must sit beside the two mean matrices and the `result/` tree and be run from that directory).
**Tools**: R 4.4.3 with ggplot2/scales (`ggsave(device = cairo_pdf)`) and base R (`cor.test(method = "spearman"/"pearson")`); settings as coded: correlations on log2(mean TPM + 1) over the 12 conditions, high-correlation threshold R ≥ 0.8, and a gene is analysed when it has mean TPM ≥ 0.5 in at least 2 conditions and both vectors have SD > 0.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details>
<summary><code>module_gene_expression_correlation.R</code> · R · 316 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/08_module_consistency/module_gene_expression_correlation.R)

```r
#!/usr/bin/env Rscript

# Compare per-gene expression patterns between Gifu and MG20 for the
# WGCNA brown and lightcyan1 modules.
#
# Filtering (on raw mean TPM):
#   1. Gifu OR MG20 has mean TPM >= 0.5 in at least 2 matched conditions.
#   2. Both Gifu and MG20 vectors have SD > 0.
# Correlations are calculated on log2(mean TPM + 1) across 12 matched
# organ/developmental conditions. Spearman is primary; Pearson is auxiliary.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(".", winslash = "/", mustWork = FALSE)
}
base_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)

gifu_file <- file.path(base_dir, "Gifu_all_samples_mean.tsv")
mg20_file <- file.path(base_dir, "MG20_all_samples_mean.tsv")
module_dir <- file.path(base_dir, "result", "06.modules_expr")
out_dir <- file.path(base_dir, "result", "13.module_gene_expression_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

r_threshold <- 0.8
tpm_threshold <- 0.5
min_conditions <- 2L

condition_order <- c(
  "flower", "leaf",
  "nodule_dpi10", "nodule_dpi21", "nodule_dpi45",
  "pod",
  "root_dpi10", "root_dpi21", "root_dpi45", "root_hpi48", "root_uni",
  "stem"
)

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing input file(s): ", paste(missing, collapse = ", "))
}

module_files <- c(
  brown = file.path(module_dir, "module_brown.csv"),
  lightcyan1 = file.path(module_dir, "module_lightcyan1.csv")
)
stop_if_missing(c(gifu_file, mg20_file, module_files))

message("Reading mean TPM matrices...")
gifu <- read.delim(gifu_file, check.names = FALSE, stringsAsFactors = FALSE)
mg20 <- read.delim(mg20_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_id" %in% names(gifu) || !"gene_id" %in% names(mg20)) {
  stop("Both mean TPM matrices must contain a gene_id column.")
}
if (anyDuplicated(gifu$gene_id) || anyDuplicated(mg20$gene_id)) {
  stop("Duplicated gene_id detected in a mean TPM matrix.")
}

rownames(gifu) <- gifu$gene_id
rownames(mg20) <- mg20$gene_id
gifu$gene_id <- NULL
mg20$gene_id <- NULL

colnames(gifu) <- sub("^Gifu_", "", colnames(gifu))
colnames(mg20) <- sub("^MG20_", "", colnames(mg20))

missing_gifu <- setdiff(condition_order, colnames(gifu))
missing_mg20 <- setdiff(condition_order, colnames(mg20))
if (length(missing_gifu) || length(missing_mg20)) {
  stop(
    "Matched conditions missing. Gifu: ", paste(missing_gifu, collapse = ", "),
    "; MG20: ", paste(missing_mg20, collapse = ", ")
  )
}

gifu <- as.matrix(gifu[, condition_order, drop = FALSE])
mg20 <- as.matrix(mg20[, condition_order, drop = FALSE])
storage.mode(gifu) <- "double"
storage.mode(mg20) <- "double"

read_module_genes <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  unique(as.character(x[[1]]))
}

module_genes <- lapply(module_files, read_module_genes)

safe_cor_test <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || sd(x) <= 0 || sd(y) <= 0) {
    return(c(r = NA_real_, p = NA_real_, n = length(x)))
  }
  z <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  c(r = unname(z$estimate), p = z$p.value, n = length(x))
}

analyse_module <- function(module_name, genes) {
  common_genes <- intersect(genes, intersect(rownames(gifu), rownames(mg20)))
  absent <- setdiff(genes, intersect(rownames(gifu), rownames(mg20)))
  if (!length(common_genes)) stop("No genes from module ", module_name, " found in both matrices.")

  g_raw <- gifu[common_genes, , drop = FALSE]
  m_raw <- mg20[common_genes, , drop = FALSE]
  n_g_ge <- rowSums(g_raw >= tpm_threshold, na.rm = TRUE)
  n_m_ge <- rowSums(m_raw >= tpm_threshold, na.rm = TRUE)
  sd_g <- apply(g_raw, 1, sd, na.rm = TRUE)
  sd_m <- apply(m_raw, 1, sd, na.rm = TRUE)
  n_complete <- rowSums(is.finite(g_raw) & is.finite(m_raw))

  pass_expression <- n_g_ge >= min_conditions | n_m_ge >= min_conditions
  pass_sd <- is.finite(sd_g) & is.finite(sd_m) & sd_g > 0 & sd_m > 0
  pass_complete <- n_complete >= 3L
  pass <- pass_expression & pass_sd & pass_complete

  status <- ifelse(
    !pass_expression, "low_expression",
    ifelse(!pass_sd, "zero_or_invalid_sd", ifelse(!pass_complete, "too_few_complete", "pass"))
  )

  result <- data.frame(
    gene_id = common_genes,
    module = module_name,
    n_Gifu_TPM_ge_0.5 = n_g_ge,
    n_MG20_TPM_ge_0.5 = n_m_ge,
    sd_Gifu_raw_TPM = sd_g,
    sd_MG20_raw_TPM = sd_m,
    n_complete_conditions = n_complete,
    filter_status = status,
    stringsAsFactors = FALSE
  )

  # Keep raw mean TPM vectors in a wide, auditable table.
  g_out <- as.data.frame(g_raw, check.names = FALSE)
  m_out <- as.data.frame(m_raw, check.names = FALSE)
  colnames(g_out) <- paste0("Gifu_", colnames(g_out), "_meanTPM")
  colnames(m_out) <- paste0("MG20_", colnames(m_out), "_meanTPM")
  result <- cbind(result, g_out, m_out)

  result$spearman_r <- NA_real_
  result$spearman_p <- NA_real_
  result$pearson_r <- NA_real_
  result$pearson_p <- NA_real_

  idx <- which(pass)
  g_log <- log2(g_raw + 1)
  m_log <- log2(m_raw + 1)
  for (i in idx) {
    sp <- safe_cor_test(g_log[i, ], m_log[i, ], "spearman")
    pe <- safe_cor_test(g_log[i, ], m_log[i, ], "pearson")
    result$spearman_r[i] <- sp["r"]
    result$spearman_p[i] <- sp["p"]
    result$pearson_r[i] <- pe["r"]
    result$pearson_p[i] <- pe["p"]
  }

  result$spearman_ge_0.8 <- !is.na(result$spearman_r) & result$spearman_r >= r_threshold
  result$pearson_ge_0.8 <- !is.na(result$pearson_r) & result$pearson_r >= r_threshold

  data.frame(
    module = module_name,
    module_genes = length(genes),
    genes_in_both_matrices = length(common_genes),
    genes_absent_from_matrix = length(absent),
    genes_pass_filter = sum(pass),
    genes_excluded = sum(!pass),
    spearman_R_ge_0.8_n = sum(result$spearman_ge_0.8),
    spearman_R_ge_0.8_pct = 100 * mean(result$spearman_ge_0.8[pass]),
    spearman_median_R = median(result$spearman_r[pass], na.rm = TRUE),
    pearson_R_ge_0.8_n = sum(result$pearson_ge_0.8),
    pearson_R_ge_0.8_pct = 100 * mean(result$pearson_ge_0.8[pass]),
    pearson_median_R = median(result$pearson_r[pass], na.rm = TRUE),
    stringsAsFactors = FALSE
  ) -> summary

  list(result = result, summary = summary)
}

analyses <- Map(analyse_module, names(module_genes), module_genes)
names(analyses) <- names(module_genes)
all_results <- do.call(rbind, lapply(analyses, `[[`, "result"))
rownames(all_results) <- NULL
summary_table <- do.call(rbind, lapply(analyses, `[[`, "summary"))
rownames(summary_table) <- NULL

passed_results <- subset(all_results, filter_status == "pass")
combined_summary <- data.frame(
  module = "brown+lightcyan1",
  module_genes = nrow(all_results),
  genes_in_both_matrices = nrow(all_results),
  genes_absent_from_matrix = 0L,
  genes_pass_filter = nrow(passed_results),
  genes_excluded = nrow(all_results) - nrow(passed_results),
  spearman_R_ge_0.8_n = sum(passed_results$spearman_r >= r_threshold, na.rm = TRUE),
  spearman_R_ge_0.8_pct = 100 * mean(passed_results$spearman_r >= r_threshold, na.rm = TRUE),
  spearman_median_R = median(passed_results$spearman_r, na.rm = TRUE),
  pearson_R_ge_0.8_n = sum(passed_results$pearson_r >= r_threshold, na.rm = TRUE),
  pearson_R_ge_0.8_pct = 100 * mean(passed_results$pearson_r >= r_threshold, na.rm = TRUE),
  pearson_median_R = median(passed_results$pearson_r, na.rm = TRUE),
  stringsAsFactors = FALSE
)
summary_table <- rbind(summary_table, combined_summary)

write.table(
  all_results,
  file.path(out_dir, "brown_lightcyan1_gene_expression_correlations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  subset(all_results, filter_status == "pass"),
  file.path(out_dir, "brown_lightcyan1_gene_expression_correlations.filtered.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  summary_table,
  file.path(out_dir, "correlation_threshold_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

plot_df <- subset(
  all_results,
  filter_status == "pass",
  select = c("gene_id", "module", "spearman_r", "pearson_r")
)

# Clean, publication-oriented distribution plot. The two modules are pooled
# for display, while the module column remains in the result table for tracing.
make_distribution_plot <- function(values, method, x_label, output_stem) {
  values <- values[is.finite(values)]
  n_high <- sum(values >= r_threshold)
  pct_high <- 100 * n_high / length(values)
  median_r <- median(values)
  annotation <- sprintf(
    "n = %s\nMedian = %.2f\nR >= %.1f: %s (%.1f%%)",
    format(length(values), big.mark = ","), median_r, r_threshold,
    format(n_high, big.mark = ","), pct_high
  )

  d <- data.frame(R = values)
  p <- ggplot(d, aes(x = R)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = 0.05, boundary = -1,
      fill = "#3C5488", color = "white", linewidth = 0.25, alpha = 0.88
    ) +
    geom_density(color = "#E64B35", linewidth = 1.15, adjust = 1) +
    geom_vline(
      xintercept = r_threshold, linetype = "22",
      color = "#B2182B", linewidth = 0.9
    ) +
    annotate(
      "label", x = -0.94, y = Inf, label = annotation,
      hjust = 0, vjust = 1.15, size = 4.7, lineheight = 1.08,
      label.size = 0.25, label.padding = unit(0.42, "lines"),
      color = "#252525", fill = scales::alpha("white", 0.92)
    ) +
    scale_x_continuous(
      limits = c(-1, 1), breaks = seq(-1, 1, 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(x = x_label, y = "Density") +
    theme_classic(base_size = 15) +
    theme(
      axis.title = element_text(size = 17, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(linewidth = 0.75, color = "black"),
      axis.ticks = element_line(linewidth = 0.65, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      plot.margin = margin(10, 14, 10, 10)
    )

  ggsave(file.path(out_dir, paste0(output_stem, ".pdf")), p,
         width = 7.2, height = 5.4, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(output_stem, ".png")), p,
         width = 7.2, height = 5.4, dpi = 400, bg = "white")
}

make_distribution_plot(
  plot_df$spearman_r, "Spearman", "Spearman correlation (rho)",
  "spearman_correlation_distribution"
)
make_distribution_plot(
  plot_df$pearson_r, "Pearson", "Pearson correlation (R)",
  "pearson_correlation_distribution"
)

# Remove the obsolete four-panel figure made by the previous script version.
unlink(file.path(out_dir, c("correlation_histogram_density.pdf",
                            "correlation_histogram_density.png")))

readme <- c(
  "Per-gene Gifu-MG20 expression correlation for WGCNA brown and lightcyan1 modules",
  "",
  paste0("Matched conditions (n=", length(condition_order), "): ", paste(condition_order, collapse = ", ")),
  "Filtering uses raw mean TPM:",
  paste0("  - Gifu OR MG20 has TPM >= ", tpm_threshold, " in at least ", min_conditions, " conditions"),
  "  - both genotype vectors have SD > 0",
  "Correlations use log2(mean TPM + 1). Spearman is primary; Pearson is auxiliary.",
  paste0("High-correlation threshold: R >= ", r_threshold, " (positive correlation only)."),
  "Gifu_root_dpi02 is not used because no matched MG20 condition exists.",
  "",
  "Files:",
  "  brown_lightcyan1_gene_expression_correlations.tsv: all module genes with TPM vectors and filter status",
  "  brown_lightcyan1_gene_expression_correlations.filtered.tsv: genes passing filters",
  "  correlation_threshold_summary.tsv: separate-module and pooled counts, percentages, and median R",
  "  spearman_correlation_distribution.pdf/png: pooled brown + lightcyan1 Spearman distribution",
  "  pearson_correlation_distribution.pdf/png: pooled brown + lightcyan1 Pearson distribution"
)
writeLines(readme, file.path(out_dir, "README.txt"))

message("Analysis complete. Output directory: ", out_dir)
print(summary_table, row.names = FALSE)
```

</details>
