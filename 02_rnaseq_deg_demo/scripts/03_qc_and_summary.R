#!/usr/bin/env Rscript

# Generate quality-control plots and a concise Markdown summary for the
# teaching RNA-seq differential-expression demo.

# ----------------------------
# 0. Locate project directories
# ----------------------------
command_args <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_args, value = TRUE)

if (length(file_argument) != 1) {
  stop("Cannot determine the script path. Run with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_argument))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1. Check required packages
# ----------------------------
required_packages <- c("ggplot2", "pheatmap")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: Rscript setup/install_r_packages.R from the repository root."
  )
}

# ----------------------------
# 2. Read input and result tables
# ----------------------------
counts <- read.csv(
  file.path(data_dir, "counts.csv"),
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  file.path(data_dir, "metadata.csv"),
  stringsAsFactors = FALSE
)

required_result_files <- c(
  "all_deseq2_results.csv",
  "deg_results.csv",
  "enrichment_results.csv",
  "normalized_counts.csv"
)

missing_result_files <- required_result_files[
  !file.exists(file.path(results_dir, required_result_files))
]

if (length(missing_result_files) > 0) {
  stop(
    "Missing result files: ",
    paste(missing_result_files, collapse = ", "),
    "\nRun scripts/01_deseq2_analysis.R and scripts/02_enrichment_demo.R first."
  )
}

metadata <- metadata[match(colnames(counts), metadata$sample), ]
rownames(metadata) <- metadata$sample

if (!identical(colnames(counts), rownames(metadata))) {
  stop("The count matrix and metadata could not be aligned.")
}

all_results <- read.csv(
  file.path(results_dir, "all_deseq2_results.csv"),
  stringsAsFactors = FALSE
)

deg_results <- read.csv(
  file.path(results_dir, "deg_results.csv"),
  stringsAsFactors = FALSE
)

enrichment_results <- read.csv(
  file.path(results_dir, "enrichment_results.csv"),
  stringsAsFactors = FALSE
)

normalized_counts <- read.csv(
  file.path(results_dir, "normalized_counts.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(normalized_counts) <- normalized_counts$gene_id
normalized_counts <- normalized_counts[, setdiff(colnames(normalized_counts), "gene_id"), drop = FALSE]
normalized_counts <- normalized_counts[, rownames(metadata), drop = FALSE]

# ----------------------------
# 3. Library-size plot
# ----------------------------
library_size_df <- data.frame(
  sample = colnames(counts),
  condition = metadata[colnames(counts), "condition"],
  total_counts = as.numeric(colSums(counts)),
  stringsAsFactors = FALSE
)

library_size_df$sample <- factor(library_size_df$sample, levels = library_size_df$sample)

library_size_plot <- ggplot2::ggplot(
  library_size_df,
  ggplot2::aes(x = sample, y = total_counts, fill = condition)
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = total_counts),
    vjust = -0.35,
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(
    values = c("control" = "#F8766D", "stress" = "#00BFC4")
  ) +
  ggplot2::labs(
    title = "RNA-seq library size check",
    subtitle = "Total simulated read counts per sample before normalization",
    x = "Sample",
    y = "Total count",
    fill = "Condition"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(figures_dir, "library_size_barplot.png"),
  library_size_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ----------------------------
# 4. Sample-correlation heatmap
# ----------------------------
log_normalized_counts <- log2(as.matrix(normalized_counts) + 1)
sample_correlation <- cor(log_normalized_counts, method = "pearson")

write.csv(
  sample_correlation,
  file.path(results_dir, "sample_correlation_matrix.csv")
)

annotation_col <- data.frame(condition = metadata$condition)
rownames(annotation_col) <- rownames(metadata)
annotation_colors <- list(
  condition = c("control" = "#F8766D", "stress" = "#00BFC4")
)

sample_pairs <- combn(colnames(sample_correlation), 2)
pair_correlations <- data.frame(
  sample_a = sample_pairs[1, ],
  sample_b = sample_pairs[2, ],
  correlation = mapply(
    function(a, b) sample_correlation[a, b],
    sample_pairs[1, ],
    sample_pairs[2, ]
  ),
  stringsAsFactors = FALSE
)
pair_correlations$comparison <- ifelse(
  metadata[pair_correlations$sample_a, "condition"] ==
    metadata[pair_correlations$sample_b, "condition"],
  "within_condition",
  "between_condition"
)

mean_within_condition_correlation <- mean(
  pair_correlations$correlation[pair_correlations$comparison == "within_condition"]
)

mean_between_condition_correlation <- mean(
  pair_correlations$correlation[pair_correlations$comparison == "between_condition"]
)

png(
  file.path(figures_dir, "sample_correlation_heatmap.png"),
  width = 1800,
  height = 1500,
  res = 220
)
pheatmap::pheatmap(
  sample_correlation,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample-to-sample correlation",
  border_color = NA,
  fontsize_number = 9
)
dev.off()


# ----------------------------
# 5. Review-friendly tables
# ----------------------------
ordered_results <- all_results[
  order(all_results$padj, -abs(all_results$log2FoldChange), na.last = TRUE),
]

top_degs_for_review <- head(
  ordered_results[ordered_results$status %in% c("Up", "Down"), ],
  15
)

write.csv(
  top_degs_for_review,
  file.path(results_dir, "top_degs_for_review.csv"),
  row.names = FALSE
)

# ----------------------------
# 6. Markdown summary
# ----------------------------
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(signif(x, digits), scientific = FALSE, trim = TRUE))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", format(signif(x, 3), scientific = TRUE, trim = TRUE))
}

make_gene_table <- function(df, max_rows = 8) {
  if (nrow(df) == 0) {
    return("_No genes passed the current threshold._")
  }

  df <- head(df, max_rows)
  rows <- sprintf(
    "| %s | %s | %s | %s |",
    df$gene_id,
    fmt_num(df$log2FoldChange, 3),
    fmt_p(df$padj),
    df$status
  )

  paste(
    c(
      "| Gene | log2FC | adjusted p-value | Status |",
      "|---|---:|---:|---|",
      rows
    ),
    collapse = "\n"
  )
}

make_enrichment_table <- function(df, max_rows = 5) {
  if (nrow(df) == 0) {
    return("_No enrichment result is available._")
  }

  df <- head(df, max_rows)
  rows <- sprintf(
    "| %s | %s | %s | %s |",
    df$gene_set,
    df$overlap_count,
    df$set_size,
    fmt_p(df$padj)
  )

  paste(
    c(
      "| Gene set | DEG overlap | Set size | adjusted p-value |",
      "|---|---:|---:|---:|",
      rows
    ),
    collapse = "\n"
  )
}

top_up <- ordered_results[
  ordered_results$status == "Up" & !is.na(ordered_results$padj),
]

top_down <- ordered_results[
  ordered_results$status == "Down" & !is.na(ordered_results$padj),
]

significant_enrichment <- enrichment_results[
  enrichment_results$padj < 0.05,
]

summary_lines <- c(
  "# RNA-seq Demo Analysis Summary",
  "",
  "## Scope",
  "",
  "This summary reports the reproducible outputs from a small simulated RNA-seq count matrix. The dataset is designed for workflow practice, not for making biological claims.",
  "",
  "## Input design",
  "",
  sprintf("- Genes tested: %s", nrow(all_results)),
  sprintf("- Samples: %s control and %s stress replicates", sum(metadata$condition == "control"), sum(metadata$condition == "stress")),
  "- Model: DESeq2 design `~ condition`, contrast `stress / control`",
  "- DEG threshold: adjusted p-value < 0.05 and absolute log2 fold change >= 1",
  "",
  "## Key results",
  "",
  sprintf("- Total teaching DEGs: %s", nrow(deg_results)),
  sprintf("- Up-regulated under stress: %s", sum(deg_results$status == "Up")),
  sprintf("- Down-regulated under stress: %s", sum(deg_results$status == "Down")),
  sprintf("- Significant simulated gene sets: %s", nrow(significant_enrichment)),
  sprintf("- Mean within-condition sample correlation: %s", fmt_num(mean_within_condition_correlation, 3)),
  sprintf("- Mean between-condition sample correlation: %s", fmt_num(mean_between_condition_correlation, 3)),
  "",
  "## Top up-regulated genes",
  "",
  make_gene_table(top_up, max_rows = 8),
  "",
  "## Top down-regulated genes",
  "",
  make_gene_table(top_down, max_rows = 8),
  "",
  "## Gene-set enrichment",
  "",
  make_enrichment_table(enrichment_results, max_rows = 5),
  "",
  "## Interpretation",
  "",
  "- PCA separates control and stress samples mainly along PC1, showing that the simulated treatment effect dominates sample variation.",
  "- The volcano plot shows a balanced set of strongly up- and down-regulated teaching genes.",
  "- The top-gene heatmap separates the two conditions and demonstrates how expression patterns can be reviewed across biological replicates.",
  "- The enrichment demo recovers the expected simulated stress-response and photosynthesis gene sets, illustrating over-representation analysis logic.",
  "",
  "## Quality-control checks added",
  "",
  "- `library_size_barplot.png` checks whether samples have comparable total counts before normalization.",
  "- `sample_correlation_heatmap.png` checks whether biological replicates are more similar to each other than to the other condition.",
  "",
  "## Limitations and next steps",
  "",
  "- Replace simulated genes with a public plant RNA-seq dataset from GEO or SRA.",
  "- Add raw-read QC with FastQC/MultiQC and transcript quantification with Salmon or STAR + featureCounts.",
  "- Use real gene identifiers and organism-specific GO/KEGG annotation.",
  "- Connect candidate genes back to literature evidence before making biological conclusions."
)

writeLines(
  summary_lines,
  con = file.path(results_dir, "analysis_summary.md"),
  useBytes = TRUE
)

cat("QC plots and analysis summary finished.\n")
