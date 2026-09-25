## =============================================================
## 00_functions.R
## Shared helper functions for the TB DEG pipeline
## (GSE222001, GSE161829, GSE99374, GSE229020)
##
## Requires: DESeq2, GEOquery, pheatmap, ggplot2, ggrepel,
##           RColorBrewer, UpSetR, matrixStats
##
## install.packages/BiocManager calls are left commented out --
##
## =============================================================

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("DESeq2","GEOquery","pheatmap","apeglm"))
#install.packages(c("ggplot2","ggrepel","RColorBrewer","UpSetR","matrixStats","pheatmap"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(UpSetR)
  library(matrixStats)
})

## -------------------------------------------------------------
## 1. FILTERING
##    Removes very low-expression genes.
##    Rule: keep genes with >= min_count reads in >= min_samples
##    samples (min_samples defaults to size of the smallest group).
## -------------------------------------------------------------
filter_low_expression <- function(counts, min_count = 10, min_samples = 3) {
  keep <- rowSums(counts >= min_count) >= min_samples
  list(filtered = counts[keep, , drop = FALSE],
       n_before = nrow(counts),
       n_after  = sum(keep))
}

plot_filtering_summary <- function(n_before, n_after, out_png, title = "Gene filtering summary") {
  df <- data.frame(stage = factor(c("Before filtering", "After filtering"),
                                   levels = c("Before filtering", "After filtering")),
                    n_genes = c(n_before, n_after))
  p <- ggplot(df, aes(x = stage, y = n_genes, fill = stage)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n_genes), vjust = -0.4, size = 4) +
    scale_fill_manual(values = c("#a6bddb", "#2166ac")) +
    labs(title = title, x = NULL, y = "Number of genes") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(out_png, p, width = 5, height = 4.5, dpi = 300)
  p
}

## -------------------------------------------------------------
## 2. NORMALIZATION (VST) + PCA
## -------------------------------------------------------------
run_vst <- function(dds) {

  vst_obj <- vst(dds, blind = TRUE)
  assay(vst_obj)
}

plot_pca <- function(vst_mat, metadata, group_col = "group", out_png,
                      title = "PCA of VST-normalized expression") {
  pca <- prcomp(t(vst_mat), scale. = FALSE)
  var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    sample = rownames(pca$x),
                    group = metadata[[group_col]][match(rownames(pca$x), rownames(metadata))])
  p <- ggplot(df, aes(PC1, PC2, color = group, label = sample)) +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
    labs(title = title,
         x = paste0("PC1 (", var_expl[1], "%)"),
         y = paste0("PC2 (", var_expl[2], "%)"),
         color = "Group") +
    theme_bw(base_size = 12)
  ggsave(out_png, p, width = 6.5, height = 5.5, dpi = 300)
  p
}

## -------------------------------------------------------------
## 3. SAMPLE CORRELATION HEATMAP
## -------------------------------------------------------------
sample_correlation <- function(vst_mat, metadata, group_col = "group",
                                out_matrix_csv, out_png,
                                title = "Sample-to-sample correlation") {
  cor_mat <- cor(vst_mat, method = "pearson")
  write.csv(cor_mat, out_matrix_csv)

  ann <- data.frame(Group = metadata[[group_col]])
  rownames(ann) <- rownames(metadata)

  png(out_png, width = 7, height = 6.5, units = "in", res = 300)
  pheatmap(cor_mat,
           annotation_col = ann,
           annotation_row = ann,
           main = title,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
  dev.off()
  cor_mat
}

## -------------------------------------------------------------
## 4. DESeq2 DEG CALLING + VOLCANO PLOT
##    contrast = c("group", "TreatmentLevel", "ReferenceLevel")
## -------------------------------------------------------------
run_deg <- function(dds, contrast, out_csv,
                     lfc_thresh = 1, padj_thresh = 0.05, shrink = TRUE) {
  res <- results(dds, contrast = contrast, alpha = padj_thresh)
  if (shrink) {
    coef_name <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
    if (coef_name %in% resultsNames(dds)) {
      res <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
    }
  }
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
  res_df$significant <- with(res_df,
                              !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) >= lfc_thresh)
  res_df <- res_df[order(res_df$padj), ]
  write.csv(res_df, out_csv, row.names = FALSE)
  res_df
}

volcano_plot <- function(res_df, out_png, title,
                          lfc_thresh = 1, padj_thresh = 0.05, top_label_n = 15) {
  df <- res_df
  df$status <- "NS"
  df$status[df$log2FoldChange >=  lfc_thresh & df$padj < padj_thresh] <- "Up"
  df$status[df$log2FoldChange <= -lfc_thresh & df$padj < padj_thresh] <- "Down"
  df$status <- factor(df$status, levels = c("Down", "NS", "Up"))
  df$negLog10Padj <- -log10(df$padj)

  top_genes <- df[df$status != "NS", ]
  top_genes <- top_genes[order(top_genes$padj), ][seq_len(min(top_label_n, nrow(top_genes))), ]

  p <- ggplot(df, aes(log2FoldChange, negLog10Padj, color = status)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c(Down = "#2166ac", NS = "grey80", Up = "#b2182b")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 3,
                              max.overlaps = 30, show.legend = FALSE, color = "black") +
    labs(title = title, x = "log2 Fold Change", y = "-log10 adjusted p-value", color = NULL) +
    theme_bw(base_size = 12)
  ggsave(out_png, p, width = 6.5, height = 5.5, dpi = 300)
  p
}

## -------------------------------------------------------------
## 5. SIGNIFICANT DEG LIST + HEATMAP
## -------------------------------------------------------------
get_significant_degs <- function(res_df, out_csv, lfc_thresh = 1, padj_thresh = 0.05) {
  sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thresh & abs(res_df$log2FoldChange) >= lfc_thresh, ]
  sig <- sig[order(sig$padj), ]
  write.csv(sig, out_csv, row.names = FALSE)
  sig
}

sig_deg_heatmap <- function(vst_mat, sig_genes, metadata, group_col = "group",
                             out_png, title = "Top significant DEGs", top_n = 50) {
  genes_use <- intersect(sig_genes, rownames(vst_mat))
  genes_use <- head(genes_use, top_n)
  mat <- vst_mat[genes_use, , drop = FALSE]
  mat_scaled <- t(scale(t(mat)))  # z-score per gene

  ann <- data.frame(Group = metadata[[group_col]])
  rownames(ann) <- rownames(metadata)

  png(out_png, width = 7.5, height = 8.5, units = "in", res = 300)
  pheatmap(mat_scaled,
           annotation_col = ann,
           show_rownames = TRUE, show_colnames = TRUE,
           main = title,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           fontsize_row = 6)
  dev.off()
}

## -------------------------------------------------------------
## 6. SHARED / UNIQUE DEGs ACROSS COMPARISONS (UpSet plot)
##    deg_lists = named list of character vectors of significant gene IDs
## -------------------------------------------------------------
shared_unique_degs <- function(deg_lists, out_csv, out_png,
                                title = "Shared and unique DEGs") {
  all_genes <- unique(unlist(deg_lists))
  membership <- sapply(deg_lists, function(x) all_genes %in% x)
  rownames(membership) <- all_genes
  out_df <- data.frame(gene = all_genes, membership * 1L)
  write.csv(out_df, out_csv, row.names = FALSE)

  upset_input <- as.data.frame(membership * 1L)
  png(out_png, width = 8, height = 5.5, units = "in", res = 300)
  print(UpSetR::upset(upset_input, sets = colnames(upset_input),
                       order.by = "freq", main.bar.color = "#2166ac",
                       sets.bar.color = "#b2182b", text.scale = 1.2))
  dev.off()
  out_df
}

## -------------------------------------------------------------
## 7. WGCNA-READY EXPORTS
##    (full VST matrix -- NOT restricted to significant DEGs -- + trait file)
## -------------------------------------------------------------
export_wgcna_expression <- function(vst_mat, out_csv) {
  write.csv(vst_mat, out_csv)
}

export_wgcna_traits <- function(metadata, trait_cols, out_csv) {
  # trait_cols: named list mapping trait name -> 0/1 (or numeric) vector
  # Example: list(TB_only = c(1,1,0,0), Comorbidity_only = c(0,0,1,1))
  traits <- as.data.frame(trait_cols)
  rownames(traits) <- rownames(metadata)
  write.csv(traits, out_csv)
}

## -------------------------------------------------------------
## 8. SAMPLE-DESIGN SCHEMATIC (simple bar of n per group)
## -------------------------------------------------------------
plot_sample_design <- function(metadata, group_col = "group", out_png,
                                title = "Sample design") {
  df <- as.data.frame(table(metadata[[group_col]]))
  colnames(df) <- c("group", "n")
  p <- ggplot(df, aes(x = group, y = n, fill = group)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4) +
    labs(title = title, x = NULL, y = "Number of samples") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(out_png, p, width = 5, height = 4.5, dpi = 300)
  p
}
