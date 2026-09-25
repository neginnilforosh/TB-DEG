## =============================================================
## GSE229020 — miRNA (NanoString nCounter) analysis  [FIXED]
## HC vs LTB (Latent), HC vs DS-TB (Drug-Susceptible), HC vs DR-TB (Drug-Resistant)
## =============================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(UpSetR)
  library(WGCNA)
})

ACC      <- "GSE229020"
BASE_DIR <- file.path("..", ACC)
LFC_TH   <- 1
PADJ_TH  <- 0.05

dir_raw   <- file.path(BASE_DIR, "01_raw_data")
dir_meta  <- file.path(BASE_DIR, "02_metadata")
dir_filt  <- file.path(BASE_DIR, "03_filtered")
dir_norm  <- file.path(BASE_DIR, "04_normalized")
dir_corr  <- file.path(BASE_DIR, "05_sample_correlation")
dir_deg   <- file.path(BASE_DIR, "06_deg_results")
dir_sig   <- file.path(BASE_DIR, "07_significant_degs")
dir_share <- file.path(BASE_DIR, "08_shared_unique_degs")
dir_wgcna <- file.path(BASE_DIR, "09_wgcna_input")
dir_fig   <- file.path(BASE_DIR, "figures")

for (d in c(dir_raw, dir_meta, dir_filt, dir_norm, dir_corr,
            dir_deg, dir_sig, dir_share, dir_wgcna, dir_fig)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- 1. Load the deposited normalized miRNA matrix ----
candidates <- list.files(dir_raw, pattern = "normalized_data", full.names = TRUE, recursive = TRUE)
candidates <- candidates[!grepl("\\.gz$", candidates)]   # prefer the already-unzipped .txt if both exist
cat("Found candidate file(s):\n"); print(candidates)
stopifnot(length(candidates) >= 1)
raw_path <- candidates[1]
full_tbl <- read.delim(raw_path, check.names = FALSE)
full_tbl <- full_tbl[, !grepl("^(V\\d+|X\\.)$|^$", names(full_tbl))]  # drop stray blank/unnamed cols

cat("\nRows before row-level cleanup:", nrow(full_tbl), "\n")
cat("Class Name breakdown (blank = the spreadsheet-export junk rows we're about to drop):\n")
print(table(full_tbl$`Class Name`, useNA = "ifany"))


full_tbl <- full_tbl[!is.na(full_tbl$`Class Name`) &
                      full_tbl$`Class Name` == "Endogenous1" &
                      !is.na(full_tbl$ID_REF) &
                      trimws(full_tbl$ID_REF) != "", ]
cat("Rows after row-level cleanup:", nrow(full_tbl), "(should be 798 real miRNAs, not 828)\n\n")

sample_cols <- setdiff(names(full_tbl),
                       c("ID_REF", "Probe Name", "miRNA_ID", "Accession", "Class Name"))
expr <- as.matrix(full_tbl[, sample_cols])
rownames(expr) <- full_tbl$ID_REF
storage.mode(expr) <- "numeric"
cat("Expression matrix:", nrow(expr), "miRNAs x", ncol(expr), "samples\n")
cat("Samples:", paste(colnames(expr), collapse = ", "), "\n")
cat("Any NA left in expr after cleanup?", anyNA(expr), "(should be FALSE)\n")

## ---- 2. Metadata -- parsed straight from the sample-name prefixes ----
metadata <- data.frame(
  sample_id = colnames(expr),
  group = dplyr::case_when(
    grepl("^HC",    colnames(expr)) ~ "HC",
    grepl("^LTB",   colnames(expr)) ~ "Latent",
    grepl("^DS-TB", colnames(expr)) ~ "DrugSusceptible",
    grepl("^DR-TB", colnames(expr)) ~ "DrugResistant",
    TRUE ~ NA_character_
  )
)
stopifnot(!any(is.na(metadata$group)))
metadata$group <- factor(metadata$group,
                         levels = c("HC", "Latent", "DrugSusceptible", "DrugResistant"))
rownames(metadata) <- metadata$sample_id
print(table(metadata$group))

write.csv(metadata, file.path(dir_meta, paste0(ACC, "_miRNA_sample_metadata.csv")), row.names = FALSE)

design_plot <- ggplot(metadata, aes(group, fill = group)) + geom_bar() +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = paste(ACC, "miRNA — sample design"), x = NULL, y = "N samples")
ggsave(file.path(dir_fig, paste0(ACC, "_miRNA_sample_design.png")), design_plot,
      width = 5, height = 4, dpi = 300)

## ---- 3. Light filtering: drop miRNAs with ~zero signal/variance across ALL samples ----
keep <- apply(expr, 1, function(x) sd(x) > 0) & rowMeans(expr) > 1
expr_filt <- expr[keep, ]
write.csv(expr_filt, file.path(dir_filt, paste0(ACC, "_filtered_miRNA_matrix.csv")))

filt_summary <- data.frame(stage = c("Before filtering", "After filtering"),
                           n = c(nrow(expr), nrow(expr_filt)))
filt_plot <- ggplot(filt_summary, aes(stage, n, fill = stage)) + geom_col() +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = paste(ACC, "miRNA filtering"), x = NULL, y = "Number of miRNAs")
ggsave(file.path(dir_fig, paste0(ACC, "_miRNA_filtering_summary.png")), filt_plot,
      width = 5, height = 4, dpi = 300)
cat("\nmiRNAs:", nrow(expr), "before filtering,", nrow(expr_filt), "after\n")

## ---- 4. log2 transform (standard for NanoString/microarray-style intensities) ----
log2_mat <- log2(expr_filt + 1)
cat("Any NaN in log2_mat?", anyNA(log2_mat), "(should be FALSE now)\n")
write.csv(log2_mat, file.path(dir_norm, paste0(ACC, "_log2_normalized_matrix.csv")))

pca <- prcomp(t(log2_mat), scale. = TRUE)
pct_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))[1:2]
pca_df <- data.frame(pca$x[, 1:2], sample_id = rownames(pca$x))
pca_df <- merge(pca_df, metadata, by = "sample_id")
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = group, label = sample_id)) +
  geom_point(size = 3) + geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = paste(ACC, "miRNA PCA"),
      x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_minimal(base_size = 12)
ggsave(file.path(dir_fig, paste0(ACC, "_miRNA_PCA.png")), pca_plot, width = 6, height = 5, dpi = 300)

## ---- 5. Sample correlation ----
sample_cor <- cor(log2_mat, method = "pearson", use = "pairwise.complete.obs")
write.csv(sample_cor, file.path(dir_corr, paste0(ACC, "_miRNA_sample_correlation_matrix.csv")))
ann <- data.frame(Group = metadata$group); rownames(ann) <- metadata$sample_id
png(file.path(dir_fig, paste0(ACC, "_miRNA_sample_correlation_heatmap.png")),
   width = 7, height = 6, units = "in", res = 300)
pheatmap(sample_cor, annotation_col = ann, annotation_row = ann,
        main = paste(ACC, "miRNA sample correlation"))
dev.off()

## ---- Outlier check  ----
avg_corr <- sort(rowMeans(sample_cor))
cat("\n--- Lowest mean sample correlation (possible outliers) ---\n")
print(head(avg_corr, 5))

## ---- 6/7. Differential expression via limma, one comparison per non-HC group ----
group <- metadata$group
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2_mat, design)

comparisons <- list(
  HC_vs_Latent          = "Latent - HC",
  HC_vs_DrugSusceptible = "DrugSusceptible - HC",
  HC_vs_DrugResistant   = "DrugResistant - HC"
)

all_sig <- list()
for (cmp_name in names(comparisons)) {
  contrast <- makeContrasts(comparisons[[cmp_name]], levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  res <- topTable(fit2, number = Inf, sort.by = "none")
  res$miRNA <- rownames(res)
  res <- res[, c("miRNA", "logFC", "P.Value", "adj.P.Val", "AveExpr", "t")]
  write.csv(res, file.path(dir_deg, paste0(ACC, "_", cmp_name, "_DE.csv")), row.names = FALSE)

  res$sig <- ifelse(res$adj.P.Val < PADJ_TH & abs(res$logFC) >= LFC_TH, "Significant", "Not significant")
  volcano <- ggplot(res, aes(logFC, -log10(P.Value), color = sig)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Significant = "firebrick", `Not significant` = "grey70")) +
    geom_vline(xintercept = c(-LFC_TH, LFC_TH), linetype = "dashed", color = "grey40") +
    labs(title = paste(ACC, "miRNA —", gsub("_", " ", cmp_name)),
        x = "log2 fold change", y = "-log10(p-value)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(dir_fig, paste0(ACC, "_", cmp_name, "_volcano.png")), volcano, width = 6, height = 5, dpi = 300)

  sig_df <- res[res$sig == "Significant", ]
  sig_df <- sig_df[order(sig_df$adj.P.Val), ]
  write.csv(sig_df, file.path(dir_sig, paste0(ACC, "_", cmp_name, "_sig_miRNAs.csv")), row.names = FALSE)
  all_sig[[cmp_name]] <- sig_df$miRNA

  if (nrow(sig_df) > 1) {
    top_n <- min(50, nrow(sig_df))
    png(file.path(dir_fig, paste0(ACC, "_", cmp_name, "_heatmap.png")),
       width = 7, height = 8, units = "in", res = 300)
    pheatmap(log2_mat[sig_df$miRNA[1:top_n], ], scale = "row", annotation_col = ann,
            main = paste(ACC, "top", top_n, "miRNAs,", gsub("_", " ", cmp_name)))
    dev.off()
  }
  cat(sprintf("%-25s : %d significant miRNAs\n", cmp_name, nrow(sig_df)))
}

## ---- 8. Shared/unique significant miRNAs across the three comparisons ----
non_empty <- all_sig[sapply(all_sig, length) > 0]
if (length(non_empty) >= 2) {
  all_mirnas <- unique(unlist(non_empty))
  membership <- as.data.frame(sapply(non_empty, function(g) as.integer(all_mirnas %in% g)))
  rownames(membership) <- all_mirnas
  write.csv(cbind(miRNA = rownames(membership), membership),
           file.path(dir_share, paste0(ACC, "_shared_unique_miRNAs.csv")), row.names = FALSE)

  png(file.path(dir_fig, paste0(ACC, "_miRNA_upset.png")), width = 8, height = 5, units = "in", res = 300)
  upset(membership, nsets = ncol(membership), order.by = "freq",
       main.bar.color = "steelblue", sets.bar.color = "darkorange")
  dev.off()
} else {
  cat("\nNOTE: fewer than 2 comparisons had any significant miRNAs -- skipping shared/unique + UpSet.\n")
}

## ---- 9. WGCNA-ready exports (full filtered matrix, not just significant miRNAs) ----
write.csv(t(log2_mat), file.path(dir_wgcna, paste0(ACC, "_WGCNA_miRNA_expression_matrix.csv")))
traits <- data.frame(
  HC              = as.integer(metadata$group == "HC"),
  Latent          = as.integer(metadata$group == "Latent"),
  DrugSusceptible = as.integer(metadata$group == "DrugSusceptible"),
  DrugResistant   = as.integer(metadata$group == "DrugResistant"),
  row.names = metadata$sample_id
)
write.csv(traits, file.path(dir_wgcna, paste0(ACC, "_WGCNA_miRNA_trait_file.csv")))

## ---- Summary ----
cat("\n========== miRNA DE COUNT SUMMARY ==========\n")
for (cmp_name in names(all_sig)) {
  cat(sprintf("%-25s : %d significant miRNAs\n", cmp_name, length(all_sig[[cmp_name]])))
}
cat("(threshold: adj.P.Val < ", PADJ_TH, ", |log2FC| >= ", LFC_TH, ")\n", sep = "")


