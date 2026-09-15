## =============================================================
## Sensitivity check: re-run HC vs LTBI DEG analysis for GSE99374
## WITHOUT the two apparent outlier samples (TU0021_CD8_LTBI,
## TP0015_CD8_LTBI) and compare against the original 4 significant genes.
## Standalone — rebuilds everything from the saved filtered counts.
## =============================================================

setwd("/Users/negin/Desktop/TB_DEG_pipeline/scripts")   # <-- edit this
source("00_functions.R")

BASE_DIR   <- file.path("..", "GSE99374")
OUTLIERS <- c("TU0021_CD8_LTBI", "TP0001_CD8_LTBI")   # samples to drop
LFC_TH     <- 1
PADJ_TH    <- 0.05

## Where the sensitivity-check outputs will go (separate from the
## original results, so nothing you already have gets overwritten)
SENS_DIR <- file.path(BASE_DIR, "sensitivity_no_outliers")
for (d in c("results", "figures")) dir.create(file.path(SENS_DIR, d), recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load what you already produced, drop the outliers ----
filt_counts <- read.csv(file.path(BASE_DIR, "03_filtered_counts", "GSE99374_filtered_count_matrix.csv"),
                         row.names = 1, check.names = FALSE)
metadata <- read.csv(file.path(BASE_DIR, "02_metadata", "GSE99374_sample_metadata.csv"))
rownames(metadata) <- metadata$sample_id

keep_samples <- setdiff(colnames(filt_counts), OUTLIERS)
filt_counts  <- filt_counts[, keep_samples]
metadata     <- metadata[keep_samples, ]
metadata$group <- factor(metadata$group, levels = c("HC", "LTBI"))

message("Samples remaining: ", ncol(filt_counts),
        " (", sum(metadata$group == "HC"), " HC / ",
        sum(metadata$group == "LTBI"), " LTBI)")

## ---- 2. Re-run DESeq2 from the filtered counts (no need to re-filter) ----
dds <- DESeqDataSetFromMatrix(filt_counts, metadata, design = ~ group)
dds <- DESeq(dds)
vst_mat <- run_vst(dds)

## ---- 3. PCA — confirm the outliers are actually gone / groups look different now ----
plot_pca(vst_mat, metadata, out_png = file.path(SENS_DIR, "figures", "GSE99374_PCA_no_outliers.png"),
          title = "GSE99374 PCA (outliers removed)")

## ---- 4. Re-run the HC vs LTBI DEG test ----
res_new <- run_deg(dds, contrast = c("group", "LTBI", "HC"),
                    out_csv = file.path(SENS_DIR, "results", "GSE99374_HC_vs_LTBI_DEG_no_outliers.csv"),
                    lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

volcano_plot(res_new, out_png = file.path(SENS_DIR, "figures", "GSE99374_volcano_no_outliers.png"),
             title = "GSE99374: HC vs LTBI (outliers removed)")

sig_new <- get_significant_degs(res_new,
                                 file.path(SENS_DIR, "results", "GSE99374_sigDEGs_no_outliers.csv"),
                                 lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

## ---- 5. THE ANSWER: compare against your original 4 genes ----
original_genes <- read.csv(file.path(BASE_DIR, "07_significant_degs", "GSE99374_HC_vs_LTBI_sigDEGs.csv"))$gene
new_genes      <- sig_new$gene

cat("\n================ SENSITIVITY CHECK RESULT ================\n")
cat("Original significant genes (with outliers):   ", paste(original_genes, collapse = ", "), "\n")
cat("New significant genes (outliers removed):      ", paste(new_genes, collapse = ", "), "\n")
cat("Still significant after removing outliers:     ", paste(intersect(original_genes, new_genes), collapse = ", "), "\n")
cat("Lost after removing outliers (were outlier-driven): ", paste(setdiff(original_genes, new_genes), collapse = ", "), "\n")
cat("New genes that only appear once outliers removed:   ", paste(setdiff(new_genes, original_genes), collapse = ", "), "\n")
cat("=============================================================\n")
