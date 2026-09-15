## =============================================================
## GSE161829 -- TBneg vs LTBI, TBneg vs ATB
## Comparisons: TBneg_vs_LTBI, TBneg_vs_ATB
##
## RUN THIS ON A MACHINE WITH INTERNET ACCESS TO NCBI GEO.
## Structure mirrors GSE222001_run.R -- see that file for inline
## explanations of every step.
## =============================================================
source("00_functions.R")
library(GEOquery)

ACC       <- "GSE161829"
BASE_DIR  <- file.path("..", ACC)
LFC_TH    <- 1
PADJ_TH   <- 0.05

dir_raw   <- file.path(BASE_DIR, "01_raw_counts")
dir_meta  <- file.path(BASE_DIR, "02_metadata")
dir_filt  <- file.path(BASE_DIR, "03_filtered_counts")
dir_norm  <- file.path(BASE_DIR, "04_normalized")
dir_corr  <- file.path(BASE_DIR, "05_sample_correlation")
dir_deg   <- file.path(BASE_DIR, "06_deg_results")
dir_sig   <- file.path(BASE_DIR, "07_significant_degs")
dir_share <- file.path(BASE_DIR, "08_shared_unique_degs")
dir_wgcna <- file.path(BASE_DIR, "09_wgcna_input")
dir_fig   <- file.path(BASE_DIR, "figures")

## STEP 1. Download ---------------------------------------------------
gse <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)
getGEOSuppFiles(ACC, baseDir = dir_raw)
# raw_counts <- read.delim(file.path(dir_raw, ACC, "<supplementary_file_name>"),
#                           row.names = 1, check.names = FALSE)

## STEP 2. Metadata -----------------------------------------------------
# Group levels expected: "TBneg", "LTBI", "ATB"
# metadata <- data.frame(
#   sample_id = rownames(pheno),
#   group     = factor(pheno$`disease state:ch1`, levels = c("TBneg", "LTBI", "ATB")),
#   replicate = seq_len(nrow(pheno))
# )
# rownames(metadata) <- metadata$sample_id
# write.csv(metadata, file.path(dir_meta, "GSE161829_sample_metadata.csv"), row.names = FALSE)
# plot_sample_design(metadata, out_png = file.path(dir_fig, "GSE161829_sample_design.png"),
#                     title = "GSE161829 sample design")

## STEP 3. Filtering ------------------------------------------------------
# raw_counts <- raw_counts[, rownames(metadata)]
# write.csv(raw_counts, file.path(dir_raw, "GSE161829_raw_count_matrix.csv"))
# filt <- filter_low_expression(raw_counts, min_count = 10, min_samples = min(table(metadata$group)))
# write.csv(filt$filtered, file.path(dir_filt, "GSE161829_filtered_count_matrix.csv"))
# plot_filtering_summary(filt$n_before, filt$n_after,
#                         out_png = file.path(dir_fig, "GSE161829_filtering_summary.png"),
#                         title = "GSE161829 gene filtering")

## STEP 4. DESeq2 + VST + PCA ----------------------------------------------
# dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
# dds <- DESeq(dds)
# vst_mat <- run_vst(dds)
# write.csv(vst_mat, file.path(dir_norm, "GSE161829_vst_normalized_matrix.csv"))
# plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, "GSE161829_PCA.png"),
#          title = "GSE161829 PCA (TBneg / LTBI / ATB)")

## STEP 5. Correlation ------------------------------------------------------
# sample_correlation(vst_mat, metadata,
#                     out_matrix_csv = file.path(dir_corr, "GSE161829_sample_correlation_matrix.csv"),
#                     out_png = file.path(dir_fig, "GSE161829_sample_correlation_heatmap.png"),
#                     title = "GSE161829 sample correlation")

## STEP 6. DEG calling -------------------------------------------------------
# res_ltbi <- run_deg(dds, contrast = c("group", "LTBI", "TBneg"),
#                      out_csv = file.path(dir_deg, "GSE161829_TBneg_vs_LTBI_DEG.csv"),
#                      lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_ltbi, out_png = file.path(dir_fig, "GSE161829_TBneg_vs_LTBI_volcano.png"),
#              title = "GSE161829: TBneg vs LTBI")
#
# res_atb <- run_deg(dds, contrast = c("group", "ATB", "TBneg"),
#                     out_csv = file.path(dir_deg, "GSE161829_TBneg_vs_ATB_DEG.csv"),
#                     lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_atb, out_png = file.path(dir_fig, "GSE161829_TBneg_vs_ATB_volcano.png"),
#              title = "GSE161829: TBneg vs ATB")

## STEP 7. Significant DEGs + heatmaps ---------------------------------------
# sig_ltbi <- get_significant_degs(res_ltbi, file.path(dir_sig, "GSE161829_TBneg_vs_LTBI_sigDEGs.csv"))
# sig_deg_heatmap(vst_mat, sig_ltbi$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE161829_TBneg_vs_LTBI_DEG_heatmap.png"),
#                  title = "GSE161829: Top DEGs, TBneg vs LTBI")
#
# sig_atb <- get_significant_degs(res_atb, file.path(dir_sig, "GSE161829_TBneg_vs_ATB_sigDEGs.csv"))
# sig_deg_heatmap(vst_mat, sig_atb$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE161829_TBneg_vs_ATB_DEG_heatmap.png"),
#                  title = "GSE161829: Top DEGs, TBneg vs ATB")

## STEP 8. Shared/unique within dataset ---------------------------------------
# shared_unique_degs(list(TBneg_vs_LTBI = sig_ltbi$gene, TBneg_vs_ATB = sig_atb$gene),
#                     out_csv = file.path(dir_share, "GSE161829_shared_unique_DEGs.csv"),
#                     out_png = file.path(dir_fig, "GSE161829_upset_within_dataset.png"),
#                     title = "GSE161829: shared/unique DEGs")

## STEP 9. WGCNA exports ---------------------------------------------------
# export_wgcna_expression(vst_mat, file.path(dir_wgcna, "GSE161829_WGCNA_expression_matrix.csv"))
# export_wgcna_traits(metadata,
#   trait_cols = list(LTBI = as.integer(metadata$group == "LTBI"),
#                      ATB  = as.integer(metadata$group == "ATB"),
#                      TBneg = as.integer(metadata$group == "TBneg")),
#   out_csv = file.path(dir_wgcna, "GSE161829_WGCNA_trait_file.csv"))
