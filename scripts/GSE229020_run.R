## =============================================================
## GSE229020 -- HC vs Latent TB, HC vs Drug Resistant TB
## Comparisons: HC_vs_LatentTB, HC_vs_DrugResistant
##
## RUN THIS ON A MACHINE WITH INTERNET ACCESS TO NCBI GEO.
## Structure mirrors GSE222001_run.R -- see that file for inline
## explanations of every step.
## =============================================================
source("00_functions.R")
library(GEOquery)

ACC       <- "GSE229020"
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

## STEP 2. Metadata -----------------------------------------------------
# Group levels expected: "HC", "LatentTB", "DrugResistant"
# metadata <- data.frame(
#   sample_id = rownames(pheno),
#   group     = factor(pheno$`disease state:ch1`,
#                       levels = c("HC", "LatentTB", "DrugResistant")),
#   replicate = seq_len(nrow(pheno))
# )
# rownames(metadata) <- metadata$sample_id
# write.csv(metadata, file.path(dir_meta, "GSE229020_sample_metadata.csv"), row.names = FALSE)
# plot_sample_design(metadata, out_png = file.path(dir_fig, "GSE229020_sample_design.png"),
#                     title = "GSE229020 sample design")

## STEP 3. Filtering ------------------------------------------------------
# raw_counts <- raw_counts[, rownames(metadata)]
# write.csv(raw_counts, file.path(dir_raw, "GSE229020_raw_count_matrix.csv"))
# filt <- filter_low_expression(raw_counts, min_count = 10, min_samples = min(table(metadata$group)))
# write.csv(filt$filtered, file.path(dir_filt, "GSE229020_filtered_count_matrix.csv"))
# plot_filtering_summary(filt$n_before, filt$n_after,
#                         out_png = file.path(dir_fig, "GSE229020_filtering_summary.png"),
#                         title = "GSE229020 gene filtering")

## STEP 4. DESeq2 + VST + PCA ----------------------------------------------
# dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
# dds <- DESeq(dds)
# vst_mat <- run_vst(dds)
# write.csv(vst_mat, file.path(dir_norm, "GSE229020_vst_normalized_matrix.csv"))
# plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, "GSE229020_PCA.png"),
#          title = "GSE229020 PCA (HC / Latent TB / Drug Resistant)")

## STEP 5. Correlation ------------------------------------------------------
# sample_correlation(vst_mat, metadata,
#                     out_matrix_csv = file.path(dir_corr, "GSE229020_sample_correlation_matrix.csv"),
#                     out_png = file.path(dir_fig, "GSE229020_sample_correlation_heatmap.png"),
#                     title = "GSE229020 sample correlation")

## STEP 6. DEG calling -------------------------------------------------------
# res_latent <- run_deg(dds, contrast = c("group", "LatentTB", "HC"),
#                        out_csv = file.path(dir_deg, "GSE229020_HC_vs_LatentTB_DEG.csv"),
#                        lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_latent, out_png = file.path(dir_fig, "GSE229020_HC_vs_LatentTB_volcano.png"),
#              title = "GSE229020: HC vs Latent TB")
#
# res_dr <- run_deg(dds, contrast = c("group", "DrugResistant", "HC"),
#                    out_csv = file.path(dir_deg, "GSE229020_HC_vs_DrugResistant_DEG.csv"),
#                    lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_dr, out_png = file.path(dir_fig, "GSE229020_HC_vs_DrugResistant_volcano.png"),
#              title = "GSE229020: HC vs Drug Resistant TB")

## STEP 7. Significant DEGs + heatmaps ---------------------------------------
# sig_latent <- get_significant_degs(res_latent, file.path(dir_sig, "GSE229020_HC_vs_LatentTB_sigDEGs.csv"))
# sig_deg_heatmap(vst_mat, sig_latent$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE229020_HC_vs_LatentTB_DEG_heatmap.png"),
#                  title = "GSE229020: Top DEGs, HC vs Latent TB")
#
# sig_dr <- get_significant_degs(res_dr, file.path(dir_sig, "GSE229020_HC_vs_DrugResistant_sigDEGs.csv"))
# sig_deg_heatmap(vst_mat, sig_dr$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE229020_HC_vs_DrugResistant_DEG_heatmap.png"),
#                  title = "GSE229020: Top DEGs, HC vs Drug Resistant")

## STEP 8. Shared/unique within dataset ---------------------------------------
# shared_unique_degs(list(HC_vs_LatentTB = sig_latent$gene, HC_vs_DrugResistant = sig_dr$gene),
#                     out_csv = file.path(dir_share, "GSE229020_shared_unique_DEGs.csv"),
#                     out_png = file.path(dir_fig, "GSE229020_upset_within_dataset.png"),
#                     title = "GSE229020: shared/unique DEGs")

## STEP 9. WGCNA exports ---------------------------------------------------
# export_wgcna_expression(vst_mat, file.path(dir_wgcna, "GSE229020_WGCNA_expression_matrix.csv"))
# export_wgcna_traits(metadata,
#   trait_cols = list(LatentTB = as.integer(metadata$group == "LatentTB"),
#                      DrugResistant = as.integer(metadata$group == "DrugResistant"),
#                      HC = as.integer(metadata$group == "HC")),
#   out_csv = file.path(dir_wgcna, "GSE229020_WGCNA_trait_file.csv"))
