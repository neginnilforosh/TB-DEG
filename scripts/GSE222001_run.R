## =============================================================
## GSE222001 -- Healthy vs Active TB, Healthy vs Latent TB
## Comparisons: Healthy_vs_Active, Healthy_vs_Latent
##
## RUN THIS ON A MACHINE WITH INTERNET ACCESS TO NCBI GEO.
## =============================================================
source("00_functions.R")
library(GEOquery)

ACC       <- "GSE222001"
BASE_DIR  <- file.path("..", ACC)     # e.g. TB_DEG_pipeline/GSE222001
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

## -------------------------------------------------------------
## STEP 1. Download raw counts + metadata from GEO
## -------------------------------------------------------------
# (a) Series matrix -> gives you sample phenotype/characteristics data
gse <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)

# (b) Supplementary raw-count file(s) -- inspect the GEO record first to
#     find the exact filename; GEO RNA-seq series usually attach a single
#     gene x sample raw-count matrix (txt/csv/xlsx) as a supplementary file.
getGEOSuppFiles(ACC, baseDir = dir_raw)
# --> after this, open the downloaded file and load it, e.g.:
# raw_counts <- read.delim(file.path(dir_raw, ACC, "GSE222001_raw_counts.txt"),
#                           row.names = 1, check.names = FALSE)

## -------------------------------------------------------------
## STEP 2. Build sample metadata
##    Inspect `pheno` (from step 1a) to identify which column holds the
##    group label (often "characteristics_ch1" / "disease state:ch1").
##    Recode into a factor with levels: "Healthy", "Active", "Latent".
## -------------------------------------------------------------
# metadata <- data.frame(
#   sample_id   = rownames(pheno),
#   group       = factor(pheno$`disease state:ch1`,
#                         levels = c("Healthy", "Active", "Latent")),
#   replicate   = seq_len(nrow(pheno)),
#   batch       = pheno$`batch:ch1`            # only if available, else drop
# )
# rownames(metadata) <- metadata$sample_id
# write.csv(metadata, file.path(dir_meta, "GSE222001_sample_metadata.csv"), row.names = FALSE)
# plot_sample_design(metadata, out_png = file.path(dir_fig, "GSE222001_sample_design.png"),
#                     title = "GSE222001 sample design")

## -------------------------------------------------------------
## STEP 3. Filter low-expression genes
## -------------------------------------------------------------
# raw_counts <- raw_counts[, rownames(metadata)]     # match column order
# write.csv(raw_counts, file.path(dir_raw, "GSE222001_raw_count_matrix.csv"))
#
# filt <- filter_low_expression(raw_counts, min_count = 10,
#                                min_samples = min(table(metadata$group)))
# write.csv(filt$filtered, file.path(dir_filt, "GSE222001_filtered_count_matrix.csv"))
# plot_filtering_summary(filt$n_before, filt$n_after,
#                         out_png = file.path(dir_fig, "GSE222001_filtering_summary.png"),
#                         title = "GSE222001 gene filtering")

## -------------------------------------------------------------
## STEP 4. DESeq2 object, VST normalization, PCA
## -------------------------------------------------------------
# dds <- DESeqDataSetFromMatrix(countData = filt$filtered,
#                                colData   = metadata,
#                                design    = ~ group)
# dds <- DESeq(dds)
#
# vst_mat <- run_vst(dds)
# write.csv(vst_mat, file.path(dir_norm, "GSE222001_vst_normalized_matrix.csv"))
# plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, "GSE222001_PCA.png"),
#          title = "GSE222001 PCA (Healthy / Active / Latent)")

## -------------------------------------------------------------
## STEP 5. Sample correlation matrix + heatmap
## -------------------------------------------------------------
# sample_correlation(vst_mat, metadata,
#                     out_matrix_csv = file.path(dir_corr, "GSE222001_sample_correlation_matrix.csv"),
#                     out_png = file.path(dir_fig, "GSE222001_sample_correlation_heatmap.png"),
#                     title = "GSE222001 sample correlation")

## -------------------------------------------------------------
## STEP 6. DEG calling -- two comparisons
## -------------------------------------------------------------
# res_active <- run_deg(dds, contrast = c("group", "Active", "Healthy"),
#                        out_csv = file.path(dir_deg, "GSE222001_Healthy_vs_Active_DEG.csv"),
#                        lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_active, out_png = file.path(dir_fig, "GSE222001_Healthy_vs_Active_volcano.png"),
#              title = "GSE222001: Healthy vs Active TB", lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
#
# res_latent <- run_deg(dds, contrast = c("group", "Latent", "Healthy"),
#                        out_csv = file.path(dir_deg, "GSE222001_Healthy_vs_Latent_DEG.csv"),
#                        lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# volcano_plot(res_latent, out_png = file.path(dir_fig, "GSE222001_Healthy_vs_Latent_volcano.png"),
#              title = "GSE222001: Healthy vs Latent TB", lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

## -------------------------------------------------------------
## STEP 7. Significant DEG lists + heatmaps
## -------------------------------------------------------------
# sig_active <- get_significant_degs(res_active,
#                 out_csv = file.path(dir_sig, "GSE222001_Healthy_vs_Active_sigDEGs.csv"),
#                 lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# sig_deg_heatmap(vst_mat, sig_active$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE222001_Healthy_vs_Active_DEG_heatmap.png"),
#                  title = "GSE222001: Top DEGs, Healthy vs Active")
#
# sig_latent <- get_significant_degs(res_latent,
#                 out_csv = file.path(dir_sig, "GSE222001_Healthy_vs_Latent_sigDEGs.csv"),
#                 lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
# sig_deg_heatmap(vst_mat, sig_latent$gene, metadata,
#                  out_png = file.path(dir_fig, "GSE222001_Healthy_vs_Latent_DEG_heatmap.png"),
#                  title = "GSE222001: Top DEGs, Healthy vs Latent")

## -------------------------------------------------------------
## STEP 8. Shared/unique DEGs between the two comparisons (UpSet)
##    NOTE: cross-dataset / cross-comorbidity UpSet plots (TB-DM, TB-HIV,
##    TB-LC, etc.) are built later, once all four datasets have been run --
##    see scripts/05_cross_dataset_upset.R
## -------------------------------------------------------------
# shared_unique_degs(list(Healthy_vs_Active = sig_active$gene,
#                          Healthy_vs_Latent = sig_latent$gene),
#                     out_csv = file.path(dir_share, "GSE222001_shared_unique_DEGs.csv"),
#                     out_png = file.path(dir_fig, "GSE222001_upset_within_dataset.png"),
#                     title = "GSE222001: shared/unique DEGs")

## -------------------------------------------------------------
## STEP 9. WGCNA-ready exports (full VST matrix, NOT DEG-restricted)
## -------------------------------------------------------------
# export_wgcna_expression(vst_mat, file.path(dir_wgcna, "GSE222001_WGCNA_expression_matrix.csv"))
# export_wgcna_traits(metadata,
#   trait_cols = list(Active = as.integer(metadata$group == "Active"),
#                      Latent = as.integer(metadata$group == "Latent"),
#                      Healthy = as.integer(metadata$group == "Healthy")),
#   out_csv = file.path(dir_wgcna, "GSE222001_WGCNA_trait_file.csv"))
