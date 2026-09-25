## =============================================================
## GSE99374 -- HC vs LTBI (single comparison)
## =============================================================
source("00_functions.R")
library(GEOquery)

ACC       <- "GSE99374"
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
metadata <- read.csv("../GSE99374/01_raw_counts/GSE99374_CD8_sample_metadata.csv")
metadata$group <- factor(metadata$group, levels = c("HC", "LTBI"))
rownames(metadata) <- metadata$sample_id
write.csv(metadata, file.path(dir_meta, "GSE99374_sample_metadata.csv"), row.names = FALSE)
plot_sample_design(metadata, out_png = file.path(dir_fig, "GSE99374_sample_design.png"))

raw_counts <- read.csv("../GSE99374/01_raw_counts/GSE99374_CD8_raw_counts.csv", row.names = 1, check.names = FALSE)
raw_counts <- raw_counts[, rownames(metadata)]

## STEP 3. Filtering (RNA-seq path) ------------------------------------------
raw_counts <- raw_counts[, rownames(metadata)]
write.csv(raw_counts, file.path(dir_raw, "GSE99374_raw_count_matrix.csv"))
filt <- filter_low_expression(raw_counts, min_count = 10, min_samples = min(table(metadata$group)))
write.csv(filt$filtered, file.path(dir_filt, "GSE99374_filtered_count_matrix.csv"))
plot_filtering_summary(filt$n_before, filt$n_after,
                        out_png = file.path(dir_fig, "GSE99374_filtering_summary.png"),
                           title = "GSE99374 gene filtering")

## STEP 4. DESeq2 + VST + PCA (RNA-seq path) ---------------------------------
dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
dds <- DESeq(dds)
vst_mat <- run_vst(dds)
write.csv(vst_mat, file.path(dir_norm, "GSE99374_vst_normalized_matrix.csv"))
plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, "GSE99374_PCA.png"),
          title = "GSE99374 PCA (HC / LTBI)")

## STEP 5. Correlation --------------------------------------------------------
sample_correlation(vst_mat, metadata,
                    out_matrix_csv = file.path(dir_corr, "GSE99374_sample_correlation_matrix.csv"),
                    out_png = file.path(dir_fig, "GSE99374_sample_correlation_heatmap.png"),
                    title = "GSE99374 sample correlation")

## STEP 6. DEG calling ---------------------------------------------------------
res_ltbi <- run_deg(dds, contrast = c("group", "LTBI", "HC"),
                     out_csv = file.path(dir_deg, "GSE99374_HC_vs_LTBI_DEG.csv"),
                     lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
volcano_plot(res_ltbi, out_png = file.path(dir_fig, "GSE99374_HC_vs_LTBI_volcano.png"),
             title = "GSE99374: HC vs LTBI")

## STEP 7. Significant DEGs + heatmap -------------------------------------------
sig_ltbi <- get_significant_degs(res_ltbi, file.path(dir_sig, "GSE99374_HC_vs_LTBI_sigDEGs.csv"))
sig_deg_heatmap(vst_mat, sig_ltbi$gene, metadata,
                 out_png = file.path(dir_fig, "GSE99374_HC_vs_LTBI_DEG_heatmap.png"),
                 title = "GSE99374: Top DEGs, HC vs LTBI")


## STEP 8. WGCNA exports ---------------------------------------------------
export_wgcna_expression(vst_mat, file.path(dir_wgcna, "GSE99374_WGCNA_expression_matrix.csv"))
export_wgcna_traits(metadata,
  trait_cols = list(LTBI = as.integer(metadata$group == "LTBI"),
                     HC   = as.integer(metadata$group == "HC")),
  out_csv = file.path(dir_wgcna, "GSE99374_WGCNA_trait_file.csv"))
