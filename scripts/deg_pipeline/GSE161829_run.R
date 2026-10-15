## GSE161829: TBneg vs LTBI and TBneg vs ATB (DESeq2).

source("00_functions.R")
library(GEOquery)
library(DESeq2)

ACC      <- "GSE161829"
BASE_DIR <- file.path("..", "..", ACC)
LFC_TH   <- 1
PADJ_TH  <- 0.05

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

for (d in c(dir_raw, dir_meta, dir_filt, dir_norm, dir_corr,
            dir_deg, dir_sig, dir_share, dir_wgcna, dir_fig)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

gse   <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)

## ---- CONFIG ----

GROUP_COL   <- "disease group:ch1"
COUNT_FILE  <- "GSE161829_RawCounts_CD4m.txt.gz"

GROUP_MAP <- c(
  "TBneg"              = "TBneg",
  "LTBI"               = "LTBI",
  "ATB_Diagnosis"      = "ATB",
  "ATB_Post_treatment" = "Exclude"
)

## ---- PART B ----

## ---- Load counts ----
count_path <- file.path(dir_raw, ACC, COUNT_FILE)
raw_counts <- read.delim(count_path, row.names = 1, check.names = FALSE)

## ---- Match count columns to GEO samples by NAME (exact one-to-one) ----
## GEO's !Sample_description records each GSM's count-column name (Sample1_CD4m ... Sample106_CD4m). Its
## order differs from pData() order at GEO positions 74-101, so renaming the columns by position mislabelled
## 7 analysed samples in an earlier run. Every check below stops the script instead of guessing.
desc_cols <- grep("^description", names(pheno), value = TRUE)
if (!length(desc_cols)) stop("GEO metadata has no !Sample_description column; cannot match count columns by name.")
gsm_col <- apply(pheno[, desc_cols, drop = FALSE], 1, function(v) {
  hit <- unique(grep("^Sample[0-9]+_CD4m$", trimws(as.character(v)), value = TRUE))
  if (length(hit) == 1) hit else NA_character_
})
names(gsm_col) <- rownames(pheno)
if (anyNA(gsm_col))
  stop("GSM(s) without exactly one count-column name in !Sample_description: ",
       paste(names(gsm_col)[is.na(gsm_col)], collapse = ", "))
if (anyDuplicated(gsm_col))
  stop("Count-column name(s) claimed by more than one GSM: ",
       paste(unique(gsm_col[duplicated(gsm_col)]), collapse = ", "))
if (anyDuplicated(colnames(raw_counts)))
  stop("Duplicated count-matrix column name(s): ",
       paste(unique(colnames(raw_counts)[duplicated(colnames(raw_counts))]), collapse = ", "))
cols_without_gsm <- setdiff(colnames(raw_counts), gsm_col)
gsm_without_col  <- names(gsm_col)[!gsm_col %in% colnames(raw_counts)]
if (length(cols_without_gsm) || length(gsm_without_col))
  stop("Count columns and GEO samples differ. Count columns with no GSM: ",
       paste(cols_without_gsm, collapse = ", "), " | GSMs with no count column: ",
       paste(gsm_without_col, collapse = ", "))
raw_counts <- raw_counts[, unname(gsm_col), drop = FALSE]      # reorder to pData() order ...
colnames(raw_counts) <- names(gsm_col)                         # ... and name each column by its GSM
if (!identical(colnames(raw_counts), rownames(pheno)))
  stop("Count columns are not in GEO sample order after matching.")

## Independent check of the group labels: each GEO title encodes the same group as the disease-group field
## (..._ATB_V1_... = ATB_Diagnosis, ..._ATB_V2_... = ATB_Post_treatment, ..._LTBI_... = LTBI, ..._TBneg_... = TBneg).
geo_group   <- as.character(pheno[[GROUP_COL]])
title_group <- ifelse(grepl("_ATB_V1", pheno$title), "ATB_Diagnosis",
               ifelse(grepl("_ATB_V2", pheno$title), "ATB_Post_treatment",
               ifelse(grepl("_LTBI_",  pheno$title), "LTBI",
               ifelse(grepl("_TBneg_", pheno$title), "TBneg", NA_character_))))
bad <- which(is.na(title_group) | title_group != geo_group)
if (length(bad))
  stop("GEO title and '", GROUP_COL, "' disagree for: ", paste(rownames(pheno)[bad], collapse = ", "))
cat("Count columns matched one-to-one to", length(gsm_col), "GEO samples via !Sample_description\n")

cat("Raw counts loaded:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

## ---- Build metadata ----
metadata <- data.frame(
  sample_id = rownames(pheno),
  group     = unname(GROUP_MAP[as.character(pheno[[GROUP_COL]])]),
  stringsAsFactors = FALSE
)
if (anyNA(metadata$group))
  stop("Unmapped GEO group label(s): ", paste(unique(geo_group[is.na(metadata$group)]), collapse = ", "))

n_before <- ncol(raw_counts)
cat("\n--- Samples before exclusion:", n_before, "---\n"); print(table(GEO_group = geo_group))
excluded <- metadata$sample_id[metadata$group == "Exclude"]     # ATB_Post_treatment, excluded by design
metadata <- metadata[metadata$group != "Exclude", ]
metadata$group <- factor(metadata$group, levels = c("TBneg", "LTBI", "ATB"))
rownames(metadata) <- metadata$sample_id

common <- intersect(colnames(raw_counts), metadata$sample_id)
silently_dropped <- setdiff(colnames(raw_counts), c(common, excluded))
missing_counts   <- setdiff(metadata$sample_id, common)
if (length(silently_dropped) || length(missing_counts))
  stop("Unexpected sample loss. Count columns dropped without an exclusion rule: ",
       paste(silently_dropped, collapse = ", "), " | metadata samples without counts: ",
       paste(missing_counts, collapse = ", "))
raw_counts <- raw_counts[, common, drop = FALSE]
metadata   <- metadata[common, , drop = FALSE]
if (ncol(raw_counts) + length(excluded) != n_before)
  stop("Sample accounting failed: ", ncol(raw_counts), " kept + ", length(excluded), " excluded != ", n_before)
if (!identical(colnames(raw_counts), rownames(metadata)))
  stop("Count columns and metadata rows are not in the same order.")
cat("Samples after excluding ATB_Post_treatment:", ncol(raw_counts), "(excluded:", length(excluded), ")\n")

cat("\n--- Final Sample Count per Group ---\n")
print(table(metadata$group))

write.csv(metadata, file.path(dir_meta, paste0(ACC, "_sample_metadata.csv")), row.names = FALSE)
plot_sample_design(metadata, out_png = file.path(dir_fig, paste0(ACC, "_sample_design.png")),
                   title = paste(ACC, "sample design"))
## ---- Filter ----
write.csv(raw_counts, file.path(dir_raw, paste0(ACC, "_raw_count_matrix.csv")))
filt <- filter_low_expression(raw_counts, min_count = 10,
                              min_samples = min(table(metadata$group)))
write.csv(filt$filtered, file.path(dir_filt, paste0(ACC, "_filtered_count_matrix.csv")))
plot_filtering_summary(filt$n_before, filt$n_after,
                       out_png = file.path(dir_fig, paste0(ACC, "_filtering_summary.png")),
                       title = paste(ACC, "gene filtering"))

## ---- DESeq2 + VST + PCA ----
if (!identical(colnames(filt$filtered), rownames(metadata)))
  stop("Filtered count matrix and metadata are not in the same sample order.")
dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
dds <- DESeq(dds)

vst_mat <- run_vst(dds)
write.csv(vst_mat, file.path(dir_norm, paste0(ACC, "_vst_normalized_matrix.csv")))
plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, paste0(ACC, "_PCA.png")),
         title = paste(ACC, "PCA (TBneg / LTBI / ATB)"))

## ---- Sample correlation ----
sample_correlation(vst_mat, metadata,
                   out_matrix_csv = file.path(dir_corr, paste0(ACC, "_sample_correlation_matrix.csv")),
                   out_png = file.path(dir_fig, paste0(ACC, "_sample_correlation_heatmap.png")),
                   title = paste(ACC, "sample correlation"))

## ---- Outlier check ----
corr_mat <- cor(vst_mat, method = "pearson")
avg_corr <- sort(rowMeans(corr_mat))
cat("\n--- Lowest mean sample correlation (possible outliers) ---\n")
print(head(avg_corr, 6))

## ---- DEGs: two comparisons ----
res_ltbi <- run_deg(dds, contrast = c("group", "LTBI", "TBneg"),
                    out_csv = file.path(dir_deg, paste0(ACC, "_TBneg_vs_LTBI_DEG.csv")),
                    lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
volcano_plot(res_ltbi, out_png = file.path(dir_fig, paste0(ACC, "_TBneg_vs_LTBI_volcano.png")),
             title = paste(ACC, ": TBneg vs LTBI"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

res_atb <- run_deg(dds, contrast = c("group", "ATB", "TBneg"),
                   out_csv = file.path(dir_deg, paste0(ACC, "_TBneg_vs_ATB_DEG.csv")),
                   lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
volcano_plot(res_atb, out_png = file.path(dir_fig, paste0(ACC, "_TBneg_vs_ATB_volcano.png")),
             title = paste(ACC, ": TBneg vs ATB"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

## ---- Significant DEGs + heatmaps ----
sig_ltbi <- get_significant_degs(res_ltbi,
              out_csv = file.path(dir_sig, paste0(ACC, "_TBneg_vs_LTBI_sigDEGs.csv")),
              lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
if (nrow(sig_ltbi) > 1) {
  sig_deg_heatmap(vst_mat, sig_ltbi$gene, metadata,
                  out_png = file.path(dir_fig, paste0(ACC, "_TBneg_vs_LTBI_DEG_heatmap.png")),
                  title = paste(ACC, ": Top DEGs, TBneg vs LTBI"))
}

sig_atb <- get_significant_degs(res_atb,
             out_csv = file.path(dir_sig, paste0(ACC, "_TBneg_vs_ATB_sigDEGs.csv")),
             lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
if (nrow(sig_atb) > 1) {
  sig_deg_heatmap(vst_mat, sig_atb$gene, metadata,
                  out_png = file.path(dir_fig, paste0(ACC, "_TBneg_vs_ATB_DEG_heatmap.png")),
                  title = paste(ACC, ": Top DEGs, TBneg vs ATB"))
}

## ---- Shared/unique between the two comparisons ----
if (nrow(sig_ltbi) > 0 && nrow(sig_atb) > 0) {
  shared_unique_degs(list(TBneg_vs_LTBI = sig_ltbi$gene, TBneg_vs_ATB = sig_atb$gene),
                     out_csv = file.path(dir_share, paste0(ACC, "_shared_unique_DEGs.csv")),
                     out_png = file.path(dir_fig, paste0(ACC, "_upset_within_dataset.png")),
                     title = paste(ACC, ": shared/unique DEGs"))
}

## ---- WGCNA exports ----
export_wgcna_expression(vst_mat, file.path(dir_wgcna, paste0(ACC, "_WGCNA_expression_matrix.csv")))
export_wgcna_traits(metadata,
  trait_cols = list(TBneg = as.integer(metadata$group == "TBneg"),
                    LTBI  = as.integer(metadata$group == "LTBI"),
                    ATB   = as.integer(metadata$group == "ATB")),
  out_csv = file.path(dir_wgcna, paste0(ACC, "_WGCNA_trait_file.csv")))

cat("\n========== DEG COUNT SUMMARY ==========\n")
cat("TBneg vs LTBI : ", nrow(sig_ltbi), " significant DEGs\n", sep = "")
cat("TBneg vs ATB  : ", nrow(sig_atb), " significant DEGs\n", sep = "")
cat("(threshold: padj < ", PADJ_TH, ", |log2FC| >= ", LFC_TH, ")\n", sep = "")
cat("=======================================\n")
