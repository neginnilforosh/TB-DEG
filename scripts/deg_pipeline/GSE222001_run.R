## =============================================================
## GSE222001 — Healthy vs Active TB, Healthy vs Latent TB
## =============================================================

source("00_functions.R")
library(GEOquery)
library(DESeq2)

ACC      <- "GSE222001"
BASE_DIR <- file.path("..", ACC)
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

## Create all folders up front so nothing fails later on a missing dir
for (d in c(dir_raw, dir_meta, dir_filt, dir_norm, dir_corr,
            dir_deg, dir_sig, dir_share, dir_wgcna, dir_fig)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


## =============================================================
## PART A 
## =============================================================

gse   <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)
getGEOSuppFiles(ACC, baseDir = dir_raw)

## --- A1. Which column holds the group label? ---
## Candidate columns (these are the ones that actually vary across samples):
vary <- names(pheno)[sapply(pheno, function(x) length(unique(x)) > 1 &&
                                                 length(unique(x)) <= 10)]
cat("\n--- Candidate group columns ---\n")
for (v in vary) cat(sprintf("  %-40s : %s\n", v, paste(unique(pheno[[v]]), collapse = " | ")))

## --- A2. What raw count file did GEO give us? ---
cat("\n--- Downloaded supplementary files ---\n")
print(list.files(file.path(dir_raw, ACC), full.names = FALSE))


## =============================================================
## CONFIG 
## =============================================================

GROUP_COL   <- "disease state:ch1"        # <-- from A1
COUNT_FILE  <- "GSE222001_count.txt" # <-- from A2

## Map the raw labels GEO uses onto clean names. 
GROUP_MAP <- c(
  "Healthy control"                = "Healthy",
  "Active tuberculosis infection"  = "Active",
  "Latent tuberculosis infection"  = "Latent"
)


## =============================================================
## PART B — runs straight through, no edits needed
## =============================================================

## ---- Load counts ----
count_path <- file.path(dir_raw, ACC, COUNT_FILE)
raw_counts <- if (grepl("\\.csv(\\.gz)?$", COUNT_FILE)) {
  read.csv(count_path, row.names = 1, check.names = FALSE)
} else if (grepl("\\.xlsx?$", COUNT_FILE)) {
  as.data.frame(readxl::read_excel(count_path)) |>
    (\(d) { rownames(d) <- d[[1]]; d[-1] })()
} else {
  read.delim(count_path, row.names = 1, check.names = FALSE)
}

# Remove featureCounts genomic annotation columns
annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
raw_counts <- raw_counts[, !colnames(raw_counts) %in% annot_cols, drop = FALSE]

cat("Raw counts loaded:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

## ---- Build metadata ----
metadata <- data.frame(
  sample_id = colnames(raw_counts),
  group     = unname(GROUP_MAP[as.character(pheno[[GROUP_COL]])]),
  stringsAsFactors = FALSE
)
stopifnot(!any(is.na(metadata$group)))
metadata$group <- factor(metadata$group, levels = c("Healthy", "Active", "Latent"))
rownames(metadata) <- metadata$sample_id

## Keep only samples present in both the counts and the metadata
common <- intersect(colnames(raw_counts), metadata$sample_id)
stopifnot(length(common) > 0)
raw_counts <- raw_counts[, common]
metadata   <- metadata[common, ]
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
dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
dds <- DESeq(dds)

vst_mat <- run_vst(dds)
write.csv(vst_mat, file.path(dir_norm, paste0(ACC, "_vst_normalized_matrix.csv")))
plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, paste0(ACC, "_PCA.png")),
         title = paste(ACC, "PCA (Healthy / Active / Latent)"))

## ---- Sample correlation ----
sample_correlation(vst_mat, metadata,
                   out_matrix_csv = file.path(dir_corr, paste0(ACC, "_sample_correlation_matrix.csv")),
                   out_png = file.path(dir_fig, paste0(ACC, "_sample_correlation_heatmap.png")),
                   title = paste(ACC, "sample correlation"))

## ---- Outlier check ----
## Ranks samples by mean correlation to all others. Anything sitting well
## below the pack is worth a sensitivity re-run before trusting the DEGs.
corr_mat <- cor(vst_mat, method = "pearson")
avg_corr <- sort(rowMeans(corr_mat))
cat("\n--- Lowest mean sample correlation (possible outliers) ---\n")
print(head(avg_corr, 6))

## ---- DEGs: two comparisons ----
res_active <- run_deg(dds, contrast = c("group", "Active", "Healthy"),
                      out_csv = file.path(dir_deg, paste0(ACC, "_Healthy_vs_Active_DEG.csv")),
                      lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
volcano_plot(res_active, out_png = file.path(dir_fig, paste0(ACC, "_Healthy_vs_Active_volcano.png")),
             title = paste(ACC, ": Healthy vs Active TB"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

res_latent <- run_deg(dds, contrast = c("group", "Latent", "Healthy"),
                      out_csv = file.path(dir_deg, paste0(ACC, "_Healthy_vs_Latent_DEG.csv")),
                      lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
volcano_plot(res_latent, out_png = file.path(dir_fig, paste0(ACC, "_Healthy_vs_Latent_volcano.png")),
             title = paste(ACC, ": Healthy vs Latent TB"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

## ---- Significant DEGs + heatmaps ----
sig_active <- get_significant_degs(res_active,
                out_csv = file.path(dir_sig, paste0(ACC, "_Healthy_vs_Active_sigDEGs.csv")),
                lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
if (nrow(sig_active) > 1) {
  sig_deg_heatmap(vst_mat, sig_active$gene, metadata,
                  out_png = file.path(dir_fig, paste0(ACC, "_Healthy_vs_Active_DEG_heatmap.png")),
                  title = paste(ACC, ": Top DEGs, Healthy vs Active"))
}

sig_latent <- get_significant_degs(res_latent,
                out_csv = file.path(dir_sig, paste0(ACC, "_Healthy_vs_Latent_sigDEGs.csv")),
                lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
if (nrow(sig_latent) > 1) {
  sig_deg_heatmap(vst_mat, sig_latent$gene, metadata,
                  out_png = file.path(dir_fig, paste0(ACC, "_Healthy_vs_Latent_DEG_heatmap.png")),
                  title = paste(ACC, ": Top DEGs, Healthy vs Latent"))
}

## ---- Shared/unique between the two comparisons ----
if (nrow(sig_active) > 0 && nrow(sig_latent) > 0) {
  shared_unique_degs(list(Healthy_vs_Active = sig_active$gene,
                          Healthy_vs_Latent = sig_latent$gene),
                     out_csv = file.path(dir_share, paste0(ACC, "_shared_unique_DEGs.csv")),
                     out_png = file.path(dir_fig, paste0(ACC, "_upset_within_dataset.png")),
                     title = paste(ACC, ": shared/unique DEGs"))
}

## ---- WGCNA exports (full VST matrix, NOT DEG-restricted) ----
export_wgcna_expression(vst_mat, file.path(dir_wgcna, paste0(ACC, "_WGCNA_expression_matrix.csv")))
export_wgcna_traits(metadata,
  trait_cols = list(Healthy = as.integer(metadata$group == "Healthy"),
                    Active  = as.integer(metadata$group == "Active"),
                    Latent  = as.integer(metadata$group == "Latent")),
  out_csv = file.path(dir_wgcna, paste0(ACC, "_WGCNA_trait_file.csv")))

## ---- The number of significants ----
cat("\n========== DEG COUNT SUMMARY ==========\n")
cat("Healthy vs Active : ", nrow(sig_active), " significant DEGs\n", sep = "")
cat("Healthy vs Latent : ", nrow(sig_latent), " significant DEGs\n", sep = "")
cat("(threshold: padj < ", PADJ_TH, ", |log2FC| >= ", LFC_TH, ")\n", sep = "")
cat("=======================================\n")

