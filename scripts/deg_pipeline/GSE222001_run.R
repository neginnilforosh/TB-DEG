## GSE222001: Healthy vs Active TB and Healthy vs Latent TB (DESeq2).

source("00_functions.R")
library(GEOquery)
library(DESeq2)

ACC      <- "GSE222001"
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

## Create output folders
for (d in c(dir_raw, dir_meta, dir_filt, dir_norm, dir_corr,
            dir_deg, dir_sig, dir_share, dir_wgcna, dir_fig)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- PART A ----

gse   <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)
getGEOSuppFiles(ACC, baseDir = dir_raw)

## --- A1. Group-label column (candidates that vary across samples): ---
vary <- names(pheno)[sapply(pheno, function(x) length(unique(x)) > 1 &&
                                                 length(unique(x)) <= 10)]
cat("\n--- Candidate group columns ---\n")
for (v in vary) cat(sprintf("  %-40s : %s\n", v, paste(unique(pheno[[v]]), collapse = " | ")))

## --- A2. What raw count file did GEO give us? ---
cat("\n--- Downloaded supplementary files ---\n")
print(list.files(file.path(dir_raw, ACC), full.names = FALSE))

## ---- CONFIG ----

GROUP_COL   <- "disease state:ch1"        # <-- from A1
COUNT_FILE  <- "GSE222001_count.txt.gz" # <-- from A2 (GEO supplies only the gzipped file; read.delim reads it directly)

## Map the raw labels GEO uses onto clean names.
GROUP_MAP <- c(
  "Healthy control"                = "Healthy",
  "Active tuberculosis infection"  = "Active",
  "Latent tuberculosis infection"  = "Latent"
)

## ---- PART B ----

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

## ---- Match count columns to GEO samples by NAME (exact one-to-one) ----
## GEO's !Sample_description holds each GSM's count-column name (TB-rep1, LTBI-rep1, NC-rep1, ...).
## The count file lists its columns in a different order from pData(), so pairing them by position
## is wrong (it swapped the Healthy and Latent labels in the existing GSE222001 outputs).
## Every check below stops the script instead of guessing.
desc_cols <- grep("^description", names(pheno), value = TRUE)
if (!length(desc_cols)) stop("GEO metadata has no !Sample_description column; cannot match count columns by name.")
gsm_name <- apply(pheno[, desc_cols, drop = FALSE], 1, function(v) {
  hit <- unique(grep("^(TB|LTBI|NC)-rep[0-9]+$", trimws(as.character(v)), value = TRUE))
  if (length(hit) == 1) hit else NA_character_
})
names(gsm_name) <- rownames(pheno)
if (anyNA(gsm_name))
  stop("GSM(s) without exactly one count-column name in !Sample_description: ",
       paste(names(gsm_name)[is.na(gsm_name)], collapse = ", "))
if (anyDuplicated(gsm_name))
  stop("Count-column name(s) claimed by more than one GSM: ",
       paste(unique(gsm_name[duplicated(gsm_name)]), collapse = ", "))
if (anyDuplicated(colnames(raw_counts)))
  stop("Duplicated count-matrix column name(s): ",
       paste(unique(colnames(raw_counts)[duplicated(colnames(raw_counts))]), collapse = ", "))
cols_without_gsm <- setdiff(colnames(raw_counts), gsm_name)
gsm_without_col  <- names(gsm_name)[!gsm_name %in% colnames(raw_counts)]
if (length(cols_without_gsm) || length(gsm_without_col))
  stop("Count columns and GEO samples differ. Count columns with no GSM: ",
       paste(cols_without_gsm, collapse = ", "), " | GSMs with no count column: ",
       paste(gsm_without_col, collapse = ", "))
col_gsm <- setNames(names(gsm_name)[match(colnames(raw_counts), gsm_name)], colnames(raw_counts))
if (anyNA(col_gsm) || anyDuplicated(col_gsm)) stop("Count-column to GSM mapping is not one-to-one.")

## Independent check of the group labels against GEO: the column-name prefix, the GEO disease state and the
## GEO title must agree (TB-repN = "Active tuberculosis infection replicate N", LTBI-repN = "Latent
## tuberculosis infection replicate N", NC-repN = "Healthy control replicate N").
PREFIX_STATE <- c(TB = "Active tuberculosis infection", LTBI = "Latent tuberculosis infection",
                  NC = "Healthy control")
geo_state <- as.character(pheno[col_gsm, GROUP_COL])
geo_title <- as.character(pheno[col_gsm, "title"])
prefix    <- sub("-rep[0-9]+$", "", colnames(raw_counts))
rep_no    <- sub("^.*-rep", "", colnames(raw_counts))
bad <- which(is.na(PREFIX_STATE[prefix]) | unname(PREFIX_STATE[prefix]) != geo_state |
             geo_title != paste(geo_state, "replicate", rep_no))
if (length(bad))
  stop("Group check failed (count-column name, GEO disease state and GEO title disagree) for: ",
       paste(colnames(raw_counts)[bad], collapse = ", "))
cat("Count columns matched one-to-one to", length(col_gsm), "GEO samples; group labels agree with GEO:\n")
print(data.frame(count_column = names(col_gsm), gsm = unname(col_gsm), geo_disease_state = geo_state,
                 row.names = NULL))

## ---- Build metadata ----
metadata <- data.frame(
  sample_id = colnames(raw_counts),
  group     = unname(GROUP_MAP[as.character(pheno[col_gsm, GROUP_COL])]),
  stringsAsFactors = FALSE
)
if (anyNA(metadata$group))
  stop("Unmapped GEO group label(s): ", paste(unique(geo_state[is.na(metadata$group)]), collapse = ", "))
metadata$group <- factor(metadata$group, levels = c("Healthy", "Active", "Latent"))
rownames(metadata) <- metadata$sample_id

## No group is excluded in this dataset: every sample must be kept
n_before <- ncol(raw_counts)
cat("\n--- Samples before exclusion:", n_before, "---\n")
common <- intersect(colnames(raw_counts), metadata$sample_id)
if (length(common) != n_before || nrow(metadata) != n_before)
  stop("Unexpected sample loss: ", n_before, " count columns, ", nrow(metadata), " metadata rows, ",
       length(common), " in both. Missing: ", paste(setdiff(colnames(raw_counts), common), collapse = ", "))
raw_counts <- raw_counts[, common, drop = FALSE]
metadata   <- metadata[common, , drop = FALSE]
if (!identical(colnames(raw_counts), rownames(metadata)))
  stop("Count columns and metadata rows are not in the same order.")
cat("Samples after exclusion:", ncol(raw_counts), "(none excluded)\n")
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
         title = paste(ACC, "PCA (Healthy / Active / Latent)"))

## ---- Sample correlation ----
sample_correlation(vst_mat, metadata,
                   out_matrix_csv = file.path(dir_corr, paste0(ACC, "_sample_correlation_matrix.csv")),
                   out_png = file.path(dir_fig, paste0(ACC, "_sample_correlation_heatmap.png")),
                   title = paste(ACC, "sample correlation"))

## ---- Outlier check: samples with low mean correlation to the others ----
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

## ---- Summary: number of significant DEGs ----
cat("\n========== DEG COUNT SUMMARY ==========\n")
cat("Healthy vs Active : ", nrow(sig_active), " significant DEGs\n", sep = "")
cat("Healthy vs Latent : ", nrow(sig_latent), " significant DEGs\n", sep = "")
cat("(threshold: padj < ", PADJ_TH, ", |log2FC| >= ", LFC_TH, ")\n", sep = "")
cat("=======================================\n")

