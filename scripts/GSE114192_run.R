## =============================================================
## GSE114192 -- Healthy_Control vs TB_Only
## =============================================================

source("00_functions.R")
library(GEOquery)
library(DESeq2)
library(stringr)

ACC       <- "GSE114192"
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

for (d in c(dir_raw, dir_meta, dir_filt, dir_norm, dir_corr,
            dir_deg, dir_sig, dir_share, dir_wgcna, dir_fig)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- Download and Load Metadata ----
cat("\nDownloading metadata from GEO...\n")
gse   <- getGEO(ACC, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pheno <- pData(gse)

## =============================================================
## CONFIGURATION
## =============================================================

GROUP_COL <- "disease state (disease_category):ch1"

GROUP_MAP <- c(
  "Healthy_Control" = "Healthy_Control",
  "TB_only"         = "TB_Only",
  "TB_DM"           = "Exclude",
  "DM_only"         = "Exclude",
  "TB_IH"           = "Exclude",
  "IH"              = "Exclude"
)

## =============================================================
## PROCESS .TAR ARCHIVE AND BUILD COUNT MATRIX
## =============================================================

tar_file <- file.path(dir_raw, ACC, paste0(ACC, "_RAW.tar"))

if (!file.exists(tar_file)) {
  cat("\nDownloading supplementary files...\n")
  getGEOSuppFiles(ACC, baseDir = dir_raw)
}

## Extract the tar archive
untar_dir <- file.path(dir_raw, ACC, "untarred")
dir.create(untar_dir, showWarnings = FALSE)
cat("\nExtracting tar archive...\n")
untar(tar_file, exdir = untar_dir)

## Read and merge all extracted count files
raw_files <- list.files(untar_dir, pattern = "\\.txt\\.gz$|\\.tsv\\.gz$|\\.counts\\.gz$|\\.txt$", full.names = TRUE)

if (length(raw_files) > 0) {
  cat("Found", length(raw_files), "count files. Merging into a single matrix...\n")
  
  count_list <- list()
  for (f in raw_files) {
    # Read each sample file (assuming 2 columns: Gene ID and Count)
    d <- read.table(f, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
    
    # Safely extract GSM ID from the filename
    gsm_id <- str_extract(basename(f), "GSM[0-9]+")
    
    if (!is.na(gsm_id)) {
      colnames(d) <- gsm_id
      count_list[[gsm_id]] <- d
    } else {
      # Fallback if GSM is not explicitly in the filename
      colnames(d) <- basename(f)
      count_list[[basename(f)]] <- d
    }
  }
  
  # Bind all columns together
  raw_counts <- do.call(cbind, count_list)
  
  # Clean up: Remove HTSeq alignment metrics if present (typically start with '__')
  raw_counts <- raw_counts[!grepl("^__", rownames(raw_counts)), ]
} else {
  stop("No count files found inside the .tar archive. Please check the archive contents.")
}

cat("Raw counts assembled successfully:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

## =============================================================
## BUILD METADATA AND FILTER GROUPS
## =============================================================

metadata <- data.frame(
  sample_id = rownames(pheno),
  group     = unname(GROUP_MAP[as.character(pheno[[GROUP_COL]])]),
  stringsAsFactors = FALSE
)
stopifnot(!any(is.na(metadata$group)))

## Remove excluded samples (e.g., DM, TB_DM)
metadata <- metadata[metadata$group != "Exclude", ]
metadata$group <- factor(metadata$group, levels = c("Healthy_Control", "TB_Only"))
rownames(metadata) <- metadata$sample_id

## Intersect to keep only samples present in both counts and metadata
common <- intersect(colnames(raw_counts), metadata$sample_id)

if (length(common) == 0) {
    cat("\nWARNING: Column names in raw_counts do not perfectly match sample_ids in metadata.\n")
    
    # Fallback: Attempt force-mapping if lengths are exactly the same
    if(ncol(raw_counts) == nrow(pheno)) {
        cat("Attempting to force-map colnames to pheno rownames...\n")
        colnames(raw_counts) <- rownames(pheno)
        common <- intersect(colnames(raw_counts), metadata$sample_id)
    } else {
        stop("Cannot resolve sample name mismatch. Manual inspection of the matrices is required.")
    }
}

raw_counts <- raw_counts[, common]
metadata   <- metadata[common, ]

cat("\n--- Final Sample Count per Group ---\n")
print(table(metadata$group))

write.csv(metadata, file.path(dir_meta, paste0(ACC, "_sample_metadata.csv")), row.names = FALSE)
plot_sample_design(metadata, out_png = file.path(dir_fig, paste0(ACC, "_sample_design.png")),
                   title = paste(ACC, "Sample Design"))

## =============================================================
## FILTERING
## =============================================================
cat("\nFiltering low expression genes...\n")
write.csv(raw_counts, file.path(dir_raw, paste0(ACC, "_raw_count_matrix.csv")))
filt <- filter_low_expression(raw_counts, min_count = 10,
                              min_samples = min(table(metadata$group)))
write.csv(filt$filtered, file.path(dir_filt, paste0(ACC, "_filtered_count_matrix.csv")))
plot_filtering_summary(filt$n_before, filt$n_after,
                       out_png = file.path(dir_fig, paste0(ACC, "_filtering_summary.png")),
                       title = paste(ACC, "Gene Filtering"))

## =============================================================
## DESEQ2 + VST NORMALIZATION + PCA
## =============================================================
cat("\nRunning DESeq2 and VST Normalization...\n")
dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
dds <- DESeq(dds)

vst_mat <- run_vst(dds)
write.csv(vst_mat, file.path(dir_norm, paste0(ACC, "_vst_normalized_matrix.csv")))
plot_pca(vst_mat, metadata, out_png = file.path(dir_fig, paste0(ACC, "_PCA.png")),
         title = paste(ACC, "PCA (Healthy Control vs TB Only)"))

## =============================================================
## SAMPLE CORRELATION & OUTLIER CHECK
## =============================================================
cat("\nGenerating Sample Correlation Heatmap...\n")
sample_correlation(vst_mat, metadata,
                   out_matrix_csv = file.path(dir_corr, paste0(ACC, "_sample_correlation_matrix.csv")),
                   out_png = file.path(dir_fig, paste0(ACC, "_sample_correlation_heatmap.png")),
                   title = paste(ACC, "Sample Correlation"))

corr_mat <- cor(vst_mat, method = "pearson")
avg_corr <- sort(rowMeans(corr_mat))
cat("\n--- Lowest Mean Sample Correlation (Possible Outliers) ---\n")
print(head(avg_corr, 6))

## =============================================================
## DIFFERENTIAL EXPRESSION ANALYSIS (DEG)
## =============================================================
cat("\nRunning DEG Analysis: TB_Only vs Healthy_Control...\n")
res_tb <- run_deg(dds, contrast = c("group", "TB_Only", "Healthy_Control"),
                  out_csv = file.path(dir_deg, paste0(ACC, "_HealthyControl_vs_TBOnly_DEG.csv")),
                  lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

volcano_plot(res_tb, out_png = file.path(dir_fig, paste0(ACC, "_HealthyControl_vs_TBOnly_volcano.png")),
             title = paste(ACC, ": Healthy Control vs TB Only"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

## =============================================================
## SIGNIFICANT DEGS & HEATMAPS
## =============================================================
sig_tb <- get_significant_degs(res_tb,
            out_csv = file.path(dir_sig, paste0(ACC, "_HealthyControl_vs_TBOnly_sigDEGs.csv")),
            lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)

if (nrow(sig_tb) > 1) {
  cat("\nGenerating Heatmap for Significant DEGs...\n")
  sig_deg_heatmap(vst_mat, sig_tb$gene, metadata,
                  out_png = file.path(dir_fig, paste0(ACC, "_HealthyControl_vs_TBOnly_DEG_heatmap.png")),
                  title = paste(ACC, ": Top DEGs, Healthy vs TB"))
}

## =============================================================
## PREPARE WGCNA INPUTS
## =============================================================
cat("\nExporting WGCNA Inputs...\n")
export_wgcna_expression(vst_mat, file.path(dir_wgcna, paste0(ACC, "_WGCNA_expression_matrix.csv")))
export_wgcna_traits(metadata,
  trait_cols = list(Healthy_Control = as.integer(metadata$group == "Healthy_Control"),
                    TB_Only         = as.integer(metadata$group == "TB_Only")),
  out_csv = file.path(dir_wgcna, paste0(ACC, "_WGCNA_trait_file.csv")))

## =============================================================
## SUMMARY
## =============================================================
cat("\n========== DEG COUNT SUMMARY ==========\n")
cat("Healthy vs TB_Only : ", nrow(sig_tb), " significant DEGs\n", sep = "")
cat("(Threshold: padj < ", PADJ_TH, ", |log2FC| >= ", LFC_TH, ")\n", sep = "")
cat("=======================================\n")
