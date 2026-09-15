## =============================================================
## Fix: generate the missing sample-correlation heatmap for GSE99374
## Standalone — reads your already-saved VST matrix + metadata,
## does not need the rest of the pipeline's session objects.
## =============================================================

setwd("/Users/negin/Desktop/TB_DEG_pipeline/scripts")   # <-- edit this
source("00_functions.R")

BASE_DIR <- file.path("..", "GSE99374")

## ---- Load what you already produced ----
vst_mat  <- as.matrix(read.csv(file.path(BASE_DIR, "04_normalized", "GSE99374_vst_normalized_matrix.csv"),
                                row.names = 1, check.names = FALSE))
metadata <- read.csv(file.path(BASE_DIR, "02_metadata", "GSE99374_sample_metadata.csv"))
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(vst_mat), ]   # enforce same order as vst_mat columns

## ---- Generate the heatmap (this ALSO re-writes the correlation CSV,
## which is fine — it'll be identical to what you already have) ----
sample_correlation(
  vst_mat, metadata, group_col = "group",
  out_matrix_csv = file.path(BASE_DIR, "05_sample_correlation", "GSE99374_sample_correlation_matrix.csv"),
  out_png        = file.path(BASE_DIR, "figures", "GSE99374_sample_correlation_heatmap.png"),
  title = "GSE99374 sample correlation"
)

message("Done — check GSE99374/figures/GSE99374_sample_correlation_heatmap.png")

