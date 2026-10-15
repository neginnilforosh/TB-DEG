## GSE114192 - TB-Diabetes framework.
## The original GSE114192_run.R maps TB_DM, DM_only, TB_IH and IH to "Exclude" and keeps only the
## 82 Healthy_Control / TB_Only samples, and its outputs are what every later step depends on.
## This script rebuilds the count matrix from the same archive WITHOUT that exclusion and writes
## everything into GSE114192/<RUN>/, so nothing already produced is overwritten.
##
## Run it from scripts/deg_pipeline (or the tutorial folder that holds 00_functions.R):
##   Rscript GSE114192_comorbidity_run.R
## Then point config.R at the new run to continue with the network pipeline.

## ---- CONFIG ----
ACC       <- "GSE114192"
RUN       <- "TB_DM"                       # output subfolder
LFC_TH    <- 1
PADJ_TH   <- 0.05
GROUP_COL <- "disease state (disease_category):ch1"

## every group is kept here; the comparisons below decide what is actually used
GROUP_MAP <- c("Healthy_Control" = "Healthy_Control",
               "TB_only"         = "TB_Only",
               "TB_DM"           = "TB_DM",
               "DM_only"         = "DM_only",
               "TB_IH"           = "Exclude",     # intermediate hyperglycaemia: small and hard to interpret
               "IH"              = "Exclude")

## case first, control second -- a positive log2FC means higher in the case group
COMPARISONS <- list(
  c("TB_DM", "DM_only"),           # MAIN: what TB adds on a diabetic background
  c("TB_DM", "TB_Only"),           # what diabetes adds on top of TB
  c("TB_DM", "Healthy_Control"),   # the full comorbid state vs healthy
  c("DM_only", "Healthy_Control")  # diabetes alone, as a reference for the others
)
## Expression ~ TB + DM + TB:DM over all four groups. The interaction term finds genes whose
## behaviour in the comorbid state is not just the sum of the two separate effects, which no
## pairwise contrast can show. It needs all four groups present with enough samples each.
RUN_INTERACTION <- TRUE
MIN_PER_GROUP   <- 3     # a group below this is dropped from a comparison
MIN_FOR_WGCNA   <- 15    # below this, a per-context WGCNA is not worth running

suppressPackageStartupMessages({
  library(DESeq2); library(GEOquery); library(stringr)
})
SCRIPT_DIR <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(local({ d <- SCRIPT_DIR
  hit <- list.files(dirname(d), pattern = "^00_functions\\.R$", recursive = TRUE, full.names = TRUE)
  if (!length(hit)) hit <- list.files(d, pattern = "^00_functions\\.R$", full.names = TRUE)
  if (!length(hit)) stop("00_functions.R not found"); hit[1] }))

BASE_DIR <- local({ d <- SCRIPT_DIR
  while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })
dir_raw  <- file.path(BASE_DIR, "01_raw_counts")          # read from the shared raw folder
OUT      <- file.path(BASE_DIR, RUN)                      # everything else written here
dirs <- list(meta = "02_metadata", filt = "03_filtered_counts", norm = "04_normalized",
             corr = "05_sample_correlation", deg = "06_deg_results", sig = "07_significant_degs",
             wgcna = "09_wgcna_input", fig = "figures")
for (d in dirs) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
P <- function(key, ...) file.path(OUT, dirs[[key]], paste0(...))
cat(">>> writing to", OUT, "-- the existing", ACC, "outputs are not touched\n")

## ---- counts from the archive (all samples, no exclusion) ----
gse   <- getGEO(ACC, GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)
tar_file <- file.path(dir_raw, ACC, paste0(ACC, "_RAW.tar"))
if (!file.exists(tar_file)) { cat(">>> downloading supplementary files...\n"); getGEOSuppFiles(ACC, baseDir = dir_raw) }
untar_dir <- file.path(dir_raw, ACC, "untarred")
dir.create(untar_dir, showWarnings = FALSE)
untar(tar_file, exdir = untar_dir)
raw_files <- list.files(untar_dir, pattern = "\\.txt\\.gz$|\\.tsv\\.gz$|\\.counts\\.gz$|\\.txt$", full.names = TRUE)
if (!length(raw_files)) stop("No count files inside the archive: ", tar_file)
count_list <- list()
for (f in raw_files) {
  d <- read.table(f, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
  id <- str_extract(basename(f), "GSM[0-9]+"); if (is.na(id)) id <- basename(f)
  colnames(d) <- id; count_list[[id]] <- d
}
raw_counts <- do.call(cbind, count_list)
raw_counts <- raw_counts[!grepl("^__", rownames(raw_counts)), ]
cat(">>> raw counts:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

## ---- metadata with the comorbid groups kept ----
metadata <- data.frame(sample_id = rownames(pheno),
                       group = unname(GROUP_MAP[as.character(pheno[[GROUP_COL]])]), stringsAsFactors = FALSE)
if (anyNA(metadata$group)) stop("Unmapped group label(s): ",
    paste(unique(as.character(pheno[[GROUP_COL]])[is.na(metadata$group)]), collapse = ", "))
metadata <- metadata[metadata$group != "Exclude", ]
rownames(metadata) <- metadata$sample_id
common <- intersect(colnames(raw_counts), metadata$sample_id)
## matched by name only; no positional fallback (pairing columns with pData() by position can mislabel samples)
if (!length(common)) stop("Sample names in the counts and the metadata do not match.")
raw_counts <- raw_counts[, common]; metadata <- metadata[common, ]
metadata$group <- factor(metadata$group)
cat(">>> groups kept:\n"); print(table(metadata$group))
write.csv(metadata, P("meta", ACC, "_", RUN, "_sample_metadata.csv"), row.names = FALSE)
plot_sample_design(metadata, out_png = P("fig", ACC, "_", RUN, "_sample_design.png"),
                   title = paste(ACC, RUN, "sample design"))

## ---- filter, normalise ----
filt <- filter_low_expression(raw_counts, min_count = 10, min_samples = min(table(metadata$group)))
write.csv(filt$filtered, P("filt", ACC, "_", RUN, "_filtered_count_matrix.csv"))
dds <- DESeqDataSetFromMatrix(filt$filtered, metadata, design = ~ group)
dds <- DESeq(dds)
vst_mat <- run_vst(dds)
plot_pca(vst_mat, metadata, out_png = P("fig", ACC, "_", RUN, "_PCA.png"), title = paste(ACC, RUN, "PCA"))
write.csv(vst_mat, P("norm", ACC, "_", RUN, "_vst_normalized_matrix.csv"))
sample_correlation(vst_mat, metadata,
                   out_matrix_csv = P("corr", ACC, "_", RUN, "_sample_correlation_matrix.csv"),
                   out_png = P("fig", ACC, "_", RUN, "_sample_correlation_heatmap.png"),
                   title = paste(ACC, RUN, "sample correlation"))

## ---- the comparisons ----
summary_rows <- list()
for (cmp in COMPARISONS) {
  case <- cmp[1]; ctrl <- cmp[2]
  if (!all(c(case, ctrl) %in% levels(metadata$group))) { cat("   skipping", case, "vs", ctrl, "(group missing)\n"); next }
  n_case <- sum(metadata$group == case); n_ctrl <- sum(metadata$group == ctrl)
  tag <- paste0(ctrl, "_vs_", case)
  cat("\n>>>", case, "(", n_case, ") vs", ctrl, "(", n_ctrl, ")\n")
  if (min(n_case, n_ctrl) < MIN_PER_GROUP) { cat("   fewer than", MIN_PER_GROUP, "samples in a group, skipped\n"); next }
  ## One shrinkage method for every contrast: apeglm needs the contrast as a model coefficient, so the
  ## control is made the reference level and only the Wald test is refitted (dispersions do not depend on
  ## the reference level). Without this, only contrasts against the first factor level were shrunk.
  dds_c <- dds
  dds_c$group <- relevel(dds_c$group, ref = ctrl)
  dds_c <- nbinomWaldTest(dds_c)
  res <- run_deg(dds_c, contrast = c("group", case, ctrl),
                 out_csv = P("deg", ACC, "_", tag, "_DEG.csv"), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
  volcano_plot(res, out_png = P("fig", ACC, "_", tag, "_volcano.png"),
               title = paste0(ACC, ": ", ctrl, " vs ", case), lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
  sig <- get_significant_degs(res, out_csv = P("sig", ACC, "_", tag, "_sigDEGs.csv"),
                              lfc_thresh = LFC_TH, padj_thresh = PADJ_TH)
  if (nrow(sig) > 1)
    sig_deg_heatmap(vst_mat, sig$gene, metadata[metadata$group %in% c(case, ctrl), ],
                    out_png = P("fig", ACC, "_", tag, "_DEG_heatmap.png"), title = paste(ACC, tag))
  summary_rows[[tag]] <- data.frame(Comparison = tag, Case = case, N_case = n_case,
                                    Control = ctrl, N_control = n_ctrl, N_significant = nrow(sig))
}
summ <- do.call(rbind, summary_rows)
write.csv(summ, P("deg", ACC, "_", RUN, "_comparison_summary.csv"), row.names = FALSE)

## ---- interaction model: Expression ~ TB + DM + TB:DM ----
if (RUN_INTERACTION) {
  four <- c("Healthy_Control", "TB_Only", "DM_only", "TB_DM")
  if (!all(four %in% levels(metadata$group))) {
    cat("\n>>> interaction model skipped: needs all four groups, missing",
        paste(setdiff(four, levels(metadata$group)), collapse = ", "), "\n")
  } else if (min(table(metadata$group)[four]) < MIN_PER_GROUP) {
    cat("\n>>> interaction model skipped: smallest group has only",
        min(table(metadata$group)[four]), "samples\n")
  } else {
    cat("\n>>> interaction model: ~ TB + DM + TB:DM\n")
    md <- metadata[metadata$group %in% four, , drop = FALSE]
    md$TB <- factor(ifelse(md$group %in% c("TB_Only", "TB_DM"), "yes", "no"), levels = c("no", "yes"))
    md$DM <- factor(ifelse(md$group %in% c("DM_only", "TB_DM"), "yes", "no"), levels = c("no", "yes"))
    print(table(md$TB, md$DM, dnn = c("TB", "DM")))
    dds_i <- DESeqDataSetFromMatrix(filt$filtered[, rownames(md)], md, design = ~ TB + DM + TB:DM)
    dds_i <- DESeq(dds_i)
    cat("    coefficients:", paste(resultsNames(dds_i), collapse = ", "), "\n")
    for (cf in c(main_TB = "TB_yes_vs_no", main_DM = "DM_yes_vs_no", interaction = "TByes.DMyes")) {
      if (!cf %in% resultsNames(dds_i)) { cat("    !! coefficient", cf, "not found, skipped\n"); next }
      r <- as.data.frame(results(dds_i, name = cf)); r$gene <- rownames(r)
      r <- r[order(r$padj), c("gene", setdiff(names(r), "gene"))]
      tag <- names(which(c(main_TB = "TB_yes_vs_no", main_DM = "DM_yes_vs_no", interaction = "TByes.DMyes") == cf))
      write.csv(r, P("deg", ACC, "_interaction_", tag, ".csv"), row.names = FALSE)
      cat("    ", tag, ": significant (padj <", PADJ_TH, "):", sum(r$padj < PADJ_TH, na.rm = TRUE), "\n")
    }
    cat("    A gene significant for the interaction behaves differently in TB-DM than TB and DM alone\n")
    cat("    would predict -- that is the point of this model.\n")
  }
}

## ---- is a per-context WGCNA worth running? ----
cat("\n>>> sample numbers per group (a per-context WGCNA needs roughly", MIN_FOR_WGCNA, "or more):\n")
print(table(metadata$group))
small <- names(which(table(metadata$group) < MIN_FOR_WGCNA))
if (length(small)) cat(">>> below that threshold:", paste(small, collapse = ", "),
                       "-- WGCNA on a comparison involving these will be unstable; say so rather than reporting modules.\n")

## ---- WGCNA inputs: one trait column per group, so any comparison can be selected later ----
export_wgcna_expression(vst_mat, P("wgcna", ACC, "_WGCNA_expression_matrix.csv"))
trait_cols <- setNames(lapply(levels(metadata$group), function(g) as.integer(metadata$group == g)),
                       levels(metadata$group))
export_wgcna_traits(metadata, trait_cols = trait_cols, out_csv = P("wgcna", ACC, "_WGCNA_trait_file.csv"))

cat("\n========== SUMMARY ==========\n"); print(summ, row.names = FALSE)
cat("\n>>> DONE. Everything in", OUT, "\n")
cat(">>> To continue with the network pipeline, add a preset to config.R, e.g.\n")
cat('    GSE114192_TBDM = list(ACC = "GSE114192", CONTROL = "TB_Only", CASE = "TB_DM",\n')
cat('                          DEG_NAME = "TBOnly_vs_TB_DM", RUN = "TB_DM",\n')
cat('                          MODULE_NAMES = character(0), MODULE_PLOT_COLORS = NULL)\n')
