## STEP 13: cell-type profile of module genes and drug targets, using the 4-week M. tuberculosis granuloma
## single-cell dataset (Broad Single Cell Portal SCP1749: macaque, 11 cell types, 2 animals).
## Data: the four files of "data files.zip" in RajarshiRay25/Single-Cell-Analysis---Tuberculosis,
## placed in external_data/SCP1749/.
## Method: pseudo-bulk CP10K per cell type; per gene, the share of its expression in each cell type
## (uniform = 9%). Lung granuloma, not blood: read it as the cell-type origin of a gene.

if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package Matrix (ships with standard R; or install.packages(\"Matrix\")).")

## ---- CONFIG ----
ACC             <- "GSE114192"
SC_SUBDIR       <- file.path("external_data", "SCP1749")
CELLTYPE_COL    <- "CellTypeAnnotations"
DETECT_MIN_CP10K <- 1
SPECIFIC_SHARE  <- 0.5
MIN_TARGETS_SC  <- 3        # a drug needs at least this many detected targets to get a cell-type profile
MIN_GENES_TEST  <- 10       # min genes per group for the module-level tests
IMMUNE_TYPES    <- c("Macrophage", "T", "Neutrophil", "Mast", "B", "pDC", "Plasma")   # for the blood-relevant drug view
DIFFUSE_FACTOR  <- 1.5      # a top cell type must hold >= this x the uniform share, otherwise the drug is labelled "diffuse"

## ---- paths ----
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()
ROOT_DIR   <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR   <- file.path(ROOT_DIR, ACC)
SC_DIR     <- file.path(ROOT_DIR, SC_SUBDIR)
dir_out    <- file.path(BASE_DIR, "18_singlecell"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
need <- file.path(SC_DIR, c("4Week_countsmatrix.mtx", "4Week_features.tsv", "4Week_barcodes.tsv", "metadata.txt"))
if (!all(file.exists(need))) stop("Single-cell files not found in ", SC_DIR, ":\n  missing: ", paste(basename(need[!file.exists(need)]), collapse = ", "),
                                  "\n  (see the header of this script for where to get them)")

## ---- helpers ----
# genes x cell-type CP10K from a sparse counts matrix and a cell-type vector
pseudobulk_cp10k <- function(X, celltype) {
  types <- names(sort(table(celltype), decreasing = TRUE))
  ind   <- Matrix::sparseMatrix(i = seq_along(celltype), j = match(celltype, types), x = 1, dims = c(length(celltype), length(types)),
                                dimnames = list(NULL, types))
  sums  <- as.matrix(X %*% ind)
  list(cp10k = sweep(sums, 2, colSums(sums), "/") * 1e4, n_cells = as.integer(table(celltype)[types]))
}
# per-gene share / top cell type / detected / specific
celltype_calls <- function(cp10k, detect_min, specific_share) {
  tot <- rowSums(cp10k)
  keep <- tot > 0
  cp <- cp10k[keep, , drop = FALSE]
  share <- cp / rowSums(cp)
  top <- max.col(share, ties.method = "first")
  out <- data.frame(Gene = rownames(cp), top_celltype = colnames(cp)[top], top_share = share[cbind(seq_len(nrow(cp)), top)],
                    top_cp10k = cp[cbind(seq_len(nrow(cp)), top)], stringsAsFactors = FALSE)
  out$detected <- out$top_cp10k >= detect_min
  out$specific <- out$detected & out$top_share >= specific_share
  list(calls = out, share = share)
}
# mean share per group (rows of `share` selected by `idx_list`)
mean_share_by_group <- function(share, groups) {
  do.call(rbind, lapply(names(groups), function(g) {
    idx <- groups[[g]]; if (length(idx) == 0) return(NULL)
    data.frame(Group = g, N_genes = length(idx), t(colMeans(share[idx, , drop = FALSE])), check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

## ---- 1. read the single-cell data ----
cat(">>> Reading the single-cell data (the 175 MB matrix takes a minute or two)...\n")
feat <- read.delim(file.path(SC_DIR, "4Week_features.tsv"), header = FALSE, stringsAsFactors = FALSE)[[1]]
bar  <- read.delim(file.path(SC_DIR, "4Week_barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)[[1]]
meta <- read.delim(file.path(SC_DIR, "metadata.txt"), stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$NAME != "TYPE", ]                       # Single Cell Portal metadata has a second "TYPE" row
if (!CELLTYPE_COL %in% names(meta)) stop("Column '", CELLTYPE_COL, "' not in metadata.txt; columns are: ", paste(names(meta), collapse = ", "))
ct <- meta[[CELLTYPE_COL]][match(bar, meta$NAME)]
if (anyNA(ct)) stop(sum(is.na(ct)), " barcodes have no cell-type annotation in metadata.txt")
X <- Matrix::readMM(file.path(SC_DIR, "4Week_countsmatrix.mtx"))
X <- as(X, "CsparseMatrix")
if (nrow(X) != length(feat) || ncol(X) != length(bar)) stop("Matrix is ", nrow(X), " x ", ncol(X), " but features/barcodes are ", length(feat), " / ", length(bar))
rownames(X) <- feat
cat(">>>", nrow(X), "genes x", ncol(X), "cells;", length(unique(ct)), "cell types:\n"); print(table(ct))

pb <- pseudobulk_cp10k(X, ct)
cp10k <- pb$cp10k
cc <- celltype_calls(cp10k, DETECT_MIN_CP10K, SPECIFIC_SHARE)
cat(">>> Genes detected (>=", DETECT_MIN_CP10K, "CP10K in their best cell type):", sum(cc$calls$detected), "of", nrow(cc$calls), "\n")

## sanity check with canonical markers: each should peak in the expected cell type
chk <- c(CD3D = "T", MS4A1 = "B", CSF3R = "Neutrophil", C1QB = "Macrophage", SFTPC = "T2P", COL1A1 = "Fibroblast", PECAM1 = "Endothelial", JCHAIN = "Plasma")
chk <- chk[names(chk) %in% cc$calls$Gene]
got <- cc$calls$top_celltype[match(names(chk), cc$calls$Gene)]
cat(">>> Marker sanity check:", paste0(names(chk), "->", got, ifelse(got == chk, "(ok)", paste0("(EXPECTED ", chk, ")")), collapse = "  "), "\n")
if (mean(got == chk) < 0.75) stop("Canonical markers do not land in the expected cell types -- wrong annotation column or wrong file pairing?")

## ---- 2. module genes -> cell types ----
gt <- read.csv(file.path(BASE_DIR, "16_final_tables", paste0(ACC, "_FINAL_GeneTable.csv")), stringsAsFactors = FALSE)
gt <- gt[!is.na(gt$Gene) & nzchar(gt$Gene) & !duplicated(gt$Gene), ]
gt$in_singlecell <- toupper(gt$Gene) %in% toupper(rownames(cp10k))
cat(">>> Module genes with a gene symbol:", nrow(gt), "| present in the macaque feature list:", sum(gt$in_singlecell),
    "(", paste(names(tapply(gt$in_singlecell, gt$Module, mean)), round(100 * tapply(gt$in_singlecell, gt$Module, mean)), "%", collapse = ", "), ")\n")

rn <- toupper(rownames(cp10k)); cp_mod <- cp10k[match(toupper(gt$Gene[gt$in_singlecell]), rn), , drop = FALSE]
rownames(cp_mod) <- gt$Gene[gt$in_singlecell]
mc <- celltype_calls(cp_mod, DETECT_MIN_CP10K, SPECIFIC_SHARE)
mg <- merge(gt[, c("Gene", "Module", "log2FC", "FDR")], mc$calls, by = "Gene")
mg$Direction <- ifelse(mg$log2FC > 0, "higher_in_TB", "lower_in_TB")
sh <- as.data.frame(mc$share); names(sh) <- paste0("share_", names(sh)); sh$Gene <- rownames(mc$share)
mg <- merge(mg, sh, by = "Gene")
mg <- mg[order(mg$Module, -mg$top_share), ]
write.csv(mg, file.path(dir_out, paste0(ACC, "_module_genes_celltype.csv")), row.names = FALSE)
cat(">>> Per module: genes assessed / detected / cell-type-specific:\n")
print(do.call(rbind, lapply(split(mg, mg$Module), function(x) data.frame(Module = x$Module[1], assessed = nrow(x), detected = sum(x$detected), specific = sum(x$specific)))), row.names = FALSE)

## ---- module x direction profile (mean share over DETECTED genes) ----
det <- mg[mg$detected, ]
shm <- as.matrix(det[, grep("^share_", names(det))]); colnames(shm) <- sub("^share_", "", colnames(shm)); rownames(shm) <- det$Gene
groups <- list()
for (m in sort(unique(det$Module))) for (dr in c("lower_in_TB", "higher_in_TB")) {
  idx <- which(det$Module == m & det$Direction == dr)
  if (length(idx) >= MIN_GENES_TEST) groups[[paste(m, dr, sep = " | ")]] <- idx
}
prof <- mean_share_by_group(shm, groups)
prof[, -(1:2)] <- round(100 * prof[, -(1:2)], 1)
write.csv(prof, file.path(dir_out, paste0(ACC, "_module_celltype_profile.csv")), row.names = FALSE)
cat("\n>>> Mean cell-type share (%) of detected genes, per module and direction (uniform expectation = ", round(100 / ncol(shm), 1), "%):\n", sep = "")
print(prof, row.names = FALSE)

## ---- lower vs higher in TB within each module, per cell type (Wilcoxon; descriptive only) ----
tests <- list()
for (m in sort(unique(det$Module))) {
  lo <- which(det$Module == m & det$Direction == "lower_in_TB"); hi <- which(det$Module == m & det$Direction == "higher_in_TB")
  if (length(lo) < MIN_GENES_TEST || length(hi) < MIN_GENES_TEST) next
  for (ctype in colnames(shm)) {
    w <- suppressWarnings(stats::wilcox.test(shm[lo, ctype], shm[hi, ctype]))
    tests[[length(tests) + 1]] <- data.frame(Module = m, CellType = ctype, N_lower = length(lo), N_higher = length(hi),
                                             Share_lower_pct = round(100 * mean(shm[lo, ctype]), 1), Share_higher_pct = round(100 * mean(shm[hi, ctype]), 1),
                                             P_wilcoxon = w$p.value, stringsAsFactors = FALSE)
  }
}
if (length(tests)) {
  tests <- do.call(rbind, tests); tests <- tests[order(tests$P_wilcoxon), ]
  write.csv(tests, file.path(dir_out, paste0(ACC, "_module_celltype_tests.csv")), row.names = FALSE)
  cat("\n>>> Strongest lower-vs-higher differences (descriptive; genes are not independent):\n"); print(utils::head(tests, 8), row.names = FALSE)
}

## ---- heatmap of the profile ----
tryCatch({
  mat <- as.matrix(prof[, -(1:2)]); rownames(mat) <- prof$Group
  grDevices::pdf(file.path(dir_out, paste0(ACC, "_module_celltype_heatmap.pdf")), width = 9, height = 5)
  op <- graphics::par(mar = c(7, 13, 3, 2))
  graphics::image(t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::colorRampPalette(c("white", "#c0392b"))(50),
                  main = "Mean cell-type share of detected module genes (%)")
  graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = rev(rownames(mat)), las = 1, cex.axis = 0.8)
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat)))
    graphics::text((j - 1) / max(ncol(mat) - 1, 1), (nrow(mat) - i) / max(nrow(mat) - 1, 1), mat[i, j], cex = 0.65)
  graphics::par(op); grDevices::dev.off()
}, error = function(e) cat("  (heatmap skipped:", conditionMessage(e), ")\n"))

## ---- 3. drug targets -> cell types ----
tg_file  <- file.path(BASE_DIR, "14_drug_targets", paste0(ACC, "_drug_targets_ALL.csv"))
fin_file <- file.path(BASE_DIR, "16_final_tables", paste0(ACC, "_FINAL_DrugTable.csv"))
if (file.exists(tg_file) && file.exists(fin_file)) {
  tg  <- read.csv(tg_file, stringsAsFactors = FALSE); fin <- read.csv(fin_file, stringsAsFactors = FALSE)
  det_genes <- cc$calls[cc$calls$detected, ]
  share_all <- cc$share[det_genes$Gene, , drop = FALSE]; rownames(share_all) <- toupper(rownames(share_all))
  imm <- intersect(IMMUNE_TYPES, colnames(share_all))
  pick_top <- function(m, factor) {                       # m: named mean-share vector summing to ~1
    o <- order(m, decreasing = TRUE)
    list(name = if (m[o[1]] >= factor / length(m)) names(m)[o[1]] else "diffuse", share = round(100 * m[o[1]], 1), order = o)
  }
  rows <- lapply(unique(fin$Drug), function(d) {
    tgs <- unique(toupper(tg$Target_Gene[tg$Drug == d])); tgs <- tgs[!is.na(tgs) & nzchar(tgs)]
    used <- intersect(tgs, rownames(share_all))
    r <- data.frame(Drug = d, N_known_targets = length(tgs), N_targets_detected_sc = length(used), stringsAsFactors = FALSE)
    r$SC_Top_CellType <- NA_character_; r$SC_Top_Share_pct <- NA_real_; r$SC_Top_Immune_CellType <- NA_character_
    r$SC_Top_Immune_Share_pct <- NA_real_; r$SC_Profile <- NA_character_
    if (length(used) >= MIN_TARGETS_SC) {
      m  <- colMeans(share_all[used, , drop = FALSE]); tp <- pick_top(m, DIFFUSE_FACTOR)
      r$SC_Top_CellType <- tp$name; r$SC_Top_Share_pct <- tp$share
      r$SC_Profile <- paste0(names(m)[tp$order[1:3]], " ", round(100 * m[tp$order[1:3]]), "%", collapse = "; ")
      si <- share_all[used, imm, drop = FALSE]; si <- si[rowSums(si) > 0, , drop = FALSE]      # re-split each target's expression among IMMUNE cell types only
      if (nrow(si) >= 1) { mi <- colMeans(si / rowSums(si)); ti <- pick_top(mi, DIFFUSE_FACTOR)
        r$SC_Top_Immune_CellType <- ti$name; r$SC_Top_Immune_Share_pct <- ti$share }
    }
    r
  })
  dr <- do.call(rbind, rows)
  write.csv(dr, file.path(dir_out, paste0(ACC, "_drug_celltype.csv")), row.names = FALSE)
  cat("\n>>> Drugs with a cell-type profile (>=", MIN_TARGETS_SC, "detected targets):", sum(!is.na(dr$SC_Top_CellType)), "of", nrow(dr), "\n")
  cat(">>> Top cell type over all 11 types:\n"); print(table(dr$SC_Top_CellType, useNA = "ifany"))
  cat(">>> Top IMMUNE cell type (immune cell types only):\n"); print(table(dr$SC_Top_Immune_CellType, useNA = "ifany"))
} else cat("\n>>> (drug files not found -- drug part skipped)\n")

cat("\n>>> DONE. Outputs in", dir_out, "\n")
cat(">>> Next: re-run WGCNA_final_tables_step10.R -- it merges the SC_* columns into the final drug table.\n")
