## =============================================================
## Test whether TARGET_MODULE (e.g. "blue": 1377 genes, no
## significant GO/KEGG/Reactome term) is an over-merging artifact.
##
## Does NOT rebuild the network/TOM from scratch -- reuses
## the pre-merge module calls already saved by step 1
## (net$unmergedColors inside <ACC>_WGCNA_workspace.rds) and just
## re-cuts the merge threshold. Fast.
## =============================================================

library(WGCNA)
options(stringsAsFactors = FALSE)

## ---- CONFIG ----
ACC             <- "GSE114192"
TRAIT_COL       <- "TB_Only"
TARGET_MODULE   <- "blue"                 # the module you're trying to split
NEW_CUT_HEIGHTS <- c(0.10, 0.15, 0.20)    # compare against the original 0.25

## ---- locate script dir (same as step 1/step 2) ----
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
SCRIPT_DIR  <- get_script_dir()
BASE_DIR    <- file.path(SCRIPT_DIR, "..", ACC)
dir_results <- file.path(BASE_DIR, "10_wgcna_results")

## ---- load the saved workspace from step 1 ----
ws_file <- file.path(dir_results, paste0(ACC, "_WGCNA_workspace.rds"))
if (!file.exists(ws_file)) stop("Missing: ", ws_file, " -- run WGCNA_TBmodules_step1.R first.")
ws <- readRDS(ws_file)
datExpr   <- ws$datExpr
datTraits <- ws$datTraits
net       <- ws$net

if (is.null(net$unmergedColors)) stop("This workspace has no net$unmergedColors -- re-run step 1 (it's part of blockwiseModules' normal output, so this shouldn't happen unless the object was saved differently).")
unmergedColors <- labels2colors(net$unmergedColors)

n_before <- sum(ws$moduleColors == TARGET_MODULE)
cat(">>> '", TARGET_MODULE, "' currently has", n_before, "genes (after the original merge).\n\n")

for (ch in NEW_CUT_HEIGHTS) {
  cat("==== mergeCutHeight =", ch, "====\n")
  merged    <- mergeCloseModules(datExpr, unmergedColors, cutHeight = ch, verbose = 0)
  newColors <- merged$colors

  old_idx <- which(ws$moduleColors == TARGET_MODULE)
  derived <- sort(unique(newColors[old_idx]))
  derived <- setdiff(derived, "grey")

  if (length(derived) <= 1) {
    cat("  '", TARGET_MODULE, "' stays as ONE module at this cutHeight (no split).\n\n")
    next
  }

  cat("  '", TARGET_MODULE, "' splits into", length(derived), "modules:", paste(derived, collapse = ", "), "\n")
  MEs_new <- orderMEs(moduleEigengenes(datExpr, newColors)$eigengenes)
  tb_numeric <- as.numeric(datTraits[[TRAIT_COL]])

  for (d in derived) {
    me_col <- paste0("ME", d)
    n_genes <- sum(newColors == d)
    if (me_col %in% names(MEs_new)) {
      r <- cor(MEs_new[[me_col]], tb_numeric, use = "p")
      p <- corPvalueStudent(r, nrow(datExpr))
      cat(sprintf("    %-12s n_genes=%-5d  cor(%s)=%6.3f  p=%.2g\n", d, n_genes, TRAIT_COL, r, p))
    }
  }
  cat("\n")
}

cat(">>> If a cutHeight here gives you a few reasonably-sized, still-TB-correlated\n")
cat("    pieces (not another single 1000+ gene blob, not 40 singleton modules),\n")
cat("    that's a good candidate mergeCutHeight to re-run WGCNA_TBmodules_step1.R\n")
cat("    with (change the mergeCutHeight value inside blockwiseModules() there) --\n")
cat("    then re-run WGCNA_TBmodules_step2_enrichment.R on the new module set.\n")
