## STEP 6: TB_UP / TB_DOWN signature for iLINCS (drugfindR prepareSignature + filterSignature).

install.packages("drugfindR", repos = c("https://cogdisreslab.r-universe.dev", "https://cran.r-project.org"))

required_pkgs <- c("drugfindR", "dplyr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "),
       '), repos = c("https://cogdisreslab.r-universe.dev", "https://cran.r-project.org"))')
}

## ---- CONFIG ----
ACC             <- "GSE114192"
LOGFC_THRESHOLD <- 0.5   # minimum |log2FC| for TB_UP / TB_DOWN (applied after prepareSignature)

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
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found
dir_deg    <- file.path(BASE_DIR, "06_deg_results")
dir_out    <- file.path(BASE_DIR, "13_ilincs_signature")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## ---- load all DEGs with padj < 0.05 (the |log2FC| filter is applied below) ----
deg_file <- file.path(dir_deg, paste0(ACC, "_HealthyControl_vs_TBOnly_DEG.csv"))
if (!file.exists(deg_file)) stop("Missing: ", deg_file)
deg <- read.csv(deg_file, stringsAsFactors = FALSE)
deg <- subset(deg, padj < 0.05)
cat(">>> Loaded", nrow(deg), "significant DEGs (before logFC filter).\n")

## ---- L1000/iLINCS needs gene SYMBOLS, not Ensembl IDs ----
have_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)
if (!have_orgdb) stop("Need clusterProfiler + org.Hs.eg.db for ENSEMBL->SYMBOL mapping (already used in step 2/4/5).")
id_map <- suppressMessages(clusterProfiler::bitr(deg$gene, fromType = "ENSEMBL", toType = "SYMBOL",
                                                  OrgDb = org.Hs.eg.db::org.Hs.eg.db))
deg_sym <- merge(deg, id_map, by.x = "gene", by.y = "ENSEMBL")
deg_sym <- deg_sym[!duplicated(deg_sym$SYMBOL), ]
cat(">>> Mapped", nrow(deg_sym), "/", nrow(deg), "to gene symbols.\n")

## ---- prepareSignature(): standardize + map to L1000 gene space ----
signature <- drugfindR::prepareSignature(deg_sym, geneColumn = "SYMBOL",
                                          logfcColumn = "log2FoldChange", pvalColumn = "padj")
cat(">>> prepareSignature() kept", nrow(signature), "/", nrow(deg_sym), "genes (some may fall outside L1000 space).\n")

## ---- filterSignature(): split into TB_UP / TB_DOWN ----
TB_UP   <- drugfindR::filterSignature(signature, direction = "up",   threshold = LOGFC_THRESHOLD)
TB_DOWN <- drugfindR::filterSignature(signature, direction = "down", threshold = LOGFC_THRESHOLD)
cat(">>> TB_UP:  ", nrow(TB_UP),   "genes\n")
cat(">>> TB_DOWN:", nrow(TB_DOWN), "genes\n")

# getConcordants() errors on missing values in the signature, so check here
for (nm in c("TB_UP", "TB_DOWN")) {
  d <- get(nm)
  if (anyNA(d)) cat("!! WARNING:", sum(!complete.cases(d)), "row(s) in", nm,
                     "have missing values -- getConcordants() will reject these later.\n")
}
if (nrow(TB_UP) < 10 || nrow(TB_DOWN) < 10) {
  cat("!! NOTE: fewer than 10 genes in one direction -- iLINCS concordance results",
      "tend to get noisy below this size; consider prop= instead of a fixed threshold if this is too thin.\n")
}

write.csv(TB_UP,   file.path(dir_out, paste0(ACC, "_TB_UP_signature.csv")),   row.names = FALSE)
write.csv(TB_DOWN, file.path(dir_out, paste0(ACC, "_TB_DOWN_signature.csv")), row.names = FALSE)
saveRDS(signature, file.path(dir_out, paste0(ACC, "_full_signature.rds")))

cat("\n>>> DONE. Saved to", dir_out, ":\n")
cat("  -", paste0(ACC, "_TB_UP_signature.csv"), "/", paste0(ACC, "_TB_DOWN_signature.csv"), "-- ready for iLINCS\n")
cat("  -", paste0(ACC, "_full_signature.rds"), "-- the full prepared signature (both directions, pre-split)\n")

