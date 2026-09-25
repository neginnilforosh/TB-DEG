## =============================================================
## STEP 4 of Ratul's roadmap: connect DEG results with WGCNA modules
##
## "For every important module, check which genes are also
##  significant TB-vs-Control DEGs. Save a table containing Gene,
##  log2FC, FDR, WGCNA module, MM/kME, and Gene Significance. The
##  genes that are both strong DEGs AND strong members of TB-related
##  modules will become the main disease genes for the next analyses."
##
## Pure base-R table join -- no Bioconductor packages needed.
## =============================================================

## ---- CONFIG ----
ACC <- "GSE114192"

## ---- locate script dir (same pattern as step 1/2) ----
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
BASE_DIR   <- file.path(SCRIPT_DIR, "..", ACC)
dir_deg    <- file.path(BASE_DIR, "06_deg_results")
dir_wgcna  <- file.path(BASE_DIR, "10_wgcna_results")

## ---- load DEG results + WGCNA module membership/gene significance ----
deg_file <- file.path(dir_deg, paste0(ACC, "_HealthyControl_vs_TBOnly_DEG.csv"))
mm_file  <- file.path(dir_wgcna, paste0(ACC, "_TBmodules_MM_GS.csv"))
if (!file.exists(deg_file)) stop("Missing: ", deg_file)
if (!file.exists(mm_file))  stop("Missing: ", mm_file, " -- run WGCNA_TBmodules_step1.R first.")

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
mm  <- read.csv(mm_file,  stringsAsFactors = FALSE)
cat("DEG table:", nrow(deg), "genes.  WGCNA TB-module table:", nrow(mm), "genes across",
    length(unique(mm$Module)), "modules (", paste(unique(mm$Module), collapse=", "), ")\n")

## ---- join on gene ID (both use ENSEMBL IDs) ----
merged <- merge(mm, deg[, c("gene", "log2FoldChange", "padj", "significant")],
                by.x = "Gene", by.y = "gene", all.x = TRUE)
if (any(is.na(merged$significant))) {
  cat("!! WARNING:", sum(is.na(merged$significant)), "module genes had no matching DEG row -- check gene ID formats.\n")
}

## ---- "main disease genes" = significant DEG AND TB-module member ----
disease_genes <- merged[!is.na(merged$significant) & merged$significant == TRUE, ]
disease_genes <- disease_genes[order(disease_genes$Module, disease_genes$padj), ]

out <- data.frame(
  Gene   = disease_genes$Gene,
  log2FC = disease_genes$log2FoldChange,
  FDR    = disease_genes$padj,
  Module = disease_genes$Module,
  MM     = disease_genes$MM,
  GS_TB  = disease_genes$GS_TB
)

out_file <- file.path(dir_wgcna, paste0(ACC, "_DiseaseGenes_DEG_WGCNA.csv"))
write.csv(out, out_file, row.names = FALSE)

cat("\n>>> Disease genes (significant DEG AND TB-module member):", nrow(out), "total\n")
print(table(out$Module))
cat("\n>>> Saved:", out_file, "\n")
cat(">>> Next (Ratul step 5): build one STRING PPI network PER MODULE from these genes\n")
