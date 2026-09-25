## =============================================================
## What this does, per the brief:
##   1. Build the WGCNA network (same as WGCNA_GSE114192_run.R)
##   2. Module-trait correlation heatmap (TB vs Control)
##   3. Select ONLY the modules significantly associated with TB
##   4. For each selected module, export: module color, N genes,
##      correlation with TB, p-value + FDR, module membership
##      (MM/kME), and gene significance (GS)
##
## Reusable for any dataset: just edit ACC and TRAIT_COL below.
## =============================================================

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
# install.packages("WGCNA")

library(WGCNA)
library(dplyr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

## ---- 1. CONFIG — edit these two lines per dataset -----------

ACC       <- "GSE114192"   # folder name, matches the rest of the repo
TRAIT_COL <- "TB_Only"     # which trait column = "has TB" (1) vs control (0)
FDR_CUTOFF      <- 0.05
MAX_SIG_MODULES <- 3   # cap: at most this many modules go on to step 3 (naming/enrichment) —
                        # picked as the strongest |correlation| among the FDR-significant ones

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }

  getwd()
}

SCRIPT_DIR <- get_script_dir()
cat("Script folder detected as:", SCRIPT_DIR, "\n")
BASE_DIR    <- file.path(SCRIPT_DIR, "..", ACC)
dir_wgcna   <- file.path(BASE_DIR, "09_wgcna_input")
dir_results <- file.path(BASE_DIR, "10_wgcna_results")
dir_fig     <- file.path(BASE_DIR, "figures")
dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)

## ---- 2. Load data (same convention as WGCNA_GSE114192_run.R) ----
cat("\n========== STEP 1: LOAD AND FORMAT DATA ==========\n")

expr_file <- file.path(dir_wgcna, paste0(ACC, "_WGCNA_expression_matrix.csv"))
if (!file.exists(expr_file)) stop("Expression matrix not found: ", expr_file)
datExpr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
if (nrow(datExpr) > ncol(datExpr)) {
  cat("Transposing expression matrix (samples as rows)...\n")
  datExpr <- as.data.frame(t(datExpr))
}
cat("Expression data:", nrow(datExpr), "samples,", ncol(datExpr), "genes.\n")

trait_file <- file.path(dir_wgcna, paste0(ACC, "_WGCNA_trait_file.csv"))
if (!file.exists(trait_file)) stop("Trait file not found: ", trait_file)
datTraits <- read.csv(trait_file, row.names = 1, check.names = FALSE)

common_samples <- intersect(rownames(datExpr), rownames(datTraits))
if (length(common_samples) == 0) stop("No matching sample IDs between expression and trait data!")
datExpr   <- datExpr[common_samples, ]
datTraits <- datTraits[common_samples, , drop = FALSE]
cat("Matched", length(common_samples), "samples.\n")

if (!(TRAIT_COL %in% colnames(datTraits))) {
  stop("TRAIT_COL '", TRAIT_COL, "' not in trait file. Available columns: ",
       paste(colnames(datTraits), collapse = ", "))
}

## ---- 3. Soft thresholding + network construction ----
cat("\n========== STEP 2: SOFT THRESHOLDING ==========\n")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 3)

pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_SoftThresholding.pdf")), width = 9, height = 5)
par(mfrow = c(1, 2)); cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower)) { cat("WARNING: no power reached R^2 cutoff; defaulting to 6.\n"); softPower <- 6 }
cat(">>> Selected soft power:", softPower, "\n")

cat("\n========== STEP 3: NETWORK CONSTRUCTION ==========\n")
net <- blockwiseModules(datExpr, power = softPower,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE, verbose = 3)

moduleColors <- labels2colors(net$colors)
cat(">>> Found", length(unique(moduleColors)), "modules (including grey/unassigned).\n")

pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_ModuleDendrogram.pdf")), width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                     "Module Colors", dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     main = paste("Gene Dendrogram and Module Colors (", ACC, ")", sep = ""))
dev.off()

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

# save the network so later steps (DEG integration, STRING export, etc.)
# don't have to rebuild it from scratch
saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs,
             datExpr = datExpr, datTraits = datTraits),
        file.path(dir_results, paste0(ACC, "_WGCNA_workspace.rds")))

## ---- 4. Module-trait correlation heatmap (all traits) ----
cat("\n========== STEP 4: MODULE-TRAIT CORRELATION ==========\n")
nSamples <- nrow(datExpr)
moduleTraitCor    <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_ModuleTrait_Heatmap.pdf")), width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),
               yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1, 1),
               main = paste("Module-Trait Relationships (", ACC, ")", sep = ""))
dev.off()

## ---- 5. Select modules significantly associated with TRAIT_COL ----
cat("\n========== STEP 5: SELECT TB-SIGNIFICANT MODULES ==========\n")
tb_cor  <- moduleTraitCor[, TRAIT_COL]
tb_pval <- moduleTraitPvalue[, TRAIT_COL]
tb_fdr  <- p.adjust(tb_pval, method = "BH")

module_names <- sub("^ME", "", names(tb_cor))
n_genes_tab  <- table(moduleColors)

module_summary <- data.frame(
  Module      = module_names,
  N_genes     = as.integer(n_genes_tab[module_names]),
  Correlation = as.numeric(tb_cor),
  P_value     = as.numeric(tb_pval),
  FDR         = as.numeric(tb_fdr)
)
module_summary <- module_summary[module_summary$Module != "grey", ]  # grey = unassigned genes
module_summary <- module_summary[order(module_summary$FDR), ]

SIG_MODULES_ALL <- module_summary$Module[module_summary$FDR < FDR_CUTOFF]
cat(">>> Modules passing FDR <", FDR_CUTOFF, ":",
    if (length(SIG_MODULES_ALL)) paste(SIG_MODULES_ALL, collapse = ", ") else "NONE",
    "(", length(SIG_MODULES_ALL), "total )\n")

# cap at MAX_SIG_MODULES, keeping the strongest EFFECT SIZE (|correlation|), not just lowest p/FDR —
# with large N, weak-but-"significant" modules (e.g. r~0.3) shouldn't outrank strong ones
ranked_sig <- module_summary[module_summary$Module %in% SIG_MODULES_ALL, ]
ranked_sig <- ranked_sig[order(-abs(ranked_sig$Correlation)), ]
SIG_MODULES <- head(ranked_sig$Module, MAX_SIG_MODULES)
cat(">>> Keeping top", MAX_SIG_MODULES, "by |correlation| ->", paste(SIG_MODULES, collapse = ", "), "\n")

module_summary$Selected <- module_summary$Module %in% SIG_MODULES  # TRUE = one of the final top-N

write.csv(module_summary,
          file.path(dir_results, paste0(ACC, "_ModuleTrait_Correlation.csv")),
          row.names = FALSE)

## ---- 6. Module Membership (MM/kME) + Gene Significance (GS) ----
cat("\n========== STEP 6: MODULE MEMBERSHIP + GENE SIGNIFICANCE ==========\n")
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM.", names(MEs))
names(MMPvalue)              <- paste0("p.MM.", names(MEs))

tb_numeric <- as.numeric(datTraits[[TRAIT_COL]])
geneTraitSignificance <- as.data.frame(cor(datExpr, tb_numeric, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.TB"
names(GSPvalue)              <- "p.GS.TB"

gene_table <- data.frame(
  Gene    = colnames(datExpr),
  Module  = moduleColors,
  GS_TB   = geneTraitSignificance$GS.TB,
  p_GS_TB = GSPvalue$p.GS.TB,
  MM      = NA_real_,
  p_MM    = NA_real_
)

# each gene's MM/kME is its correlation with ITS OWN module's eigengene
for (m in unique(moduleColors)) {
  idx    <- which(moduleColors == m)
  col_mm <- paste0("MM.ME", m)
  col_p  <- paste0("p.MM.ME", m)
  if (col_mm %in% names(geneModuleMembership)) {
    gene_table$MM[idx]   <- geneModuleMembership[idx, col_mm]
    gene_table$p_MM[idx] <- MMPvalue[idx, col_p]
  }
}

gene_table_sig <- gene_table[gene_table$Module %in% SIG_MODULES, ]
gene_table_sig <- gene_table_sig[order(gene_table_sig$Module, -abs(gene_table_sig$MM)), ]

write.csv(gene_table_sig,
          file.path(dir_results, paste0(ACC, "_TBmodules_MM_GS.csv")),
          row.names = FALSE)

write.csv(data.frame(Gene = colnames(datExpr), ModuleColor = moduleColors),
          file.path(dir_results, paste0(ACC, "_Gene_Module_Assignment.csv")),
          row.names = FALSE)

cat("\n>>> DONE. New/updated files in", dir_results, ":\n")
cat("  -", paste0(ACC, "_ModuleTrait_Correlation.csv"), " <- 1 row/module: N_genes, correlation, p, FDR\n")
cat("  -", paste0(ACC, "_TBmodules_MM_GS.csv"), "      <- 1 row/gene, TB-significant modules only: MM/kME + GS\n")
cat("  -", paste0(ACC, "_Gene_Module_Assignment.csv"), " <- 1 row/gene, all modules (unchanged format)\n")
cat("  -", paste0(ACC, "_WGCNA_workspace.rds"), "       <- net/MEs/datExpr/datTraits for the next step (DEG integration)\n")

