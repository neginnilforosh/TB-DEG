## STEP 17: integrate the GSE165709 results with the earlier work.
## Reads the per-cohort RNA and ATAC results from step 16 and answers, for each gene:
##   - does it respond to Mtb in HC, in PLWH, or in both?
##   - is the response supported at the expression level, the chromatin level, or both?
##   - was it already important in the previous analysis (module, DEG overlap, priority tier)?
## Outputs in external_results/GSE165709/integration/. Base R only.

## ---- CONFIG ----
COHORTS   <- c("HC", "PLWH")       # first = reference (normal response)
PADJ_CUT  <- 0.05
LFC_CUT   <- 0                     # raise to e.g. 0.5 for a stricter definition of "responds"

## ---- paths ----
get_script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
EXT      <- file.path(ROOT_DIR, "external_results", "GSE165709")
dir_out  <- file.path(EXT, "integration"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) { if (!file.exists(f)) stop("Missing: ", f, "\n  Run GSE165709_differential_step16.R for both assays first."); read.csv(f, stringsAsFactors = FALSE) }

## ---- 1. RNA: one row per gene per cohort ----
rna <- lapply(setNames(COHORTS, COHORTS), function(co)
  rd(file.path(EXT, "RNA", paste0("GSE165709_RNA_INF_vs_NEG_in_", co, ".csv"))))
rna_diff <- rd(file.path(EXT, "RNA", paste0("GSE165709_RNA_response_", COHORTS[2], "_vs_", COHORTS[1], ".csv")))
## use the gene symbol when step 16 produced one (the RNA ids are "ENSG..._SYMBOL"), so that the
## RNA, the ATAC annotation and the earlier gene tables are all keyed the same way
key <- function(x) if ("SYMBOL" %in% names(x)) x$SYMBOL else x$feature
for (co in names(rna)) rna[[co]]$gene_key <- key(rna[[co]])
rna_diff$gene_key <- key(rna_diff)
genes <- sort(unique(unlist(lapply(rna, function(x) x$gene_key))))
cat(">>> RNA: ", length(genes), " genes tested | significant: ",
    paste(COHORTS, vapply(rna, function(x) sum(x$padj < PADJ_CUT, na.rm = TRUE), 1L), collapse = ", "), "\n", sep = "")

G <- data.frame(Gene = genes, stringsAsFactors = FALSE)
for (co in COHORTS) {
  m <- match(G$Gene, rna[[co]]$gene_key)
  G[[paste0("RNA_log2FC_", co)]] <- rna[[co]]$log2FoldChange[m]
  G[[paste0("RNA_padj_", co)]]   <- rna[[co]]$padj[m]
  G[[paste0("RNA_responds_", co)]] <- !is.na(rna[[co]]$padj[m]) & rna[[co]]$padj[m] < PADJ_CUT &
                                       abs(rna[[co]]$log2FoldChange[m]) > LFC_CUT
}
m <- match(G$Gene, rna_diff$gene_key)
G$RNA_response_diff_log2FC <- rna_diff$log2FoldChange[m]
G$RNA_response_diff_padj   <- rna_diff$padj[m]
G$RNA_response_diff_pvalue <- rna_diff$pvalue[m]

## ---- 2. ATAC: a gene counts as responding if any annotated peak responds ----
atac_genes <- function(co) {
  f <- file.path(EXT, "ATAC", paste0("GSE165709_ATAC_INF_vs_NEG_in_", co, "_annotated.csv"))
  a <- rd(f); a <- a[!is.na(a$SYMBOL), ]
  sig <- a[!is.na(a$padj) & a$padj < PADJ_CUT & abs(a$log2FoldChange) > LFC_CUT, ]
  sig <- sig[order(sig$padj), ]; best <- sig[!duplicated(sig$SYMBOL), ]
  list(tested = unique(a$SYMBOL), sig = best)
}
at <- lapply(setNames(COHORTS, COHORTS), atac_genes)
cat(">>> ATAC: genes with a significant peak | ",
    paste(COHORTS, vapply(at, function(x) nrow(x$sig), 1L), collapse = ", "), "\n", sep = "")
for (co in COHORTS) {
  m <- match(G$Gene, at[[co]]$sig$SYMBOL)
  G[[paste0("ATAC_log2FC_", co)]]   <- at[[co]]$sig$log2FoldChange[m]
  G[[paste0("ATAC_padj_", co)]]     <- at[[co]]$sig$padj[m]
  G[[paste0("ATAC_annotation_", co)]] <- at[[co]]$sig$annotation[m]
  G[[paste0("ATAC_responds_", co)]] <- G$Gene %in% at[[co]]$sig$SYMBOL
  G[[paste0("ATAC_tested_", co)]]   <- G$Gene %in% at[[co]]$tested
}

## ---- 3. support across layers and the blunted response ----
for (co in COHORTS) {
  r <- G[[paste0("RNA_responds_", co)]]; a <- G[[paste0("ATAC_responds_", co)]]
  G[[paste0("Support_", co)]] <- ifelse(r & a, "RNA+ATAC", ifelse(r, "RNA only", ifelse(a, "ATAC only", "none")))
}
ref <- COHORTS[1]; alt <- COHORTS[2]
G$Response_pattern <- with(G, ifelse(G[[paste0("RNA_responds_", ref)]] & G[[paste0("RNA_responds_", alt)]], "both cohorts",
                             ifelse(G[[paste0("RNA_responds_", ref)]], paste(ref, "only"),
                             ifelse(G[[paste0("RNA_responds_", alt)]], paste(alt, "only"), "neither"))))
## "Blunted": responds in the reference cohort and not in the other one, with a smaller fold change there.
## The formal interaction test (RNA_response_diff_padj) is reported separately because it is much less
## powerful: it tests a difference of two log fold changes, so its standard error is about 1.4x larger
## than either single contrast, and in this dataset it yields no FDR-significant gene at all. A gene that
## is significant in one cohort and not the other is NOT by itself evidence that the responses differ --
## hence the two columns: Blunted (the descriptive pattern) and Blunted_interaction_support (the test).
smaller <- !is.na(G[[paste0("RNA_log2FC_", alt)]]) &
  abs(G[[paste0("RNA_log2FC_", alt)]]) < abs(G[[paste0("RNA_log2FC_", ref)]])
G$Blunted <- G[[paste0("RNA_responds_", ref)]] & !G[[paste0("RNA_responds_", alt)]] & smaller
G$Blunted_interaction_p    <- G$RNA_response_diff_pvalue
G$Blunted_interaction_padj <- G$RNA_response_diff_padj
G$Blunted_interaction_support <- ifelse(!is.na(G$RNA_response_diff_padj) & G$RNA_response_diff_padj < PADJ_CUT, "FDR",
                                 ifelse(!is.na(G$RNA_response_diff_pvalue) & G$RNA_response_diff_pvalue < 0.05, "nominal only", "none"))

## ---- 4. link to the earlier analysis ----
gp <- rd(file.path(run_out_dir(file.path(ROOT_DIR, ACC)), "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
m <- match(G$Gene, gp$Gene)
G$Prev_Module <- gp$Module[m]; G$Prev_Priority <- gp$Final_Priority[m]
G$Prev_DEG_overlap <- gp$DEG_WGCNA_overlap[m]; G$Prev_log2FC <- gp$log2FC[m]
G$In_previous_modules <- !is.na(G$Prev_Module)
G <- G[order(-(G$Response_pattern == "both cohorts"), G$RNA_padj_HC), ]
write.csv(G, file.path(dir_out, "GSE165709_gene_integration.csv"), row.names = FALSE)

cat("\n>>> response pattern (RNA):\n"); print(table(G$Response_pattern))
cat("\n>>> layer support:\n"); for (co in COHORTS) { cat("  ", co, ": "); print(table(G[[paste0("Support_", co)]])) }
cat("\n>>> genes responding in", ref, "but not in", alt, ", with a smaller fold change:", sum(G$Blunted, na.rm = TRUE), "\n")
cat("    of these, the formal interaction test supports the difference for",
    sum(G$Blunted & G$Blunted_interaction_support == "FDR", na.rm = TRUE), "genes at FDR <", PADJ_CUT, "and",
    sum(G$Blunted & G$Blunted_interaction_support == "nominal only", na.rm = TRUE), "at nominal p < 0.05.\n")
cat("    (A gene significant in one cohort but not the other is a weaker claim than a tested difference;\n")
cat("     use Blunted to generate candidates and Blunted_interaction_support to grade them.)\n")

## ---- 5. do the previously important genes respond more often? ----
test <- function(sel, label, col) {
  u <- !is.na(G$Prev_Module)                       # only genes that the earlier analysis covered
  a <- sum(G[[col]][u & sel], na.rm = TRUE); n <- sum(u & sel)
  b <- sum(G[[col]][u & !sel], na.rm = TRUE); nb <- sum(u & !sel)
  if (n == 0 || nb == 0) return(NULL)
  p <- stats::fisher.test(matrix(c(a, n - a, b, nb - b), 2), alternative = "greater")$p.value
  data.frame(Layer = col, Group = label, N = n, Responding = a, Pct = round(100 * a / n, 1),
             Pct_rest = round(100 * b / nb, 1), P_fisher = signif(p, 3))
}
tests <- list()
for (col in c(paste0("RNA_responds_", COHORTS), paste0("ATAC_responds_", COHORTS))) {
  tests <- c(tests, list(test(G$Prev_Priority == "High", "High priority", col),
                         test(G$Prev_DEG_overlap %in% TRUE, "DEG-overlap genes", col)))
  for (mo in sort(unique(na.omit(G$Prev_Module)))) tests <- c(tests, list(test(G$Prev_Module == mo, paste("module", mo), col)))
}
tt <- do.call(rbind, tests)
if (!is.null(tt)) {
  write.csv(tt, file.path(dir_out, "GSE165709_previous_gene_overlap_tests.csv"), row.names = FALSE)
  cat("\n>>> do the previously important genes respond more often than the rest?\n"); print(tt, row.names = FALSE)
}

## genes supported in every layer: previously important, and responding in RNA and ATAC
core <- G[G$In_previous_modules & G$Support_HC == "RNA+ATAC", ]
write.csv(core, file.path(dir_out, "GSE165709_core_supported_genes.csv"), row.names = FALSE)
cat("\n>>> previously-important genes with both RNA and ATAC support in", ref, ":", nrow(core), "\n")
if (nrow(core)) print(utils::head(core[, c("Gene", "Prev_Module", "Prev_Priority", "RNA_log2FC_HC", "ATAC_log2FC_HC", "Response_pattern")], 10), row.names = FALSE)

cat("\n>>> DONE. Outputs in", dir_out, "\n")
cat(">>> Note: these are alveolar macrophages infected in vitro, while the earlier work is whole blood --\n")
cat("    agreement between them is meaningful, but a gene can respond in one setting and not the other.\n")
