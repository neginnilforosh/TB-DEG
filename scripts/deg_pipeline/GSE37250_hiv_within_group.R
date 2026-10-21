## GSE37250: where does the weaker interferon response to TB in people with HIV come from?
##
## The TB x HIV interaction says that the difference between active and latent TB is smaller in people with HIV.
## That can happen in two ways, which the interaction alone cannot separate:
##   (a) interferon genes are LOWER in active TB with HIV than in active TB without HIV (a lower peak), or
##   (b) interferon genes are already HIGHER in latent TB with HIV than without (HIV raises the baseline),
## so the two direct comparisons below are run with the same model as the main analysis
## (limma, ~ country + group, BeadChip as a blocking factor):
##   ATB_HIVneg -> ATB_HIVpos     and     LTBI_HIVneg -> LTBI_HIVpos
## Positive NES = higher in the HIV-positive group. Run after GSE37250_run.R:
##   Rscript scripts/deg_pipeline/GSE37250_hiv_within_group.R
## Outputs in GSE37250/06_deg_results/ (files named *_HIV_within_*.csv, not *_DEG.csv, so they stay out of the
## network pipeline).

ACC <- "GSE37250"; PADJ <- 0.05; LFC <- 1; SEED <- 42
for (p in c("limma")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, ACC)
x <- as.matrix(read.csv(file.path(B, "04_normalized", paste0(ACC, "_vst_normalized_matrix.csv")), row.names = 1, check.names = FALSE))
S <- read.csv(file.path(B, "02_metadata", paste0(ACC, "_sample_metadata.csv")), stringsAsFactors = FALSE)
S <- S[match(colnames(x), S$sample_id), ]; stopifnot(identical(S$sample_id, colnames(x)))

run <- function(ctrl, case) {
  k <- S$group %in% c(ctrl, case); s <- S[k, ]; s$g <- factor(s$group, levels = c(ctrl, case))
  d <- model.matrix(~ country + g, s); y <- x[, k]
  dc <- limma::duplicateCorrelation(y, d, block = s$chip)
  fit <- limma::eBayes(limma::lmFit(y, d, block = s$chip, correlation = dc$consensus.correlation))
  t <- limma::topTable(fit, coef = ncol(d), number = Inf, sort.by = "none")
  out <- data.frame(gene = rownames(t), baseMean = t$AveExpr, log2FoldChange = t$logFC, pvalue = t$P.Value, padj = t$adj.P.Val)
  out <- out[order(out$pvalue), ]
  write.csv(out, file.path(B, "06_deg_results", paste0(ACC, "_HIV_within_", ctrl, "_vs_", case, ".csv")), row.names = FALSE)
  cat(sprintf("  %-26s n = %3d vs %3d | padj < %.2f: %5d | |log2FC| >= %g: %4d (up %d, down %d)\n", paste(ctrl, "->", case),
              sum(s$g == ctrl), sum(s$g == case), PADJ, sum(out$padj < PADJ), LFC, sum(out$padj < PADJ & abs(out$log2FoldChange) >= LFC),
              sum(out$padj < PADJ & out$log2FoldChange >= LFC), sum(out$padj < PADJ & out$log2FoldChange <= -LFC)))
  out
}
cat("========== HIV positive vs negative, within the same TB state ==========\n")
res <- list(ATB = run("ATB_HIVneg", "ATB_HIVpos"), LTBI = run("LTBI_HIVneg", "LTBI_HIVpos"))

if (all(vapply(c("clusterProfiler", "msigdbr", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(rownames(x), "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  gsea <- function(d) { d$stat <- sign(d$log2FoldChange) * -log10(pmax(d$pvalue, 1e-300))
    d <- merge(d[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g, pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0)) }
  ifn <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
  tab <- do.call(rbind, lapply(names(res), function(nm) { g <- gsea(res[[nm]])
    write.csv(g[, c("ID", "setSize", "NES", "pvalue", "p.adjust")], file.path(B, "06_deg_results", paste0(ACC, "_HIV_within_", nm, "_GSEA_HALLMARK.csv")), row.names = FALSE)
    cbind(Within = nm, g[g$ID %in% ifn, c("ID", "NES", "p.adjust")]) }))
  cat("\n========== interferon Hallmark sets, HIV positive vs negative ==========\n"); print(tab, row.names = FALSE)
  cat("\nReading: a NEGATIVE NES within active TB means a lower interferon peak with HIV (a);\n",
      "a POSITIVE NES within latent TB means HIV raises the baseline (b). Both can be true at once.\n", sep = "")
} else cat("(clusterProfiler/msigdbr not installed: GSEA skipped)\n")
cat(">>> DONE\n")
