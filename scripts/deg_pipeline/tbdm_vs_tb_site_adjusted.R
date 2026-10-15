## TB with diabetes vs TB alone, all samples, adjusted for recruitment site (GSE114192).
##
## Both groups were recruited at all four sites, so site can be modelled as a covariate without being confounded
## with the comparison; this is the main TB-DM comparison agreed for the analysis. The script fits
##   ~ site + group   on all TB-only and TB-DM samples
## and compares the result with the unadjusted run (~ group), reports the significant genes, and runs GSEA on the
## Hallmark interferon gene sets.
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/tbdm_vs_tb_site_adjusted.R
## Outputs in GSE114192/TB_DM/tb_vs_tbdm_site_adjusted/.

## ---- CONFIG ----
PADJ <- 0.05; LFC <- 1; SEED <- 42

for (p in c("DESeq2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(DESeq2))
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, "GSE114192", "TB_DM"); OUT <- file.path(B, "tb_vs_tbdm_site_adjusted"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

## ---- data ----
counts <- rd(B, "03_filtered_counts", "GSE114192_TB_DM_filtered_count_matrix.csv"); rownames(counts) <- counts[, 1]
counts <- round(as.matrix(counts[, -1])); storage.mode(counts) <- "integer"
meta <- rd(B, "02_metadata", "GSE114192_TB_DM_sample_metadata.csv")
geo  <- read.delim(file.path(B, "00_design", "GSE114192_geo_series_matrix_metadata.tsv"), stringsAsFactors = FALSE)
meta$site <- geo$field_site_GEO[match(meta$sample_id, geo$gsm)]
meta <- meta[match(colnames(counts), meta$sample_id), ]
if (!identical(meta$sample_id, colnames(counts)) || anyNA(meta$site)) stop("Counts, metadata and site do not match by sample name.")
keep <- meta$group %in% c("TB_Only", "TB_DM")
m <- droplevels(meta[keep, ]); m$group <- factor(m$group, levels = c("TB_Only", "TB_DM")); m$site <- factor(m$site)
cat(">>> samples by group and site:\n"); print(table(m$group, m$site))
write.csv(as.data.frame.matrix(table(m$group, m$site)), file.path(OUT, "GSE114192_TB_vs_TBDM_site_by_group.csv"))

## ---- model ----
dds <- DESeq(DESeqDataSetFromMatrix(counts[, keep], m, ~ site + group), quiet = TRUE)
coef <- "group_TB_DM_vs_TB_Only"
res <- results(dds, name = coef)
shrink <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
res_s <- lfcShrink(dds, coef = coef, type = shrink, res = res, quiet = TRUE)
out <- data.frame(gene = rownames(res_s), baseMean = res_s$baseMean, log2FoldChange = res_s$log2FoldChange,
                  lfcSE = res_s$lfcSE, pvalue = res_s$pvalue, padj = res_s$padj, LFC_type = paste0(shrink, " shrunk"))
out <- out[order(out$padj), ]
write.csv(out, file.path(OUT, "GSE114192_TB_Only_vs_TB_DM_site_adjusted_DEG.csv"), row.names = FALSE)
sig <- out[!is.na(out$padj) & out$padj < PADJ & abs(out$log2FoldChange) >= LFC, ]
write.csv(sig, file.path(OUT, "GSE114192_TB_Only_vs_TB_DM_site_adjusted_sigDEGs.csv"), row.names = FALSE)
cat("\n>>> TB-DM vs TB-only, ~ site + group:", nrow(sig), "significant genes (padj <", PADJ, ", |log2FC| >=", LFC, ");",
    sum(out$padj < PADJ, na.rm = TRUE), "with padj <", PADJ, "at any fold change\n")

## ---- compare with the unadjusted run ----
ref <- rd(B, "06_deg_results", "GSE114192_TB_Only_vs_TB_DM_DEG.csv"); ref <- ref[!is.na(ref$padj), ]; rownames(ref) <- ref$gene
ref_sig <- ref$gene[ref$padj < PADJ & abs(ref$log2FoldChange) >= LFC]
x <- out[!is.na(out$padj), ]; rownames(x) <- x$gene; g <- intersect(rownames(x), rownames(ref)); rs <- intersect(ref_sig, g)
cmp <- data.frame(Unadjusted_significant = length(ref_sig), Adjusted_significant = nrow(sig),
  Kept = length(intersect(ref_sig, sig$gene)), Same_direction_pct = round(100 * mean(sign(x[rs, "log2FoldChange"]) == sign(ref[rs, "log2FoldChange"])), 1),
  LFC_spearman = round(stats::cor(x[g, "log2FoldChange"], ref[g, "log2FoldChange"], method = "spearman"), 3))
write.csv(cmp, file.path(OUT, "GSE114192_TB_vs_TBDM_adjusted_vs_unadjusted.csv"), row.names = FALSE)
cat("\n========== site-adjusted vs unadjusted ==========\n"); print(cmp, row.names = FALSE)
if (nrow(sig)) {
  sym <- if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, sig$gene, "SYMBOL", "ENSEMBL")) else sig$gene
  cat("\n    significant genes:", paste(ifelse(is.na(sym), sig$gene, sym), sprintf("(%+.1f)", sig$log2FoldChange), collapse = ", "), "\n")
}

## ---- interferon pathways (same ranking as the other GSEA runs) ----
if (all(vapply(c("clusterProfiler", "msigdbr", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(out$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  run_gsea <- function(d) {
    d <- d[!is.na(d$pvalue) & !is.na(d$log2FoldChange), ]; d$stat <- sign(d$log2FoldChange) * -log10(pmax(d$pvalue, 1e-300))
    d <- merge(d[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g,
                                                        pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0))
  }
  ga <- run_gsea(out); gu <- run_gsea(ref)
  write.csv(ga[, c("ID", "setSize", "NES", "pvalue", "p.adjust")], file.path(OUT, "GSE114192_TB_vs_TBDM_site_adjusted_GSEA_HALLMARK.csv"), row.names = FALSE)
  ifn <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
  tab <- merge(gu[gu$ID %in% ifn, c("ID", "NES", "p.adjust")], ga[ga$ID %in% ifn, c("ID", "NES", "p.adjust")], by = "ID", suffixes = c("_unadjusted", "_site_adjusted"))
  cat("\n========== interferon Hallmark sets, TB-DM vs TB-only ==========\n"); print(tab, row.names = FALSE)
  top <- ga[ga$p.adjust < 0.05, ]; top <- top[order(top$NES), c("ID", "NES", "p.adjust")]
  cat("\n    all Hallmark sets with padj < 0.05 (site-adjusted):\n"); print(top, row.names = FALSE)
} else cat("\n(clusterProfiler/msigdbr/org.Hs.eg.db not installed: GSEA skipped)\n")
cat("\n>>> DONE. Outputs in", OUT, "\n")
