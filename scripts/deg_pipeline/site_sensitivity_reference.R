## Site sensitivity analysis for the reference comparison (GSE114192, healthy vs TB).
##
## Why: healthy controls were recruited only in South Africa and Romania, TB patients in all four countries
## (also Indonesia and Peru). Part of the TB-vs-healthy difference could therefore be a country difference.
## This script redoes the differential expression with country handled in three ways and checks whether the
## main results survive:
##   A  reference          all 82 samples, ~ group                     (the published analysis, read from disk)
##   B  SA+RO only         57 samples, ~ group                         (drop the countries without controls)
##   C  SA+RO, adjusted    57 samples, ~ site + group
##   D  all, adjusted      82 samples, ~ site + group                  (group effect informed by SA+RO only)
## plus a direct test of the concern: within TB patients, how different are Indonesia/Peru from SA/RO?
## Reported for each design: log2FC agreement with A, significant genes kept, direction, the High-priority
## genes, and (if clusterProfiler + msigdbr are installed) the interferon Hallmark NES.
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/site_sensitivity_reference.R
## Outputs in GSE114192/sensitivity_site/. Nothing else is changed.

## ---- CONFIG ----
PADJ <- 0.05; LFC <- 1; SEED <- 42
SAME_SITES <- c("South_Africa", "Romania")

for (p in c("DESeq2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(DESeq2))
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, "GSE114192"); OUT <- file.path(B, "sensitivity_site"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

## ---- data ----
counts <- as.matrix(rd(B, "03_filtered_counts", "GSE114192_filtered_count_matrix.csv")); rownames(counts) <- counts[, 1]
counts <- counts[, -1, drop = FALSE]; storage.mode(counts) <- "numeric"; counts <- round(counts)
meta <- rd(B, "02_metadata", "GSE114192_sample_metadata.csv")
geo  <- read.delim(file.path(B, "TB_DM", "00_design", "GSE114192_geo_series_matrix_metadata.tsv"), stringsAsFactors = FALSE)
meta$site <- geo$field_site_GEO[match(meta$sample_id, geo$gsm)]
if (anyNA(meta$site)) stop("No site for: ", paste(meta$sample_id[is.na(meta$site)], collapse = ", "))
meta <- meta[match(colnames(counts), meta$sample_id), ]
if (!identical(meta$sample_id, colnames(counts))) stop("Counts and metadata do not match by sample name.")
meta$group <- factor(meta$group, levels = c("Healthy_Control", "TB_Only")); meta$site <- factor(meta$site)
cat(">>> samples by group and site:\n"); print(table(meta$group, meta$site))

ref <- rd(B, "06_deg_results", "GSE114192_HealthyControl_vs_TBOnly_DEG.csv")
ref <- ref[!is.na(ref$padj), ]; rownames(ref) <- ref$gene
ref_sig <- ref$gene[ref$padj < PADJ & abs(ref$log2FoldChange) >= LFC]
pri <- rd(B, "19_consensus", "GSE114192_GenePriority_all_modules.csv")
high <- pri$Ensembl[pri$Final_Priority == "High"]
cat(">>> reference: ", length(ref_sig), " significant genes, ", length(high), " High-priority genes\n", sep = "")

run_de <- function(keep, design, coef = "group_TB_Only_vs_Healthy_Control") {
  m <- droplevels(meta[keep, ]); dds <- DESeqDataSetFromMatrix(counts[, keep], m, design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = coef)
  shrink <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
  res <- lfcShrink(dds, coef = coef, type = shrink, res = res, quiet = TRUE)
  out <- data.frame(gene = rownames(res), log2FoldChange = res$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
  attr(out, "shrink") <- shrink; out
}
compare <- function(x, label, n) {
  shrink <- if (is.null(attr(x, "shrink"))) "apeglm (reference)" else attr(x, "shrink")
  x <- x[!is.na(x$padj), ]; rownames(x) <- x$gene; g <- intersect(rownames(x), rownames(ref))
  sig <- x$gene[x$padj < PADJ & abs(x$log2FoldChange) >= LFC]
  kept <- intersect(ref_sig, sig)
  same_dir_ref <- mean(sign(x[intersect(ref_sig, g), "log2FoldChange"]) == sign(ref[intersect(ref_sig, g), "log2FoldChange"]))
  h <- intersect(high, g)
  data.frame(Design = label, Samples = n,
    LFC_spearman_vs_reference = round(stats::cor(x[g, "log2FoldChange"], ref[g, "log2FoldChange"], method = "spearman"), 3),
    Significant = length(sig), Reference_sig_kept = length(kept), Pct_reference_sig_kept = round(100 * length(kept) / length(ref_sig), 1),
    Pct_reference_sig_same_direction = round(100 * same_dir_ref, 1),
    High_genes_padj_0.05 = sum(x[h, "padj"] < PADJ), High_genes_same_direction = sum(sign(x[h, "log2FoldChange"]) == sign(ref[h, "log2FoldChange"])),
    High_genes_tested = length(h), LFC_shrinkage = shrink)
}

sa_ro <- meta$site %in% SAME_SITES
cat(">>> running B (SA+RO, ~ group) ...\n");        B_res <- run_de(sa_ro, ~ group)
cat(">>> running C (SA+RO, ~ site + group) ...\n"); C_res <- run_de(sa_ro, ~ site + group)
cat(">>> running D (all, ~ site + group) ...\n");   D_res <- run_de(rep(TRUE, nrow(meta)), ~ site + group)
for (nm in c("B", "C", "D")) write.csv(get(paste0(nm, "_res")), file.path(OUT, paste0("GSE114192_site_design_", nm, "_DEG.csv")), row.names = FALSE)
summ <- rbind(compare(ref, "A reference: all, ~ group", nrow(meta)),
              compare(B_res, "B SA+RO only, ~ group", sum(sa_ro)),
              compare(C_res, "C SA+RO only, ~ site + group", sum(sa_ro)),
              compare(D_res, "D all, ~ site + group", nrow(meta)))
write.csv(summ, file.path(OUT, "GSE114192_site_sensitivity_summary.csv"), row.names = FALSE)
cat("\n========== does the reference result survive country handling? ==========\n"); print(summ, row.names = FALSE)

## ---- the concern itself: within TB, Indonesia/Peru vs South Africa/Romania ----
tb <- meta$group == "TB_Only"; meta$region <- factor(ifelse(meta$site %in% SAME_SITES, "SA_RO", "ID_PE"), levels = c("SA_RO", "ID_PE"))
cat("\n>>> within TB patients: Indonesia/Peru (", sum(tb & meta$region == "ID_PE"), ") vs SA/RO (", sum(tb & meta$region == "SA_RO"), ")\n", sep = "")
m <- droplevels(meta[tb, ]); dds <- DESeq(DESeqDataSetFromMatrix(counts[, tb], m, ~ region), quiet = TRUE)
r <- as.data.frame(results(dds, name = "region_ID_PE_vs_SA_RO")); r$gene <- rownames(r)
write.csv(r, file.path(OUT, "GSE114192_TB_only_region_DEG.csv"), row.names = FALSE)
reg_sig <- r$gene[!is.na(r$padj) & r$padj < PADJ & abs(r$log2FoldChange) >= LFC]
cat("    genes differing by region within TB (padj < 0.05, |log2FC| >= 1):", length(reg_sig), "\n")
cat("    of the", length(ref_sig), "reference DEGs:", length(intersect(reg_sig, ref_sig)), "| of the", length(high), "High genes:",
    length(intersect(reg_sig, high)), "\n")

## ---- interferon pathways (optional) ----
if (all(vapply(c("clusterProfiler", "msigdbr", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(ref$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  nes <- function(x) {
    x <- x[!is.na(x$pvalue), ]; x$stat <- sign(x$log2FoldChange) * -log10(pmax(x$pvalue, 1e-300))
    d <- merge(x[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); g <- clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g,
                                                pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0)
    g <- as.data.frame(g); g[g$ID %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"), c("ID", "NES", "p.adjust")]
  }
  ifn <- do.call(rbind, lapply(list(A = ref, B = B_res, C = C_res, D = D_res), nes))
  ifn$Design <- rep(c("A", "B", "C", "D"), each = 2)
  write.csv(ifn, file.path(OUT, "GSE114192_site_sensitivity_interferon_GSEA.csv"), row.names = FALSE)
  cat("\n========== interferon Hallmark NES by design ==========\n"); print(ifn[, c("Design", "ID", "NES", "p.adjust")], row.names = FALSE)
} else cat("\n(clusterProfiler/msigdbr/org.Hs.eg.db not installed: interferon GSEA skipped)\n")
cat("\n>>> DONE. Outputs in", OUT, "\n")
