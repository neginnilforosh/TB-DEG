## Site sensitivity analysis for the TB-diabetes analysis (GSE114192).
##
## Healthy and diabetes-only participants come from South Africa and Romania only; TB and TB-DM from all four
## sites. Two results depend on that design and are checked here:
##  1. the main contrast TB-DM vs DM-only (which genes and High-priority genes it finds)
##  2. the TB x DM interaction, whose interferon gene sets were negative ("diabetes blunts the interferon
##     response to TB")
## Each is redone (B) on the South Africa + Romania samples only, (C) the same with site as a covariate, and
## (D) on all samples with site as a covariate, and compared with the published run (A). Within TB-DM patients
## the two regions are also compared directly.
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/site_sensitivity_tbdm.R
## Outputs in GSE114192/TB_DM/sensitivity_site/. Nothing else is changed.

## ---- CONFIG ----
PADJ <- 0.05; LFC <- 1; SEED <- 42
SAME_SITES <- c("South_Africa", "Romania")

for (p in c("DESeq2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(DESeq2))
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, "GSE114192", "TB_DM"); OUT <- file.path(B, "sensitivity_site"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

## ---- data ----
counts <- rd(B, "03_filtered_counts", "GSE114192_TB_DM_filtered_count_matrix.csv"); rownames(counts) <- counts[, 1]
counts <- round(as.matrix(counts[, -1])); storage.mode(counts) <- "integer"
meta <- rd(B, "02_metadata", "GSE114192_TB_DM_sample_metadata.csv")
geo  <- read.delim(file.path(B, "00_design", "GSE114192_geo_series_matrix_metadata.tsv"), stringsAsFactors = FALSE)
meta$site <- geo$field_site_GEO[match(meta$sample_id, geo$gsm)]
meta <- meta[match(colnames(counts), meta$sample_id), ]
if (!identical(meta$sample_id, colnames(counts)) || anyNA(meta$site)) stop("Counts, metadata and site do not match by sample name.")
meta$site <- factor(meta$site); meta$region <- factor(ifelse(meta$site %in% SAME_SITES, "SA_RO", "ID_PE"), levels = c("SA_RO", "ID_PE"))
meta$TB <- factor(ifelse(meta$group %in% c("TB_Only", "TB_DM"), "yes", "no"), levels = c("no", "yes"))
meta$DM <- factor(ifelse(meta$group %in% c("DM_only", "TB_DM"), "yes", "no"), levels = c("no", "yes"))
cat(">>> samples by group and site:\n"); print(table(meta$group, meta$site))

shrink_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
fit <- function(keep, design) { m <- droplevels(meta[keep, ]); DESeq(DESeqDataSetFromMatrix(counts[, keep], m, design), quiet = TRUE) }
tab <- function(dds, coef, shrink = TRUE) {
  r <- results(dds, name = coef)
  if (shrink) r <- lfcShrink(dds, coef = coef, type = shrink_type, res = r, quiet = TRUE)
  data.frame(gene = rownames(r), log2FoldChange = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
}
sa_ro <- meta$site %in% SAME_SITES

## ================= 1. main contrast: TB-DM vs DM-only =================
ref <- rd(B, "06_deg_results", "GSE114192_DM_only_vs_TB_DM_DEG.csv"); ref <- ref[!is.na(ref$padj), ]; rownames(ref) <- ref$gene
ref_sig <- ref$gene[ref$padj < PADJ & abs(ref$log2FoldChange) >= LFC]
pri <- rd(B, "19_consensus", "GSE114192_GenePriority_all_modules.csv"); high <- pri$Ensembl[pri$Final_Priority == "High"]
pair <- meta$group %in% c("DM_only", "TB_DM"); meta$group2 <- factor(ifelse(meta$group == "TB_DM", "TB_DM", "DM_only"), levels = c("DM_only", "TB_DM"))
coef_pair <- "group2_TB_DM_vs_DM_only"
cat("\n>>> main contrast | reference:", length(ref_sig), "significant genes,", length(high), "High genes\n")
designs <- list(B = list(keep = pair & sa_ro, f = ~ group2,        lab = "B SA+RO only, ~ group"),
                C = list(keep = pair & sa_ro, f = ~ site + group2, lab = "C SA+RO only, ~ site + group"),
                D = list(keep = pair,         f = ~ site + group2, lab = "D all, ~ site + group"))
compare <- function(x, label, n) {
  x <- x[!is.na(x$padj), ]; rownames(x) <- x$gene; g <- intersect(rownames(x), rownames(ref)); h <- intersect(high, g)
  sig <- x$gene[x$padj < PADJ & abs(x$log2FoldChange) >= LFC]; rs <- intersect(ref_sig, g)
  data.frame(Design = label, Samples = n,
    LFC_spearman_vs_reference = round(stats::cor(x[g, "log2FoldChange"], ref[g, "log2FoldChange"], method = "spearman"), 3),
    Significant = length(sig), Pct_reference_sig_kept = round(100 * length(intersect(ref_sig, sig)) / length(ref_sig), 1),
    Pct_reference_sig_same_direction = round(100 * mean(sign(x[rs, "log2FoldChange"]) == sign(ref[rs, "log2FoldChange"])), 1),
    High_genes_padj_0.05 = sum(x[h, "padj"] < PADJ, na.rm = TRUE), High_genes_same_direction = sum(sign(x[h, "log2FoldChange"]) == sign(ref[h, "log2FoldChange"])),
    High_genes_tested = length(h))
}
summ <- list(compare(ref, "A reference: all, ~ group", sum(pair)))
for (nm in names(designs)) {
  d <- designs[[nm]]; cat(">>> running", d$lab, "...\n")
  x <- tab(fit(d$keep, d$f), coef_pair); write.csv(x, file.path(OUT, paste0("GSE114192_TBDM_vs_DM_design_", nm, "_DEG.csv")), row.names = FALSE)
  summ[[length(summ) + 1]] <- compare(x, d$lab, sum(d$keep))
}
summ <- do.call(rbind, summ); summ$LFC_shrinkage <- c("as published", rep(shrink_type, 3))
write.csv(summ, file.path(OUT, "GSE114192_TBDM_site_sensitivity_summary.csv"), row.names = FALSE)
cat("\n========== main contrast: does it survive site handling? ==========\n"); print(summ, row.names = FALSE)

## within TB-DM: Indonesia/Peru vs South Africa/Romania
tbdm <- meta$group == "TB_DM"
r <- tab(fit(tbdm, ~ region), "region_ID_PE_vs_SA_RO", shrink = FALSE)
write.csv(r, file.path(OUT, "GSE114192_TBDM_only_region_DEG.csv"), row.names = FALSE)
reg_sig <- r$gene[!is.na(r$padj) & r$padj < PADJ & abs(r$log2FoldChange) >= LFC]
cat("\n>>> within TB-DM patients, ID/PE (", sum(tbdm & meta$region == "ID_PE"), ") vs SA/RO (", sum(tbdm & meta$region == "SA_RO"),
    "): ", length(reg_sig), " genes differ | of the reference DEGs: ", length(intersect(reg_sig, ref_sig)),
    " | of the High genes: ", length(intersect(reg_sig, high)), "\n", sep = "")

## ================= 2. interaction model =================
four <- meta$group %in% c("Healthy_Control", "TB_Only", "DM_only", "TB_DM")
idesigns <- list(I2 = list(keep = four & sa_ro, f = ~ TB + DM + TB:DM,        lab = "I2 SA+RO only"),
                 I3 = list(keep = four & sa_ro, f = ~ site + TB + DM + TB:DM, lab = "I3 SA+RO only, + site"),
                 I4 = list(keep = four,         f = ~ site + TB + DM + TB:DM, lab = "I4 all, + site"))
ires <- list(I1 = rd(B, "06_deg_results", "GSE114192_interaction_interaction.csv"))
for (nm in names(idesigns)) {
  d <- idesigns[[nm]]; cat(">>> interaction:", d$lab, "...\n")
  ires[[nm]] <- tab(fit(d$keep, d$f), "TByes.DMyes", shrink = FALSE)
  write.csv(ires[[nm]], file.path(OUT, paste0("GSE114192_interaction_", nm, ".csv")), row.names = FALSE)
}
labs <- c(I1 = "I1 reference: all, no site", vapply(idesigns, `[[`, "", "lab"))
isumm <- data.frame(Design = labs[names(ires)], Interaction_genes_padj_0.05 = vapply(ires, function(x) sum(x$padj < PADJ, na.rm = TRUE), 0L))

## interferon gene sets on the interaction coefficient (same ranking as the published interaction GSEA)
if (all(vapply(c("clusterProfiler", "msigdbr", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(ires$I1$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  ifn <- function(x) {
    x <- x[!is.na(x$pvalue) & !is.na(x$log2FoldChange), ]; x$stat <- sign(x$log2FoldChange) * -log10(pmax(x$pvalue, 1e-300))
    d <- merge(x[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); g <- as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g,
                                                             pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0))
    s <- function(id, col) { v <- g[g$ID == id, col]; if (length(v)) v else NA }
    c(IFNG_NES = s("HALLMARK_INTERFERON_GAMMA_RESPONSE", "NES"), IFNG_padj = s("HALLMARK_INTERFERON_GAMMA_RESPONSE", "p.adjust"),
      IFNA_NES = s("HALLMARK_INTERFERON_ALPHA_RESPONSE", "NES"), IFNA_padj = s("HALLMARK_INTERFERON_ALPHA_RESPONSE", "p.adjust"))
  }
  isumm <- cbind(isumm, do.call(rbind, lapply(ires, ifn)))
} else cat("(clusterProfiler/msigdbr/org.Hs.eg.db not installed: interferon GSEA on the interaction skipped)\n")
write.csv(isumm, file.path(OUT, "GSE114192_interaction_site_sensitivity.csv"), row.names = FALSE)
cat("\n========== TB x DM interaction: does the interferon blunting survive site handling? ==========\n"); print(isumm, row.names = FALSE)
cat("\nA negative interferon NES that stays significant in I2-I4 supports the blunting finding; if it disappears once\n",
    "site is handled, the finding was driven by where the TB-DM patients were recruited.\n", sep = "")
cat(">>> DONE. Outputs in", OUT, "\n")
