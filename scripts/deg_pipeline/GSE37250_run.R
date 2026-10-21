## GSE37250 (whole blood, Illumina HumanHT-12 v4): expression stage of the TB-HIV context.
##
## Run after GSE37250_design_check.R and GSE37250_chip_mapping.R, which record the sample groups and the BeadChip of
## every array. Steps:
##   1. raw signal from the non-normalised file, columns renamed to GEO samples through the expression-based
##      array-to-sample links (never by column order)
##   2. log2 and quantile normalisation; probes kept if detected (detection p < 0.05) in at least as many samples
##      as the smallest group; probes mapped to Ensembl genes, keeping the most highly expressed probe per gene
##   3. limma, ~ country + group, with BeadChip as a blocking factor (duplicateCorrelation):
##        main       active vs latent TB in people with HIV
##        secondary  active TB vs other diseases in people with HIV; active vs latent TB without HIV
##        interaction (active vs latent TB) x (HIV positive vs negative): does HIV change the response to TB?
##   4. outputs in the same layout and column names as the RNA-seq datasets, so the network steps run unchanged
## Arrays compress fold changes compared with RNA-seq; the same thresholds as elsewhere (padj < 0.05, |log2FC| >= 1)
## are kept for comparability, and the threshold-free analyses (GSEA, evidence ranking) are unaffected.
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/GSE37250_run.R

ACC <- "GSE37250"
PADJ <- 0.05; LFC <- 1; DET_P <- 0.05; SEED <- 42
for (p in c("limma")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, ACC); RAW <- file.path(ROOT, "external_data", ACC); DES <- file.path(B, "00_design")
P <- function(...) { f <- file.path(B, ...); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); f }
rd <- function(f) { if (!file.exists(f)) stop("Missing: ", f, " (run the design-check and chip-mapping scripts first)"); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

## ---- samples ----
S <- rd(file.path(DES, paste0(ACC, "_sample_table.csv")))
M <- rd(file.path(DES, paste0(ACC, "_array_to_sample.csv"))); M <- M[M$accepted, ]
dis <- c("active tuberculosis" = "ATB", "latent TB infection" = "LTBI", "other disease" = "OD")
S$disease <- unname(dis[S$disease.state]); S$hiv <- ifelse(grepl("positive", S$hiv.status, ignore.case = TRUE), "HIVpos", "HIVneg")
S$group <- paste(S$disease, S$hiv, sep = "_"); S$country <- gsub(" ", "_", S$geographical.region)
if (anyNA(S$disease)) stop("Unexpected disease labels: ", paste(unique(S$disease.state[is.na(S$disease)]), collapse = ", "))
cat(">>> samples by group:\n"); print(table(S$group, S$country))

## ---- raw signal, renamed by the expression-based links ----
f <- list.files(RAW, pattern = "non.?normali[sz]ed", full.names = TRUE)[1]
R <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
probe <- R[[1]]
sig <- as.matrix(R[, grep("AVG_Signal", names(R))]); colnames(sig) <- sub("\\.AVG_Signal$", "", colnames(sig))
det <- as.matrix(R[, grep("Detection", names(R))]); colnames(det) <- sub("\\.Detection.*$", "", colnames(det))
rownames(sig) <- rownames(det) <- probe
gsm <- M$gsm[match(colnames(sig), M$array)]
if (anyNA(gsm) || anyDuplicated(gsm)) stop("Not every array has a unique accepted sample link.")
colnames(sig) <- gsm; colnames(det) <- M$gsm[match(colnames(det), M$array)]
S <- S[match(colnames(sig), S$gsm), ]; det <- det[, colnames(sig)]
stopifnot(identical(S$gsm, colnames(sig)))

## ---- normalisation and probe filter ----
x <- limma::normalizeBetweenArrays(log2(pmax(sig, 1)), method = "quantile")
n_min <- min(table(S$group))
keep <- rowSums(det < DET_P) >= n_min
x <- x[keep, ]; cat(">>> probes detected in at least", n_min, "samples:", sum(keep), "of", length(keep), "\n")

## ---- probes to Ensembl genes ----
ann_f <- file.path(DES, paste0(ACC, "_probe_annotation.csv"))
if (!file.exists(ann_f)) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("Missing package: GEOquery")
  gse <- GEOquery::getGEO(filename = file.path(RAW, paste0(ACC, "_series_matrix.txt.gz")), GSEMatrix = TRUE, getGPL = TRUE)
  fd <- Biobase::fData(gse)
  ec <- intersect(c("Entrez_Gene_ID", "ENTREZ_GENE_ID", "Entrez"), names(fd))[1]; sc <- intersect(c("Symbol", "ILMN_Gene", "SYMBOL"), names(fd))[1]
  ann <- data.frame(probe = rownames(fd), entrez = if (!is.na(ec)) as.character(fd[[ec]]) else NA, symbol = if (!is.na(sc)) as.character(fd[[sc]]) else NA)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Missing package: org.Hs.eg.db")
  e2 <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(stats::na.omit(ann$entrez[nzchar(ann$entrez)])), "ENSEMBL", "ENTREZID", multiVals = "first"))
  s2 <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(stats::na.omit(ann$symbol[nzchar(ann$symbol)])), "ENSEMBL", "SYMBOL", multiVals = "first"))
  ann$ensembl <- unname(e2[ann$entrez]); miss <- is.na(ann$ensembl); ann$ensembl[miss] <- unname(s2[ann$symbol[miss]])
  write.csv(ann, ann_f, row.names = FALSE)
}
ann <- read.csv(ann_f, stringsAsFactors = FALSE)
ens <- ann$ensembl[match(rownames(x), ann$probe)]
ok <- !is.na(ens) & nzchar(ens); x <- x[ok, ]; ens <- ens[ok]
o <- order(-rowMeans(x)); x <- x[o, ]; ens <- ens[o]; x <- x[!duplicated(ens), ]; rownames(x) <- ens[!duplicated(ens)]
cat(">>> genes after mapping to Ensembl (most expressed probe per gene):", nrow(x), "\n")

## ---- write the shared inputs ----
meta <- data.frame(sample_id = S$gsm, group = S$group, disease = S$disease, hiv = S$hiv, country = S$country, chip = S$chip)
write.csv(meta, P("02_metadata", paste0(ACC, "_sample_metadata.csv")), row.names = FALSE)
write.csv(x, P("04_normalized", paste0(ACC, "_vst_normalized_matrix.csv")))      # name kept for the pipeline; values are normalised log2 array intensities
write.csv(x, P("09_wgcna_input", paste0(ACC, "_WGCNA_expression_matrix.csv")))
tr <- sapply(sort(unique(S$group)), function(g) as.integer(S$group == g)); rownames(tr) <- S$gsm
write.csv(tr, P("09_wgcna_input", paste0(ACC, "_WGCNA_trait_file.csv")))

## ---- limma ----
de <- function(ctrl, case) {
  k <- S$group %in% c(ctrl, case); s <- S[k, ]; s$g <- factor(s$group, levels = c(ctrl, case))
  d <- model.matrix(~ country + g, s); y <- x[, k]
  dc <- limma::duplicateCorrelation(y, d, block = s$chip)
  fit <- limma::eBayes(limma::lmFit(y, d, block = s$chip, correlation = dc$consensus.correlation))
  t <- limma::topTable(fit, coef = ncol(d), number = Inf, sort.by = "none")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, ncol(d)]
  out <- data.frame(gene = rownames(t), baseMean = t$AveExpr, log2FoldChange = t$logFC, lfcSE = se, pvalue = t$P.Value, padj = t$adj.P.Val)
  out$significant <- out$padj < PADJ & abs(out$log2FoldChange) >= LFC
  out <- out[order(out$pvalue), ]
  name <- paste0(ctrl, "_vs_", case)
  write.csv(out, P("06_deg_results", paste0(ACC, "_", name, "_DEG.csv")), row.names = FALSE)
  write.csv(out[out$significant, ], P("07_significant_degs", paste0(ACC, "_", name, "_sigDEGs.csv")), row.names = FALSE)
  cat(sprintf("  %-28s n = %3d vs %3d | chip correlation %.3f | padj < %.2f: %5d | DEGs (|log2FC| >= %g): %4d (up %d, down %d)\n",
              name, sum(s$g == ctrl), sum(s$g == case), dc$consensus.correlation, PADJ, sum(out$padj < PADJ), LFC,
              sum(out$significant), sum(out$significant & out$log2FoldChange > 0), sum(out$significant & out$log2FoldChange < 0)))
  out
}
cat("\n========== differential expression (limma, ~ country + group, BeadChip as block) ==========\n")
main <- de("LTBI_HIVpos", "ATB_HIVpos")
invisible(de("OD_HIVpos", "ATB_HIVpos")); invisible(de("LTBI_HIVneg", "ATB_HIVneg"))

## ---- TB x HIV interaction ----
k <- S$disease %in% c("ATB", "LTBI"); s <- S[k, ]
s$TB <- factor(ifelse(s$disease == "ATB", "yes", "no"), levels = c("no", "yes")); s$HIV <- factor(ifelse(s$hiv == "HIVpos", "yes", "no"), levels = c("no", "yes"))
d <- model.matrix(~ country + TB * HIV, s); y <- x[, k]
dc <- limma::duplicateCorrelation(y, d, block = s$chip)
fit <- limma::eBayes(limma::lmFit(y, d, block = s$chip, correlation = dc$consensus.correlation))
ti <- limma::topTable(fit, coef = "TByes:HIVyes", number = Inf, sort.by = "none")
inter <- data.frame(gene = rownames(ti), baseMean = ti$AveExpr, log2FoldChange = ti$logFC, pvalue = ti$P.Value, padj = ti$adj.P.Val)
inter <- inter[order(inter$pvalue), ]
write.csv(inter, P("06_deg_results", paste0(ACC, "_interaction_TBxHIV.csv")), row.names = FALSE)   # not *_DEG.csv on purpose
cat(sprintf("\n  TB x HIV interaction (%d samples): %d genes with padj < %.2f\n", nrow(s), sum(inter$padj < PADJ), PADJ))

## ---- interferon gene sets: main comparison and interaction ----
if (all(vapply(c("clusterProfiler", "msigdbr", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(rownames(x), "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  gsea <- function(d) { d$stat <- sign(d$log2FoldChange) * -log10(pmax(d$pvalue, 1e-300))
    d <- merge(d[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g, pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0)) }
  gi <- gsea(inter); write.csv(gi[, c("ID", "setSize", "NES", "pvalue", "p.adjust")], P("06_deg_results", paste0(ACC, "_interaction_TBxHIV_GSEA_HALLMARK.csv")), row.names = FALSE)
  gm <- gsea(main)
  ifn <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
  tab <- merge(gm[gm$ID %in% ifn, c("ID", "NES", "p.adjust")], gi[gi$ID %in% ifn, c("ID", "NES", "p.adjust")], by = "ID", suffixes = c("_ATBvsLTBI_HIVpos", "_TBxHIV"))
  cat("\n========== interferon Hallmark sets ==========\n"); print(tab, row.names = FALSE)
  cat("  (a negative interaction NES means TB raises these genes less in people with HIV)\n")
} else cat("\n(clusterProfiler/msigdbr not installed: interferon GSEA skipped)\n")

cat("\n>>> DONE. Inputs for the network steps are in", B, "\n")
cat(">>> Add to scripts/wgcna_network/config.R, inside PRESETS:\n")
cat('    GSE37250_TBHIV = list(ACC = "GSE37250", CONTROL = "LTBI_HIVpos", CASE = "ATB_HIVpos",\n',
    '                         DEG_NAME = "LTBI_HIVpos_vs_ATB_HIVpos", RUN = "TB_HIV",\n',
    '                         MODULE_NAMES = character(0), MODULE_PLOT_COLORS = NULL),\n', sep = "")
