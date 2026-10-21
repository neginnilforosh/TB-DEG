## Validation of the TB-HIV context in two independent whole-blood RNA-seq cohorts of people with HIV.
##
##   GSE107104  Uganda, advanced HIV: TB-HIV vs HIV only, ART-naive participants only (the 2 on ART are left out).
##              GEO provides normalised counts only, so log2(x + 1) values are analysed with limma-trend.
##   GSE162164  India, advanced HIV: TB-HIV vs HIV only. Raw counts are analysed with DESeq2 (apeglm shrinkage),
##              adjusted for sex.
## Samples are linked to their groups by name only: through the GSM accession, the sample ID shown in the GEO title,
## or the per-sample supplementary file name; if none works the script stops rather than guessing.
## The discovery comparison (GSE37250) is active vs latent TB in people with HIV, while these cohorts compare TB-HIV
## with HIV without TB, so agreement is read as consistency of the TB signal in HIV, not as an exact replication.
## For each cohort:
##   - direction and significance of the prioritised TB-HIV genes (evidence-family tiers) and of the 8 core genes
##   - signature score (mean z of up genes minus mean z of down genes) AUC, bootstrap CI, and 1,000 random gene sets
##   - Hallmark GSEA (interferon sets reported)
## Run after the TB-HIV pipeline:  Rscript scripts/deg_pipeline/tbhiv_validation.R
## Outputs in GSE37250/TB_HIV/23_hiv_validation/.

CORE <- c("GBP1", "PARP9", "FCGR1A", "SERPING1", "AIM2", "CASP1", "TNFSF10", "FBXO6")
PADJ <- 0.05; SEED <- 42; N_RANDOM <- 1000; N_BOOT <- 1000
for (p in c("limma", "DESeq2", "readxl", "GEOquery", "Biobase", "org.Hs.eg.db")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT, "GSE37250", "TB_HIV", "23_hiv_validation"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ev <- read.csv(file.path(ROOT, "GSE37250", "TB_HIV", "20_evidence_families", "GSE37250_Gene_evidence_families.csv"), stringsAsFactors = FALSE)
ev <- ev[!is.na(ev$log2FC), ]

## ---------- helpers ----------
pheno <- function(acc) {
  gse <- GEOquery::getGEO(filename = file.path(ROOT, "external_data", acc, paste0(acc, "_series_matrix.txt.gz")), GSEMatrix = TRUE, getGPL = FALSE)
  Biobase::pData(gse)
}
## link expression columns to GSM accessions by name: GSM in the column, ID in brackets in the title, or the column
## name inside the per-sample supplementary file name
## sample IDs are compared after normalisation: split into letter and number parts, drop separators and leading
## zeros, so "IR1-0010" and "IR-1-010" are the same ID, while "IR1-001" and "IR1-010" stay different
norm_id <- function(x) vapply(regmatches(toupper(x), gregexpr("[A-Z]+|[0-9]+", toupper(x))), function(p) {
  p <- ifelse(grepl("^[0-9]+$", p), suppressWarnings(as.character(as.integer(p))), p); paste(p, collapse = "_") }, "")
link_columns <- function(cols, pd) {
  sup <- apply(pd[, grep("^supplementary_file", names(pd)), drop = FALSE], 1, paste, collapse = " ")
  tid <- sub(".*\\[(.*)\\].*", "\\1", pd$title); tid[tid == pd$title] <- NA
  ntid <- ifelse(is.na(tid), NA, norm_id(ifelse(is.na(tid), "", tid)))
  vapply(cols, function(cn) {
    j <- which(pd$geo_accession == cn); if (length(j) == 1) return(pd$geo_accession[j])
    j <- which(!is.na(ntid) & ntid == norm_id(cn)); if (length(j) == 1) return(pd$geo_accession[j])
    j <- which(vapply(sup, function(s) grepl(cn, s, fixed = TRUE), TRUE)); if (length(j) == 1) return(pd$geo_accession[j])
    j <- which(vapply(pd$title, function(t) grepl(cn, t, fixed = TRUE), TRUE)); if (length(j) == 1) return(pd$geo_accession[j])
    NA_character_ }, "")
}
to_ensembl <- function(ids) {
  ids <- sub("\\..*$", "", ids)
  if (mean(grepl("^ENSG", ids)) > 0.5) return(ids)
  type <- if (mean(grepl("^[0-9]+$", ids)) > 0.5) "ENTREZID" else "SYMBOL"
  unname(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ids, "ENSEMBL", type, multiVals = "first")))
}
collapse <- function(m, ens) { ok <- !is.na(ens); m <- m[ok, , drop = FALSE]; ens <- ens[ok]
  o <- order(-rowMeans(m)); m <- m[o, , drop = FALSE]; ens <- ens[o]; m <- m[!duplicated(ens), , drop = FALSE]; rownames(m) <- ens[!duplicated(ens)]; m }
auc <- function(s, y) { r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

validate <- function(acc, de, logx, y) {
  cat("\n  -- replication of the prioritised TB-HIV genes --\n")
  m <- merge(ev[, c("Ensembl", "Gene", "log2FC", "Evidence_tier")], de, by.x = "Ensembl", by.y = "gene")
  m$same <- sign(m$log2FC) == sign(m$log2FoldChange); m$rep <- m$same & m$padj < PADJ
  tiers <- do.call(rbind, lapply(c("High", "Medium", "Low"), function(t) { s <- m[m$Evidence_tier == t, ]; if (!nrow(s)) return(NULL)
    data.frame(Dataset = acc, Tier = t, N = nrow(s), Pct_same_direction = round(100 * mean(s$same), 1),
               P_direction = signif(stats::binom.test(sum(s$same), nrow(s), 0.5, alternative = "greater")$p.value, 3),
               Pct_reproduced = round(100 * mean(s$rep), 1)) }))
  h <- m$rep[m$Evidence_tier == "High"]; l <- m$rep[m$Evidence_tier == "Low"]
  tiers$P_High_vs_Low <- NA
  if (length(h) && length(l)) tiers$P_High_vs_Low[tiers$Tier == "High"] <- signif(stats::fisher.test(matrix(c(sum(h), sum(!h), sum(l), sum(!l)), 2), alternative = "greater")$p.value, 3)
  print(tiers, row.names = FALSE)
  core <- m[m$Gene %in% CORE, c("Gene", "log2FC", "log2FoldChange", "padj")]; names(core)[2:3] <- c("log2FC_discovery", "log2FC_validation")
  cat("\n  core genes:\n"); print(core[order(core$padj), ], row.names = FALSE)
  write.csv(m, file.path(OUT, paste0(acc, "_gene_replication.csv")), row.names = FALSE)

  cat("\n  -- signature score (High genes, directions from the discovery cohort) --\n")
  hg <- ev[ev$Evidence_tier == "High", ]; hg <- hg[hg$Ensembl %in% rownames(logx), ]
  z <- t(scale(t(logx))); z[is.na(z)] <- 0
  up <- hg$Ensembl[hg$log2FC > 0]; dn <- hg$Ensembl[hg$log2FC < 0]
  score <- (if (length(up)) colMeans(z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(z[dn, , drop = FALSE]) else 0)
  a <- auc(score, y); set.seed(SEED)
  bt <- replicate(N_BOOT, { i1 <- sample(which(y == 1), replace = TRUE); i0 <- sample(which(y == 0), replace = TRUE); auc(score[c(i1, i0)], y[c(i1, i0)]) })
  pool <- setdiff(rownames(z), hg$Ensembl)
  null <- replicate(N_RANDOM, { s <- colMeans(z[sample(pool, nrow(hg)), , drop = FALSE]); max(auc(s, y), auc(-s, y)) })
  res <- data.frame(Dataset = acc, Genes = nrow(hg), AUC = round(a, 3), CI_low = round(stats::quantile(bt, .025), 3), CI_high = round(stats::quantile(bt, .975), 3),
                    Random_median = round(stats::median(null), 3), P_vs_random = signif((sum(null >= a) + 1) / (N_RANDOM + 1), 3))
  print(res, row.names = FALSE)

  if (all(vapply(c("clusterProfiler", "msigdbr"), requireNamespace, logical(1), quietly = TRUE))) {
    hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
    map <- suppressMessages(clusterProfiler::bitr(de$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    d <- de; d$stat <- sign(d$log2FoldChange) * -log10(pmax(d$pvalue, 1e-300))
    d <- merge(d[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
    set.seed(SEED); g <- as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g, pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0))
    write.csv(g[, c("ID", "setSize", "NES", "pvalue", "p.adjust")], file.path(OUT, paste0(acc, "_GSEA_HALLMARK.csv")), row.names = FALSE)
    cat("\n  -- Hallmark GSEA, TB-HIV vs HIV (interferon and the five strongest) --\n")
    print(g[g$ID %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE") | rank(g$p.adjust) <= 5, c("ID", "NES", "p.adjust")], row.names = FALSE)
  }
  list(tiers = tiers, auc = res)
}

all_tiers <- list(); all_auc <- list()

## ---------- GSE107104 ----------
cat("\n================ GSE107104 (Uganda, ART-naive) ================\n")
pd <- pheno("GSE107104")
grp <- ifelse(grepl("TB", pd$`disease status:ch1`), "TB_HIV", "HIV")
naive <- grepl("naive|naïve", pd$`treatment:ch1`, ignore.case = TRUE)
x <- as.data.frame(readxl::read_excel(file.path(ROOT, "external_data", "GSE107104", "GSE107104_normalizedCounts.xlsx")))
idc <- names(x)[!vapply(x, is.numeric, TRUE)][1]; ids <- x[[idc]]; x <- as.matrix(x[, vapply(x, is.numeric, TRUE)])
gsm <- link_columns(colnames(x), pd)
cat(">>> columns linked to GEO samples by name:", sum(!is.na(gsm)), "of", ncol(x), "| id column:", idc, "| example ids:", paste(utils::head(ids, 3), collapse = ", "), "\n")
if (anyNA(gsm)) stop("Unlinked columns: ", paste(colnames(x)[is.na(gsm)], collapse = ", "))
colnames(x) <- gsm; k <- match(gsm, pd$geo_accession)
keep <- naive[k]; y <- as.integer(grp[k] == "TB_HIV")[keep]; x <- x[, keep]
cat(">>> ART-naive samples: TB-HIV", sum(y == 1), "| HIV only", sum(y == 0), "\n")
lx <- collapse(log2(x + 1), to_ensembl(ids)); lx <- lx[rowMeans(lx) > 1, ]
des <- model.matrix(~ factor(y)); fit <- limma::eBayes(limma::lmFit(lx, des), trend = TRUE)
t <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
de <- data.frame(gene = rownames(t), log2FoldChange = t$logFC, pvalue = t$P.Value, padj = t$adj.P.Val)
write.csv(de[order(de$pvalue), ], file.path(OUT, "GSE107104_TBHIV_vs_HIV_limma.csv"), row.names = FALSE)
cat(">>> genes:", nrow(lx), "| padj < 0.05:", sum(de$padj < PADJ), "\n")
r <- validate("GSE107104", de, lx, y); all_tiers[[1]] <- r$tiers; all_auc[[1]] <- r$auc

## ---------- GSE162164 ----------
cat("\n================ GSE162164 (India) ================\n")
pd <- pheno("GSE162164")
grp <- ifelse(grepl("TB", pd$title), "TB_HIV", "HIV")
f <- file.path(ROOT, "external_data", "GSE162164", "GSE162164_counts.txt.gz")
cnt <- as.matrix(utils::read.table(f, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
gsm <- link_columns(colnames(cnt), pd)
cat(">>> columns linked to GEO samples by name:", sum(!is.na(gsm)), "of", ncol(cnt), "| example ids:", paste(utils::head(rownames(cnt), 3), collapse = ", "), "\n")
## The combined table uses laboratory IDs (RICC_...) that GEO never ties to samples, but GEO also provides one count
## file per sample, named by its GSM accession. Each combined column is linked to the per-sample file that holds
## exactly the same counts: identity of the data, not column order, decides the link.
tarf <- file.path(ROOT, "external_data", "GSE162164", "GSE162164_RAW.tar")
if (anyNA(gsm) && file.exists(tarf)) {
  td <- file.path(tempdir(), "GSE162164_RAW"); dir.create(td, showWarnings = FALSE); utils::untar(tarf, exdir = td)
  files <- list.files(td, pattern = "^GSM", full.names = TRUE)
  per <- lapply(files, function(f) { d <- utils::read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE)
    num <- which(vapply(d, function(v) all(!is.na(suppressWarnings(as.numeric(v[-1])))), TRUE)); num <- num[num > 1]
    if (!length(num)) return(NULL)
    v <- suppressWarnings(as.numeric(d[[num[length(num)]]])); names(v) <- d[[1]]; v[!is.na(v)] })
  names(per) <- sub("_.*$", "", basename(files))
  match_one <- function(j) {
    hit <- vapply(per, function(v) { if (is.null(v)) return(FALSE); g <- intersect(names(v), rownames(cnt)); length(g) > 1000 && all(v[g] == cnt[g, j]) }, TRUE)
    if (sum(hit) == 1) names(per)[hit] else NA_character_ }
  by_content <- vapply(seq_len(ncol(cnt)), match_one, "")
  cat(">>> columns linked by identical counts in the per-sample GEO files:", sum(!is.na(by_content)), "of", ncol(cnt), "\n")
  gsm[is.na(gsm)] <- by_content[is.na(gsm)]
  if (anyDuplicated(stats::na.omit(gsm))) stop("Two columns matched the same sample; the link is ambiguous.")
  write.csv(data.frame(column = colnames(cnt), gsm = gsm, title = pd$title[match(gsm, pd$geo_accession)]),
            file.path(OUT, "GSE162164_column_to_sample.csv"), row.names = FALSE)
}
if (anyNA(gsm)) stop("Unlinked columns: ", paste(colnames(cnt)[is.na(gsm)], collapse = ", "))
colnames(cnt) <- gsm; k <- match(gsm, pd$geo_accession)
sx <- pd$`Sex:ch1`[k]; if (is.null(sx)) sx <- rep("NA", length(k))
y <- as.integer(grp[k] == "TB_HIV"); cat(">>> TB-HIV", sum(y == 1), "| HIV only", sum(y == 0), "\n")
ens <- to_ensembl(rownames(cnt)); ok <- !is.na(ens); cnt <- rowsum(cnt[ok, ], ens[ok])
cnt <- cnt[rowSums(cnt >= 10) >= min(table(y)), ]
md <- data.frame(row.names = colnames(cnt), y = factor(y), sex = factor(sx))
dds <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(round(cnt), md, if (length(unique(sx)) > 1) ~ sex + y else ~ y), quiet = TRUE)
rs <- DESeq2::lfcShrink(dds, coef = "y_1_vs_0", type = if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal", quiet = TRUE)
de <- data.frame(gene = rownames(rs), log2FoldChange = rs$log2FoldChange, pvalue = rs$pvalue, padj = rs$padj); de <- de[!is.na(de$padj), ]
write.csv(de[order(de$pvalue), ], file.path(OUT, "GSE162164_TBHIV_vs_HIV_DESeq2.csv"), row.names = FALSE)
cat(">>> genes:", nrow(cnt), "| padj < 0.05:", sum(de$padj < PADJ), "\n")
vs <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
r <- validate("GSE162164", de, vs, y); all_tiers[[2]] <- r$tiers; all_auc[[2]] <- r$auc

write.csv(do.call(rbind, all_tiers), file.path(OUT, "TBHIV_validation_replication_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_auc), file.path(OUT, "TBHIV_validation_signature_AUC.csv"), row.names = FALSE)
cat("\n>>> DONE. Outputs in", OUT, "\n")
