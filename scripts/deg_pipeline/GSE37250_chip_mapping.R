## GSE37250: recover the BeadChip (technical batch) of every sample.
##
## The non-normalised file names its columns by BeadChip barcode and position (e.g. 5483347021_A), but GEO does not
## say which barcode belongs to which sample, and matching by column order is exactly the mistake that mislabelled
## two cohorts earlier. Instead, each raw array is linked to a GEO sample by its expression profile: the raw signal
## is log2-transformed and correlated (Spearman, over all probes) with the processed values of every GEO sample in
## the series matrix. A link is accepted only when it is mutual (the array's best sample has that array as its best
## array) and clearly better than the second-best match. The chip is then the barcode before the underscore; an
## HT-12 BeadChip holds 12 arrays, so a chip is the natural technical batch.
## Run after GSE37250_design_check.R:  Rscript scripts/deg_pipeline/GSE37250_chip_mapping.R
## Outputs in GSE37250/00_design/.

ACC <- "GSE37250"
MIN_GAP <- 0.01        # best correlation must beat the second best by at least this much
for (p in c("GEOquery", "Biobase")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT, ACC, "00_design"); RAW <- file.path(ROOT, "external_data", ACC)
S <- read.csv(file.path(OUT, paste0(ACC, "_sample_table.csv")), stringsAsFactors = FALSE)

## ---- processed values per GEO sample ----
gse <- GEOquery::getGEO(filename = file.path(RAW, paste0(ACC, "_series_matrix.txt.gz")), GSEMatrix = TRUE)
P <- Biobase::exprs(gse)
if (max(P, na.rm = TRUE) > 100) P <- log2(pmax(P, 1))
cat(">>> series matrix:", nrow(P), "probes x", ncol(P), "samples\n")

## ---- raw signal per array ----
f <- list.files(RAW, pattern = "non.?normali[sz]ed", full.names = TRUE)[1]
R <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
id <- R[[1]]; sig <- grep("AVG_Signal", names(R), value = TRUE)
R <- as.matrix(R[, sig]); rownames(R) <- id; colnames(R) <- sub("\\.AVG_Signal$", "", sig)
R <- log2(pmax(R, 1))
cat(">>> raw file:", nrow(R), "probes x", ncol(R), "arrays\n")

## ---- correlate on shared probes (ranks, so the different normalisations do not matter) ----
pr <- intersect(rownames(R), rownames(P)); cat(">>> shared probes:", length(pr), "\n")
rk <- function(m) apply(m, 2, rank)
C <- stats::cor(rk(R[pr, ]), rk(P[pr, ]))          # arrays x samples, Spearman
best_s <- colnames(C)[apply(C, 1, which.max)]
best_a <- rownames(C)[apply(C, 2, which.max)]; names(best_a) <- colnames(C)
gap <- apply(C, 1, function(r) { s <- sort(r, decreasing = TRUE); s[1] - s[2] })
mutual <- best_a[best_s] == rownames(C)
M <- data.frame(array = rownames(C), gsm = best_s, rho = round(apply(C, 1, max), 4), gap = round(gap, 4),
                mutual = mutual, accepted = mutual & gap >= MIN_GAP)
M$chip <- sub("_.*$", "", M$array); M$position <- sub("^[^_]*_", "", M$array)
cat(">>> accepted links:", sum(M$accepted), "of", nrow(M), "| median best rho", round(stats::median(M$rho), 3),
    "| median gap to second best", round(stats::median(M$gap), 3), "\n")
if (any(duplicated(M$gsm[M$accepted]))) stop("Two arrays matched the same sample; the mapping is ambiguous.")
write.csv(M, file.path(OUT, paste0(ACC, "_array_to_sample.csv")), row.names = FALSE)

S$chip <- M$chip[match(S$gsm, M$gsm[M$accepted])]
write.csv(S, file.path(OUT, paste0(ACC, "_sample_table.csv")), row.names = FALSE)
cat(">>> samples with a chip:", sum(!is.na(S$chip)), "of", nrow(S), "| chips:", length(unique(stats::na.omit(S$chip))), "\n")

## ---- chips against groups ----
S$group <- paste(S$disease.state, "|", S$hiv.status)
per_chip <- table(S$chip)
cat("\n>>> arrays per chip:", paste(names(table(per_chip)), table(per_chip), sep = " x", collapse = ", "), "\n")
mixed <- tapply(S$group, S$chip, function(g) length(unique(g)))
cat(">>> chips holding more than one group:", sum(mixed > 1), "of", length(mixed), "\n")
ctry <- tapply(S$geographical.region, S$chip, function(g) length(unique(g)))
cat(">>> chips holding samples from both countries:", sum(ctry > 1), "of", length(ctry), "\n")
for (lab in c("HIV positive", "HIV negative")) {
  s <- S[!is.na(S$chip) & grepl(lab, S$hiv.status) & S$disease.state %in% c("active tuberculosis", "latent TB infection"), ]
  if (nrow(s) < 10) next
  tb <- table(s$chip, s$disease.state); x <- stats::chisq.test(tb, simulate.p.value = TRUE, B = 5000)
  shared <- sum(rowSums(tb > 0) == 2)
  cat(sprintf("  %s, active vs latent TB: %d chips, %d contain both groups; chip-group association p = %.3g\n",
              lab, nrow(tb), shared, x$p.value))
}
write.csv(as.data.frame.matrix(table(S$chip, S$group)), file.path(OUT, paste0(ACC, "_chip_by_group.csv")))
cat("\n>>> DONE. Outputs in", OUT, "\n")
