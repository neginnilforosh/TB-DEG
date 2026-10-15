## Metabolomics of TB progression in the GC6-74 cohort (ST000974 plasma, ST000975 serum, ST000976 diluted plasma).
##
## Cases (progressors) are household contacts who developed active TB during follow-up; controls (non-progressors)
## stayed TB-free. Each case was matched to about four controls (CASE_GROUP = matched set), samples were taken at
## enrolment (BL) and 6 and 18 months later while participants were TB-free, DONOR_ID links repeated samples of one
## person, and TIME_TO_TB gives the months from a sample to the case's diagnosis. Each study is analysed on its own.
##
##  A. Main analysis, baseline: progressors vs non-progressors at BL, compared within matched sets
##     (~ matched set + group); sensitivity without matching (~ site + sex + age + group).
##  B. Closeness to diagnosis: all pre-diagnosis samples; case samples split into < 5 and >= 5 months before
##     diagnosis and each compared with control samples (~ group + time point + site + sex + age), with repeated
##     samples of the same person handled as correlated (limma duplicateCorrelation, block = DONOR_ID).
##  The kynurenine/tryptophan ratio, a read-out of IDO1 activity, is added as an extra feature when both are measured.
## Samples taken after diagnosis (TIME_TO_TB < 0) are excluded. ST000976 is small and reported as descriptive only.
## Needs scripts/wgcna_network/metabolomics_metadata_check.R to have been run (it writes the sample table).
## Run from anywhere in the repository:  Rscript scripts/wgcna_network/metabolomics_progression.R
## Outputs in GSE114192/21_metabolomics/progression/.

## ---- CONFIG ----
STUDIES   <- c(ST000974 = "Plasma", ST000975 = "Serum", ST000976 = "Plasma diluted with RPMI")
DESCRIPTIVE_ONLY <- c("ST000976")
MIN_FRAC  <- 0.8      # keep metabolites measured (non-missing, > 0) in at least 80% of the analysed samples
NEAR_MONTHS <- 5      # "< 5 months before diagnosis" as in the original study
FDR <- 0.05
AMINO <- c("tryptophan", "kynurenine", "histidine", "leucine", "isoleucine", "valine", "methionine", "phenylalanine",
           "tyrosine", "threonine", "serine", "glycine", "alanine", "asparagine", "glutamine", "glutamate", "aspartate",
           "lysine", "arginine", "proline", "kynurenine/tryptophan ratio")

## smoking sensitivity: MET_SMOKING=1 adds the SMOKE field as a covariate everywhere and writes to progression_smoking/
SMOKING <- identical(Sys.getenv("MET_SMOKING"), "1")
EXTRA <- if (SMOKING) "SMOKE" else character(0)
for (p in c("limma", "httr", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
MET  <- file.path(ROOT, "GSE114192", "21_metabolomics")
OUT  <- file.path(MET, if (SMOKING) "progression_smoking" else "progression"); dir.create(OUT, showWarnings = FALSE)
CACHE <- file.path(MET, "progression")   # abundance matrices are shared between the two runs
if (SMOKING) cat(">>> smoking sensitivity run: SMOKE added as a covariate\n")
st_file <- file.path(MET, "metadata_check", "GC6-74_metabolomics_sample_table.csv")
if (!file.exists(st_file)) stop("Run metabolomics_metadata_check.R first (missing ", st_file, ").")
S_ALL <- read.csv(st_file, stringsAsFactors = FALSE)

## ---- abundance matrix from the REST service (cached) ----
get_matrix <- function(st) {
  dir.create(CACHE, showWarnings = FALSE); cache <- file.path(CACHE, paste0(st, "_abundance_matrix.csv"))
  if (file.exists(cache)) { m <- as.matrix(read.csv(cache, row.names = 1, check.names = FALSE)); return(m) }
  r <- httr::GET(paste0("https://www.metabolomicsworkbench.org/rest/study/study_id/", st, "/data"), httr::timeout(300))
  if (httr::status_code(r) != 200) stop("Could not download data for ", st)
  ## parse without simplification: one record per metabolite, with a DATA field of sample -> value.
  ## Records can have empty fields (e.g. no RefMet name), which breaks data-frame conversion, so build it by hand.
  recs <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  if (!is.null(recs$DATA)) recs <- list(recs)
  val <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1]])
  nm <- vapply(recs, function(x) { n <- val(x$metabolite_name); if (is.na(n) || !nzchar(n)) val(x$refmet_name) else n }, "")
  vals <- lapply(recs, function(x) { d <- x$DATA; if (is.null(d)) return(numeric(0))
    v <- vapply(d, function(z) if (is.null(z) || !length(z) || identical(z, "")) NA_real_ else suppressWarnings(as.numeric(z[[1]])), 0); v })
  samples <- unique(unlist(lapply(vals, names)))
  if (!length(samples)) stop("No per-sample values found for ", st, "; record fields: ", paste(names(recs[[1]]), collapse = ", "))
  m <- t(vapply(vals, function(v) { out <- rep(NA_real_, length(samples)); names(out) <- samples; out[names(v)] <- v; out }, numeric(length(samples))))
  ok <- !is.na(nm) & nzchar(nm); m <- m[ok, , drop = FALSE]; rownames(m) <- make.unique(nm[ok]); colnames(m) <- samples
  cat("  downloaded", nrow(m), "metabolites x", ncol(m), "samples\n")
  write.csv(m, cache); m
}

fit_limma <- function(x, design, coefs, block = NULL) {
  cor <- NULL
  if (!is.null(block)) { dc <- limma::duplicateCorrelation(x, design, block = block); cor <- dc$consensus.correlation }
  fit <- limma::eBayes(limma::lmFit(x, design, block = block, correlation = cor))
  list(fit = fit, cor = cor, tabs = lapply(setNames(coefs, coefs), function(k) {
    t <- limma::topTable(fit, coef = k, number = Inf, sort.by = "none"); data.frame(metabolite = rownames(t), log2FC = t$logFC, pvalue = t$P.Value, padj = t$adj.P.Val) }))
}
key_rows <- function(t, study, analysis, contrast) {
  k <- t[tolower(t$metabolite) %in% AMINO, ]
  if (!nrow(k)) return(NULL)
  data.frame(Study = study, Analysis = analysis, Contrast = contrast, Metabolite = k$metabolite,
             log2FC = round(k$log2FC, 3), pvalue = signif(k$pvalue, 3), padj = signif(k$padj, 3))
}

summary_rows <- list(); key <- list()
for (st in names(STUDIES)) {
  cat("\n==========", st, "(", STUDIES[[st]], ")", if (st %in% DESCRIPTIVE_ONLY) "[descriptive only]" else "", "==========\n")
  s <- S_ALL[S_ALL$study == st, ]; s$TIME_TO_TB <- suppressWarnings(as.numeric(s$TIME_TO_TB)); s$AGE <- suppressWarnings(as.numeric(s$AGE))
  m <- tryCatch(get_matrix(st), error = function(e) { cat("  !!", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) next
  s <- s[s$sample %in% colnames(m), ]
  post <- s$GROUP == "case" & !is.na(s$TIME_TO_TB) & s$TIME_TO_TB < 0
  if (any(post)) cat("  excluded", sum(post), "case sample(s) taken after diagnosis (TIME_TO_TB < 0)\n")
  s <- s[!post, ]
  if (SMOKING) { ok <- !is.na(s$SMOKE) & nzchar(s$SMOKE); cat("  smoking status:", paste(names(table(s$SMOKE[ok])), table(s$SMOKE[ok]), collapse = ", "),
                 "| samples without it dropped:", sum(!ok), "\n"); s <- s[ok, ] }
  x <- log2(m[, s$sample, drop = FALSE]); x[!is.finite(x)] <- NA
  ## kynurenine / tryptophan ratio (log2 scale: difference)
  ik <- which(tolower(rownames(x)) == "kynurenine"); it <- which(tolower(rownames(x)) == "tryptophan")
  if (length(ik) == 1 && length(it) == 1) x <- rbind(x, "kynurenine/tryptophan ratio" = x[ik, ] - x[it, ])
  cat("  samples:", ncol(x), "| metabolites:", nrow(x), "| participants:", length(unique(s$subject)), "\n")

  ## ---------- A. baseline, within matched sets ----------
  b <- s[s$TIMEPOINT == "BL", ]
  sets <- tapply(b$GROUP, b$CASE_GROUP, function(g) all(c("case", "control") %in% g)); b <- b[b$CASE_GROUP %in% names(sets)[sets], ]
  xb <- x[, b$sample, drop = FALSE]; xb <- xb[rowMeans(!is.na(xb)) >= MIN_FRAC, , drop = FALSE]
  b$GROUP <- factor(b$GROUP, levels = c("control", "case"))
  cat("  A. baseline:", sum(b$GROUP == "case"), "cases,", sum(b$GROUP == "control"), "controls in", length(unique(b$CASE_GROUP)),
      "matched sets |", nrow(xb), "metabolites\n")
  if (sum(b$GROUP == "case") >= 3) {
    rA <- fit_limma(xb, model.matrix(stats::as.formula(paste("~ factor(CASE_GROUP) +", paste(c(EXTRA, "GROUP"), collapse = " + "))), b), "GROUPcase")$tabs[[1]]
    cv0 <- c("REGION", "GENDER", "AGE", EXTRA); covs <- cv0[vapply(cv0, function(v) length(unique(stats::na.omit(b[[v]]))) > 1, TRUE)]
    okc <- stats::complete.cases(b[, covs, drop = FALSE])
    rA2 <- fit_limma(xb[, okc], model.matrix(stats::as.formula(paste("~", paste(c(covs, "GROUP"), collapse = " + "))), b[okc, ]), "GROUPcase")$tabs[[1]]
    write.csv(rA[order(rA$pvalue), ], file.path(OUT, paste0(st, "_A_baseline_matched.csv")), row.names = FALSE)
    write.csv(rA2[order(rA2$pvalue), ], file.path(OUT, paste0(st, "_A_baseline_unmatched_adjusted.csv")), row.names = FALSE)
    cat("     significant (FDR <", FDR, "): matched", sum(rA$padj < FDR, na.rm = TRUE), "| unmatched + site/sex/age", sum(rA2$padj < FDR, na.rm = TRUE), "\n")
    top <- utils::head(rA[order(rA$pvalue), ], 5); cat("     top (matched):", paste(sprintf("%s %+.2f (p %.1e)", top$metabolite, top$log2FC, top$pvalue), collapse = "; "), "\n")
    key[[length(key) + 1]] <- key_rows(rA, st, "A baseline", "progressor vs non-progressor (matched sets)")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(Study = st, Analysis = "A baseline (matched)", Cases = sum(b$GROUP == "case"),
      Controls = sum(b$GROUP == "control"), Metabolites = nrow(xb), Significant_FDR = sum(rA$padj < FDR, na.rm = TRUE))
  }

  ## ---------- B. closeness to diagnosis, repeated samples handled ----------
  s$grp3 <- ifelse(s$GROUP == "control", "control", ifelse(!is.na(s$TIME_TO_TB) & s$TIME_TO_TB < NEAR_MONTHS, "near", "far"))
  s <- s[!(s$GROUP == "case" & is.na(s$TIME_TO_TB)), ]
  s$grp3 <- factor(s$grp3, levels = c("control", "far", "near"))
  cv1 <- c("TIMEPOINT", "REGION", "GENDER", "AGE", EXTRA); covs <- cv1[vapply(cv1, function(v) length(unique(stats::na.omit(s[[v]]))) > 1, TRUE)]
  s <- s[stats::complete.cases(s[, covs, drop = FALSE]), ]
  xs <- x[, s$sample, drop = FALSE]; xs <- xs[rowMeans(!is.na(xs)) >= MIN_FRAC, , drop = FALSE]
  cat("  B. closeness to diagnosis: case samples <", NEAR_MONTHS, "months:", sum(s$grp3 == "near"), "| >=", NEAR_MONTHS, "months:",
      sum(s$grp3 == "far"), "| control samples:", sum(s$grp3 == "control"), "\n")
  if (sum(s$grp3 == "near") >= 3 && sum(s$grp3 == "far") >= 3) {
    d <- model.matrix(stats::as.formula(paste("~ 0 + grp3 +", paste(covs, collapse = " + "))), s)
    colnames(d) <- make.names(colnames(d))
    dc <- limma::duplicateCorrelation(xs, d, block = s$subject)
    fit <- limma::lmFit(xs, d, block = s$subject, correlation = dc$consensus.correlation)
    cm <- limma::makeContrasts(near_vs_control = grp3near - grp3control, far_vs_control = grp3far - grp3control,
                               near_vs_far = grp3near - grp3far, levels = d)
    f2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
    cat("     within-person correlation:", round(dc$consensus.correlation, 3), "\n")
    for (k in colnames(cm)) {
      t <- limma::topTable(f2, coef = k, number = Inf, sort.by = "none")
      t <- data.frame(metabolite = rownames(t), log2FC = t$logFC, pvalue = t$P.Value, padj = t$adj.P.Val)
      write.csv(t[order(t$pvalue), ], file.path(OUT, paste0(st, "_B_", k, ".csv")), row.names = FALSE)
      cat(sprintf("     %-16s significant (FDR < %.2f): %d\n", k, FDR, sum(t$padj < FDR, na.rm = TRUE)))
      key[[length(key) + 1]] <- key_rows(t, st, "B closeness to diagnosis", k)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Study = st, Analysis = paste("B", k),
        Cases = sum(s$grp3 == ifelse(k == "far_vs_control", "far", "near")), Controls = if (k == "near_vs_far") sum(s$grp3 == "far") else sum(s$grp3 == "control"),
        Metabolites = nrow(xs), Significant_FDR = sum(t$padj < FDR, na.rm = TRUE))
    }
  } else cat("     too few case samples in one of the two windows; skipped\n")
}

if (length(summary_rows)) {
  S <- do.call(rbind, summary_rows); S$Note <- ifelse(S$Study %in% DESCRIPTIVE_ONLY, "descriptive only (small study)", "")
  write.csv(S, file.path(OUT, "GC6-74_progression_summary.csv"), row.names = FALSE)
  cat("\n========== summary ==========\n"); print(S, row.names = FALSE)
}
if (length(key)) {
  K <- do.call(rbind, key); write.csv(K, file.path(OUT, "GC6-74_amino_acids_and_kynurenine.csv"), row.names = FALSE)
  cat("\n========== tryptophan, kynurenine and their ratio ==========\n")
  print(K[tolower(K$Metabolite) %in% c("tryptophan", "kynurenine", "kynurenine/tryptophan ratio"), ], row.names = FALSE)
  if (SMOKING) { sm <- list.files(OUT, pattern = "_A_baseline_matched\\.csv$", full.names = TRUE)
    for (f in sm) { t <- read.csv(f); t <- t[grepl("naphthyl|cotinine|hydroxycotinine", t$metabolite, ignore.case = TRUE), ]
      if (nrow(t)) { cat("\n  smoking markers,", basename(f), "\n"); print(t, row.names = FALSE) } } }
}
cat("\nThe three studies share some participants (ST000974-ST000975 and ST000975-ST000976), so they are not independent\n",
    "replicates; read agreement between them as consistency across sample types, not as validation.\n", sep = "")
cat(">>> DONE. Outputs in", OUT, "\n")
