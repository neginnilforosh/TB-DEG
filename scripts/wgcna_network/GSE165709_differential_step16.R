## STEP 16: GSE165709 (alveolar macrophages) - response to Mtb infection, per cohort.
## Set ASSAY to "RNA" (GSE165708) or "ATAC" (GSE165703) and run it twice.
## Three results per assay: the infection effect in HC, the infection effect in PLWH, and the
## difference between the two responses (the genes/regions that respond in HC but not in PLWH).
## ATAC peaks are additionally annotated to genes. Outputs in external_results/GSE165709/.

## ---- CONFIG ----
ASSAY        <- "ATAC"             # "RNA" or "ATAC"
COUNT_FILE   <- list(RNA  = "GSE165708_non-normalized_estimated_counts_matrix.txt",
                     ATAC = "GSE165703_ATAC_FeatureCount_raw_quantification_matrix.txt")
KEEP_COHORTS <- c("HC", "PLWH")   # PrEP (treatment) is left out of the main analysis
CASE_LEVEL   <- "INF"; CTRL_LEVEL <- "NEG"
MIN_COUNT    <- 10
PADJ_CUT     <- 0.05
TSS_REGION   <- c(-3000, 3000)

req <- c("DESeq2")
if (ASSAY == "ATAC") req <- c(req, "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomicRanges", "IRanges")
ok <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing package(s): ", paste(req[!ok], collapse = ", "),
                   "\n  BiocManager::install(c(", paste0('"', req[!ok], '"', collapse = ", "), "))")

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
ROOT_DIR <- local({ d <- SCRIPT_DIR
  while (length(list.files(d, pattern = "^GSE[0-9]+$")) == 0 && dirname(d) != d) d <- dirname(d); d })
dir_out <- file.path(ROOT_DIR, "external_results", "GSE165709", ASSAY)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
count_path <- file.path(ROOT_DIR, "external_data", "GSE165709", COUNT_FILE[[ASSAY]])
if (!file.exists(count_path)) stop("Count matrix not found: ", count_path,
                                   "\n  Download it from GEO into external_data/GSE165709/ (that folder is git-ignored).")

## ---- 1. counts and sample table ----
counts <- read.delim(count_path, header = TRUE, row.names = 1, check.names = FALSE)   # both files are tab-separated
cat(">>> ", ASSAY, ": ", nrow(counts), " features x ", ncol(counts), " samples\n", sep = "")
samples <- colnames(counts)
## The two subseries name their samples differently, so each is parsed on its own terms.
## ATAC  "AMC12_Mtb_INF": prefix gives the cohort (AMC = HC, AMP = PrEP, AMV = PLWH).
## RNA   "HC.AMC2.MTB.1": the cohort is written out, and the last field is 1 = Mtb-challenged,
##       0 = non-challenged. Note that AMC covers BOTH HC (AMC2-22) and PrEP (AMC23-37) here,
##       so the prefix must NOT be used to set the cohort for the RNA file.
meta <- if (ASSAY == "ATAC") {
  d <- sub("_Mtb_(INF|NEG)$", "", samples)
  data.frame(sample = samples,
             cohort = c(AMC = "HC", AMP = "PrEP", AMV = "PLWH")[sub("[0-9_].*$", "", d)],
             infection = ifelse(grepl("_Mtb_INF$", samples), CASE_LEVEL, CTRL_LEVEL),
             donor = sub("_(2|LO|DUP2)$", "", d),       # technical variants of the same donor
             stringsAsFactors = FALSE)
} else {
  p <- strsplit(samples, ".", fixed = TRUE)
  data.frame(sample = samples,
             cohort = vapply(p, `[`, "", 1),
             infection = ifelse(vapply(p, function(x) x[length(x)], "") == "1", CASE_LEVEL, CTRL_LEVEL),
             donor = vapply(p, `[`, "", 2), stringsAsFactors = FALSE)
}
if (anyNA(meta$cohort)) stop("Unrecognised sample prefix: ",
                             paste(unique(samples[is.na(meta$cohort)]), collapse = ", "))
cat(">>> cohorts in the file:", paste(names(table(meta$cohort)), table(meta$cohort), collapse = ", "), "\n")
cat(">>> infection groups:", paste(names(table(meta$infection)), table(meta$infection), collapse = ", "), "\n")
if (ASSAY == "RNA") cat(">>> (RNA: last name field 1 = Mtb-challenged, 0 = non-challenged)\n")
if (!all(KEEP_COHORTS %in% meta$cohort)) stop("Cohort(s) not found: ",
    paste(setdiff(KEEP_COHORTS, meta$cohort), collapse = ", "), " -- check the sample-name patterns above.")

dropped <- setdiff(unique(meta$cohort), KEEP_COHORTS)
if (length(dropped)) cat(">>> leaving out:", paste(dropped, collapse = ", "),
                         "(", sum(!meta$cohort %in% KEEP_COHORTS), "samples )\n")
meta <- meta[meta$cohort %in% KEEP_COHORTS, ]

## a paired design needs both conditions from the same donor
paired_donors <- names(which(tapply(meta$infection, meta$donor, function(x) length(unique(x))) == 2))
if (length(paired_donors) < length(unique(meta$donor)))
  cat(">>> dropping", sum(!meta$donor %in% paired_donors), "samples whose donor has no INF/NEG pair\n")
meta <- meta[meta$donor %in% paired_donors, ]
meta$cohort <- factor(meta$cohort, levels = KEEP_COHORTS)
meta$infection <- factor(meta$infection, levels = c(CTRL_LEVEL, CASE_LEVEL))

## donor index WITHIN cohort: the recipe for individuals nested in groups (DESeq2 vignette),
## which lets one model give both the per-cohort effects and the difference between them
meta <- meta[order(meta$cohort, meta$donor), ]
meta$donor.n <- factor(unlist(lapply(split(meta$donor, meta$cohort),
                                     function(d) match(d, unique(d)))))
rownames(meta) <- meta$sample
counts <- counts[, meta$sample, drop = FALSE]
cat(">>> analysed:", nrow(meta), "samples |",
    paste(names(table(meta$cohort)), table(meta$cohort), collapse = ", "), "| donors:",
    paste(names(table(meta$cohort[!duplicated(meta$donor)])), table(meta$cohort[!duplicated(meta$donor)]), collapse = ", "), "\n")
write.csv(meta, file.path(dir_out, paste0("GSE165709_", ASSAY, "_sample_metadata.csv")), row.names = FALSE)

## ---- 2. model ----
design <- ~ cohort + cohort:donor.n + cohort:infection
mm <- model.matrix(design, meta)
mm <- mm[, colSums(is.na(mm)) == 0 & apply(mm, 2, function(x) !all(x == 0)), drop = FALSE]
if (qr(mm)$rank < ncol(mm)) {            # unbalanced donor numbers leave empty columns
  keep <- qr(mm)$pivot[seq_len(qr(mm)$rank)]
  cat(">>> dropping", ncol(mm) - length(keep), "redundant design columns:",
      paste(setdiff(colnames(mm), colnames(mm)[keep]), collapse = ", "), "\n")
  mm <- mm[, keep, drop = FALSE]
}
## the RNA matrix holds fractional estimated counts, so it is rounded for DESeq2
dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(counts)), meta, design = ~ 1)
dds <- dds[rowSums(DESeq2::counts(dds)) >= MIN_COUNT, ]
cat(">>>", nrow(dds), "features kept after filtering\n")
dds <- DESeq2::DESeq(dds, full = mm, betaPrior = FALSE)
cat(">>> coefficients:", paste(DESeq2::resultsNames(dds), collapse = ", "), "\n")

coef_of <- function(co) grep(paste0("cohort", co, ".*infection"), DESeq2::resultsNames(dds), value = TRUE)[1]
save_res <- function(r, tag) {
  d <- as.data.frame(r); d$feature <- rownames(d)
  ## RNA row names are "ENSG00000000003_TSPAN6": split them so the symbol can be matched later
  if (ASSAY == "RNA" && all(grepl("^ENSG[0-9]+_", utils::head(d$feature, 20)))) {
    d$ensembl <- sub("_.*$", "", d$feature); d$SYMBOL <- sub("^[^_]*_", "", d$feature)
  }
  d <- d[order(d$padj), c("feature", setdiff(names(d), "feature"))]
  write.csv(d, file.path(dir_out, paste0("GSE165709_", ASSAY, "_", tag, ".csv")), row.names = FALSE)
  cat("   ", tag, ": significant (padj <", PADJ_CUT, "):", sum(d$padj < PADJ_CUT, na.rm = TRUE), "\n")
  d
}
cat(">>> infection effect per cohort\n")
res <- list()
for (co in KEEP_COHORTS) {
  cf <- coef_of(co)
  if (is.na(cf)) { cat("    !! no infection coefficient for", co, "\n"); next }
  res[[co]] <- save_res(DESeq2::results(dds, name = cf), paste0("INF_vs_NEG_in_", co))
}
## difference between the two responses: responds in HC but not (or less) in PLWH
if (length(KEEP_COHORTS) == 2 && all(!is.na(vapply(KEEP_COHORTS, coef_of, character(1))))) {
  cat(">>> difference between the two responses\n")
  res$diff <- save_res(DESeq2::results(dds, contrast = list(coef_of(KEEP_COHORTS[2]), coef_of(KEEP_COHORTS[1]))),
                       paste0("response_", KEEP_COHORTS[2], "_vs_", KEEP_COHORTS[1]))
}

## ---- 3. annotate ATAC peaks to genes ----
if (ASSAY == "ATAC") {
  cat(">>> annotating peaks\n")
  for (tag in names(res)) {
    pk <- res[[tag]]$feature
    gr <- GenomicRanges::GRanges(seqnames = sub(":.*", "", pk),
          ranges = IRanges::IRanges(start = as.numeric(sub(".*:(.*)-.*", "\\1", pk)),
                                    end   = as.numeric(sub(".*-(.*)", "\\1", pk))))
    names(gr) <- pk
    an <- as.data.frame(ChIPseeker::annotatePeak(gr,
            TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
            tssRegion = TSS_REGION, annoDb = "org.Hs.eg.db", verbose = FALSE))
    if (nrow(an) != length(gr)) stop("annotatePeak returned ", nrow(an), " rows for ", length(gr), " peaks.")
    an$feature <- names(gr)                       # the id travels with the peak; never match on coordinates
    out <- merge(res[[tag]], an[, intersect(c("feature", "annotation", "distanceToTSS", "SYMBOL"), names(an))],
                 by = "feature", all.x = TRUE)
    out <- out[order(out$padj), ]
    nm <- if (tag == "diff") paste0("response_", KEEP_COHORTS[2], "_vs_", KEEP_COHORTS[1]) else paste0("INF_vs_NEG_in_", tag)
    write.csv(out, file.path(dir_out, paste0("GSE165709_ATAC_", nm, "_annotated.csv")), row.names = FALSE)
    cat("    ", nm, ":", sum(!is.na(out$SYMBOL)), "of", nrow(out), "peaks annotated |",
        length(unique(out$SYMBOL[!is.na(out$padj) & out$padj < PADJ_CUT & !is.na(out$SYMBOL)])), "genes with a significant peak\n")
  }
}

cat("\n>>> DONE. Outputs in", dir_out, "\n")
cat(">>> Run this script once with ASSAY <- \"RNA\" and once with ASSAY <- \"ATAC\", then run GSE165709_integration_step17.R.\n")
