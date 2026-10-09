## STEP 18: is there an ATAC signal at the level of a GENE SET, even though almost no single peak
## survives FDR? Per-peak testing penalises every one of the ~116,000 peaks separately; asking whether
## the peaks of a gene set are shifted as a group pools that evidence instead of splitting it.
##
## Three tests per gene set (High priority, DEG-overlap, each module), per cohort:
##   1. Wilcoxon on the per-peak Wald statistic: are the set's peaks shifted up or down as a group?
##   2. Wilcoxon on |statistic|: are they moved more than the rest, in either direction?
##   3. Sample-label permutation of the set membership (keeps the peak-per-gene structure), which is
##      the honest null here because peaks of the same gene are correlated and Wilcoxon assumes they
##      are not. The permutation p-value is the one to quote; Wilcoxon is shown for comparison.
## Outputs in external_results/GSE165709/atac_geneset/.  Base R only.

## ---- CONFIG ----
COHORTS    <- c("HC", "PLWH")
N_PERM     <- 2000
MIN_GENES  <- 10        # a set needs at least this many genes with a peak
TSS_ONLY   <- FALSE     # TRUE = use only promoter peaks (cleaner link, far fewer peaks)

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
EXT     <- file.path(ROOT_DIR, "external_results", "GSE165709")
dir_out <- file.path(EXT, "atac_geneset"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) { if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE) }

## ---- gene sets from the previous analysis ----
gp <- rd(file.path(run_out_dir(file.path(ROOT_DIR, ACC)), "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
cat(">>> previous gene table:", nrow(gp), "genes |", ACC, "\n")
SETS <- list("High priority" = gp$Gene[gp$Final_Priority == "High"],
             "High or Medium" = gp$Gene[gp$Final_Priority %in% c("High", "Medium")],
             "DEG-overlap"    = gp$Gene[gp$DEG_WGCNA_overlap %in% TRUE])
for (mo in sort(unique(gp$Module))) SETS[[paste("module", mo)]] <- gp$Gene[gp$Module == mo]
SETS <- lapply(SETS, function(x) unique(x[!is.na(x) & nzchar(x)]))

## ---- test one cohort ----
run_cohort <- function(co) {
  f <- file.path(EXT, "ATAC", paste0("GSE165709_ATAC_INF_vs_NEG_in_", co, "_annotated.csv"))
  a <- rd(f)
  a <- a[!is.na(a$SYMBOL) & !is.na(a$stat), ]
  if (TSS_ONLY) a <- a[grepl("^Promoter", a$annotation), ]
  cat("\n>>>", co, ":", nrow(a), "annotated peaks with a test statistic |",
      length(unique(a$SYMBOL)), "genes | significant peaks (padj<0.05):", sum(a$padj < 0.05, na.rm = TRUE), "\n")

  genes_all <- unique(a$SYMBOL)
  by_gene <- split(seq_len(nrow(a)), a$SYMBOL)          # peaks belonging to each gene
  out <- list()
  for (nm in names(SETS)) {
    inset <- intersect(SETS[[nm]], genes_all)
    if (length(inset) < MIN_GENES) { cat("   ", nm, ": only", length(inset), "genes with a peak, skipped\n"); next }
    idx  <- unlist(by_gene[inset], use.names = FALSE)
    rest <- setdiff(seq_len(nrow(a)), idx)
    s_in <- a$stat[idx]; s_out <- a$stat[rest]

    w_signed <- suppressWarnings(stats::wilcox.test(s_in, s_out))$p.value
    w_abs    <- suppressWarnings(stats::wilcox.test(abs(s_in), abs(s_out)))$p.value
    obs      <- mean(abs(s_in))

    ## permutation: resample whole GENES, so the number of peaks per gene stays realistic
    n_g <- length(inset)
    set.seed(42)
    null <- vapply(seq_len(N_PERM), function(i) {
      g <- sample(genes_all, n_g)
      mean(abs(a$stat[unlist(by_gene[g], use.names = FALSE)]))
    }, numeric(1))
    p_perm <- (sum(null >= obs) + 1) / (N_PERM + 1)

    out[[length(out) + 1]] <- data.frame(
      Cohort = co, Set = nm, N_genes_in_set = length(SETS[[nm]]), N_genes_with_peak = n_g, N_peaks = length(idx),
      Mean_abs_stat = round(obs, 3), Mean_abs_stat_rest = round(mean(abs(s_out)), 3),
      Null_mean = round(mean(null), 3), P_permutation = signif(p_perm, 3),
      P_wilcoxon_signed = signif(w_signed, 3), P_wilcoxon_abs = signif(w_abs, 3),
      Mean_log2FC = round(mean(a$log2FoldChange[idx], na.rm = TRUE), 3), stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

res <- do.call(rbind, lapply(COHORTS, run_cohort))
if (is.null(res)) stop("No gene set had enough genes with an annotated peak.")
res$FDR_permutation <- signif(stats::p.adjust(res$P_permutation, "BH"), 3)
res <- res[order(res$Cohort, res$P_permutation), ]
write.csv(res, file.path(dir_out, paste0("GSE165709_ATAC_geneset_tests", if (TSS_ONLY) "_promoterOnly" else "", ".csv")), row.names = FALSE)

cat("\n========== gene-set level ATAC results ==========\n")
print(res[, c("Cohort", "Set", "N_genes_with_peak", "N_peaks", "Mean_abs_stat", "Null_mean",
              "P_permutation", "FDR_permutation", "P_wilcoxon_abs")], row.names = FALSE)
cat("\nMean_abs_stat is the average |Wald statistic| of the set's peaks; Null_mean is what random gene\n")
cat("sets of the same size give. P_permutation is the one to quote -- the Wilcoxon columns treat peaks\n")
cat("of the same gene as independent, which they are not, so they are optimistic.\n")
if (all(res$FDR_permutation >= 0.05))
  cat("\n>>> No gene set is shifted beyond chance. The ATAC layer carries no usable signal here, at peak\n    level or set level, and that is the honest result to report.\n")
cat("\n>>> DONE. Output in", dir_out, "\n")
