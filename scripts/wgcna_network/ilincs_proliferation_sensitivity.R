## Sensitivity analysis: does the GSE161829 drug ranking disagree with the pilot because of proliferation genes?
##
## About half of the GSE161829 TB_UP signature is cell-division genes (CDK1, TOP2A, CCNB1, ...), which the
## whole-blood pilot signature does not contain. Compounds that block cell division reverse such a signature
## regardless of TB. This script removes the proliferation genes, queries iLINCS again with what is left,
## and compares three UP-only rankings with the pilot's UP-only ranking:
##   full GSE161829 signature   vs pilot
##   proliferation-free version vs pilot
## Two comparisons per pair: (1) per LINCS signature, do the two queries score it in the same direction
## (Spearman of similarity over shared signatures); (2) per named compound, overlap of the top-N reversers
## against a permutation null. If removing proliferation genes makes the agreement clearly better, the
## earlier non-replication of drugs is explained by signature composition rather than by noise.
## Run from scripts/wgcna_network (needs internet once; the new query is cached).
## Outputs in GSE161829/TBneg_vs_ATB/13_ilincs_signature/sensitivity_no_proliferation/.

## ---- CONFIG ----
ILINCS_LIBRARY <- "CP"
TOP_N   <- c(50, 100, 200)
N_PERM  <- 2000
SEED    <- 42
## proliferation genes: Hallmark G2M checkpoint, E2F targets and mitotic spindle (msigdbr), plus CDKN1B
PROLIF_SETS <- c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_MITOTIC_SPINDLE")

for (p in c("drugfindR", "msigdbr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
sig_dir   <- file.path(ROOT, "GSE161829", "TBneg_vs_ATB", "13_ilincs_signature")
pilot_dir <- file.path(ROOT, "GSE114192", "13_ilincs_signature")
out_dir   <- file.path(sig_dir, "sensitivity_no_proliferation"); dir.create(out_dir, showWarnings = FALSE)

## ---- 1. proliferation-free signature ----
up <- read.csv(file.path(sig_dir, "GSE161829_TB_UP_signature.csv"), stringsAsFactors = FALSE)
hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
prolif <- unique(c(hs$gene_symbol[hs$gs_name %in% PROLIF_SETS], "CDKN1B"))
is_prolif <- up$Name_GeneSymbol %in% prolif
cat(">>> TB_UP:", nrow(up), "genes |", sum(is_prolif), "proliferation genes removed:",
    paste(up$Name_GeneSymbol[is_prolif], collapse = ", "), "\n")
cat("    kept:", paste(up$Name_GeneSymbol[!is_prolif], collapse = ", "), "\n")
up_clean <- up[!is_prolif, ]
write.csv(up_clean, file.path(out_dir, "GSE161829_TB_UP_signature_no_proliferation.csv"), row.names = FALSE)
if (nrow(up_clean) < 5) stop("Fewer than 5 genes left after removing proliferation genes; the query would be meaningless.")

## ---- 2. iLINCS query (cached) ----
raw_clean <- file.path(out_dir, "GSE161829_concordants_UP_no_proliferation_raw.csv")
if (!file.exists(raw_clean)) {
  cat(">>> Querying iLINCS with", nrow(up_clean), "genes (needs internet)...\n")
  res <- drugfindR::getConcordants(up_clean, ilincsLibrary = ILINCS_LIBRARY)
  write.csv(res, raw_clean, row.names = FALSE)
}
rd <- function(f) { if (!file.exists(f)) stop("Missing: ", f); x <- read.csv(f, stringsAsFactors = FALSE); names(x) <- tolower(names(x)); x }
pilot <- rd(file.path(pilot_dir, "GSE114192_concordants_UP_raw.csv"))
full  <- rd(file.path(sig_dir,   "GSE161829_concordants_UP_raw.csv"))
clean <- rd(raw_clean)

## ---- 3. comparisons ----
per_signature <- function(a, b) {        # same LINCS signature scored by both queries
  m <- merge(a[, c("signatureid", "similarity")], b[, c("signatureid", "similarity")], by = "signatureid")
  c(n = nrow(m), rho = suppressWarnings(stats::cor(m$similarity.x, m$similarity.y, method = "spearman")))
}
per_drug <- function(x) {                # best (most negative) score per named compound
  x <- x[!grepl("^BRD-", x$treatment), ]
  s <- tapply(x$similarity, toupper(x$treatment), min); sort(s)
}
top_overlap <- function(a, b, n) {
  da <- per_drug(a); db <- per_drug(b); univ <- intersect(names(da), names(db))
  ta <- utils::head(names(da)[names(da) %in% univ], n); tb <- utils::head(names(db)[names(db) %in% univ], n)
  obs <- length(intersect(ta, tb))
  set.seed(SEED); null <- replicate(N_PERM, length(intersect(sample(univ, n), sample(univ, n))))
  c(shared = obs, expected = round(mean(null), 1), p = signif((sum(null >= obs) + 1) / (N_PERM + 1), 3),
    rho_all = round(suppressWarnings(stats::cor(rank(da[univ]), rank(db[univ]), method = "spearman")), 3))
}
rows <- list()
for (cmp in list(list("Full GSE161829 signature", full), list("Proliferation genes removed", clean))) {
  ps <- per_signature(pilot, cmp[[2]])
  for (n in TOP_N) {
    o <- top_overlap(pilot, cmp[[2]], n)
    rows[[length(rows) + 1]] <- data.frame(GSE161829_signature = cmp[[1]], Top_N = n,
      Shared_with_pilot = o[["shared"]], Expected_by_chance = o[["expected"]], P = o[["p"]],
      Drug_rank_rho = o[["rho_all"]], Signatures_compared = ps[["n"]],
      Per_signature_rho = round(ps[["rho"]], 3))
  }
}
R <- do.call(rbind, rows)
write.csv(R, file.path(out_dir, "GSE161829_vs_pilot_drug_agreement.csv"), row.names = FALSE)
cat("\n========== agreement with the pilot's drug ranking (UP-only queries) ==========\n"); print(R, row.names = FALSE)
cat("\nPer_signature_rho > 0 means the two disease signatures push the same LINCS signatures the same way.\n",
    "If removing proliferation genes raises it and the top-N overlap, the earlier disagreement was driven by\n",
    "the proliferation component of the GSE161829 signature.\n", sep = "")
cat(">>> DONE. Outputs in", out_dir, "\n")
