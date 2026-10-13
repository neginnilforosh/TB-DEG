## GSEA on the coefficients of the interaction model  Expression ~ TB + DM + TB:DM  (GSE114192).
##
## Why: comparing two separate GSEA runs (TB-DM vs TB: interferon DOWN; DM vs healthy: interferon UP)
## suggests that diabetes changes the interferon response to TB, but two separate tests are not a test
## of the difference. The TB:DM coefficient is that difference -- (TB-DM minus DM) minus (TB minus
## healthy) -- so a gene set enriched on it is the formal evidence. A negative NES for an interferon set
## means TB raises that set LESS in people with diabetes than in people without.
##
## Uses the same ranking (sign(log2FC) * -log10 p), gene-set collections, size limits and seed as
## WGCNA_GSEA_step11.R, so the numbers are comparable. Run from scripts/wgcna_network.
## Outputs in GSE114192/TB_DM/17_gsea_interaction/.

required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "msigdbr")
ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing: ", paste(required_pkgs[!ok], collapse = ", "), " (same packages as step 11)")

## ---- CONFIG ----
ACC        <- "GSE114192"
RUN        <- "TB_DM"
COEFS      <- c(interaction = "interaction", main_TB = "main_TB", main_DM = "main_DM")  # file suffixes
MIN_GS <- 15; MAX_GS <- 500; PADJ_SIG <- 0.05; SEED <- 42
COLLECTIONS <- list(HALLMARK = list(collection = "H"),
                    REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"))

## ---- paths ----
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, ACC)) && dirname(d) != d) d <- dirname(d); d })
in_dir  <- file.path(ROOT, ACC, RUN, "06_deg_results")
out_dir <- file.path(ROOT, ACC, RUN, "17_gsea_interaction"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- gene sets (once) ----
term2gene <- do.call(rbind, lapply(names(COLLECTIONS), function(nm) {
  x <- do.call(msigdbr::msigdbr, c(list(species = "Homo sapiens"), COLLECTIONS[[nm]]))
  data.frame(term = x$gs_name, gene = x$gene_symbol, Collection = nm, stringsAsFactors = FALSE)
}))
coll_of <- unique(term2gene[, c("term", "Collection")])

## same ranking as step 11
build_ranking <- function(deg) {
  deg <- deg[!is.na(deg$pvalue) & !is.na(deg$log2FoldChange), ]
  deg$stat <- sign(deg$log2FoldChange) * -log10(pmax(deg$pvalue, 1e-300))
  map <- suppressMessages(clusterProfiler::bitr(deg$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  d <- merge(deg[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL")
  d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$SYMBOL), ]
  sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE)
}

all_ifn <- list()
for (nm in names(COEFS)) {
  f <- file.path(in_dir, paste0(ACC, "_interaction_", COEFS[[nm]], ".csv"))
  if (!file.exists(f)) { cat("!! missing:", f, "\n"); next }
  r <- build_ranking(read.csv(f, stringsAsFactors = FALSE))
  cat("\n==========", nm, "| ranked genes:", length(r), "==========\n")
  set.seed(SEED)
  g <- clusterProfiler::GSEA(r, TERM2GENE = term2gene[, c("term", "gene")], minGSSize = MIN_GS, maxGSSize = MAX_GS,
                             pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
  res <- as.data.frame(g)
  if (!nrow(res)) { cat("  no gene set tested\n"); next }
  res <- merge(res, coll_of, by.x = "ID", by.y = "term", all.x = TRUE)
  res <- res[order(res$p.adjust), c("Collection", "ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")]
  write.csv(res, file.path(out_dir, paste0(ACC, "_GSEA_", nm, ".csv")), row.names = FALSE)
  cat("  significant (padj <", PADJ_SIG, "):", sum(res$p.adjust < PADJ_SIG), "\n")
  print(utils::head(res[res$p.adjust < PADJ_SIG, c("Collection", "ID", "NES", "p.adjust")], 8), row.names = FALSE)
  ifn <- res[grepl("INTERFERON", res$ID), c("ID", "NES", "pvalue", "p.adjust")]
  if (nrow(ifn)) all_ifn[[nm]] <- cbind(Coefficient = nm, ifn)
}

## ---- the question this script exists for ----
if (length(all_ifn)) {
  ifn <- do.call(rbind, all_ifn)
  write.csv(ifn, file.path(out_dir, paste0(ACC, "_GSEA_interferon_by_coefficient.csv")), row.names = FALSE)
  cat("\n========== interferon gene sets on each coefficient ==========\n"); print(ifn, row.names = FALSE)
  x <- ifn[ifn$Coefficient == "interaction", ]
  if (nrow(x)) {
    sig_dn <- x$p.adjust < PADJ_SIG & x$NES < 0
    cat("\n>>> Interaction term:",
        if (any(sig_dn)) paste0("interferon set(s) significantly NEGATIVE (", paste(x$ID[sig_dn], collapse = ", "),
                                ") -- formal support that diabetes blunts the interferon response to TB.")
        else "no interferon set reaches padj <0.05 with a negative NES -- the blunting seen in the two separate\n    comparisons is NOT confirmed by the formal test; report it as a pattern only.", "\n")
  }
}
cat("\n>>> DONE. Outputs in", out_dir, "\n")
