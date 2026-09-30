## =============================================================
## STEP 11: GSEA on the full ranked DEG list, per comparison
## Why GSEA (and not another over-representation tool) for the blue module: blue is
## ~1,400 genes that mostly shift a little and together (961 lower / 416 higher in TB),
## which a hard cutoff + ORA blurs. GSEA uses every gene, needs no cutoff, and reports a
## direction (NES > 0 = higher in TB_Only / the case group, NES < 0 = lower).
##
## What it does, per comparison found in <dataset>/06_deg_results/*_DEG.csv (DESeq2-style):
##   1. ranking statistic = sign(log2FC) * -log10(pvalue)        (all genes, no filtering)
##   2. clusterProfiler::GSEA against MSigDB sets: Hallmark, Reactome, GO BP and the
##      C8 cell-type signature sets (the direct test of the "blue = lymphocyte signal" idea)
##   3. saves every result (unfiltered) + a significant-only table
##   4. for GSE114192 only: a term x WGCNA-module table -- how many of each significant
##      pathway's leading-edge genes sit in green / purple / blue
## The miRNA (GSE229020) tables are skipped: gene sets are for genes, not miRNAs.
## =============================================================

required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "msigdbr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  msigdbr needs its data package too:  install.packages(c(\"msigdbr\", \"msigdbdf\"),",
       " repos = c(\"https://igordot.r-universe.dev\", \"https://cloud.r-project.org\"))")
}

## ---- CONFIG ----
MODULE_DATASET <- "GSE114192"   # the dataset whose WGCNA modules get overlaid on the GSEA results
MIN_GS <- 15; MAX_GS <- 500     # gene-set size limits (genes present in the ranked list)
PADJ_SIG <- 0.05
COLLECTIONS <- list(            # name -> arguments for msigdbr()
  HALLMARK = list(collection = "H"),
  REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"),
  GOBP     = list(collection = "C5", subcollection = "GO:BP"),
  CELLTYPE = list(collection = "C8")
)

## ---- locate this script's folder and the repo root (works from scripts/ or scripts/wgcna_network/) ----
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()
ROOT_DIR <- local({ d <- SCRIPT_DIR
  while (length(list.files(d, pattern = "^GSE[0-9]+$")) == 0 && dirname(d) != d) d <- dirname(d); d })
cat(">>> Repo root:", ROOT_DIR, "\n")

## =============================================================
## helpers (pure R)
## =============================================================
# DEG table (ENSEMBL ids) + mapping table -> named, decreasing ranking vector keyed by SYMBOL
build_ranking <- function(deg, id_map) {
  deg <- deg[!is.na(deg$pvalue) & !is.na(deg$log2FoldChange), ]
  deg$stat <- sign(deg$log2FoldChange) * -log10(pmax(deg$pvalue, 1e-300))
  d <- merge(deg[, c("gene", "stat")], id_map, by.x = "gene", by.y = "ENSEMBL")
  d <- d[!is.na(d$SYMBOL) & nzchar(d$SYMBOL), ]
  d <- d[order(-abs(d$stat)), ]
  d <- d[!duplicated(d$SYMBOL), ]                     # one gene per symbol: keep its strongest signal
  r <- setNames(d$stat, d$SYMBOL)
  sort(r, decreasing = TRUE)
}

# significant GSEA terms -> per-module leading-edge counts (module_map: data.frame SYMBOL, Module)
leading_edge_by_module <- function(gsea_df, module_map) {
  out <- list()
  mod_size <- table(module_map$Module)
  for (i in seq_len(nrow(gsea_df))) {
    le <- strsplit(as.character(gsea_df$core_enrichment[i]), "/", fixed = TRUE)[[1]]
    if (length(le) == 0) next
    mods <- module_map$Module[match(le, module_map$SYMBOL)]
    tab  <- table(factor(mods, levels = names(mod_size)))
    for (m in names(mod_size)) {
      if (tab[[m]] == 0) next
      out[[length(out) + 1]] <- data.frame(
        Collection = gsea_df$Collection[i], Term = gsea_df$Description[i], NES = gsea_df$NES[i],
        p.adjust = gsea_df$p.adjust[i], n_leading_edge = length(le), Module = m,
        n_in_module = as.integer(tab[[m]]), pct_of_leading_edge = round(100 * tab[[m]] / length(le), 1),
        pct_of_module = round(100 * tab[[m]] / mod_size[[m]], 2), stringsAsFactors = FALSE)
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

## =============================================================
## gene sets (once)
## =============================================================
cat(">>> Loading MSigDB gene sets...\n")
term2gene <- list()
for (nm in names(COLLECTIONS)) {
  t2g <- tryCatch({
    x <- do.call(msigdbr::msigdbr, c(list(species = "Homo sapiens"), COLLECTIONS[[nm]]))
    unique(data.frame(term = x$gs_name, gene = x$gene_symbol, stringsAsFactors = FALSE))
  }, error = function(e) { cat("  !!", nm, "failed to load:", conditionMessage(e), "\n"); NULL })
  if (!is.null(t2g) && nrow(t2g) > 0) { term2gene[[nm]] <- t2g; cat("  ", nm, ":", length(unique(t2g$term)), "sets\n") }
}
if (length(term2gene) == 0) stop("No gene-set collection could be loaded -- check the msigdbr/msigdbdf install (see PACKAGES.md).")

## =============================================================
## per-comparison GSEA
## =============================================================
deg_files <- list.files(ROOT_DIR, pattern = "_DEG\\.csv$", recursive = TRUE, full.names = TRUE)
deg_files <- deg_files[grepl("/06_deg_results/", deg_files) & !grepl("sensitivity", deg_files)]
cat(">>> Found", length(deg_files), "DEG tables.\n")

for (f in deg_files) {
  hdr <- names(read.csv(f, nrows = 1))
  if (!all(c("gene", "log2FoldChange", "pvalue") %in% hdr)) { cat("  skipping (not a DESeq2-style table):", basename(f), "\n"); next }
  ds   <- sub("^(GSE[0-9]+)_.*$", "\\1", basename(f))
  comp <- sub("_DEG\\.csv$", "", basename(f))
  out_dir <- file.path(ROOT_DIR, ds, "17_gsea"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\n==========", comp, "==========\n")

  deg <- read.csv(f, stringsAsFactors = FALSE)
  id_map <- suppressMessages(clusterProfiler::bitr(unique(deg$gene), fromType = "ENSEMBL", toType = "SYMBOL",
                                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  ranking <- build_ranking(deg, id_map)
  cat("  ranked genes:", length(ranking), "(of", nrow(deg), "in the table)\n")

  all_res <- list()
  for (nm in names(term2gene)) {
    res <- tryCatch(clusterProfiler::GSEA(geneList = ranking, TERM2GENE = term2gene[[nm]],
                                          minGSSize = MIN_GS, maxGSSize = MAX_GS, pvalueCutoff = 1,
                                          pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE),
                    error = function(e) { cat("  !!", nm, "GSEA failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) { cat("  ", nm, ": no gene sets passed the size filter\n"); next }
    df$Collection <- nm
    write.csv(df, file.path(out_dir, paste0(comp, "_GSEA_", nm, ".csv")), row.names = FALSE)
    all_res[[nm]] <- df
    cat(sprintf("   %-9s tested %4d sets | significant (padj<%.2f): %3d  (up in case: %d, down: %d)\n", nm, nrow(df), PADJ_SIG,
                sum(df$p.adjust < PADJ_SIG), sum(df$p.adjust < PADJ_SIG & df$NES > 0), sum(df$p.adjust < PADJ_SIG & df$NES < 0)))
  }
  if (length(all_res) == 0) next
  combined <- do.call(rbind, all_res)
  sig <- combined[combined$p.adjust < PADJ_SIG, ]
  sig <- sig[order(sig$p.adjust), ]
  write.csv(sig[, setdiff(names(sig), c("leading_edge"))], file.path(out_dir, paste0(comp, "_GSEA_significant.csv")), row.names = FALSE)
  cat("  top 5 significant:\n"); print(utils::head(sig[, c("Collection", "Description", "NES", "p.adjust")], 5), row.names = FALSE)

  ## ---- module overlay (only for the dataset that has WGCNA modules) ----
  mm_file <- file.path(ROOT_DIR, ds, "10_wgcna_results", paste0(ds, "_TBmodules_MM_GS.csv"))
  if (ds == MODULE_DATASET && file.exists(mm_file) && nrow(sig) > 0) {
    mm  <- read.csv(mm_file, stringsAsFactors = FALSE)
    map <- merge(mm[, c("Gene", "Module")], id_map, by.x = "Gene", by.y = "ENSEMBL")
    map <- unique(map[!is.na(map$SYMBOL), c("SYMBOL", "Module")])
    map <- map[!duplicated(map$SYMBOL), ]
    tm <- leading_edge_by_module(sig, map)
    if (!is.null(tm)) {
      tm <- tm[order(tm$Module, -tm$n_in_module), ]
      write.csv(tm, file.path(out_dir, paste0(comp, "_GSEA_term_by_module.csv")), row.names = FALSE)
      cat("  module overlay saved (", nrow(tm), "term x module rows). Blue's strongest pathways:\n")
      b <- tm[tm$Module == "blue", ]; print(utils::head(b[, c("Collection", "Term", "NES", "n_in_module", "pct_of_leading_edge")], 6), row.names = FALSE)
    }
  }
}

cat("\n>>> DONE. Results in <dataset>/17_gsea/. How to read them:\n")
cat("    NES > 0 = pathway higher in the case group; NES < 0 = lower. The CELLTYPE collection is the\n")
cat("    test of the blue-module lymphocyte idea: T-cell signature sets with NES < 0 whose leading edge\n")
cat("    falls mostly in 'blue' (see *_GSEA_term_by_module.csv) support it; nothing there argues against it.\n")
cat(">>> Note: 6 comparisons x 4 collections; the fgsea permutations make the first run take a few minutes.\n")

