## =============================================================
## Input:  <ACC>_TBmodules_MM_GS.csv  (from WGCNA_TBmodules_step1.R —
##          already just the Selected top-N modules, e.g. green/purple/blue)
## Output, per module: one CSV each for GO Biological Process, KEGG,
##          Reactome, and Hallmark (MSigDB H), plus one
##          <ACC>_TopTerms_Summary.csv to actually pick the name from.
##
## =============================================================

# One-time setup (Bioconductor + CRAN):
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))
#
# msigdbr changed its data delivery in v10+: it now needs the separate
# 'msigdbdf' package (gene-set data moved out of CRAN for size reasons).
# Install BOTH from the maintainer's r-universe repo:
#install.packages(c("msigdbr", "msigdbdf"),
#                  repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))
# If that's out of date by the time you run this, check https://igordot.github.io/msigdbr/


required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "dplyr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing/broken package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Try: BiocManager::install(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "), "))",
       "\n  Then confirm each loads cleanly with library(<name>) BEFORE rerunning this script.")
}

have_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
if (!have_msigdbr) {
  cat("!! msigdbr not installed -> Hallmark enrichment will be SKIPPED.\n")
  cat("   Install with: install.packages(c('msigdbr','msigdbdf'), repos = c('https://igordot.r-universe.dev','https://cloud.r-project.org'))\n")
}

## ---- CONFIG ----
ACC         <- "GSE114192"
PADJ_CUTOFF <- 0.05   # for the "top term" summary only; full tables are saved unfiltered

## ---- locate this script's own folder  ----
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
SCRIPT_DIR  <- get_script_dir()
cat("Script folder detected as:", SCRIPT_DIR, "\n")
BASE_DIR    <- file.path(SCRIPT_DIR, "..", ACC)
dir_results <- file.path(BASE_DIR, "10_wgcna_results")
dir_enrich  <- file.path(BASE_DIR, "11_module_enrichment")
dir.create(dir_enrich, recursive = TRUE, showWarnings = FALSE)

## ---- load the Selected module gene lists from step 1 ----
mm_gs_file <- file.path(dir_results, paste0(ACC, "_TBmodules_MM_GS.csv"))
if (!file.exists(mm_gs_file)) stop("Run WGCNA_TBmodules_step1.R first — missing: ", mm_gs_file)
genes_tab <- read.csv(mm_gs_file, stringsAsFactors = FALSE)
SELECTED_MODULES <- unique(genes_tab$Module)
cat(">>> Modules to annotate:", paste(SELECTED_MODULES, collapse = ", "), "\n")

## background = every gene that went into the WGCNA network (fair test),
## NOT the whole genome
bg_file <- file.path(dir_results, paste0(ACC, "_Gene_Module_Assignment.csv"))
if (!file.exists(bg_file)) stop("Missing: ", bg_file)
bg_ensembl <- read.csv(bg_file, stringsAsFactors = FALSE)$Gene

## ---- ID mapping: ENSEMBL -> ENTREZ + SYMBOL (KEGG/Reactome/Hallmark need Entrez) ----
cat(">>> Mapping", length(bg_ensembl), "background genes ENSEMBL -> ENTREZ/SYMBOL...\n")
id_map <- suppressMessages(
  clusterProfiler::bitr(bg_ensembl, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"),
                         OrgDb = org.Hs.eg.db::org.Hs.eg.db)
)
id_map <- id_map[!duplicated(id_map$ENSEMBL), ]
cat("    mapped", nrow(id_map), "/", length(bg_ensembl),
    "(", round(100 * nrow(id_map) / length(bg_ensembl), 1), "%)\n")
bg_entrez <- unique(na.omit(id_map$ENTREZID))

hallmark_sets <- NULL
if (have_msigdbr) {
  hallmark_sets <- dplyr::distinct(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    gs_name, entrez_gene
  )
  cat(">>> Loaded", length(unique(hallmark_sets$gs_name)), "Hallmark gene sets.\n")
}

safe_enrich <- function(fn, label) {
  tryCatch(as.data.frame(fn()), error = function(e) {
    cat("    !!", label, "failed:", conditionMessage(e), "\n")
    data.frame()
  })
}

top1 <- function(df, mod, label) {
  if (nrow(df) == 0 || !("p.adjust" %in% names(df)) || min(df$p.adjust, na.rm = TRUE) > PADJ_CUTOFF) {
    return(data.frame(Module = mod, Source = label, Term = NA, p.adjust = NA, Count = NA))
  }
  d <- df[df$p.adjust < PADJ_CUTOFF, ]
  d <- d[order(d$p.adjust), ]
  data.frame(Module = mod, Source = label, Term = d$Description[1],
             p.adjust = d$p.adjust[1], Count = d$Count[1])
}

run_one_module <- function(mod) {
  cat("\n----", mod, "----\n")
  ens    <- genes_tab$Gene[genes_tab$Module == mod]
  entrez <- unique(na.omit(id_map$ENTREZID[id_map$ENSEMBL %in% ens]))
  cat(length(ens), "genes ->", length(entrez), "with an Entrez ID\n")
  if (length(entrez) < 3) {
    cat("    too few mapped genes, skipping this module.\n")
    return(NULL)
  }

  go_df <- safe_enrich(function() clusterProfiler::enrichGO(gene = entrez, universe = bg_entrez,
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                       ont = "BP", pAdjustMethod = "BH",
                                       pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "GO_BP")
  kegg_df <- safe_enrich(function() clusterProfiler::enrichKEGG(gene = entrez, universe = bg_entrez, organism = "hsa",
                                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), "KEGG")
  reac_df <- safe_enrich(function() ReactomePA::enrichPathway(gene = entrez, universe = bg_entrez, organism = "human",
                                              pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                              readable = TRUE), "Reactome")
  hm_df <- data.frame()
  if (have_msigdbr) {
    hm_df <- safe_enrich(function() clusterProfiler::enricher(gene = entrez, universe = bg_entrez, TERM2GENE = hallmark_sets,
                                         pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), "Hallmark")
  }

  write.csv(go_df,   file.path(dir_enrich, paste0(ACC, "_", mod, "_GO_BP.csv")),    row.names = FALSE)
  write.csv(kegg_df, file.path(dir_enrich, paste0(ACC, "_", mod, "_KEGG.csv")),     row.names = FALSE)
  write.csv(reac_df, file.path(dir_enrich, paste0(ACC, "_", mod, "_Reactome.csv")), row.names = FALSE)
  write.csv(hm_df,   file.path(dir_enrich, paste0(ACC, "_", mod, "_Hallmark.csv")), row.names = FALSE)

  rbind(top1(go_df, mod, "GO_BP"), top1(kegg_df, mod, "KEGG"),
        top1(reac_df, mod, "Reactome"), top1(hm_df, mod, "Hallmark"))
}

top_summary <- do.call(rbind, lapply(SELECTED_MODULES, run_one_module))
write.csv(top_summary, file.path(dir_enrich, paste0(ACC, "_TopTerms_Summary.csv")), row.names = FALSE)

cat("\n>>> DONE. Per-module GO_BP/KEGG/Reactome/Hallmark tables +",
    paste0(ACC, "_TopTerms_Summary.csv"), "in", dir_enrich, "\n")
