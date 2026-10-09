## Shared/unique significant DEGs across the gene-level datasets (UpSet plot).
## Run after the GSE*_run.R scripts, from scripts/deg_pipeline. GSE229020 (miRNA panel) is excluded:
## miRNA identifiers cannot overlap gene symbols.
##
## The four datasets store gene identifiers in three different formats:
##   GSE222001, GSE114192  plain Ensembl          ENSG00000057704
##   GSE161829             Ensembl.version|symbol ENSG00000138496.16|PARP9
##   GSE99374              gene symbol            HIST1H2AG
## Comparing them as written would report zero overlap for GSE161829 and GSE99374, so everything is
## converted to a gene SYMBOL first (GSE99374 only has symbols, so symbols are the common currency).
## Genes that cannot be mapped are counted and reported rather than silently dropped.

source("00_functions.R")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE))
  stop("Need org.Hs.eg.db and AnnotationDbi: BiocManager::install(c(\"org.Hs.eg.db\", \"AnnotationDbi\"))")

ROOT <- file.path("..", "..")
FILES <- list(
  GSE222001_Healthy_vs_Active = "GSE222001/07_significant_degs/GSE222001_Healthy_vs_Active_sigDEGs.csv",
  GSE222001_Healthy_vs_Latent = "GSE222001/07_significant_degs/GSE222001_Healthy_vs_Latent_sigDEGs.csv",
  GSE161829_TBneg_vs_LTBI     = "GSE161829/07_significant_degs/GSE161829_TBneg_vs_LTBI_sigDEGs.csv",
  GSE161829_TBneg_vs_ATB      = "GSE161829/07_significant_degs/GSE161829_TBneg_vs_ATB_sigDEGs.csv",
  GSE99374_HC_vs_LTBI         = "GSE99374/07_significant_degs/GSE99374_HC_vs_LTBI_sigDEGs.csv",
  GSE114192_HC_vs_TBOnly      = "GSE114192/07_significant_degs/GSE114192_HealthyControl_vs_TBOnly_sigDEGs.csv")

## ---- read and harmonise ----
to_symbol <- function(ids) {
  ids <- trimws(as.character(ids))
  sym <- sub("^.*\\|", "", ids[grepl("\\|", ids)])          # "ENSG...|PARP9" -> "PARP9"
  ens <- sub("\\.[0-9]+$", "", ids[!grepl("\\|", ids) & grepl("^ENSG", ids)])
  oth <- ids[!grepl("\\|", ids) & !grepl("^ENSG", ids)]      # already a symbol
  mapped <- character(0)
  if (length(ens)) {
    m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens),
                                                keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
    mapped <- unname(m[!is.na(m)])
  }
  list(symbols = unique(c(sym, mapped, oth)), n_in = length(unique(ids)),
       n_unmapped = length(unique(ens)) - length(mapped))
}

deg_lists <- list(); report <- list()
for (nm in names(FILES)) {
  f <- file.path(ROOT, FILES[[nm]])
  if (!file.exists(f)) { cat("!! missing, skipped:", FILES[[nm]], "\n"); next }
  g <- read.csv(f, stringsAsFactors = FALSE)$gene
  r <- to_symbol(g)
  report[[nm]] <- data.frame(Comparison = nm, N_input = r$n_in, N_symbols = length(r$symbols),
                             N_unmapped = r$n_unmapped, stringsAsFactors = FALSE)
  if (length(r$symbols) == 0) { cat("   ", nm, ": no significant genes, left out of the plot\n"); next }
  deg_lists[[nm]] <- r$symbols
}
rep_df <- do.call(rbind, report)
cat("\n>>> gene identifiers after conversion to symbols:\n"); print(rep_df, row.names = FALSE)
if (length(deg_lists) < 2) stop("Fewer than two comparisons have significant genes; nothing to compare.")

out_dir <- file.path(ROOT, "cross_dataset_comparison")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
write.csv(rep_df, file.path(out_dir, "cross_dataset_id_conversion_report.csv"), row.names = FALSE)

## ---- shared / unique ----
shared_unique_degs(deg_lists,
                   out_csv = file.path(out_dir, "cross_dataset_shared_unique_DEGs.csv"),
                   out_png = file.path(out_dir, "figures", "cross_dataset_upset.png"),
                   title = "Shared and unique DEGs across TB datasets")

## ---- plain-language summary next to the plot ----
tab <- read.csv(file.path(out_dir, "cross_dataset_shared_unique_DEGs.csv"), stringsAsFactors = FALSE)
cols <- setdiff(names(tab), "gene")
tab$N_datasets <- rowSums(tab[, cols, drop = FALSE])
write.csv(tab[order(-tab$N_datasets), ], file.path(out_dir, "cross_dataset_shared_unique_DEGs.csv"), row.names = FALSE)
cat("\n>>> genes by how many comparisons they are significant in:\n"); print(table(tab$N_datasets))
shared <- tab$gene[tab$N_datasets >= 2]
cat("\n>>> significant in two or more comparisons:", length(shared), "genes\n")
if (length(shared)) cat("   ", paste(utils::head(shared, 25), collapse = ", "),
                        if (length(shared) > 25) paste0(" ... (+", length(shared) - 25, " more)") else "", "\n")
cat("\n>>> DONE. Files in", normalizePath(out_dir, mustWork = FALSE), "\n")
cat(">>> Note: the comparisons differ greatly in size (from 0 to ~300 significant genes), so a small\n")
cat("    overlap mostly reflects that, not a biological disagreement.\n")
