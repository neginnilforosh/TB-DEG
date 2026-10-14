## STEP 20: bridge the metabolomics layer to the gene layer through ENZYMES, not pathway membership.
##
## Why: joining the two layers on shared KEGG pathways gave zero overlap, for a structural reason --
## the gene modules enrich immune-signalling pathways, which hold almost no metabolites, while the
## metabolites map to metabolic pathways, which hold almost no immune genes. The direct biochemical link
## is the enzyme: a metabolite is made or used by specific enzymes, and those enzymes are genes whose
## expression we measured. Tryptophan -> IDO1 (an interferon-inducible enzyme) is the obvious example.
##
## For every metabolite with a KEGG id (step 19): KEGG compound -> EC numbers -> human genes. Then:
##   1. a bridge table: metabolite change next to the expression change, module and priority of each enzyme
##   2. a test: are enzymes of the SIGNIFICANT metabolites more often differentially expressed / in the
##      disease modules than enzymes of the other MEASURED metabolites? (enzymes vs enzymes, so the test is
##      not just "metabolic genes differ from immune genes")
## Caveat printed with the results: the metabolites come from a plasma progression cohort (ST000974,
## M06) and the genes from a blood transcriptome (this run), so this is consistency across cohorts,
## not a paired measurement.
## Run from scripts/wgcna_network after step 19. Outputs in <run>/22_metabolite_gene_bridge/.

## ---- CONFIG ----
STUDY_ID     <- "ST000974"
MET_PADJ     <- 0.05     # "significant" metabolites
MET_P_BROAD  <- 0.01     # broader set (nominal), reported alongside -- there are few FDR hits
GENE_PADJ    <- 0.05; GENE_LFC <- 1
KEGG_SLEEP   <- 0.34     # KEGG asks for <= 3 requests per second

for (p in c("httr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

## ---- paths ----
get_script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR <- file.path(ROOT_DIR, ACC); OUT_DIR <- run_out_dir(BASE_DIR)
dir_out  <- file.path(OUT_DIR, "22_metabolite_gene_bridge"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) { if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE) }

met <- rd(file.path(in_dir(BASE_DIR, "21_metabolomics"), paste0(STUDY_ID, "_differential_metabolites_with_KEGG.csv")))
met <- met[!is.na(met$kegg_id) & grepl("^C[0-9]{5}$", met$kegg_id), ]
met <- met[order(met$pvalue), ]; met <- met[!duplicated(met$kegg_id), ]
cat(">>>", nrow(met), "measured metabolites with a KEGG compound id |",
    sum(met$padj < MET_PADJ), "significant (padj <", MET_PADJ, ") |", sum(met$pvalue < MET_P_BROAD), "with p <", MET_P_BROAD, "\n")

## ---- KEGG, cached ----
cache_file <- file.path(dir_out, "kegg_cache.rds")
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
kegg_get <- function(path) {                     # returns a 2-column data frame (from, to) or NULL
  if (!is.null(cache[[path]])) return(cache[[path]])
  Sys.sleep(KEGG_SLEEP)
  r <- tryCatch(httr::GET(paste0("https://rest.kegg.jp/", path), httr::timeout(60)), error = function(e) NULL)
  if (is.null(r) || !httr::status_code(r) %in% c(200, 404)) stop("KEGG request failed: ", path, " -- re-run to retry (results so far are cached)")
  txt <- if (httr::status_code(r) == 404) "" else httr::content(r, "text", encoding = "UTF-8")
  lines <- strsplit(trimws(txt), "\n")[[1]]; lines <- lines[nzchar(lines)]
  out <- if (length(lines)) do.call(rbind, lapply(strsplit(lines, "\t"), function(p) data.frame(from = p[1], to = p[2]))) else
    data.frame(from = character(0), to = character(0))
  cache[[path]] <<- out; saveRDS(cache, cache_file); out
}
batched <- function(ids, prefix, target) {       # KEGG accepts up to 10 ids joined by '+'
  ids <- sort(unique(ids)); if (!length(ids)) return(data.frame(from = character(0), to = character(0)))
  do.call(rbind, lapply(split(ids, ceiling(seq_along(ids) / 10)), function(b)
    kegg_get(paste0("link/", target, "/", paste0(prefix, b, collapse = "+")))))
}

cat(">>> KEGG: compound -> enzyme (EC)...\n")
c2e <- batched(met$kegg_id, "cpd:", "enzyme")
c2e$kegg_id <- sub("^cpd:", "", c2e$from); c2e$EC <- sub("^ec:", "", c2e$to)
cat("    ", length(unique(c2e$kegg_id)), "compounds have at least one enzyme |", length(unique(c2e$EC)), "EC numbers\n")
cat(">>> KEGG: enzyme -> human gene...\n")
e2g <- batched(unique(c2e$EC), "ec:", "hsa")
e2g$EC <- sub("^ec:", "", e2g$from); e2g$ENTREZID <- sub("^hsa:", "", e2g$to)
bridge <- merge(c2e[, c("kegg_id", "EC")], e2g[, c("EC", "ENTREZID")], by = "EC")
cat("    ", length(unique(bridge$ENTREZID)), "human enzyme genes linked to measured metabolites\n")
if (!nrow(bridge)) stop("No enzyme genes found -- nothing to bridge.")

## ---- genes: ids, expression, modules ----
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Missing package: org.Hs.eg.db")
ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(bridge$ENTREZID),
                        keytype = "ENTREZID", columns = c("SYMBOL", "ENSEMBL")))
bridge <- merge(bridge, ann, by = "ENTREZID")

deg <- rd(file.path(in_dir(BASE_DIR, "06_deg_results"), paste0(ACC, "_", DEG_NAME, "_DEG.csv")))
deg <- deg[!is.na(deg$padj), c("gene", "log2FoldChange", "padj")]
deg$DEG <- deg$padj < GENE_PADJ & abs(deg$log2FoldChange) >= GENE_LFC
pri <- file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv"))
pri <- if (file.exists(pri)) read.csv(pri, stringsAsFactors = FALSE)[, c("Ensembl", "Module", "Module_Name", "Final_Priority")] else NULL

bridge <- merge(bridge, deg, by.x = "ENSEMBL", by.y = "gene")              # keep enzymes that were measured
if (!is.null(pri)) bridge <- merge(bridge, pri, by.x = "ENSEMBL", by.y = "Ensembl", all.x = TRUE)
bridge$In_module <- if (is.null(pri)) NA else !is.na(bridge$Module)
bridge <- merge(met[, c("metabolite", "kegg_id", "log2FC", "pvalue", "padj")], bridge, by = "kegg_id",
                suffixes = c("_metabolite", "_gene"))
names(bridge) <- sub("^log2FC$", "Metabolite_log2FC", names(bridge))
names(bridge)[names(bridge) == "log2FoldChange"] <- "Gene_log2FC"
names(bridge)[names(bridge) == "pvalue"] <- "Metabolite_p"
names(bridge)[names(bridge) == "padj_metabolite"] <- "Metabolite_padj"
names(bridge)[names(bridge) == "padj_gene"] <- "Gene_padj"
bridge$Metabolite_significant <- bridge$Metabolite_padj < MET_PADJ
bridge <- bridge[order(bridge$Metabolite_p, bridge$Gene_padj), ]
keep <- c("metabolite", "kegg_id", "Metabolite_log2FC", "Metabolite_p", "Metabolite_padj", "Metabolite_significant",
          "EC", "SYMBOL", "ENSEMBL", "Gene_log2FC", "Gene_padj", "DEG", "In_module", "Module", "Module_Name", "Final_Priority")
bridge <- unique(bridge[, intersect(keep, names(bridge))])
write.csv(bridge, file.path(dir_out, paste0(STUDY_ID, "_", ACC, "_metabolite_enzyme_bridge.csv")), row.names = FALSE)

## ---- what the significant metabolites connect to ----
sig <- bridge[bridge$Metabolite_significant, ]
cat("\n========== significant metabolites and their enzyme genes ==========\n")
if (!nrow(sig)) cat("  none of the significant metabolites has a measured enzyme gene\n") else
  for (m in unique(sig$metabolite)) {
    s <- sig[sig$metabolite == m, ]
    cat(sprintf("  %s (%s, log2FC %+.2f): ", m, s$kegg_id[1], s$Metabolite_log2FC[1]))
    cat(paste0(s$SYMBOL, ifelse(s$DEG, sprintf("[DEG %+.1f]", s$Gene_log2FC), ""),
               ifelse(!is.na(s$Module), paste0("{", s$Module, "}"), ""), collapse = ", "), "\n")
  }

## ---- test: enzymes of changed metabolites vs enzymes of the other measured metabolites ----
gene_level <- function(changed) {
  ## one row per enzyme gene; "changed" = linked to at least one changed metabolite
  ch <- unique(bridge$ENSEMBL[changed]); all_g <- unique(bridge[, c("ENSEMBL", "DEG", "In_module")])
  all_g$Linked_to_changed <- all_g$ENSEMBL %in% ch; all_g
}
tests <- list()
for (set in list(list(nm = paste0("padj<", MET_PADJ), v = bridge$Metabolite_padj < MET_PADJ),
                 list(nm = paste0("p<", MET_P_BROAD), v = bridge$Metabolite_p < MET_P_BROAD))) {
  gl <- gene_level(set$v)
  for (out in c("DEG", "In_module")) {
    if (all(is.na(gl[[out]]))) next
    tb <- table(factor(gl$Linked_to_changed, c(TRUE, FALSE)), factor(gl[[out]], c(TRUE, FALSE)))
    if (sum(tb[1, ]) == 0) next
    ft <- stats::fisher.test(tb, alternative = "greater")
    tests[[length(tests) + 1]] <- data.frame(Metabolite_set = set$nm, Outcome = out,
      Enzymes_linked_to_changed = sum(tb[1, ]), Of_which_outcome = tb[1, 1],
      Pct_changed = round(100 * tb[1, 1] / sum(tb[1, ]), 1),
      Enzymes_other = sum(tb[2, ]), Pct_other = round(100 * tb[2, 1] / max(1, sum(tb[2, ])), 1),
      Odds_ratio = unname(round(ft$estimate, 2)), P_value = signif(ft$p.value, 3))
  }
}
if (length(tests)) {
  tt <- do.call(rbind, tests)
  write.csv(tt, file.path(dir_out, paste0(STUDY_ID, "_", ACC, "_enzyme_enrichment_tests.csv")), row.names = FALSE)
  cat("\n========== are enzymes of changed metabolites more often changed themselves? ==========\n")
  print(tt, row.names = FALSE)
  cat("(background = enzymes of the other measured metabolites, so this compares enzymes with enzymes)\n")
}
cat("\nCaveat: metabolites from a plasma progression cohort (", STUDY_ID, ", M06) and genes from this run's blood\n",
    "transcriptome are different people -- read this as cross-cohort consistency, not a paired measurement.\n", sep = "")
cat(">>> DONE. Outputs in", dir_out, "\n")
