## =============================================================
## STEP 10 (FINAL) of Ratul's roadmap: two master tables.
##
## "The first should be a gene table containing: Gene, log2FC, FDR,
##  Module, Module biological name, MM/kME, Gene Significance, PPI
##  degree, Betweenness, RWR score, Diffusion score, and Pathways.
##  The second should be a drug table containing: Drug, iLINCS
##  reversal score, Known targets, Network proximity, RWR score,
##  Modules affected, and Pathways affected."
##
## Pure local data wrangling -- everything needed was already produced
## by steps 1-9, no new network calls.
## =============================================================

required_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse=", "),
                       " (already installed for steps 2/4/5/6/9)")

## ---- CONFIG ----
ACC <- "GSE114192"
## Fill in blue's name here once you have one (Hallmark run, or a manual call) --
## everything below carries NA through gracefully until then.
MODULE_NAMES <- c(green = "Interferon_Response", purple = "Ribosome_Biogenesis_Chromatin", blue = NA)


## ---- name harmonization helper (same as step 9) ----
## Network node names (STRING preferredName) can differ from the current org.Hs.eg.db SYMBOL
## for the same gene (e.g. RIGI vs DDX58). Joining PPI/RWR columns on SYMBOL alone would
## leave exactly those genes -- often important hubs -- with empty PPI/RWR cells, so:
## exact SYMBOL match first, ALIAS as a guarded fallback (see below).
ensembl_to_network_name <- function(ens_ids, network_names, annot = NULL) {
  ens_ids <- unique(ens_ids)
  if (is.null(annot)) {
    annot <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ens_ids,
                                                     keytype = "ENSEMBL", columns = c("SYMBOL", "ALIAS")))
  }
  sym_by_ens   <- split(annot$SYMBOL, annot$ENSEMBL)
  alias_by_ens <- split(annot$ALIAS,  annot$ENSEMBL)
  syms    <- lapply(ens_ids, function(e) unique(stats::na.omit(sym_by_ens[[e]])))
  aliases <- lapply(ens_ids, function(e) unique(stats::na.omit(alias_by_ens[[e]])))
  first_hit <- function(x) { h <- intersect(x, network_names); if (length(h)) h[1] else NA_character_ }
  exact  <- vapply(syms,    first_hit, character(1))   # 1) exact SYMBOL match
  viaal  <- vapply(aliases, first_hit, character(1))   # 2) ALIAS match (fallback only)
  nm <- exact
  use_alias <- is.na(exact) & !is.na(viaal)
  ## An alias hit is only trusted when it can't be a false friend: the network node must NOT already
  ## be claimed by another gene's exact symbol, and must not be claimed by 2+ different genes' aliases.
  claimed      <- unique(stats::na.omit(exact))
  alias_counts <- table(viaal[use_alias])
  unique_alias <- names(alias_counts)[alias_counts == 1]
  ok <- use_alias & !(viaal %in% claimed) & (viaal %in% unique_alias)
  nm[ok] <- viaal[ok]
  first_sym <- vapply(syms, function(s) if (length(s)) s[1] else NA_character_, character(1))
  data.frame(ENSEMBL = ens_ids, SYMBOL = unname(first_sym), NetName = unname(nm), stringsAsFactors = FALSE)
}


## ---- drug-name key (used to join DGIdb drug names onto iLINCS compound names) ----
## DGIdb names are UPPERCASE and often carry salt/hydrate suffixes ("AMLODIPINE BESYLATE",
## "AFURESERTIB HYDROCHLORIDE"); iLINCS uses mixed case and the bare name ("Amlodipine").
## An exact string match therefore found only 32 of 395 drugs; upper-casing gets 297, and
## also stripping trailing salt/hydrate words gets 367 (measured on the real GSE114192 files).
SALT_WORDS <- paste0("(DIHYDROCHLORIDE|HYDROCHLORIDE|HYDROBROMIDE|HCL|MESYLATE|MESILATE|BESYLATE|BESILATE|",
                     "MALEATE|FUMARATE|CITRATE|SULFATE|SULPHATE|PHOSPHATE|TARTRATE|ACETATE|SUCCINATE|",
                     "TOSYLATE|TOSILATE|LACTATE|GLUCONATE|NITRATE|SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|",
                     "CHLORIDE|BROMIDE|ANHYDROUS|MONOHYDRATE|DIHYDRATE|TRIHYDRATE|HEMIHYDRATE|HYDRATE)")
drug_key <- function(x) {
  x <- toupper(trimws(x))
  repeat {
    y <- trimws(sub(paste0("[[:space:],()-]+", SALT_WORDS, "$"), "", x))
    if (identical(y, x)) break
    x <- y
  }
  x
}

## ---- locate script dir ----
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
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found
dir_deg    <- file.path(BASE_DIR, "06_deg_results")
dir_wgcna  <- file.path(BASE_DIR, "10_wgcna_results")
dir_enrich <- file.path(BASE_DIR, "11_module_enrichment")
dir_string <- file.path(BASE_DIR, "12_string_ppi")
dir_sig    <- file.path(BASE_DIR, "13_ilincs_signature")
dir_dt     <- file.path(BASE_DIR, "14_drug_targets")
dir_net    <- file.path(BASE_DIR, "15_network_algorithms")
dir_out    <- file.path(BASE_DIR, "16_final_tables")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## =============================================================
## GENE TABLE
## =============================================================
cat(">>> Building gene table...\n")
deg    <- read.csv(file.path(dir_deg, paste0(ACC, "_HealthyControl_vs_TBOnly_DEG.csv")), stringsAsFactors = FALSE)
mm_gs  <- read.csv(file.path(dir_wgcna, paste0(ACC, "_TBmodules_MM_GS.csv")), stringsAsFactors = FALSE)  # Ensembl IDs
cent   <- read.csv(file.path(dir_string, paste0(ACC, "_COMBINED_centrality.csv")), stringsAsFactors = FALSE)  # network (STRING) names
scores <- read.csv(file.path(dir_net, paste0(ACC, "_GeneScores_RWR_Diffusion.csv")), stringsAsFactors = FALSE)  # network (STRING) names

## Ensembl -> (current symbol, network node name)
name_map <- ensembl_to_network_name(mm_gs$Gene, union(cent$Gene, scores$Gene))
cat(">>> ", sum(!is.na(name_map$NetName)), "/", nrow(name_map), "module genes matched to a node in the combined network",
    "(", sum(is.na(name_map$NetName)), "have no STRING node -- their PPI/RWR cells stay empty).\n")

gene_table <- merge(mm_gs, name_map, by.x = "Gene", by.y = "ENSEMBL", all.x = TRUE)   # adds SYMBOL, NetName
gene_table <- merge(gene_table, deg[, c("gene", "log2FoldChange", "padj")], by.x = "Gene", by.y = "gene", all.x = TRUE)
gene_table <- merge(gene_table, cent[, c("Gene", "Degree", "Betweenness")],
                    by.x = "NetName", by.y = "Gene", all.x = TRUE)
gene_table <- merge(gene_table, scores[, c("Gene", "RWR_score", "Diffusion_score")],
                    by.x = "NetName", by.y = "Gene", all.x = TRUE)
gene_table <- gene_table[!duplicated(gene_table$Gene), ]   # Gene here = Ensembl ID; one row per gene
gene_table$Module_Name <- MODULE_NAMES[gene_table$Module]

## Pathway per gene: is this gene among the members of its module's TOP GO_BP term?
## (enrichGO(readable=TRUE) result has a "geneID" column: "/"-separated org.Hs.eg.db symbols,
##  so this compares against SYMBOL, not the network name)
gene_table$Top_Pathway <- NA_character_
for (mod in unique(gene_table$Module)) {
  go_file <- file.path(dir_enrich, paste0(ACC, "_", mod, "_GO_BP.csv"))
  if (!file.exists(go_file)) next
  go <- read.csv(go_file, stringsAsFactors = FALSE)
  go <- go[order(go$p.adjust), ]
  if (nrow(go) == 0 || is.na(go$p.adjust[1]) || go$p.adjust[1] > 0.05) next  # no significant term for this module (e.g. blue)
  top_term_genes <- strsplit(go$geneID[1], "/")[[1]]
  idx <- gene_table$Module == mod & gene_table$SYMBOL %in% top_term_genes
  gene_table$Top_Pathway[idx] <- go$Description[1]
}

gene_table_out <- data.frame(
  Gene = gene_table$SYMBOL, Network_Name = gene_table$NetName, Ensembl = gene_table$Gene, log2FC = gene_table$log2FoldChange,
  FDR = gene_table$padj, Module = gene_table$Module, Module_Name = gene_table$Module_Name,
  MM = gene_table$MM, Gene_Significance = gene_table$GS_TB,
  PPI_Degree = gene_table$Degree, PPI_Betweenness = gene_table$Betweenness,
  RWR_score = gene_table$RWR_score, Diffusion_score = gene_table$Diffusion_score,
  Top_Pathway = gene_table$Top_Pathway
)
gene_table_out <- gene_table_out[order(gene_table_out$Module, -gene_table_out$RWR_score), ]
write.csv(gene_table_out, file.path(dir_out, paste0(ACC, "_FINAL_GeneTable.csv")), row.names = FALSE)
cat(">>> Gene table:", nrow(gene_table_out), "genes.\n")

## =============================================================
## DRUG TABLE
## =============================================================
cat("\n>>> Building drug table...\n")
cand  <- read.csv(file.path(dir_sig, paste0(ACC, "_iLINCS_candidate_compounds.csv")), stringsAsFactors = FALSE)
dt    <- read.csv(file.path(dir_dt, paste0(ACC, "_Drug_Target_TBnetwork.csv")), stringsAsFactors = FALSE)
prox  <- read.csv(file.path(dir_net, paste0(ACC, "_NetworkProximity.csv")), stringsAsFactors = FALSE)
cand$.key <- drug_key(cand$Target)

drugs <- unique(dt$Drug)  # scope = drugs that have at least one target IN the TB network
rows <- list()
for (drug in drugs) {
  d_sub <- dt[dt$Drug == drug, ]
  targets <- unique(d_sub$Drug_Target)
  modules_affected <- unique(d_sub$TB_Module)
  pathways_affected <- unique(na.omit(MODULE_NAMES[modules_affected]))

  sim_rows <- cand[cand$.key == drug_key(drug), ]
  reversal_score <- if (nrow(sim_rows)) min(sim_rows$Similarity) else NA  # most negative = strongest reversal

  prox_row <- prox[prox$Drug == drug, ]
  proximity_z <- if (nrow(prox_row)) prox_row$Z_score[1] else NA
  proximity_p <- if (nrow(prox_row) && "P_empirical" %in% names(prox_row)) prox_row$P_empirical[1] else NA   # present from the fixed step 9 on

  target_scores <- scores$RWR_score[scores$Gene %in% targets]
  rwr_drug_score <- if (length(target_scores)) mean(target_scores, na.rm = TRUE) else NA

  rows[[length(rows) + 1]] <- data.frame(
    Drug = drug, iLINCS_Reversal_Score = reversal_score, Known_Targets = paste(targets, collapse = ";"),
    Network_Proximity_Z = proximity_z, Network_Proximity_P = proximity_p, RWR_score = rwr_drug_score,
    Modules_Affected = paste(modules_affected, collapse = ";"),
    Pathways_Affected = if (length(pathways_affected)) paste(pathways_affected, collapse = ";") else NA_character_
  )
}
drug_table_out <- do.call(rbind, rows)
drug_table_out <- drug_table_out[order(drug_table_out$iLINCS_Reversal_Score), ]  # strongest reversal first

## ---- optional: per-drug pathways from Enrichr on the drug's known targets (made by step 12) ----
enr_file <- file.path(BASE_DIR, "17_pathway_annotation", paste0(ACC, "_DrugPathways_Enrichr_summary.csv"))
if (file.exists(enr_file)) {
  enr <- read.csv(enr_file, stringsAsFactors = FALSE)
  drug_table_out <- merge(drug_table_out, enr, by = "Drug", all.x = TRUE, sort = FALSE)
  drug_table_out <- drug_table_out[order(drug_table_out$iLINCS_Reversal_Score), ]
  cat(">>> Merged Enrichr pathway summary:", sum(!is.na(drug_table_out$Enrichr_N_targets_used)), "drugs have Enrichr results.\n")
} else {
  cat(">>> (No Enrichr summary yet -- run WGCNA_Enrichr_drugs_step12.R, then re-run this script to add the pathway columns.)\n")
}

## ---- optional: cell-type profile of each drug's targets (made by step 13, macaque TB-granuloma single-cell data) ----
sc_file <- file.path(BASE_DIR, "18_singlecell", paste0(ACC, "_drug_celltype.csv"))
if (file.exists(sc_file)) {
  sc <- read.csv(sc_file, stringsAsFactors = FALSE)
  sc <- sc[, c("Drug", "N_targets_detected_sc", "SC_Top_CellType", "SC_Top_Share_pct", "SC_Top_Immune_CellType", "SC_Top_Immune_Share_pct", "SC_Profile")]
  drug_table_out <- merge(drug_table_out, sc, by = "Drug", all.x = TRUE, sort = FALSE)
  drug_table_out <- drug_table_out[order(drug_table_out$iLINCS_Reversal_Score), ]
  cat(">>> Merged single-cell cell-type profile:", sum(!is.na(drug_table_out$SC_Top_CellType)), "drugs have one.\n")
} else {
  cat(">>> (No single-cell profile yet -- run WGCNA_singlecell_step13.R, then re-run this script.)\n")
}
write.csv(drug_table_out, file.path(dir_out, paste0(ACC, "_FINAL_DrugTable.csv")), row.names = FALSE)
cat(">>> Drug table:", nrow(drug_table_out), "drugs.\n")
cat(">>>", sum(!is.na(drug_table_out$iLINCS_Reversal_Score)), "/", nrow(drug_table_out),
    "have an iLINCS reversal score;", sum(is.finite(drug_table_out$Network_Proximity_Z)), "have a finite proximity Z.\n")
if (mean(is.na(drug_table_out$iLINCS_Reversal_Score)) > 0.3)
  cat("!! More than 30% of drugs have no iLINCS score -- check the drug-name join before trusting the ranking.\n")

cat("\n>>> DONE. Both final tables in", dir_out, ":\n")
cat("  -", paste0(ACC, "_FINAL_GeneTable.csv"))
cat("\n  -", paste0(ACC, "_FINAL_DrugTable.csv"), "\n")
cat("\n>>> Reminder: blue's Module_Name/Pathways_Affected are NA until it gets a real\n")
cat("    biological name (Hallmark run, or manual review) -- everything else already\n")
cat("    carries through correctly regardless.\n")
