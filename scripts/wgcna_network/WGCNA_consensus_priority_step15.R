## STEP 15: consensus prioritisation. Genes and drugs are ranked by agreement across several independent
## criteria, not by one score. Reads the step 10/9/3 outputs of the active run; base R only.
## Outputs in <run>/19_consensus/.

## ---- CONFIG ----
TOP_FRACTION   <- 0.20    # "top-ranked" = best 20% of a criterion
MERGE_RWR_DIFF <- TRUE    # RWR and diffusion are near-duplicates (correlation printed below): count them as ONE vote
HIGH_N         <- 4       # consensus tiers (number of criteria met)
MED_N          <- 2
PROX_Z_CUT     <- -1      # drug: targets closer to the disease genes than average
MIN_TARGETS    <- 3       # drug: targets inside the TB network
P_CUT          <- 0.05    # drug: proximity and module-enrichment p-values
## drug-target cleaning before counting targets (set both to c()/FALSE to reproduce the original pilot numbers)
EXCLUDE_SCREEN_TARGETS <- c("NFE2L2")   # linked to >100 compounds by screening-assay records, not mechanism
EXCLUDE_ADME_TARGETS   <- TRUE          # drug-metabolism / transport genes (CYP, UGT, SULT, SLC, ABC)
ADME_PATTERN <- "^(CYP[0-9]|UGT[0-9]|SULT[0-9]|SLC[0-9]|ABC[A-G][0-9])"

## ---- paths ----
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
source(file.path(SCRIPT_DIR, "config.R"))
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })
OUT_DIR  <- run_out_dir(BASE_DIR)
dir_out  <- file.path(OUT_DIR, "19_consensus"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE) }

## ---- helper: flag the best TOP_FRACTION of a criterion (NA never counts) ----
top_flag <- function(x, higher_is_better = TRUE) {
  v <- if (higher_is_better) x else -x
  thr <- stats::quantile(v, 1 - TOP_FRACTION, na.rm = TRUE)
  !is.na(v) & v >= thr
}

## =============================================================
## GENES
## =============================================================
cat("\n========== GENES ==========\n")
g  <- rd(OUT_DIR, "16_final_tables", paste0(ACC, "_FINAL_GeneTable.csv"))
dg <- rd(OUT_DIR, "10_wgcna_results", paste0(ACC, "_DiseaseGenes_DEG_WGCNA.csv"))
g$DEG_WGCNA_overlap <- g$Ensembl %in% dg$Gene
g$In_PPI_network    <- !is.na(g$PPI_Degree)
cat(">>> ", nrow(g), "genes in the selected modules (", paste(names(table(g$Module)), table(g$Module), collapse = ", "),
    "); of these", sum(g$DEG_WGCNA_overlap), "are also significant DEGs.\n")
cat(">>> ", sum(!g$In_PPI_network), "genes have no node in the STRING network, so they have no PPI/RWR/diffusion score",
    "and can meet at most the two WGCNA criteria -- column In_PPI_network keeps that visible.\n")

## the five criteria; RWR and diffusion are reported separately but may count as one vote
crit <- list(Top_MM = top_flag(abs(g$MM)), Top_GS = top_flag(abs(g$Gene_Significance)),
             Top_PPI_Degree = top_flag(g$PPI_Degree), Top_RWR = top_flag(g$RWR_score),
             Top_Diffusion = top_flag(g$Diffusion_score))
for (n in names(crit)) g[[n]] <- crit[[n]]
## per-metric rank (1 = best; ties share the best rank; genes without a value stay NA), so each
## Top_* decision can be traced back to where the gene actually sits on that metric
rank_desc <- function(x) rank(-x, ties.method = "min", na.last = "keep")
g$Rank_MM             <- rank_desc(abs(g$MM))
g$Rank_GS             <- rank_desc(abs(g$Gene_Significance))
g$Rank_PPI_Degree     <- rank_desc(g$PPI_Degree)
g$Rank_PPI_Betweenness <- rank_desc(g$PPI_Betweenness)   # reported only; not a consensus criterion
g$Rank_RWR            <- rank_desc(g$RWR_score)
g$Rank_Diffusion      <- rank_desc(g$Diffusion_score)
rho <- suppressWarnings(stats::cor(g$RWR_score, g$Diffusion_score, method = "spearman", use = "complete.obs"))
cat(">>> RWR vs diffusion rank correlation:", round(rho, 2),
    if (MERGE_RWR_DIFF) " -- counted as ONE criterion (MERGE_RWR_DIFF = TRUE)\n" else " -- counted separately\n")
votes <- if (MERGE_RWR_DIFF) {
  cbind(crit$Top_MM, crit$Top_GS, crit$Top_PPI_Degree, crit$Top_RWR | crit$Top_Diffusion)
} else do.call(cbind, crit)
g$N_criteria <- rowSums(votes)
g$Criteria_met <- apply(as.data.frame(crit), 1, function(r) paste(sub("^Top_", "", names(crit))[r], collapse = "+"))
g$Final_Priority <- ifelse(g$N_criteria >= HIGH_N, "High", ifelse(g$N_criteria >= MED_N, "Medium", "Low"))
g$Direction_in_TB <- ifelse(g$log2FC > 0, "up", "down")

cols <- c("Gene", "Ensembl", "Module", "Module_Name", "log2FC", "FDR", "Direction_in_TB", "MM", "Gene_Significance",
          "PPI_Degree", "PPI_Betweenness", "RWR_score", "Diffusion_score",
          "Rank_MM", "Rank_GS", "Rank_PPI_Degree", "Rank_PPI_Betweenness", "Rank_RWR", "Rank_Diffusion",
          "Top_MM", "Top_GS", "Top_PPI_Degree", "Top_RWR", "Top_Diffusion", "N_criteria", "Criteria_met", "Final_Priority", "DEG_WGCNA_overlap", "In_PPI_network")
gene_out <- g[order(-g$N_criteria, g$FDR), intersect(cols, names(g))]
write.csv(gene_out, file.path(dir_out, paste0(ACC, "_GenePriority_all_modules.csv")), row.names = FALSE)
write.csv(gene_out[gene_out$DEG_WGCNA_overlap, ], file.path(dir_out, paste0(ACC, "_GenePriority_DEG_overlap.csv")), row.names = FALSE)
cat(">>> Priority (all module genes):", paste(names(table(gene_out$Final_Priority)), table(gene_out$Final_Priority), collapse = ", "), "\n")
sub188 <- gene_out[gene_out$DEG_WGCNA_overlap, ]
cat(">>> Priority (DEG-overlap genes only):", paste(names(table(sub188$Final_Priority)), table(sub188$Final_Priority), collapse = ", "), "\n")
cat(">>> High-priority genes per module:", paste(names(table(gene_out$Module[gene_out$Final_Priority == "High"])),
                                                  table(gene_out$Module[gene_out$Final_Priority == "High"]), collapse = ", "), "\n")
cat("\ntop 10 by consensus:\n"); print(utils::head(gene_out[, c("Gene", "Module", "log2FC", "FDR", "N_criteria", "Final_Priority", "DEG_WGCNA_overlap")], 10), row.names = FALSE)

## =============================================================
## DRUGS
## =============================================================
cat("\n========== DRUGS ==========\n")
d  <- rd(OUT_DIR, "16_final_tables", paste0(ACC, "_FINAL_DrugTable.csv"))
mo <- rd(OUT_DIR, "15_network_algorithms", paste0(ACC, "_ModuleOverlapScoring.csv"))

clean_targets <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(EXCLUDE_SCREEN_TARGETS)) x <- setdiff(x, EXCLUDE_SCREEN_TARGETS)
  if (isTRUE(EXCLUDE_ADME_TARGETS)) x <- x[!grepl(ADME_PATTERN, x)]
  x
}
raw_targets <- strsplit(d$Known_Targets, ";")
targets <- lapply(raw_targets, clean_targets)
n_raw <- vapply(raw_targets, function(x) length(unique(x[nzchar(x)])), integer(1))
d$N_TB_Network_Targets <- vapply(targets, length, integer(1))
cat(">>> target cleaning: excluded", if (length(EXCLUDE_SCREEN_TARGETS)) paste(EXCLUDE_SCREEN_TARGETS, collapse = ",") else "none",
    "| ADME genes excluded:", isTRUE(EXCLUDE_ADME_TARGETS), "\n")
cat("    drugs whose network targets all came from excluded genes:", sum(n_raw > 0 & d$N_TB_Network_Targets == 0),
    "| drugs with >=", MIN_TARGETS, "targets: before", sum(n_raw >= MIN_TARGETS), "after", sum(d$N_TB_Network_Targets >= MIN_TARGETS), "\n")
d$Known_Targets_used <- vapply(targets, paste, character(1), collapse = ";")
## informational only: share of the drug's (cleaned) targets that fall in the disease network
all_t <- file.path(OUT_DIR, "14_drug_targets", paste0(ACC, "_drug_targets_ALL.csv"))
if (file.exists(all_t)) {
  at <- read.csv(all_t, stringsAsFactors = FALSE)
  tot <- vapply(split(at$Target_Gene, at$Drug), function(x) length(clean_targets(x)), integer(1))
  d$N_targets_total <- unname(tot[d$Drug])
  d$Target_fraction_in_network <- round(d$N_TB_Network_Targets / d$N_targets_total, 3)
}
high_genes <- gene_out$Gene[gene_out$Final_Priority == "High"]
d$N_HighPriority_Targets <- vapply(targets, function(x) length(intersect(unique(x), high_genes)), integer(1))
d$HighPriority_Targets   <- vapply(targets, function(x) paste(intersect(unique(x), high_genes), collapse = ";"), character(1))
best <- stats::aggregate(P_value ~ Drug, mo, min); names(best)[2] <- "Module_Enrichment_P"
mods <- stats::aggregate(Module ~ Drug, mo[mo$P_value < P_CUT, , drop = FALSE], function(x) paste(sort(unique(x)), collapse = ";"))
names(mods)[2] <- "Modules_Enriched"
d <- merge(merge(d, best, by = "Drug", all.x = TRUE), mods, by = "Drug", all.x = TRUE)

dcrit <- list(
  Strong_reversal   = top_flag(d$iLINCS_Reversal_Score, higher_is_better = FALSE),   # most negative 20%
  Negative_prox_Z   = !is.na(d$Network_Proximity_Z) & d$Network_Proximity_Z <= PROX_Z_CUT,
  Significant_prox_P= !is.na(d$Network_Proximity_P) & d$Network_Proximity_P < P_CUT,
  Targets_in_TBnet  = d$N_TB_Network_Targets >= MIN_TARGETS,
  Module_enriched   = !is.na(d$Module_Enrichment_P) & d$Module_Enrichment_P < P_CUT,
  Hits_high_gene    = d$N_HighPriority_Targets >= 1,
  Top_target_RWR    = top_flag(d$RWR_score))   # targets well connected to the disease seeds (mean RWR of the drug's targets)
## drug ranks for the continuous criteria (1 = best); threshold criteria (proximity Z/P, counts) are in their own columns
d$Rank_iLINCS_reversal <- rank(d$iLINCS_Reversal_Score, ties.method = "min", na.last = "keep")   # most negative first
d$Rank_Proximity_Z     <- rank(d$Network_Proximity_Z,   ties.method = "min", na.last = "keep")   # most negative first
d$Rank_Target_RWR      <- rank(-d$RWR_score,            ties.method = "min", na.last = "keep")
for (n in names(dcrit)) d[[n]] <- dcrit[[n]]
d$N_criteria <- rowSums(do.call(cbind, dcrit))
d$Criteria_met <- apply(as.data.frame(dcrit), 1, function(r) paste(names(dcrit)[r], collapse = "+"))
d$Final_Priority <- ifelse(d$N_criteria >= HIGH_N, "High", ifelse(d$N_criteria >= MED_N, "Medium", "Low"))
cat(">>> drugs meeting each criterion (of", nrow(d), "):",
    paste(names(dcrit), vapply(dcrit, sum, integer(1)), collapse = ", "), "\n")

dcols <- c("Drug", "iLINCS_Reversal_Score", "Rank_iLINCS_reversal", "Network_Proximity_Z", "Rank_Proximity_Z", "Network_Proximity_P", "N_TB_Network_Targets",
           "Known_Targets", "Known_Targets_used", "N_targets_total", "Target_fraction_in_network", "Modules_Affected", "Modules_Enriched", "Module_Enrichment_P",
           "N_HighPriority_Targets", "HighPriority_Targets", "RWR_score", "Rank_Target_RWR", "Enrichr_KEGG", "Enrichr_Reactome",
           "SC_Top_Immune_CellType", names(dcrit), "N_criteria", "Criteria_met", "Final_Priority")
drug_out <- d[order(-d$N_criteria, d$iLINCS_Reversal_Score), intersect(dcols, names(d))]
write.csv(drug_out, file.path(dir_out, paste0(ACC, "_DrugPriority.csv")), row.names = FALSE)
cat(">>> Priority:", paste(names(table(drug_out$Final_Priority)), table(drug_out$Final_Priority), collapse = ", "), "\n")
cat("\ntop 10 by consensus:\n")
print(utils::head(drug_out[, c("Drug", "iLINCS_Reversal_Score", "Network_Proximity_Z", "N_TB_Network_Targets",
                               "N_HighPriority_Targets", "N_criteria", "Final_Priority")], 10), row.names = FALSE)

## ---- caveats, written next to the tables so the numbers are not read alone ----
n_padj <- sum(stats::p.adjust(d$Network_Proximity_P, "BH") < P_CUT, na.rm = TRUE)
notes <- c(
  paste0("Consensus prioritisation for ", ACC, " (", CONTROL, " vs ", CASE, "), ", format(Sys.Date()), "."),
  paste0("Top-ranked = best ", round(100 * TOP_FRACTION), "% of a criterion. High = ", HIGH_N,
         "+ criteria, Medium = ", MED_N, "-", HIGH_N - 1, ", Low = below that."), "",
  "Read with these caveats:",
  paste0("- RWR and diffusion rank-correlate ", round(rho, 2), ", so they are not independent evidence; they ",
         if (MERGE_RWR_DIFF) "count as one criterion here." else "count separately here, which inflates N_criteria."),
  paste0("- ", sum(!g$In_PPI_network), " of ", nrow(g), " module genes have no STRING node and can reach at most ",
         "two criteria; a Low tier there means 'not assessable', not 'unimportant' (see In_PPI_network)."),
  paste0("- Network proximity is weak in this run: ", sum(dcrit$Significant_prox_P), " of ", nrow(d),
         " drugs reach p < ", P_CUT, " and ", n_padj, " survive FDR correction. Treat it as a tie-breaker, not evidence."),
  paste0("- Module enrichment reaches p < ", P_CUT, " for ", sum(dcrit$Module_enriched), " drugs; most drugs have ",
         "too few targets for this test."),
  paste0("- Drug targets: ", if (length(EXCLUDE_SCREEN_TARGETS)) paste(EXCLUDE_SCREEN_TARGETS, collapse = ", ") else "no screening-assay genes",
         " excluded", if (isTRUE(EXCLUDE_ADME_TARGETS)) " and drug-metabolism/transport genes (CYP, UGT, SULT, SLC, ABC) excluded" else "",
         " before counting; ", sum(n_raw > 0 & d$N_TB_Network_Targets == 0), " drugs lost all their network targets this way."),
  ## module labels come from this run's MODULE_NAMES in config.R, so they are never another run's labels
  if (length(MODULE_NAMES))
    paste0("- Module labels (from config.R): ", paste(names(MODULE_NAMES), "=", MODULE_NAMES, collapse = "; "), ".")
  else "- Module labels: none set in config.R (MODULE_NAMES is empty), so modules are identified by WGCNA colour only.")
writeLines(notes, file.path(dir_out, paste0(ACC, "_CONSENSUS_NOTES.txt")))
cat("\n>>> DONE. Tables and ", paste0(ACC, "_CONSENSUS_NOTES.txt"), " in ", dir_out, "\n", sep = "")
