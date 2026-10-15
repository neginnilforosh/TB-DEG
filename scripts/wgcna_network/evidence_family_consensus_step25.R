## STEP 25: evidence-family consensus for genes and drugs.
##
## Why: step 15 counts yes/no criteria, and several criteria measure the same thing (degree, RWR and diffusion all
## come from one PPI network; proximity Z and P from one test). Close to a threshold a gene or drug flips tier for
## no real reason (one proximity Z moved from -1.00 to -1.03 and a drug became High). Here the evidence is grouped
## into independent FAMILIES, each family gives one continuous score between 0 and 1 (a percentile within the run),
## and the final score is the mean over the families a gene or drug has. Tiers are a readable summary of that score,
## and every run reports how much the ranking moves when one family is left out or a threshold changes.
##
## Gene families                                           Drug families
##   Expression     |log2FC| and -log10 FDR                  Reversal        iLINCS reversal strength
##   Co-expression  |module membership| and |gene signif.|   Target network  targets in the disease network: count and share
##   Network        degree, betweenness, RWR/diffusion       Target genes    best gene evidence score among the targets, target RWR
##   Replication    same direction and significance in       Module overlap  module enrichment of the targets
##                  independent active-TB cohorts            (Proximity      sensitivity only: never significant after FDR)
##   Infection      same direction and significance when
##                  macrophages are infected (GSE165709)
## Run after steps 15 and 21:  Rscript evidence_family_consensus_step25.R   (TB_RUN=... for other runs)
## Outputs in <run>/20_evidence_families/. Nothing else is changed.

## ---- CONFIG ----
MIN_FAMILIES_GENE <- 3      # a gene needs a score in at least this many families
MIN_FAMILIES_DRUG <- 3      # a drug needs reversal evidence plus at least this many families in total
STRONG    <- 0.80           # a family "supports" a gene/drug when its score is in the top 20%
HIGH_RULE <- list(strong_frac = 0.75, min_score = 0.75)   # proposed: strong in >= 75% of the families (genes 4 of 5, drugs 3 of 4)
MED_RULE  <- list(min_strong = 2, min_score = 0.60)    # change here; everything is recomputed
TOP_N_GENES <- 50; TOP_N_DRUGS <- 20
CAP_LOG10 <- 10             # cap -log10(p) so a handful of extreme p-values do not dominate
ACTIVE_PHENOTYPE <- "Active"

get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "20_evidence_families"); dir.create(dir_out, showWarnings = FALSE)
rd <- function(f, need = TRUE) { if (!file.exists(f)) { if (need) stop("Missing: ", f) else return(NULL) }
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

pct <- function(x) { r <- rep(NA_real_, length(x)); ok <- !is.na(x); if (sum(ok) > 1) r[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1) else r[ok] <- 1; r }
nlog <- function(p) pmin(-log10(pmax(p, 1e-300)), CAP_LOG10)
rmean <- function(...) { m <- cbind(...); out <- rowMeans(m, na.rm = TRUE); out[rowSums(!is.na(m)) == 0] <- NA; out }
jac <- function(a, b) if (!length(union(a, b))) NA else length(intersect(a, b)) / length(union(a, b))

## ======================= GENES =======================
g <- rd(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
g$Fam_Expression   <- rmean(pct(abs(g$log2FC)), pct(nlog(g$FDR)))
g$Fam_Coexpression <- rmean(pct(abs(g$MM)), pct(abs(g$Gene_Significance)))
g$Fam_Network      <- rmean(pct(g$PPI_Degree), pct(g$PPI_Betweenness), pct(rmean(pct(g$RWR_score), pct(g$Diffusion_score))))
g$Fam_Network[!g$In_PPI_network %in% TRUE] <- NA

## replication: independent active-TB cohorts (step 21). Agreement = sign(match) x -log10 padj, averaged over cohorts
ev <- rd(file.path(OUT_DIR, "23_external_validation", paste0(ACC, "_external_validation_long.csv")), need = FALSE)
if (!is.null(ev)) {
  ev <- ev[ev$Phenotype == ACTIVE_PHENOTYPE & !is.na(ev$Ext_log2FC) & !is.na(ev$Ext_padj), ]
  ev$agree <- ifelse(sign(ev$Ext_log2FC) == sign(ev$Discovery_log2FC), 1, -1) * nlog(ev$Ext_padj)
  a <- tapply(ev$agree, ev$Ensembl, mean); g$Fam_Replication <- pct(unname(a[g$Ensembl]))
} else g$Fam_Replication <- NA_real_

## infection response: GSE165709 macrophages from healthy donors, infected vs uninfected
inf <- rd(file.path(ROOT_DIR, "external_results", "GSE165709", "integration", "GSE165709_gene_integration.csv"), need = FALSE)
if (!is.null(inf)) {
  inf <- inf[!is.na(inf$RNA_log2FC_HC) & !is.na(inf$RNA_padj_HC), ]; inf <- inf[!duplicated(inf$Gene), ]
  m <- match(g$Gene, inf$Gene)
  agree <- ifelse(sign(inf$RNA_log2FC_HC[m]) == sign(g$log2FC), 1, -1) * nlog(inf$RNA_padj_HC[m])
  g$Fam_Infection <- pct(agree)
} else g$Fam_Infection <- NA_real_

gfam <- c("Fam_Expression", "Fam_Coexpression", "Fam_Network", "Fam_Replication", "Fam_Infection")
gfam <- gfam[vapply(gfam, function(f) any(!is.na(g[[f]])), logical(1))]
score_genes <- function(fams, strong = STRONG) {
  M <- as.matrix(g[, fams, drop = FALSE]); n_av <- rowSums(!is.na(M))
  s <- rowMeans(M, na.rm = TRUE); s[n_av < min(MIN_FAMILIES_GENE, length(fams))] <- NA
  list(score = s, n_available = n_av, n_strong = rowSums(M >= strong, na.rm = TRUE))
}
tier <- function(score, n_strong, n_fam) ifelse(is.na(score), "Insufficient evidence",
  ifelse(n_strong >= ceiling(HIGH_RULE$strong_frac * n_fam) & score >= HIGH_RULE$min_score, "High",
  ifelse(n_strong >= MED_RULE$min_strong | score >= MED_RULE$min_score, "Medium", "Low")))
main <- score_genes(gfam)
g$Families_available <- main$n_available; g$Families_strong <- main$n_strong
g$Evidence_score <- round(main$score, 4); g$Evidence_rank <- rank(-g$Evidence_score, ties.method = "min", na.last = "keep")
g$Evidence_tier <- tier(g$Evidence_score, g$Families_strong, length(gfam))

## sensitivity: leave one family out, and the "strong" threshold at 0.7 / 0.9
variants <- c(list(main = gfam), setNames(lapply(gfam, function(f) setdiff(gfam, f)), paste0("without_", sub("Fam_", "", gfam))))
top_of <- function(s, n) { o <- order(-s, na.last = NA); g$Ensembl[utils::head(o, n)] }
main_top <- top_of(g$Evidence_score, TOP_N_GENES); main_high <- g$Ensembl[g$Evidence_tier == "High"]
ranks <- sapply(variants, function(v) rank(-score_genes(v)$score, ties.method = "min", na.last = "keep"))
sens <- do.call(rbind, lapply(names(variants), function(v) {
  sv <- score_genes(variants[[v]]); t <- tier(sv$score, sv$n_strong, length(variants[[v]]))
  data.frame(Variant = v, Families = paste(sub("Fam_", "", variants[[v]]), collapse = "+"),
             Spearman_with_main = round(suppressWarnings(stats::cor(sv$score, g$Evidence_score, method = "spearman", use = "pairwise")), 3),
             TopN_overlap_with_main = sum(top_of(sv$score, TOP_N_GENES) %in% main_top), N_High = sum(t == "High"),
             High_Jaccard_with_main = round(jac(g$Ensembl[t == "High"], main_high), 3))
}))
for (thr in c(0.70, 0.90)) { sv <- score_genes(gfam, thr); t <- tier(sv$score, sv$n_strong, length(gfam))
  sens <- rbind(sens, data.frame(Variant = paste0("strong_threshold_", thr), Families = paste(sub("Fam_", "", gfam), collapse = "+"),
    Spearman_with_main = 1, TopN_overlap_with_main = TOP_N_GENES, N_High = sum(t == "High"),
    High_Jaccard_with_main = round(jac(g$Ensembl[t == "High"], main_high), 3))) }
g$Rank_best_across_variants  <- apply(ranks, 1, function(r) suppressWarnings(min(r, na.rm = TRUE)))
g$Rank_worst_across_variants <- apply(ranks, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
g$Stable_topN <- apply(ranks, 1, function(r) all(!is.na(r) & r <= TOP_N_GENES))
g$Previous_tier_step15 <- g$Final_Priority

keep <- c("Gene", "Ensembl", "Module", "Module_Name", "log2FC", "FDR", gfam, "Families_available", "Families_strong",
          "Evidence_score", "Evidence_rank", "Rank_best_across_variants", "Rank_worst_across_variants", "Stable_topN",
          "Evidence_tier", "Previous_tier_step15", "DEG_WGCNA_overlap")
go <- g[order(g$Evidence_rank, na.last = TRUE), keep]
write.csv(go, file.path(dir_out, paste0(ACC, "_Gene_evidence_families.csv")), row.names = FALSE)
write.csv(sens, file.path(dir_out, paste0(ACC, "_Gene_evidence_sensitivity.csv")), row.names = FALSE)
cat(">>> run", run_id, "| GENES:", nrow(g), "| families:", paste(sub("Fam_", "", gfam), collapse = ", "), "\n")
cat("    tiers:", paste(names(table(g$Evidence_tier)), table(g$Evidence_tier), collapse = ", "),
    "| stable in the top", TOP_N_GENES, "in every variant:", sum(g$Stable_topN), "\n")
print(table(Step15 = g$Final_Priority, Families = g$Evidence_tier))
cat("\n  top 15 genes:\n"); print(utils::head(go[, c("Gene", "Module", "Evidence_score", "Families_strong", "Rank_best_across_variants",
                                                       "Rank_worst_across_variants", "Evidence_tier", "Previous_tier_step15")], 15), row.names = FALSE)
cat("\n  sensitivity:\n"); print(sens, row.names = FALSE)

## ======================= DRUGS =======================
d <- rd(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_DrugPriority.csv")))
d$Fam_Reversal <- pct(-d$iLINCS_Reversal_Score)
d$Fam_Target_network <- rmean(pct(d$N_TB_Network_Targets), pct(d$Target_fraction_in_network))
d$Fam_Target_network[!(d$N_TB_Network_Targets > 0)] <- NA
gene_score <- setNames(g$Evidence_score, toupper(g$Gene))
best_target <- vapply(strsplit(ifelse(is.na(d$Known_Targets_used), "", d$Known_Targets_used), ";"), function(t) {
  s <- gene_score[toupper(t)]; if (all(is.na(s))) NA_real_ else max(s, na.rm = TRUE) }, 0)
d$Best_target_gene_score <- round(best_target, 4)
d$Fam_Target_genes <- rmean(pct(best_target), pct(d$RWR_score))
d$Fam_Module_overlap <- pct(nlog(d$Module_Enrichment_P))
d$Fam_Proximity <- pct(-d$Network_Proximity_Z)
dfam <- c("Fam_Reversal", "Fam_Target_network", "Fam_Target_genes", "Fam_Module_overlap")
score_drugs <- function(fams, strong = STRONG) {
  M <- as.matrix(d[, fams, drop = FALSE]); n_av <- rowSums(!is.na(M)); s <- rowMeans(M, na.rm = TRUE)
  ok <- n_av >= min(MIN_FAMILIES_DRUG, length(fams)) & (!"Fam_Reversal" %in% fams | !is.na(d$Fam_Reversal)); s[!ok] <- NA
  list(score = s, n_available = n_av, n_strong = rowSums(M >= strong, na.rm = TRUE))
}
dm <- score_drugs(dfam)
d$Families_available <- dm$n_available; d$Families_strong <- dm$n_strong
d$Evidence_score <- round(dm$score, 4); d$Evidence_rank <- rank(-d$Evidence_score, ties.method = "min", na.last = "keep")
d$Evidence_tier <- tier(d$Evidence_score, d$Families_strong, length(dfam))
dtop <- function(s, n) { o <- order(-s, na.last = NA); d$Drug[utils::head(o, n)] }
d_main_top <- dtop(d$Evidence_score, TOP_N_DRUGS)
dvariants <- c(list(main = dfam, with_Proximity = c(dfam, "Fam_Proximity")),
               setNames(lapply(dfam, function(f) setdiff(dfam, f)), paste0("without_", sub("Fam_", "", dfam))))
dranks <- sapply(dvariants, function(v) rank(-score_drugs(v)$score, ties.method = "min", na.last = "keep"))
dsens <- do.call(rbind, lapply(names(dvariants), function(v) { sv <- score_drugs(dvariants[[v]])
  data.frame(Variant = v, Families = paste(sub("Fam_", "", dvariants[[v]]), collapse = "+"), Drugs_scored = sum(!is.na(sv$score)),
             Spearman_with_main = round(suppressWarnings(stats::cor(sv$score, d$Evidence_score, method = "spearman", use = "pairwise")), 3),
             TopN_overlap_with_main = sum(dtop(sv$score, TOP_N_DRUGS) %in% d_main_top)) }))
d$Rank_best_across_variants  <- apply(dranks, 1, function(r) suppressWarnings(min(r, na.rm = TRUE)))
d$Rank_worst_across_variants <- apply(dranks, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
d$Stable_topN <- apply(dranks, 1, function(r) all(!is.na(r) & r <= TOP_N_DRUGS))
d$Rank_best_across_variants[!is.finite(d$Rank_best_across_variants)] <- NA; d$Rank_worst_across_variants[!is.finite(d$Rank_worst_across_variants)] <- NA
d$Previous_tier_step15 <- d$Final_Priority
dkeep <- c("Drug", "iLINCS_Reversal_Score", "N_TB_Network_Targets", "Target_fraction_in_network", "Known_Targets_used",
           "Best_target_gene_score", "Module_Enrichment_P", "Network_Proximity_Z", dfam, "Fam_Proximity", "Families_available",
           "Families_strong", "Evidence_score", "Evidence_rank", "Rank_best_across_variants", "Rank_worst_across_variants",
           "Stable_topN", "Evidence_tier", "Previous_tier_step15")
do_ <- d[order(d$Evidence_rank, na.last = TRUE), dkeep]
write.csv(do_, file.path(dir_out, paste0(ACC, "_Drug_evidence_families.csv")), row.names = FALSE)
write.csv(dsens, file.path(dir_out, paste0(ACC, "_Drug_evidence_sensitivity.csv")), row.names = FALSE)
cat("\n>>> DRUGS:", nrow(d), "| scored (reversal + at least", MIN_FAMILIES_DRUG, "families):", sum(!is.na(d$Evidence_score)),
    "| tiers:", paste(names(table(d$Evidence_tier)), table(d$Evidence_tier), collapse = ", "),
    "| stable in the top", TOP_N_DRUGS, "in every variant:", sum(d$Stable_topN), "\n")
print(table(Step15 = d$Final_Priority, Families = d$Evidence_tier))
cat("\n  top 15 drugs:\n"); print(utils::head(do_[, c("Drug", "Evidence_score", "Families_strong", "Rank_best_across_variants",
                                                       "Rank_worst_across_variants", "Evidence_tier", "Previous_tier_step15")], 15), row.names = FALSE)
cat("\n  sensitivity:\n"); print(dsens, row.names = FALSE)

writeLines(c(paste0("Evidence-family consensus for ", ACC, " (", CONTROL, " vs ", CASE, ")."),
  paste0("Gene families: ", paste(sub("Fam_", "", gfam), collapse = ", "), ". Drug families: ", paste(sub("Fam_", "", dfam), collapse = ", "),
         "; proximity is reported and tested as a sensitivity variant only."),
  "Each family score is a within-run percentile (0 = weakest, 1 = strongest); the evidence score is the mean over available families.",
  paste0("A family supports a gene or drug when its score is >= ", STRONG, ". Proposed tiers: High = support from at least ",
         round(100 * HIGH_RULE$strong_frac), "% of the families (genes ", ceiling(HIGH_RULE$strong_frac * length(gfam)), " of ", length(gfam),
         ", drugs ", ceiling(HIGH_RULE$strong_frac * length(dfam)), " of ", length(dfam), ") and score >= ", HIGH_RULE$min_score, "; Medium = at least ", MED_RULE$min_strong,
         " supporting families or score >= ", MED_RULE$min_score, "; otherwise Low. Genes need ", MIN_FAMILIES_GENE,
         " families with data, drugs need reversal evidence and ", MIN_FAMILIES_DRUG, " families."),
  paste0("Rank_best/worst_across_variants: range of the rank when each family is left out in turn (and, for drugs, when proximity is added). ",
         "Stable_topN = in the top ", TOP_N_GENES, " genes / ", TOP_N_DRUGS, " drugs in every variant.")),
  file.path(dir_out, paste0(ACC, "_EVIDENCE_FAMILIES_NOTES.txt")))
cat("\n>>> DONE. Outputs in", dir_out, "\n")
