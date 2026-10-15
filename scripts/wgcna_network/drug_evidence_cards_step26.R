## STEP 26: evidence cards for candidate drugs.
##
## The drug ranking is less stable than the gene ranking, so instead of a fixed shortlist each candidate gets a card
## with the evidence that matters for interpreting it:
##   1. reversal of BOTH halves of the disease signature: best (most negative) iLINCS similarity against the
##      up-regulated genes and against the down-regulated genes
##   2. consistency across LINCS experiments: share of the compound's signatures (cell lines, doses, times) that
##      reverse the signature, and the number of cell lines in which it does
##   3. mechanism: targets in the disease network, which of them are High-priority genes, which modules they sit in,
##      and whether the targets point to cell-division machinery (antiproliferative drugs reverse proliferation
##      signatures regardless of TB)
##   4. stability: evidence score, rank range across evidence-family variants (step 25)
##   5. network proximity, as supporting information only (no drug passes FDR)
## Candidates = drugs stable in the top 20 in every variant, or High in step 25, or in the top CARD_TOP by score.
## Run after step 25:  Rscript drug_evidence_cards_step26.R     (TB_RUN=... for other runs)
## Outputs in <run>/26_drug_evidence_cards/.

## ---- CONFIG ----
CARD_TOP  <- 25        # also card the best 25 by evidence score
STRONG_SIM <- -0.321   # same strength threshold as the iLINCS consensus step
CELLCYCLE_TARGET <- "^(CDK[0-9]+|PLK[0-9]|AURK[ABC]|TOP1|TOP2A|TOP2B|TYMS|RRM[12]|DHFR|KIF11|WEE1|CHEK[12]|TTK|BUB1B?|CDC25[ABC]|HDAC[1-3])$"

get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
source(file.path(SCRIPT_DIR, "drug_name_key.R"))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "26_drug_evidence_cards"); dir.create(dir_out, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) NULL else read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

d <- rd(OUT_DIR, "20_evidence_families", paste0(ACC, "_Drug_evidence_families.csv"))
g <- rd(OUT_DIR, "20_evidence_families", paste0(ACC, "_Gene_evidence_families.csv"))
if (is.null(d) || is.null(g)) stop("Run step 25 (evidence_family_consensus_step25.R) first.")
up <- rd(OUT_DIR, "13_ilincs_signature", paste0(ACC, "_concordants_UP_raw.csv"))
dn <- rd(OUT_DIR, "13_ilincs_signature", paste0(ACC, "_concordants_DOWN_raw.csv"))
if (is.null(up)) stop("No iLINCS UP results for this run.")
prox <- rd(OUT_DIR, "15_network_algorithms", paste0(ACC, "_NetworkProximity_withFDR.csv"))
## all known targets (not only those inside the disease network) -- used to recognise antiproliferative mechanisms
all_t <- rd(OUT_DIR, "14_drug_targets", paste0(ACC, "_drug_targets_ALL.csv"))
all_by_key <- if (is.null(all_t)) list() else split(all_t$Target_Gene, drug_key(all_t$Drug))

## ---- candidates ----
cand <- d[(d$Stable_topN %in% TRUE) | d$Evidence_tier == "High" | (!is.na(d$Evidence_rank) & d$Evidence_rank <= CARD_TOP), ]
cand <- cand[order(cand$Evidence_rank), ]
cat(">>> run", run_id, "|", nrow(cand), "candidate drugs (stable top 20, High, or top", CARD_TOP, "by evidence score)\n")

## ---- iLINCS per direction ----
per_dir <- function(raw, key) {
  if (is.null(raw)) return(c(best = NA, n_sig = 0, frac_rev = NA, n_lines_rev = 0))
  r <- raw[drug_key(raw$treatment) == key, ]
  if (!nrow(r)) return(c(best = NA, n_sig = 0, frac_rev = NA, n_lines_rev = 0))
  c(best = min(r$similarity), n_sig = nrow(r), frac_rev = mean(r$similarity < 0),
    n_lines_rev = length(unique(r$cellline[r$similarity <= STRONG_SIM])))
}
gene_tier <- setNames(g$Evidence_tier, toupper(g$Gene)); gene_mod <- setNames(g$Module, toupper(g$Gene))
cards <- do.call(rbind, lapply(seq_len(nrow(cand)), function(i) {
  x <- cand[i, ]; k <- drug_key(x$Drug)
  u <- per_dir(up, k); w <- per_dir(dn, k)
  tg <- if (is.na(x$Known_Targets_used) || !nzchar(x$Known_Targets_used)) character(0) else strsplit(x$Known_Targets_used, ";")[[1]]
  hi <- tg[gene_tier[toupper(tg)] %in% "High"]
  mods <- sort(unique(stats::na.omit(gene_mod[toupper(tg)])))
  at <- unique(c(tg, all_by_key[[k]])); cc <- sort(unique(at[grepl(CELLCYCLE_TARGET, at)]))
  pz <- if (!is.null(prox) && "Drug" %in% names(prox)) prox[prox$Drug == x$Drug, ] else NULL
  data.frame(Drug = x$Drug, Evidence_rank = x$Evidence_rank, Evidence_score = x$Evidence_score, Evidence_tier = x$Evidence_tier,
    Stable_top20 = x$Stable_topN, Rank_range = paste(x$Rank_best_across_variants, x$Rank_worst_across_variants, sep = "-"),
    Best_sim_UP = round(u[["best"]], 3), Best_sim_DOWN = round(w[["best"]], 3),
    Reverses_both = !is.na(u[["best"]]) && !is.na(w[["best"]]) && u[["best"]] <= STRONG_SIM && w[["best"]] <= STRONG_SIM,
    LINCS_signatures_UP = u[["n_sig"]], Share_reversing_UP = round(u[["frac_rev"]], 2), Cell_lines_reversing_UP = u[["n_lines_rev"]],
    LINCS_signatures_DOWN = w[["n_sig"]], Share_reversing_DOWN = round(w[["frac_rev"]], 2), Cell_lines_reversing_DOWN = w[["n_lines_rev"]],
    Targets_in_network = length(tg), Targets = paste(tg, collapse = ";"), High_priority_targets = paste(hi, collapse = ";"),
    Target_modules = paste(mods, collapse = ";"), Cell_division_targets = paste(cc, collapse = ";"),
    Proximity_Z = if (!is.null(pz) && nrow(pz)) round(pz$Z_score[1], 2) else round(x$Network_Proximity_Z, 2),
    Proximity_FDR = if (!is.null(pz) && nrow(pz) && "FDR_BH" %in% names(pz)) signif(pz$FDR_BH[1], 2) else NA,
    stringsAsFactors = FALSE)
}))
## a short reading per drug
cards$Reading <- with(cards, paste0(
  ifelse(Reverses_both, "reverses both halves of the signature", ifelse(!is.na(Best_sim_UP) & Best_sim_UP <= STRONG_SIM, "reverses the up-regulated genes only", "weak reversal")),
  ifelse(!is.na(Share_reversing_UP) & Share_reversing_UP >= 0.6, paste0(", consistently (", round(100 * Share_reversing_UP), "% of signatures, ", Cell_lines_reversing_UP, " cell lines)"),
         ifelse(!is.na(Share_reversing_UP), paste0(", in ", round(100 * Share_reversing_UP), "% of signatures"), "")),
  ifelse(nzchar(High_priority_targets), paste0("; hits High-priority genes ", gsub(";", ", ", High_priority_targets)), ""),
  ifelse(nzchar(Cell_division_targets), paste0("; cell-division targets (", gsub(";", ", ", Cell_division_targets), ")"), ""),
  ifelse(Stable_top20 %in% TRUE, "; stable across evidence variants", "")))
write.csv(cards, file.path(dir_out, paste0(ACC, "_drug_evidence_cards.csv")), row.names = FALSE)

## ---- readable version ----
md <- c(paste0("# Drug evidence cards: ", ACC, " (", CASE, " vs ", CONTROL, ")"), "",
  "Candidates are drugs that stay in the top 20 in every evidence-family variant, are High in the evidence-family ranking, or are in the",
  paste0("top ", CARD_TOP, " by evidence score. Similarities are iLINCS concordance scores (negative = reversal; strong at <= ", STRONG_SIM, ")."),
  "Network proximity is shown as supporting information only: no drug passes FDR.", "",
  "| Drug | Rank (range) | Tier | Reverses up / down | Share of signatures reversing (up) | High-priority targets | Cell-division targets | Reading |",
  "| --- | --- | --- | --- | --- | --- | --- | --- |",
  sprintf("| %s | %s (%s) | %s%s | %s / %s | %s | %s | %s | %s |", cards$Drug, cards$Evidence_rank, cards$Rank_range, cards$Evidence_tier,
          ifelse(cards$Stable_top20 %in% TRUE, ", stable", ""), cards$Best_sim_UP, ifelse(is.na(cards$Best_sim_DOWN), "not queried", cards$Best_sim_DOWN),
          ifelse(is.na(cards$Share_reversing_UP), "", paste0(round(100 * cards$Share_reversing_UP), "%")),
          gsub(";", ", ", cards$High_priority_targets), gsub(";", ", ", cards$Cell_division_targets), cards$Reading))
writeLines(md, file.path(dir_out, paste0(ACC, "_drug_evidence_cards.md")))

cat("    reverse both halves:", sum(cards$Reverses_both), "| consistent (>= 60% of signatures):", sum(cards$Share_reversing_UP >= 0.6, na.rm = TRUE),
    "| hit High-priority genes:", sum(nzchar(cards$High_priority_targets)), "| cell-division targets:", sum(nzchar(cards$Cell_division_targets)), "\n\n")
best <- cards[cards$Reverses_both & cards$Share_reversing_UP >= 0.6 & nzchar(cards$High_priority_targets), ]
cat(">>> strongest-supported (both halves reversed, consistent across signatures, hits High-priority genes):",
    if (nrow(best)) paste0(best$Drug, ifelse(nzchar(best$Cell_division_targets), " [has cell-division targets]", ""), collapse = ", ") else "none", "\n")
cat("    (cell-division targets are counted among ALL known targets; multi-kinase inhibitors list many kinases from screening panels,\n",
    "     so the flag marks a possible antiproliferative mechanism, not a proven one)\n\n", sep = "")
print(cards[, c("Drug", "Evidence_rank", "Rank_range", "Best_sim_UP", "Best_sim_DOWN", "Share_reversing_UP", "High_priority_targets", "Cell_division_targets")], row.names = FALSE)
cat("\n>>> DONE. Outputs in", dir_out, "\n")
