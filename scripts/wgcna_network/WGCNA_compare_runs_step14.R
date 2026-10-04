## STEP 14: compare two runs (default: pilot GSE114192_TB vs replication GSE161829_ATB) at the level of
## DEGs, WGCNA modules, GSEA pathways and drug candidates. Base R only; reads the outputs of steps 1-13.
## Drugs are compared on UP-only iLINCS results for both runs, rebuilt from the raw concordants with the
## same rule as drugfindR::consensusConcordants (|similarity| >= cutoff, strongest signature per compound).

RUN_A  <- "GSE114192_TB"
RUN_B  <- "GSE161829_ATB"
SIM_CUTOFF <- 0.321
SIG_PADJ <- 0.05; SIG_LFC <- 1          # DEG rule used in both runs
GSEA_SETS <- c("HALLMARK", "CELLTYPE", "REACTOME")

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
source(file.path(SCRIPT_DIR, "config.R"))            # PRESETS, clean_gene_ids()
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (length(list.files(d, pattern = "^GSE[0-9]+$")) == 0 && dirname(d) != d) d <- dirname(d); d })
A <- PRESETS[[RUN_A]]; B <- PRESETS[[RUN_B]]
run_dir  <- function(p) if (nzchar(p$RUN)) file.path(ROOT_DIR, p$ACC, p$RUN) else file.path(ROOT_DIR, p$ACC)
gsea_file <- function(p, coll) {                     # dataset-level 17_gsea first, then the run folder
  fn <- paste0(p$ACC, "_", p$DEG_NAME, "_GSEA_", coll, ".csv")
  for (d in c(file.path(ROOT_DIR, p$ACC, "17_gsea"), file.path(run_dir(p), "17_gsea"))) if (file.exists(file.path(d, fn))) return(file.path(d, fn))
  file.path(ROOT_DIR, p$ACC, "17_gsea", fn)
}
out_dir  <- file.path(ROOT_DIR, "comparisons", paste0(RUN_A, "_vs_", RUN_B))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
summary_lines <- c(paste("Comparison:", RUN_A, "vs", RUN_B), paste("Date:", format(Sys.Date())), "")
note <- function(...) { line <- paste0(...); cat(line, "\n"); summary_lines <<- c(summary_lines, line) }

## ---- helpers ----
drug_key <- function(x) {                             # same key as step 10: case-insensitive, salt/hydrate words removed
  salts <- "(DIHYDROCHLORIDE|HYDROCHLORIDE|HYDROBROMIDE|HCL|MESYLATE|MESILATE|BESYLATE|BESILATE|MALEATE|FUMARATE|CITRATE|SULFATE|SULPHATE|PHOSPHATE|TARTRATE|ACETATE|SUCCINATE|TOSYLATE|TOSILATE|LACTATE|GLUCONATE|NITRATE|SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|CHLORIDE|BROMIDE|ANHYDROUS|MONOHYDRATE|DIHYDRATE|TRIHYDRATE|HEMIHYDRATE|HYDRATE)"
  x <- toupper(trimws(x))
  repeat { y <- trimws(sub(paste0("[[:space:],()-]+", salts, "$"), "", x)); if (identical(y, x)) break; x <- y }
  x
}
consensus_unpaired <- function(raw, cutoff) {         # drugfindR rule: |sim| >= cutoff, keep the strongest signature per compound
  raw <- raw[abs(raw$similarity) >= cutoff, ]
  best <- stats::ave(abs(raw$similarity), raw$treatment, FUN = max)
  raw <- raw[abs(raw$similarity) == best, ]
  raw <- raw[!duplicated(raw$treatment), ]
  data.frame(Compound = raw$treatment, Similarity = raw$similarity, stringsAsFactors = FALSE)
}
fisher_over <- function(a, b, universe) {             # one-sided overlap test inside a shared universe
  a <- intersect(a, universe); b <- intersect(b, universe)
  k <- length(intersect(a, b))
  p <- stats::fisher.test(matrix(c(k, length(a) - k, length(b) - k, length(universe) - length(union(a, b))), 2), alternative = "greater")$p.value
  list(k = k, nA = length(a), nB = length(b), N = length(universe), p = p, expected = length(a) * length(b) / length(universe))
}

## =============================================================
## 1. DEGs
## =============================================================
cat("\n========== 1. DEGs ==========\n")
read_deg <- function(p) { d <- rd(ROOT_DIR, p$ACC, "06_deg_results", paste0(p$ACC, "_", p$DEG_NAME, "_DEG.csv"))
  d$gene <- clean_gene_ids(d$gene); d <- d[order(d$padj), ]; d[!duplicated(d$gene), c("gene", "log2FoldChange", "padj")] }
da <- read_deg(A); db <- read_deg(B)
dd <- merge(da, db, by = "gene", suffixes = c("_A", "_B"))
rho <- suppressWarnings(stats::cor(dd$log2FoldChange_A, dd$log2FoldChange_B, method = "spearman", use = "complete.obs"))
sigA <- dd$gene[!is.na(dd$padj_A) & dd$padj_A < SIG_PADJ & abs(dd$log2FoldChange_A) >= SIG_LFC]
sigB <- dd$gene[!is.na(dd$padj_B) & dd$padj_B < SIG_PADJ & abs(dd$log2FoldChange_B) >= SIG_LFC]
ov <- fisher_over(sigA, sigB, dd$gene)
both <- dd[dd$gene %in% intersect(sigA, sigB), ]
same_dir <- sum(sign(both$log2FoldChange_A) == sign(both$log2FoldChange_B))
note("DEGs: ", nrow(dd), " shared genes; Spearman correlation of log2FC = ", round(rho, 2))
note("Significant (padj < ", SIG_PADJ, ", |log2FC| >= ", SIG_LFC, "): ", RUN_A, " ", length(sigA), ", ", RUN_B, " ", length(sigB),
     "; in both ", ov$k, " (expected by chance ", round(ov$expected, 1), ", Fisher p = ", signif(ov$p, 2), "); same direction: ", same_dir, "/", ov$k)
write.csv(both[order(both$padj_B), ], file.path(out_dir, "DEGs_significant_in_both.csv"), row.names = FALSE)

## =============================================================
## 2. Modules
## =============================================================
cat("\n========== 2. WGCNA modules ==========\n")
modA <- rd(run_dir(A), "10_wgcna_results", paste0(A$ACC, "_Gene_Module_Assignment.csv")); modA$Gene <- clean_gene_ids(modA$Gene)
modB <- rd(run_dir(B), "10_wgcna_results", paste0(B$ACC, "_Gene_Module_Assignment.csv")); modB$Gene <- clean_gene_ids(modB$Gene)
selA <- with(rd(run_dir(A), "10_wgcna_results", paste0(A$ACC, "_ModuleTrait_Correlation.csv")), Module[Selected %in% TRUE])
selB <- with(rd(run_dir(B), "10_wgcna_results", paste0(B$ACC, "_ModuleTrait_Correlation.csv")), Module[Selected %in% TRUE])
U <- intersect(modA$Gene, modB$Gene)
rows <- list()
for (ma in selA) for (mb in setdiff(unique(modB$ModuleColor), "grey")) {
  f <- fisher_over(modA$Gene[modA$ModuleColor == ma], modB$Gene[modB$ModuleColor == mb], U)
  if (f$k == 0) next
  shared <- intersect(modA$Gene[modA$ModuleColor == ma], modB$Gene[modB$ModuleColor == mb])
  sd <- dd[dd$gene %in% shared, ]
  rows[[length(rows) + 1]] <- data.frame(Module_A = ma, Module_B = mb, Selected_B = mb %in% selB, Size_A = f$nA, Size_B = f$nB,
    Overlap = f$k, Expected = round(f$expected, 1), Jaccard = round(f$k / (f$nA + f$nB - f$k), 3), P_fisher = f$p,
    Same_direction_pct = if (nrow(sd)) round(100 * mean(sign(sd$log2FoldChange_A) == sign(sd$log2FoldChange_B)), 1) else NA)
}
mo <- do.call(rbind, rows); mo$FDR <- stats::p.adjust(mo$P_fisher, "BH"); mo <- mo[order(mo$Module_A, mo$P_fisher), ]
write.csv(mo, file.path(out_dir, "module_overlap.csv"), row.names = FALSE)
note("Modules (universe = ", length(U), " genes in both networks). Best match in ", RUN_B, " for each selected ", RUN_A, " module:")
for (ma in selA) { x <- mo[mo$Module_A == ma, ][1, ]
  note("  ", ma, " (", x$Size_A, " genes) -> ", x$Module_B, if (x$Selected_B) " [selected]" else "", ": ", x$Overlap, " shared (expected ",
       x$Expected, ", FDR ", signif(x$FDR, 2), "), same direction ", x$Same_direction_pct, "%") }

## =============================================================
## 3. Pathways (GSEA)
## =============================================================
cat("\n========== 3. GSEA ==========\n")
for (coll in GSEA_SETS) {
  fa <- gsea_file(A, coll); fb <- gsea_file(B, coll)
  if (!file.exists(fa) || !file.exists(fb)) { note("GSEA ", coll, ": file missing, skipped"); next }
  ga <- read.csv(fa)[, c("ID", "NES", "p.adjust")]; gb <- read.csv(fb)[, c("ID", "NES", "p.adjust")]
  g <- merge(ga, gb, by = "ID", suffixes = c("_A", "_B"))
  g$Status <- with(g, ifelse(p.adjust_A < 0.05 & p.adjust_B < 0.05, ifelse(sign(NES_A) == sign(NES_B), "replicated", "opposite"),
                      ifelse(p.adjust_A < 0.05, paste("only", RUN_A), ifelse(p.adjust_B < 0.05, paste("only", RUN_B), "neither"))))
  g <- g[order(g$Status != "replicated", -abs(g$NES_A + g$NES_B)), ]
  write.csv(g, file.path(out_dir, paste0("GSEA_", coll, "_comparison.csv")), row.names = FALSE)
  rep <- g[g$Status == "replicated", ]
  note("GSEA ", coll, ": ", nrow(g), " sets in both; NES correlation (Spearman) ", round(stats::cor(g$NES_A, g$NES_B, method = "spearman"), 2),
       "; replicated (FDR < 0.05 in both, same sign) ", nrow(rep), ", opposite ", sum(g$Status == "opposite"))
  if (coll == "HALLMARK" && nrow(rep)) note("  replicated Hallmark: ", paste0(sub("^HALLMARK_", "", rep$ID), " (", sprintf("%+.2f/%+.2f", rep$NES_A, rep$NES_B), ")", collapse = ", "))
}

## =============================================================
## 4. Drugs
## =============================================================
cat("\n========== 4. Drugs ==========\n")
rawA <- rd(run_dir(A), "13_ilincs_signature", paste0(A$ACC, "_concordants_UP_raw.csv"))
rawB <- rd(run_dir(B), "13_ilincs_signature", paste0(B$ACC, "_concordants_UP_raw.csv"))
ca <- consensus_unpaired(rawA, SIM_CUTOFF); cb <- consensus_unpaired(rawB, SIM_CUTOFF)
ca$Key <- drug_key(ca$Compound); cb$Key <- drug_key(cb$Compound)
## self-check: the rebuilt UP-only result for B must equal step 7's output when B was run UP-only
modeB <- file.path(run_dir(B), "13_ilincs_signature", paste0(B$ACC, "_consensus_mode.txt"))
if (file.exists(modeB) && grepl("UP only", readLines(modeB, warn = FALSE)[1])) {
  candB <- rd(run_dir(B), "13_ilincs_signature", paste0(B$ACC, "_iLINCS_candidate_compounds.csv"))
  ok <- setequal(candB$Target, cb$Compound[cb$Similarity < 0])
  note("Self-check: rebuilt UP-only reversal list for ", RUN_B, " matches step 7 output: ", ok)
}
universe <- intersect(unique(drug_key(rawA$treatment)), unique(drug_key(rawB$treatment)))
revA <- unique(ca$Key[ca$Similarity < 0]); revB <- unique(cb$Key[cb$Similarity < 0])
f <- fisher_over(revA, revB, universe)
note("Drugs (UP-only signatures, compounds tested in both = ", length(universe), "): reversal in ", RUN_A, " ", f$nA, ", in ", RUN_B, " ", f$nB,
     "; in both ", f$k, " (expected ", round(f$expected), ", Fisher p = ", signif(f$p, 2), ")")
## agreement of the drug rankings themselves (not only the large reversal lists)
ka <- ca[!duplicated(ca$Key), ]; kb <- cb[!duplicated(cb$Key), ]
kk <- intersect(ka$Key, kb$Key)
sa_all <- setNames(ka$Similarity, ka$Key)[kk]; sb_all <- setNames(kb$Similarity, kb$Key)[kk]
note("  rank agreement over ", length(kk), " compounds scored in both: Spearman of similarity ", round(stats::cor(sa_all, sb_all, method = "spearman"), 2))
for (N in c(50, 100)) { ta <- names(sort(sa_all))[1:N]; tb <- names(sort(sb_all))[1:N]; o <- length(intersect(ta, tb))
  note("  top-", N, " strongest reversers: ", o, " shared (expected ", round(N * N / length(kk), 1), ", p = ",
       signif(stats::phyper(o - 1, N, length(kk) - N, N, lower.tail = FALSE), 2), ")") }
ms <- merge(rawA[, c("signatureid", "similarity")], rawB[, c("signatureid", "similarity")], by = "signatureid", suffixes = c("_A", "_B"))
note("  same LINCS signatures scored in both: ", nrow(ms), "; Spearman of similarity ", round(stats::cor(ms$similarity_A, ms$similarity_B, method = "spearman"), 2))
upA <- rd(run_dir(A), "13_ilincs_signature", paste0(A$ACC, "_TB_UP_signature.csv"))$Name_GeneSymbol
upB <- rd(run_dir(B), "13_ilincs_signature", paste0(B$ACC, "_TB_UP_signature.csv"))$Name_GeneSymbol
note("  TB_UP signature genes: ", length(upA), " vs ", length(upB), "; shared: ", paste(sort(intersect(upA, upB)), collapse = ", "))
sa <- stats::aggregate(Similarity ~ Key, ca[ca$Similarity < 0, ], min); sb <- stats::aggregate(Similarity ~ Key, cb[cb$Similarity < 0, ], min)
rep_drugs <- merge(sa, sb, by = "Key", suffixes = c("_A", "_B"))
rep_drugs$Mean_similarity <- (rep_drugs$Similarity_A + rep_drugs$Similarity_B) / 2
rep_drugs <- rep_drugs[order(rep_drugs$Mean_similarity), ]
write.csv(rep_drugs, file.path(out_dir, "drugs_reversal_in_both_UPonly.csv"), row.names = FALSE)

## final-table level: drugs that also have targets in BOTH TB networks
fa <- rd(run_dir(A), "16_final_tables", paste0(A$ACC, "_FINAL_DrugTable.csv")); fb <- rd(run_dir(B), "16_final_tables", paste0(B$ACC, "_FINAL_DrugTable.csv"))
fa$Key <- drug_key(fa$Drug); fb$Key <- drug_key(fb$Drug)
keepA <- intersect(c("Key", "Drug", "iLINCS_Reversal_Score", "Known_Targets", "Modules_Affected", "Enrichr_KEGG", "Enrichr_Reactome", "SC_Top_Immune_CellType"), names(fa))
keepB <- intersect(c("Key", "iLINCS_Reversal_Score", "Known_Targets", "Modules_Affected", "Enrichr_KEGG"), names(fb))
ft <- merge(fa[, keepA], fb[, keepB], by = "Key", suffixes = c("_A", "_B"))
ft$Reversal_UPonly_both <- ft$Key %in% rep_drugs$Key
ft <- ft[order(!ft$Reversal_UPonly_both, ft$iLINCS_Reversal_Score_A + ft$iLINCS_Reversal_Score_B), ]
write.csv(ft, file.path(out_dir, "drugs_in_both_final_tables.csv"), row.names = FALSE)
note("Drugs in both final drug tables (target in both TB networks): ", nrow(ft), "; of these reversing in both (UP-only): ", sum(ft$Reversal_UPonly_both))

writeLines(summary_lines, file.path(out_dir, "SUMMARY.txt"))

## ---- README.md files, so GitHub shows a readable page for each comparison and an index ----
files_md <- c(
  "| File | Content |", "|---|---|",
  "| `SUMMARY.txt` | all numbers above, plain text |",
  "| `DEGs_significant_in_both.csv` | genes significant in both runs, with log2FC and FDR in each |",
  "| `module_overlap.csv` | overlap of every selected module of run A with every module of run B |",
  "| `GSEA_HALLMARK_comparison.csv`, `GSEA_CELLTYPE_comparison.csv`, `GSEA_REACTOME_comparison.csv` | NES and FDR of each gene set in both runs, with a replicated / opposite / only-one label |",
  "| `drugs_reversal_in_both_UPonly.csv` | compounds that reverse the TB_UP signature in both runs |",
  "| `drugs_in_both_final_tables.csv` | drugs present in both final drug tables (targets in both TB networks) |")
body <- sub("^  ", "    ", summary_lines[-(1:3)])
md <- c(paste0("# ", RUN_A, " vs ", RUN_B), "",
        paste0("Generated by `scripts/wgcna_network/WGCNA_compare_runs_step14.R` on ", format(Sys.Date()), "."), "",
        "## Results", "", "```", body, "```", "", "## Files", "", files_md)
writeLines(md, file.path(out_dir, "README.md"))
idx_dirs <- sort(list.dirs(file.path(ROOT_DIR, "comparisons"), recursive = FALSE, full.names = FALSE))
idx <- c("# Comparisons between runs", "",
         "Each folder compares two runs of the network pipeline (DEGs, WGCNA modules, GSEA pathways, drug candidates).",
         "Open a folder to see its summary.", "",
         paste0("- [", idx_dirs, "](", idx_dirs, "/)"))
writeLines(idx, file.path(ROOT_DIR, "comparisons", "README.md"))
cat("\n>>> DONE. Outputs in", out_dir, "\n")
