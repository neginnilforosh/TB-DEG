## STEP 21: external validation of the prioritised genes (Ratul, 10 Oct: "use the other active-TB/LTBI
## datasets mainly as an external validation and phenotype-characterisation layer ... test whether the
## final prioritised genes are reproduced in independent active-TB datasets, also altered in LTBI, more
## characteristic of active TB than latent TB. Do not make simple cross-dataset DEG overlap the main analysis.")
##
## For the active run's consensus genes (step 15), looks each gene up in every external comparison and asks:
##   1. direction:  does it move the same way as in the discovery run?   (sign test against 50%)
##   2. reproduced: same direction AND padj < 0.05 there?                (High vs Low tier, Fisher)
##   3. active vs latent: within a dataset that has both, is the change bigger in active TB than in LTBI?
## Outputs in <run>/23_external_validation/. Base R + ggplot2.

## ---- CONFIG ----
## comparisons: dataset, DEG file stem, phenotype, note. log2FC sign must mean "higher in disease".
EXTERNAL <- list(
  list(ds = "GSE161829", stem = "TBneg_vs_ATB",       pheno = "Active", note = "whole blood, n = 20 vs 25"),
  list(ds = "GSE222001", stem = "Healthy_vs_Active",  pheno = "Active", note = "n = 5 vs 5 -- low power"),
  list(ds = "GSE161829", stem = "TBneg_vs_LTBI",      pheno = "LTBI",   note = "whole blood, n = 20 vs 40"),
  list(ds = "GSE222001", stem = "Healthy_vs_Latent",  pheno = "LTBI",   note = "n = 5 vs 5 -- low power"),
  list(ds = "GSE99374",  stem = "HC_vs_LTBI",         pheno = "LTBI",   note = "sorted CD8 T cells, not whole blood"))
PADJ <- 0.05
TIERS <- c("High", "Medium", "Low")

suppressPackageStartupMessages(library(ggplot2))
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "23_external_validation"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
FIG_DPI <- 600
save_gg <- function(f, p, w, h) { ggsave(f, p, width = w, height = h, dpi = FIG_DPI); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h) }

pri <- read.csv(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")), stringsAsFactors = FALSE)
pri <- pri[!is.na(pri$log2FC), ]
pri$Disc_sign <- sign(pri$log2FC)
cat(">>> discovery run:", run_id, "|", nrow(pri), "genes |", paste(names(table(pri$Final_Priority)), table(pri$Final_Priority), collapse = ", "), "\n")

## one key per external row: Ensembl without version when there is one, otherwise the symbol
key_of <- function(id) {
  ens <- ifelse(grepl("^ENSG", id), sub("\\..*$", "", sub("\\|.*$", "", id)), NA)
  sym <- ifelse(grepl("\\|", id), sub("^.*\\|", "", id), ifelse(grepl("^ENSG", id), NA, id))
  data.frame(ens = ens, sym = sym, stringsAsFactors = FALSE)
}
long <- list()
for (e in EXTERNAL) {
  f <- file.path(ROOT_DIR, e$ds, "06_deg_results", paste0(e$ds, "_", e$stem, "_DEG.csv"))
  if (!file.exists(f)) { cat("!! missing, skipped:", f, "\n"); next }
  x <- read.csv(f, stringsAsFactors = FALSE); k <- key_of(x$gene); x <- cbind(x, k)
  by_ens <- x[!is.na(x$ens), ]; by_ens <- by_ens[!duplicated(by_ens$ens), ]
  by_sym <- x[!is.na(x$sym), ];  by_sym <- by_sym[order(by_sym$pvalue), ]; by_sym <- by_sym[!duplicated(by_sym$sym), ]
  m <- match(pri$Ensembl, by_ens$ens); use <- by_ens[m, ]
  miss <- is.na(m) & !is.na(pri$Gene); ms <- match(pri$Gene[miss], by_sym$sym)
  use[miss, ] <- by_sym[ms, ]
  comp <- paste(e$ds, e$stem, sep = ":")
  long[[comp]] <- data.frame(Gene = pri$Gene, Ensembl = pri$Ensembl, Module = pri$Module,
    Final_Priority = pri$Final_Priority, Discovery_log2FC = pri$log2FC, Comparison = comp, Phenotype = e$pheno,
    Ext_log2FC = use$log2FoldChange, Ext_padj = use$padj, stringsAsFactors = FALSE)
  cat(sprintf("  %-32s %-6s matched %5d / %d genes  (%s)\n", comp, e$pheno, sum(!is.na(use$log2FoldChange)), nrow(pri), e$note))
}
L <- do.call(rbind, long)
L$Same_direction <- ifelse(is.na(L$Ext_log2FC), NA, sign(L$Ext_log2FC) == sign(L$Discovery_log2FC))
L$Reproduced <- !is.na(L$Ext_padj) & L$Ext_padj < PADJ & L$Same_direction %in% TRUE
write.csv(L, file.path(dir_out, paste0(ACC, "_external_validation_long.csv")), row.names = FALSE)

## ---- 1 + 2: per comparison and tier ----
summ <- do.call(rbind, lapply(split(L, L$Comparison), function(d) do.call(rbind, lapply(TIERS, function(t) {
  s <- d[d$Final_Priority == t & !is.na(d$Same_direction), ]
  if (!nrow(s)) return(NULL)
  bt <- stats::binom.test(sum(s$Same_direction), nrow(s), 0.5, alternative = "greater")
  data.frame(Comparison = d$Comparison[1], Phenotype = d$Phenotype[1], Tier = t, N_tested = nrow(s),
             Pct_same_direction = round(100 * mean(s$Same_direction), 1), P_direction = signif(bt$p.value, 3),
             N_reproduced = sum(s$Reproduced), Pct_reproduced = round(100 * mean(s$Reproduced), 1))
}))))
## High vs Low reproduction (is the consensus tier informative, beyond "these genes are TB genes"?)
summ$P_High_vs_Low_reproduced <- NA
for (comp in unique(summ$Comparison)) {
  d <- L[L$Comparison == comp & !is.na(L$Same_direction), ]
  h <- d$Reproduced[d$Final_Priority == "High"]; l <- d$Reproduced[d$Final_Priority == "Low"]
  if (length(h) && length(l)) {
    p <- stats::fisher.test(matrix(c(sum(h), sum(!h), sum(l), sum(!l)), 2), alternative = "greater")$p.value
    summ$P_High_vs_Low_reproduced[summ$Comparison == comp & summ$Tier == "High"] <- signif(p, 3)
  }
}
write.csv(summ, file.path(dir_out, paste0(ACC, "_external_validation_summary.csv")), row.names = FALSE)
cat("\n========== reproduction of the prioritised genes ==========\n"); print(summ, row.names = FALSE)

## ---- 3: active vs latent within the same dataset ----
al <- list()
for (ds in unique(sapply(EXTERNAL, `[[`, "ds"))) {
  a <- L[grepl(paste0("^", ds, ":"), L$Comparison) & L$Phenotype == "Active", ]
  b <- L[grepl(paste0("^", ds, ":"), L$Comparison) & L$Phenotype == "LTBI", ]
  if (!nrow(a) || !nrow(b)) next
  m <- merge(a[, c("Ensembl", "Gene", "Final_Priority", "Discovery_log2FC", "Ext_log2FC", "Reproduced")],
             b[, c("Ensembl", "Ext_log2FC", "Ext_padj")], by = "Ensembl", suffixes = c("_active", "_LTBI"))
  m <- m[!is.na(m$Ext_log2FC_active) & !is.na(m$Ext_log2FC_LTBI), ]
  ## signed in the discovery direction, so "bigger" means a bigger change the TB way
  m$Shift_active <- m$Ext_log2FC_active * sign(m$Discovery_log2FC)
  m$Shift_LTBI   <- m$Ext_log2FC_LTBI   * sign(m$Discovery_log2FC)
  m$Active_specific <- m$Reproduced & !(m$Ext_padj < PADJ & sign(m$Ext_log2FC_LTBI) == sign(m$Discovery_log2FC))
  m$Dataset <- ds
  for (t in TIERS) {
    s <- m[m$Final_Priority == t, ]; if (nrow(s) < 3) next
    w <- stats::wilcox.test(s$Shift_active, s$Shift_LTBI, paired = TRUE, alternative = "greater")
    al[[paste(ds, t)]] <- data.frame(Dataset = ds, Tier = t, N = nrow(s),
      Median_shift_active = round(stats::median(s$Shift_active), 3), Median_shift_LTBI = round(stats::median(s$Shift_LTBI), 3),
      Pct_bigger_in_active = round(100 * mean(s$Shift_active > s$Shift_LTBI), 1), P_paired = signif(w$p.value, 3),
      N_active_specific = sum(s$Active_specific, na.rm = TRUE))
  }
  write.csv(m, file.path(dir_out, paste0(ACC, "_active_vs_latent_", ds, ".csv")), row.names = FALSE)
}
if (length(al)) {
  al <- do.call(rbind, al)
  write.csv(al, file.path(dir_out, paste0(ACC, "_active_vs_latent_summary.csv")), row.names = FALSE)
  cat("\n========== active TB vs LTBI (same dataset, change in the discovery direction) ==========\n"); print(al, row.names = FALSE)
}

## ---- figure: High genes across comparisons ----
H <- L[L$Final_Priority == "High" & !is.na(L$Gene), ]
if (nrow(H)) {
  H$Comparison <- factor(H$Comparison, levels = unique(L$Comparison))
  ord <- H[H$Comparison == levels(H$Comparison)[1], ]; ord <- ord$Gene[order(ord$Discovery_log2FC)]
  H$Gene <- factor(H$Gene, levels = unique(c(ord, H$Gene)))
  H$lab <- ifelse(H$Reproduced, "*", "")
  lim <- max(abs(H$Ext_log2FC), na.rm = TRUE)
  p <- ggplot(H, aes(Comparison, Gene, fill = Ext_log2FC)) + geom_tile(colour = "white") +
    geom_text(aes(label = lab), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim), na.value = "grey85", name = "log2FC") +
    labs(title = paste0(ACC, " High-priority genes in external datasets"), subtitle = "* = same direction and padj < 0.05",
         x = NULL, y = NULL) + theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_gg(file.path(dir_out, paste0(ACC, "_external_validation_heatmap.png")), p, 6.5, max(4, 0.16 * nlevels(H$Gene) + 1.5))
}
cat("\n>>> DONE. Outputs in", dir_out, "\n")
