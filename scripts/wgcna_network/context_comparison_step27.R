## STEP 27: compare the clinical contexts once each has been analysed on its own.
##
## Each context (TB vs healthy, TB with diabetes vs diabetes, later TB with HIV) has its own differential expression,
## modules, network and evidence-family ranking. This step lines them up gene by gene to answer:
##   1. which genes respond to TB in every context (core) and which only in one (context-specific),
##      keeping the direction of change
##   2. whether a gene keeps its place in the network across contexts (module, centrality, evidence rank) or is
##      "rewired": prioritised in one context but not in another, or moved between modules
##   3. which Hallmark pathways are shared, and which drugs are supported in more than one context
## Add a context by adding its run id to CONTEXTS; nothing else changes.
## Run from scripts/wgcna_network after step 25 (and step 11 GSEA) for every context:
##     Rscript context_comparison_step27.R
## Outputs in comparisons/contexts/.

## ---- CONFIG ----
CONTEXTS <- c(TB = "GSE114192_TB", TB_DM = "GSE114192_TBDM_vs_DM")   # add TB_HIV = "<run id>" when available
PADJ <- 0.05; LFC <- 1
TOP_RANK <- 50          # "prioritised" = in the top 50 by evidence score
REWIRE_GAP <- 200       # evidence rank difference that counts as a large move

get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
Sys.setenv(TB_RUN = CONTEXTS[[1]])
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, "GSE114192", "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT_DIR, "comparisons", "contexts"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL

## ---- load each context ----
ctx <- lapply(names(CONTEXTS), function(nm) {
  p <- PRESETS[[CONTEXTS[[nm]]]]; base <- file.path(ROOT_DIR, p$ACC); run <- if (nzchar(p$RUN)) file.path(base, p$RUN) else base
  deg_dir <- if (dir.exists(file.path(run, "06_deg_results"))) file.path(run, "06_deg_results") else file.path(base, "06_deg_results")
  deg <- rd(file.path(deg_dir, paste0(p$ACC, "_", p$DEG_NAME, "_DEG.csv")))
  ev  <- rd(file.path(run, "20_evidence_families", paste0(p$ACC, "_Gene_evidence_families.csv")))
  gp  <- rd(file.path(run, "19_consensus", paste0(p$ACC, "_GenePriority_all_modules.csv")))
  dr  <- rd(file.path(run, "20_evidence_families", paste0(p$ACC, "_Drug_evidence_families.csv")))
  gs  <- rd(file.path(base, "17_gsea", paste0(p$ACC, "_", p$DEG_NAME, "_GSEA_HALLMARK.csv")))
  am  <- rd(file.path(run, "10_wgcna_results", paste0(p$ACC, "_Gene_Module_Assignment.csv")))
  miss <- c(deg = is.null(deg), evidence = is.null(ev), consensus = is.null(gp))
  if (any(miss)) stop("Context ", nm, ": missing ", paste(names(miss)[miss], collapse = ", "), " (run its pipeline through step 25).")
  list(name = nm, label = paste(p$CASE, "vs", p$CONTROL), deg = deg, ev = ev, gp = gp, dr = dr, gs = gs, am = am, mods = p$MODULE_NAMES)
}); names(ctx) <- names(CONTEXTS)
cat(">>> contexts:", paste(sprintf("%s (%s)", names(ctx), vapply(ctx, `[[`, "", "label")), collapse = "; "), "\n")

## ---- 1. differential expression across contexts ----
key <- Reduce(intersect, lapply(ctx, function(x) x$deg$gene[!is.na(x$deg$padj)]))
D <- data.frame(Ensembl = key)
for (nm in names(ctx)) { x <- ctx[[nm]]$deg; m <- match(key, x$gene)
  D[[paste0("log2FC_", nm)]] <- x$log2FoldChange[m]; D[[paste0("padj_", nm)]] <- x$padj[m]
  D[[paste0("DEG_", nm)]] <- !is.na(x$padj[m]) & x$padj[m] < PADJ & abs(x$log2FoldChange[m]) >= LFC }
degc <- paste0("DEG_", names(ctx)); lfcc <- paste0("log2FC_", names(ctx))
D$N_contexts_DEG <- rowSums(D[, degc, drop = FALSE])
D$Same_direction <- apply(sign(D[, lfcc, drop = FALSE]), 1, function(s) length(unique(s)) == 1)
D$DEG_class <- ifelse(D$N_contexts_DEG == length(ctx) & D$Same_direction, "core (all contexts, same direction)",
               ifelse(D$N_contexts_DEG == length(ctx), "all contexts, opposite direction",
               ifelse(D$N_contexts_DEG >= 1, paste("specific:", apply(D[, degc, drop = FALSE], 1, function(r) paste(names(ctx)[r], collapse = "+"))), "not DE")))
cat("\n========== 1. differential expression ==========\n")
for (nm in names(ctx)) cat(sprintf("  %-6s %4d DEGs\n", nm, sum(D[[paste0("DEG_", nm)]])))
print(table(D$DEG_class[D$DEG_class != "not DE"]))
if (length(ctx) == 2) {
  a <- D[[degc[1]]]; b <- D[[degc[2]]]; ft <- stats::fisher.test(table(a, b), alternative = "greater")
  cat(sprintf("  overlap %d (expected %.1f, Fisher p %.1e) | log2FC Spearman over all %d genes %.2f\n", sum(a & b),
              sum(a) * sum(b) / nrow(D), ft$p.value, nrow(D), stats::cor(D[[lfcc[1]]], D[[lfcc[2]]], method = "spearman")))
}

## ---- 2. prioritised genes and network rewiring ----
G <- NULL
for (nm in names(ctx)) {
  e <- ctx[[nm]]$ev; g <- ctx[[nm]]$gp
  x <- merge(e[, c("Ensembl", "Gene", "Module", "Module_Name", "Evidence_score", "Evidence_rank", "Evidence_tier", "Stable_topN")],
             g[, intersect(c("Ensembl", "MM", "PPI_Degree", "PPI_Betweenness", "Rank_RWR"), names(g))], by = "Ensembl", all.x = TRUE)
  names(x)[-(1:2)] <- paste0(names(x)[-(1:2)], "_", nm)
  G <- if (is.null(G)) x else merge(G, x, by = c("Ensembl", "Gene"), all = TRUE)
}
rk <- paste0("Evidence_rank_", names(ctx)); tier <- paste0("Evidence_tier_", names(ctx)); mod <- paste0("Module_", names(ctx))
G$In_all_networks <- stats::complete.cases(G[, rk, drop = FALSE])
G$N_contexts_top <- rowSums(sapply(rk, function(k) !is.na(G[[k]]) & G[[k]] <= TOP_RANK))
G$N_contexts_High <- rowSums(sapply(tier, function(k) G[[k]] %in% "High"))
G$Rank_gap <- apply(G[, rk, drop = FALSE], 1, function(r) if (all(is.na(r))) NA else diff(range(r, na.rm = TRUE)))
## where does each gene sit in every context, including modules that were not selected, and is it changed there?
for (nm in names(ctx)) {
  am <- ctx[[nm]]$am; G[[paste0("Module_all_", nm)]] <- if (is.null(am)) NA else am$ModuleColor[match(G$Ensembl, am$Gene)]
  x <- ctx[[nm]]$deg; m <- match(G$Ensembl, x$gene)
  G[[paste0("log2FC_", nm)]] <- x$log2FoldChange[m]; G[[paste0("padj_", nm)]] <- x$padj[m]
}
changed_everywhere_same_way <- function(i) {
  lf <- unlist(G[i, paste0("log2FC_", names(ctx))]); pa <- unlist(G[i, paste0("padj_", names(ctx))])
  all(!is.na(pa) & pa < PADJ) && length(unique(sign(lf))) == 1
}
G$Changed_in_all_same_direction <- vapply(seq_len(nrow(G)), changed_everywhere_same_way, TRUE)
G$Class <- ifelse(G$N_contexts_High == length(ctx), "core: High in every context",
           ifelse(G$N_contexts_High >= 1 & !G$In_all_networks & G$Changed_in_all_same_direction,
                  "module rewiring: changed the same way everywhere, but outside the selected modules of another context",
           ifelse(G$N_contexts_High >= 1 & !G$In_all_networks, "context-specific: High in one context, not changed in another",
           ifelse(G$N_contexts_High >= 1 & G$Rank_gap >= REWIRE_GAP, "rewired: High in one context, far down in another",
           ifelse(G$N_contexts_High >= 1, "High in some contexts, close in the others", "not prioritised")))))
G <- G[order(-G$N_contexts_High, apply(G[, rk, drop = FALSE], 1, function(r) suppressWarnings(min(r, na.rm = TRUE)))), ]
write.csv(G, file.path(OUT, "context_gene_comparison_and_rewiring.csv"), row.names = FALSE)
cat("\n========== 2. prioritised genes across contexts ==========\n"); print(table(G$Class[G$Class != "not prioritised"]))
core <- G[G$Class == "core: High in every context", ]
cat("  core genes (", nrow(core), "):", paste(core$Gene, collapse = ", "), "\n")
rew <- G[G$Class == "rewired: High in one context, far down in another", ]
if (nrow(rew)) { cat("  rewired examples:\n"); print(utils::head(rew[, c("Gene", rk, mod)], 12), row.names = FALSE) }
mr <- G[grepl("^module rewiring", G$Class), ]
if (nrow(mr)) for (nm in names(ctx)) { s_ <- mr[mr[[paste0("Evidence_tier_", nm)]] %in% "High", ]
  if (!nrow(s_)) next
  other <- setdiff(names(ctx), nm)
  cat(sprintf("  module rewiring, High in %s (%d) -- module in %s: %s\n", nm, nrow(s_), paste(other, collapse = "/"),
      paste(sprintf("%s[%s]", s_$Gene, apply(s_[, paste0("Module_all_", other), drop = FALSE], 1, paste, collapse = "/")), collapse = ", "))) }
spec <- G[grepl("^context-specific", G$Class), ]
for (nm in names(ctx)) { s_ <- spec$Gene[spec[[paste0("Evidence_tier_", nm)]] %in% "High"]
  cat(sprintf("  truly specific to %s (%d): %s\n", nm, length(s_), paste(utils::head(s_, 25), collapse = ", "))) }
## module correspondence: shared genes between each pair of selected modules
if (length(ctx) >= 2) {
  pairs <- utils::combn(names(ctx), 2, simplify = FALSE); MC <- list()
  for (p in pairs) {
    a <- ctx[[p[1]]]$ev; b <- ctx[[p[2]]]$ev; univ <- intersect(a$Ensembl, b$Ensembl)
    for (ma in unique(a$Module)) for (mb in unique(b$Module)) {
      A <- a$Ensembl[a$Module == ma]; B <- b$Ensembl[b$Module == mb]; k <- length(intersect(A, B))
      exp <- length(intersect(A, univ)) * length(intersect(B, univ)) / max(1, length(univ))
      MC[[length(MC) + 1]] <- data.frame(Context_A = p[1], Module_A = ma, Name_A = unname(ctx[[p[1]]]$mods[ma]),
        Context_B = p[2], Module_B = mb, Name_B = unname(ctx[[p[2]]]$mods[mb]), Shared_genes = k, Expected = round(exp, 1))
    }
  }
  MC <- do.call(rbind, MC); MC$Enrichment <- round(MC$Shared_genes / pmax(MC$Expected, 0.1), 2)
  write.csv(MC, file.path(OUT, "context_module_correspondence.csv"), row.names = FALSE)
  cat("\n  module correspondence (shared genes / expected):\n"); print(MC[MC$Shared_genes >= 10, c("Context_A", "Module_A", "Name_A", "Context_B", "Module_B", "Name_B", "Shared_genes", "Expected")], row.names = FALSE)
}

## ---- 3. pathways ----
if (all(vapply(ctx, function(x) !is.null(x$gs), TRUE))) {
  P <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), lapply(names(ctx), function(nm) {
    x <- ctx[[nm]]$gs[, c("ID", "NES", "p.adjust")]; names(x)[-1] <- paste0(c("NES_", "padj_"), nm); x }))
  sig <- sapply(names(ctx), function(nm) !is.na(P[[paste0("padj_", nm)]]) & P[[paste0("padj_", nm)]] < PADJ)
  P$Shared_same_sign <- rowSums(sig) == length(ctx) & apply(sign(P[, paste0("NES_", names(ctx)), drop = FALSE]), 1, function(s) length(unique(s)) == 1)
  write.csv(P, file.path(OUT, "context_hallmark_pathways.csv"), row.names = FALSE)
  cat("\n========== 3. Hallmark pathways ==========\n  significant in every context with the same sign:", sum(P$Shared_same_sign), "\n")
  print(P[P$Shared_same_sign, c("ID", paste0("NES_", names(ctx)))][order(-abs(P$NES_TB[P$Shared_same_sign])), ], row.names = FALSE)
}

## ---- 4. drugs ----
if (all(vapply(ctx, function(x) !is.null(x$dr), TRUE))) {
  R <- Reduce(function(a, b) merge(a, b, by = "Drug", all = TRUE), lapply(names(ctx), function(nm) {
    x <- ctx[[nm]]$dr[, c("Drug", "Evidence_rank", "Evidence_tier", "Stable_topN")]; names(x)[-1] <- paste0(names(x)[-1], "_", nm); x }))
  R$N_contexts_stable <- rowSums(sapply(names(ctx), function(nm) R[[paste0("Stable_topN_", nm)]] %in% TRUE))
  R$N_contexts_top25 <- rowSums(sapply(names(ctx), function(nm) !is.na(R[[paste0("Evidence_rank_", nm)]]) & R[[paste0("Evidence_rank_", nm)]] <= 25))
  R <- R[order(-R$N_contexts_top25, -R$N_contexts_stable), ]
  write.csv(R, file.path(OUT, "context_drug_comparison.csv"), row.names = FALSE)
  both <- R[R$N_contexts_top25 == length(ctx), ]
  cat("\n========== 4. drugs ==========\n  in the top 25 of every context (", nrow(both), "):", paste(both$Drug, collapse = ", "), "\n")
  st <- R$Drug[R$N_contexts_stable == length(ctx)]; cat("  stable top 20 in every context (", length(st), "):", paste(st, collapse = ", "), "\n")
}
write.csv(D, file.path(OUT, "context_DEG_comparison.csv"), row.names = FALSE)
cat("\n>>> DONE. Outputs in", OUT, "\n")
