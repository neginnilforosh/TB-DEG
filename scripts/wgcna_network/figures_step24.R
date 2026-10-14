## STEP 24: manuscript figure set for one run (Ratul, 10 Oct: the A-G figure list, "generate the figures
## directly from the analysis scripts", 600 DPI raster + vector PDF, consistent labels and sizes).
##
## Reads only the tables earlier steps saved -- nothing is recomputed -- and draws the figures that did not
## exist yet. Each figure is independent: if its input is missing it is skipped with a message, so the
## script works on any run (pilot, TB-DM, later TB-HIV). Run it once per run:
##     Rscript figures_step24.R                                   # active run
##     TB_RUN=GSE114192_TBDM_vs_DM Rscript figures_step24.R
## Outputs in <run>/figures_manuscript/ (PNG 600 dpi + PDF of every figure).

## ---- CONFIG ----
FIG_DPI <- 600
TOP_N   <- 20            # bars in "top" plots
TOP_TERMS <- 8           # terms per module in enrichment dot plots
BASE_SIZE <- 9

for (p in c("ggplot2", "igraph")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages({ library(ggplot2) })
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR <- file.path(ROOT_DIR, ACC); OUT_DIR <- run_out_dir(BASE_DIR)
FIG <- file.path(OUT_DIR, "figures_manuscript"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
P <- function(...) file.path(OUT_DIR, ...)
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
theme_set(theme_bw(base_size = BASE_SIZE) + theme(panel.grid.minor = element_blank()))
mod_label <- function(m) ifelse(!is.na(MODULE_NAMES[m]) & nzchar(MODULE_NAMES[m]), paste0(m, " (", MODULE_NAMES[m], ")"), m)
CMP <- paste0(CASE, " vs ", CONTROL)

made <- character(0); skipped <- character(0)
save_gg <- function(name, p, w, h) {
  f <- file.path(FIG, paste0(ACC, "_", name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = FIG_DPI, bg = "white"); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h, bg = "white")
  made <<- c(made, name)
}
save_base <- function(name, w, h, draw) {
  f <- file.path(FIG, paste0(ACC, "_", name, ".png"))
  png(f, width = w, height = h, units = "in", res = FIG_DPI); draw(); dev.off()
  pdf(sub("\\.png$", ".pdf", f), width = w, height = h); draw(); dev.off()
  made <<- c(made, name)
}
fig <- function(name, expr) tryCatch(expr, error = function(e) { skipped <<- c(skipped, paste0(name, ": ", conditionMessage(e))) })
need <- function(x, what) if (is.null(x) || !NROW(x)) stop("missing input (", what, ")")
cat(">>> figures for", run_id, "->", FIG, "\n")

## ======================= B. DEG =======================
fig("B_DEG_counts_bar", {
  fs <- list.files(file.path(in_dir(BASE_DIR, "06_deg_results")), pattern = "_DEG\\.csv$", full.names = TRUE)
  fs <- fs[!grepl("interaction", fs)]; need(fs, "06_deg_results")
  d <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f); s <- !is.na(x$padj) & x$padj < 0.05 & abs(x$log2FoldChange) >= 1
    data.frame(Comparison = sub("_DEG\\.csv$", "", sub(paste0("^", ACC, "_"), "", basename(f))),
               Direction = c("Up", "Down"), N = c(sum(s & x$log2FoldChange > 0), sum(s & x$log2FoldChange < 0))) }))
  d$Comparison <- gsub("_vs_", " vs ", d$Comparison)
  save_gg("B_DEG_counts_bar", ggplot(d, aes(Comparison, N, fill = Direction)) + geom_col(position = "dodge") +
    geom_text(aes(label = N), position = position_dodge(0.9), vjust = -0.3, size = 2.6) +
    scale_fill_manual(values = c(Up = "#b2182b", Down = "#2166ac")) + coord_cartesian(clip = "off") +
    labs(title = "Significant DEGs (padj < 0.05, |log2FC| >= 1)", x = NULL, y = "Genes") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)), 5.5, 4)
})

## ======================= A. QC: sample clustering / outliers =======================
fig("A_sample_clustering_outliers", {
  base <- dirname(in_dir(BASE_DIR, "02_metadata"))
  vf <- list.files(file.path(base, "04_normalized"), pattern = "_vst_normalized_matrix\\.csv$", full.names = TRUE)[1]
  mf <- list.files(file.path(base, "02_metadata"), pattern = "_sample_metadata\\.csv$", full.names = TRUE)[1]
  if (is.na(vf) || is.na(mf)) stop("missing input (vst matrix / metadata)")
  x <- as.matrix(read.csv(vf, row.names = 1, check.names = FALSE)); md <- read.csv(mf, stringsAsFactors = FALSE)
  grp <- setNames(md$group, md$sample_id)[colnames(x)]
  x <- x[order(-apply(x, 1, stats::var))[seq_len(min(2000, nrow(x)))], ]          # top-variance genes
  cc <- stats::cor(x); hc <- stats::hclust(stats::as.dist(1 - cc), "average")
  mc <- (colSums(cc) - 1) / (ncol(cc) - 1); thr <- stats::quantile(mc, .25) - 1.5 * stats::IQR(mc)
  out <- names(mc)[mc < thr]
  write.csv(data.frame(Sample = names(mc), Group = grp, Mean_correlation = round(mc, 4), Outlier = mc < thr),
            file.path(FIG, paste0(ACC, "_A_sample_outliers.csv")), row.names = FALSE)
  cols <- setNames(grDevices::hcl.colors(length(unique(grp)), "Dark 3"), unique(grp))
  save_base("A_sample_clustering_outliers", max(7, ncol(x) * 0.07), 5, function() {
    op <- graphics::par(mar = c(6, 4, 3, 1)); on.exit(graphics::par(op))
    dd <- stats::as.dendrogram(hc)
    dd <- stats::dendrapply(dd, function(n) { if (is.leaf(n)) { l <- attr(n, "label")
      attr(n, "nodePar") <- list(pch = NA, lab.col = if (l %in% out) "red" else cols[grp[l]], lab.cex = 0.45) }; n })
    plot(dd, main = paste0("Sample clustering (1 - Pearson, top-variance genes); red = outlier (", length(out), ")"), ylab = "1 - r")
    graphics::legend("topright", legend = names(cols), text.col = cols, bty = "n", cex = 0.7, inset = c(0, 0.08))
  })
})

## ======================= C. WGCNA =======================
ws <- if (file.exists(P("10_wgcna_results", paste0(ACC, "_WGCNA_workspace.rds")))) readRDS(P("10_wgcna_results", paste0(ACC, "_WGCNA_workspace.rds"))) else NULL
fig("C_eigengene_clustering", {
  need(ws, "WGCNA workspace"); me <- ws$MEs; me <- me[, colnames(me) != "MEgrey", drop = FALSE]
  hc <- stats::hclust(stats::as.dist(1 - stats::cor(me)), "average")
  save_base("C_eigengene_clustering", 7, 4.5, function() {
    plot(hc, main = "Module eigengene clustering (merge height 0.25 marked)", xlab = "", sub = "", ylab = "1 - r", cex = 0.7)
    graphics::abline(h = 0.25, col = "red", lty = 2) })
})
mmgs <- rd(P("10_wgcna_results", paste0(ACC, "_TBmodules_MM_GS.csv")))
fig("C_MM_vs_GS", {
  need(mmgs, "TBmodules_MM_GS"); mmgs$ModuleLab <- mod_label(mmgs$Module)
  r <- do.call(rbind, lapply(split(mmgs, mmgs$ModuleLab), function(d) { ct <- stats::cor.test(abs(d$MM), abs(d$GS_TB))
    data.frame(ModuleLab = d$ModuleLab[1], lab = sprintf("r = %.2f, p = %.1e", ct$estimate, ct$p.value)) }))
  save_gg("C_MM_vs_GS", ggplot(mmgs, aes(abs(MM), abs(GS_TB))) + geom_point(aes(colour = Module), size = 0.5, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.3, size = 2.6) +
    scale_colour_identity() + facet_wrap(~ ModuleLab) +
    labs(title = "Module membership vs gene significance", x = "|MM|", y = paste0("|GS| (", CMP, ")")), 7.5, 3.2)
})

## ======================= D. enrichment =======================
enr <- do.call(rbind, lapply(c("GO_BP", "KEGG", "Reactome"), function(db) do.call(rbind, lapply(names(MODULE_NAMES), function(m) {
  d <- rd(P("11_module_enrichment", paste0(ACC, "_", m, "_", db, ".csv")))
  if (is.null(d) || !nrow(d) || !"p.adjust" %in% names(d)) return(NULL)
  d <- d[d$p.adjust < 0.05, ]; if (!nrow(d)) return(NULL)
  gr <- vapply(strsplit(d$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 0)
  data.frame(DB = db, Module = m, Description = d$Description, GeneRatio = gr, Count = d$Count, padj = d$p.adjust)
}))))
for (db in c("GO_BP", "KEGG", "Reactome")) fig(paste0("D_dotplot_", db), {
  d <- enr[!is.null(enr) & enr$DB == db, , drop = FALSE]; need(d, paste(db, "terms with padj < 0.05"))
  d <- do.call(rbind, lapply(split(d, d$Module), function(z) utils::head(z[order(z$padj), ], TOP_TERMS)))
  d$Description <- substr(d$Description, 1, 60); d$ModuleLab <- mod_label(d$Module)
  d$Term <- factor(paste(d$Description), levels = unique(d$Description[order(d$Module, -d$padj)]))
  save_gg(paste0("D_dotplot_", db), ggplot(d, aes(GeneRatio, Term, size = Count, colour = -log10(padj))) + geom_point() +
    facet_grid(ModuleLab ~ ., scales = "free_y", space = "free_y") + scale_colour_viridis_c() +
    labs(title = paste0(db, " enrichment of the selected modules"), y = NULL) +
    theme(strip.text.y = element_text(angle = 0)), 7.5, max(3, 0.18 * nrow(d) + 1.2))
})
fig("D_module_pathway_heatmap", {
  need(enr, "enrichment terms"); top <- unique(unlist(lapply(split(enr, enr$Module), function(z) utils::head(z$Description[order(z$padj)], 6))))
  d <- expand.grid(Description = top, Module = names(MODULE_NAMES), stringsAsFactors = FALSE)
  best <- stats::aggregate(padj ~ Description + Module, enr[enr$Description %in% top, ], min)
  d <- merge(d, best, all.x = TRUE); d$score <- ifelse(is.na(d$padj), 0, -log10(d$padj)); d$ModuleLab <- mod_label(d$Module)
  d$Description <- substr(d$Description, 1, 60)
  save_gg("D_module_pathway_heatmap", ggplot(d, aes(ModuleLab, Description, fill = score)) + geom_tile(colour = "white") +
    scale_fill_gradient(low = "white", high = "#b2182b", name = "-log10 padj") +
    labs(title = "Top enriched terms across modules (GO BP / KEGG / Reactome)", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)), 7.5, max(3, 0.2 * length(top) + 1.5))
})

## ======================= E. network =======================
cen <- rd(P("12_string_ppi", paste0(ACC, "_COMBINED_centrality.csv")))
gs  <- rd(P("15_network_algorithms", paste0(ACC, "_GeneScores_RWR_Diffusion.csv")))
gp  <- rd(P("19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
topbar <- function(name, d, val, title, xlab) fig(name, {
  need(d, name); d <- utils::head(d[order(-d[[val]]), ], TOP_N); d$Gene <- factor(d$Gene, levels = rev(d$Gene))
  save_gg(name, ggplot(d, aes(.data[[val]], Gene, fill = Module)) + geom_col() + scale_fill_identity(guide = "legend",
    labels = function(x) mod_label(x)) + labs(title = title, x = xlab, y = NULL, fill = "Module"), 5.5, 4.8)
})
topbar("E_top_PPI_hubs", cen, "Degree", "Top PPI hubs (combined network)", "Degree")
topbar("E_top_RWR_genes", gs, "RWR_score", "Top genes by random walk with restart", "RWR score")
topbar("E_top_diffusion_genes", gs, "Diffusion_score", "Top genes by network diffusion", "Diffusion score")

net_edges <- rd(P("12_string_ppi", paste0(ACC, "_COMBINED_STRING_network.csv")))
comm <- rd(P("15_network_algorithms", paste0(ACC, "_Communities.csv")))
fig("E_network_labelled", {
  need(net_edges, "combined network"); need(cen, "centrality")
  g <- igraph::graph_from_data_frame(unique(net_edges[, c("preferredName_A", "preferredName_B")]), directed = FALSE)
  g <- igraph::simplify(g)
  cmp <- igraph::components(g); g <- igraph::induced_subgraph(g, which(cmp$membership == which.max(cmp$csize)))  # small detached pieces squash the layout
  set.seed(42); lay <- igraph::layout_with_fr(g)
  v <- data.frame(Gene = igraph::V(g)$name, x = lay[, 1], y = lay[, 2])
  v <- merge(v, cen[, c("Gene", "Module")], all.x = TRUE)
  if (!is.null(gs)) v <- merge(v, gs[, c("Gene", "RWR_score")], all.x = TRUE) else v$RWR_score <- 1
  hi <- if (!is.null(gp)) gp$Gene[gp$Final_Priority == "High"] else character(0)
  el <- igraph::as_edgelist(g); e <- data.frame(x = lay[match(el[, 1], igraph::V(g)$name), 1], y = lay[match(el[, 1], igraph::V(g)$name), 2],
    xend = lay[match(el[, 2], igraph::V(g)$name), 1], yend = lay[match(el[, 2], igraph::V(g)$name), 2])
  lab <- v[v$Gene %in% hi, ]
  p <- ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey85", linewidth = 0.1) +
    geom_point(data = v, aes(x, y, size = RWR_score, colour = Module), alpha = 0.8) + scale_colour_identity(guide = "legend", labels = function(x) mod_label(x)) +
    scale_size_continuous(range = c(0.3, 3.5), name = "RWR score") +
    labs(title = "Combined disease PPI network", subtitle = "node size = RWR score; labels = High-priority genes", colour = "Module") +
    theme_void(base_size = BASE_SIZE)
  if (nrow(lab)) p <- p + (if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = lab, aes(x, y, label = Gene), size = 2, max.overlaps = 40, segment.size = 0.2) else
    geom_text(data = lab, aes(x, y, label = Gene), size = 2, vjust = -0.8))
  save_gg("E_network_labelled", p, 8, 7)
  ## communities: same layout, coloured by community (largest 8)
  if (!is.null(comm)) {
    v2 <- merge(v, comm[, c("Gene", "Community")], all.x = TRUE)
    big <- names(sort(table(v2$Community), decreasing = TRUE))[1:min(8, length(unique(v2$Community)))]
    v2$Comm <- ifelse(v2$Community %in% big, paste("C", v2$Community), "other")
    save_gg("E_network_communities", ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey88", linewidth = 0.1) +
      geom_point(data = v2, aes(x, y, colour = Comm), size = 0.7) + labs(title = "Network communities (largest 8)", colour = "Community") +
      theme_void(base_size = BASE_SIZE), 7.5, 6.5)
  }
})

fig("E_consensus_genes", {
  need(gp, "consensus genes"); d <- gp[gp$Final_Priority == "High", ]; need(d, "High genes")
  crit <- grep("^Top_", names(d), value = TRUE)
  L <- do.call(rbind, lapply(crit, function(cname) data.frame(Gene = d$Gene, Module = d$Module, Criterion = sub("^Top_", "", cname), Met = d[[cname]])))
  L$Gene <- factor(L$Gene, levels = d$Gene[order(d$Module, -d$N_criteria)])
  save_gg("E_consensus_genes", ggplot(L, aes(Criterion, Gene, fill = ifelse(Met %in% TRUE, Module, "white"))) + geom_tile(colour = "grey80") +
    scale_fill_identity() + labs(title = "High-priority genes: criteria met (filled = top 20%)", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 30, hjust = 1)), 4.5, max(4, 0.13 * nrow(d) + 1.2))
})

## ======================= F. GSEA =======================
fig("F_GSEA_dotplot", {
  d <- do.call(rbind, lapply(c("HALLMARK", "REACTOME", "CELLTYPE"), function(cl) {
    x <- rd(file.path(BASE_DIR, "17_gsea", paste0(ACC, "_", DEG_NAME, "_GSEA_", cl, ".csv")))
    if (is.null(x)) return(NULL); x <- x[x$p.adjust < 0.05, ]; if (!nrow(x)) return(NULL)
    x <- rbind(utils::head(x[order(-x$NES), ], 8), utils::head(x[order(x$NES), ], 8)); x <- x[!duplicated(x$ID), ]
    data.frame(Collection = cl, Term = substr(sub("^(HALLMARK|REACTOME)_", "", x$Description), 1, 55), NES = x$NES, padj = x$p.adjust, Size = x$setSize) }))
  need(d, "17_gsea tables"); d$Term <- factor(d$Term, levels = unique(d$Term[order(d$NES)]))
  save_gg("F_GSEA_dotplot", ggplot(d, aes(NES, Term, size = Size, colour = -log10(padj))) + geom_point() + geom_vline(xintercept = 0, linetype = 2) +
    facet_grid(Collection ~ ., scales = "free_y", space = "free_y") + scale_colour_viridis_c() +
    labs(title = paste0("GSEA, ", CMP, " (padj < 0.05)"), y = NULL) + theme(strip.text.y = element_text(angle = 0)), 7.5, max(4, 0.15 * nrow(d) + 1.5))
})

## ======================= G. drugs =======================
## use the same cleaned target lists as the consensus (step 15 drops screening-only and ADME targets)
clean_dt <- function(dt, dp) {
  if (is.null(dt) || is.null(dp) || !"Known_Targets_used" %in% names(dp)) return(dt)
  used <- strsplit(ifelse(is.na(dp$Known_Targets_used), "", dp$Known_Targets_used), ";"); names(used) <- dp$Drug
  dt[mapply(function(dr, tg) tg %in% used[[dr]], dt$Drug, dt$Drug_Target), , drop = FALSE]
}
dp <- rd(P("19_consensus", paste0(ACC, "_DrugPriority.csv")))
fig("G_iLINCS_top_reversal", {
  x <- rd(P("13_ilincs_signature", paste0(ACC, "_iLINCS_candidate_compounds.csv"))); need(x, "iLINCS candidates")
  x <- x[!grepl("^BRD-", x$Target), ]; x <- stats::aggregate(Similarity ~ Target, x, min)
  x <- utils::head(x[order(x$Similarity), ], TOP_N); x$Target <- factor(x$Target, levels = rev(x$Target))
  save_gg("G_iLINCS_top_reversal", ggplot(x, aes(Similarity, Target)) + geom_col(fill = "#2166ac") +
    labs(title = "Strongest reversal of the disease signature (iLINCS, named compounds)", x = "Concordance (more negative = stronger reversal)", y = NULL), 6, 4.8)
})
fig("G_network_proximity", {
  x <- rd(P("15_network_algorithms", paste0(ACC, "_NetworkProximity.csv"))); need(x, "network proximity")
  x$lab <- ifelse(rank(x$P_empirical, ties.method = "first") <= 10, x$Drug, "")
  p <- ggplot(x, aes(Z_score, -log10(P_empirical), size = N_targets_in_network)) + geom_point(alpha = 0.5, colour = "grey30") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") + geom_vline(xintercept = 0, linetype = 3) +
    labs(title = "Network proximity of drug targets to the disease network", subtitle = "red line = empirical p 0.05 (no FDR correction)",
         x = "Proximity Z (negative = closer than random)", y = "-log10 empirical p", size = "Targets in network")
  if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(aes(label = lab), size = 2.2, max.overlaps = 30)
  save_gg("G_network_proximity", p, 6.5, 5)
})
dt <- clean_dt(rd(P("14_drug_targets", paste0(ACC, "_Drug_Target_TBnetwork.csv"))), dp)
fig("G_drug_module_heatmap", {
  need(dp, "drug priority"); need(dt, "drug targets")
  keep <- dp$Drug[dp$Final_Priority %in% c("High", "Medium")]; d <- dt[dt$Drug %in% keep, ]; need(d, "High/Medium drug targets")
  d <- stats::aggregate(Drug_Target ~ Drug + TB_Module, d, function(z) length(unique(z)))
  d$TB_Module <- mod_label(d$TB_Module); d$Drug <- factor(d$Drug, levels = rev(intersect(dp$Drug[order(-dp$N_criteria)], d$Drug)))
  save_gg("G_drug_module_heatmap", ggplot(d, aes(TB_Module, Drug, fill = Drug_Target)) + geom_tile(colour = "white") + geom_text(aes(label = Drug_Target), size = 2) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15", name = "Targets") + labs(title = "High/Medium drugs: targets per module", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = 5.5), axis.text.x = element_text(angle = 20, hjust = 1)), 5.5, max(4, 0.1 * length(unique(d$Drug)) + 1.5))
})
fig("G_drug_target_network", {
  need(dp, "drug priority"); need(dt, "drug targets")
  hd <- dp$Drug[dp$Final_Priority == "High"]; d <- unique(dt[dt$Drug %in% hd, c("Drug", "Drug_Target", "TB_Module")]); need(d, "High drug targets")
  g <- igraph::graph_from_data_frame(d[, 1:2], directed = FALSE); set.seed(42); lay <- igraph::layout_with_fr(g)
  v <- data.frame(name = igraph::V(g)$name, x = lay[, 1], y = lay[, 2]); v$type <- ifelse(v$name %in% d$Drug, "Drug", "Target")
  v$col <- ifelse(v$type == "Drug", "grey20", d$TB_Module[match(v$name, d$Drug_Target)])
  el <- igraph::as_edgelist(g); e <- data.frame(x = v$x[match(el[, 1], v$name)], y = v$y[match(el[, 1], v$name)], xend = v$x[match(el[, 2], v$name)], yend = v$y[match(el[, 2], v$name)])
  save_gg("G_drug_target_network", ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey70", linewidth = 0.3) +
    geom_point(data = v, aes(x, y, shape = type, colour = col), size = 3) + scale_colour_identity() + scale_shape_manual(values = c(Drug = 15, Target = 16)) +
    geom_text(data = v, aes(x, y, label = name), size = 2, vjust = -1.1) + labs(title = "High-priority drugs and their disease-network targets", subtitle = "target colour = module", shape = NULL) +
    theme_void(base_size = BASE_SIZE), 7.5, 6.5)
})
fig("G_consensus_drugs", {
  need(dp, "drug priority"); d <- dp[dp$Final_Priority %in% c("High", "Medium"), ]; need(d, "High/Medium drugs")
  crit <- c("Strong_reversal", "Negative_prox_Z", "Significant_prox_P", "Targets_in_TBnet", "Module_enriched", "Hits_high_gene", "Top_target_RWR")
  crit <- intersect(crit, names(d)); d <- utils::head(d[order(-d$N_criteria, d$iLINCS_Reversal_Score), ], 40)
  L <- do.call(rbind, lapply(crit, function(cn) data.frame(Drug = d$Drug, Criterion = cn, Met = d[[cn]] %in% TRUE, Tier = d$Final_Priority)))
  L$Drug <- factor(L$Drug, levels = rev(d$Drug))
  save_gg("G_consensus_drugs", ggplot(L, aes(Criterion, Drug, fill = ifelse(Met, Tier, "none"))) + geom_tile(colour = "grey80") +
    scale_fill_manual(values = c(High = "#a50f15", Medium = "#fc9272", none = "white"), name = "Tier (met)") +
    labs(title = "Drug consensus: criteria met (top 40 High/Medium)", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6)), 5.5, max(4, 0.13 * nrow(d) + 1.5))
})

cat("\n>>> made", length(made), "figures (PNG 600 dpi + PDF):\n   ", paste(made, collapse = "\n    "), "\n")
if (length(skipped)) cat("\n>>> skipped", length(skipped), ":\n   ", paste(skipped, collapse = "\n    "), "\n")
cat(">>> DONE. Folder:", FIG, "\n")
