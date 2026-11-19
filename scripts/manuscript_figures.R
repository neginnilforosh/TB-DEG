

for (p in c("ggplot2", "pheatmap", "UpSetR", "ggrepel", "scales"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT, "manuscript_figures"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(ROOT, ...); if (!file.exists(f)) { message("  (missing: ", f, ")"); return(NULL) }
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }

CTX <- list(
  TB     = list(run = "GSE114192",          acc = "GSE114192", label = "TB vs healthy"),
  TB_DM  = list(run = "GSE114192/TB_DM",    acc = "GSE114192", label = "TB with diabetes vs diabetes"),
  TB_HIV = list(run = "GSE37250/TB_HIV",    acc = "GSE37250",  label = "Active vs latent TB, HIV-positive"))
COL <- c(up = "#B2182B", down = "#2166AC", ns = "grey75")
theme_ms <- theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10), legend.key.size = unit(3.5, "mm"),
                                                  strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_gg <- function(p, name, w = 85, h = 70) {
  wrap <- function(x, n) if (is.null(x)) x else paste(strwrap(x, n), collapse = "\n")
  p$labels$title <- wrap(p$labels$title, floor(w / 1.9)); p$labels$subtitle <- wrap(p$labels$subtitle, floor(w / 1.6))
  for (ext in c("png", "pdf")) ggsave(file.path(OUT, paste0(name, ".", ext)), p, width = w, height = h, units = "mm", dpi = 600)
  message("  saved ", name)
}
save_base <- function(fun, name, w = 85, h = 70) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "mm", res = 600); fun(); dev.off()
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w / 25.4, height = h / 25.4, onefile = FALSE); fun(); dev.off()
  message("  saved ", name)
}
fam_cols <- c(Fam_Expression = "Expression", Fam_Coexpression = "Co-expression", Fam_Network = "Network",
              Fam_Replication = "Replication", Fam_Infection = "Infection response")

## ---- 2A volcano of the reference comparison, labelled with gene symbols ----
message("== TB volcano")
dv <- rd("GSE114192", "06_deg_results", "GSE114192_HealthyControl_vs_TBOnly_DEG.csv")
if (!is.null(dv)) {
  ens <- sub("\\..*$", "", dv$gene); sym <- setNames(rep(NA_character_, nrow(dv)), ens)
  if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE))
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  evv <- rd("GSE114192", "20_evidence_families", "GSE114192_Gene_evidence_families.csv")
  if (!is.null(evv)) { m <- setNames(evv$Gene, evv$Ensembl); fill <- is.na(sym) & ens %in% names(m); sym[fill] <- m[ens[fill]] }
  dv$symbol <- unname(sym); dv <- dv[!is.na(dv$padj), ]
  dv$dir <- ifelse(dv$padj < 0.05 & dv$log2FoldChange >= 1, "up", ifelse(dv$padj < 0.05 & dv$log2FoldChange <= -1, "down", "ns"))
  key <- c("STAT1", "GBP1", "GBP5", "BATF2", "FCGR1A", "ANKRD22", "SERPING1", "CASP1", "AIM2", "IFI27")
  top <- order(dv$padj)[seq_len(min(12, nrow(dv)))]
  dv$Label <- ifelse((seq_len(nrow(dv)) %in% top | dv$symbol %in% key) & dv$dir != "ns" & !is.na(dv$symbol), dv$symbol, "")
  p <- ggplot(dv, aes(log2FoldChange, -log10(padj), colour = dir)) + geom_point(size = 0.5, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = "grey50") + geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_text_repel(aes(label = Label), size = 2.1, colour = "black", max.overlaps = 40, segment.size = 0.2) + scale_colour_manual(values = COL, guide = "none") +
    labs(title = "TB vs healthy controls (GSE114192)", subtitle = sprintf("%d up, %d down (FDR < 0.05, |log2FC| >= 1)", sum(dv$dir == "up"), sum(dv$dir == "down")),
         x = "log2 fold change", y = "-log10 adjusted p") + theme_ms
  save_gg(p, "Fig2A_volcano_TB")
}

for (k in names(CTX)) {
  c0 <- CTX[[k]]; message("== ", k)
  ## ---- 2C module-trait ----
  mt <- rd(c0$run, "10_wgcna_results", paste0(c0$acc, "_ModuleTrait_Correlation.csv"))
  if (!is.null(mt)) {
    mt <- mt[mt$Module != "grey", ]; mt$Module <- factor(mt$Module, levels = mt$Module[order(mt$Correlation)])
    mt$sig <- ifelse(mt$FDR < 0.001, "***", ifelse(mt$FDR < 0.01, "**", ifelse(mt$FDR < 0.05, "*", "")))
    p <- ggplot(mt, aes(Correlation, Module, fill = Correlation)) + geom_col(width = 0.75, colour = ifelse(mt$Selected, "black", NA)) +
      geom_text(aes(label = sig, hjust = ifelse(Correlation > 0, -0.2, 1.2)), size = 2.6) +
      scale_fill_gradient2(low = COL["down"], high = COL["up"], mid = "white", limits = c(-1, 1)) +
      scale_x_continuous(limits = c(-1, 1)) + labs(title = paste0("Module-trait correlation: ", c0$label), x = "Correlation with disease (eigengene)", y = NULL) +
      theme_ms + theme(legend.position = "none")
    save_gg(p, paste0("Fig2C_module_trait_", k), 85, 80)
  }
  ## ---- 2D dendrogram ----
  ws <- file.path(ROOT, c0$run, "10_wgcna_results", paste0(c0$acc, "_WGCNA_workspace.rds"))
  if (file.exists(ws)) {
    w <- readRDS(ws); net <- if (!is.null(w$net)) w$net else w
    if (!is.null(net$dendrograms) && !is.null(net$colors)) {
      tree <- net$dendrograms[[1]]; cols <- net$colors[net$blockGenes[[1]]]
      if (is.numeric(cols) && requireNamespace("WGCNA", quietly = TRUE)) cols <- WGCNA::labels2colors(cols)
      save_base(function() {
        layout(matrix(1:2), heights = c(4, 1)); par(mar = c(0, 3, 2, 1), cex = 0.6)
        plot(tree, labels = FALSE, hang = 0.03, main = paste0("Gene dendrogram: ", c0$label), xlab = "", sub = "", ylab = "Height")
        par(mar = c(1, 3, 0, 1)); image(matrix(seq_along(tree$order)), col = cols[tree$order], axes = FALSE)
        mtext("Module", side = 2, las = 1, line = 0.5, cex = 0.6)
      }, paste0("Fig2D_dendrogram_", k), 170, 80)
    } else message("  (workspace has no dendrogram)")
  }
  ## ---- 3B degree distribution, 3D score agreement ----
  ce <- rd(c0$run, "12_string_ppi", paste0(c0$acc, "_COMBINED_centrality.csv"))
  sc <- rd(c0$run, "15_network_algorithms", paste0(c0$acc, "_GeneScores_RWR_Diffusion.csv"))
  if (!is.null(ce)) {
    dd <- as.data.frame(table(ce$Degree)); names(dd) <- c("Degree", "N"); dd$Degree <- as.numeric(as.character(dd$Degree)); dd <- dd[dd$Degree > 0, ]
    p <- ggplot(dd, aes(Degree, N)) + geom_point(size = 0.9, alpha = 0.7) + scale_x_log10() + scale_y_log10() +
      geom_smooth(method = "lm", se = FALSE, colour = COL["up"], linewidth = 0.5, formula = y ~ x) +
      labs(title = paste0("Degree distribution: ", c0$label), subtitle = sprintf("%d proteins", nrow(ce)), x = "Degree (log)", y = "Number of proteins (log)") + theme_ms
    save_gg(p, paste0("Fig3B_degree_distribution_", k))
  }
  if (!is.null(ce) && !is.null(sc)) {
    m <- merge(ce[, c("Gene", "Degree", "Betweenness", "Closeness", "Eigenvector")], sc[, c("Gene", "RWR_score", "Diffusion_score")], by = "Gene")
    names(m) <- c("Gene", "Degree", "Betweenness", "Closeness", "Eigenvector", "RWR", "Diffusion")
    cm <- stats::cor(m[, -1], method = "spearman", use = "pairwise.complete.obs")
    save_base(function() pheatmap::pheatmap(cm, display_numbers = TRUE, number_format = "%.2f", fontsize = 7, fontsize_number = 6,
                                            color = grDevices::colorRampPalette(c("white", "#FDDBC7", COL["up"]))(50), breaks = seq(0, 1, length.out = 51),
                                            main = paste0("Network score agreement (Spearman)\n", c0$label), border_color = NA),
              paste0("Fig3D_score_agreement_", k), 120, 95)
    rk <- m; for (v in names(m)[-1]) rk[[v]] <- rank(-m[[v]], ties.method = "min")
    top <- rk[order(rowMeans(rk[, -1])), ][1:20, ]; rownames(top) <- top$Gene
    save_base(function() pheatmap::pheatmap(as.matrix(top[, -1]), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
                                            number_format = "%.0f", fontsize = 7, fontsize_number = 5.5, border_color = NA,
                                            color = grDevices::colorRampPalette(c(COL["up"], "#FDDBC7", "white"))(50),
                                            main = paste0("Ranks of the top 20 genes by each score\n", c0$label)),
              paste0("Fig3D_top_gene_ranks_", k), 110, 110)
  }
  ## ---- 4A UpSet of strong evidence families, 4B top genes ----
  ev <- rd(c0$run, "20_evidence_families", paste0(c0$acc, "_Gene_evidence_families.csv"))
  if (!is.null(ev)) {
    ev$Gene[is.na(ev$Gene) | ev$Gene == ""] <- ev$Ensembl[is.na(ev$Gene) | ev$Gene == ""]
    fc <- intersect(names(fam_cols), names(ev))
    sets <- lapply(setNames(fc, fam_cols[fc]), function(f) ev$Gene[!is.na(ev[[f]]) & ev[[f]] >= 0.8])
    sets <- sets[vapply(sets, length, 1) > 0]
    save_base(function() print(UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 25, order.by = "freq",
                                             mainbar.y.label = "Genes supported by this combination", sets.x.label = "Genes per family",
                                             text.scale = 1.1, main.bar.color = "grey30", sets.bar.color = COL["up"])),
              paste0("Fig4A_upset_evidence_families_", k), 150, 95)
    top <- ev[order(ev$Evidence_rank), ][1:30, ]
    mat <- as.matrix(top[, fc]); colnames(mat) <- fam_cols[fc]; rownames(mat) <- paste0(top$Gene, ifelse(top$Stable_topN %in% TRUE, " *", ""))
    ann <- data.frame(Tier = top$Evidence_tier, Module = top$Module, row.names = rownames(mat))
    modcol <- setNames(unique(top$Module), unique(top$Module)); modcol[] <- ifelse(modcol %in% grDevices::colors(), modcol, "grey50")
    save_base(function() pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f",
                                            fontsize = 7, fontsize_number = 5.5, border_color = NA, na_col = "grey92", annotation_row = ann,
                                            annotation_colors = list(Tier = c(High = COL[["up"]], Medium = "#F4A582", Low = "grey80"), Module = modcol),
                                            color = grDevices::colorRampPalette(c("white", "#FDDBC7", COL["up"]))(50), breaks = seq(0, 1, length.out = 51),
                                            main = paste0("Top 30 genes: evidence-family scores\n", c0$label, " (* stable in the top 50)")),
              paste0("Fig4B_top_genes_evidence_", k), 130, 150)
  }
  ## ---- 6D reversal vs proximity ----
  dr <- rd(c0$run, "20_evidence_families", paste0(c0$acc, "_Drug_evidence_families.csv"))
  if (!is.null(dr) && "Network_Proximity_Z" %in% names(dr)) {
    d <- dr[!is.na(dr$Network_Proximity_Z) & !is.na(dr$iLINCS_Reversal_Score), ]
    nice <- function(x) ifelse(grepl("[0-9]", x), x, tools::toTitleCase(tolower(x)))
    d$Label <- ifelse(d$Stable_topN %in% TRUE | d$Evidence_rank <= 8, nice(d$Drug), "")
    p <- ggplot(d, aes(iLINCS_Reversal_Score, Network_Proximity_Z)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
      geom_point(aes(colour = Evidence_tier), size = 1.1, alpha = 0.8) +
      geom_text_repel(aes(label = Label), size = 2.2, max.overlaps = 30, segment.size = 0.2) +
      scale_colour_manual(values = c(High = COL[["up"]], Medium = "#F4A582", Low = "grey70")) +
      labs(title = paste0("Signature reversal vs network proximity: ", c0$label),
           subtitle = "No compound is significantly closer than random after FDR correction",
           x = "iLINCS similarity (more negative = stronger reversal)", y = "Network proximity Z (more negative = closer)", colour = "Evidence tier") + theme_ms
    save_gg(p, paste0("Fig6D_reversal_vs_proximity_", k), 110, 85)
  }
}

## ---- 4D GSEA curves for the reference (needs clusterProfiler, msigdbr, enrichplot) ----
if (all(vapply(c("clusterProfiler", "msigdbr", "enrichplot", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))) {
  de <- rd("GSE114192", "06_deg_results", "GSE114192_HealthyControl_vs_TBOnly_DEG.csv")
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H"); t2g <- unique(hs[, c("gs_name", "gene_symbol")])
  map <- suppressMessages(clusterProfiler::bitr(de$gene, "ENSEMBL", "SYMBOL", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  de <- de[!is.na(de$pvalue), ]; de$stat <- sign(de$log2FoldChange) * -log10(pmax(de$pvalue, 1e-300))
  de <- merge(de[, c("gene", "stat")], map, by.x = "gene", by.y = "ENSEMBL"); de <- de[order(-abs(de$stat)), ]; de <- de[!duplicated(de$SYMBOL), ]
  ## same gene-set size limits as WGCNA_GSEA_step11.R (15-500), so the curves are drawn from the same test; the NES
  ## and FDR in the titles are read from the step-11 Hallmark table, which is the result reported in the text
  set.seed(42); g <- clusterProfiler::GSEA(sort(setNames(de$stat, de$SYMBOL), decreasing = TRUE), TERM2GENE = t2g, minGSSize = 15, maxGSSize = 500,
                                          pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0)
  hm <- rd("GSE114192", "17_gsea", "GSE114192_HealthyControl_vs_TBOnly_GSEA_HALLMARK.csv")
  for (s in c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")) {
    r <- hm[hm$ID == s, ]
    p <- enrichplot::gseaplot2(g, geneSetID = s, title = sprintf("%s (NES %+.2f, FDR %.1e)", gsub("_", " ", sub("HALLMARK_", "", s)), r$NES, r$p.adjust), base_size = 8)
    save_gg(p, paste0("Fig4D_GSEA_", sub("HALLMARK_", "", s), "_TB"), 110, 85)
  }
} else message("(4D GSEA curves skipped: clusterProfiler/msigdbr/enrichplot not installed)")

## ---- 7A/7B TB-HIV volcano and PCA ----
message("== TB-HIV expression")
dh <- rd("GSE37250", "06_deg_results", "GSE37250_LTBI_HIVpos_vs_ATB_HIVpos_DEG.csv")
ev_h <- rd("GSE37250/TB_HIV", "20_evidence_families", "GSE37250_Gene_evidence_families.csv")
if (!is.null(dh)) {
  dh$dir <- ifelse(dh$padj < 0.05 & dh$log2FoldChange >= 1, "up", ifelse(dh$padj < 0.05 & dh$log2FoldChange <= -1, "down", "ns"))
  sym <- if (!is.null(ev_h)) setNames(ev_h$Gene, ev_h$Ensembl) else character(0)
  dh$Label <- ifelse(dh$gene %in% names(sym) & rank(dh$pvalue) <= 15, sym[dh$gene], "")
  p <- ggplot(dh, aes(log2FoldChange, -log10(padj), colour = dir)) + geom_point(size = 0.5, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = "grey50") + geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_text_repel(aes(label = Label), size = 2.1, colour = "black", max.overlaps = 30, segment.size = 0.2) + scale_colour_manual(values = COL, guide = "none") +
    labs(title = "Active vs latent TB in people with HIV (GSE37250)", subtitle = sprintf("%d up, %d down (FDR < 0.05, |log2FC| >= 1)", sum(dh$dir == "up"), sum(dh$dir == "down")),
         x = "log2 fold change", y = "-log10 adjusted p") + theme_ms
  save_gg(p, "Fig7A_TBHIV_volcano")
}
xh <- file.path(ROOT, "GSE37250", "04_normalized", "GSE37250_vst_normalized_matrix.csv")
mh <- rd("GSE37250", "02_metadata", "GSE37250_sample_metadata.csv")
if (file.exists(xh) && !is.null(mh)) {
  x <- as.matrix(utils::read.csv(xh, row.names = 1, check.names = FALSE)); mh <- mh[match(colnames(x), mh$sample_id), ]
  keep <- mh$group %in% c("ATB_HIVpos", "LTBI_HIVpos"); v <- x[, keep]; v <- v[order(-apply(v, 1, stats::var))[1:2000], ]
  pc <- stats::prcomp(t(v), scale. = FALSE); ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], Group = ifelse(mh$group[keep] == "ATB_HIVpos", "Active TB", "Latent TB"), Country = mh$country[keep])
  p <- ggplot(df, aes(PC1, PC2, colour = Group, shape = Country)) + geom_point(size = 1, alpha = 0.8) +
    scale_colour_manual(values = c("Active TB" = COL[["up"]], "Latent TB" = COL[["down"]])) +
    labs(title = "People with HIV (GSE37250): 2,000 most variable genes", x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) + theme_ms
  save_gg(p, "Fig7B_TBHIV_PCA", 95, 75)
}
## ---- 7C interferon in the interaction models ----
idm <- rd("GSE114192/TB_DM", "sensitivity_site", "GSE114192_interaction_site_sensitivity.csv")
ihv <- rd("GSE37250", "06_deg_results", "GSE37250_interaction_TBxHIV_GSEA_HALLMARK.csv")
wa  <- rd("GSE37250", "06_deg_results", "GSE37250_HIV_within_ATB_GSEA_HALLMARK.csv")
wl  <- rd("GSE37250", "06_deg_results", "GSE37250_HIV_within_LTBI_GSEA_HALLMARK.csv")
rows <- list()
if (!is.null(idm)) for (i in seq_len(nrow(idm))) rows[[length(rows) + 1]] <- data.frame(
  Analysis = paste("TB x diabetes:", sub("^I[0-9] ", "", idm$Design[i])), Set = c("IFN-gamma", "IFN-alpha"), NES = c(idm$IFNG_NES[i], idm$IFNA_NES[i]), FDR = c(idm$IFNG_padj[i], idm$IFNA_padj[i]))
add_g <- function(g, lab) if (!is.null(g)) { s <- g[g$ID %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"), ]
  rows[[length(rows) + 1]] <<- data.frame(Analysis = lab, Set = ifelse(grepl("GAMMA", s$ID), "IFN-gamma", "IFN-alpha"), NES = s$NES, FDR = s$p.adjust) }
add_g(ihv, "TB x HIV: interaction"); add_g(wa, "HIV+ vs HIV-, within active TB"); add_g(wl, "HIV+ vs HIV-, within latent TB")
if (length(rows)) {
  d <- do.call(rbind, rows); d$Analysis <- factor(d$Analysis, levels = rev(unique(d$Analysis)))
  p <- ggplot(d, aes(NES, Analysis, fill = NES)) + geom_col(width = 0.7) + facet_wrap(~ Set) + geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_gradient2(low = COL["down"], high = COL["up"], mid = "white", guide = "none") +
    labs(title = "Interferon response gene sets in the comorbidity models", subtitle = "All shown NES have FDR < 0.05 unless the bar is pale",
         x = "Normalised enrichment score", y = NULL) + theme_ms + aes(alpha = FDR < 0.05) + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none")
  save_gg(p, "Fig7C_interferon_comorbidity_models", 170, 80)
}

## ---- 8A-8C comparison across contexts ----
message("== comparison across contexts")
gc <- rd("comparisons", "contexts", "context_gene_comparison_and_rewiring.csv")
dc <- rd("comparisons", "contexts", "context_DEG_comparison.csv")
if (!is.null(gc) && !is.null(dc)) {
  core <- gc[grepl("^core", gc$Class), c("Gene", "Ensembl")]
  m <- merge(core, dc, by = "Ensembl"); mat <- as.matrix(m[, c("log2FC_TB", "log2FC_TB_DM", "log2FC_TB_HIV")]); rownames(mat) <- m$Gene
  colnames(mat) <- c("TB", "TB with diabetes", "TB with HIV"); lim <- max(abs(mat), na.rm = TRUE)
  save_base(function() pheatmap::pheatmap(mat, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize = 8, fontsize_number = 6.5,
                                          border_color = NA, color = grDevices::colorRampPalette(c(COL["down"], "white", COL["up"]))(51), breaks = seq(-lim, lim, length.out = 52),
                                          main = "Core genes: log2 fold change by context"), "Fig8A_core_genes_across_contexts", 110, 95)
}
hp <- rd("comparisons", "contexts", "context_hallmark_pathways.csv")
if (!is.null(hp)) {
  s <- hp[hp$Shared_same_sign %in% TRUE, ]; mat <- as.matrix(s[, c("NES_TB", "NES_TB_DM", "NES_TB_HIV")])
  rownames(mat) <- tools::toTitleCase(tolower(gsub("_", " ", sub("HALLMARK_", "", s$ID)))); colnames(mat) <- c("TB", "TB with diabetes", "TB with HIV")
  mat <- mat[order(-mat[, 1]), , drop = FALSE]
  save_base(function() pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize = 7, fontsize_number = 6,
                                          border_color = NA, color = grDevices::colorRampPalette(c(COL["down"], "white", COL["up"]))(51), breaks = seq(-3, 3, length.out = 52),
                                          main = "Hallmark sets significant in all contexts (NES)"), "Fig8B_hallmark_across_contexts", 130, 115)
}
mc <- rd("comparisons", "contexts", "context_module_correspondence.csv")
if (!is.null(mc)) {
  for (pr in unique(paste(mc$Context_A, mc$Context_B, sep = "|"))) {
    s <- mc[paste(mc$Context_A, mc$Context_B, sep = "|") == pr, ]; a <- strsplit(pr, "\\|")[[1]]
    s$A <- paste0(s$Module_A, " (", gsub("_", " ", s$Name_A), ")"); s$B <- paste0(s$Module_B, " (", gsub("_", " ", s$Name_B), ")")
    p <- ggplot(s, aes(B, A, fill = log2(pmax(Enrichment, 0.05)))) + geom_tile(colour = "white") + geom_text(aes(label = Shared_genes), size = 2.4) +
      scale_fill_gradient2(low = COL["down"], high = COL["up"], mid = "white", name = "log2 obs/exp") +
      labs(title = paste("Shared genes between selected modules:", a[1], "vs", a[2]), x = a[2], y = a[1]) + theme_ms + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_gg(p, paste0("Fig8C_module_correspondence_", a[1], "_vs_", a[2]), 130, 75)
  }
}

## ---- 9A-9C supporting omics ----
message("== supporting omics")
mac <- rd("external_results", "GSE165709", "integration", "GSE165709_gene_integration.csv")
ev1 <- rd("GSE114192", "20_evidence_families", "GSE114192_Gene_evidence_families.csv")
if (!is.null(mac) && !is.null(ev1)) {
  hg <- ev1$Gene[ev1$Evidence_tier == "High"]; d <- mac[mac$Gene %in% hg, ]
  d$Status <- ifelse(d$RNA_responds_HC & d$RNA_responds_PLWH, "Both", ifelse(d$RNA_responds_HC, "Healthy donors only", ifelse(d$RNA_responds_PLWH, "HIV only", "Neither")))
  p <- ggplot(d, aes(RNA_log2FC_HC, RNA_log2FC_PLWH, colour = Status)) + geom_abline(linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) + geom_point(size = 1.2, alpha = 0.85) +
    geom_text_repel(aes(label = ifelse(Gene %in% c("STAT1", "GBP1", "GBP5", "CASP1", "AIM2", "IDO1", "CXCL10", "BATF2", "IRF7"), Gene, "")), size = 2.2, colour = "black", segment.size = 0.2) +
    scale_colour_manual(values = c(Both = COL[["up"]], `Healthy donors only` = "#F4A582", `HIV only` = "#92C5DE", Neither = "grey70")) +
    labs(title = "Response of High-priority TB genes to M. tuberculosis in alveolar macrophages",
         x = "log2 fold change, infected vs uninfected (healthy donors)", y = "log2 fold change (people with HIV)", colour = "Significant in") + theme_ms
  save_gg(p, "Fig9A_macrophage_infection_response", 120, 90)
}
aa <- rd("GSE114192", "21_metabolomics", "progression", "GC6-74_amino_acids_and_kynurenine.csv")
if (!is.null(aa)) {
  d <- aa[aa$Study %in% c("ST000974", "ST000975") & aa$Metabolite %in% c("tryptophan", "kynurenine", "kynurenine/tryptophan ratio"), ]
  d$Window <- c("A baseline" = "Enrolment", "B closeness to diagnosis" = NA)[d$Analysis]
  d$Window[d$Contrast == "far_vs_control"] <- ">= 5 months before diagnosis"; d$Window[d$Contrast == "near_vs_control"] <- "< 5 months before diagnosis"
  d <- d[!is.na(d$Window), ]; d$Window <- factor(d$Window, levels = c("Enrolment", ">= 5 months before diagnosis", "< 5 months before diagnosis"))
  d$Study <- c(ST000974 = "Plasma (ST000974)", ST000975 = "Serum (ST000975)")[d$Study]
  d$star <- ifelse(d$padj < 0.05, "*", ifelse(d$pvalue < 0.05, "\u00b7", ""))
  p <- ggplot(d, aes(Window, log2FC, fill = Metabolite)) + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = star, y = log2FC + sign(log2FC) * 0.04), position = position_dodge(0.8), size = 3) + geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_wrap(~ Study) + scale_fill_manual(values = c(tryptophan = COL[["down"]], kynurenine = "#F4A582", `kynurenine/tryptophan ratio` = COL[["up"]])) +
    labs(title = "Future TB cases vs non-progressors (GC6-74)", subtitle = "* FDR < 0.05; \u00b7 nominal p < 0.05", x = NULL, y = "log2 fold change") +
    theme_ms + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
  save_gg(p, "Fig9B_kynurenine_tryptophan", 170, 90)
}
for (k in names(CTX)) {
  pr <- rd(CTX[[k]]$run, "18_singlecell", paste0(CTX[[k]]$acc, "_module_celltype_profile.csv"))
  if (is.null(pr)) next
  mat <- as.matrix(pr[, -(1:2)]); rownames(mat) <- paste0(pr$Group, " (n=", pr$N_genes, ")")
  save_base(function() pheatmap::pheatmap(mat, cluster_rows = FALSE, display_numbers = TRUE, number_format = "%.0f", fontsize = 7, fontsize_number = 5.5, border_color = NA,
                                          color = grDevices::colorRampPalette(c("white", "#FDDBC7", COL["up"]))(50),
                                          main = paste0("Module-gene expression by cell type (%)\n", CTX[[k]]$label, "; even spread = 9%")),
            paste0("Fig9C_celltype_origin_", k), 170, 75)
}
message(">>> DONE. Panels in ", OUT)
