## =============================================================

## Uses STRING's REST API directly (https://string-db.org/help/api/) --
## no manual upload needed, no STRINGdb Bioconductor package required.
## =============================================================

# install.packages(c("httr", "igraph"))   # both plain CRAN, no Bioconductor

required_pkgs <- c("httr", "igraph")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "), "))")
}

## ---- CONFIG ----
ACC             <- "GSE114192"
SPECIES         <- 9606     # human
REQUIRED_SCORE  <- 400      # STRING "medium confidence" default (0-1000); use 700+ for high-confidence only
USE_FULL_MODULE <- TRUE     # TRUE  = every gene in the module (what Ratul's wording asks for)
                             # FALSE = only the step-3 "disease genes" (DEG + module overlap) subset
N_LABELS_FULL   <- 15        # how many top hubs to label on the full-network plot (rest stay unlabeled dots)
N_HUBS_SUBNET   <- 30        # size of the separate, fully-labeled "hub-only" zoomed-in plot

## ---- locate script dir  ----
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
BASE_DIR   <- file.path(SCRIPT_DIR, "..", ACC)
dir_wgcna  <- file.path(BASE_DIR, "10_wgcna_results")
dir_string <- file.path(BASE_DIR, "12_string_ppi")
dir.create(dir_string, recursive = TRUE, showWarnings = FALSE)

## ---- load module gene lists ----
if (USE_FULL_MODULE) {
  all_genes   <- read.csv(file.path(dir_wgcna, paste0(ACC, "_Gene_Module_Assignment.csv")), stringsAsFactors = FALSE)
  mod_summary <- read.csv(file.path(dir_wgcna, paste0(ACC, "_ModuleTrait_Correlation.csv")), stringsAsFactors = FALSE)
  sel <- mod_summary$Module[mod_summary$Selected == TRUE]
  gene_lists <- setNames(lapply(sel, function(m) all_genes$Gene[all_genes$ModuleColor == m]), sel)
} else {
  dg <- read.csv(file.path(dir_wgcna, paste0(ACC, "_DiseaseGenes_DEG_WGCNA.csv")), stringsAsFactors = FALSE)
  gene_lists <- split(dg$Gene, dg$Module)
}
for (m in names(gene_lists)) cat(m, ":", length(gene_lists[[m]]), "genes\n")

## ---- ENSEMBL -> SYMBOL (STRING resolves symbols far more reliably than raw ENSG) ----
have_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)
ens_to_symbol <- function(ens_ids) {
  if (have_orgdb) {
    m <- suppressMessages(clusterProfiler::bitr(ens_ids, fromType = "ENSEMBL", toType = "SYMBOL",
                                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    return(unique(na.omit(m$SYMBOL)))
  }
  cat("  (org.Hs.eg.db/clusterProfiler not available -- sending raw Ensembl IDs to STRING, which usually still works but is less reliable)\n")
  unique(ens_ids)
}

## ---- STRING REST API helpers ----
string_get_ids <- function(symbols) {
  r <- httr::POST("https://string-db.org/api/tsv/get_string_ids",
                   body = list(identifiers = paste(symbols, collapse = "\r"),
                               species = SPECIES, limit = 1, echo_query = 1,
                               caller_identity = "negin_tb_thesis_wgcna"),
                   encode = "form", httr::timeout(120))
  httr::stop_for_status(r)
  read.delim(text = httr::content(r, "text", encoding = "UTF-8"), stringsAsFactors = FALSE)
}
string_get_network <- function(string_ids) {
  r <- httr::POST("https://string-db.org/api/tsv/network",
                   body = list(identifiers = paste(string_ids, collapse = "\r"),
                               species = SPECIES, required_score = REQUIRED_SCORE,
                               caller_identity = "negin_tb_thesis_wgcna"),
                   encode = "form", httr::timeout(180))
  httr::stop_for_status(r)
  read.delim(text = httr::content(r, "text", encoding = "UTF-8"), stringsAsFactors = FALSE)
}

## ---- run per module ----
for (mod in names(gene_lists)) {
  cat("\n====", mod, "(", length(gene_lists[[mod]]), "genes ) ====\n")
  net_file  <- file.path(dir_string, paste0(ACC, "_", mod, "_STRING_network.csv"))
  cent_file <- file.path(dir_string, paste0(ACC, "_", mod, "_centrality.csv"))

  symbols <- ens_to_symbol(gene_lists[[mod]])
  cat("  ", length(symbols), "gene symbols to query\n")

  ids_df <- tryCatch(string_get_ids(symbols), error = function(e) {
    cat("  !! get_string_ids failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(ids_df) || nrow(ids_df) == 0) { cat("  skipping (no STRING IDs resolved)\n"); next }
  cat("  get_string_ids columns:", paste(colnames(ids_df), collapse = ", "), "\n")
  string_ids <- unique(ids_df$stringId)
  cat("  resolved", length(string_ids), "/", length(symbols), "to STRING IDs\n")

  net_df <- tryCatch(string_get_network(string_ids), error = function(e) {
    cat("  !! network fetch failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(net_df) || nrow(net_df) == 0) { cat("  skipping (no interactions returned)\n"); next }
  write.csv(net_df, net_file, row.names = FALSE)
  cat("  network columns:", paste(colnames(net_df), collapse = ", "), "\n")
  cat("  ", nrow(net_df), "interactions ->", net_file, "\n")

  ## ---- build graph + the 4 centrality measures  ----
  g <- igraph::simplify(igraph::graph_from_data_frame(
    net_df[, c("preferredName_A", "preferredName_B")], directed = FALSE))

  cent <- data.frame(
    Gene        = igraph::V(g)$name,
    Degree      = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness   = igraph::closeness(g, normalized = TRUE),
    Eigenvector = igraph::eigen_centrality(g)$vector
  )
  cent <- cent[order(-cent$Degree, -cent$Eigenvector), ]
  write.csv(cent, cent_file, row.names = FALSE)

  cat("  top 5 hubs by degree:\n")
  print(utils::head(cent, 5))

  ## ---- network diagrams ----
  tryCatch({
    top_labels <- utils::head(cent$Gene, N_LABELS_FULL)

    igraph::V(g)$size        <- 2 + 8 * sqrt(igraph::degree(g) / max(igraph::degree(g)))
    igraph::V(g)$label       <- ifelse(igraph::V(g)$name %in% top_labels, igraph::V(g)$name, NA)
    igraph::V(g)$label.cex   <- 1.0
    igraph::V(g)$label.color <- "black"
    igraph::V(g)$color       <- ifelse(igraph::V(g)$name %in% top_labels, "#e74c3c", "#a9c4e0")
    igraph::V(g)$frame.color <- NA

    set.seed(42)
    layout_full <- igraph::layout_with_fr(g, niter = 2000)  # more iterations -> nodes spread out more, less crowding

    ## push each label OUTWARD, away from the graph's center, instead of centering it on
    ## the node
    centroid <- colMeans(layout_full)
    ang <- atan2(layout_full[, 2] - centroid[2], layout_full[, 1] - centroid[1])
    igraph::V(g)$label.degree <- ang
    igraph::V(g)$label.dist   <- ifelse(!is.na(igraph::V(g)$label), 1.4, 0)

    full_plot_file <- file.path(dir_string, paste0(ACC, "_", mod, "_STRING_network_full.pdf"))
    grDevices::pdf(full_plot_file, width = 18, height = 18, pointsize = 14)
    plot(g, layout = layout_full,
         edge.color = grDevices::adjustcolor("grey70", alpha.f = 0.35), edge.width = 0.4,
         vertex.label.font = 2,
         main = paste0(mod, " module STRING network (", ACC, ") -- top ", N_LABELS_FULL, " hubs labeled, n=", igraph::vcount(g)))
    grDevices::dev.off()
    cat("  network plot ->", full_plot_file, "\n")

    ## fully-labeled view of just the top hub genes + their edges among each other
    n_sub <- min(N_HUBS_SUBNET, igraph::vcount(g))
    hub_genes_sub <- utils::head(cent$Gene, n_sub)
    g_sub <- igraph::induced_subgraph(g, vids = which(igraph::V(g)$name %in% hub_genes_sub))
    igraph::V(g_sub)$label       <- igraph::V(g_sub)$name
    igraph::V(g_sub)$label.cex   <- 0.9
    igraph::V(g_sub)$label.color <- "black"
    igraph::V(g_sub)$size        <- 8 + 14 * sqrt(igraph::degree(g_sub) / max(1, max(igraph::degree(g_sub))))
    igraph::V(g_sub)$color       <- "#e74c3c"
    igraph::V(g_sub)$frame.color <- NA

    set.seed(42)
    layout_sub <- igraph::layout_with_fr(g_sub, niter = 2000)
    hub_plot_file <- file.path(dir_string, paste0(ACC, "_", mod, "_STRING_hubs_top", n_sub, ".pdf"))
    grDevices::pdf(hub_plot_file, width = 12, height = 12, pointsize = 16)
    plot(g_sub, layout = layout_sub, edge.color = "grey55", edge.width = 0.8,
         main = paste0(mod, ": top ", n_sub, " hub genes (", ACC, ")"))
    grDevices::dev.off()
    cat("  hub-subnetwork plot ->", hub_plot_file, "\n")
  }, error = function(e) cat("  !! plotting failed:", conditionMessage(e), "(tables above are still saved fine)\n"))
}

cat("\n>>> DONE. Per-module STRING interaction tables + centrality tables in", dir_string, "\n")
