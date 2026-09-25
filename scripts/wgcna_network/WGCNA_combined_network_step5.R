## =============================================================
## STEP 6 of Ratul's roadmap: ONE combined TB PPI network
##
## "After analyzing each module separately, combine the important
##  genes from all significant modules into one TB PPI network.
##  Keep the module identity of every gene." -- used later for
##  RWR, network proximity, and diffusion.
##
## Re-queries STRING on the UNION of all modules' genes (not just
## pooling the 3 separate per-module edge lists from step 5) so
## cross-module ("bridge") interactions are actually captured --
## those are exactly what a combined network is for.
## =============================================================

required_pkgs <- c("httr", "igraph")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "), "))")
}

## ---- CONFIG (keep consistent with step 5) ----
ACC             <- "GSE114192"
SPECIES         <- 9606
REQUIRED_SCORE  <- 400
USE_FULL_MODULE <- TRUE
N_LABELS_FULL   <- 20
MODULE_COLORS   <- c(green = "#2ecc71", purple = "#9b59b6", blue = "#3498db")  # plot colors, edit if you rename blue

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
BASE_DIR   <- file.path(SCRIPT_DIR, "..", ACC)
dir_wgcna  <- file.path(BASE_DIR, "10_wgcna_results")
dir_string <- file.path(BASE_DIR, "12_string_ppi")
dir.create(dir_string, recursive = TRUE, showWarnings = FALSE)

## ---- load module gene lists (same logic/choice as step 5) ----
if (USE_FULL_MODULE) {
  all_genes   <- read.csv(file.path(dir_wgcna, paste0(ACC, "_Gene_Module_Assignment.csv")), stringsAsFactors = FALSE)
  mod_summary <- read.csv(file.path(dir_wgcna, paste0(ACC, "_ModuleTrait_Correlation.csv")), stringsAsFactors = FALSE)
  sel <- mod_summary$Module[mod_summary$Selected == TRUE]
  gene_lists <- setNames(lapply(sel, function(m) all_genes$Gene[all_genes$ModuleColor == m]), sel)
} else {
  dg <- read.csv(file.path(dir_wgcna, paste0(ACC, "_DiseaseGenes_DEG_WGCNA.csv")), stringsAsFactors = FALSE)
  gene_lists <- split(dg$Gene, dg$Module)
}

have_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)
ens_to_symbol <- function(ens_ids) {
  if (have_orgdb) {
    m <- suppressMessages(clusterProfiler::bitr(ens_ids, fromType = "ENSEMBL", toType = "SYMBOL",
                                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    return(unique(na.omit(m$SYMBOL)))
  }
  unique(ens_ids)
}

## ---- tag every gene symbol with its module of origin ----
symbol_module_map <- do.call(rbind, lapply(names(gene_lists), function(m) {
  syms <- ens_to_symbol(gene_lists[[m]])
  data.frame(Symbol = syms, Module = m, stringsAsFactors = FALSE)
}))
dupes <- symbol_module_map$Symbol[duplicated(symbol_module_map$Symbol)]
if (length(dupes) > 0) {
  cat("!! ", length(dupes), "symbol(s) mapped to >1 module (ID-mapping collision) -- keeping first occurrence:",
      paste(unique(dupes), collapse = ", "), "\n")
  symbol_module_map <- symbol_module_map[!duplicated(symbol_module_map$Symbol), ]
}
cat(">>> Combined set:", nrow(symbol_module_map), "genes across", length(unique(symbol_module_map$Module)),
    "modules (", paste(table(symbol_module_map$Module), names(table(symbol_module_map$Module)), collapse = ", "), ")\n")

## ---- STRING REST API helpers (same as step 5) ----
string_get_ids <- function(symbols) {
  r <- httr::POST("https://string-db.org/api/tsv/get_string_ids",
                   body = list(identifiers = paste(symbols, collapse = "\r"),
                               species = SPECIES, limit = 1, echo_query = 1,
                               caller_identity = "negin_tb_thesis_wgcna"),
                   encode = "form", httr::timeout(180))
  httr::stop_for_status(r)
  read.delim(text = httr::content(r, "text", encoding = "UTF-8"), stringsAsFactors = FALSE)
}
string_get_network <- function(string_ids) {
  r <- httr::POST("https://string-db.org/api/tsv/network",
                   body = list(identifiers = paste(string_ids, collapse = "\r"),
                               species = SPECIES, required_score = REQUIRED_SCORE,
                               caller_identity = "negin_tb_thesis_wgcna"),
                   encode = "form", httr::timeout(240))
  httr::stop_for_status(r)
  read.delim(text = httr::content(r, "text", encoding = "UTF-8"), stringsAsFactors = FALSE)
}

## ---- one fresh STRING query on the COMBINED gene set ----
## (this is what actually reveals cross-module "bridge" interactions --
##  simply pooling step 5's 3 separate network CSVs would miss those entirely)
cat("\n>>> Querying STRING on the combined set (may take longer than any single module)...\n")
ids_df <- string_get_ids(symbol_module_map$Symbol)
cat("  resolved", length(unique(ids_df$stringId)), "/", nrow(symbol_module_map), "to STRING IDs\n")

## Authoritative StringID -> Module map, built from OUR OWN query (queryItem = the exact
## symbol we submitted, already tied to a module in symbol_module_map) -- never by matching
## STRING's own preferredName afterward. (queryItem and STRING's preferredName can legitimately
## differ for the same protein, e.g. via a synonym -- that mismatch is what silently dropped
## DDX58 and 9 other genes to "unmapped" last run, even though DDX58 was genuinely in green.)
id_module_map <- merge(ids_df[, c("queryItem", "stringId")], symbol_module_map,
                        by.x = "queryItem", by.y = "Symbol")
id_module_map <- id_module_map[!duplicated(id_module_map$stringId), ]

net_df <- string_get_network(unique(ids_df$stringId))
combined_net_file <- file.path(dir_string, paste0(ACC, "_COMBINED_STRING_network.csv"))
write.csv(net_df, combined_net_file, row.names = FALSE)
cat("  ", nrow(net_df), "interactions ->", combined_net_file, "\n")

## ---- build graph using STRING IDs as node identity (stable, guaranteed to match id_module_map) ----
g <- igraph::simplify(igraph::graph_from_data_frame(
  net_df[, c("stringId_A", "stringId_B")], directed = FALSE))

mod_lookup  <- setNames(id_module_map$Module, id_module_map$stringId)
name_lookup <- setNames(c(net_df$preferredName_A, net_df$preferredName_B),
                         c(net_df$stringId_A, net_df$stringId_B))
name_lookup <- name_lookup[!duplicated(names(name_lookup))]

igraph::V(g)$module      <- mod_lookup[igraph::V(g)$name]
igraph::V(g)$displayName <- ifelse(igraph::V(g)$name %in% names(name_lookup),
                                    name_lookup[igraph::V(g)$name], igraph::V(g)$name)
n_unmapped <- sum(is.na(igraph::V(g)$module))
if (n_unmapped > 0) cat("  !!", n_unmapped, "node(s) still unmapped after the fix -- inspect these:",
                         paste(utils::head(igraph::V(g)$displayName[is.na(igraph::V(g)$module)], 10), collapse = ", "), "\n")
igraph::V(g)$module[is.na(igraph::V(g)$module)] <- "unmapped"

## ---- flag cross-module ("bridge") edges: these are new information step 5 couldn't see ----
em      <- igraph::ends(g, igraph::E(g), names = TRUE)  # STRING IDs (stable join key)
mod_a   <- igraph::V(g)$module[match(em[, 1], igraph::V(g)$name)]
mod_b   <- igraph::V(g)$module[match(em[, 2], igraph::V(g)$name)]
is_bridge <- mod_a != mod_b
cat("\n>>> Cross-module bridge edges:", sum(is_bridge), "/", length(is_bridge),
    "(", round(100 * mean(is_bridge), 1), "% )\n")
bridge_table <- data.frame(GeneA = name_lookup[em[is_bridge, 1]], ModuleA = mod_a[is_bridge],
                            GeneB = name_lookup[em[is_bridge, 2]], ModuleB = mod_b[is_bridge])
bridge_file <- file.path(dir_string, paste0(ACC, "_COMBINED_bridge_edges.csv"))
write.csv(bridge_table, bridge_file, row.names = FALSE)
cat("  ->", bridge_file, "\n")

## ---- centrality on the COMBINED network (same 4 measures, now network-wide) ----
cent <- data.frame(
  Gene        = igraph::V(g)$displayName,
  Module      = igraph::V(g)$module,
  Degree      = igraph::degree(g),
  Betweenness = igraph::betweenness(g, normalized = TRUE),
  Closeness   = igraph::closeness(g, normalized = TRUE),
  Eigenvector = igraph::eigen_centrality(g)$vector
)
cent <- cent[order(-cent$Degree, -cent$Eigenvector), ]
cent_file <- file.path(dir_string, paste0(ACC, "_COMBINED_centrality.csv"))
write.csv(cent, cent_file, row.names = FALSE)
cat("\ntop 10 hubs overall (any module):\n")
print(utils::head(cent, 10))
cat("\ntop 5 hubs WITH THE MOST cross-module bridges (interesting for network medicine):\n")
bridge_deg <- table(c(bridge_table$GeneA, bridge_table$GeneB))
if (length(bridge_deg) > 0) print(utils::head(sort(bridge_deg, decreasing = TRUE), 5))

## ---- plot: nodes colored by module of origin ----
tryCatch({
  top_labels <- utils::head(cent$Gene, N_LABELS_FULL)  # already displayName
  igraph::V(g)$size        <- 2 + 8 * sqrt(igraph::degree(g) / max(igraph::degree(g)))
  igraph::V(g)$label       <- ifelse(igraph::V(g)$displayName %in% top_labels, igraph::V(g)$displayName, NA)
  igraph::V(g)$label.cex   <- 1.1
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$color       <- MODULE_COLORS[igraph::V(g)$module]
  igraph::V(g)$color[is.na(igraph::V(g)$color)] <- "grey50"
  igraph::V(g)$frame.color <- NA
  igraph::E(g)$color <- ifelse(is_bridge, grDevices::adjustcolor("red", alpha.f = 0.5),
                                grDevices::adjustcolor("grey75", alpha.f = 0.25))
  igraph::E(g)$width <- ifelse(is_bridge, 1.2, 0.3)

  set.seed(42)
  layout_combined <- igraph::layout_with_fr(g, niter = 3000)  # this graph is the biggest one, needs the most spreading-out

  centroid <- colMeans(layout_combined)
  ang <- atan2(layout_combined[, 2] - centroid[2], layout_combined[, 1] - centroid[1])
  igraph::V(g)$label.degree <- ang
  igraph::V(g)$label.dist   <- ifelse(!is.na(igraph::V(g)$label), 1.6, 0)

  plot_file <- file.path(dir_string, paste0(ACC, "_COMBINED_STRING_network.pdf"))
  grDevices::pdf(plot_file, width = 22, height = 22, pointsize = 16)
  plot(g, layout = layout_combined, vertex.label.font = 2,
       main = paste0("Combined TB PPI network (", ACC, ") -- colored by module, red edges = cross-module bridges"))
  graphics::legend("bottomleft", legend = names(MODULE_COLORS), col = MODULE_COLORS, pch = 19, bty = "n", cex = 1.5)
  grDevices::dev.off()
  cat("\ncombined network plot ->", plot_file, "\n")
}, error = function(e) cat("!! plotting failed:", conditionMessage(e), "(tables above are still saved fine)\n"))

cat("\n>>> DONE. Combined network + bridge edges + centrality in", dir_string, "\n")
cat(">>> Next (Ratul step 7): prepare the TB_UP / TB_DOWN expression signature for iLINCS.\n")
