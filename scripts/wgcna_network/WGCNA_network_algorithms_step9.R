## STEP 9: RWR, diffusion, network proximity, community detection and module-overlap scoring
## on the combined TB network. RWR uses a column-normalized transition matrix, diffusion a
## symmetric-normalized one (Vanunu et al. 2010).

required_pkgs <- c("igraph")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing package: igraph. Install with: install.packages(\"igraph\")")

## ---- CONFIG ----
RESTART_R    <- 0.5     # RWR restart probability (0.3-0.5 is standard)
DIFFUSE_A    <- 0.5     # diffusion keep-vs-spread balance
MAX_ITER     <- 200
TOLERANCE    <- 1e-10
N_PERMUTATIONS <- 1000  # for the network-proximity null distribution

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
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found
OUT_DIR <- run_out_dir(BASE_DIR)
dir_string <- file.path(OUT_DIR, "12_string_ppi")
dir_wgcna  <- file.path(OUT_DIR, "10_wgcna_results")
dir_dt     <- file.path(OUT_DIR, "14_drug_targets")
dir_out    <- file.path(OUT_DIR, "15_network_algorithms")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## ---- rebuild the combined network graph ----
net_df  <- read.csv(file.path(dir_string, paste0(ACC, "_COMBINED_STRING_network.csv")), stringsAsFactors = FALSE)
cent_df <- read.csv(file.path(dir_string, paste0(ACC, "_COMBINED_centrality.csv")), stringsAsFactors = FALSE)  # Gene(=displayName), Module
dg      <- read.csv(file.path(dir_wgcna, paste0(ACC, "_DiseaseGenes_DEG_WGCNA.csv")), stringsAsFactors = FALSE)  # 188 seed genes (Ensembl IDs)
dt      <- read.csv(file.path(dir_dt, paste0(ACC, "_Drug_Target_TBnetwork.csv")), stringsAsFactors = FALSE)      # from step 8

g <- igraph::simplify(igraph::graph_from_data_frame(net_df[, c("stringId_A", "stringId_B")], directed = FALSE))
name_lookup <- setNames(c(net_df$preferredName_A, net_df$preferredName_B),
                         c(net_df$stringId_A, net_df$stringId_B))
name_lookup <- name_lookup[!duplicated(names(name_lookup))]
igraph::V(g)$displayName <- ifelse(igraph::V(g)$name %in% names(name_lookup),
                                    name_lookup[igraph::V(g)$name], igraph::V(g)$name)
n_nodes <- igraph::vcount(g)
cat(">>> Network:", n_nodes, "nodes,", igraph::ecount(g), "edges.\n")

## ---- seed genes = the 188 disease genes ----
## matched to network nodes by name: exact symbol first, alias as fallback (e.g. RIGI = DDX58)
ensembl_to_network_name <- function(ens_ids, network_names, annot = NULL) {
  ens_ids <- unique(ens_ids)
  if (is.null(annot)) {
    annot <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ens_ids,
                                                     keytype = "ENSEMBL", columns = c("SYMBOL", "ALIAS")))
  }
  sym_by_ens   <- split(annot$SYMBOL, annot$ENSEMBL)
  alias_by_ens <- split(annot$ALIAS,  annot$ENSEMBL)
  syms    <- lapply(ens_ids, function(e) unique(stats::na.omit(sym_by_ens[[e]])))
  aliases <- lapply(ens_ids, function(e) unique(stats::na.omit(alias_by_ens[[e]])))
  first_hit <- function(x) { h <- intersect(x, network_names); if (length(h)) h[1] else NA_character_ }
  exact  <- vapply(syms,    first_hit, character(1))   # 1) exact SYMBOL match
  viaal  <- vapply(aliases, first_hit, character(1))   # 2) ALIAS match (fallback only)
  nm <- exact
  use_alias <- is.na(exact) & !is.na(viaal)
  ## an alias hit is used only if the node is not claimed by another gene's symbol or by 2+ aliases
  claimed      <- unique(stats::na.omit(exact))
  alias_counts <- table(viaal[use_alias])
  unique_alias <- names(alias_counts)[alias_counts == 1]
  ok <- use_alias & !(viaal %in% claimed) & (viaal %in% unique_alias)
  nm[ok] <- viaal[ok]
  first_sym <- vapply(syms, function(s) if (length(s)) s[1] else NA_character_, character(1))
  data.frame(ENSEMBL = ens_ids, SYMBOL = unname(first_sym), NetName = unname(nm), stringsAsFactors = FALSE)
}
if (!(requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE))) {
  stop("Need org.Hs.eg.db + AnnotationDbi.")
}
seed_map   <- ensembl_to_network_name(dg$Gene, igraph::V(g)$displayName)
seed_names <- unique(stats::na.omit(seed_map$NetName))
seed_idx   <- which(igraph::V(g)$displayName %in% seed_names)
cat(">>> Seed (disease) genes matched to a network node:", sum(!is.na(seed_map$NetName)), "/", nrow(seed_map),
    "(the rest simply aren't in the STRING-derived network)\n")
if (length(seed_idx) == 0) stop("No seed genes matched the network -- check gene ID/symbol consistency before continuing.")

adj <- igraph::as_adjacency_matrix(g, sparse = FALSE)
node_names <- igraph::V(g)$name  # STRING IDs, stable row/col order for everything below

## ---- 1. RANDOM WALK WITH RESTART (RWR) ----
cat("\n========== 1. RWR ==========\n")
deg_vec <- igraph::degree(g)
W <- sweep(adj, 2, pmax(colSums(adj), 1e-12), "/")  # column-normalize -> transition matrix

p0 <- rep(0, n_nodes); p0[seed_idx] <- 1 / length(seed_idx)
p <- p0
for (i in seq_len(MAX_ITER)) {
  p_new <- (1 - RESTART_R) * as.vector(W %*% p) + RESTART_R * p0
  if (sum(abs(p_new - p)) < TOLERANCE) { p <- p_new; break }
  p <- p_new
}
cat(">>> RWR converged after", i, "iterations.\n")
rwr_score <- setNames(p, node_names)

## ---- 2. Diffusion (symmetric normalization) ----
cat("\n========== 2. Diffusion ==========\n")
d_sqrt_inv <- 1 / sqrt(pmax(deg_vec, 1e-12))
W_sym <- sweep(sweep(adj, 1, d_sqrt_inv, "*"), 2, d_sqrt_inv, "*")  # D^-1/2 A D^-1/2

p0b <- p0  # same seed vector
q <- p0b
for (i in seq_len(MAX_ITER)) {
  q_new <- DIFFUSE_A * as.vector(W_sym %*% q) + (1 - DIFFUSE_A) * p0b
  if (sum(abs(q_new - q)) < TOLERANCE) { q <- q_new; break }
  q <- q_new
}
cat(">>> Diffusion converged after", i, "iterations.\n")
diffusion_score <- setNames(q, node_names)

## ---- 3. Network proximity: closest distance from drug targets to seeds (Guney et al. 2016) ----
## unreachable nodes excluded; degree-matched null; empirical p-value reported with the Z-score
cat("\n========== 3. Network Proximity (per drug) ==========\n")
dist_mat <- igraph::distances(g)  # full shortest-path matrix, reused for every drug (compute once)
N_DEGREE_BINS <- 20

nearest_seed_dist <- apply(dist_mat[, seed_idx, drop = FALSE], 1, min)
reachable <- is.finite(nearest_seed_dist)
cat(">>>", sum(!reachable), "of", n_nodes, "nodes have no path to any seed (small disconnected components) -- excluded from proximity.\n")

closest_distance <- function(target_idx, seed_idx) {
  if (length(target_idx) == 0) return(NA_real_)
  d <- apply(dist_mat[target_idx, seed_idx, drop = FALSE], 1, min)
  d <- d[is.finite(d)]
  if (length(d) == 0) return(NA_real_)
  mean(d)
}

## degree bins cut on degree values, so equal-degree nodes share a bin
deg_breaks <- unique(stats::quantile(deg_vec, probs = seq(0, 1, length.out = N_DEGREE_BINS + 1), type = 1))
bin_id <- as.integer(cut(deg_vec, breaks = deg_breaks, include.lowest = TRUE, labels = FALSE))
pool_by_bin <- split(which(reachable), bin_id[reachable])
sample_degree_matched <- function(target_idx) {
  vapply(target_idx, function(t) {
    p <- pool_by_bin[[as.character(bin_id[t])]]
    if (is.null(p) || length(p) == 0) p <- which(reachable)
    p[sample.int(length(p), 1)]
  }, integer(1))
}

set.seed(42)
proximity_results <- list()
drugs <- unique(dt$Drug)
for (drug in drugs) {
  targets    <- unique(dt$Drug_Target[dt$Drug == drug])
  target_idx <- which(igraph::V(g)$displayName %in% targets)
  if (length(target_idx) == 0) next
  reach_idx  <- target_idx[reachable[target_idx]]
  if (length(reach_idx) == 0) {
    proximity_results[[drug]] <- data.frame(Drug = drug, N_targets_in_network = length(target_idx),
        N_targets_reachable = 0, Observed_distance = NA_real_, Null_mean = NA_real_, Null_sd = NA_real_,
        Z_score = NA_real_, P_empirical = NA_real_)
    next
  }
  d_obs  <- closest_distance(reach_idx, seed_idx)
  null_d <- vapply(seq_len(N_PERMUTATIONS), function(i) closest_distance(sample_degree_matched(reach_idx), seed_idx), numeric(1))
  null_sd <- stats::sd(null_d, na.rm = TRUE)
  z <- if (is.finite(null_sd) && null_sd > 0) (d_obs - mean(null_d, na.rm = TRUE)) / null_sd else NA_real_
  p_emp <- (sum(null_d <= d_obs, na.rm = TRUE) + 1) / (sum(!is.na(null_d)) + 1)  # smaller = closer than chance

  proximity_results[[drug]] <- data.frame(Drug = drug, N_targets_in_network = length(target_idx),
      N_targets_reachable = length(reach_idx), Observed_distance = d_obs,
      Null_mean = mean(null_d, na.rm = TRUE), Null_sd = null_sd, Z_score = z, P_empirical = p_emp)
}
proximity_table <- do.call(rbind, proximity_results)
proximity_table <- proximity_table[order(proximity_table$Z_score), ]  # most negative = closest to disease = best
write.csv(proximity_table, file.path(dir_out, paste0(ACC, "_NetworkProximity.csv")), row.names = FALSE)
cat(">>> Network proximity computed for", nrow(proximity_table), "/", length(drugs), "drugs;",
    sum(is.finite(proximity_table$Z_score)), "have a finite Z-score.\n")
cat(">>> Most disease-proximal (most negative Z, i.e. closer than degree-matched random targets would be):\n")
print(utils::head(proximity_table, 5))

## ---- 4. COMMUNITY DETECTION ----
cat("\n========== 4. Community Detection ==========\n")
communities <- igraph::cluster_louvain(g)
comm_table <- data.frame(Gene = igraph::V(g)$displayName, Community = igraph::membership(communities),
                          Module = cent_df$Module[match(igraph::V(g)$displayName, cent_df$Gene)])
write.csv(comm_table, file.path(dir_out, paste0(ACC, "_Communities.csv")), row.names = FALSE)
cat(">>> Found", length(unique(comm_table$Community)), "communities. Modularity:", round(igraph::modularity(communities), 3), "\n")
cat(">>> Community x WGCNA-module cross-tab:\n")
print(table(comm_table$Community, comm_table$Module, useNA = "ifany"))

## ---- 5. Module-overlap scoring (hypergeometric, per drug and module) ----
cat("\n========== 5. Module-Overlap Scoring ==========\n")
module_sizes <- table(cent_df$Module)
overlap_rows <- list()
for (drug in drugs) {
  d_sub <- dt[dt$Drug == drug, ]
  n_drug_targets <- length(unique(d_sub$Drug_Target))
  for (mod in names(module_sizes)) {
    n_hit <- length(unique(d_sub$Drug_Target[d_sub$TB_Module == mod]))
    if (n_hit == 0) next
    # P(>= n_hit targets in the module) under the hypergeometric null
    p_val <- stats::phyper(n_hit - 1, module_sizes[[mod]], n_nodes - module_sizes[[mod]], n_drug_targets, lower.tail = FALSE)
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      Drug = drug, Module = mod, N_targets_in_module = n_hit,
      N_drug_targets_total = n_drug_targets, Module_size = module_sizes[[mod]], P_value = p_val
    )
  }
}
overlap_table <- do.call(rbind, overlap_rows)
overlap_table <- overlap_table[order(overlap_table$P_value), ]
write.csv(overlap_table, file.path(dir_out, paste0(ACC, "_ModuleOverlapScoring.csv")), row.names = FALSE)
cat(">>> Module-overlap scored for", length(unique(overlap_table$Drug)), "drugs across", length(module_sizes), "modules.\n")

## ---- per-gene RWR + diffusion scores ----
gene_scores <- data.frame(Gene = igraph::V(g)$displayName, Module = cent_df$Module[match(igraph::V(g)$displayName, cent_df$Gene)],
                           RWR_score = rwr_score[node_names], Diffusion_score = diffusion_score[node_names])
gene_scores <- gene_scores[order(-gene_scores$RWR_score), ]
write.csv(gene_scores, file.path(dir_out, paste0(ACC, "_GeneScores_RWR_Diffusion.csv")), row.names = FALSE)

cat("\n>>> DONE. All 5 outputs in", dir_out, ":\n")
cat("  -", paste0(ACC, "_GeneScores_RWR_Diffusion.csv"))
cat("\n  -", paste0(ACC, "_NetworkProximity.csv"), "(per drug)")
cat("\n  -", paste0(ACC, "_Communities.csv"))
cat("\n  -", paste0(ACC, "_ModuleOverlapScoring.csv"), "(per drug x module)")
