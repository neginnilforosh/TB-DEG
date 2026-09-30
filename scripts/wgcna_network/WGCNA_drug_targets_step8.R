## =============================================================
## Target lookup: DGIdb v5's public GraphQL API (dgidb.org/api/graphql,
## no key needed).
## =============================================================

required_pkgs <- c("httr", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "), "))")
}

## ---- CONFIG ----
ACC <- "GSE114192"
DGIDB_URL <- "https://dgidb.org/api/graphql"

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
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found
dir_sig    <- file.path(BASE_DIR, "13_ilincs_signature")
dir_string <- file.path(BASE_DIR, "12_string_ppi")
dir_out    <- file.path(BASE_DIR, "14_drug_targets")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## ---- load candidate compounds (step 7) + combined network (step 5) ----
cand_file <- file.path(dir_sig, paste0(ACC, "_iLINCS_candidate_compounds.csv"))
net_file  <- file.path(dir_string, paste0(ACC, "_COMBINED_centrality.csv"))
if (!file.exists(cand_file)) stop("Missing: ", cand_file, " -- run WGCNA_iLINCS_query_step7.R first.")
if (!file.exists(net_file))  stop("Missing: ", net_file,  " -- run WGCNA_combined_network_step5.R first.")

candidates <- read.csv(cand_file, stringsAsFactors = FALSE)
tb_network <- read.csv(net_file, stringsAsFactors = FALSE)  # Gene, Module, Degree, Betweenness, Closeness, Eigenvector

drug_col <- intersect(c("treatment", "compound", "Compound", "Treatment", "name","Target"), names(candidates))[1]
if (is.na(drug_col)) stop("Couldn't find a compound-name column in ", cand_file,
                          " -- columns are: ", paste(names(candidates), collapse = ", "))
drug_names <- unique(candidates[[drug_col]])
cat(">>> Looking up targets for", length(drug_names), "candidate compounds via DGIdb...\n")

## ---- DGIdb GraphQL query (drug -> target genes), by symmetry with the confirmed gene-side query ----
dgidb_query <- '
query GetDrugInteractions($drugs: [String!]!) {
  drugs(names: $drugs) {
    nodes {
      name
      conceptId
      interactions {
        gene { name longName }
        interactionTypes { type }
        interactionScore
        sources { fullName }
      }
    }
  }
}'

resp <- httr::POST(DGIDB_URL, encode = "json", httr::timeout(120),
                    body = list(query = dgidb_query, variables = list(drugs = as.list(drug_names))))

if (httr::status_code(resp) != 200 || !is.null(httr::content(resp, "parsed")$errors)) {
  cat("!! DGIdb query failed or returned GraphQL errors. Response:\n")
  print(httr::content(resp, "text", encoding = "UTF-8"))
  cat("\n>>> Fetching DGIdb's own schema so you can see the real field name instead of guessing:\n")
  schema_q <- '{ __schema { queryType { fields { name } } } }'
  schema_resp <- httr::POST(DGIDB_URL, encode = "json", body = list(query = schema_q))
  print(httr::content(schema_resp, "parsed"))
  stop("Fix the query field name above (probably just 'drugs' -> whatever the schema printout shows) and re-run.")
}

parsed <- httr::content(resp, "parsed")
nodes <- parsed$data$drugs$nodes
cat(">>> DGIdb matched", length(nodes), "/", length(drug_names), "compounds.\n")

## ---- flatten to one row per Drug-Target pair ----
rows <- list()
for (n in nodes) {
  for (interaction in n$interactions) {
    rows[[length(rows) + 1]] <- data.frame(
      Drug        = n$name,
      Target_Gene = interaction$gene$name,
      Interaction_Types = paste(sapply(interaction$interactionTypes, function(x) x$type), collapse = ";"),
      Sources     = paste(sapply(interaction$sources, function(x) x$fullName), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
drug_target_table <- if (length(rows)) do.call(rbind, rows) else
  data.frame(Drug=character(), Target_Gene=character(), Interaction_Types=character(), Sources=character())
cat(">>> Total Drug-Target pairs from DGIdb:", nrow(drug_target_table), "\n")
write.csv(drug_target_table, file.path(dir_out, paste0(ACC, "_drug_targets_ALL.csv")), row.names = FALSE)

## ---- keep only targets that are actually IN our combined TB PPI network ----

merged <- merge(drug_target_table, tb_network, by.x = "Target_Gene", by.y = "Gene")
merged <- merged[order(merged$Drug, -merged$Degree), ]
cat(">>> Of those,", nrow(merged), "pairs have a target gene present in the combined TB network",
    "(", length(unique(merged$Drug)), "/", length(drug_names), "compounds have at least one hit).\n")

final_table <- data.frame(
  Drug            = merged$Drug,
  Drug_Target     = merged$Target_Gene,
  TB_Module       = merged$Module,
  Target_Degree_in_TB_network = merged$Degree,
  Interaction_Types = merged$Interaction_Types,
  Sources         = merged$Sources
)
out_file <- file.path(dir_out, paste0(ACC, "_Drug_Target_TBnetwork.csv"))
write.csv(final_table, out_file, row.names = FALSE)

cat("\n>>> DONE. Saved to", dir_out, ":\n")
cat("  -", paste0(ACC, "_drug_targets_ALL.csv"), "        <- every DGIdb hit, unfiltered\n")
cat("  -", paste0(ACC, "_Drug_Target_TBnetwork.csv"), "   <- only targets present in the combined TB network (Drug -> Target -> Module)\n")
cat("\n>>> Compounds with NO target landing in the TB network at all:\n")
no_hit <- setdiff(drug_names, unique(merged$Drug))
if (length(no_hit)) print(no_hit) else cat("  (none -- every compound had at least one hit)\n")
cat("\n>>> Note: 'Biological pathway' labels are only available for the green/purple modules (Interferon_Response / Ribosome_Biogenesis-Chromatin).\n")
cat("    blue still has no enrichment name, so its targets show Module='blue' with no pathway label.\n")


