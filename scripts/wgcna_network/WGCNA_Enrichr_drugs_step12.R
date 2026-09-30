## STEP 12: Enrichr (GO BP, KEGG, Reactome, WikiPathways, DSigDB) on each drug's full DGIdb target list.
## Drugs with < MIN_TARGETS targets are skipped and listed; DSigDB gives similar drug signatures,
## not pathways. Results are cached per drug, so an interrupted run can simply be restarted.

if (!requireNamespace("enrichR", quietly = TRUE)) stop("Missing package enrichR. Install with: install.packages(\"enrichR\")")
## enrichR sets its connection options only when attached with library()
library(enrichR)
if (!isTRUE(getOption("enrichR.live"))) stop("Enrichr website is not reachable from this machine (check internet/VPN); nothing was run.")

## ---- CONFIG ----
ACC          <- "GSE114192"
MIN_TARGETS  <- 3       # skip drugs with fewer known targets than this
PADJ_KEEP    <- 0.05
N_KEEP       <- 10      # top terms kept per drug per library in the long table
PAUSE_SEC    <- 0.5     # politeness delay between drugs
## newest Enrichr library matching each pattern (library names carry a year)
LIB_PATTERNS <- c(GO_BP        = "^GO_Biological_Process_20[0-9]{2}$",
                  KEGG         = "^KEGG_20[0-9]{2}_Human$",
                  Reactome     = "^Reactome_(Pathways_)?20[0-9]{2}$",
                  WikiPathways = "^WikiPathways?_20[0-9]{2}_Human$",
                  DSigDB       = "^DSigDB$")

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
BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found
dir_dt  <- file.path(BASE_DIR, "14_drug_targets")
dir_fin <- file.path(BASE_DIR, "16_final_tables")
dir_out <- file.path(BASE_DIR, "17_pathway_annotation")
dir_cache <- file.path(dir_out, "cache")
dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)

## ---- helpers ----
# newest library whose name matches the pattern (year = the 4 digits in the name)
pick_library <- function(available, pattern) {
  hit <- available[grepl(pattern, available)]
  if (length(hit) == 0) return(NA_character_)
  yr <- suppressWarnings(as.integer(sub(".*_(20[0-9]{2}).*", "\\1", hit)))
  yr[is.na(yr)] <- 0L
  hit[order(yr, hit, decreasing = TRUE)][1]
}
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
# one enrichr() result (list of data.frames) -> long table of significant terms
tidy_enrichr <- function(res, drug, lib_of, padj_keep, n_keep) {
  out <- list()
  for (lib in names(res)) {
    df <- res[[lib]]
    if (!is.data.frame(df) || nrow(df) == 0) next
    df <- df[!is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_keep, , drop = FALSE]
    if (nrow(df) == 0) next
    df <- utils::head(df[order(df$Adjusted.P.value), ], n_keep)
    out[[lib]] <- data.frame(Drug = drug, Family = lib_of[[lib]], Library = lib, Term = df$Term, Overlap = df$Overlap,
                             Adjusted.P.value = df$Adjusted.P.value, Odds.Ratio = df$Odds.Ratio,
                             Combined.Score = df$Combined.Score, Genes = df$Genes, stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else NULL
}
# long table -> one row per drug, best term per family as "Term [adj p]"
summarise_drugs <- function(long, drugs_run, n_targets, families) {
  s <- data.frame(Drug = drugs_run, Enrichr_N_targets_used = unname(n_targets[drugs_run]), stringsAsFactors = FALSE)
  for (fam in families) {
    col <- paste0("Enrichr_", fam)
    s[[col]] <- NA_character_
    if (is.null(long)) next
    sub <- long[long$Family == fam, ]
    if (nrow(sub) == 0) next
    best <- sub[order(sub$Drug, sub$Adjusted.P.value), ]
    best <- best[!duplicated(best$Drug), ]
    lab <- paste0(best$Term, " [", signif(best$Adjusted.P.value, 2), "]")
    s[[col]][match(best$Drug, s$Drug)] <- lab
  }
  s
}

## ---- inputs ----
tg  <- read.csv(file.path(dir_dt, paste0(ACC, "_drug_targets_ALL.csv")), stringsAsFactors = FALSE)
fin <- read.csv(file.path(dir_fin, paste0(ACC, "_FINAL_DrugTable.csv")), stringsAsFactors = FALSE)
drugs <- unique(fin$Drug)
targets <- lapply(setNames(drugs, drugs), function(d) { g <- unique(tg$Target_Gene[tg$Drug == d]); g[!is.na(g) & nzchar(g)] })
n_targets <- vapply(targets, length, integer(1))
run_drugs  <- names(n_targets)[n_targets >= MIN_TARGETS]
skipped    <- data.frame(Drug = names(n_targets)[n_targets < MIN_TARGETS], N_known_targets = n_targets[n_targets < MIN_TARGETS], stringsAsFactors = FALSE)
cat(">>>", length(drugs), "drugs in the final table;", length(run_drugs), "have >=", MIN_TARGETS, "known targets and will be run;",
    nrow(skipped), "skipped (too few targets).\n")
write.csv(skipped, file.path(dir_out, paste0(ACC, "_DrugPathways_skipped.csv")), row.names = FALSE)

## ---- which Enrichr libraries ----
dbs_df <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) NULL)
if (is.null(dbs_df) || !("libraryName" %in% names(dbs_df))) stop("Could not get the list of Enrichr libraries -- the website may be down; try again later.")
available <- dbs_df$libraryName
libs <- vapply(LIB_PATTERNS, function(p) pick_library(available, p), character(1))
cat(">>> Libraries:\n"); for (f in names(libs)) cat(sprintf("    %-13s %s\n", f, ifelse(is.na(libs[[f]]), "!! NOT FOUND (skipped) -- check names with enrichR::listEnrichrDbs()", libs[[f]])))
libs <- libs[!is.na(libs)]
if (length(libs) == 0) stop("None of the expected Enrichr libraries were found -- the library naming may have changed.")
lib_of <- setNames(names(libs), libs)                 # library name -> family

## ---- per-drug Enrichr (cached, resumable) ----
cat(">>> Running Enrichr for", length(run_drugs), "drugs (cached to", dir_cache, ")...\n")
long_list <- list(); n_new <- 0; n_fail <- 0
for (i in seq_along(run_drugs)) {
  d <- run_drugs[i]; cf <- file.path(dir_cache, paste0(safe_name(d), ".rds"))
  if (file.exists(cf)) { long_list[[d]] <- readRDS(cf); next }
  res <- NULL
  for (attempt in 1:2) {
    res <- tryCatch(suppressMessages(enrichR::enrichr(targets[[d]], unname(libs))), error = function(e) NULL)
    if (!is.null(res)) break
    Sys.sleep(5)
  }
  if (is.null(res)) { n_fail <- n_fail + 1; cat("  !! Enrichr failed for", d, "(will be retried on the next run)\n"); next }
  tidy <- tidy_enrichr(res, d, lib_of, PADJ_KEEP, N_KEEP)
  saveRDS(if (is.null(tidy)) data.frame() else tidy, cf)          # cache even "nothing significant", so it isn't re-queried
  long_list[[d]] <- if (is.null(tidy)) data.frame() else tidy
  n_new <- n_new + 1
  if (n_new %% 25 == 0) cat("   ...", i, "/", length(run_drugs), "drugs processed\n")
  Sys.sleep(PAUSE_SEC)
}
long <- do.call(rbind, long_list[vapply(long_list, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))])
cat(">>> New Enrichr calls this run:", n_new, "| failed:", n_fail, "| drugs with >=1 significant term:",
    if (is.null(long)) 0 else length(unique(long$Drug)), "\n")

done_drugs <- names(long_list)
summ <- summarise_drugs(long, done_drugs, n_targets, names(libs))
write.csv(summ, file.path(dir_out, paste0(ACC, "_DrugPathways_Enrichr_summary.csv")), row.names = FALSE)
if (!is.null(long)) write.csv(long, file.path(dir_out, paste0(ACC, "_DrugPathways_Enrichr_long.csv")), row.names = FALSE)

cat("\n>>> DONE. In", dir_out, ":\n")
cat("  -", paste0(ACC, "_DrugPathways_Enrichr_summary.csv"), "(one row per drug; step 10 merges this into the final drug table)\n")
cat("  -", paste0(ACC, "_DrugPathways_Enrichr_long.csv"), "(every significant term, with the overlapping target genes)\n")
cat("  -", paste0(ACC, "_DrugPathways_skipped.csv"), "(drugs with too few known targets)\n")
cat(">>> Next: re-run WGCNA_final_tables_step10.R so the drug table picks the summary up.\n")
