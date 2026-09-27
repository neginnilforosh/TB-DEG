## =============================================================

## Follows drugfindR's own documented paired workflow exactly:
##   getConcordants(TB_UP,   ilincsLibrary="CP")
##   getConcordants(TB_DOWN, ilincsLibrary="CP")
##   consensusConcordants(concordants_up, concordants_down, paired=TRUE)
## =============================================================

required_pkgs <- c("drugfindR", "dplyr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "),
       '), repos = c("https://cogdisreslab.r-universe.dev", "https://cran.r-project.org"))')
}

## ---- CONFIG ----
ACC                  <- "GSE114192"
ILINCS_LIBRARY       <- "CP"     # Chemical Perturbagen library (drug repurposing)
SIMILARITY_CUTOFF    <- 0.321    # drugfindR's own default for consensusConcordants()
REVERSAL_ONLY        <- TRUE     # keep only NEGATIVE similarity (true reversal), drop mimetic (positive) hits

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
dir_sig    <- file.path(BASE_DIR, "13_ilincs_signature")
dir_out    <- dir_sig

## ---- load the TB_UP / TB_DOWN signature from step 6 ----
up_file   <- file.path(dir_sig, paste0(ACC, "_TB_UP_signature.csv"))
down_file <- file.path(dir_sig, paste0(ACC, "_TB_DOWN_signature.csv"))
if (!file.exists(up_file) || !file.exists(down_file)) {
  stop("Missing TB_UP/TB_DOWN files -- run WGCNA_iLINCS_signature_step6.R first.")
}
TB_UP   <- read.csv(up_file,   stringsAsFactors = FALSE)
TB_DOWN <- read.csv(down_file, stringsAsFactors = FALSE)
cat(">>> TB_UP:", nrow(TB_UP), "genes.  TB_DOWN:", nrow(TB_DOWN), "genes.\n")
cat(">>> Querying iLINCS (", ILINCS_LIBRARY, "library ) -- this makes live network calls",
    "and can take a while for each direction...\n")

## ---- query iLINCS separately for each direction (drugfindR's documented pattern) ----
concordants_up   <- tryCatch(drugfindR::getConcordants(TB_UP,   ilincsLibrary = ILINCS_LIBRARY),
                              error = function(e) { cat("!! getConcordants(TB_UP) failed:", conditionMessage(e), "\n"); NULL })
concordants_down <- tryCatch(drugfindR::getConcordants(TB_DOWN, ilincsLibrary = ILINCS_LIBRARY),
                              error = function(e) { cat("!! getConcordants(TB_DOWN) failed:", conditionMessage(e), "\n"); NULL })

if (is.null(concordants_up) || is.null(concordants_down)) {
  stop("At least one getConcordants() call failed -- see message above. Common causes: ",
       "no internet from this machine to ilincs.org, or ilincs.org temporarily down. ",
       "Nothing else in this script ran, so nothing was overwritten.")
}
cat(">>> Got", nrow(concordants_up), "concordants for TB_UP,", nrow(concordants_down), "for TB_DOWN.\n")

write.csv(concordants_up,   file.path(dir_out, paste0(ACC, "_concordants_UP_raw.csv")),   row.names = FALSE)
write.csv(concordants_down, file.path(dir_out, paste0(ACC, "_concordants_DOWN_raw.csv")), row.names = FALSE)

## ---- paired consensus (drugfindR's documented combination step) ----
consensus <- drugfindR::consensusConcordants(concordants_up, concordants_down,
                                              paired = TRUE, cutoff = SIMILARITY_CUTOFF)
cat(">>> Consensus (|similarity| >=", SIMILARITY_CUTOFF, "):", nrow(consensus), "compounds.\n")

sim_col <- if ("similarity" %in% names(consensus)) "similarity" else
           if ("Similarity" %in% names(consensus)) "Similarity" else NA
if (is.na(sim_col)) {
  cat("!! WARNING: couldn't find a similarity column by the expected name -- columns are:",
      paste(names(consensus), collapse = ", "), "-- check REVERSAL_ONLY filtering manually.\n")
} else if (REVERSAL_ONLY) {
  n_before <- nrow(consensus)
  consensus <- consensus[consensus[[sim_col]] < 0, ]
  consensus <- consensus[order(consensus[[sim_col]]), ]  # most negative (strongest reversal) first
  cat(">>> Kept", nrow(consensus), "/", n_before, "reversal (negative-similarity) compounds.\n")
}

out_file <- file.path(dir_out, paste0(ACC, "_iLINCS_candidate_compounds.csv"))
write.csv(consensus, out_file, row.names = FALSE)

cat("\n>>> DONE. Raw per-direction results + final candidate table saved to", dir_out, ":\n")
cat("  -", paste0(ACC, "_concordants_UP_raw.csv"), "/", paste0(ACC, "_concordants_DOWN_raw.csv"), "\n")
cat("  -", paste0(ACC, "_iLINCS_candidate_compounds.csv"), " <- ranked reversal candidates, strongest first\n")
