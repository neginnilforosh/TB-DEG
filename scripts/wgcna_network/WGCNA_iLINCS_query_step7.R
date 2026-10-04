## STEP 7: iLINCS compounds that reverse the TB signature (drugfindR getConcordants on UP and DOWN,
## then paired consensusConcordants). Only negative similarity (reversal) is kept.

required_pkgs <- c("drugfindR", "dplyr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) {
  stop("Missing package(s): ", paste(required_pkgs[!pkg_ok], collapse = ", "),
       "\n  Install with: install.packages(c(", paste0('"', required_pkgs[!pkg_ok], '"', collapse = ", "),
       '), repos = c("https://cogdisreslab.r-universe.dev", "https://cran.r-project.org"))')
}

## ---- CONFIG ----
ILINCS_LIBRARY       <- "CP"     # Chemical Perturbagen library (drug repurposing)
SIMILARITY_CUTOFF    <- 0.321    # drugfindR's own default for consensusConcordants()
REVERSAL_ONLY        <- TRUE
MIN_GENES_DIR        <- 5        # a direction with fewer genes is not queried (too small for iLINCS)     # keep only NEGATIVE similarity (true reversal), drop mimetic (positive) hits

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
dir_sig    <- file.path(OUT_DIR, "13_ilincs_signature")
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

## ---- query iLINCS separately for each direction ----
use_up <- nrow(TB_UP) >= MIN_GENES_DIR; use_down <- nrow(TB_DOWN) >= MIN_GENES_DIR
if (!use_up && !use_down) stop("TB_UP and TB_DOWN both have fewer than ", MIN_GENES_DIR, " genes: signature too small for iLINCS.")
if (!use_up)   cat("!! TB_UP has only", nrow(TB_UP), "genes: not queried, TB_DOWN only.\n")
if (!use_down) cat("!! TB_DOWN has only", nrow(TB_DOWN), "genes: not queried, TB_UP only.\n")
query <- function(sig, label) tryCatch(drugfindR::getConcordants(sig, ilincsLibrary = ILINCS_LIBRARY),
                                       error = function(e) { cat("!! getConcordants(", label, ") failed:", conditionMessage(e), "\n"); NULL })
concordants_up   <- if (use_up)   query(TB_UP, "TB_UP") else NULL
concordants_down <- if (use_down) query(TB_DOWN, "TB_DOWN") else NULL
if ((use_up && is.null(concordants_up)) || (use_down && is.null(concordants_down))) {
  stop("A getConcordants() call failed (see above): check the internet connection to ilincs.org. Nothing was written.")
}
if (!is.null(concordants_up))   write.csv(concordants_up,   file.path(dir_out, paste0(ACC, "_concordants_UP_raw.csv")),   row.names = FALSE)
if (!is.null(concordants_down)) write.csv(concordants_down, file.path(dir_out, paste0(ACC, "_concordants_DOWN_raw.csv")), row.names = FALSE)
cat(">>> Concordants: TB_UP", if (is.null(concordants_up)) "not queried" else nrow(concordants_up),
    "| TB_DOWN", if (is.null(concordants_down)) "not queried" else nrow(concordants_down), "\n")

## ---- consensus: paired if both directions returned results, otherwise the single available one ----
dirs <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(UP = concordants_up, DOWN = concordants_down))
if (length(dirs) == 0) stop("iLINCS returned no concordant signatures.")
mode <- if (length(dirs) == 2) "paired (TB_UP + TB_DOWN)" else paste0("single direction (", names(dirs), " only)")
consensus <- if (length(dirs) == 2) {
  drugfindR::consensusConcordants(dirs$UP, dirs$DOWN, paired = TRUE, cutoff = SIMILARITY_CUTOFF)
} else {
  drugfindR::consensusConcordants(dirs[[1]], paired = FALSE, cutoff = SIMILARITY_CUTOFF)
}
writeLines(mode, file.path(dir_out, paste0(ACC, "_consensus_mode.txt")))
cat(">>> Consensus mode:", mode, "|", nrow(consensus), "compounds with |similarity| >=", SIMILARITY_CUTOFF, "\n")

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
