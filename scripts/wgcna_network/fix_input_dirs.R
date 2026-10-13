## Make the network scripts read their INPUTS from the run folder when it has them.
##
## Until now every run took its DEG-pipeline inputs (06_deg_results, 09_wgcna_input) from the dataset
## folder and only wrote OUTPUTS to the run folder. That was right for the pilot and for GSE161829,
## whose inputs live in the dataset folder, but the TB-DM run's inputs were written by
## GSE114192_comorbidity_run.R into GSE114192/TB_DM/, so step 1 picked up the pilot's 82-sample
## matrix and trait file instead.
##
## This adds in_dir() to config.R -- the run folder's copy if it exists, otherwise the dataset
## folder's -- and uses it in the four scripts that read DEG-pipeline outputs (steps 1, 3, 6, 10).
## The pilot (RUN = "") and GSE161829 (no inputs inside its run folder) resolve exactly as before.
##
## Run once from scripts/wgcna_network:   Rscript fix_input_dirs.R     (safe to run twice)

here <- getwd()
if (!file.exists(file.path(here, "config.R"))) stop("Run this from the folder that holds config.R")

## ---- 1. config.R: add in_dir() right after run_out_dir() ----
cfg <- readLines("config.R", warn = FALSE)
if (any(grepl("^in_dir <- function", cfg))) {
  cat("config.R: in_dir() already present\n")
} else {
  at <- grep("^run_out_dir <- function", cfg)
  if (length(at) != 1) stop("Could not find run_out_dir() in config.R")
  ## run_out_dir may span one or several lines; insert after the line that closes it
  end <- at
  if (!grepl("\\}\\s*$", cfg[at]) && grepl("\\{", cfg[at])) {
    while (end < length(cfg) && !grepl("^\\}", cfg[end])) end <- end + 1
  }
  add <- c("## inputs: the run folder's copy when it has one (e.g. GSE114192/TB_DM/09_wgcna_input),",
           "## otherwise the dataset folder's -- so existing runs resolve exactly as before",
           "in_dir <- function(base, sub) {",
           "  r <- file.path(run_out_dir(base), sub)",
           "  if (dir.exists(r)) r else file.path(base, sub)",
           "}")
  cfg <- append(cfg, add, after = end)
  writeLines(cfg, "config.R")
  cat("config.R: in_dir() added\n")
}

## ---- 2. the four scripts ----
targets <- c(WGCNA_TBmodules_step1.R        = "09_wgcna_input",
             WGCNA_DEG_integration_step3.R  = "06_deg_results",
             WGCNA_iLINCS_signature_step6.R = "06_deg_results",
             WGCNA_final_tables_step10.R    = "06_deg_results")
for (f in names(targets)) {
  if (!file.exists(f)) { cat("  skipped, not found:", f, "\n"); next }
  sub <- targets[[f]]
  x <- readLines(f, warn = FALSE)
  old <- paste0('file.path(BASE_DIR, "', sub, '")')
  new <- paste0('in_dir(BASE_DIR, "', sub, '")')
  hits <- grepl(old, x, fixed = TRUE)
  if (!any(hits)) {
    cat(sprintf("  %-32s already patched\n", f)); next
  }
  x[hits] <- gsub(old, new, x[hits], fixed = TRUE)
  writeLines(x, f)
  cat(sprintf("  %-32s patched (%d line)\n", f, sum(hits)))
}

## ---- 3. check ----
for (f in c("config.R", names(targets))) if (file.exists(f)) invisible(parse(f))
cat("\nAll files parse. Run step 1 again.\n")
