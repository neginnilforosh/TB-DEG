## Run ONCE from the repo root:   Rscript fix_script_paths.R scripts/wgcna_network
##
## Why: the scripts computed the data folder as  <script folder>/../<ACC>.  That only works when
## the scripts sit directly in scripts/.  After moving them to scripts/wgcna_network/ (as on
## GitHub) they look in scripts/GSE114192/ and stop with "cannot open file". This rewrites that one
## line in every script so it walks UP from the script's folder until it finds the real dataset folder (one that contains
## 02_metadata/, so a stray empty folder from an earlier failed run can't fool it) --
## so it works from scripts/, scripts/wgcna_network/, or anywhere deeper.

args <- commandArgs(trailingOnly = TRUE)
target_dir <- if (length(args) >= 1) args[1] else "scripts/wgcna_network"
if (!dir.exists(target_dir)) stop("Folder not found: ", target_dir, " (run this from the repo root)")

old_pat  <- '^(\\s*)BASE_DIR\\s*<-\\s*file\\.path\\(SCRIPT_DIR,\\s*"\\.\\.",\\s*ACC\\)\\s*$'
new_line <- 'BASE_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })  # walk up until the real dataset folder (has 02_metadata/) is found'

files <- list.files(target_dir, pattern = "\\.R$", full.names = TRUE)
n_patched <- 0
for (f in files) {
  x   <- readLines(f, warn = FALSE)
  hit <- grepl(old_pat, x, perl = TRUE)
  if (any(hit)) {
    x[hit] <- new_line
    writeLines(x, f)
    n_patched <- n_patched + 1
    cat("patched          ", basename(f), "\n")
  } else if (any(grepl("walk up until the", x, fixed = TRUE))) {
    cat("already patched  ", basename(f), "\n")
  } else {
    cat("no BASE_DIR line ", basename(f), "\n")
  }
}
cat("\nDone:", n_patched, "of", length(files), "scripts patched.\n")
