## TB-HIV validation datasets: what is there before deciding how to analyse them.
##
## For GSE107104 and GSE162164 (whole-blood RNA-seq, people with HIV with and without active TB):
##   - sample characteristics from the series matrix: TB status, ART, CD4 and anything else recorded
##   - the supplementary files on GEO; expression-like files are downloaded and inspected
##   - for each expression table: gene identifiers, number of columns, whether column names match GEO sample
##     titles or accessions (samples are only ever linked by name), and whether the values look like raw counts
##     (non-negative integers) or normalised values, which decides between DESeq2/voom and limma on log values
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/tbhiv_validation_check.R
## Outputs in <accession>/00_design/; downloads in external_data/ (not tracked).

ACCS <- c("GSE107104", "GSE162164")
for (p in c("GEOquery", "Biobase")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
options(timeout = 3600)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })

read_any <- function(f) {
  if (grepl("\\.xlsx?$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) { cat("    (install readxl to read", basename(f), ")\n"); return(NULL) }
    return(as.data.frame(readxl::read_excel(f)))
  }
  sep <- if (grepl("\\.csv(\\.gz)?$", f, ignore.case = TRUE)) "," else "\t"
  tryCatch(utils::read.table(f, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "", fill = TRUE),
           error = function(e) { cat("    could not read:", conditionMessage(e), "\n"); NULL })
}

for (acc in ACCS) {
  cat("\n================", acc, "================\n")
  RAW <- file.path(ROOT, "external_data", acc); OUT <- file.path(ROOT, acc, "00_design")
  dir.create(RAW, recursive = TRUE, showWarnings = FALSE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  gse <- tryCatch(GEOquery::getGEO(acc, GSEMatrix = TRUE, destdir = RAW), error = function(e) { cat("  series matrix failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(gse)) next
  pd <- do.call(rbind, lapply(gse, Biobase::pData))
  chars <- grep("^characteristics_ch1", names(pd), value = TRUE)
  S <- data.frame(gsm = pd$geo_accession, title = pd$title, source = pd$source_name_ch1, stringsAsFactors = FALSE)
  for (cc in chars) for (i in seq_len(nrow(pd))) { v <- pd[[cc]][i]; if (is.na(v) || !nzchar(v)) next
    k <- make.names(trimws(sub(":.*$", "", v))); S[i, k] <- trimws(sub("^[^:]*:", "", v)) }
  write.csv(S, file.path(OUT, paste0(acc, "_sample_table.csv")), row.names = FALSE)
  cat(">>>", nrow(S), "samples | platform", paste(unique(pd$platform_id), collapse = ", "), "\n")
  cat(">>> characteristics:", paste(setdiff(names(S), c("gsm", "title", "source")), collapse = ", "), "\n")
  for (k in setdiff(names(S), c("gsm", "title"))) {
    v <- S[[k]]; if (length(unique(v)) <= 12) { cat("  ", k, ":\n"); print(table(v, useNA = "ifany")) }
  }
  cat("  example titles:", paste(utils::head(S$title, 6), collapse = " | "), "\n")

  sup <- tryCatch(GEOquery::getGEOSuppFiles(acc, baseDir = RAW, makeDirectory = FALSE, fetch_files = FALSE), error = function(e) NULL)
  if (is.null(sup) || !nrow(sup)) { cat(">>> no supplementary files listed\n"); next }
  cat(">>> supplementary files:\n"); print(sup$fname)
  want <- sup$fname[grepl("count|express|fpkm|tpm|norm|matrix|txt|csv|tsv|xls", sup$fname, ignore.case = TRUE) & !grepl("\\.tar$|RAW", sup$fname)]
  for (fn in want) {
    dest <- file.path(RAW, fn)
    if (!file.exists(dest)) tryCatch(utils::download.file(sup$url[sup$fname == fn], dest, mode = "wb", quiet = TRUE),
                                     error = function(e) cat("  download failed:", fn, "\n"))
    if (!file.exists(dest)) next
    cat("\n  --", fn, "(", round(file.size(dest) / 1e6, 1), "MB )\n")
    d <- read_any(dest); if (is.null(d)) next
    cat("    ", nrow(d), "rows x", ncol(d), "columns | first columns:", paste(utils::head(names(d), 5), collapse = " | "), "\n")
    cat("     first identifiers:", paste(utils::head(as.character(d[[1]]), 4), collapse = ", "), "\n")
    num <- vapply(d, is.numeric, TRUE); m <- as.matrix(d[, num, drop = FALSE])
    if (ncol(m)) {
      ints <- all(abs(m - round(m)) < 1e-8, na.rm = TRUE); neg <- any(m < 0, na.rm = TRUE)
      cat(sprintf("     numeric columns: %d | all integers: %s | negative values: %s | max %.1f | median %.2f\n",
                  ncol(m), ints, neg, max(m, na.rm = TRUE), stats::median(m, na.rm = TRUE)))
      cat("     looks like:", if (ints && !neg) "RAW COUNTS" else if (!neg && max(m, na.rm = TRUE) > 100) "normalised counts (not integers)" else "log-scale or other values", "\n")
      hits_t <- sum(vapply(colnames(m), function(cn) any(vapply(S$title, function(t) grepl(t, cn, fixed = TRUE) || grepl(cn, t, fixed = TRUE), TRUE)), TRUE))
      hits_g <- sum(colnames(m) %in% S$gsm)
      cat("     columns matching GEO titles:", hits_t, "| matching GSM accessions:", hits_g, "of", ncol(m), "\n")
      cat("     column names:", paste(utils::head(colnames(m), 8), collapse = " | "), "\n")
    }
  }
}
cat("\n>>> DONE\n")
