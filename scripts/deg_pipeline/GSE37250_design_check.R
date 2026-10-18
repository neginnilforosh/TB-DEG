## GSE37250 design check, before the TB-HIV analysis.
##
## Questions to answer from the data GEO provides:
##   1. how many samples per disease group, HIV status and country
##   2. whether any technical batch information exists (BeadChip barcode, array position, scan date, plate, run),
##      and whether it is balanced across the groups that will be compared
##   3. which clinical variables are available (age, sex, CD4 count, antiretroviral treatment, ...)
## Sources: the series matrix (sample characteristics) and the header of the non-normalised data file, whose
## column names often carry the BeadChip barcode. Samples are always matched by name, never by position.
## Run from anywhere in the repository:  Rscript scripts/deg_pipeline/GSE37250_design_check.R
## Outputs in GSE37250/00_design/.

ACC <- "GSE37250"
for (p in c("GEOquery", "Biobase")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, " (BiocManager::install('", p, "'))")
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT, ACC, "00_design"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RAW <- file.path(ROOT, "external_data", ACC); dir.create(RAW, recursive = TRUE, showWarnings = FALSE)

## ---- 1. sample characteristics ----
options(timeout = 3600)   # the series matrix is about 100 MB; the default 2-minute limit is too short
gz_ok <- function(f) file.exists(f) && identical(suppressWarnings(system2("gzip", c("-t", shQuote(f)), stdout = FALSE, stderr = FALSE)), 0L)   # gzip -t catches truncated downloads
sm <- file.path(RAW, paste0(ACC, "_series_matrix.txt.gz"))
if (file.exists(sm) && !gz_ok(sm)) { cat(">>> removing an incomplete download:", basename(sm), "\n"); file.remove(sm) }
gse <- if (file.exists(sm)) list(GEOquery::getGEO(filename = sm, GSEMatrix = TRUE)) else GEOquery::getGEO(ACC, GSEMatrix = TRUE, destdir = RAW)
pd <- do.call(rbind, lapply(gse, Biobase::pData))
cat(">>>", nrow(pd), "samples on", length(unique(pd$platform_id)), "platform(s):", paste(unique(pd$platform_id), collapse = ", "), "\n")
## "key: value" characteristics -> one column per key
chars <- grep("^characteristics_ch1", names(pd), value = TRUE)
kv <- do.call(rbind, lapply(seq_len(nrow(pd)), function(i) {
  v <- unlist(pd[i, chars]); v <- v[nzchar(v) & !is.na(v)]
  data.frame(gsm = pd$geo_accession[i], key = trimws(sub(":.*$", "", v)), value = trimws(sub("^[^:]*:", "", v)))
}))
keys <- unique(kv$key)
S <- data.frame(gsm = pd$geo_accession, title = pd$title, stringsAsFactors = FALSE)
for (k in keys) S[[make.names(k)]] <- kv$value[kv$key == k][match(S$gsm, kv$gsm[kv$key == k])]
for (f in intersect(c("description", "source_name_ch1", "scan_protocol", "hyb_protocol", "submission_date", "last_update_date"), names(pd)))
  S[[f]] <- pd[[f]]
cat(">>> characteristics found:", paste(keys, collapse = " | "), "\n")

find_col <- function(pattern) { h <- grep(pattern, names(S), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
c_dis <- find_col("disease"); c_hiv <- find_col("hiv"); c_geo <- find_col("geograph|region|country|site")
cat(">>> using: disease =", c_dis, "| HIV =", c_hiv, "| country =", c_geo, "\n")

## ---- 2. technical batch: characteristics, descriptions and the non-normalised file header ----
batch_cols <- grep("chip|sentrix|barcode|batch|array|slide|scan|plate|run|date|hyb", names(S), ignore.case = TRUE, value = TRUE)
batch_cols <- setdiff(batch_cols, c("submission_date", "last_update_date"))   # GEO bookkeeping, not lab batches
cat(">>> batch-like fields in the sample characteristics:", if (length(batch_cols)) paste(batch_cols, collapse = ", ") else "none", "\n")
hdr <- NULL
f_raw <- list.files(RAW, pattern = "non.?normali[sz]ed|raw", ignore.case = TRUE, full.names = TRUE)
f_raw <- f_raw[!grepl("series_matrix", f_raw) & vapply(f_raw, function(f) !grepl("\\.gz$", f) || gz_ok(f), TRUE)]
if (!length(f_raw)) {
  sup <- tryCatch(GEOquery::getGEOSuppFiles(ACC, baseDir = RAW, makeDirectory = FALSE, filter_regex = "non.?normali[sz]ed|raw", fetch_files = TRUE),
                  error = function(e) { cat("  supplementary download failed:", conditionMessage(e), "\n"); NULL })
  f_raw <- list.files(RAW, pattern = "non.?normali[sz]ed|raw", ignore.case = TRUE, full.names = TRUE)
  f_raw <- f_raw[!grepl("series_matrix", f_raw)]
}
if (length(f_raw)) {
  con <- if (grepl("\\.gz$", f_raw[1])) gzfile(f_raw[1]) else file(f_raw[1]); first <- readLines(con, n = 10); close(con)
  hl <- grep("\t", first)[1]; hdr <- strsplit(first[hl], "\t")[[1]]
  cat(">>> non-normalised file", basename(f_raw[1]), "has", length(hdr), "columns; examples:", paste(utils::head(hdr, 6), collapse = " | "), "\n")
  sig <- grep("pval|detection", hdr, ignore.case = TRUE, invert = TRUE, value = TRUE)[-1]
  ## chip barcode pattern: 10-digit Sentrix barcode followed by _A.._L or _R0xC0y
  bc <- regmatches(sig, regexpr("[0-9]{9,10}_[A-L]|[0-9]{9,12}_R[0-9]{2}C[0-9]{2}", sig))
  if (length(bc)) {
    cat("    Sentrix barcodes found in", length(bc), "column names\n")
    ## link a column to a sample only if its title or GSM accession appears in the column name (never by position)
    link <- function(i) { j <- which(vapply(sig, function(x) grepl(S$title[i], x, fixed = TRUE) || grepl(S$gsm[i], x, fixed = TRUE), TRUE))
      if (length(j) == 1) sig[j] else NA_character_ }
    S$raw_column <- vapply(seq_len(nrow(S)), link, "")
    S$chip <- vapply(S$raw_column, function(x) { if (is.na(x)) return(NA_character_)
      m <- regmatches(x, regexpr("[0-9]{9,12}(?=_)", x, perl = TRUE)); if (length(m)) m else NA_character_ }, "")
    cat("    linked to samples by title or accession:", sum(!is.na(S$raw_column)), "of", nrow(S),
        if (sum(!is.na(S$raw_column)) == 0) "-- column names cannot be tied to samples, so the chip is not used" else "", "\n")
  } else cat("    no BeadChip barcodes in the column names\n")
}
write.csv(S, file.path(OUT, paste0(ACC, "_sample_table.csv")), row.names = FALSE)

## ---- 3. tables ----
tab <- function(...) { t <- table(..., useNA = "ifany"); print(t); invisible(t) }
if (!is.na(c_dis) && !is.na(c_hiv)) {
  S$group <- paste(S[[c_dis]], "| HIV", S[[c_hiv]])
  cat("\n========== group x country ==========\n"); t1 <- if (!is.na(c_geo)) tab(S$group, S[[c_geo]]) else tab(S$group)
  write.csv(as.data.frame.matrix(t1), file.path(OUT, paste0(ACC, "_group_by_country.csv")))
  for (b in c(batch_cols, if ("chip" %in% names(S)) "chip")) {
    if (length(unique(stats::na.omit(S[[b]]))) < 2) next
    cat("\n========== group x", b, "==========\n"); t2 <- tab(S$group, S[[b]])
    write.csv(as.data.frame.matrix(t2), file.path(OUT, paste0(ACC, "_group_by_", make.names(b), ".csv")))
    hivpos <- grepl("pos|\\+|yes", S[[c_hiv]], ignore.case = TRUE)
    if (any(hivpos)) { x <- stats::chisq.test(table(S$group[hivpos], S[[b]][hivpos]), simulate.p.value = TRUE, B = 2000)
      cat(sprintf("  association of %s with group among HIV-positive samples: chi-square p = %.3g\n", b, x$p.value)) }
  }
}
cat("\n========== clinical variables available ==========\n")
clin <- grep("age|sex|gender|cd4|art|antiretro|treat|bmi|smok|diabet", names(S), ignore.case = TRUE, value = TRUE)
if (!length(clin)) cat("  none beyond disease, HIV status and country\n") else for (k in clin)
  cat(sprintf("  %-25s %d of %d filled | examples: %s\n", k, sum(!is.na(S[[k]]) & nzchar(S[[k]])), nrow(S), paste(utils::head(unique(S[[k]]), 5), collapse = ", ")))
cat("\n>>> DONE. Outputs in", OUT, "(raw downloads in external_data/, which is not tracked)\n")
