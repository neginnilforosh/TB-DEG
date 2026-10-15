## Metabolomics metadata check (ST000974 plasma, ST000975 serum, ST000976 plasma diluted with RPMI).
##
## Before choosing the main metabolomics analysis we need to know, for each study:
##   - whether participant identifiers exist, so repeated samples from one person (BL, M06, M18) can be linked
##   - whether time to TB diagnosis is recorded (e.g. < 5 months vs >= 5 months before diagnosis)
##   - how many samples there are per group (progressor / non-progressor), time point and site
##   - whether the three studies share samples or participants
## The script downloads each study's mwTab file from the Metabolomics Workbench REST service and reads the
## SUBJECT_SAMPLE_FACTORS section (subject, sample, factors, additional sample data).
## Run from anywhere in the repository:  Rscript scripts/wgcna_network/metabolomics_metadata_check.R
## Outputs in GSE114192/21_metabolomics/metadata_check/.

STUDIES <- c(ST000974 = "Plasma", ST000975 = "Serum", ST000976 = "Plasma diluted with RPMI")
API <- "https://www.metabolomicsworkbench.org/rest/study"
for (p in c("httr", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
OUT <- file.path(ROOT, "GSE114192", "21_metabolomics", "metadata_check"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

get_text <- function(url) {
  r <- httr::GET(url, httr::timeout(120))
  if (httr::status_code(r) != 200) stop("Request failed (", httr::status_code(r), "): ", url)
  httr::content(r, "text", encoding = "UTF-8")
}
analysis_ids <- function(st) {
  j <- jsonlite::fromJSON(get_text(paste0(API, "/study_id/", st, "/analysis")), simplifyVector = FALSE)
  if (!is.null(j$analysis_id)) j <- list(j)
  unique(vapply(j, function(x) x$analysis_id, ""))
}
## SUBJECT_SAMPLE_FACTORS lines: tag \t subject \t sample \t factors \t additional data
parse_ssf <- function(txt) {
  l <- grep("^SUBJECT_SAMPLE_FACTORS", strsplit(txt, "\r?\n")[[1]], value = TRUE)
  l <- l[!grepl("SUBJECT\\(optional\\)", l)]
  p <- strsplit(l, "\t")
  d <- data.frame(subject = trimws(vapply(p, function(x) if (length(x) >= 2) x[2] else "", "")),
                  sample  = trimws(vapply(p, function(x) if (length(x) >= 3) x[3] else "", "")),
                  factors = trimws(vapply(p, function(x) if (length(x) >= 4) x[4] else "", "")),
                  additional = trimws(vapply(p, function(x) if (length(x) >= 5) paste(x[5:length(x)], collapse = " ") else "", "")),
                  stringsAsFactors = FALSE)
  f <- lapply(strsplit(d$factors, "\\s*\\|\\s*"), function(x) { k <- sub(":.*$", "", x); v <- sub("^[^:]*:", "", x); setNames(trimws(v), trimws(k)) })
  keys <- unique(unlist(lapply(f, names)))
  for (k in keys) d[[k]] <- vapply(f, function(x) if (k %in% names(x)) x[[k]] else NA_character_, "")
  ## additional sample data: "KEY=value; KEY=value" -> one column per key
  a <- lapply(strsplit(d$additional, "\\s*;\\s*"), function(x) { x <- x[grepl("=", x)]
    setNames(trimws(sub("^[^=]*=", "", x)), trimws(sub("=.*$", "", x))) })
  for (k in setdiff(unique(unlist(lapply(a, names))), names(d))) d[[k]] <- vapply(a, function(x) if (k %in% names(x)) x[[k]] else NA_character_, "")
  d
}

all <- list(); summ <- list()
for (st in names(STUDIES)) {
  cat("\n==========", st, "(", STUDIES[[st]], ") ==========\n")
  ids <- tryCatch(analysis_ids(st), error = function(e) { cat("  could not list analyses:", conditionMessage(e), "\n"); character(0) })
  if (!length(ids)) next
  txt <- get_text(paste0(API, "/analysis_id/", ids[1], "/mwtab/txt"))
  writeLines(txt, file.path(OUT, paste0(st, "_", ids[1], "_mwtab.txt")))
  d <- parse_ssf(txt); d$study <- st; d$sample_type <- STUDIES[[st]]
  cat("  analysis", ids[1], "|", nrow(d), "samples | factor fields:", paste(setdiff(names(d), c("subject", "sample", "factors", "additional", "study", "sample_type")), collapse = ", "), "\n")
  if (!any(nzchar(d$subject) & d$subject != "-") && "DONOR_ID" %in% names(d)) d$subject <- d$DONOR_ID   # IDs sit in the additional data
  has_subj <- any(!is.na(d$subject) & nzchar(d$subject) & d$subject != "-")
  cat("  participant identifiers:", if (has_subj) paste(length(unique(d$subject)), "distinct for", nrow(d), "samples") else "NONE", "\n")
  add <- unique(unlist(lapply(strsplit(d$additional[nzchar(d$additional)], ";\\s*"), function(x) sub("=.*$", "", x))))
  cat("  additional sample data fields:", if (length(add)) paste(add, collapse = ", ") else "none", "\n")
  ttd <- grep("time|month|diagnos|progress|days", c(names(d), add), ignore.case = TRUE, value = TRUE)
  cat("  fields that could hold time to diagnosis:", if (length(ttd)) paste(unique(ttd), collapse = ", ") else "none found", "\n")
  if (all(c("GROUP", "TIMEPOINT") %in% names(d))) { cat("  samples by group and time point:\n"); print(table(d$GROUP, d$TIMEPOINT)) }
  if (all(c("GROUP", "REGION") %in% names(d))) { cat("  samples by group and site:\n"); print(table(d$GROUP, d$REGION)) }
  if (has_subj && "TIMEPOINT" %in% names(d)) {
    per <- table(d$subject)
    cat("  participants:", length(per), "| with more than one sample:", sum(per > 1), "| samples per participant:",
        paste(names(table(per)), table(per), sep = " x", collapse = ", "), "\n")
    tp <- tapply(d$TIMEPOINT, d$subject, function(x) paste(sort(unique(x)), collapse = "+"))
    cat("  time-point patterns per participant:\n"); print(sort(table(tp), decreasing = TRUE))
    if ("GROUP" %in% names(d)) { cat("  participants by group:\n"); print(table(tapply(d$GROUP, d$subject, `[`, 1))) }
  }
  if ("CASE_GROUP" %in% names(d)) { cat("  CASE_GROUP by group and time point:\n"); print(table(d$CASE_GROUP, paste(d$GROUP, d$TIMEPOINT), useNA = "ifany")) }
  if ("TIME_TO_TB" %in% names(d)) {
    tt <- suppressWarnings(as.numeric(d$TIME_TO_TB))
    cat("  TIME_TO_TB values (cases): numeric", sum(!is.na(tt) & d$GROUP == "case"), "of", sum(d$GROUP == "case"),
        "| examples:", paste(utils::head(unique(d$TIME_TO_TB[d$GROUP == "case"]), 8), collapse = ", "), "\n")
    if (any(!is.na(tt))) { cat("  TIME_TO_TB summary by time point (cases):\n"); print(tapply(tt[d$GROUP == "case"], d$TIMEPOINT[d$GROUP == "case"], summary)) }
    cat("  TIME_TO_TB examples (controls):", paste(utils::head(unique(d$TIME_TO_TB[d$GROUP == "control"]), 5), collapse = ", "), "\n")
  }
  summ[[st]] <- data.frame(Study = st, Sample_type = STUDIES[[st]], Analysis = ids[1], Samples = nrow(d),
    Participant_IDs = if (has_subj) length(unique(d$subject)) else 0, Time_to_diagnosis_field = if (length(ttd)) paste(unique(ttd), collapse = ";") else "none")
  all[[st]] <- d
}
if (length(all)) {
  cols <- c("study", "sample_type", "subject", "sample", "GROUP", "TIMEPOINT", "REGION", "GENDER", "CASE_GROUP", "TIME_TO_TB", "SAMPLE_TYPE", "AGE", "SMOKE")
  A <- do.call(rbind, lapply(all, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[, cols] }))
  write.csv(A, file.path(OUT, "GC6-74_metabolomics_sample_table.csv"), row.names = FALSE)
  write.csv(do.call(rbind, summ), file.path(OUT, "GC6-74_metabolomics_metadata_summary.csv"), row.names = FALSE)
  cat("\n========== overlap between studies ==========\n")
  st <- names(all)
  for (i in seq_along(st)) for (j in seq_along(st)) if (i < j) {
    s <- length(intersect(all[[i]]$sample, all[[j]]$sample))
    u <- length(setdiff(intersect(all[[i]]$subject, all[[j]]$subject), c("", "-", NA)))
    cat(sprintf("  %s vs %s: %d shared sample IDs, %d shared participant IDs\n", st[i], st[j], s, u))
  }
}
cat("\n>>> DONE. Outputs in", OUT, "\n")
