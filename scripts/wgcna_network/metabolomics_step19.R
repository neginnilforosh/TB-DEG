## STEP 19: metabolomics layer.
## Downloads a TB study from Metabolomics Workbench, runs differential abundance between two groups
## (the metabolite equivalent of the DEG step), maps the metabolites to KEGG pathways, and joins them
## to the gene side through shared pathways. Outputs in <run>/21_metabolomics/.
##
## FINDING A STUDY (do this first, the script needs a study ID):
##   https://www.metabolomicsworkbench.org/data/browse.php  -> search "tuberculosis"
##   Pick a study with two clear groups (TB vs control) and a reasonable sample size, then set
##   STUDY_ID below. Run with CHECK_ONLY = TRUE first: it prints the study summary and the sample
##   groups without analysing anything, so you can confirm the group labels before committing.

## ---- CONFIG ----
STUDY_ID    <- "ST000974"          # e.g. "ST001234" -- required
CHECK_ONLY  <- FALSE
CASE_LABEL  <- "case"
CTRL_LABEL  <- "control"        # the factor value for control
FACTOR_NAME <- ""          # which factor column to use; "" = the first one with >1 value
## Metabolomics Workbench packs several variables into one string, e.g.
##   "GROUP:case | TIMEPOINT:BL | REGION:MRC | GENDER:F"
## Set SPLIT_FACTORS = TRUE to break that into real columns, then name the one holding the groups.
SPLIT_FACTORS <- TRUE
GROUP_FIELD   <- "GROUP"         # the field inside the compound string that holds case/control
SUBSET <- c(TIMEPOINT = "M06") # keep only these samples; c() = keep everything
COVARIATES    <- c("REGION", "GENDER") # adjusted for, if limma is installed
PADJ_CUT    <- 0.05
LOG2        <- TRUE        # log2-transform before testing (abundances are usually right-skewed)
MIN_FRAC    <- 0.5         # a metabolite must be measured in at least this fraction of samples

req <- c("httr", "jsonlite")
ok <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing package(s): ", paste(req[!ok], collapse = ", "),
                   "\n  install.packages(c(", paste0('"', req[!ok], '"', collapse = ", "), "))")

## ---- paths ----
get_script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "21_metabolomics"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

if (!nzchar(STUDY_ID)) stop("Set STUDY_ID first -- see the header for how to find a TB study.")
API <- "https://www.metabolomicsworkbench.org/rest/study/study_id"
fetch <- function(what) {
  r <- tryCatch(httr::GET(paste(API, STUDY_ID, what, sep = "/"), httr::timeout(120)), error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200)
    stop("Could not reach Metabolomics Workbench for ", STUDY_ID, "/", what,
         " -- check the study ID and your internet connection.")
  txt <- httr::content(r, "text", encoding = "UTF-8")
  out <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (is.null(out)) stop("Unexpected (non-JSON) answer for ", what, ": ", substr(txt, 1, 200))
  out
}
## the API returns either one record or a named list of records; make it a data frame either way
as_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && length(x) && is.list(x[[1]])) return(do.call(rbind, lapply(x, as.data.frame, stringsAsFactors = FALSE)))
  as.data.frame(x, stringsAsFactors = FALSE)
}

## ---- 1. what is in this study? ----
summ <- fetch("summary")
cat(">>> ", STUDY_ID, ": ", as.character(summ$study_title)[1], "\n", sep = "")
cat("    institute:", as.character(summ$institute)[1], "| samples:", as.character(summ$number_of_samples)[1], "\n")

fac <- as_df(fetch("factors"))
## "A:x | B:y" -> columns A and B
if (SPLIT_FACTORS && "factors" %in% names(fac)) {
  parts <- strsplit(as.character(fac$factors), "\\s*\\|\\s*")
  keys  <- unique(unlist(lapply(parts, function(p) sub(":.*$", "", p))))
  for (k in keys) fac[[k]] <- vapply(parts, function(p) {
    hit <- p[sub(":.*$", "", p) == k]; if (length(hit)) sub("^[^:]*:", "", hit[1]) else NA_character_ }, character(1))
  cat(">>> split the compound factor into:", paste(keys, collapse = ", "), "\n")
}
fac_cols <- setdiff(names(fac), c("study_id", "local_sample_id", "sample_id", "subject_id", "Sample ID", "factors", "raw_data", "mb_sample_id"))
cat("\n>>> sample factors available:", paste(fac_cols, collapse = ", "), "\n")
## show every factor column, so a column with only one value (e.g. sample_source = "Blood (plasma)")
## does not send you back for another run; the group labels usually sit in "factors"
for (fc in fac_cols) {
  v <- as.character(fac[[fc]]); n <- length(unique(v))
  cat("\n--- ", fc, " (", n, " distinct value", if (n != 1) "s" else "", ") ---\n", sep = "")
  if (n == 1) cat("    all samples: ", v[1], "  <- cannot split a comparison on this column\n", sep = "")
  else if (n > 30) { cat("    too many values to list; first few:\n"); print(utils::head(sort(table(v), decreasing = TRUE), 10)) }
  else print(sort(table(v), decreasing = TRUE))
}
use_factor <- if (nzchar(FACTOR_NAME)) FACTOR_NAME else if (GROUP_FIELD %in% names(fac)) GROUP_FIELD else {
  multi <- fac_cols[vapply(fac_cols, function(c) length(unique(as.character(fac[[c]]))) > 1, TRUE)]
  if (length(multi)) multi[1] else fac_cols[1]
}
if (!use_factor %in% names(fac)) stop("Factor '", use_factor, "' not found. Available: ", paste(fac_cols, collapse = ", "))
cat("\n>>> using factor column: '", use_factor, "'", if (!nzchar(FACTOR_NAME)) " (first one with more than one value; set FACTOR_NAME to choose another)" else "", "\n", sep = "")
if (length(unique(as.character(fac[[use_factor]]))) < 2)
  cat("!! This column has a single value, so no comparison can be made -- set FACTOR_NAME to another column.\n")
write.csv(fac, file.path(dir_out, paste0(STUDY_ID, "_sample_factors.csv")), row.names = FALSE)

if (CHECK_ONLY) {
  cat("\n>>> CHECK_ONLY is TRUE. Copy the two group labels above into CASE_LABEL and CTRL_LABEL,\n")
  cat("    set CHECK_ONLY <- FALSE, and run again.\n")
  quit(save = "no")
}
if (!nzchar(CASE_LABEL) || !nzchar(CTRL_LABEL)) stop("Set CASE_LABEL and CTRL_LABEL (see the check run).")

## ---- 2. abundance matrix ----
dat <- fetch("data")
dd  <- as_df(dat)
## The API returns one row per metabolite. The per-sample values arrive either as a nested DATA field
## or already flattened into columns named "DATA.<sample>", alongside text columns such as
## analysis_summary and refmet_name -- those must not be mistaken for samples.
META_COLS <- c("study_id", "analysis_id", "analysis_summary", "metabolite_name", "metabolite_id",
               "refmet_name", "units", "DATA")
if ("DATA" %in% names(dd) && is.list(dd$DATA) && !is.data.frame(dd$DATA)) {
  mat <- do.call(rbind, lapply(dd$DATA, function(x) unlist(x)))
} else {
  samp_cols <- grep("^DATA\\.", names(dd), value = TRUE)
  if (!length(samp_cols)) samp_cols <- setdiff(names(dd), META_COLS)
  mat <- as.matrix(dd[, samp_cols, drop = FALSE])
  colnames(mat) <- sub("^DATA\\.", "", colnames(mat))
}
name_src <- intersect(c("metabolite_name", "refmet_name"), names(dd))[1]
if (is.na(name_src)) stop("No metabolite name column; columns are: ", paste(names(dd), collapse = ", "))
rownames(mat) <- make.unique(as.character(dd[[name_src]]))
mode(mat) <- "numeric"
if (all(is.na(mat))) stop("All abundances are missing after conversion -- the data shape is not what this script expects.")
cat("\n>>> abundance matrix:", nrow(mat), "metabolites x", ncol(mat), "samples\n")

id_col <- intersect(c("local_sample_id", "sample_id", "Sample ID"), names(fac))[1]
if (is.na(id_col)) stop("No sample-id column in the factor table; columns are: ", paste(names(fac), collapse = ", "))
grp <- setNames(as.character(fac[[use_factor]]), as.character(fac[[id_col]]))
## keep only the requested slice (e.g. one timepoint)
if (length(SUBSET)) {
  keep_s <- rep(TRUE, nrow(fac))
  for (k in names(SUBSET)) {
    if (!k %in% names(fac)) stop("SUBSET field '", k, "' is not a factor column.")
    keep_s <- keep_s & as.character(fac[[k]]) == SUBSET[[k]]
  }
  cat(">>> subset", paste(names(SUBSET), SUBSET, sep = "=", collapse = ", "), ":", sum(keep_s), "of", nrow(fac), "samples\n")
  fac <- fac[keep_s, , drop = FALSE]
  grp <- grp[as.character(fac[[id_col]])]
}
common <- intersect(colnames(mat), names(grp))
if (!length(common)) stop("Sample names in the data and the factor table do not match.\n  data: ",
                          paste(utils::head(colnames(mat), 3), collapse = ", "), "\n  factors: ",
                          paste(utils::head(names(grp), 3), collapse = ", "))
mat <- mat[, common, drop = FALSE]; g <- grp[common]
keep <- g %in% c(CASE_LABEL, CTRL_LABEL)
if (sum(keep) < 6) stop("Only ", sum(keep), " samples match the two labels -- check the spelling.")
mat <- mat[, keep, drop = FALSE]; g <- factor(g[keep], levels = c(CTRL_LABEL, CASE_LABEL))
cat(">>> analysed:", ncol(mat), "samples |", paste(names(table(g)), table(g), collapse = ", "), "\n")

## ---- 3. differential abundance ----
mat <- mat[rowMeans(!is.na(mat) & mat > 0) >= MIN_FRAC, , drop = FALSE]
cat(">>>", nrow(mat), "metabolites kept (measured in at least", round(100 * MIN_FRAC), "% of samples)\n")
x <- if (LOG2) log2(mat + 1) else mat
covs <- intersect(COVARIATES, names(fac))
cov_df <- if (length(covs)) fac[match(colnames(x), as.character(fac[[id_col]])), covs, drop = FALSE] else NULL
use_limma <- length(covs) > 0 && requireNamespace("limma", quietly = TRUE) &&
             all(vapply(cov_df, function(v) length(unique(v[!is.na(v)])) > 1, TRUE))
if (length(covs) && !use_limma)
  cat("!! Covariates", paste(covs, collapse = ", "), "not used (limma missing, or a covariate has one level);",
      "falling back to an unadjusted test.\n")
if (use_limma) {
  cat(">>> testing with limma, adjusting for:", paste(covs, collapse = ", "), "\n")
  design <- stats::model.matrix(stats::as.formula(paste("~", paste(c(covs, "grp_"), collapse = " + "))),
                                data = cbind(cov_df, grp_ = g))
  fit <- limma::eBayes(limma::lmFit(x, design))
  coef_nm <- grep("^grp_", colnames(design), value = TRUE)[1]
  tt <- limma::topTable(fit, coef = coef_nm, number = Inf, sort.by = "none")
  res <- data.frame(metabolite = rownames(x),
                    mean_case = rowMeans(x[, g == CASE_LABEL, drop = FALSE], na.rm = TRUE),
                    mean_control = rowMeans(x[, g == CTRL_LABEL, drop = FALSE], na.rm = TRUE),
                    log2FC = tt$logFC, pvalue = tt$P.Value, stringsAsFactors = FALSE)
} else
res <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
  a <- x[i, g == CASE_LABEL]; b <- x[i, g == CTRL_LABEL]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) return(NULL)
  w <- suppressWarnings(stats::wilcox.test(a, b))          # non-parametric: abundances are rarely normal
  data.frame(metabolite = rownames(x)[i], mean_case = mean(a), mean_control = mean(b),
             log2FC = mean(a) - mean(b), pvalue = w$p.value, stringsAsFactors = FALSE)
}))
res$padj <- stats::p.adjust(res$pvalue, "BH")
res <- res[order(res$padj, res$pvalue), ]
write.csv(res, file.path(dir_out, paste0(STUDY_ID, "_differential_metabolites.csv")), row.names = FALSE)
sig <- res[!is.na(res$padj) & res$padj < PADJ_CUT, ]
cat(">>> significant metabolites (padj <", PADJ_CUT, "):", nrow(sig), "of", nrow(res), "\n")
if (nrow(sig)) print(utils::head(sig[, c("metabolite", "log2FC", "padj")], 10), row.names = FALSE)

## ---- 4. metabolite -> KEGG pathway ----
## The study's own metabolite table carries database ids where the submitters provided them.
mets <- tryCatch(as_df(fetch("metabolites")), error = function(e) NULL)
kegg_col <- if (!is.null(mets)) intersect(c("kegg_id", "KEGG", "kegg"), names(mets))[1] else NA
map <- NULL
if (!is.na(kegg_col)) {
  name_col <- intersect(c("metabolite_name", "metabolite"), names(mets))[1]
  map <- unique(mets[, c(name_col, kegg_col)]); names(map) <- c("metabolite", "kegg_id")
  cat("\n>>> KEGG ids taken from the study's own metabolite table\n")
} else {
  ## This study has no KEGG ids, only names. RefMet was tried first but its REST answer carries no
  ## KEGG field at all (only refmet_name, formula, mass, classes, refmet_id, pubchem_cid), so the
  ## names are matched against KEGG's own compound index instead.
  cat("\n>>> No KEGG ids in the study table; matching names against KEGG (one request each,\n")
  cat("    so this takes a few minutes)...\n")
  norm_nm <- function(x) tolower(trimws(gsub("\\*+$", "", x)))
  ## KEGG's search is substring-based, so "tryptophan" also returns CMP- and 5-hydroxy- forms.
  ## Only an exact synonym match counts, and the L- form wins when both exist (L-Histidine over
  ## generic Histidine), because that is the form the human pathways use.
  kegg_pick <- function(query, lines) {
    if (!length(lines) || !nzchar(paste(lines, collapse = ""))) return(NA_character_)
    parts <- strsplit(lines, "\t")
    parts <- parts[vapply(parts, length, 1L) >= 2]
    if (!length(parts)) return(NA_character_)
    ids  <- vapply(parts, function(p) sub("^cpd:", "", p[1]), character(1))
    syns <- lapply(parts, function(p) trimws(strsplit(p[2], ";")[[1]]))
    q <- norm_nm(query)
    hit <- which(vapply(syns, function(s) any(norm_nm(s) == q) || any(norm_nm(s) == paste0("l-", q)), TRUE))
    if (!length(hit)) return(NA_character_)
    first_L <- vapply(syns[hit], function(s) grepl("^l-", norm_nm(s[1])), TRUE)
    if (any(first_L)) ids[hit[which(first_L)[1]]] else ids[hit[order(ids[hit])][1]]
  }
  kegg_one <- function(nm) {
    u <- paste0("https://rest.kegg.jp/find/compound/", utils::URLencode(gsub("\\*+$", "", nm), reserved = TRUE))
    r <- tryCatch(httr::GET(u, httr::timeout(30)), error = function(e) NULL)
    if (is.null(r)) return("REQUEST_FAILED")                 # kept distinct from a genuine no-match
    if (httr::status_code(r) != 200) return("REQUEST_FAILED")
    txt <- httr::content(r, "text", encoding = "UTF-8")
    if (!nzchar(trimws(txt))) return(NA_character_)
    kegg_pick(nm, strsplit(txt, "\n")[[1]])
  }
  cache_file <- file.path(dir_out, paste0(STUDY_ID, "_kegg_name_cache.csv"))
  cache <- if (file.exists(cache_file)) read.csv(cache_file, stringsAsFactors = FALSE) else
           data.frame(metabolite = character(0), kegg_id = character(0), stringsAsFactors = FALSE)
  cache <- cache[!(cache$kegg_id %in% "REQUEST_FAILED"), , drop = FALSE]   # retry anything that errored
  todo <- setdiff(res$metabolite, cache$metabolite)
  if (length(todo)) {
    cat("    ", length(todo), "names to look up (", nrow(cache), "already cached)\n")
    got <- vapply(seq_along(todo), function(i) {
      if (i %% 50 == 0) cat("     ...", i, "/", length(todo), "\n")
      Sys.sleep(0.34)                                        # KEGG asks for no more than 3 requests a second
      kegg_one(todo[i])
    }, character(1))
    cache <- rbind(cache, data.frame(metabolite = todo, kegg_id = got, stringsAsFactors = FALSE))
    write.csv(cache, cache_file, row.names = FALSE)
  }
  n_fail <- sum(cache$kegg_id %in% "REQUEST_FAILED")
  if (n_fail) cat("!!", n_fail, "lookups failed to reach KEGG; re-run to retry just those.\n")
  cache$kegg_id[cache$kegg_id %in% "REQUEST_FAILED"] <- NA_character_
  map <- cache
  n_hit <- sum(!is.na(map$kegg_id) & nzchar(map$kegg_id))
  cat(">>> matched", n_hit, "of", nrow(map), "metabolite names to a KEGG compound\n")
  if (n_hit == 0) cat("!! None matched -- the pathway join cannot be made for this study.\n")
}
if (!is.null(map) && any(!is.na(map$kegg_id) & nzchar(map$kegg_id))) {
  res2 <- merge(res, map, by = "metabolite", all.x = TRUE)
  write.csv(res2, file.path(dir_out, paste0(STUDY_ID, "_differential_metabolites_with_KEGG.csv")), row.names = FALSE)
  sig_kegg <- unique(res2$kegg_id[res2$metabolite %in% sig$metabolite & !is.na(res2$kegg_id) & nzchar(res2$kegg_id)])
  cat("\n>>> significant metabolites with a KEGG id:", length(sig_kegg), "\n")

  ## KEGG compound -> pathway, straight from the KEGG REST service
  link <- tryCatch(read.delim("https://rest.kegg.jp/link/pathway/compound", header = FALSE,
                              col.names = c("cpd", "path"), stringsAsFactors = FALSE), error = function(e) NULL)
  nm <- tryCatch(read.delim("https://rest.kegg.jp/list/pathway/hsa", header = FALSE,
                            col.names = c("path", "name"), stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(link) || is.null(nm)) {
    cat("!! Could not reach rest.kegg.jp; skipping the pathway join.\n")
  } else {
    link$cpd  <- sub("^cpd:", "", link$cpd)
    link$pid  <- sub("^path:map", "", link$path)
    nm$pid    <- sub("^(path:)?hsa", "", nm$path)
    met_paths <- merge(link[link$cpd %in% sig_kegg, ], nm, by = "pid")
    met_paths <- unique(met_paths[, c("pid", "name", "cpd")])
    agg <- stats::aggregate(cpd ~ pid + name, met_paths, function(x) paste(sort(unique(x)), collapse = ";"))
    agg$N_metabolites <- vapply(strsplit(agg$cpd, ";"), length, 1L)
    agg <- agg[order(-agg$N_metabolites), ]
    write.csv(agg, file.path(dir_out, paste0(STUDY_ID, "_metabolite_pathways.csv")), row.names = FALSE)
    cat(">>> pathways containing a significant metabolite:", nrow(agg), "\n")
    if (nrow(agg)) print(utils::head(agg[, c("name", "N_metabolites")], 10), row.names = FALSE)

    ## ---- 5. shared pathways with the gene side ----
    gene_kegg <- list.files(file.path(OUT_DIR, "11_module_enrichment"), pattern = "_KEGG\\.csv$", full.names = TRUE)
    if (!length(gene_kegg)) {
      cat("\n!! No module KEGG enrichment found; run the module enrichment step first to join the layers.\n")
    } else {
      gp <- do.call(rbind, lapply(gene_kegg, function(f) {
        d <- read.csv(f, stringsAsFactors = FALSE)
        if (!nrow(d) || !("p.adjust" %in% names(d))) return(NULL)
        d <- d[d$p.adjust < PADJ_CUT, , drop = FALSE]
        if (!nrow(d)) return(NULL)
        data.frame(Module = sub(".*_([a-z]+)_KEGG\\.csv$", "\\1", basename(f)),
                   pid = sub("^hsa", "", d$ID), pathway = d$Description,
                   gene_padj = d$p.adjust, genes = d$geneID, stringsAsFactors = FALSE)
      }))
      if (is.null(gp)) {
        cat("\n!! No module has a significant KEGG pathway, so no shared pathway can be found.\n")
      } else {
        shared <- merge(gp, agg[, c("pid", "name", "cpd", "N_metabolites")], by = "pid")
        shared <- shared[order(shared$gene_padj), ]
        write.csv(shared, file.path(dir_out, paste0(STUDY_ID, "_shared_pathways_genes_metabolites.csv")), row.names = FALSE)
        cat("\n>>> pathways significant on BOTH layers:", nrow(shared), "\n")
        if (nrow(shared)) print(utils::head(shared[, c("pathway", "Module", "gene_padj", "N_metabolites")], 10), row.names = FALSE)
        cat("    These are the links between the transcriptomics and metabolomics layers.\n")
      }
    }
  }
}
cat("\n>>> DONE. Outputs in", dir_out, "\n")
