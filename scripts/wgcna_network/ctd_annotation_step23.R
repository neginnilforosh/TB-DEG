## STEP 23: CTD annotation (Ratul, 10 Oct: "use CTD as a supportive knowledge-annotation resource. For the
## final important genes, annotate known associations with tuberculosis, diabetes, HIV, relevant
## immune/inflammatory diseases, chemicals/drugs. Use CTD as an annotation/validation layer rather than
## as the primary method.")
##
## Queries the CTD Batch Query service for the run's High + Medium consensus genes (curated gene-disease
## associations and chemical-gene interactions), caches every answer, and writes one row per gene with
## yes/no labels and the evidence type. If the service cannot be reached, download the files by hand
## (see MANUAL_* below) and the script uses those instead.
## Outputs in <run>/25_ctd_annotation/.

## ---- CONFIG ----
TIERS <- c("High", "Medium")
DISEASE_LABELS <- c(                                   # label = regular expression on the CTD disease name
  Tuberculosis = "Tuberculosis",
  Diabetes     = "Diabetes Mellitus",
  HIV          = "HIV Infections|Acquired Immunodeficiency Syndrome|HIV",
  Immune_inflammatory = "Inflammation|Arthritis|Lupus|Autoimmune|Sepsis|Psoriasis|Inflammatory Bowel|Crohn|Colitis|Asthma")
BATCH <- 100                                           # genes per request
CTD_URL <- "https://ctdbase.org/tools/batchQuery.go"
## manual fallback: CTD Batch Query web page, input type "Gene", paste the gene list, choose the report,
## download as CSV, and put the files here (leave "" to use the online query)
MANUAL_DISEASES  <- ""   # e.g. "CTD_batch_diseases_curated.csv"
MANUAL_CHEMICALS <- ""   # e.g. "CTD_batch_cgixns.csv"

for (p in c("httr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "25_ctd_annotation"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

pri <- read.csv(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")), stringsAsFactors = FALSE)
genes <- unique(pri$Gene[pri$Final_Priority %in% TIERS & !is.na(pri$Gene) & nzchar(pri$Gene)])
cat(">>> run", run_id, "|", length(genes), paste(TIERS, collapse = "+"), "genes to annotate\n")

## CTD csv: header line starts with "# Input"; column names vary slightly between reports
read_ctd <- function(txt) {
  lines <- strsplit(txt, "\r?\n")[[1]]; lines <- lines[nzchar(lines)]
  h <- grep("^#?\\s*Input", lines)[1]
  if (is.na(h)) return(NULL)
  lines[h] <- sub("^#\\s*", "", lines[h])
  d <- utils::read.csv(text = paste(lines[h:length(lines)], collapse = "\n"), stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- gsub("[^A-Za-z]", "", names(d)); d
}
cache_file <- file.path(dir_out, "ctd_cache.rds")
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
ctd_query <- function(report) {
  out <- list()
  for (b in split(sort(genes), ceiling(seq_along(genes) / BATCH))) {
    key <- paste(report, paste(b, collapse = "|"))
    if (is.null(cache[[key]])) {
      Sys.sleep(1)
      r <- tryCatch(httr::POST(CTD_URL, body = list(inputType = "gene", inputTerms = paste(b, collapse = "|"),
                    report = report, format = "csv"), encode = "form", httr::timeout(120)), error = function(e) NULL)
      if (is.null(r) || httr::status_code(r) != 200)
        stop("CTD request failed (", report, "). Re-run to retry -- finished batches are cached -- or use the MANUAL_* files.")
      d <- read_ctd(httr::content(r, "text", encoding = "UTF-8"))
      if (is.null(d)) stop("CTD returned something that is not a CSV report; first characters:\n",
                           substr(httr::content(r, "text", encoding = "UTF-8"), 1, 300))
      cache[[key]] <<- d; saveRDS(cache, cache_file)
    }
    out[[length(out) + 1]] <- cache[[key]]
  }
  do.call(rbind, out)
}
get_report <- function(report, manual) {
  if (nzchar(manual)) { f <- file.path(dir_out, manual); if (!file.exists(f)) stop("Manual file not found: ", f)
    return(read_ctd(paste(readLines(f, warn = FALSE), collapse = "\n"))) }
  ctd_query(report)
}

## ---- diseases ----
dis <- get_report("diseases_curated", MANUAL_DISEASES)
need <- c("GeneSymbol", "DiseaseName")
if (!all(need %in% names(dis))) stop("Disease report lacks columns ", paste(setdiff(need, names(dis)), collapse = ", "),
                                     "; it has: ", paste(names(dis), collapse = ", "))
ev_col <- intersect(c("DirectEvidence"), names(dis))
cat(">>> CTD curated gene-disease rows:", nrow(dis), "for", length(unique(dis$GeneSymbol)), "genes\n")

ann <- data.frame(Gene = genes, stringsAsFactors = FALSE)
for (lab in names(DISEASE_LABELS)) {
  hit <- dis[grepl(DISEASE_LABELS[[lab]], dis$DiseaseName, ignore.case = TRUE), ]
  ann[[paste0("CTD_", lab)]] <- ann$Gene %in% hit$GeneSymbol
  ann[[paste0("CTD_", lab, "_diseases")]] <- vapply(ann$Gene, function(g)
    paste(utils::head(unique(hit$DiseaseName[hit$GeneSymbol == g]), 5), collapse = "; "), "")
  if (length(ev_col)) ann[[paste0("CTD_", lab, "_evidence")]] <- vapply(ann$Gene, function(g)
    paste(unique(stats::na.omit(hit[[ev_col]][hit$GeneSymbol == g & nzchar(hit[[ev_col]])])), collapse = ";"), "")
}
ann$CTD_any_disease <- ann$Gene %in% dis$GeneSymbol
ann$CTD_N_diseases  <- vapply(ann$Gene, function(g) length(unique(dis$DiseaseName[dis$GeneSymbol == g])), 0L)

## ---- chemicals ----
chem <- tryCatch(get_report("cgixns", MANUAL_CHEMICALS), error = function(e) { cat("!! chemicals skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(chem) && all(c("GeneSymbol", "ChemicalName") %in% names(chem))) {
  if ("Organism" %in% names(chem)) chem <- chem[grepl("Homo sapiens", chem$Organism), ]
  ann$CTD_N_chemicals <- vapply(ann$Gene, function(g) length(unique(chem$ChemicalName[chem$GeneSymbol == g])), 0L)
  ann$CTD_top_chemicals <- vapply(ann$Gene, function(g) {
    t <- sort(table(chem$ChemicalName[chem$GeneSymbol == g]), decreasing = TRUE); paste(utils::head(names(t), 5), collapse = "; ") }, "")
}

out <- merge(pri[pri$Gene %in% genes, c("Gene", "Ensembl", "Module", "Module_Name", "log2FC", "FDR", "Final_Priority")],
             ann, by = "Gene", all.x = TRUE)
out <- out[order(match(out$Final_Priority, TIERS), -out$CTD_Tuberculosis, out$Gene), ]
write.csv(out, file.path(dir_out, paste0(ACC, "_CTD_annotation.csv")), row.names = FALSE)

## ---- summary: is the High tier enriched for known TB genes compared with Medium? ----
cat("\n========== CTD labels by tier ==========\n")
sm <- do.call(rbind, lapply(TIERS, function(t) { s <- out[out$Final_Priority == t, ]
  data.frame(Tier = t, N = nrow(s), t(sapply(names(DISEASE_LABELS), function(l) round(100 * mean(s[[paste0("CTD_", l)]], na.rm = TRUE), 1)))) }))
names(sm)[-(1:2)] <- paste0("Pct_", names(DISEASE_LABELS)); print(sm, row.names = FALSE)
write.csv(sm, file.path(dir_out, paste0(ACC, "_CTD_summary.csv")), row.names = FALSE)
cat("\nKnown TB genes in the High tier:", paste(out$Gene[out$Final_Priority == "High" & out$CTD_Tuberculosis %in% TRUE], collapse = ", "), "\n")
cat("High-tier genes with NO curated TB association (potentially new):",
    paste(out$Gene[out$Final_Priority == "High" & out$CTD_Tuberculosis %in% FALSE], collapse = ", "), "\n")
cat("\nNote: CTD 'curated' associations come from the literature, so well-studied genes collect more labels.\n",
    "A missing label means 'not curated', not 'not involved'.\n", sep = "")
cat(">>> DONE. Outputs in", dir_out, "\n")
