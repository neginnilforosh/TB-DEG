# Reproducibility record

Generated with R 4.5.0 on macOS 27.0.1.

Package versions: `package_versions.csv`. Full session: `sessionInfo.txt`.

## Annotation and database versions

| Resource | Version |
|---|---|
| org.Hs.eg.db | 3.22.0 |
| msigdbr | 26.1.1 |
| msigdbdf | not installed |
| STRINGdb | not installed |
| clusterProfiler | 4.18.4 |
| ChIPseeker | 1.46.1 |
| TxDb.Hsapiens.UCSC.hg38.knownGene | 3.22.0 |
| drugfindR | 1.0.0 |
| WGCNA | 1.74 |
| DESeq2 | 1.50.2 |

org.Hs.eg.db source dates: EGSOURCEDATE = 2025-Sep24; GOSOURCEDATE = 2025-07-22; ENSOURCEDATE = 2025-Sep03

Web resources queried at run time (versions are those live on the run date): STRING, iLINCS,
Enrichr, DGIdb, KEGG REST, Metabolomics Workbench. Query dates are not stored in the public outputs;
see reproducibility/external_resource_versions.csv for which resources are version-pinned.

## Run presets (scripts/wgcna_network/config.R)

```r
PRESETS <- list(
  GSE114192_TB = list(                       # pilot: Healthy_Control vs TB_Only
    ACC = "GSE114192", CONTROL = "Healthy_Control", CASE = "TB_Only",
    DEG_NAME = "HealthyControl_vs_TBOnly",     # <ACC>_<DEG_NAME>_DEG.csv in 06_deg_results
    RUN = "",                                  # "" = outputs directly in the dataset folder
    MODULE_NAMES = c(green = "Interferon_Response", purple = "Ribosome_Biogenesis_Chromatin", blue = "Mixed_CellComposition"),
    MODULE_PLOT_COLORS = c(green = "#2ecc71", purple = "#9b59b6", blue = "#3498db")),
  
  GSE161829_ATB = list(                      # replication: TBneg vs active TB (LTBI samples excluded)
    ACC = "GSE161829", CONTROL = "TBneg", CASE = "ATB",
    DEG_NAME = "TBneg_vs_ATB",
    RUN = "TBneg_vs_ATB",                      # outputs in GSE161829/TBneg_vs_ATB/
    MODULE_NAMES = c(midnightblue = "M1_Cell_Cycle_Proliferation",
                greenyellow = "M2_Mixed_with_IFN_cluster",
                 salmon       = "M3_OxPhos_weak"),              # fill in after step 2
    MODULE_PLOT_COLORS = NULL) ,                # NULL = use the WGCNA colour names
  
  GSE114192_TBDM_vs_DM = list(
    ACC = "GSE114192", CONTROL = "DM_only", CASE = "TB_DM",
    DEG_NAME = "DM_only_vs_TB_DM",
    RUN = "TB_DM",
    MODULE_NAMES = c(black = "M1_Interferon_Complement",
                 greenyellow = "M2_TLR_NFkB_Inflammation",
                 green = "M3_Mixed_CellComposition"),
    MODULE_PLOT_COLORS = NULL),
	
  GSE114192_TBDM_vs_TB = list(
	  ACC = "GSE114192", CONTROL = "TB_Only", CASE = "TB_DM",
    DEG_NAME = "TB_Only_vs_TB_DM", 
    RUN = "TB_DM",
    MODULE_NAMES = character(0), 
    MODULE_PLOT_COLORS = NULL)

)
```

## Random seeds

| Script | Line | Call |
|---|---|---|
| scripts/deg_pipeline/site_sensitivity_reference.R | 102 | `set.seed(SEED); g <- clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g,` |
| scripts/deg_pipeline/site_sensitivity_tbdm.R | 108 | `set.seed(SEED); g <- as.data.frame(clusterProfiler::GSEA(sort(setNames(d$stat, d$SYMBOL), decreasing = TRUE), TERM2GENE = t2g,` |
| scripts/wgcna_network/figures_step24.R | 159 | `set.seed(42); lay <- igraph::layout_with_fr(g)` |
| scripts/wgcna_network/figures_step24.R | 248 | `g <- igraph::graph_from_data_frame(d[, 1:2], directed = FALSE); set.seed(42); lay <- igraph::layout_with_fr(g)` |
| scripts/wgcna_network/GSE114192_interaction_gsea.R | 55 | `set.seed(SEED)` |
| scripts/wgcna_network/GSE165709_atac_geneset_step18.R | 70 | `set.seed(42)` |
| scripts/wgcna_network/ilincs_proliferation_sensitivity.R | 69 | `set.seed(SEED); null <- replicate(N_PERM, length(intersect(sample(univ, n), sample(univ, n))))` |
| scripts/wgcna_network/ml_validation_step22.R | 58 | `boot_ci <- function(score, y) { set.seed(SEED); i1 <- which(y == 1); i0 <- which(y == 0)` |
| scripts/wgcna_network/ml_validation_step22.R | 80 | `set.seed(SEED)` |
| scripts/wgcna_network/ml_validation_step22.R | 91 | `set.seed(SEED); pool <- setdiff(common, genes)` |
| scripts/wgcna_network/ml_validation_step22.R | 128 | `set.seed(SEED)` |
| scripts/wgcna_network/WGCNA_combined_network_step5.R | 172 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_GSEA_step11.R | 102 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_network_algorithms_step9.R | 151 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_STRING_step4.R | 132 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_STRING_step4.R | 161 | `set.seed(42)` |

## Parameters by script

### scripts/deg_pipeline/GSE114192_comorbidity_run.R

```r
ACC       <- "GSE114192"
RUN       <- "TB_DM"                       # output subfolder
LFC_TH    <- 1
PADJ_TH   <- 0.05
GROUP_COL <- "disease state (disease_category):ch1"
## every group is kept here; the comparisons below decide what is actually used
GROUP_MAP <- c("Healthy_Control" = "Healthy_Control",
               "TB_only"         = "TB_Only",
               "TB_DM"           = "TB_DM",
               "DM_only"         = "DM_only",
               "TB_IH"           = "Exclude",     # intermediate hyperglycaemia: small and hard to interpret
               "IH"              = "Exclude")
## case first, control second -- a positive log2FC means higher in the case group
COMPARISONS <- list(
  c("TB_DM", "DM_only"),           # MAIN: what TB adds on a diabetic background
  c("TB_DM", "TB_Only"),           # what diabetes adds on top of TB
  c("TB_DM", "Healthy_Control"),   # the full comorbid state vs healthy
  c("DM_only", "Healthy_Control")  # diabetes alone, as a reference for the others
)
## Expression ~ TB + DM + TB:DM over all four groups. The interaction term finds genes whose
## behaviour in the comorbid state is not just the sum of the two separate effects, which no
## pairwise contrast can show. It needs all four groups present with enough samples each.
RUN_INTERACTION <- TRUE
MIN_PER_GROUP   <- 3     # a group below this is dropped from a comparison
MIN_FOR_WGCNA   <- 15    # below this, a per-context WGCNA is not worth running
suppressPackageStartupMessages({
  library(DESeq2); library(GEOquery); library(stringr)
})
SCRIPT_DIR <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(local({ d <- SCRIPT_DIR
  hit <- list.files(dirname(d), pattern = "^00_functions\\.R$", recursive = TRUE, full.names = TRUE)
  if (!length(hit)) hit <- list.files(d, pattern = "^00_functions\\.R$", full.names = TRUE)
  if (!length(hit)) stop("00_functions.R not found"); hit[1] }))
BASE_DIR <- local({ d <- SCRIPT_DIR
  while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); file.path(d, ACC) })
dir_raw  <- file.path(BASE_DIR, "01_raw_counts")          # read from the shared raw folder
OUT      <- file.path(BASE_DIR, RUN)                      # everything else written here
dirs <- list(meta = "02_metadata", filt = "03_filtered_counts", norm = "04_normalized",
             corr = "05_sample_correlation", deg = "06_deg_results", sig = "07_significant_degs",
             wgcna = "09_wgcna_input", fig = "figures")
for (d in dirs) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
P <- function(key, ...) file.path(OUT, dirs[[key]], paste0(...))
cat(">>> writing to", OUT, "-- the existing", ACC, "outputs are not touched\n")
```

### scripts/deg_pipeline/GSE114192_run.R

```r
GROUP_COL <- "disease state (disease_category):ch1"
GROUP_MAP <- c(
  "Healthy_Control" = "Healthy_Control",
  "TB_only"         = "TB_Only",
  "TB_DM"           = "Exclude",
  "DM_only"         = "Exclude",
  "TB_IH"           = "Exclude",
  "IH"              = "Exclude"
)
```

### scripts/deg_pipeline/GSE161829_run.R

```r
GROUP_COL   <- "disease group:ch1"
COUNT_FILE  <- "GSE161829_RawCounts_CD4m.txt.gz"
GROUP_MAP <- c(
  "TBneg"              = "TBneg",
  "LTBI"               = "LTBI",
  "ATB_Diagnosis"      = "ATB",
  "ATB_Post_treatment" = "Exclude"
)
```

### scripts/deg_pipeline/GSE222001_run.R

```r
GROUP_COL   <- "disease state:ch1"        # <-- from A1
COUNT_FILE  <- "GSE222001_count.txt.gz" # <-- from A2 (GEO supplies only the gzipped file; read.delim reads it directly)
## Map the raw labels GEO uses onto clean names.
GROUP_MAP <- c(
  "Healthy control"                = "Healthy",
  "Active tuberculosis infection"  = "Active",
  "Latent tuberculosis infection"  = "Latent"
)
```

### scripts/deg_pipeline/site_sensitivity_reference.R

```r
PADJ <- 0.05; LFC <- 1; SEED <- 42
SAME_SITES <- c("South_Africa", "Romania")
for (p in c("DESeq2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(DESeq2))
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, "GSE114192"); OUT <- file.path(B, "sensitivity_site"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
```

### scripts/deg_pipeline/site_sensitivity_tbdm.R

```r
PADJ <- 0.05; LFC <- 1; SEED <- 42
SAME_SITES <- c("South_Africa", "Romania")
for (p in c("DESeq2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(DESeq2))
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
B <- file.path(ROOT, "GSE114192", "TB_DM"); OUT <- file.path(B, "sensitivity_site"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("Missing: ", f); read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
```

### scripts/wgcna_network/ctd_annotation_step23.R

```r
TIERS <- c("High", "Medium")
DISEASE_LABELS <- c(                                   # label = regular expression on the CTD disease name
  Tuberculosis = "Tuberculosis",
  Mycobacterial = "Mycobacteri|Tuberculosis|BCG",      # wider: also inherited susceptibility to mycobacterial disease
  Diabetes     = "Diabetes Mellitus",
  HIV          = "HIV Infections|Acquired Immunodeficiency Syndrome|HIV",
  Immune_inflammatory = "Inflammation|Arthritis|Lupus|Autoimmune|Sepsis|Psoriasis|Inflammatory Bowel|Crohn|Colitis|Asthma")
BATCH <- 100                                           # genes per request
CTD_URL <- "https://ctdbase.org/tools/batchQuery.go"
## manual fallback: CTD Batch Query web page, input type "Gene", paste the gene list, choose the report,
## download as CSV, and put the files here (leave "" to use the online query)
MANUAL_DISEASES  <- ""   # e.g. "CTD_batch_diseases_curated.csv"
MANUAL_CHEMICALS <- ""   # e.g. "CTD_batch_cgixns.csv"
## RECOMMENDED: CTD blocks scripted requests (they are redirected to a browser check), so download the two
## bulk files once with a browser from https://ctdbase.org/downloads/ and put them in external_data/CTD/:
##   "Gene-disease associations" -> CTD_curated_genes_diseases.csv.gz   (curated only, small)
##   "Chemical-gene interactions" -> CTD_chem_gene_ixns.csv.gz          (large; optional)
## They are then used for every run (pilot, TB-DM, TB-HIV) and take precedence over the online query.
BULK_DISEASES  <- "CTD_curated_genes_diseases.csv.gz"
BULK_CHEMICALS <- "CTD_chem_gene_ixns.csv.gz"
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
## CTD bulk files: comment lines start with '#'; the column names are on the line after "# Fields:"
read_bulk <- function(path, fallback) {
  con <- gzfile(path); head <- readLines(con, n = 60); close(con)
  com <- grep("^#", head); hdr <- NULL
  fi <- grep("^#\\s*Fields:", head)
  if (length(fi) && fi[1] < length(head)) hdr <- strsplit(sub("^#\\s*", "", head[fi[1] + 1]), ",")[[1]]
  if (is.null(hdr) || length(hdr) < 3) hdr <- fallback
  skip <- if (length(com)) max(com) else 0
  d <- if (requireNamespace("data.table", quietly = TRUE))
    as.data.frame(data.table::fread(path, skip = skip, header = FALSE, sep = ",", quote = "\"", fill = TRUE, showProgress = FALSE))
  else utils::read.csv(gzfile(path), skip = skip, header = FALSE, stringsAsFactors = FALSE, quote = "\"")
  names(d)[seq_len(min(length(hdr), ncol(d)))] <- hdr[seq_len(min(length(hdr), ncol(d)))]
  names(d) <- gsub("[^A-Za-z]", "", names(d))
  d <- d[d$GeneSymbol %in% genes, , drop = FALSE]
  ## the full (non-curated) gene-disease file also holds inferred rows; keep curated ones only
  if ("DirectEvidence" %in% names(d)) d <- d[!is.na(d$DirectEvidence) & nzchar(d$DirectEvidence), , drop = FALSE]
  d
}
BULK_DIR <- file.path(ROOT_DIR, "external_data", "CTD")
get_report <- function(report, manual) {
  bulk <- file.path(BULK_DIR, if (report == "diseases_curated") BULK_DISEASES else BULK_CHEMICALS)
  if (file.exists(bulk)) {
    cat(">>> using CTD bulk file", basename(bulk), "\n")
    return(read_bulk(bulk, if (report == "diseases_curated")
      c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceChemicalName", "InferenceScore", "OmimIDs", "PubMedIDs") else
      c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID", "GeneForms", "Organism", "OrganismID", "Interaction", "InteractionActions", "PubMedIDs")))
  }
  if (nzchar(manual)) { f <- file.path(dir_out, manual); if (!file.exists(f)) stop("Manual file not found: ", f)
    return(read_ctd(paste(readLines(f, warn = FALSE), collapse = "\n"))) }
  ctd_query(report)
}
```

### scripts/wgcna_network/evidence_family_consensus_step25.R

```r
MIN_FAMILIES_GENE <- 3      # a gene needs a score in at least this many families
MIN_FAMILIES_DRUG <- 3      # a drug needs reversal evidence plus at least this many families in total
STRONG    <- 0.80           # a family "supports" a gene/drug when its score is in the top 20%
HIGH_RULE <- list(strong_frac = 0.75, min_score = 0.75)   # proposed: strong in >= 75% of the families (genes 4 of 5, drugs 3 of 4)
MED_RULE  <- list(min_strong = 2, min_score = 0.60)    # change here; everything is recomputed
TOP_N_GENES <- 50; TOP_N_DRUGS <- 20
CAP_LOG10 <- 10             # cap -log10(p) so a handful of extreme p-values do not dominate
ACTIVE_PHENOTYPE <- "Active"
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "20_evidence_families"); dir.create(dir_out, showWarnings = FALSE)
rd <- function(f, need = TRUE) { if (!file.exists(f)) { if (need) stop("Missing: ", f) else return(NULL) }
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
pct <- function(x) { r <- rep(NA_real_, length(x)); ok <- !is.na(x); if (sum(ok) > 1) r[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1) else r[ok] <- 1; r }
nlog <- function(p) pmin(-log10(pmax(p, 1e-300)), CAP_LOG10)
rmean <- function(...) { m <- cbind(...); out <- rowMeans(m, na.rm = TRUE); out[rowSums(!is.na(m)) == 0] <- NA; out }
jac <- function(a, b) if (!length(union(a, b))) NA else length(intersect(a, b)) / length(union(a, b))
## ======================= GENES =======================
g <- rd(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
g$Fam_Expression   <- rmean(pct(abs(g$log2FC)), pct(nlog(g$FDR)))
g$Fam_Coexpression <- rmean(pct(abs(g$MM)), pct(abs(g$Gene_Significance)))
g$Fam_Network      <- rmean(pct(g$PPI_Degree), pct(g$PPI_Betweenness), pct(rmean(pct(g$RWR_score), pct(g$Diffusion_score))))
g$Fam_Network[!g$In_PPI_network %in% TRUE] <- NA
## replication: independent active-TB cohorts (step 21). Agreement = sign(match) x -log10 padj, averaged over cohorts
ev <- rd(file.path(OUT_DIR, "23_external_validation", paste0(ACC, "_external_validation_long.csv")), need = FALSE)
if (!is.null(ev)) {
  ev <- ev[ev$Phenotype == ACTIVE_PHENOTYPE & !is.na(ev$Ext_log2FC) & !is.na(ev$Ext_padj), ]
  ev$agree <- ifelse(sign(ev$Ext_log2FC) == sign(ev$Discovery_log2FC), 1, -1) * nlog(ev$Ext_padj)
  a <- tapply(ev$agree, ev$Ensembl, mean); g$Fam_Replication <- pct(unname(a[g$Ensembl]))
} else g$Fam_Replication <- NA_real_
## infection response: GSE165709 macrophages from healthy donors, infected vs uninfected
inf <- rd(file.path(ROOT_DIR, "external_results", "GSE165709", "integration", "GSE165709_gene_integration.csv"), need = FALSE)
if (!is.null(inf)) {
  inf <- inf[!is.na(inf$RNA_log2FC_HC) & !is.na(inf$RNA_padj_HC), ]; inf <- inf[!duplicated(inf$Gene), ]
  m <- match(g$Gene, inf$Gene)
  agree <- ifelse(sign(inf$RNA_log2FC_HC[m]) == sign(g$log2FC), 1, -1) * nlog(inf$RNA_padj_HC[m])
  g$Fam_Infection <- pct(agree)
} else g$Fam_Infection <- NA_real_
gfam <- c("Fam_Expression", "Fam_Coexpression", "Fam_Network", "Fam_Replication", "Fam_Infection")
gfam <- gfam[vapply(gfam, function(f) any(!is.na(g[[f]])), logical(1))]
score_genes <- function(fams, strong = STRONG) {
  M <- as.matrix(g[, fams, drop = FALSE]); n_av <- rowSums(!is.na(M))
  s <- rowMeans(M, na.rm = TRUE); s[n_av < min(MIN_FAMILIES_GENE, length(fams))] <- NA
  list(score = s, n_available = n_av, n_strong = rowSums(M >= strong, na.rm = TRUE))
}
tier <- function(score, n_strong, n_fam) ifelse(is.na(score), "Insufficient evidence",
  ifelse(n_strong >= ceiling(HIGH_RULE$strong_frac * n_fam) & score >= HIGH_RULE$min_score, "High",
  ifelse(n_strong >= MED_RULE$min_strong | score >= MED_RULE$min_score, "Medium", "Low")))
main <- score_genes(gfam)
g$Families_available <- main$n_available; g$Families_strong <- main$n_strong
g$Evidence_score <- round(main$score, 4); g$Evidence_rank <- rank(-g$Evidence_score, ties.method = "min", na.last = "keep")
g$Evidence_tier <- tier(g$Evidence_score, g$Families_strong, length(gfam))
## sensitivity: leave one family out, and the "strong" threshold at 0.7 / 0.9
variants <- c(list(main = gfam), setNames(lapply(gfam, function(f) setdiff(gfam, f)), paste0("without_", sub("Fam_", "", gfam))))
top_of <- function(s, n) { o <- order(-s, na.last = NA); g$Ensembl[utils::head(o, n)] }
main_top <- top_of(g$Evidence_score, TOP_N_GENES); main_high <- g$Ensembl[g$Evidence_tier == "High"]
ranks <- sapply(variants, function(v) rank(-score_genes(v)$score, ties.method = "min", na.last = "keep"))
sens <- do.call(rbind, lapply(names(variants), function(v) {
  sv <- score_genes(variants[[v]]); t <- tier(sv$score, sv$n_strong, length(variants[[v]]))
  data.frame(Variant = v, Families = paste(sub("Fam_", "", variants[[v]]), collapse = "+"),
             Spearman_with_main = round(suppressWarnings(stats::cor(sv$score, g$Evidence_score, method = "spearman", use = "pairwise")), 3),
             TopN_overlap_with_main = sum(top_of(sv$score, TOP_N_GENES) %in% main_top), N_High = sum(t == "High"),
             High_Jaccard_with_main = round(jac(g$Ensembl[t == "High"], main_high), 3))
}))
for (thr in c(0.70, 0.90)) { sv <- score_genes(gfam, thr); t <- tier(sv$score, sv$n_strong, length(gfam))
  sens <- rbind(sens, data.frame(Variant = paste0("strong_threshold_", thr), Families = paste(sub("Fam_", "", gfam), collapse = "+"),
    Spearman_with_main = 1, TopN_overlap_with_main = TOP_N_GENES, N_High = sum(t == "High"),
    High_Jaccard_with_main = round(jac(g$Ensembl[t == "High"], main_high), 3))) }
g$Rank_best_across_variants  <- apply(ranks, 1, function(r) suppressWarnings(min(r, na.rm = TRUE)))
g$Rank_worst_across_variants <- apply(ranks, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
g$Stable_topN <- apply(ranks, 1, function(r) all(!is.na(r) & r <= TOP_N_GENES))
g$Previous_tier_step15 <- g$Final_Priority
keep <- c("Gene", "Ensembl", "Module", "Module_Name", "log2FC", "FDR", gfam, "Families_available", "Families_strong",
          "Evidence_score", "Evidence_rank", "Rank_best_across_variants", "Rank_worst_across_variants", "Stable_topN",
          "Evidence_tier", "Previous_tier_step15", "DEG_WGCNA_overlap")
go <- g[order(g$Evidence_rank, na.last = TRUE), keep]
write.csv(go, file.path(dir_out, paste0(ACC, "_Gene_evidence_families.csv")), row.names = FALSE)
write.csv(sens, file.path(dir_out, paste0(ACC, "_Gene_evidence_sensitivity.csv")), row.names = FALSE)
cat(">>> run", run_id, "| GENES:", nrow(g), "| families:", paste(sub("Fam_", "", gfam), collapse = ", "), "\n")
cat("    tiers:", paste(names(table(g$Evidence_tier)), table(g$Evidence_tier), collapse = ", "),
    "| stable in the top", TOP_N_GENES, "in every variant:", sum(g$Stable_topN), "\n")
print(table(Step15 = g$Final_Priority, Families = g$Evidence_tier))
cat("\n  top 15 genes:\n"); print(utils::head(go[, c("Gene", "Module", "Evidence_score", "Families_strong", "Rank_best_across_variants",
                                                       "Rank_worst_across_variants", "Evidence_tier", "Previous_tier_step15")], 15), row.names = FALSE)
cat("\n  sensitivity:\n"); print(sens, row.names = FALSE)
## ======================= DRUGS =======================
d <- rd(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_DrugPriority.csv")))
d$Fam_Reversal <- pct(-d$iLINCS_Reversal_Score)
d$Fam_Target_network <- rmean(pct(d$N_TB_Network_Targets), pct(d$Target_fraction_in_network))
d$Fam_Target_network[!(d$N_TB_Network_Targets > 0)] <- NA
gene_score <- setNames(g$Evidence_score, toupper(g$Gene))
best_target <- vapply(strsplit(ifelse(is.na(d$Known_Targets_used), "", d$Known_Targets_used), ";"), function(t) {
  s <- gene_score[toupper(t)]; if (all(is.na(s))) NA_real_ else max(s, na.rm = TRUE) }, 0)
d$Best_target_gene_score <- round(best_target, 4)
d$Fam_Target_genes <- rmean(pct(best_target), pct(d$RWR_score))
d$Fam_Module_overlap <- pct(nlog(d$Module_Enrichment_P))
d$Fam_Proximity <- pct(-d$Network_Proximity_Z)
dfam <- c("Fam_Reversal", "Fam_Target_network", "Fam_Target_genes", "Fam_Module_overlap")
score_drugs <- function(fams, strong = STRONG) {
  M <- as.matrix(d[, fams, drop = FALSE]); n_av <- rowSums(!is.na(M)); s <- rowMeans(M, na.rm = TRUE)
  ok <- n_av >= min(MIN_FAMILIES_DRUG, length(fams)) & (!"Fam_Reversal" %in% fams | !is.na(d$Fam_Reversal)); s[!ok] <- NA
  list(score = s, n_available = n_av, n_strong = rowSums(M >= strong, na.rm = TRUE))
}
dm <- score_drugs(dfam)
d$Families_available <- dm$n_available; d$Families_strong <- dm$n_strong
d$Evidence_score <- round(dm$score, 4); d$Evidence_rank <- rank(-d$Evidence_score, ties.method = "min", na.last = "keep")
d$Evidence_tier <- tier(d$Evidence_score, d$Families_strong, length(dfam))
dtop <- function(s, n) { o <- order(-s, na.last = NA); d$Drug[utils::head(o, n)] }
d_main_top <- dtop(d$Evidence_score, TOP_N_DRUGS)
dvariants <- c(list(main = dfam, with_Proximity = c(dfam, "Fam_Proximity")),
               setNames(lapply(dfam, function(f) setdiff(dfam, f)), paste0("without_", sub("Fam_", "", dfam))))
dranks <- sapply(dvariants, function(v) rank(-score_drugs(v)$score, ties.method = "min", na.last = "keep"))
dsens <- do.call(rbind, lapply(names(dvariants), function(v) { sv <- score_drugs(dvariants[[v]])
  data.frame(Variant = v, Families = paste(sub("Fam_", "", dvariants[[v]]), collapse = "+"), Drugs_scored = sum(!is.na(sv$score)),
             Spearman_with_main = round(suppressWarnings(stats::cor(sv$score, d$Evidence_score, method = "spearman", use = "pairwise")), 3),
             TopN_overlap_with_main = sum(dtop(sv$score, TOP_N_DRUGS) %in% d_main_top)) }))
d$Rank_best_across_variants  <- apply(dranks, 1, function(r) suppressWarnings(min(r, na.rm = TRUE)))
d$Rank_worst_across_variants <- apply(dranks, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
d$Stable_topN <- apply(dranks, 1, function(r) all(!is.na(r) & r <= TOP_N_DRUGS))
d$Rank_best_across_variants[!is.finite(d$Rank_best_across_variants)] <- NA; d$Rank_worst_across_variants[!is.finite(d$Rank_worst_across_variants)] <- NA
d$Previous_tier_step15 <- d$Final_Priority
dkeep <- c("Drug", "iLINCS_Reversal_Score", "N_TB_Network_Targets", "Target_fraction_in_network", "Known_Targets_used",
           "Best_target_gene_score", "Module_Enrichment_P", "Network_Proximity_Z", dfam, "Fam_Proximity", "Families_available",
           "Families_strong", "Evidence_score", "Evidence_rank", "Rank_best_across_variants", "Rank_worst_across_variants",
           "Stable_topN", "Evidence_tier", "Previous_tier_step15")
do_ <- d[order(d$Evidence_rank, na.last = TRUE), dkeep]
write.csv(do_, file.path(dir_out, paste0(ACC, "_Drug_evidence_families.csv")), row.names = FALSE)
write.csv(dsens, file.path(dir_out, paste0(ACC, "_Drug_evidence_sensitivity.csv")), row.names = FALSE)
cat("\n>>> DRUGS:", nrow(d), "| scored (reversal + at least", MIN_FAMILIES_DRUG, "families):", sum(!is.na(d$Evidence_score)),
    "| tiers:", paste(names(table(d$Evidence_tier)), table(d$Evidence_tier), collapse = ", "),
    "| stable in the top", TOP_N_DRUGS, "in every variant:", sum(d$Stable_topN), "\n")
print(table(Step15 = d$Final_Priority, Families = d$Evidence_tier))
cat("\n  top 15 drugs:\n"); print(utils::head(do_[, c("Drug", "Evidence_score", "Families_strong", "Rank_best_across_variants",
                                                       "Rank_worst_across_variants", "Evidence_tier", "Previous_tier_step15")], 15), row.names = FALSE)
cat("\n  sensitivity:\n"); print(dsens, row.names = FALSE)
writeLines(c(paste0("Evidence-family consensus for ", ACC, " (", CONTROL, " vs ", CASE, ")."),
  paste0("Gene families: ", paste(sub("Fam_", "", gfam), collapse = ", "), ". Drug families: ", paste(sub("Fam_", "", dfam), collapse = ", "),
         "; proximity is reported and tested as a sensitivity variant only."),
  "Each family score is a within-run percentile (0 = weakest, 1 = strongest); the evidence score is the mean over available families.",
  paste0("A family supports a gene or drug when its score is >= ", STRONG, ". Proposed tiers: High = support from at least ",
         round(100 * HIGH_RULE$strong_frac), "% of the families (genes ", ceiling(HIGH_RULE$strong_frac * length(gfam)), " of ", length(gfam),
         ", drugs ", ceiling(HIGH_RULE$strong_frac * length(dfam)), " of ", length(dfam), ") and score >= ", HIGH_RULE$min_score, "; Medium = at least ", MED_RULE$min_strong,
         " supporting families or score >= ", MED_RULE$min_score, "; otherwise Low. Genes need ", MIN_FAMILIES_GENE,
         " families with data, drugs need reversal evidence and ", MIN_FAMILIES_DRUG, " families."),
  paste0("Rank_best/worst_across_variants: range of the rank when each family is left out in turn (and, for drugs, when proximity is added). ",
         "Stable_topN = in the top ", TOP_N_GENES, " genes / ", TOP_N_DRUGS, " drugs in every variant.")),
  file.path(dir_out, paste0(ACC, "_EVIDENCE_FAMILIES_NOTES.txt")))
cat("\n>>> DONE. Outputs in", dir_out, "\n")
```

### scripts/wgcna_network/external_validation_step21.R

```r
## comparisons: dataset, DEG file stem, phenotype, note. log2FC sign must mean "higher in disease".
EXTERNAL <- list(
  list(ds = "GSE161829", stem = "TBneg_vs_ATB",       pheno = "Active", note = "CD4 memory T cells (GSE161829_RawCounts_CD4m), not whole blood, n = 20 vs 25"),
  list(ds = "GSE222001", stem = "Healthy_vs_Active",  pheno = "Active", note = "n = 5 vs 5 -- low power"),
  list(ds = "GSE161829", stem = "TBneg_vs_LTBI",      pheno = "LTBI",   note = "CD4 memory T cells (GSE161829_RawCounts_CD4m), not whole blood, n = 20 vs 40"),
  list(ds = "GSE222001", stem = "Healthy_vs_Latent",  pheno = "LTBI",   note = "n = 5 vs 5 -- low power"),
  list(ds = "GSE99374",  stem = "HC_vs_LTBI",         pheno = "LTBI",   note = "sorted CD8 T cells, not whole blood"))
PADJ <- 0.05
TIERS <- c("High", "Medium", "Low")
suppressPackageStartupMessages(library(ggplot2))
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
OUT_DIR  <- run_out_dir(file.path(ROOT_DIR, ACC))
dir_out  <- file.path(OUT_DIR, "23_external_validation"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
FIG_DPI <- 600
save_gg <- function(f, p, w, h) { ggsave(f, p, width = w, height = h, dpi = FIG_DPI); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h) }
pri <- read.csv(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")), stringsAsFactors = FALSE)
pri <- pri[!is.na(pri$log2FC), ]
pri$Disc_sign <- sign(pri$log2FC)
cat(">>> discovery run:", run_id, "|", nrow(pri), "genes |", paste(names(table(pri$Final_Priority)), table(pri$Final_Priority), collapse = ", "), "\n")
## one key per external row: Ensembl without version when there is one, otherwise the symbol
key_of <- function(id) {
  ens <- ifelse(grepl("^ENSG", id), sub("\\..*$", "", sub("\\|.*$", "", id)), NA)
  sym <- ifelse(grepl("\\|", id), sub("^.*\\|", "", id), ifelse(grepl("^ENSG", id), NA, id))
  data.frame(ens = ens, sym = sym, stringsAsFactors = FALSE)
}
long <- list()
for (e in EXTERNAL) {
  f <- file.path(ROOT_DIR, e$ds, "06_deg_results", paste0(e$ds, "_", e$stem, "_DEG.csv"))
  if (!file.exists(f)) { cat("!! missing, skipped:", f, "\n"); next }
  x <- read.csv(f, stringsAsFactors = FALSE); k <- key_of(x$gene); x <- cbind(x, k)
  by_ens <- x[!is.na(x$ens), ]; by_ens <- by_ens[!duplicated(by_ens$ens), ]
  by_sym <- x[!is.na(x$sym), ];  by_sym <- by_sym[order(by_sym$pvalue), ]; by_sym <- by_sym[!duplicated(by_sym$sym), ]
  m <- match(pri$Ensembl, by_ens$ens); use <- by_ens[m, ]
  miss <- is.na(m) & !is.na(pri$Gene); ms <- match(pri$Gene[miss], by_sym$sym)
  use[miss, ] <- by_sym[ms, ]
  comp <- paste(e$ds, e$stem, sep = ":")
  long[[comp]] <- data.frame(Gene = pri$Gene, Ensembl = pri$Ensembl, Module = pri$Module,
    Final_Priority = pri$Final_Priority, Discovery_log2FC = pri$log2FC, Comparison = comp, Phenotype = e$pheno,
    Ext_log2FC = use$log2FoldChange, Ext_padj = use$padj, stringsAsFactors = FALSE)
  cat(sprintf("  %-32s %-6s matched %5d / %d genes  (%s)\n", comp, e$pheno, sum(!is.na(use$log2FoldChange)), nrow(pri), e$note))
}
L <- do.call(rbind, long)
L$Same_direction <- ifelse(is.na(L$Ext_log2FC), NA, sign(L$Ext_log2FC) == sign(L$Discovery_log2FC))
L$Reproduced <- !is.na(L$Ext_padj) & L$Ext_padj < PADJ & L$Same_direction %in% TRUE
write.csv(L, file.path(dir_out, paste0(ACC, "_external_validation_long.csv")), row.names = FALSE)
```

### scripts/wgcna_network/figures_step24.R

```r
FIG_DPI <- 600
TOP_N   <- 20            # bars in "top" plots
TOP_TERMS <- 8           # terms per module in enrichment dot plots
BASE_SIZE <- 9
for (p in c("ggplot2", "igraph")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages({ library(ggplot2) })
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR <- file.path(ROOT_DIR, ACC); OUT_DIR <- run_out_dir(BASE_DIR)
FIG <- file.path(OUT_DIR, "figures_manuscript"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
P <- function(...) file.path(OUT_DIR, ...)
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
theme_set(theme_bw(base_size = BASE_SIZE) + theme(panel.grid.minor = element_blank()))
mod_label <- function(m) ifelse(!is.na(MODULE_NAMES[m]) & nzchar(MODULE_NAMES[m]), paste0(m, " (", MODULE_NAMES[m], ")"), m)
CMP <- paste0(CASE, " vs ", CONTROL)
made <- character(0); skipped <- character(0)
save_gg <- function(name, p, w, h) {
  f <- file.path(FIG, paste0(ACC, "_", name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = FIG_DPI, bg = "white"); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h, bg = "white")
  made <<- c(made, name)
}
save_base <- function(name, w, h, draw) {
  f <- file.path(FIG, paste0(ACC, "_", name, ".png"))
  png(f, width = w, height = h, units = "in", res = FIG_DPI); draw(); dev.off()
  pdf(sub("\\.png$", ".pdf", f), width = w, height = h); draw(); dev.off()
  made <<- c(made, name)
}
fig <- function(name, expr) tryCatch(expr, error = function(e) { skipped <<- c(skipped, paste0(name, ": ", conditionMessage(e))) })
need <- function(x, what) if (is.null(x) || !NROW(x)) stop("missing input (", what, ")")
cat(">>> figures for", run_id, "->", FIG, "\n")
## ======================= B. DEG =======================
fig("B_DEG_counts_bar", {
  fs <- list.files(file.path(in_dir(BASE_DIR, "06_deg_results")), pattern = "_DEG\\.csv$", full.names = TRUE)
  fs <- fs[!grepl("interaction", fs)]; need(fs, "06_deg_results")
  d <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f); s <- !is.na(x$padj) & x$padj < 0.05 & abs(x$log2FoldChange) >= 1
    data.frame(Comparison = sub("_DEG\\.csv$", "", sub(paste0("^", ACC, "_"), "", basename(f))),
               Direction = c("Up", "Down"), N = c(sum(s & x$log2FoldChange > 0), sum(s & x$log2FoldChange < 0))) }))
  d$Comparison <- gsub("_vs_", " vs ", d$Comparison)
  save_gg("B_DEG_counts_bar", ggplot(d, aes(Comparison, N, fill = Direction)) + geom_col(position = "dodge") +
    geom_text(aes(label = N), position = position_dodge(0.9), vjust = -0.3, size = 2.6) +
    scale_fill_manual(values = c(Up = "#b2182b", Down = "#2166ac")) + coord_cartesian(clip = "off") +
    labs(title = "Significant DEGs (padj < 0.05, |log2FC| >= 1)", x = NULL, y = "Genes") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)), 5.5, 4)
})
## ======================= A. QC: sample clustering / outliers =======================
fig("A_sample_clustering_outliers", {
  base <- dirname(in_dir(BASE_DIR, "02_metadata"))
  vf <- list.files(file.path(base, "04_normalized"), pattern = "_vst_normalized_matrix\\.csv$", full.names = TRUE)[1]
  mf <- list.files(file.path(base, "02_metadata"), pattern = "_sample_metadata\\.csv$", full.names = TRUE)[1]
  if (is.na(vf) || is.na(mf)) stop("missing input (vst matrix / metadata)")
  x <- as.matrix(read.csv(vf, row.names = 1, check.names = FALSE)); md <- read.csv(mf, stringsAsFactors = FALSE)
  grp <- setNames(md$group, md$sample_id)[colnames(x)]
  x <- x[order(-apply(x, 1, stats::var))[seq_len(min(2000, nrow(x)))], ]          # top-variance genes
  cc <- stats::cor(x); hc <- stats::hclust(stats::as.dist(1 - cc), "average")
  mc <- (colSums(cc) - 1) / (ncol(cc) - 1); thr <- stats::quantile(mc, .25) - 1.5 * stats::IQR(mc)
  out <- names(mc)[mc < thr]
  write.csv(data.frame(Sample = names(mc), Group = grp, Mean_correlation = round(mc, 4), Outlier = mc < thr),
            file.path(FIG, paste0(ACC, "_A_sample_outliers.csv")), row.names = FALSE)
  cols <- setNames(grDevices::hcl.colors(length(unique(grp)), "Dark 3"), unique(grp))
  save_base("A_sample_clustering_outliers", max(7, ncol(x) * 0.07), 5, function() {
    op <- graphics::par(mar = c(6, 4, 3, 1)); on.exit(graphics::par(op))
    dd <- stats::as.dendrogram(hc)
    dd <- stats::dendrapply(dd, function(n) { if (is.leaf(n)) { l <- attr(n, "label")
      attr(n, "nodePar") <- list(pch = NA, lab.col = if (l %in% out) "red" else cols[grp[l]], lab.cex = 0.45) }; n })
    plot(dd, main = paste0("Sample clustering (1 - Pearson, top-variance genes); red = outlier (", length(out), ")"), ylab = "1 - r")
    graphics::legend("topright", legend = names(cols), text.col = cols, bty = "n", cex = 0.7, inset = c(0, 0.08))
  })
})
## ======================= C. WGCNA =======================
ws <- if (file.exists(P("10_wgcna_results", paste0(ACC, "_WGCNA_workspace.rds")))) readRDS(P("10_wgcna_results", paste0(ACC, "_WGCNA_workspace.rds"))) else NULL
fig("C_eigengene_clustering", {
  need(ws, "WGCNA workspace"); me <- ws$MEs; me <- me[, colnames(me) != "MEgrey", drop = FALSE]
  hc <- stats::hclust(stats::as.dist(1 - stats::cor(me)), "average")
  save_base("C_eigengene_clustering", 7, 4.5, function() {
    plot(hc, main = "Module eigengene clustering (merge height 0.25 marked)", xlab = "", sub = "", ylab = "1 - r", cex = 0.7)
    graphics::abline(h = 0.25, col = "red", lty = 2) })
})
mmgs <- rd(P("10_wgcna_results", paste0(ACC, "_TBmodules_MM_GS.csv")))
fig("C_MM_vs_GS", {
  need(mmgs, "TBmodules_MM_GS"); mmgs$ModuleLab <- mod_label(mmgs$Module)
  r <- do.call(rbind, lapply(split(mmgs, mmgs$ModuleLab), function(d) { ct <- stats::cor.test(abs(d$MM), abs(d$GS_TB))
    data.frame(ModuleLab = d$ModuleLab[1], lab = sprintf("r = %.2f, p = %.1e", ct$estimate, ct$p.value)) }))
  save_gg("C_MM_vs_GS", ggplot(mmgs, aes(abs(MM), abs(GS_TB))) + geom_point(aes(colour = Module), size = 0.5, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.3, size = 2.6) +
    scale_colour_identity() + facet_wrap(~ ModuleLab) +
    labs(title = "Module membership vs gene significance", x = "|MM|", y = paste0("|GS| (", CMP, ")")), 7.5, 3.2)
})
## ======================= D. enrichment =======================
enr <- do.call(rbind, lapply(c("GO_BP", "KEGG", "Reactome"), function(db) do.call(rbind, lapply(names(MODULE_NAMES), function(m) {
  d <- rd(P("11_module_enrichment", paste0(ACC, "_", m, "_", db, ".csv")))
  if (is.null(d) || !nrow(d) || !"p.adjust" %in% names(d)) return(NULL)
  d <- d[d$p.adjust < 0.05, ]; if (!nrow(d)) return(NULL)
  gr <- vapply(strsplit(d$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 0)
  data.frame(DB = db, Module = m, Description = d$Description, GeneRatio = gr, Count = d$Count, padj = d$p.adjust)
}))))
for (db in c("GO_BP", "KEGG", "Reactome")) fig(paste0("D_dotplot_", db), {
  d <- enr[!is.null(enr) & enr$DB == db, , drop = FALSE]; need(d, paste(db, "terms with padj < 0.05"))
  d <- do.call(rbind, lapply(split(d, d$Module), function(z) utils::head(z[order(z$padj), ], TOP_TERMS)))
  d$Description <- substr(d$Description, 1, 60); d$ModuleLab <- mod_label(d$Module)
  d$Term <- factor(paste(d$Description), levels = unique(d$Description[order(d$Module, -d$padj)]))
  save_gg(paste0("D_dotplot_", db), ggplot(d, aes(GeneRatio, Term, size = Count, colour = -log10(padj))) + geom_point() +
    facet_grid(ModuleLab ~ ., scales = "free_y", space = "free_y") + scale_colour_viridis_c() +
    labs(title = paste0(db, " enrichment of the selected modules"), y = NULL) +
    theme(strip.text.y = element_text(angle = 0)), 7.5, max(3, 0.18 * nrow(d) + 1.2))
})
fig("D_module_pathway_heatmap", {
  need(enr, "enrichment terms"); top <- unique(unlist(lapply(split(enr, enr$Module), function(z) utils::head(z$Description[order(z$padj)], 6))))
  d <- expand.grid(Description = top, Module = names(MODULE_NAMES), stringsAsFactors = FALSE)
  best <- stats::aggregate(padj ~ Description + Module, enr[enr$Description %in% top, ], min)
  d <- merge(d, best, all.x = TRUE); d$score <- ifelse(is.na(d$padj), 0, -log10(d$padj)); d$ModuleLab <- mod_label(d$Module)
  d$Description <- substr(d$Description, 1, 60)
  save_gg("D_module_pathway_heatmap", ggplot(d, aes(ModuleLab, Description, fill = score)) + geom_tile(colour = "white") +
    scale_fill_gradient(low = "white", high = "#b2182b", name = "-log10 padj") +
    labs(title = "Top enriched terms across modules (GO BP / KEGG / Reactome)", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)), 7.5, max(3, 0.2 * length(top) + 1.5))
})
## ======================= E. network =======================
cen <- rd(P("12_string_ppi", paste0(ACC, "_COMBINED_centrality.csv")))
gs  <- rd(P("15_network_algorithms", paste0(ACC, "_GeneScores_RWR_Diffusion.csv")))
gp  <- rd(P("19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")))
topbar <- function(name, d, val, title, xlab) fig(name, {
  need(d, name); d <- utils::head(d[order(-d[[val]]), ], TOP_N); d$Gene <- factor(d$Gene, levels = rev(d$Gene))
  save_gg(name, ggplot(d, aes(.data[[val]], Gene, fill = Module)) + geom_col() + scale_fill_identity(guide = "legend",
    labels = function(x) mod_label(x)) + labs(title = title, x = xlab, y = NULL, fill = "Module"), 5.5, 4.8)
})
topbar("E_top_PPI_hubs", cen, "Degree", "Top PPI hubs (combined network)", "Degree")
topbar("E_top_RWR_genes", gs, "RWR_score", "Top genes by random walk with restart", "RWR score")
topbar("E_top_diffusion_genes", gs, "Diffusion_score", "Top genes by network diffusion", "Diffusion score")
net_edges <- rd(P("12_string_ppi", paste0(ACC, "_COMBINED_STRING_network.csv")))
comm <- rd(P("15_network_algorithms", paste0(ACC, "_Communities.csv")))
fig("E_network_labelled", {
  need(net_edges, "combined network"); need(cen, "centrality")
  g <- igraph::graph_from_data_frame(unique(net_edges[, c("preferredName_A", "preferredName_B")]), directed = FALSE)
  g <- igraph::simplify(g)
  cmp <- igraph::components(g); g <- igraph::induced_subgraph(g, which(cmp$membership == which.max(cmp$csize)))  # small detached pieces squash the layout
  set.seed(42); lay <- igraph::layout_with_fr(g)
  v <- data.frame(Gene = igraph::V(g)$name, x = lay[, 1], y = lay[, 2])
  v <- merge(v, cen[, c("Gene", "Module")], all.x = TRUE)
  if (!is.null(gs)) v <- merge(v, gs[, c("Gene", "RWR_score")], all.x = TRUE) else v$RWR_score <- 1
  hi <- if (!is.null(gp)) gp$Gene[gp$Final_Priority == "High"] else character(0)
  el <- igraph::as_edgelist(g); e <- data.frame(x = lay[match(el[, 1], igraph::V(g)$name), 1], y = lay[match(el[, 1], igraph::V(g)$name), 2],
    xend = lay[match(el[, 2], igraph::V(g)$name), 1], yend = lay[match(el[, 2], igraph::V(g)$name), 2])
  lab <- v[v$Gene %in% hi, ]
  p <- ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey85", linewidth = 0.1) +
    geom_point(data = v, aes(x, y, size = RWR_score, colour = Module), alpha = 0.8) + scale_colour_identity(guide = "legend", labels = function(x) mod_label(x)) +
    scale_size_continuous(range = c(0.3, 3.5), name = "RWR score") +
    labs(title = "Combined disease PPI network", subtitle = "node size = RWR score; labels = High-priority genes", colour = "Module") +
    theme_void(base_size = BASE_SIZE)
  if (nrow(lab)) p <- p + (if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = lab, aes(x, y, label = Gene), size = 2, max.overlaps = 40, segment.size = 0.2) else
    geom_text(data = lab, aes(x, y, label = Gene), size = 2, vjust = -0.8))
  save_gg("E_network_labelled", p, 8, 7)
  ## communities: same layout, coloured by community (largest 8)
  if (!is.null(comm)) {
    v2 <- merge(v, comm[, c("Gene", "Community")], all.x = TRUE)
    big <- names(sort(table(v2$Community), decreasing = TRUE))[1:min(8, length(unique(v2$Community)))]
    v2$Comm <- ifelse(v2$Community %in% big, paste("C", v2$Community), "other")
    save_gg("E_network_communities", ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey88", linewidth = 0.1) +
      geom_point(data = v2, aes(x, y, colour = Comm), size = 0.7) + labs(title = "Network communities (largest 8)", colour = "Community") +
      theme_void(base_size = BASE_SIZE), 7.5, 6.5)
  }
})
fig("E_consensus_genes", {
  need(gp, "consensus genes"); d <- gp[gp$Final_Priority == "High", ]; need(d, "High genes")
  crit <- grep("^Top_", names(d), value = TRUE)
  L <- do.call(rbind, lapply(crit, function(cname) data.frame(Gene = d$Gene, Module = d$Module, Criterion = sub("^Top_", "", cname), Met = d[[cname]])))
  L$Gene <- factor(L$Gene, levels = d$Gene[order(d$Module, -d$N_criteria)])
  save_gg("E_consensus_genes", ggplot(L, aes(Criterion, Gene, fill = ifelse(Met %in% TRUE, Module, "white"))) + geom_tile(colour = "grey80") +
    scale_fill_identity() + labs(title = "High-priority genes: criteria met (filled = top 20%)", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 30, hjust = 1)), 4.5, max(4, 0.13 * nrow(d) + 1.2))
})
## ======================= F. GSEA =======================
fig("F_GSEA_dotplot", {
  d <- do.call(rbind, lapply(c("HALLMARK", "REACTOME", "CELLTYPE"), function(cl) {
    x <- rd(file.path(BASE_DIR, "17_gsea", paste0(ACC, "_", DEG_NAME, "_GSEA_", cl, ".csv")))
    if (is.null(x)) return(NULL); x <- x[x$p.adjust < 0.05, ]; if (!nrow(x)) return(NULL)
    x <- rbind(utils::head(x[order(-x$NES), ], 8), utils::head(x[order(x$NES), ], 8)); x <- x[!duplicated(x$ID), ]
    data.frame(Collection = cl, Term = substr(sub("^(HALLMARK|REACTOME)_", "", x$Description), 1, 55), NES = x$NES, padj = x$p.adjust, Size = x$setSize) }))
  need(d, "17_gsea tables"); d$Term <- factor(d$Term, levels = unique(d$Term[order(d$NES)]))
  save_gg("F_GSEA_dotplot", ggplot(d, aes(NES, Term, size = Size, colour = -log10(padj))) + geom_point() + geom_vline(xintercept = 0, linetype = 2) +
    facet_grid(Collection ~ ., scales = "free_y", space = "free_y") + scale_colour_viridis_c() +
    labs(title = paste0("GSEA, ", CMP, " (padj < 0.05)"), y = NULL) + theme(strip.text.y = element_text(angle = 0)), 7.5, max(4, 0.15 * nrow(d) + 1.5))
})
## ======================= G. drugs =======================
## use the same cleaned target lists as the consensus (step 15 drops screening-only and ADME targets)
clean_dt <- function(dt, dp) {
  if (is.null(dt) || is.null(dp) || !"Known_Targets_used" %in% names(dp)) return(dt)
  used <- strsplit(ifelse(is.na(dp$Known_Targets_used), "", dp$Known_Targets_used), ";"); names(used) <- dp$Drug
  dt[mapply(function(dr, tg) tg %in% used[[dr]], dt$Drug, dt$Drug_Target), , drop = FALSE]
}
dp <- rd(P("19_consensus", paste0(ACC, "_DrugPriority.csv")))
fig("G_iLINCS_top_reversal", {
  x <- rd(P("13_ilincs_signature", paste0(ACC, "_iLINCS_candidate_compounds.csv"))); need(x, "iLINCS candidates")
  x <- x[!grepl("^BRD-", x$Target), ]; x <- stats::aggregate(Similarity ~ Target, x, min)
  x <- utils::head(x[order(x$Similarity), ], TOP_N); x$Target <- factor(x$Target, levels = rev(x$Target))
  save_gg("G_iLINCS_top_reversal", ggplot(x, aes(Similarity, Target)) + geom_col(fill = "#2166ac") +
    labs(title = "Strongest reversal of the disease signature (iLINCS, named compounds)", x = "Concordance (more negative = stronger reversal)", y = NULL), 6, 4.8)
})
fig("G_network_proximity", {
  x <- rd(P("15_network_algorithms", paste0(ACC, "_NetworkProximity.csv"))); need(x, "network proximity")
  x$lab <- ifelse(rank(x$P_empirical, ties.method = "first") <= 10, x$Drug, "")
  p <- ggplot(x, aes(Z_score, -log10(P_empirical), size = N_targets_in_network)) + geom_point(alpha = 0.5, colour = "grey30") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") + geom_vline(xintercept = 0, linetype = 3) +
    labs(title = "Network proximity of drug targets to the disease network", subtitle = "red line = empirical p 0.05 (no FDR correction)",
         x = "Proximity Z (negative = closer than random)", y = "-log10 empirical p", size = "Targets in network")
  if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(aes(label = lab), size = 2.2, max.overlaps = 30)
  save_gg("G_network_proximity", p, 6.5, 5)
})
dt <- clean_dt(rd(P("14_drug_targets", paste0(ACC, "_Drug_Target_TBnetwork.csv"))), dp)
fig("G_drug_module_heatmap", {
  need(dp, "drug priority"); need(dt, "drug targets")
  keep <- dp$Drug[dp$Final_Priority %in% c("High", "Medium")]; d <- dt[dt$Drug %in% keep, ]; need(d, "High/Medium drug targets")
  d <- stats::aggregate(Drug_Target ~ Drug + TB_Module, d, function(z) length(unique(z)))
  d$TB_Module <- mod_label(d$TB_Module); d$Drug <- factor(d$Drug, levels = rev(intersect(dp$Drug[order(-dp$N_criteria)], d$Drug)))
  save_gg("G_drug_module_heatmap", ggplot(d, aes(TB_Module, Drug, fill = Drug_Target)) + geom_tile(colour = "white") + geom_text(aes(label = Drug_Target), size = 2) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15", name = "Targets") + labs(title = "High/Medium drugs: targets per module", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = 5.5), axis.text.x = element_text(angle = 20, hjust = 1)), 5.5, max(4, 0.1 * length(unique(d$Drug)) + 1.5))
})
fig("G_drug_target_network", {
  need(dp, "drug priority"); need(dt, "drug targets")
  hd <- dp$Drug[dp$Final_Priority == "High"]; d <- unique(dt[dt$Drug %in% hd, c("Drug", "Drug_Target", "TB_Module")]); need(d, "High drug targets")
  g <- igraph::graph_from_data_frame(d[, 1:2], directed = FALSE); set.seed(42); lay <- igraph::layout_with_fr(g)
  v <- data.frame(name = igraph::V(g)$name, x = lay[, 1], y = lay[, 2]); v$type <- ifelse(v$name %in% d$Drug, "Drug", "Target")
  v$col <- ifelse(v$type == "Drug", "grey20", d$TB_Module[match(v$name, d$Drug_Target)])
  el <- igraph::as_edgelist(g); e <- data.frame(x = v$x[match(el[, 1], v$name)], y = v$y[match(el[, 1], v$name)], xend = v$x[match(el[, 2], v$name)], yend = v$y[match(el[, 2], v$name)])
  save_gg("G_drug_target_network", ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey70", linewidth = 0.3) +
    geom_point(data = v, aes(x, y, shape = type, colour = col), size = 3) + scale_colour_identity() + scale_shape_manual(values = c(Drug = 15, Target = 16)) +
    geom_text(data = v, aes(x, y, label = name), size = 2, vjust = -1.1) + labs(title = "High-priority drugs and their disease-network targets", subtitle = "target colour = module", shape = NULL) +
    theme_void(base_size = BASE_SIZE), 7.5, 6.5)
})
fig("G_consensus_drugs", {
  need(dp, "drug priority"); d <- dp[dp$Final_Priority %in% c("High", "Medium"), ]; need(d, "High/Medium drugs")
  crit <- c("Strong_reversal", "Negative_prox_Z", "Significant_prox_P", "Targets_in_TBnet", "Module_enriched", "Hits_high_gene", "Top_target_RWR")
  crit <- intersect(crit, names(d)); d <- utils::head(d[order(-d$N_criteria, d$iLINCS_Reversal_Score), ], 40)
  L <- do.call(rbind, lapply(crit, function(cn) data.frame(Drug = d$Drug, Criterion = cn, Met = d[[cn]] %in% TRUE, Tier = d$Final_Priority)))
  L$Drug <- factor(L$Drug, levels = rev(d$Drug))
  save_gg("G_consensus_drugs", ggplot(L, aes(Criterion, Drug, fill = ifelse(Met, Tier, "none"))) + geom_tile(colour = "grey80") +
    scale_fill_manual(values = c(High = "#a50f15", Medium = "#fc9272", none = "white"), name = "Tier (met)") +
    labs(title = "Drug consensus: criteria met (top 40 High/Medium)", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6)), 5.5, max(4, 0.13 * nrow(d) + 1.5))
})
cat("\n>>> made", length(made), "figures (PNG 600 dpi + PDF):\n   ", paste(made, collapse = "\n    "), "\n")
if (length(skipped)) cat("\n>>> skipped", length(skipped), ":\n   ", paste(skipped, collapse = "\n    "), "\n")
cat(">>> DONE. Folder:", FIG, "\n")
```

### scripts/wgcna_network/GSE114192_interaction_gsea.R

```r
ACC        <- "GSE114192"
RUN        <- "TB_DM"
COEFS      <- c(interaction = "interaction", main_TB = "main_TB", main_DM = "main_DM")  # file suffixes
MIN_GS <- 15; MAX_GS <- 500; PADJ_SIG <- 0.05; SEED <- 42
COLLECTIONS <- list(HALLMARK = list(collection = "H"),
                    REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"))
```

### scripts/wgcna_network/GSE165709_atac_geneset_step18.R

```r
COHORTS    <- c("HC", "PLWH")
N_PERM     <- 2000
MIN_GENES  <- 10        # a set needs at least this many genes with a peak
TSS_ONLY   <- FALSE     # TRUE = use only promoter peaks (cleaner link, far fewer peaks)
```

### scripts/wgcna_network/GSE165709_differential_step16.R

```r
ASSAY        <- "ATAC"             # "RNA" or "ATAC"
COUNT_FILE   <- list(RNA  = "GSE165708_non-normalized_estimated_counts_matrix.txt",
                     ATAC = "GSE165703_ATAC_FeatureCount_raw_quantification_matrix.txt")
KEEP_COHORTS <- c("HC", "PLWH")   # PrEP (treatment) is left out of the main analysis
CASE_LEVEL   <- "INF"; CTRL_LEVEL <- "NEG"
MIN_COUNT    <- 10
PADJ_CUT     <- 0.05
TSS_REGION   <- c(-3000, 3000)
req <- c("DESeq2")
if (ASSAY == "ATAC") req <- c(req, "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomicRanges", "IRanges")
ok <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing package(s): ", paste(req[!ok], collapse = ", "),
                   "\n  BiocManager::install(c(", paste0('"', req[!ok], '"', collapse = ", "), "))")
```

### scripts/wgcna_network/GSE165709_integration_step17.R

```r
COHORTS   <- c("HC", "PLWH")       # first = reference (normal response)
PADJ_CUT  <- 0.05
LFC_CUT   <- 0                     # raise to e.g. 0.5 for a stricter definition of "responds"
```

### scripts/wgcna_network/ilincs_proliferation_sensitivity.R

```r
ILINCS_LIBRARY <- "CP"
TOP_N   <- c(50, 100, 200)
N_PERM  <- 2000
SEED    <- 42
## proliferation genes: Hallmark G2M checkpoint, E2F targets and mitotic spindle (msigdbr), plus CDKN1B
PROLIF_SETS <- c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_MITOTIC_SPINDLE")
for (p in c("drugfindR", "msigdbr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })
sig_dir   <- file.path(ROOT, "GSE161829", "TBneg_vs_ATB", "13_ilincs_signature")
pilot_dir <- file.path(ROOT, "GSE114192", "13_ilincs_signature")
out_dir   <- file.path(sig_dir, "sensitivity_no_proliferation"); dir.create(out_dir, showWarnings = FALSE)
```

### scripts/wgcna_network/metabolomics_bridge_step20.R

```r
STUDY_ID     <- "ST000974"
MET_PADJ     <- 0.05     # "significant" metabolites
MET_P_BROAD  <- 0.01     # broader set (nominal), reported alongside -- there are few FDR hits
GENE_PADJ    <- 0.05; GENE_LFC <- 1
KEGG_SLEEP   <- 0.34     # KEGG asks for <= 3 requests per second
for (p in c("httr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
```

### scripts/wgcna_network/metabolomics_step19.R

```r
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
```

### scripts/wgcna_network/ml_validation_step22.R

```r
TIERS     <- c("High")              # genes used as features; c("High","Medium") for a larger panel
TESTS <- list(                       # dataset, group A (control), group B (case), label
  list(ds = "GSE161829", a = "TBneg",   b = "ATB",    lab = "GSE161829 ATB vs TBneg"),
  list(ds = "GSE161829", a = "LTBI",    b = "ATB",    lab = "GSE161829 ATB vs LTBI"),
  list(ds = "GSE161829", a = "TBneg",   b = "LTBI",   lab = "GSE161829 LTBI vs TBneg"),
  list(ds = "GSE222001", a = "Healthy", b = "Active", lab = "GSE222001 Active vs Healthy (n=5+5)"),
  list(ds = "GSE222001", a = "Latent",  b = "Active", lab = "GSE222001 Active vs Latent (n=5+5)"))
SEED <- 42; N_BOOT <- 1000; N_RANDOM <- 1000; ALPHA <- 0.5
N_OUTER <- 5; N_REPEATS <- 10        # nested CV: 5 outer folds, repeated 10 times
for (p in c("glmnet", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(ggplot2))
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR <- file.path(ROOT_DIR, ACC); OUT_DIR <- run_out_dir(BASE_DIR)
dir_out  <- file.path(OUT_DIR, "24_ml_validation"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
FIG_DPI <- 600
save_gg <- function(f, p, w, h) { ggsave(f, p, width = w, height = h, dpi = FIG_DPI); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h) }
```

### scripts/wgcna_network/WGCNA_combined_network_step5.R

```r
SPECIES         <- 9606
REQUIRED_SCORE  <- 400
USE_FULL_MODULE <- TRUE
N_LABELS_FULL   <- 20
```

### scripts/wgcna_network/WGCNA_consensus_priority_step15.R

```r
TOP_FRACTION   <- 0.20    # "top-ranked" = best 20% of a criterion
MERGE_RWR_DIFF <- TRUE    # RWR and diffusion are near-duplicates (correlation printed below): count them as ONE vote
HIGH_N         <- 4       # consensus tiers (number of criteria met)
MED_N          <- 2
PROX_Z_CUT     <- -1      # drug: targets closer to the disease genes than average
MIN_TARGETS    <- 3       # drug: targets inside the TB network
P_CUT          <- 0.05    # drug: proximity and module-enrichment p-values
## drug-target cleaning before counting targets (set both to c()/FALSE to reproduce the original pilot numbers)
EXCLUDE_SCREEN_TARGETS <- c("NFE2L2")   # linked to >100 compounds by screening-assay records, not mechanism
EXCLUDE_ADME_TARGETS   <- TRUE          # drug-metabolism / transport genes (CYP, UGT, SULT, SLC, ABC)
ADME_PATTERN <- "^(CYP[0-9]|UGT[0-9]|SULT[0-9]|SLC[0-9]|ABC[A-G][0-9])"
```

### scripts/wgcna_network/WGCNA_drug_targets_step8.R

```r
DGIDB_URL <- "https://dgidb.org/api/graphql"
```

### scripts/wgcna_network/WGCNA_Enrichr_drugs_step12.R

```r
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
```

### scripts/wgcna_network/WGCNA_GSEA_step11.R

```r
MIN_GS <- 15; MAX_GS <- 500     # gene-set size limits (genes present in the ranked list)
PADJ_SIG <- 0.05
COLLECTIONS <- list(            # name -> arguments for msigdbr()
  HALLMARK = list(collection = "H"),
  REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"),
  GOBP     = list(collection = "C5", subcollection = "GO:BP"),
  CELLTYPE = list(collection = "C8")
)
```

### scripts/wgcna_network/WGCNA_iLINCS_query_step7.R

```r
ILINCS_LIBRARY       <- "CP"     # Chemical Perturbagen library (drug repurposing)
SIMILARITY_CUTOFF    <- 0.321    # drugfindR's own default for consensusConcordants()
REVERSAL_ONLY        <- TRUE
MIN_GENES_DIR        <- 5        # a direction with fewer genes is not queried (too small for iLINCS)     # keep only NEGATIVE similarity (true reversal), drop mimetic (positive) hits
```

### scripts/wgcna_network/WGCNA_iLINCS_signature_step6.R

```r
LOGFC_THRESHOLD <- 0.5   # minimum |log2FC| for TB_UP / TB_DOWN (applied after prepareSignature)
```

### scripts/wgcna_network/WGCNA_network_algorithms_step9.R

```r
RESTART_R    <- 0.5     # RWR restart probability (0.3-0.5 is standard)
DIFFUSE_A    <- 0.5     # diffusion keep-vs-spread balance
MAX_ITER     <- 200
TOLERANCE    <- 1e-10
N_PERMUTATIONS <- 1000  # for the network-proximity null distribution
```

### scripts/wgcna_network/WGCNA_singlecell_step13.R

```r
SC_SUBDIR       <- file.path("external_data", "SCP1749")
CELLTYPE_COL    <- "CellTypeAnnotations"
DETECT_MIN_CP10K <- 1
SPECIFIC_SHARE  <- 0.5
MIN_TARGETS_SC  <- 3        # a drug needs at least this many detected targets to get a cell-type profile
MIN_GENES_TEST  <- 10       # min genes per group for the module-level tests
IMMUNE_TYPES    <- c("Macrophage", "T", "Neutrophil", "Mast", "B", "pDC", "Plasma")   # for the blood-relevant drug view
DIFFUSE_FACTOR  <- 1.5      # a top cell type must hold >= this x the uniform share, otherwise the drug is labelled "diffuse"
```

### scripts/wgcna_network/WGCNA_STRING_step4.R

```r
SPECIES         <- 9606     # human
REQUIRED_SCORE  <- 400      # STRING "medium confidence" default (0-1000); use 700+ for high-confidence only
USE_FULL_MODULE <- TRUE     # TRUE  = every gene in the module (default)
                             # FALSE = only the step-3 "disease genes" (DEG + module overlap) subset
N_LABELS_FULL   <- 15        # how many top hubs to label on the full-network plot (rest stay unlabeled dots)
N_HUBS_SUBNET   <- 30        # size of the separate, fully-labeled "hub-only" zoomed-in plot
```

### scripts/wgcna_network/WGCNA_TBmodules_step2_enrichment.R

```r
PADJ_CUTOFF <- 0.05   # for the "top term" summary only; full tables are saved unfiltered
```

