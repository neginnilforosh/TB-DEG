# Reproducibility record

Generated 2026-10-14 01:55 with R version 4.5.0 (2025-04-11) on macOS 27.0.1.

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
Enrichr, DGIdb, KEGG REST, Metabolomics Workbench. See the run dates in each output folder.

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
    MODULE_NAMES = character(0),               # fill in after step 2
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
| scripts/wgcna_network/GSE114192_interaction_gsea.R | 55 | `set.seed(SEED)` |
| scripts/wgcna_network/GSE165709_atac_geneset_step18.R | 70 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_combined_network_step5.R | 172 | `set.seed(42)` |
| scripts/wgcna_network/WGCNA_GSEA_step11.R | 92 | `set.seed(42)` |
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
COUNT_FILE  <- "GSE222001_count.txt" # <-- from A2
## Map the raw labels GEO uses onto clean names.
GROUP_MAP <- c(
  "Healthy control"                = "Healthy",
  "Active tuberculosis infection"  = "Active",
  "Latent tuberculosis infection"  = "Latent"
)
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

