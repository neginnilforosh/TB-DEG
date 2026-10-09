## Configuration shared by all wgcna_network scripts.
## Pick the run here (or set the environment variable TB_RUN, e.g. TB_RUN=GSE161829_ATB Rscript step1.R).

ACTIVE_RUN <- "GSE114192_TB"
##ACTIVE_RUN <- "GSE161829_ATB"
PRESETS <- list(
  GSE114192_TB = list(                       # pilot: Healthy_Control vs TB_Only
    ACC = "GSE114192", CONTROL = "Healthy_Control", CASE = "TB_Only",
    DEG_NAME = "HealthyControl_vs_TBOnly",     # <ACC>_<DEG_NAME>_DEG.csv in 06_deg_results
    RUN = "",                                  # "" = outputs directly in the dataset folder
    MODULE_NAMES = c(green = "Interferon_Response", purple = "Ribosome_Biogenesis_Chromatin", blue = NA),
    MODULE_PLOT_COLORS = c(green = "#2ecc71", purple = "#9b59b6", blue = "#3498db")),
  GSE161829_ATB = list(                      # replication: TBneg vs active TB (LTBI samples excluded)
    ACC = "GSE161829", CONTROL = "TBneg", CASE = "ATB",
    DEG_NAME = "TBneg_vs_ATB",
    RUN = "TBneg_vs_ATB",                      # outputs in GSE161829/TBneg_vs_ATB/
    MODULE_NAMES = character(0),               # fill in after step 2
    MODULE_PLOT_COLORS = NULL)                 # NULL = use the WGCNA colour names
)

## ---- apply the active preset ----
run_id <- Sys.getenv("TB_RUN", ACTIVE_RUN)
if (!run_id %in% names(PRESETS)) stop("Unknown run '", run_id, "'. Available: ", paste(names(PRESETS), collapse = ", "))
cfg <- PRESETS[[run_id]]
ACC <- cfg$ACC; CONTROL <- cfg$CONTROL; CASE <- cfg$CASE; TRAIT_COL <- CASE
DEG_NAME <- cfg$DEG_NAME; RUN <- cfg$RUN
MODULE_NAMES <- cfg$MODULE_NAMES; MODULE_PLOT_COLORS <- cfg$MODULE_PLOT_COLORS
DEG_FILE_NAME <- paste0(ACC, "_", DEG_NAME, "_DEG.csv")
cat(">>> Run:", run_id, "|", ACC, ":", CONTROL, "vs", CASE, if (nzchar(RUN)) paste0("| outputs in ", ACC, "/", RUN) else "", "\n")

## output folder of the run (inputs 02/06/07/09 always stay in the dataset folder)
run_out_dir <- function(base_dir) if (nzchar(RUN)) file.path(base_dir, RUN) else base_dir

## "ENSG00000138496.16|PARP9" -> "ENSG00000138496"; plain Ensembl IDs are left unchanged
clean_gene_ids <- function(x) sub("\\.[0-9]+$", "", sub("\\|.*$", "", x))

