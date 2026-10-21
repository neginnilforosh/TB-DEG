## Configuration shared by all wgcna_network scripts.
## Pick the run here (or set the environment variable TB_RUN, e.g. TB_RUN=GSE161829_ATB Rscript step1.R).

##ACTIVE_RUN <- "GSE114192_TBDM_vs_DM"
ACTIVE_RUN <- "GSE114192_TB"
##ACTIVE_RUN <- "GSE161829_ATB"
PRESETS <- list(
  GSE114192_TB = list( 
    ACC = "GSE114192",
    CONTROL = "Healthy_Control",
    CASE = "TB_Only",
    DEG_NAME = "HealthyControl_vs_TBOnly",
    RUN = "",
    MODULE_NAMES = c(green = "Interferon_Response", purple = 							"Ribosome_Biogenesis_Chromatin", blue = 											"Mixed_CellComposition"),
    MODULE_PLOT_COLORS = c(green = "#2ecc71", purple = "#9b59b6"
    								, blue = "#3498db")),
  
  GSE161829_ATB = list(
    ACC = "GSE161829",
    CONTROL = "TBneg",
    CASE = "ATB",
    DEG_NAME = "TBneg_vs_ATB",
    RUN = "TBneg_vs_ATB",                
    MODULE_NAMES = c(midnightblue = "M1_Cell_Cycle_Proliferation",
                greenyellow = "M2_Mixed_with_IFN_cluster",
                 salmon       = "M3_OxPhos_weak"),        
    MODULE_PLOT_COLORS = NULL) ,               
  
  GSE114192_TBDM_vs_DM = list(
    ACC = "GSE114192", CONTROL = "DM_only", CASE = "TB_DM",
    DEG_NAME = "DM_only_vs_TB_DM",
    RUN = "TB_DM",
    MODULE_NAMES = c(black = "M1_Interferon_Complement",
                 greenyellow = "M2_TLR_NFkB_Inflammation",
                 green = "M3_Mixed_CellComposition"),
    MODULE_PLOT_COLORS = NULL),
	
  GSE114192_TBDM_vs_TB = list(
	 ACC = "GSE114192",
	 CONTROL = "TB_Only",
	 CASE = "TB_DM",
    DEG_NAME = "TB_Only_vs_TB_DM", 
    RUN = "TB_DM",
    MODULE_NAMES = character(0), 
    MODULE_PLOT_COLORS = NULL),
    
  GSE37250_TBHIV = list(
   ACC = "GSE37250",
   CONTROL = "LTBI_HIVpos",
   CASE = "ATB_HIVpos",
   DEG_NAME = "LTBI_HIVpos_vs_ATB_HIVpos",
   RUN = "TB_HIV",
   MODULE_NAMES = c(yellow = "M1_Inflammatory_Neutrophil",
                 turquoise = "M2_Ribosome_Biogenesis_MYC",
                 purple = "M3_Interferon_Response"),
   MODULE_PLOT_COLORS = NULL)

)
## ---- apply the active preset ----
run_id <- trimws(Sys.getenv("TB_RUN", ACTIVE_RUN))
if (!run_id %in% names(PRESETS)) stop("Unknown run '", run_id, "'. Available: ", paste(names(PRESETS), collapse = ", "))
cfg <- PRESETS[[run_id]]
ACC <- cfg$ACC; CONTROL <- cfg$CONTROL; CASE <- cfg$CASE; TRAIT_COL <- CASE
DEG_NAME <- cfg$DEG_NAME; RUN <- cfg$RUN
MODULE_NAMES <- cfg$MODULE_NAMES; MODULE_PLOT_COLORS <- cfg$MODULE_PLOT_COLORS
DEG_FILE_NAME <- paste0(ACC, "_", DEG_NAME, "_DEG.csv")
cat(">>> Run:", run_id, "|", ACC, ":", CONTROL, "vs", CASE, if (nzchar(RUN)) paste0("| outputs in ", ACC, "/", RUN) else "", "\n")

## output folder of the run (inputs 02/06/07/09 always stay in the dataset folder)
run_out_dir <- function(base_dir) if (nzchar(RUN)) file.path(base_dir, RUN) else base_dir
## inputs: the run folder's copy when it has one (e.g. GSE114192/TB_DM/09_wgcna_input),
## otherwise the dataset folder's -- so existing runs resolve exactly as before
in_dir <- function(base, sub) {
  r <- file.path(run_out_dir(base), sub)
  if (dir.exists(r)) r else file.path(base, sub)
}

## "ENSG00000138496.16|PARP9" -> "ENSG00000138496"; plain Ensembl IDs are left unchanged
clean_gene_ids <- function(x) sub("\\.[0-9]+$", "", sub("\\|.*$", "", x))

