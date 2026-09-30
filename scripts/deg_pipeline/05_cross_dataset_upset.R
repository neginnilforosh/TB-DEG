## =============================================================
## 05_cross_dataset_upset.R
## Shared/unique DEG comparison across the gene-level datasets
## (GSE222001, GSE161829, GSE99374, GSE114192). Run after the *_run.R
## scripts have produced their 07_significant_degs/*.csv files.
## =============================================================
source("00_functions.R")

read_gene_col <- function(path) read.csv(path, stringsAsFactors = FALSE)$gene

deg_lists <- list(
  GSE222001_Healthy_vs_Active  = read_gene_col("../../GSE222001/07_significant_degs/GSE222001_Healthy_vs_Active_sigDEGs.csv"),
  GSE222001_Healthy_vs_Latent  = read_gene_col("../../GSE222001/07_significant_degs/GSE222001_Healthy_vs_Latent_sigDEGs.csv"),
  GSE161829_TBneg_vs_LTBI      = read_gene_col("../../GSE161829/07_significant_degs/GSE161829_TBneg_vs_LTBI_sigDEGs.csv"),
  GSE161829_TBneg_vs_ATB       = read_gene_col("../../GSE161829/07_significant_degs/GSE161829_TBneg_vs_ATB_sigDEGs.csv"),
  GSE99374_HC_vs_LTBI          = read_gene_col("../../GSE99374/07_significant_degs/GSE99374_HC_vs_LTBI_sigDEGs.csv"),
  GSE114192_HC_vs_TBOnly       = read_gene_col("../../GSE114192/07_significant_degs/GSE114192_HealthyControl_vs_TBOnly_sigDEGs.csv")
  # GSE229020 left out: NanoString miRNA panel -- its tables list miRNAs, not genes
)

out_dir <- file.path("..", "..", "cross_dataset_comparison")
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)

shared_unique_degs(deg_lists,
                    out_csv = file.path(out_dir, "cross_dataset_shared_unique_DEGs.csv"),
                    out_png = file.path(out_dir, "figures", "cross_dataset_upset.png"),
                    title = "Shared and unique DEGs across TB datasets/comorbidities")


