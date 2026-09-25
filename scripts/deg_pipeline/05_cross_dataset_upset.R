## =============================================================
## 05_cross_dataset_upset.R
## Shared/unique DEG comparison ACROSS the four datasets/comorbidity
## contexts, run after all four *_run.R scripts have produced their
## 07_significant_degs/*.csv files.
## =============================================================
source("00_functions.R")

read_gene_col <- function(path) read.csv(path, stringsAsFactors = FALSE)$gene

deg_lists <- list(
  GSE222001_Healthy_vs_Active  = read_gene_col("../GSE222001/07_significant_degs/GSE222001_Healthy_vs_Active_sigDEGs.csv"),
  GSE222001_Healthy_vs_Latent  = read_gene_col("../GSE222001/07_significant_degs/GSE222001_Healthy_vs_Latent_sigDEGs.csv"),
  GSE161829_TBneg_vs_LTBI      = read_gene_col("../GSE161829/07_significant_degs/GSE161829_TBneg_vs_LTBI_sigDEGs.csv"),
  GSE161829_TBneg_vs_ATB       = read_gene_col("../GSE161829/07_significant_degs/GSE161829_TBneg_vs_ATB_sigDEGs.csv"),
  GSE99374_HC_vs_LTBI          = read_gene_col("../GSE99374/07_significant_degs/GSE99374_HC_vs_LTBI_sigDEGs.csv"),
  GSE229020_HC_vs_LatentTB     = read_gene_col("../GSE229020/07_significant_degs/GSE229020_HC_vs_LatentTB_sigDEGs.csv"),
  GSE229020_HC_vs_DrugResistant = read_gene_col("../GSE229020/07_significant_degs/GSE229020_HC_vs_DrugResistant_sigDEGs.csv")
)

out_dir <- file.path("..", "cross_dataset_comparison")
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)

shared_unique_degs(deg_lists,
                    out_csv = file.path(out_dir, "cross_dataset_shared_unique_DEGs.csv"),
                    out_png = file.path(out_dir, "figures", "cross_dataset_upset.png"),
                    title = "Shared and unique DEGs across TB datasets/comorbidities")


