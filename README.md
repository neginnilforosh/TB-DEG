# TB DEG Pipeline — GSE222001, GSE161829, GSE99374, GSE229020

## What this is
A complete, ready-to-run R/Bioconductor pipeline (DESeq2-based) that takes
each dataset from raw counts through DEG calling, producing every file and
figure requested, organized in one folder per dataset:

```
TB_DEG_pipeline/
├── scripts/
│   ├── 00_functions.R              <- shared functions (filtering, VST/PCA,
│   │                                   correlation heatmap, DEG+volcano,
│   │                                   sig-DEG heatmap, UpSet, WGCNA export)
│   ├── GSE222001_run.R
│   ├── GSE161829_run.R
│   ├── GSE99374_run.R
│   ├── GSE229020_run.R
│   └── 05_cross_dataset_upset.R    <- shared/unique DEGs across datasets
├── GSE222001/
│   ├── 01_raw_counts/
│   ├── 02_metadata/
│   ├── 03_filtered_counts/
│   ├── 04_normalized/
│   ├── 05_sample_correlation/
│   ├── 06_deg_results/
│   ├── 07_significant_degs/
│   ├── 08_shared_unique_degs/
│   ├── 09_wgcna_input/
│   └── figures/
├── GSE161829/  (same subfolders)
├── GSE99374/   (same subfolders)
└── GSE229020/  (same subfolders)
```


## How to actually run this

copy this folder to any machine with internet access and R installed (RStudio, a university
cluster, Google Colab-with-R, etc.), install the packages listed at the
top of `00_functions.R`, then run each `GSE*_run.R` script in `scripts/`
top to bottom. Each commented-out block corresponds 1:1 to one deliverable
in the brief (raw matrix → metadata → filtered matrix → VST/PCA →
correlation heatmap → DEG/volcano → sig-DEG/heatmap → shared-unique/UpSet
→ WGCNA exports). Uncomment as you confirm each dataset's exact column
names from GEO (see note below).


one manual check is needed per dataset before the metadata
step: GEO's `characteristics_ch1` field naming differs between series
(e.g. `"disease state:ch1"` vs `"diagnosis:ch1"` vs `"group:ch1"`), and
each dataset needs its raw-count supplementary file identified by hand
(filename shown on the dataset's GEO page). The scripts flag exactly
where to plug that in.

## Design notes
- Filtering: keep genes with ≥10 reads in ≥ (size of smallest group) samples.
- Normalization: DESeq2 VST (blind=TRUE) for QC/PCA/correlation/WGCNA input.
- DEG calling: DESeq2 Wald test per pairwise contrast vs the control group,
  apeglm LFC shrinkage, significance = padj < 0.05 & |log2FC| ≥ 1.
- WGCNA expression export is the FULL VST matrix (not DEG-restricted), per
  the brief.
- Cross-dataset shared/unique DEGs and the UpSet plot are handled in
  `scripts/05_cross_dataset_upset.R`, after all four datasets are run —
  edit the comparison list there once comorbidity metadata (TB-DM, TB-HIV,
  TB-LC, drug-resistant, etc.) is attached to each dataset's trait file.
