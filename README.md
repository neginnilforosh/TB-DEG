# Network medicine for tuberculosis

This repository asks a simple question: **which host genes drive tuberculosis, and which existing drugs might push them back towards normal?**

It starts from blood gene-expression data of people with TB and healthy controls, finds the genes and gene networks that change, ranks them by several independent kinds of evidence, and then looks for approved or experimental compounds whose effects run opposite to the disease. The same analysis is being repeated for TB with diabetes and TB with HIV, so the three contexts can be compared at the end.

Everything here is computational. A drug that ranks highly is a hypothesis worth testing, not a treatment.

---

## What the analysis does

1. **Differential expression.** Compare patients with controls and find the genes that go up or down (DESeq2 for RNA-seq, limma for arrays and miRNA panels).
2. **Co-expression modules.** Group genes that rise and fall together across all samples (WGCNA), keep the groups most associated with disease, and name them by the pathways they are enriched for.
3. **Protein networks.** Place the module genes on the STRING protein-interaction network and measure how central each one is (degree, betweenness, random walk with restart, network diffusion).
4. **Gene ranking.** Rank every gene on each measure separately, then flag the genes that score well on several of them at once.
5. **Drug search.** Look for compounds whose expression signature reverses the disease signature (iLINCS), map their known targets onto the disease network (DGIdb), and rank them by how many independent lines of evidence agree.
6. **Supporting layers.** Pathway analysis of the whole ranked gene list (GSEA), cell-type mapping (single-cell data), chromatin accessibility (ATAC-seq), plasma metabolites, literature labels (CTD), and checks in independent cohorts.

## Datasets

| Dataset | Sample type and platform | Groups | Used for |
|---|---|---|---|
| GSE114192 | Whole blood, RNA-seq | Healthy 36, TB 46, diabetes 52, TB + diabetes 61 | Main analysis (TB) and TB-diabetes |
| GSE50834 | PBMC, Illumina arrays | HIV 23, HIV + TB 21 | TB-HIV (planned) |
| GSE161829 | Memory CD4+ T cells, RNA-seq | No TB 20, latent TB 40, active TB 25 | Independent check |
| GSE222001 | NK cells, RNA-seq | Healthy 5, latent TB 5, active TB 5 | Independent check (small) |
| GSE99374 | CD8+ T cells, RNA-seq | Healthy 19, latent TB 21 | Early analysis |
| GSE229020 | miRNA panel (NanoString) | Healthy 4, latent 6, drug-sensitive 6, drug-resistant 4 | miRNA differential expression |
| GSE165709 | Macrophages, RNA-seq + ATAC-seq | Healthy and people living with HIV, infected vs uninfected | Do the TB genes respond to infection, and is their chromatin opened? |
| ST000974 | Plasma metabolomics | TB cases vs controls | Metabolite layer |
| SCP1749 | Single-cell RNA, lung granuloma | — | Which cell types express the module genes |

## What we have found so far (TB vs healthy)

- **Interferon signalling is the strongest signal.** The top pathways are the interferon-gamma and interferon-alpha responses, and the most connected gene in the disease network is STAT1.
- **The blood changes in composition, not only in gene activity.** One large module mixes lymphocyte genes that go down with neutrophil and macrophage genes that go up. Single-cell data and cell-type gene sets agree with this reading.
- **The core TB genes respond directly to infection.** Most of the top-ranked genes also respond when lung macrophages are infected with *M. tuberculosis* in the lab (GSE165709), and several of them, mostly interferon-driven, respond less in cells from people living with HIV.
- **Interferon links genes to metabolites.** Tryptophan is lower in the plasma of TB cases (ST000974), while IDO1, KYNU and WARS1, the enzymes that use it, are higher in TB blood. All three are switched on by interferon.
- **Drug predictions are the least stable part.** Ranking compounds by signature reversal changes noticeably between similar inputs, and network proximity does not separate drugs from random after correcting for multiple tests. Drug candidates should therefore be read together with their other evidence, not on their own.

## Where things stand

| Part | Status |
|---|---|
| TB vs healthy | Complete |
| Supporting layers (GSEA, single-cell, ATAC-seq, metabolomics, CTD) | Complete |
| TB + diabetes | Full pipeline has run; results are preliminary because healthy and diabetes-only samples come from fewer recruitment sites than the TB samples, which has to be accounted for |
| Independent cohorts (GSE161829, GSE222001) | Sample labels corrected; all analyses redone. Pathways and the interferon module replicate; drug candidates do not |
| TB + HIV | Not started. The dataset's two groups were processed at different times, so the best way to use it is still being decided |
| Comparison of TB, TB-diabetes and TB-HIV | Waits for TB-HIV |

## Repository layout

```
scripts/
  deg_pipeline/      one script per dataset: download, filter, normalise, differential expression
  wgcna_network/     the network analysis, run step by step (step1 ... step24); config.R picks the dataset
  helpers/           small Python tools for annotated tables and checks
  record_reproducibility.R
GSE114192/           main analysis; TB_DM/ holds the TB-diabetes analysis in the same structure
GSE161829/ GSE222001/ GSE99374/ GSE229020/   other cohorts
GSE50834/            TB-HIV dataset audit
external_results/    GSE165709 RNA-seq and ATAC-seq results
reproducibility/     package versions, parameters, seeds, external resources
```

Inside each dataset the folders are numbered in the order they are produced: `01_raw_counts` → `02_metadata` → … → `09_wgcna_input` for the expression stage, then `10_wgcna_results` → … → `19_consensus` for the network stage, and `21`–`25` for the supporting layers. Figures are in `figures/`; the full figure set for each analysis is in `figures_manuscript/` (600 dpi PNG and PDF).

The two final tables most people will want are:

- `GSE114192/19_consensus/GSE114192_GenePriority_all_modules.csv` — every module gene with its expression change, module, network scores, the rank on each measure and the overall priority
- `GSE114192/19_consensus/GSE114192_DrugPriority.csv` — every candidate compound with its reversal score, targets, modules hit and overall priority

## How to run

1. Install the packages listed in `PACKAGES.md`.
2. **Expression stage:** set the working directory to `scripts/deg_pipeline` and run the `GSE*_run.R` script for a dataset. It downloads the data from GEO, matches samples to their labels by name, and writes folders `01`–`09`.
3. **Network stage:** choose the analysis in `scripts/wgcna_network/config.R` (or set `TB_RUN`, for example `TB_RUN=GSE114192_TB Rscript WGCNA_TBmodules_step1.R`) and run the steps in order. Steps 4, 5, 7, 8 and 12 query web services and need an internet connection.
4. **Single-cell data:** download the SCP1749 four-week granuloma matrix from [this repository](https://github.com/RajarshiRay25/Single-Cell-Analysis---Tuberculosis) into `external_data/SCP1749/`.
5. **CTD:** download `CTD_curated_genes_diseases.csv.gz` (and optionally `CTD_chem_gene_ixns.csv.gz`) from ctdbase.org/downloads into `external_data/CTD/`.

Large raw files and downloads stay out of Git; the scripts fetch them again.

## Limitations worth knowing

- Bulk blood data cannot tell a change in cell numbers from a change in gene activity inside cells.
- The independent cohorts use sorted cell types (CD4, CD8, NK cells), not whole blood, so a gene that does not replicate there may simply not be expressed in those cells.
- Web services (STRING, iLINCS, DGIdb, Enrichr, KEGG) return their current version at the time of the query; `reproducibility/` records what was used.
- The single-cell reference is macaque lung tissue from two animals and is used only to say which cell type a gene belongs to.
- CTD contains very few curated TB genes, so the absence of a TB label there does not mean a gene is new.
