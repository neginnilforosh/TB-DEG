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
| GSE165709 | Alveolar macrophages (cell type as recorded in the analysis scripts), RNA-seq + ATAC-seq | Healthy donors and people living with HIV; cells infected with *M. tuberculosis* in the lab vs uninfected | Do the TB genes respond to infection, and is their chromatin opened? |
| ST000974 | Plasma metabolomics (Metabolomics Workbench) | Samples labelled case (92) and control (274) by the study; how the case group is described here is still open (decision D7) | Metabolite layer (supporting) |
| SCP1749 | Single-cell RNA, lung granuloma | — | Which cell types express the module genes |

## What we have found so far (TB vs healthy)

- **Interferon signalling is the strongest signal.** The top pathways are the interferon-gamma and interferon-alpha responses (GSEA NES +2.72 and +2.69). In the combined disease network STAT1 has the highest eigenvector centrality; MYC has the most direct interactions.
- **The blood changes in composition, not only in gene activity.** One large module mixes lymphocyte genes that go down with neutrophil and macrophage genes that go up. Single-cell data and cell-type gene sets agree with this reading.
- **The core TB genes respond directly to infection.** 44 of the 47 High-priority genes also change significantly when alveolar macrophages from healthy donors are infected with *M. tuberculosis* in the lab, and 32 of 47 in cells from people living with HIV (GSE165709). Some genes appear to respond less in the cells from people living with HIV, but a formal test of the difference between the two groups found no gene or chromatin region significant after multiple-testing correction, so this is a descriptive observation only.
- **Interferon may link genes to metabolites.** In ST000974, tryptophan is lower in plasma samples labelled as cases than in controls (log2 fold change −0.16, adjusted p = 0.006), while IDO1, KYNU and WARS1, enzymes that use tryptophan, are higher in TB blood in GSE114192 (log2 fold change 1.2–1.4, adjusted p < 1e-5). All three are induced by interferon. The two datasets come from different people, and the ST000974 case definition is still open (D7), so this is a supporting observation, not a finding on its own.
- **Drug predictions are the least stable part.** Ranking compounds by signature reversal changes noticeably between similar inputs, and network proximity does not separate drugs from random after correcting for multiple tests. Drug candidates should therefore be read together with their other evidence, not on their own.

## Where things stand

| Part | Status |
|---|---|
| TB vs healthy | Expression, modules, networks and drug evidence are complete. The module-selection rule (D3) and the final gene and drug priority rules (D4, D5) are still open, so the priority tables are provisional |
| Supporting layers (GSEA, single-cell, ATAC-seq, metabolomics, CTD) | Run, not yet validated. The ATAC-seq peak-to-gene annotation method still has to be documented, the ST000974 case definition is open (D7), and the CTD annotation has not been independently re-checked |
| TB + diabetes | Full pipeline has run; results are preliminary because healthy and diabetes-only samples come from fewer recruitment sites than the TB samples. A site sensitivity analysis (`GSE114192/TB_DM/sensitivity_site/`) finds the TB-DM vs DM changes largely stable (52–54% of significant genes kept with smaller samples, at least 99.3% with the same direction, all 38 High genes kept), but the TB × DM interaction is not: 4 genes without a site term, 1 in South Africa and Romania only, 0 once site is modelled. The final design is still open (D1, D13) |
| Independent cohorts (GSE161829, GSE222001) | Sample labels corrected; all analyses redone. Pathways and the interferon module replicate; drug candidates do not |
| TB + HIV | Not started. The dataset's two groups were processed at different times, so the best way to use it is still being decided (D10). It stays blocked until the open pre-HIV decisions and checks are complete |
| Comparison of TB, TB-diabetes and TB-HIV | Waits for TB-HIV |
| Code corrections and reproducibility updates (drug-name matching, proximity P_BH, removal of DSigDB self-hits, nested-CV column labels, GSEA warning log, package checks) | In the scripts and tested; the stored drug, proximity, Enrichr, consensus and ML-label outputs were produced before them and have not been regenerated yet |


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
comparisons/         reference vs replication run comparison (step 14)
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
- **Site confounding in the main TB vs healthy analysis (GSE114192).** Healthy controls were recruited only in South Africa (24) and Romania (12), while the TB patients come from all four sites (Indonesia 14, Peru 11, Romania 10, South Africa 11); chi-square 29.2 (df 3), Cramér's V 0.60, from the counts in `GSE114192/TB_DM/00_design/site_by_group.csv`. Some TB vs healthy differences may therefore partly reflect site. The main analysis does not adjust for site, and its modules, networks and machine-learning estimates inherit this limitation. A sensitivity analysis (`GSE114192/sensitivity_site/`) supports the main results: with site in the model, or with only South Africa and Romania, 52–57% of the 300 significant genes stay significant at the smaller sample size, 99.7% keep their direction, all 47 High-priority genes stay significant in the same direction, and the interferon-gamma NES is 2.72–2.82.
- **Machine-learning estimates.** Discovery AUCs computed on the same samples used to choose the genes are optimistic (signature score 0.975). Nested cross-validation, which re-selects the genes inside every training fold (5 folds × 10 repeats), gives 0.967 for both models; the interval stored with it is the spread across repeats, not a confidence interval (the current output file still labels it `CI_low`/`CI_high`; the corrected script labels it `Repeat_AUC_q025`/`Repeat_AUC_q975` from the next run on). The external checks use small cohorts of sorted cells: GSE161829 active TB vs no TB, AUC 0.72 (bootstrap 95% CI 0.57–0.86; not better than random gene sets, p = 0.08) and active vs latent TB 0.74 (0.60–0.86); GSE222001 has 5 people per group (0.88, 0.60–1.00). None of these estimates accounts for site (`GSE114192/24_ml_validation/`).
- **Network proximity is a negative result.** None of the 395 candidate drugs is closer to the disease network than random target sets after multiple-testing correction (FDR < 0.05; `GSE114192/15_network_algorithms/GSE114192_NetworkProximity_withFDR.csv`). Proximity is therefore not used as validated evidence; how it should enter the drug priority tiers is still open (D5).
- **Enrichr results are annotation.** Each drug's known targets are tested for pathway enrichment (corrected within each drug and library, against the genome background). An enriched pathway shows which pathways the targets belong to; it does not show that the drug reverses that pathway, and enrichment across many drugs is not biological validation. The current annotation tables still include each drug's own DSigDB entry (a circular match); the corrected script removes these from the next run on (checked offline against the saved results, not yet in a full online run).
- **The TB × DM interaction is exploratory.** Its few significant genes disappear once recruitment site is accounted for (see above), and the diabetes groups differ in treatment and body-mass index, which cannot be adjusted for.
- **GSEA gene mapping.** Before GSEA, Ensembl IDs are mapped to gene symbols; about 1% of IDs map to several symbols, and each symbol then carries the same statistic (the "ties" warnings). The current method is reproducible: re-running it gives the stored results exactly. In an exploratory comparison (reference and GSE161829 TB comparisons, not a validated alternative), two other mapping policies left Hallmark results almost unchanged (0-2 gene sets crossing padj 0.05, no direction changes) but moved 11-30 Reactome gene sets across padj 0.05, with NES rank correlation of 0.998 or more. Reactome results close to the 0.05 cut-off should be read with this in mind. The mapping policy is unchanged; whether to change it is still open.
- **Random draws in step 9.** Network proximity uses 1,000 degree-matched random target sets per drug, drawn in sequence from one seed, and community detection (Louvain) is not separately seeded. When the drug list changes, proximity values move within Monte Carlo error (in tests, Z-scores correlated at 0.998 or more with the previous run; nominal p < 0.05 counts moved, e.g. 4 to 2 in the reference) and community labels change. Options for making these draws stable per drug have been identified; none has been adopted.
- **Web-service versions.** All existing STRING network outputs were produced through the unversioned STRING endpoint, so the STRING version behind them is not known. The scripts now pin STRING 12.5; this applies to new runs only. iLINCS, DGIdb, Enrichr and KEGG are not version-pinned; `reproducibility/external_resource_versions.csv` records what is known.
- The single-cell reference is macaque lung tissue from two animals and is used only to say which cell type a gene belongs to.
- CTD contains very few curated TB genes, so the absence of a TB label there does not mean a gene is new.
