# Network medicine for tuberculosis

This repository asks a simple question: **which host genes drive tuberculosis, and which existing drugs might push them back towards normal?**

It starts from blood gene-expression data of people with TB and healthy controls, finds the genes and gene networks that change, ranks them by several independent kinds of evidence, and then looks for approved or experimental compounds whose effects run opposite to the disease. The same analysis is being repeated for TB with diabetes and TB with HIV, so the three contexts can be compared at the end.

Everything here is computational. A drug that ranks highly is a hypothesis worth testing, not a treatment.

---

## What the analysis does

1. **Differential expression.** Compare patients with controls and find the genes that go up or down (DESeq2 for RNA-seq, limma for arrays and miRNA panels).
2. **Co-expression modules.** Group genes that rise and fall together across all samples (WGCNA), keep the groups most associated with disease, and name them by the pathways they are enriched for.
3. **Protein networks.** Place the module genes on the STRING protein-interaction network and measure how central each one is (degree, betweenness, random walk with restart, network diffusion).
4. **Gene ranking.** Group the evidence into independent families (expression, co-expression, network, replication in other cohorts, response to infection), give each family one score, combine them, and check how much the ranking moves when one family is left out.
5. **Drug search.** Look for compounds whose expression signature reverses the disease signature (iLINCS), map their known targets onto the disease network (DGIdb), rank them by evidence families in the same way, and examine the best candidates one by one in evidence cards.
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
| ST000974, ST000975, ST000976 | Plasma, serum and diluted plasma metabolomics from the GC6-74 household-contact cohort (Metabolomics Workbench) | Contacts who developed active TB within two years (progressors) vs contacts who stayed healthy, matched about 1:4 and sampled while still TB-free; partly the same participants across the three studies | Metabolite layer (supporting) |
| SCP1749 | Single-cell RNA, lung granuloma | — | Which cell types express the module genes |

## What we have found so far (TB vs healthy)

- **Interferon signalling is the strongest signal.** The top pathways are the interferon-gamma and interferon-alpha responses (GSEA NES +2.72 and +2.69). In the combined disease network STAT1 has the highest eigenvector centrality; MYC has the most direct interactions.
- **The blood changes in composition, not only in gene activity.** One large module mixes lymphocyte genes that go down with neutrophil and macrophage genes that go up. Single-cell data and cell-type gene sets agree with this reading.
- **The core TB genes respond directly to infection.** 44 of the 47 High-priority genes also change significantly when alveolar macrophages from healthy donors are infected with *M. tuberculosis* in the lab, and 32 of 47 in cells from people living with HIV (GSE165709). Some genes appear to respond less in the cells from people living with HIV, but a formal test of the difference between the two groups found no gene or chromatin region significant after multiple-testing correction, so this is a descriptive observation only.
- **Interferon may link genes to metabolites.** In the GC6-74 cohort the kynurenine/tryptophan ratio, a read-out of the interferon-induced enzyme IDO1, did not differ between future TB cases and controls at enrolment or more than five months before diagnosis, but rose within five months of diagnosis in plasma and in serum (adjusted p 0.12 and 0.03), with a matching fall in tryptophan. IDO1, KYNU and WARS1 are higher in active-TB blood in GSE114192. The metabolite and gene data come from different people, so this is a supporting observation; sample numbers close to diagnosis are small.
- **Drug predictions are the least stable part.** Ranking compounds by signature reversal changes noticeably between similar inputs, and network proximity does not separate drugs from random after correcting for multiple tests. Only a few compounds reverse both halves of the signature consistently across LINCS experiments and act on a high-priority gene, and many top candidates also have cell-division targets, so part of the drug ranking reflects antiproliferative effects rather than TB biology.
- **TB with and without diabetes share an interferon-gamma core.** 25 genes are high-priority in both contexts, nearly all interferon-gamma related (STAT1, GBP1, GBP5, CASP1, CD274, ...). The type I interferon genes (ISG15, IFIT3, IFI6, OASL, ...) still go up in TB with diabetes but form a separate small module there instead of joining the core, which fits the weaker interferon-alpha response in diabetes.

## Where things stand

| Part | Status |
|---|---|
| TB vs healthy | Complete. Genes and drugs are ranked by evidence families (`20_evidence_families/`), with the criteria count of step 15 kept for comparison. A site sensitivity analysis (`GSE114192/sensitivity_site/`) keeps all 47 High genes of the criteria count significant with unchanged direction and the interferon responses as the top pathways |
| Supporting layers (GSEA, single-cell, ATAC-seq, metabolomics, CTD) | Complete. The metabolomics result does not change when smoking is added as a covariate (serum) |
| TB + diabetes | Complete. The main comparison is TB with diabetes vs TB alone on all samples with recruitment site as a covariate (`GSE114192/TB_DM/tb_vs_tbdm_site_adjusted/`); both groups come from all four sites. The interferon-gamma and interferon-alpha responses are weaker in TB with diabetes (NES −1.84 and −2.22), and the TB × diabetes interaction, fitted on the two sites that have all four groups and adjusted for site, gives the same result (`GSE114192/TB_DM/sensitivity_site/`). Neutrophil-granule and red-cell genes are higher in TB with diabetes, which may reflect blood cell composition |
| Independent cohorts (GSE161829, GSE222001) | Sample labels corrected; all analyses redone. Pathways and the interferon module replicate; drug candidates do not |
| TB + HIV | Not started. The first candidate (GSE50834) has disease group and processing batch almost completely overlapping, so alternative datasets (GSE37250, GSE107104, GSE162164) have been compared and the choice is being agreed |
| Comparison across contexts | TB vs TB with diabetes complete (`comparisons/contexts/`); TB with HIV will be added when its dataset is chosen |


## Repository layout

```
scripts/
  deg_pipeline/      one script per dataset: download, filter, normalise, differential expression
  wgcna_network/     the network analysis, run step by step (step1 ... step27); config.R picks the dataset
  helpers/           small Python tools for annotated tables and checks
  record_reproducibility.R
GSE114192/           main analysis; TB_DM/ holds the TB-diabetes analysis in the same structure
GSE161829/ GSE222001/ GSE99374/ GSE229020/   other cohorts
GSE50834/            TB-HIV dataset audit
external_results/    GSE165709 RNA-seq and ATAC-seq results
reproducibility/     package versions, parameters, seeds, external resources
comparisons/         reference vs replication run (step 14) and comparison across clinical contexts (contexts/, step 27)
```

Inside each dataset the folders are numbered in the order they are produced: `01_raw_counts` → `02_metadata` → … → `09_wgcna_input` for the expression stage, `10_wgcna_results` → … → `19_consensus` for the network stage, `20_evidence_families` for the final ranking, and `21`–`26` for the supporting layers, validation and drug evidence cards. Figures are in `figures/`; the full figure set for each analysis is in `figures_manuscript/` (600 dpi PNG and PDF).

The final tables most people will want are:

- `GSE114192/20_evidence_families/GSE114192_Gene_evidence_families.csv` — every module gene with its family scores, evidence score, rank, rank range across variants and tier
- `GSE114192/20_evidence_families/GSE114192_Drug_evidence_families.csv` and `GSE114192/26_drug_evidence_cards/` — every candidate compound ranked the same way, and a card for each of the best ones
- `GSE114192/19_consensus/` — the earlier criteria-count ranking, kept for comparison

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
- **Network proximity is a negative result.** None of the 395 candidate drugs is closer to the disease network than random target sets after multiple-testing correction (FDR < 0.05; `GSE114192/15_network_algorithms/GSE114192_NetworkProximity_withFDR.csv`). Proximity is therefore not used as validated evidence; how it should enter the drug priority tiers is still being decided.
- **Enrichr results are annotation.** Each drug's known targets are tested for pathway enrichment (corrected within each drug and library, against the genome background). An enriched pathway shows which pathways the targets belong to; it does not show that the drug reverses that pathway, and enrichment across many drugs is not biological validation. The current annotation tables still include each drug's own DSigDB entry (a circular match); the corrected script removes these from the next run on (checked offline against the saved results, not yet in a full online run).
- **The TB × DM interaction is a gene-set result.** No single gene survives once recruitment site is accounted for; the interferon gene sets do (see above). The diabetes groups also differ in treatment and body-mass index, which cannot be adjusted for.
- **GSEA gene mapping.** Before GSEA, Ensembl IDs are mapped to gene symbols; about 1% of IDs map to several symbols, and each symbol then carries the same statistic (the "ties" warnings). The current method is reproducible: re-running it gives the stored results exactly. In an exploratory comparison (reference and GSE161829 TB comparisons, not a validated alternative), two other mapping policies left Hallmark results almost unchanged (0-2 gene sets crossing padj 0.05, no direction changes) but moved 11-30 Reactome gene sets across padj 0.05, with NES rank correlation of 0.998 or more. Reactome results close to the 0.05 cut-off should be read with this in mind. The mapping policy is unchanged; whether to change it is still open.
- **Random draws in step 9.** Network proximity uses 1,000 degree-matched random target sets per drug, drawn in sequence from one seed, and community detection (Louvain) is not separately seeded. When the drug list changes, proximity values move within Monte Carlo error (in tests, Z-scores correlated at 0.998 or more with the previous run; nominal p < 0.05 counts moved, e.g. 4 to 2 in the reference) and community labels change. Options for making these draws stable per drug have been identified; none has been adopted.
- **Web-service versions.** All existing STRING network outputs were produced through the unversioned STRING endpoint, so the STRING version behind them is not known. The scripts now pin STRING 12.5; this applies to new runs only. iLINCS, DGIdb, Enrichr and KEGG are not version-pinned; `reproducibility/external_resource_versions.csv` records what is known.
- The single-cell reference is macaque lung tissue from two animals and is used only to say which cell type a gene belongs to.
- CTD contains very few curated TB genes, so the absence of a TB label there does not mean a gene is new.
