# TB Systems-Biology Pipeline — GSE99374, GSE222001, GSE161829, GSE229020, GSE114192

Part of a network-medicine thesis project: computational drug-repurposing
approaches for tuberculosis (drug-resistant, latent, drug-susceptible TB).

## What this is

Two connected stages:

1. **DEG pipeline** (all 5 datasets): raw counts → DESeq2/limma DEG calling.
2. **Network-medicine pipeline**, configured per run in `scripts/wgcna_network/config.R`
   (pilot: GSE114192, Healthy_Control vs TB_Only, n = 82; replication: GSE161829, TBneg vs ATB, n = 45):
   WGCNA modules → enrichment → DEG×WGCNA integration → STRING PPI → combined network →
   iLINCS drug-reversal signature → drug targets (DGIdb) → network algorithms →
   final gene and drug tables, annotated with GSEA, per-drug pathways (Enrichr) and
   single-cell cell-type profiles.

GSE114192 was chosen as the pilot because it is the largest clean binary
TB-vs-Control comparison of the 5 datasets. GSE161829 (TBneg vs active TB) is the
replication run; step 14 compares the two. GSE222001, GSE99374 and GSE229020 are at the
DEG stage (GSEA is available for the gene-level ones).

```
TB-DEG/
├── scripts/
│   ├── deg_pipeline/                 <- Stage 1 (run with working directory = scripts/deg_pipeline)
│   │   ├── 00_functions.R            shared: filtering, VST/PCA, DEG + volcano, etc.
│   │   ├── GSE99374_run.R, GSE222001_run.R, GSE161829_run.R, GSE114192_run.R
│   │   ├── GSE229020_miRNA_run.R     NanoString miRNA panel, not RNA-seq — uses limma
│   │   ├── sensitivity_check_no_outliers.R
│   │   └── 05_cross_dataset_upset.R  shared/unique DEGs across the gene-level datasets
│   └── wgcna_network/                <- Stage 2, run in step order; works from any working directory
│       ├── config.R                           run presets (dataset, groups, DEG table, output folder)
│       ├── WGCNA_TBmodules_step1.R            module–trait correlation, top-3 modules, MM/GS
│       ├── WGCNA_remerge_test.R               diagnostic only: alternative merge thresholds
│       ├── WGCNA_TBmodules_step2_enrichment.R GO BP / KEGG / Reactome per module
│       ├── WGCNA_DEG_integration_step3.R      DEG × WGCNA -> "main disease genes"
│       ├── WGCNA_STRING_step4.R               STRING PPI + centrality, per module
│       ├── WGCNA_combined_network_step5.R     all modules in one network + bridge edges
│       ├── WGCNA_iLINCS_signature_step6.R     TB_UP / TB_DOWN signature (drugfindR)
│       ├── WGCNA_iLINCS_query_step7.R         iLINCS reversal compounds (drugfindR)
│       ├── WGCNA_drug_targets_step8.R         DGIdb targets, mapped onto the TB network
│       ├── WGCNA_network_algorithms_step9.R   RWR, diffusion, proximity, communities, module overlap
│       ├── WGCNA_final_tables_step10.R        FINAL gene table + FINAL drug table (merges 12 and 13 when present)
│       ├── WGCNA_GSEA_step11.R                GSEA on the full ranked DEG list, per comparison
│       ├── WGCNA_Enrichr_drugs_step12.R       Enrichr on each drug's known targets (resumable cache)
│       ├── WGCNA_singlecell_step13.R          cell-type profile of module genes and drug targets
│       └── WGCNA_compare_runs_step14.R        compares two runs: DEGs, modules, GSEA, drugs
├── GSE114192/ GSE222001/ GSE161829/ GSE99374/ GSE229020/
├── comparisons/                  step-14 output, one folder per pair of runs
├── external_data/                (NOT in git — see "Single-cell data" below)
├── PACKAGES.md                   every package, install commands, known problems
└── .gitignore
```

Per-dataset folders: `01_raw_counts` → `02_metadata` → `03_filtered_counts` → `04_normalized` →
`05_sample_correlation` → `06_deg_results` → `07_significant_degs` → [`08_shared_unique_degs`] →
`09_wgcna_input` → `figures/`. GSE114192 adds `10_wgcna_results` → `11_module_enrichment` →
`12_string_ppi` → `13_ilincs_signature` → `14_drug_targets` → `15_network_algorithms` →
`16_final_tables` → `17_gsea` → `17_pathway_annotation` → `18_singlecell`. The replication run writes
the same 10–18 folders into `GSE161829/TBneg_vs_ATB/`; GSEA results (step 11) are always in
`<dataset>/17_gsea/`, one set of files per comparison.

Raw GEO downloads (`01_raw_counts/**.gz` etc.) are not tracked; the `*_run.R` scripts re-download them.

## How to run

**Stage 1.** Install packages from `PACKAGES.md`, set the working directory to `scripts/deg_pipeline`
(RStudio: Session → Set Working Directory → To Source File Location) and run each `GSE*_run.R`.
GEO's `characteristics_ch1` naming differs between series, and each raw-count supplementary file was
identified by hand; the scripts mark where. DEG calls: DESeq2 Wald test vs control, apeglm shrinkage,
significance = padj < 0.05 & |log2FC| ≥ 1 (limma for the GSE229020 miRNA panel). WGCNA input is the full
VST matrix.

**Stage 2.** Choose the run in `scripts/wgcna_network/config.R` (`ACTIVE_RUN`, or the environment
variable `TB_RUN`, e.g. `TB_RUN=GSE114192_TB Rscript WGCNA_TBmodules_step1.R`). Each preset sets the
dataset, the control and case groups (other samples, e.g. LTBI, are left out), the DEG table and the
output folder; a new dataset or comparison only needs a new preset. Run the scripts with `Rscript` or
RStudio's Source button so a script stops at its first error. Then run the steps in order (the remerge
test is optional; step 14 compares two finished runs). Steps 4, 5,
7, 8 and 12 call web services (STRING, iLINCS, DGIdb, Enrichr), so they need internet and can be slow;
step 12 caches every drug and can simply be re-run after an interruption. After steps 12 and 13,
re-run step 10 so the final drug table picks up their columns.

**Single-cell data (step 13).** Uses the 4-week *M. tuberculosis* granuloma dataset (Broad Single Cell
Portal SCP1749) as packaged in
[RajarshiRay25/Single-Cell-Analysis---Tuberculosis](https://github.com/RajarshiRay25/Single-Cell-Analysis---Tuberculosis).
Download `data files.zip` from that repository and put `4Week_countsmatrix.mtx`, `4Week_features.tsv`,
`4Week_barcodes.tsv` and `metadata.txt` into `external_data/SCP1749/` (175 MB, git-ignored).

## Current results (GSE114192)

- **Modules** (top 3 by |correlation| with TB): **green** r = +0.73, interferon/antiviral response;
  **purple** r = −0.70, ribosome biogenesis / chromatin; **blue** r = −0.63, 1,377 genes, no significant
  over-representation.
- **Blue explained by GSEA + single-cell:** the network is unsigned, and blue mixes lymphocyte genes that
  are lower in TB (naive T-cell sets, e.g. HAY_BONE_MARROW_NAIVE_T_CELL NES −2.59) with myeloid genes
  that are higher in TB. In the granuloma single-cell data, blue genes lower in TB are expressed mostly
  in T/B cells (mean share T 18.7 %, B 11.3 %) and blue genes higher in TB mostly in macrophages and
  neutrophils (13.8 %, 11.8 %), against 9.1 % expected. Bulk data cannot tell fewer cells from lower
  expression per cell.
- **Strongest pathways (GSEA Hallmark):** IFN-γ response NES +2.72, IFN-α response +2.69.
- **Disease genes:** 188 (143 green, 45 blue, 0 purple). **Hubs:** STAT1 (green), PARP1 (purple), MYC (blue).
- **Drugs:** 4,251 iLINCS reversal compounds → 395 with a DGIdb target in the TB network. Final drug
  table: 367 with a reversal score, 352 with Enrichr pathways (43 have < 3 known targets), 244 with a
  cell-type profile.
- **Network proximity is not informative here:** no drug's targets are closer to the disease genes than
  degree-matched random targets (mean Z ≈ 0). It stays in the table but is not used for ranking.

## Replication (GSE161829, TBneg vs ATB) — step 14

- **Pathways replicate:** 14 Hallmark sets are significant with the same sign in both runs, led by
  IFN-γ response (NES +2.72 / +2.49) and IFN-α response (+2.69 / +2.51), with IL6-JAK-STAT3, complement,
  allograft rejection, oxidative phosphorylation and mTORC1; 86 cell-type and 76 Reactome sets also
  replicate. E2F targets, G2M checkpoint, MYC targets V2 and mitotic spindle go in opposite directions
  (lower in the pilot, higher in GSE161829).
- **Genes and modules replicate only weakly:** log2FC correlation 0.22; 5 genes significant in both
  (0.6 expected by chance), all in the same direction; module overlaps are small.
- **Drug candidates do not replicate:** reversal hits overlap no more than expected by chance (350 vs
  380), and the same LINCS signatures score in opposite directions (Spearman −0.42). The two TB_UP
  signatures share only STAT1, PSME2 and UBE2L6, and about half of the 18-gene GSE161829 signature is
  cell-cycle genes (the pilot's has none), which most likely drives the disagreement.
- GSE161829 has very few L1000 genes going down, so its iLINCS query used TB_UP only; the comparison
  therefore uses UP-only results for both runs.

## Known limitations

- Replication so far covers one second cohort; drug candidates are not yet reproducible across cohorts.
- DGIdb target lists include drug-metabolising enzymes and transporters (CYPs, SLC, ABC), so roughly
  15 % of the top KEGG/Reactome labels per drug are ADME-type rather than mechanism-of-action.
- The single-cell reference is macaque lung granuloma (2 animals, coarse labels), not human blood;
  read it as "which immune cell type a gene belongs to", not as blood composition. Its p-values are
  descriptive only.
