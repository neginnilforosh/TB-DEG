# TB Systems-Biology Pipeline — GSE99374, GSE222001, GSE161829, GSE229020, GSE114192

Part of a network-medicine thesis project: computational drug-repurposing
approaches for tuberculosis (drug-resistant, latent, drug-susceptible TB).

## What this is

Two connected stages:

1. **DEG pipeline** (all 5 datasets): raw counts → DESeq2/limma DEG calling.
2. **WGCNA / network-medicine pipeline** (GSE114192 pilot only, for now):
   co-expression modules → enrichment-based naming → DEG×WGCNA integration →
   STRING PPI → combined multi-module network. Feeds into the later
   drug-repurposing steps (iLINCS, network proximity/RWR/diffusion) —
   not yet in this repo.

GSE114192 (Healthy_Control vs TB_Only, n=82) was chosen as the WGCNA pilot:
it's the largest clean binary TB-vs-Control comparison of the 5 datasets.
The other four remain at the DEG stage until the pilot's approach is
validated .

```
TB_DEG_pipeline/
├── scripts/
│   ├── deg_pipeline/                 <- Stage 1: DESeq2/limma DEG calling
│   │   ├── 00_functions.R            (shared: filtering, VST/PCA, DEG+volcano, etc.)
│   │   ├── GSE99374_run.R
│   │   ├── GSE222001_run.R
│   │   ├── GSE161829_run.R
│   │   ├── GSE229020_miRNA_run.R     (NanoString miRNA panel, not RNA-seq — uses limma)
│   │   ├── GSE114192_run.R
│   │   ├── sensitivity_check_no_outliers.R
│   │   └── 05_cross_dataset_upset.R  (shared/unique DEGs across datasets)
│   └── wgcna_network/                <- Stage 2: WGCNA -> ... -> combined network (GSE114192)
│       ├── WGCNA_TBmodules_step1.R           module-trait correlation, top-N module selection, MM/GS
│       ├── WGCNA_remerge_test.R              diagnostic: test alternate merge thresholds on one module
│       ├── WGCNA_TBmodules_step2_enrichment.R GO BP / KEGG / Reactome / Hallmark per module
│       ├── WGCNA_DEG_integration_step3.R     DEG x WGCNA -> "main disease genes" table
│       ├── WGCNA_STRING_step4.R              STRING PPI + centrality, per module
│       └── WGCNA_combined_network_step5.R    all modules combined into one network + bridge edges
├── GSE114192/    (see folder layout below, includes 10_wgcna_results/, 11_module_enrichment/, 12_string_ppi/)
├── GSE222001/
├── GSE161829/
├── GSE99374/     (+ sensitivity_no_outliers/ subfolder: outlier-removed re-run)
├── GSE229020/    (miRNA, NanoString nCounter panel — not gene-level RNA-seq)
├── PACKAGES.md   <- every package needed, install commands, known gotchas
└── .gitignore
```

Per-dataset folder layout (same subfolders across all 5):
`01_raw_counts` (or `01_raw_data`) → `02_metadata` → `03_filtered_counts` →
`04_normalized` → `05_sample_correlation` → `06_deg_results` →
`07_significant_degs` → [`08_shared_unique_degs`, only where a dataset has
more than one comparison] → `09_wgcna_input` → `figures/`. GSE114192 adds:
`10_wgcna_results/` → `11_module_enrichment/` → `12_string_ppi/`.

## How to run Stage 1 (DEG pipeline)

Install packages from `PACKAGES.md`, then run each `GSE*_run.R` script in
`scripts/deg_pipeline/` top to bottom. One manual check per dataset before
the metadata step: GEO's `characteristics_ch1` field naming differs between
series (e.g. `"disease state:ch1"` vs `"diagnosis:ch1"` vs `"group:ch1"`),
and each dataset's raw-count supplementary file has to be identified by
hand from its GEO page. The scripts flag exactly where to plug that in.

Design notes:
- Filtering: keep genes with ≥10 reads in ≥ (size of smallest group) samples.
- Normalization: DESeq2 VST (blind=TRUE) for QC/PCA/correlation/WGCNA input.
- DEG calling: DESeq2 Wald test per pairwise contrast vs control, apeglm LFC
  shrinkage, significance = padj < 0.05 & |log2FC| ≥ 1 (limma for GSE229020's
  miRNA panel instead of DESeq2, since it isn't RNA-seq).
- WGCNA expression export is the FULL VST matrix (not DEG-restricted).

## How to run Stage 2 (WGCNA / network-medicine pipeline, GSE114192)

Run the scripts in `scripts/wgcna_network/` in the order listed above. Each
one reads the previous step's output from `GSE114192/10_wgcna_results/` (or
`11_module_enrichment/`, `12_string_ppi/`) and is safe to re-run — it
overwrites just its own files. `WGCNA_remerge_test.R` is a one-off
diagnostic (checks whether an odd module is an over-merging artifact), not
part of the main sequential flow.

Current status / open items:
- 3 TB-significant modules selected (top 3 by |correlation|, not just FDR):
  **green** (r=+0.73) = clean, strong interferon/antiviral response
  signature (GO/KEGG/Reactome agree); **purple** (r=-0.70) = ribosome
  biogenesis/rRNA processing + chromatin/SUMOylation (less consistent);
  **blue** (r=-0.63, largest at 1377 genes) = no significant enrichment in
  GO/KEGG/Reactome despite ruling out over-merging as the cause — still
  open, Hallmark not yet tried.
- DEG × WGCNA integration: 188 "main disease genes" (143 green, 45 blue, 0
  purple — purple's module-level correlation doesn't come from any single
  strong DEG).
- STRING hubs: STAT1 (green), PARP1 (purple), MYC (blue).
- Combined network: substantial cross-module connectivity, especially
  blue↔purple and blue↔green.
- Not yet done: iLINCS TB_UP/TB_DOWN signature, network proximity/RWR/
  diffusion, drug-target mapping.

See `PACKAGES.md` for every package this needs and the specific network/
installation issues hit along the way (Bioconductor timeouts, msigdbr's
new `msigdbdf` dependency, KEGG/STRING API flakiness).
