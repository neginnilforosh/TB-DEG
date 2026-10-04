# Environment setup

R version used during development: 4.3+ (also tested clean on 4.5.x).
Run everything on a machine with internet access — GEO downloads, and
later, live calls to KEGG/STRING, need it.

## Install everything in one go

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor
BiocManager::install(c(
  "DESeq2", "GEOquery", "apeglm", "limma", "WGCNA",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA"
))

# CRAN
install.packages(c(
  "pheatmap", "ggplot2", "ggrepel", "RColorBrewer", "UpSetR",
  "matrixStats", "dplyr", "httr", "jsonlite", "igraph", "enrichR", "Matrix"
))

# drugfindR (iLINCS, steps 6-7) is distributed via r-universe
install.packages("drugfindR", repos = c("https://cogdisreslab.r-universe.dev", "https://cloud.r-project.org"))

# msigdbr changed its data delivery in v10+: the gene-set data itself moved
# out of the CRAN package into a separate 'msigdbdf' package (too big for
# CRAN). Install BOTH from the maintainer's r-universe repo:
install.packages(c("msigdbr", "msigdbdf"),
                  repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))
```

## Known gotchas (all hit at least once this project)

- **`BiocManager::install()` times out on `bioconductor.org/config.yaml`**
  ("Bioconductor version cannot be validated; no internet connection?"): this
  is a network issue on your machine/network, not a code problem. Usually
  transient — retry in a few minutes. If it persists: test
  `download.file("https://cran.r-project.org", "test.html")` — if that also
  fails, it's general connectivity, not Bioconductor specifically; try a
  different network (mobile hotspot, VPN off/on). On a university HPC
  cluster, outbound internet is often blocked on purpose — check for a
  preloaded module (`module avail r`, `module avail bioconductor`) instead
  of installing yourself.
- **`library(clusterProfiler)` errors with "no package called..."**: means
  the install above never actually completed (usually because of the
  Bioconductor connectivity issue above). Run `library(clusterProfiler)`
  alone and check for a real error before re-running any pipeline script.
- **KEGG lookups inside `enrichKEGG()` time out**: it fetches
  `rest.kegg.jp` live every time (not bundled offline). Usually transient —
  just retry. If it keeps failing, `options(timeout = 300)` before
  running.
- **`Could not resolve host: datasetStatistics`** from enrichR: the package sets its connection
  options only when attached with `library(enrichR)`; calling `enrichR::listEnrichrDbs()` without
  attaching leaves them empty. Step 12 attaches it itself.
- **`object '.Random.seed' not found` in GSEA**: a fresh `Rscript` session has no random-number state yet;
  step 11 calls `set.seed(42)` before GSEA.
- **`entrez_gene` not found (msigdbr)**: msigdbr >= 10 renamed it to `ncbi_gene`; step 2 accepts both.
- **iLINCS / DGIdb / Enrichr calls fail or time out**: these are live web services. Retry later;
  step 12 keeps a per-drug cache in `17_pathway_annotation/cache/`, so a re-run only redoes failed drugs.
- **STRING API (`string-db.org`) calls slow/fail for large gene sets**:
  normal for 1000+ genes (the `blue` module, for example) — just slower,
  not broken. Each script uses `httr::timeout()` generously for this.

## Which script needs which packages

| Script(s) | Packages |
|---|---|
| `00_functions.R`, `GSE*_run.R` | DESeq2, GEOquery, pheatmap, ggplot2, ggrepel, RColorBrewer, UpSetR, matrixStats, apeglm |
| `GSE229020_miRNA_run.R` | limma, ggplot2, ggrepel, pheatmap, UpSetR, WGCNA |
| `WGCNA_TBmodules_step1.R`, `WGCNA_remerge_test.R` | WGCNA |
| `WGCNA_TBmodules_step2_enrichment.R` | clusterProfiler, org.Hs.eg.db, ReactomePA, dplyr |
| `WGCNA_DEG_integration_step3.R` | none beyond base R |
| `WGCNA_STRING_step4.R`, `WGCNA_combined_network_step5.R` | httr, igraph, clusterProfiler + org.Hs.eg.db (ENSEMBL→SYMBOL) |
| `WGCNA_iLINCS_signature_step6.R`, `WGCNA_iLINCS_query_step7.R` | drugfindR, dplyr, clusterProfiler + org.Hs.eg.db |
| `WGCNA_drug_targets_step8.R` | httr, jsonlite |
| `WGCNA_network_algorithms_step9.R` | igraph, org.Hs.eg.db, AnnotationDbi |
| `WGCNA_final_tables_step10.R` | org.Hs.eg.db, AnnotationDbi |
| `WGCNA_GSEA_step11.R` | clusterProfiler, org.Hs.eg.db, msigdbr + msigdbdf |
| `WGCNA_Enrichr_drugs_step12.R` | enrichR |
| `WGCNA_singlecell_step13.R` | Matrix (+ the SCP1749 files in `external_data/SCP1749/`, see README) |
| `WGCNA_compare_runs_step14.R` | none beyond base R |
