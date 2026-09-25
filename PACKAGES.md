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
  "matrixStats", "dplyr", "httr", "igraph"
))

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
- **STRING API (`string-db.org`) calls slow/fail for large gene sets**:
  normal for 1000+ genes (the `blue` module, for example) — just slower,
  not broken. Each script uses `httr::timeout()` generously for this.

## Which script needs which packages

| Script(s) | Packages |
|---|---|
| `00_functions.R`, `GSE*_run.R` | DESeq2, GEOquery, pheatmap, ggplot2, ggrepel, RColorBrewer, UpSetR, matrixStats, apeglm |
| `GSE229020_miRNA_run.R` | limma, ggplot2, ggrepel, pheatmap, UpSetR, WGCNA |
| `WGCNA_GSE114192_run.R`, `WGCNA_TBmodules_step1.R`, `WGCNA_remerge_test.R` | WGCNA |
| `WGCNA_TBmodules_step2_enrichment.R` | clusterProfiler, org.Hs.eg.db, ReactomePA, dplyr, (optional) msigdbr + msigdbdf |
| `WGCNA_DEG_integration_step3.R` | none beyond base R |
| `WGCNA_STRING_step4.R`, `WGCNA_combined_network_step5.R` | httr, igraph, (optional) clusterProfiler + org.Hs.eg.db for ENSEMBL→SYMBOL mapping |
