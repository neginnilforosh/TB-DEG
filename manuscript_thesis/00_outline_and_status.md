# Manuscript / thesis outline and status (checkpoint 2026-10-14)

Working title (thesis): *Computational network medicine approaches for tuberculosis treatment in patients diagnosed with drug-resistant, latent, and drug-susceptible tuberculosis.*

Status legend:
- **VERIFIED**: number read from a result file at this checkpoint.
- **PRELIMINARY**: exists, but depends on pending decisions.
- **BLOCKED**: needs a decision or rerun.
- **NOT STARTED**

All numbers below were read from repository files on 2026-10-14. None were recomputed by inference.

## 1. Introduction — NOT STARTED
TB burden and comorbidities (DM, HIV); host-directed therapy and drug repurposing; network medicine. Citations must be added from the literature; none are drafted here.

## 2. Methods

### 2.1 Reference cohort (TB vs HC) — VERIFIED
- GSE114192 (Eckold et al. 2021, PMID 32533832), whole-blood RNA-seq, GPL18573, HTSeq counts.
- Healthy_Control n = 36, TB_Only n = 46 (82 samples). 13,939 genes after low-expression filtering.
- DESeq2 `~ group`, Wald test, apeglm-shrunk log2FC, BH-FDR. Significance: padj < 0.05 and |log2FC| ≥ 1.
- WGCNA (unsigned, minModuleSize 30, mergeCutHeight 0.25). Soft power reconstructed as 16 from the saved figure, so it is **unverified**; log it on any rerun.
- Module selection: BH-FDR < 0.05, then top 3 by |r|. **D3 pending.**
- STRING (score ≥ 0.4); centralities; Louvain; RWR (restart 0.5); diffusion.
- Network proximity: closest distance, 1,000 degree-matched permutations, seed 42.
- iLINCS via drugfindR; GSEA with sign(LFC) × −log10(p) ranking.
- Versions: `reproducibility/package_versions.csv` (R 4.5.0, DESeq2 1.50.2, WGCNA 1.74, limma 3.66.0, clusterProfiler 4.18.4, drugfindR 1.0.0, igraph 2.3.4, msigdbr 26.1.1). apeglm, glmnet and data.table versions are not recorded.

### 2.2 TB-DM — PRELIMINARY / BLOCKED
Design record: `GSE114192/TB_DM/00_*`. Pending D1 (site), D2 (sex), D3, D4 and D13. Contrasts C1–C5 are in `00_contrast_plan.csv`.

### 2.3 TB-HIV (GSE50834) — NOT STARTED
Readiness audit done: Illumina HT-12 v4, PBMC, 23 HIV vs 21 HIV/TB arrays, 43 patients. Group is confounded with processing period. **D10 pending.**

### 2.4 Validation layers — BLOCKED (label errors)
- GSE161829 (CD4 memory T cells): 7 of 85 samples mislabelled.
- GSE222001 (NK cells, 5 per group): Healthy and Latent swapped.
- GSE99374 (CD8 T cells): OK.
- GSE165709 (RNA + ATAC, HC/PLWH): supporting layer.
- ST000974: case definition unverified.
- CTD: annotation only.

## 3. Results

### 3.1 Reference: TB vs HC — VERIFIED
- **DEGs:** 300 significant (284 up, 16 down in TB; all with padj ≤ 0.0053).
- **WGCNA:** 18 merged modules (52 before merging; 6,320 grey genes). Selected modules:

| Module | Genes | r with TB | FDR |
|---|---|---|---|
| green | 356 | +0.73 | 2.0e-13 |
| purple | 213 | −0.70 | 1.3e-12 |
| blue | 1,377 | −0.63 | 1.1e-9 |

- **DEG–WGCNA overlap:** 188 genes. MM and GS reproduced exactly from the workspace.
- **GSEA Hallmark:** 28 sets with padj < 0.05 (20 positive, 8 negative). Top sets:

| Hallmark set | NES | padj |
|---|---|---|
| IFN-gamma response | +2.72 | 1e-37 |
| IFN-alpha response | +2.69 | 1e-24 |
| IL6-JAK-STAT3 | +2.10 | — |
| Inflammatory response | +2.05 | — |

- **PPI:** combined STRING network of 1,626 nodes and 14,726 edge rows.
- **Consensus genes:** 47 High, 373 Medium, 1,526 Low (1,088 of the Low genes meet 0 criteria). High-set Jaccard ranges 0.45–0.47 for the 10%/30% cut-offs, so report the High set as threshold-dependent.
- **Network proximity:** negative result. 395 drugs tested, 4 nominal p < 0.05, 0 with FDR < 0.05 (minimum FDR 1.0).
- **Consensus drugs:** 8 High under the current rule, but 6 of the 8 depend on proximity votes. Under D5 options B2 or C, only ipriflavone and gossypol remain High. Drugs are computational hypotheses only.

### 3.2 External validation — BLOCKED until D11 reruns
Existing direction-concordance numbers must not be reported, because the GSE161829 and GSE222001 labels are wrong.

### 3.3 GSE165709 — VERIFIED descriptive
| Layer | HC: INF vs NEG | PLWH: INF vs NEG | Formal HC/PLWH difference |
|---|---|---|---|
| RNA genes (padj < 0.05) | 2,259 | 669 | 0 |
| ATAC regions (padj < 0.05) | 238 | 0 | 0 |

44 of 47 High genes respond in HC RNA. Five genes have RNA+ATAC support, in HC only.

### 3.4 TB-DM — PRELIMINARY (not for reporting)
Current preliminary numbers:
- 295 DEGs TB-DM vs DM; 31 DEGs TB-DM vs TB.
- 3 selected modules, 1,563 genes.
- 38 High genes, 12 High drugs (0 proximity results with FDR < 0.05).

### 3.5 TB-HIV, cross-context integration — NOT STARTED (by instruction)

## 4. Discussion — NOT STARTED
Points to cover:
- the interferon signature
- cell-composition versus regulation
- the negative proximity finding
- site and batch confounding
- cell-type mismatch of the validation sets
- limitations

## 5. Figures
See `reproducibility/checkpoint_2026-10-14/FIGURE_INVENTORY.csv`.
- No PNG contains an embedded DPI tag.
- Scripts specify 600 dpi, except GSE229020, which uses 300.
- 31 raster figures lack a PDF/SVG companion.
