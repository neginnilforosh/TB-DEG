# TB-DEG project status and action checklist

Audit date: 14 Oct 2026 (Europe/Helsinki). Repository: `~/Desktop/TB-DEG`, branch `main`, HEAD `f2b3a52`.
Audit type: read-only. No scientific analysis was run to produce this checklist.

**Status labels.** `[x] DONE — VERIFIED` · `[~] PARTIALLY COMPLETE` · `[ ] NOT STARTED` ·
`[!] BLOCKED — DECISION REQUIRED` · `[?] NOT VERIFIED`

**Priority.** P0 = integrity risk, wrong sample labels, or threat to reported results.
P1 = decision or correction blocking a core analysis. P2 = reproducibility, documentation, figure quality.
P3 = optional or later thesis work.

**Owner.** Me = Negin. Ratul = supervisor decision. Later step = an analysis stage that cannot start yet.

**Evidence basis.** A task is `[x]` only when a file in the repository supports it now. Checks run earlier on
14 Oct whose output was not saved as a file are marked `[?]`, even if they passed at the time.

---

## 1. Current repository state and recent cleanup

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| R-01 | [x] DONE — VERIFIED | Record Git state | `git rev-parse HEAD` = `f2b3a52`; local `origin/main` ref = `f2b3a52` (not fetched); 0 staged; 2 modified tracked files (`README.md`, `scripts/wgcna_network/external_validation_step21.R`); 54 untracked files in 26 status entries | Nothing | — | P2 | — | — |
| R-02 | [x] DONE — VERIFIED | Duplicate TB-DM GSEA files removed | `GSE114192/TB_DM/17_gsea/` absent; 21 authoritative TB-DM GSEA files present and unchanged vs HEAD in `GSE114192/17_gsea/` | Nothing | — | P2 | — | — |
| R-03 | [x] DONE — VERIFIED | Finder duplicates removed | `GSE114192/01_raw_counts/GSE114192/untarred/`: 249 files, all canonical `GSM*_RSEQ*.txt.gz`, 0 `* 2`/`* 3` copies | Nothing | — | P2 | — | — |
| R-04 | [x] DONE — VERIFIED | Cleanup docs contain the intended information | `README.md` has "Current status", "Open decisions (pending Ratul)" (D1–D13) and the merged items in "Known limitations"; `GSE114192/TB_DM/00_dataset_decisions.txt` §11–§12 | Older README sections still carry pre-correction statements (see REF-15) | — | P2 | — | — |
| R-05 | [x] DONE — VERIFIED | Reference and other tracked scientific outputs unchanged | `git diff HEAD` is empty for `GSE114192/`, `GSE161829/`, `GSE222001/`, `GSE99374/`, `GSE229020/`, `comparisons/`, `cross_dataset_comparison/`, `external_results/`, `reproducibility/`; no file changed since the post-cleanup checksum snapshot | Nothing | — | P0 | — | — |
| R-06 | [x] DONE — VERIFIED | Sample-label audit files still exist | `GSE161829/00_checkpoint_2026-10-14/GSE161829_count_column_to_GSM_audit.csv` (106 rows); `GSE222001/00_checkpoint_2026-10-14/GSE222001_sample_label_audit.csv` (15 rows) | Keep until corrected reruns are verified | D11 | P0 | — | — |
| R-07 | [x] DONE — VERIFIED | Evidence for open decisions still exists | `GSE114192/TB_DM/00_dataset_decisions.txt`, `00_contrast_plan.csv`, `00_sample_manifest.csv`, `00_design/` (6 files); module-selection CSVs (2); `*/19_consensus/diagnostics_2026-10-14/` (5 + 5); `reproducibility/checkpoint_2026-10-14/` (6 files); `GSE50834/` (2) | Keep until D3–D5 and D10 are resolved | D3–D5, D10 | P1 | — | — |
| R-08 | [~] PARTIALLY COMPLETE | Remove temporary or redundant files | Remaining: `scripts/checkpoint/__pycache__/rds_min.cpython-310.pyc` (cache); `reproducibility/checkpoint_2026-10-14/RESULTS_MANIFEST.csv` lists files deleted in the cleanup (stale snapshot); 4 dated folder names (`00_checkpoint_2026-10-14`, `checkpoint_2026-10-14`, `diagnostics_2026-10-14`); ignored local files (`.DS_Store` ×8, `.RData`, `.Rhistory`, `.Rapp.history`) | Decide on the cache, the stale manifest and dated folder names once D3–D5/D11 are closed | Rules 10 and 14 of the cleanup (keep until decisions) | P2 | Revisit after D11 and D3–D5 | Me |
| R-09 | [!] BLOCKED — DECISION REQUIRED | Uncommitted local changes | 2 modified + 54 untracked files exist only locally | Decide what to commit; whether to keep the comment/label edit in `external_validation_step21.R` | Your approval | P2 | Review `git diff` and the untracked list | Me |
| R-10 | [x] DONE — VERIFIED | Tracked raw files and KEGG caches kept as decided | `GSE161829/01_raw_counts/GSE161829/GSE161829_raw_counts.txt`, `GSE99374/01_raw_counts/GSE99374/GSE99374_RawCounts_CD8.xlsx`, `GSE114192/22_metabolite_gene_bridge/kegg_cache.rds`, `GSE114192/TB_DM/22_metabolite_gene_bridge/kegg_cache.rds` are still tracked | None now. `GSE161829_raw_counts.txt` holds the original `SampleN_CD4m` column names needed for the D11 fix | — | P2 | — | — |
| R-11 | [?] NOT VERIFIED | Evidence lost in cleanup | Full `DECISIONS_FOR_RATUL.md` and `DIAGNOSTICS_NOTES.md` were deleted after merging; copies exist only as chat attachments, not on disk | Save the copies somewhere you control if you need the full proposal text | — | P2 | Save the two chat attachments outside the repo | Me |

---

## 2. Completed reference analysis (GSE114192, TB vs healthy control)

Scientific status: **computationally complete; interpretation partly unresolved** (D3–D5, validation layers).
Do not rerun or overwrite.

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| REF-01 | [x] DONE — VERIFIED | Sample selection and metadata | `GSE114192/02_metadata/GSE114192_sample_metadata.csv` (HC 36, TB_Only 46); count-matrix columns are GSM IDs (`GSE114192/01_raw_counts/GSE114192_raw_count_matrix.csv`, 82 columns, no duplicates) | Nothing | — | P1 | — | — |
| REF-02 | [x] DONE — VERIFIED | DE: filtering, normalisation, DESeq2, shrinkage | `scripts/deg_pipeline/GSE114192_run.R`; `03_filtered_counts/`, `04_normalized/`, `06_deg_results/GSE114192_HealthyControl_vs_TBOnly_DEG.csv` (13,939 genes, apeglm-shrunk); 300 significant (284 up, 16 down) | Nothing | — | P1 | — | — |
| REF-03 | [~] PARTIALLY COMPLETE | DEG tables carry ID, symbol, effect, significance, direction, contrast | Original tables: Ensembl ID, log2FC, p, padj only. Annotated copies add symbol, direction, comparison, LFC type: `06_deg_results/..._DEG_annotated.csv`, `07_significant_degs/..._sigDEGs_annotated.csv` (generated by `scripts/checkpoint/annotate_deg_tables.py`) | 421 of 13,939 genes have no symbol; the Wald statistic is absent from shrunk tables; the R pipeline does not write these columns itself | — | P2 | Build the columns into the R pipeline on the next approved rerun | Later step |
| REF-04 | [~] PARTIALLY COMPLETE | WGCNA: module–trait, merging, eigengenes | `10_wgcna_results/GSE114192_ModuleTrait_Correlation.csv`, `..._Gene_Module_Assignment.csv`, `..._WGCNA_workspace.rds`; eigengenes and merging summary exported from the workspace (52 → 18 modules) | Soft power not saved; reconstructed ≈16 from the figure (`reproducibility/checkpoint_2026-10-14/WGCNA_soft_threshold_reconstructed_from_figure.csv`, approximate) | — | P2 | Log the power on any future run | Later step |
| REF-05 | [!] BLOCKED — DECISION REQUIRED | WGCNA module selection | Rule: FDR < 0.05, then top 3 by \|r\| (`WGCNA_TBmodules_step1.R`). 12 modules pass FDR; 3 selected (green, purple, blue). Comparison: `10_wgcna_results/checkpoint_2026-10-14/module_selection_comparison.csv` | Choose the rule | D3 | P1 | Present the comparison table to Ratul | Ratul |
| REF-06 | [?] NOT VERIFIED | MM/GS values correct | `10_wgcna_results/GSE114192_TBmodules_MM_GS.csv` (1,946 genes). Recomputed and matched on 14 Oct, but that check was not saved as a script | Re-check with a saved script if needed for the thesis | — | P2 | Optional: save an MM/GS check script | Me |
| REF-07 | [x] DONE — VERIFIED | GSEA | `17_gsea/GSE114192_HealthyControl_vs_TBOnly_GSEA_{HALLMARK,REACTOME,GOBP,CELLTYPE}.csv`; Hallmark: 28 sets padj < 0.05 (IFN-γ NES +2.72, IFN-α +2.69) | Nothing | — | P2 | — | — |
| REF-08 | [~] PARTIALLY COMPLETE | STRING PPI, centralities, communities, RWR, diffusion | `12_string_ppi/` (combined network 1,626 nodes), `15_network_algorithms/`; seeds recorded in `reproducibility/PARAMETERS.md` | STRING version not recorded (queried live; `STRINGdb not installed`) | — | P2 | Record STRING version/date on next run | Later step |
| REF-09 | [!] BLOCKED — DECISION REQUIRED | Consensus gene ranking | `19_consensus/GSE114192_GenePriority_all_modules.csv` (all required columns present; 47 High / 373 Medium / 1,526 Low). Diagnostics: 1,088 "Low" genes meet 0 criteria; High-set Jaccard 0.45–0.47 at 10%/30% cut-offs (`19_consensus/diagnostics_2026-10-14/`) | Tier scheme decision | D4 | P1 | Ratul decides; ranking recomputation would need your approval | Ratul |
| REF-10 | [!] BLOCKED — DECISION REQUIRED | Drug ranking and network proximity | `19_consensus/GSE114192_DrugPriority.csv`; proximity with BH-FDR: `15_network_algorithms/GSE114192_NetworkProximity_withFDR.csv` (395 drugs, 4 nominal p < 0.05, 0 FDR < 0.05). 6 of 8 High drugs rely on proximity votes | Decide the role of proximity; README "not used for ranking" text contradicts step 15 | D5 | P1 | Ratul decides | Ratul |
| REF-11 | [!] BLOCKED — DECISION REQUIRED | External validation (step 21) | `23_external_validation/` uses GSE161829 and GSE222001, both with label errors | Corrected reruns | D11 | P0 | See VAL-01/VAL-02 | Ratul → Me |
| REF-12 | [!] BLOCKED — DECISION REQUIRED | ML validation (step 22) | `24_ml_validation/GSE114192_ml_validation_AUC.csv` | In-sample discovery AUCs; external rows affected by label errors | D11, D12 | P0 | See ML-01–ML-03 | Ratul → Me |
| REF-13 | [?] NOT VERIFIED | Metabolomics (ST000974) | `21_metabolomics/`, `22_metabolite_gene_bridge/`; M06 subset, 26 case vs 85 control, 9/699 metabolites padj < 0.05 | Meaning of "case" not verified; no run log showing whether REGION/GENDER adjustment ran | D7 | P1 | See SUP-01 | Me |
| REF-14 | [~] PARTIALLY COMPLETE | Single-cell support (SCP1749) | `18_singlecell/`; macaque lung granuloma, 2 animals | Supporting evidence only; stated in README | — | P3 | — | — |
| REF-15 | [~] PARTIALLY COMPLETE | Limitations and negative results reported honestly | README "Known limitations" now covers label errors, ML leakage, proximity FDR, soft power, GSE165709, ST000974, CTD, DPI | README "Current results" and "Replication" sections still state pre-correction findings; status block flags this but the text itself is unchanged | D11 | P2 | Rewrite those sections after the D11 reruns | Me |
| REF-16 | [x] DONE — VERIFIED | Reference outputs preserved | See R-05 | Do not rerun or overwrite without explicit approval | — | P0 | — | — |
| REF-17 | [~] PARTIALLY COMPLETE | Latent positional-mapping risk | `GSE114192_run.R` (lines 113–121) and `GSE114192_comorbidity_run.R` (≈ line 92) fall back to positional renaming if names fail to match. Not triggered now | Replace the fallback with an error (script change) | Your approval | P2 | Include in the D11 fix | Me |

---

## 3. TB-DM (GSE114192/TB_DM)

Scientific status: **preliminary**. Do not rerun or finalise until D1–D5 and D13 are decided.

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| DM-01 | [x] DONE — VERIFIED | Four-group selection and exclusions | `GSE114192/TB_DM/02_metadata/GSE114192_TB_DM_sample_metadata.csv` (TB_DM 61, DM_only 52, TB_Only 46, HC 36 = 195); `00_sample_manifest.csv` (249 rows, inclusion reasons; IH 10 and TB-IH 44 excluded) | Nothing | — | P1 | — | — |
| DM-02 | [x] DONE — VERIFIED | Manifest, contrast plan, decision record present | `00_sample_manifest.csv`, `00_contrast_plan.csv` (C1–C5), `00_dataset_decisions.txt` | Nothing | — | P1 | — | — |
| DM-03 | [!] BLOCKED — DECISION REQUIRED | Site confounding and identifiability | `00_design/site_by_group.csv`, `design_identifiability.csv`, `site_group_association.txt` (Cramér's V 0.34; HC/DM only at SA+RO) | Choose the site design | D1 | P1 | Ratul decides | Ratul |
| DM-04 | [!] BLOCKED — DECISION REQUIRED | Sex covariate | Not in GEO; options in `00_dataset_decisions.txt` §4 | Choose an option | D2 | P1 | Ratul decides | Ratul |
| DM-05 | [~] PARTIALLY COMPLETE | Primary: TB-DM vs DM-only | `06_deg_results/GSE114192_DM_only_vs_TB_DM_DEG.csv`; 295 significant (`GSE114192_TB_DM_comparison_summary.csv`); `~ group`, apeglm-shrunk | Site-confounded; preliminary | D1, D2, D13 | P1 | Rerun after decisions | Later step |
| DM-06 | [~] PARTIALLY COMPLETE | Secondary: TB-DM vs TB-only | `06_deg_results/GSE114192_TB_Only_vs_TB_DM_DEG.csv`; 31 significant; unshrunk | Preliminary | D1, D2, D13 | P1 | Rerun after decisions | Later step |
| DM-07 | [~] PARTIALLY COMPLETE | Interaction model (exploratory) | `06_deg_results/GSE114192_interaction_{interaction,main_TB,main_DM}.csv`; `17_gsea_interaction/` (4 files); labelled exploratory in `00_dataset_decisions.txt` §8 | Only SA+RO contain all four groups; rerun with site handling | D1, D13 | P1 | Rerun after decisions | Later step |
| DM-08 | [!] BLOCKED — DECISION REQUIRED | Consistent shrinkage across contrasts | Only `DM_only_vs_TB_DM` is apeglm-shrunk; others unshrunk (`00_dataset_decisions.txt` §9) | Apply one rule in the rerun | D13 | P1 | Include in rerun plan | Later step |
| DM-09 | [~] PARTIALLY COMPLETE | WGCNA soft power, module rule, eigengenes | Eigengenes and merging summary exported (46 → 29 modules); soft power not saved; R² never reaches 0.85, fallback power 6 most likely used (reconstructed, approximate) | Log power; decide module rule (18 modules pass FDR, top 3 kept; brown r 0.49 excluded) | D3, D13 | P1 | Ratul decides D3 | Ratul |
| DM-10 | [~] PARTIALLY COMPLETE | GSEA | 4 TB-DM comparisons in `GSE114192/17_gsea/` (authoritative location) | Preliminary | D13 | P2 | — | Later step |
| DM-11 | [!] BLOCKED — DECISION REQUIRED | Network, consensus genes, drugs, proximity | `12_string_ppi/`, `15_network_algorithms/` (+ `_withFDR.csv`: 432 drugs, 3 nominal, 0 FDR), `19_consensus/` (38 High genes, 12 High drugs; 10 of 12 rely on proximity) | Preliminary; rules undecided | D3, D4, D5 | P1 | Ratul decides | Ratul |
| DM-12 | [~] PARTIALLY COMPLETE | Downstream TB-DM layers | `TB_DM/23_external_validation/`, `24_ml_validation/`, `25_ctd_annotation/`, `22_metabolite_gene_bridge/`, `figures/`, `figures_manuscript/` | Built on preliminary results and on mislabelled validation sets | D11, D13 | P1 | Regenerate after rerun | Later step |
| DM-13 | [x] DONE — VERIFIED | Preliminary status labelled; decisions linked | `00_dataset_decisions.txt` §11–§12; README status table | Nothing | — | P2 | — | — |
| DM-14 | [!] BLOCKED — DECISION REQUIRED | Approved TB-DM rerun | Plan: `00_contrast_plan.csv` | Needs D1–D5 and an R environment | D1–D5, D13 | P1 | Ask Ratul for D1–D5 | Ratul |

---

## 4. TB-HIV (GSE50834)

Scientific status: **not started**. Do not begin DE before D10.

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| HIV-01 | [x] DONE — VERIFIED | Platform and material | `GSE50834/00_dataset_audit.txt`: GPL10558 (Illumina HumanHT-12 V4), PBMC, 44 arrays (HIV only 23, HIV/TB 21) | Nothing | — | P1 | — | — |
| HIV-02 | [x] DONE — VERIFIED | Patient / repeat-array relationships | `GSE50834/00_sample_table_from_GEO.csv`: 43 patients; patient 2 has two arrays (GSM1230903/04) | Decide how to handle the repeat | D10 | P1 | Include in D10 | Ratul |
| HIV-03 | [~] PARTIALLY COMPLETE | Name-based sample matching | Barcode-to-GSM map from GEO `!Sample_description` is in the sample table | No run script yet | D10 | P1 | Use this mapping when the script is written | Later step |
| HIV-04 | [!] BLOCKED — DECISION REQUIRED | Group vs processing-batch confounding | 34 of 44 arrays: scan period predicts group; only 2007–08 has both groups (5 vs 5) | Decide whether and how to proceed | D10 | P1 | Ratul decides | Ratul |
| HIV-05 | [ ] NOT STARTED | Search for an alternative or additional TB-HIV dataset | None | Search GEO/ArrayExpress for HIV vs HIV/TB with balanced processing | D10 option C | P1 | Do the search if Ratul chooses it | Me |
| HIV-06 | [~] PARTIALLY COMPLETE | Platform-appropriate method | Proposal in `00_dataset_audit.txt`: limma neqc on non-normalised data, detection filter, probe→gene rule | Raw files (`GSE50834_non-normalized.txt.gz`, `GSE50834_RAW.tar`) not downloaded | D10 | P1 | Download after D10 | Me |
| HIV-07 | [ ] NOT STARTED | Differential expression (TB-HIV vs HIV) | None | Full pipeline | D10 | P1 | Wait | Later step |
| HIV-08 | [ ] NOT STARTED | WGCNA, networks, ranking, drugs for TB-HIV | None | Full pipeline | D10, HIV-07, D3–D5 | P1 | Wait | Later step |

---

## 5. External validation and replication

| ID | Status | Dataset / task | Material, platform, samples | Label verification | Dependents | What remains | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|---|
| VAL-01 | [!] BLOCKED — DECISION REQUIRED | GSE161829 | Memory CD4+ T cells (cryopreserved PBMC), SmartSeq2, GPL24676; analysed 85 (TBneg 20, LTBI 40, ATB 25); 21 post-treatment arrays excluded | **Error:** `GSE161829_run.R` lines 49–51 rename columns by position; 28 columns mismatched, 7 wrong group labels (`GSE161829_count_column_to_GSM_audit.csv`) | `GSE161829/03`–`09`, `17_gsea/`, `figures/`; `GSE161829/TBneg_vs_ATB/`; `comparisons/`; `cross_dataset_comparison/`; `23`/`24` in both GSE114192 runs | Name-based fix and rerun | P0 | Ask Ratul to approve D11 | Ratul → Me |
| VAL-02 | [!] BLOCKED — DECISION REQUIRED | GSE222001 | NK cells, HiSeq X Ten; 15 (5 per group) | **Error:** `GSE222001_run.R` lines 77–80 assign labels by position; Healthy and Latent swapped (`GSE222001_sample_label_audit.csv`) | `GSE222001/03`–`09`, `17_gsea/`, `figures/`; `cross_dataset_comparison/`; `23`/`24` in both GSE114192 runs | Name-based fix and rerun | P0 | Ask Ratul to approve D11 | Ratul → Me |
| VAL-03 | [~] PARTIALLY COMPLETE | GSE99374 | Sorted CD8+ T cells; 40 (HC 19, LTBI 21); 4 significant DEGs; outlier sensitivity in `GSE99374/sensitivity_no_outliers/` | Labels match author sample-name suffixes (40/40); not cross-checked against GEO characteristics | Step 21 LTBI rows | GEO cross-check (optional) | P2 | Optional check | Me |
| VAL-04 | [?] NOT VERIFIED | GSE229020 | NanoString nCounter miRNA panel; 20 (HC 4, Latent 6, DrugSusceptible 6, DrugResistant 4); 1/0/0 significant miRNAs | Labels parsed from column-name prefixes in `GSE229020_miRNA_run.R`; not checked against GEO; biological material not verified | Not used in steps 21–22 | Verify material and labels against GEO | P2 | Check GEO record | Me |
| VAL-05 | [?] NOT VERIFIED | GSE165709 (RNA + ATAC) | Results in `external_results/GSE165709/`; RNA HC 16+16, PLWH 20+20; ATAC HC 14+14, PLWH 17+17 | Sample labels not re-checked against GEO in this audit | Integration step 17, ATAC gene-set step 18 | Verify labels; compare methods with supervisor's `external_data/GSE165709/ATAC1.Rmd` (present, git-ignored) | P2 | Read `ATAC1.Rmd` against step 16 | Me |
| VAL-06 | [!] BLOCKED — DECISION REQUIRED | Replication (step 14) | `comparisons/GSE114192_TB_vs_GSE161829_ATB/` | Depends on GSE161829 labels | VAL-01 | Rerun after D11 | P0 | Wait for D11 | Later step |
| VAL-07 | [!] BLOCKED — DECISION REQUIRED | Cross-dataset UpSet | `cross_dataset_comparison/` (`05_cross_dataset_upset.R`) | Uses GSE161829 and GSE222001 significant lists | VAL-01, VAL-02 | Rerun after D11 | P0 | Wait for D11 | Later step |
| VAL-08 | [?] NOT VERIFIED | GSE161829 significant-gene count discrepancy | `GSE161829/07_significant_degs/GSE161829_TBneg_vs_ATB_sigDEGs.csv` has 50 genes; `comparisons/.../SUMMARY.txt` reports 39 | Likely due to the 12,074 shared-gene universe in step 14; not confirmed | VAL-01 | Recheck after rerun | P2 | — | Later step |

---

## 6. Machine learning and supporting evidence

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| ML-01 | [x] DONE — VERIFIED | Leakage identified | `scripts/wgcna_network/ml_validation_step22.R`: features from the High tier chosen on all 82 samples; lambda.min chosen on the reported folds. AUCs 0.975 / 0.894 are in-sample (documented in README) | Report as apparent performance only | — | P0 | — | — |
| ML-02 | [ ] NOT STARTED | Nested cross-validation | No nested-CV code anywhere in `scripts/` | Implement after approval | D12 | P1 | Ratul decides D12 | Ratul |
| ML-03 | [!] BLOCKED — DECISION REQUIRED | Out-of-sample (external) AUCs | `24_ml_validation/` rows for GSE161829 and GSE222001 | Invalid until labels are corrected | D11 | P0 | Wait for D11 | Later step |
| SUP-01 | [?] NOT VERIFIED | ST000974 case definition | `21_metabolomics/ST000974_sample_factors.csv` (GROUP case/control; BL/M06/M18; AHRI/MRC/SUN) | Read the Metabolomics Workbench study page by hand (automated access is blocked); record what "case" means | D7 | P1 | Open the ST000974 page | Me |
| SUP-02 | [?] NOT VERIFIED | ST000974 covariate adjustment actually applied | `metabolomics_step19.R` adjusts for REGION/GENDER only if limma is installed; no run log saved | Confirm from a log or rerun (approval needed) | SUP-01 | P2 | Check after SUP-01 | Me |
| SUP-03 | [x] DONE — VERIFIED | GSE165709 described by its actual evidence | README: HC 2,259 genes / 238 regions; PLWH 669 / 0; response difference 0 / 0; 5 RNA+ATAC genes in HC only | Nothing | — | P2 | — | — |
| SUP-04 | [x] DONE — VERIFIED | CTD used as annotation only | `25_ctd_annotation/`; README note | CTD release not recorded | — | P2 | Record version on next run | Later step |
| SUP-05 | [~] PARTIALLY COMPLETE | Metabolite–gene bridge | `22_metabolite_gene_bridge/` (both runs) | Depends on ST000974 case definition | SUP-01 | P2 | — | Later step |
| SUP-06 | [ ] NOT STARTED | Mechanism-of-action / latent-space drug prediction | No script found | Later-stage work | D8 | P3 | Wait | Later step |
| SUP-07 | [ ] NOT STARTED | H3K27ac ChIP-seq layer | None | Later-stage work | D7 | P3 | Wait | Later step |

---

## 7. Cross-context integration

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| INT-01 | [!] BLOCKED — DECISION REQUIRED | Independent analyses complete | Reference complete; TB-DM preliminary; TB-HIV not started | TB-DM rerun, TB-HIV pipeline | D1–D5, D10, D13 | P1 | Wait | Later step |
| INT-02 | [x] DONE — VERIFIED | Raw count matrices not combined across studies | `WGCNA_compare_runs_step14.R` merges result tables by gene ID only; `05_cross_dataset_upset.R` uses significant-gene lists | Keep this rule | — | P0 | — | — |
| INT-03 | [ ] NOT STARTED | Cross-context DEG table (TB, TB-DM, TB-HIV) with direction pattern | None | Requires all three contexts | INT-01 | P1 | Wait | Later step |
| INT-04 | [ ] NOT STARTED | Core/shared network and context-specific rewiring | None | Requires all three networks | INT-01 | P1 | Wait | Later step |
| INT-05 | [~] PARTIALLY COMPLETE | Preserve direction; separate untested from non-significant | Annotated DEG tables have `direction` and `padj_status`; `cross_dataset_comparison/` direction handling not verified | Apply in INT-03 | INT-03 | P2 | — | Later step |

---

## 8. Manuscript, thesis, and final deliverables

| ID | Status | Task | Evidence | What remains | Dependencies / blockers | Pri | Next action | Owner |
|---|---|---|---|---|---|---|---|---|
| MS-01 | [~] PARTIALLY COMPLETE | Outline | `manuscript_thesis/00_outline_and_status.md` (draft; under your review; numbers taken from repository files on 14 Oct) | Your review; decide whether it belongs in the repository | D9 | P2 | Review the outline | Me |
| MS-02 | [x] DONE — VERIFIED | Final gene and drug tables (reference) | `19_consensus/GSE114192_GenePriority_all_modules.csv` (Gene, Module, log2FC, FDR, MM, GS, PPI degree, betweenness, RWR, diffusion, priority); `19_consensus/GSE114192_DrugPriority.csv` (reversal score, proximity Z/P, network targets, modules, enrichment P, high-priority targets, priority) | Tier content depends on D4/D5 | D4, D5 | P1 | — | — |
| MS-03 | [?] NOT VERIFIED | Figure resolution (600 dpi) and formats | 92 PNG and 94 PDF under `GSE*/`, `comparisons/`, `cross_dataset_comparison/`; no PNG has an embedded DPI tag; no TIFF; 31 raster figures lack a PDF/SVG (`reproducibility/checkpoint_2026-10-14/FIGURE_INVENTORY.csv`); scripts set 600 dpi except `GSE229020_miRNA_run.R` (300) | Export TIFF/PDF with DPI tags from the scripts on the next approved run | Approved reruns | P2 | Plan the export settings | Later step |
| MS-04 | [~] PARTIALLY COMPLETE | Figures traceable to scripts | `scripts/wgcna_network/figures_step24.R` → `figures_manuscript/` (41 files per run); per-step figures in `figures/` | TB-DM and validation figures will change after reruns | D11, D13 | P2 | — | Later step |
| MS-05 | [ ] NOT STARTED | Methods text | Outline only | Write after decisions | D1–D5 | P2 | Wait | Me |
| MS-06 | [ ] NOT STARTED | Introduction and discussion | None | Write | — | P3 | — | Me |
| MS-07 | [~] PARTIALLY COMPLETE | Limitations | README "Known limitations" | Move to thesis text | — | P2 | — | Me |
| MS-08 | [!] BLOCKED — DECISION REQUIRED | Manuscript vs thesis plan and scope | Ratul's 13 Oct email points to manuscript and thesis drafting | Decide plan and scope | D6, D9 | P1 | Discuss with Ratul | Ratul |
| MS-09 | [~] PARTIALLY COMPLETE | Software and database versions | `reproducibility/package_versions.csv` (30 packages), `sessionInfo.txt` (R 4.5.0), `PARAMETERS.md` (regenerated by `scripts/record_reproducibility.R`) | apeglm, glmnet, data.table versions missing; STRING, DGIdb, KEGG, CTD releases not recorded | — | P2 | Add on next `record_reproducibility.R` run | Later step |

---

## 9. Supervisor decisions required

| ID | Decision | Blocks | Evidence |
|---|---|---|---|
| D1 | TB-DM site handling | DM-03, DM-05–07, DM-14 | `GSE114192/TB_DM/00_dataset_decisions.txt` §3; `00_design/` |
| D2 | Sex covariate for TB-DM | DM-04–06, DM-14 | same file §4 |
| D3 | WGCNA module-selection rule | REF-05, DM-09, HIV-08 | `*/10_wgcna_results/checkpoint_2026-10-14/module_selection_comparison.csv` |
| D4 | Consensus gene tier scheme | REF-09, DM-11 | `*/19_consensus/diagnostics_2026-10-14/gene_*` |
| D5 | Role of network proximity in drug tiers | REF-10, DM-11 | `*/19_consensus/diagnostics_2026-10-14/drug_criteria_provenance_and_DRAFT_rules.csv` |
| D6 | Thesis scope (drug-resistant / latent TB) | MS-08 | — |
| D7 | Role of metabolomics, scRNA, GSE165709, ChIP-seq | REF-13, SUP-01, SUP-07 | README "Known limitations" |
| D8 | MoA / latent-space prediction | SUP-06 | — |
| D9 | Manuscript vs thesis plan | MS-01, MS-08 | `manuscript_thesis/00_outline_and_status.md` |
| D10 | GSE50834 suitability (batch confounding) | HIV-02–HIV-08, INT-01 | `GSE50834/00_dataset_audit.txt` |
| D11 | Name-based matching and reruns for GSE161829 and GSE222001 | REF-11, REF-12, VAL-01/02/06/07, ML-03 | the two audit CSVs |
| D12 | ML validation design (nested CV) | ML-02 | `ml_validation_step22.R` |
| D13 | Approve the TB-DM rerun | DM-05–DM-14 | `00_contrast_plan.csv` |

Additional question for D3–D5: should the chosen rules also be applied to the reference rankings?
That would recompute rankings only (not DE), but it changes a completed analysis and needs your approval.

---

## 10. Ordered next steps and roadmap

### Five most important immediate actions (in order)
1. **Send Ratul the D11 request** (P0): both label errors, their audit CSVs and the dependency list (VAL-01, VAL-02). Owner: Me.
2. **Ask Ratul for D1–D5 together** (P1), with the evidence files in section 9, because they unblock the TB-DM rerun and settle the ranking rules for all contexts. Owner: Me → Ratul.
3. **Ask Ratul for D10** (P1): GSE50834 confounding, with options A/B/C from `GSE50834/00_dataset_audit.txt`. Owner: Me → Ratul.
4. **Verify the ST000974 case definition by hand** (SUP-01, P1): no approval needed; it is a reading task. Owner: Me.
5. **Review `manuscript_thesis/00_outline_and_status.md` and the uncommitted changes** (MS-01, R-09, P2): decide what to keep and commit. Owner: Me.

**Single best next action:** send Ratul the D11 request (item 1). It is the only P0 item that needs his decision, and it unblocks the correction of every affected validation result. Not executed here.

### Tasks that must wait for Ratul
D1–D13 as listed in section 9; specifically DM-03–DM-14, REF-05, REF-09, REF-10, HIV-04–HIV-08, ML-02, MS-08.

### Tasks that need your approval before changing scientific outputs
- Name-based fixes to `GSE161829_run.R` and `GSE222001_run.R`, and reruns of their dependents (D11).
- Replacing the positional fallback in `GSE114192_run.R` / `GSE114192_comorbidity_run.R` (REF-17).
- Any recomputation of reference rankings under new D3–D5 rules.
- The TB-DM rerun (D13), written to a new run folder rather than over the preliminary outputs.
- Nested-CV ML (D12).

### Tasks that must not be repeated
- The TB vs healthy-control reference DE, WGCNA, GSEA, networks and drug analyses (REF-01–REF-08).
- The cleanup (R-02, R-03): already completed and verified.
- The GSE50834 metadata audit (HIV-01, HIV-02).

### Optional tasks that can wait (P3)
SUP-06 (MoA prediction), SUP-07 (ChIP-seq), MS-06 (introduction/discussion), REF-14 (single-cell refinements).

### Completed vs incomplete vs unverified
- **Completed and verified:** R-01–R-07, R-10; REF-01, REF-02, REF-07, REF-16; DM-01, DM-02, DM-13; HIV-01, HIV-02; ML-01; SUP-03, SUP-04; INT-02; MS-02.
- **Incomplete or blocked:** all `[~]`, `[ ]` and `[!]` items above, including the TB-DM pipeline, the whole TB-HIV pipeline, all validation results using GSE161829 or GSE222001, nested-CV ML, and cross-context integration.
- **Unverified claims:** REF-06 (MM/GS check not saved), REF-13/SUP-01/SUP-02 (ST000974), VAL-04 (GSE229020), VAL-05 (GSE165709 labels), VAL-08 (count discrepancy), MS-03 (figure DPI), R-11 (full text of deleted notes kept only outside the repo).
