#!/usr/bin/env python3
"""Additional file 1, one workbook with Supplementary Tables S1-S11 and the extended Tables 3-6.
Run after manuscript_tables.py (reads manuscript_tables/Main_tables.xlsx and Supplementary_tables.xlsx) and adds what those
two files do not hold: gene symbols for every gene in S1, ranking stability (S4), GO biological process and
cell-type GSEA (S5), the full iLINCS screen (S6) and software, versions and parameters (S11). A contents sheet
describes every sheet and names its source file. Every result value is read from a file; the parameter rows of
S11 restate the settings of the step scripts (soft-thresholding powers are read from the *_WGCNA_soft_threshold.csv tables).
Output: manuscript_tables/Additional_file_1.xlsx"""
import sys, re
from pathlib import Path
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment, PatternFill

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
for f in ("Main_tables.xlsx", "Supplementary_tables.xlsx"):
    if not (ROOT / "manuscript_tables" / f).exists(): sys.exit(f"Missing manuscript_tables/{f}: run scripts/manuscript_tables.py first.")
OUT = ROOT / "manuscript_tables" / "Additional_file_1.xlsx"
sup = pd.read_excel(ROOT / "manuscript_tables/Supplementary_tables.xlsx", sheet_name=None)
main = pd.read_excel(ROOT / "manuscript_tables/Main_tables.xlsx", sheet_name=None)
CTX = {"TB": ("GSE114192", "GSE114192", "HealthyControl_vs_TBOnly", "TB vs healthy controls (GSE114192)"),
       "TB-DM": ("GSE114192/TB_DM", "GSE114192", "DM_only_vs_TB_DM", "TB with diabetes vs diabetes alone (GSE114192)"),
       "TB-HIV": ("GSE37250/TB_HIV", "GSE37250", "LTBI_HIVpos_vs_ATB_HIVpos", "Active vs latent TB in people with HIV (GSE37250)")}
rd = lambda *p: pd.read_csv(ROOT.joinpath(*p))

# ---- gene symbol map: GSE114192 annotated DEG table + GSE37250 probe annotation
m1 = rd("GSE114192/06_deg_results/GSE114192_HealthyControl_vs_TBOnly_DEG_annotated.csv")[["gene", "gene_symbol"]].dropna()
m2 = rd("GSE37250/00_design/GSE37250_probe_annotation.csv")[["ensembl", "symbol"]].dropna().rename(columns={"ensembl": "gene", "symbol": "gene_symbol"})
sym = dict(zip(m2.gene, m2.gene_symbol)); sym.update(dict(zip(m1.gene, m1.gene_symbol)))

S = {}; DESC = []
def add(name, df, desc, src):
    S[name] = df; DESC.append((name, desc, len(df), src))

for k, (run, acc, deg, lab) in CTX.items():
    d = sup[f"S1_{k}"].drop(columns=["Gene"], errors="ignore").copy()
    d.insert(1, "Gene_symbol", d.gene.map(sym))
    d = d.rename(columns={"gene": "Ensembl"})
    if k == "TB-HIV":
        d = d.rename(columns={"baseMean": "Mean_log2_expression", "lfcSE": "SE"})
        note = "limma (~ country + group, BeadChip block); log2FoldChange = active minus latent TB; Mean_log2_expression = average log2 signal"
    else:
        note = "DESeq2 with apeglm-shrunken log2 fold changes; baseMean = mean normalised count"
    add(f"S1_{k}", d, f"Differential expression, {lab}. All tested genes. significant = padj < 0.05 and |log2FC| ≥ 1. {note}.",
        f"{run}/06_deg_results or {acc}/06_deg_results/{acc}_{deg}_DEG.csv")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S2_{k}_nodes", sup[f"S2_{k}_nodes"], f"Combined STRING network of the selected modules, nodes with module and centrality (betweenness and closeness normalised), {lab}.", f"{run}/12_string_ppi")
    e = sup[f"S2_{k}_edges"].rename(columns={"preferredName_A": "Protein_A", "preferredName_B": "Protein_B", "score": "STRING_combined_score"})
    add(f"S2_{k}_edges", e, f"Combined STRING network, edges (combined score ≥ 0.4), {lab}.", f"{run}/12_string_ppi/{acc}_COMBINED_STRING_network.csv")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S3_{k}", sup[f"S3_{k}"], f"Random walk with restart (restart 0.5) and network diffusion (α = 0.5) scores seeded with the disease genes, and Louvain community, {lab}.", f"{run}/15_network_algorithms")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S4_{k}", sup[f"S4_{k}"], f"Evidence-family scores of all module genes (Fam_* = within-run percentile, 0 to 1; empty = no data for that family), evidence score, rank, rank range across leave-one-family-out variants and tier, {lab}.", f"{run}/20_evidence_families/{acc}_Gene_evidence_families.csv")
for k, (run, acc, deg, lab) in CTX.items():
    v = rd(run, "20_evidence_families", f"{acc}_Gene_evidence_sensitivity.csv")
    add(f"S4_{k}_sensitivity", v, f"Stability of the gene ranking: Spearman correlation, top-50 overlap and High-tier Jaccard index of each variant against the main ranking, {lab}.", f"{run}/20_evidence_families/{acc}_Gene_evidence_sensitivity.csv")
for k, (run, acc, deg, lab) in CTX.items():
    parts = []
    for coll, suf in [("Hallmark", "HALLMARK"), ("Reactome", "REACTOME"), ("GO biological process", "GOBP"), ("Cell type (C8)", "CELLTYPE")]:
        f = ROOT / acc / "17_gsea" / f"{acc}_{deg}_GSEA_{suf}.csv"
        if f.exists():
            x = pd.read_csv(f).drop(columns=["Collection"], errors="ignore"); x.insert(0, "Collection", coll); parts.append(x)
    g = pd.concat(parts, ignore_index=True)
    keep = [c for c in ["Collection", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment"] if c in g.columns]
    add(f"S5_{k}", g[keep], f"Pre-ranked GSEA (ranking sign(log2FC) × −log10 p; gene sets of 15 to 500 genes; seed 42; BH adjustment), all tested sets of four collections, {lab}. core_enrichment = leading-edge genes.", f"{acc}/17_gsea/{acc}_{deg}_GSEA_*.csv")
for k, (run, acc, deg, lab) in CTX.items():
    up = rd(run, "13_ilincs_signature", f"{acc}_concordants_UP_raw.csv"); up.insert(0, "Query", "Up-regulated genes")
    dn = rd(run, "13_ilincs_signature", f"{acc}_concordants_DOWN_raw.csv") if (ROOT / run / "13_ilincs_signature" / f"{acc}_concordants_DOWN_raw.csv").exists() else None
    if dn is not None: dn.insert(0, "Query", "Down-regulated genes")
    allc = pd.concat([x for x in (up, dn) if x is not None], ignore_index=True)
    add(f"S6_{k}_screen", allc, f"iLINCS chemical perturbagen screen, every LINCS signature returned for each half of the signature (similarity < 0 = reversal; strong at ≤ −0.321), {lab}.", f"{run}/13_ilincs_signature/{acc}_concordants_*_raw.csv")
    add(f"S6_{k}_compounds", sup[f"S6_{k}"], f"Consensus reversal per compound after the drugfindR paired consensus (strongest reversing signature per compound), {lab}.", f"{run}/13_ilincs_signature/{acc}_iLINCS_candidate_compounds.csv")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S7_{k}_targets", sup[f"S7_{k}_targets"], f"DGIdb drug-target interactions of the reversing compounds (before removal of drug-metabolism and transport genes), with the name used for the query, {lab}.", f"{run}/14_drug_targets/{acc}_drug_targets_ALL.csv")
    add(f"S7_{k}_proximity", sup[f"S7_{k}_proximity"], f"Network proximity (closest measure) of each compound's targets to the disease genes, against 1,000 degree-matched random target sets: Z-score, empirical p and BH-adjusted p, {lab}.", f"{run}/15_network_algorithms/{acc}_NetworkProximity*.csv")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S8_{k}", sup[f"S8_{k}"], f"All pairs of the 25 best-ranked compounds: target overlap, network distances and separation, community overlap, fingerprint Tanimoto similarity, pathway overlap and complementarity label, {lab}.", f"drug_complementarity/{k.replace('-', '_')}/{k.replace('-', '_')}_drug_pairs.csv")
add("S9_macrophage_genes", sup["S9_macrophage_genes"], "GSE165709 alveolar macrophages, M. tuberculosis-infected vs uninfected cells of the same donor: RNA-seq and ATAC-seq response per gene in HIV-negative donors (HC) and people with HIV on ART (PLWH), response difference, and link to the TB gene ranking (Prev_*).", "external_results/GSE165709/integration/GSE165709_gene_integration.csv")
add("S9_chromatin_genesets", sup["S9_chromatin_genesets"], "GSE165709 ATAC-seq, gene-set tests: mean |Wald statistic| of peaks near each gene set against 2,000 random gene sets.", "external_results/GSE165709/atac_geneset/GSE165709_ATAC_geneset_tests.csv")
for k, (run, acc, deg, lab) in CTX.items():
    add(f"S10_{k}_celltype", sup[f"S10_{k}_celltype"], f"Share (%) of module-gene expression by cell type in SCP1749 macaque lung granulomas (pseudo-bulk mean expression; an even spread is about 9%), by module and direction, {lab}.", f"{run}/18_singlecell/{acc}_module_celltype_profile.csv")
add("S10_metabolomics_summary", sup["S10_metabolomics_summary"], "GC6-74 metabolomics (ST000974 plasma, ST000975 serum, ST000976 diluted plasma): analyses, numbers of cases and control samples, metabolites tested and metabolites significant at FDR < 0.05.", "GSE114192/21_metabolomics/progression/GC6-74_progression_summary.csv")
add("S10_metabolomics_key", sup["S10_metabolomics_key"], "GC6-74 metabolomics: amino acids, kynurenine and the kynurenine/tryptophan ratio in future TB cases vs matched non-progressors, by study and analysis (limma; BH adjustment).", "GSE114192/21_metabolomics/progression/GC6-74_amino_acids_and_kynurenine.csv")

# ---- S11 software, versions and parameters
pv = rd("reproducibility/package_versions.csv")
use = {"DESeq2": "Differential expression (RNA-seq)", "apeglm": "Fold-change shrinkage", "limma": "Differential expression (arrays, miRNA, metabolomics)",
       "WGCNA": "Co-expression modules", "clusterProfiler": "Over-representation analysis and GSEA", "ReactomePA": "Reactome over-representation",
       "msigdbr": "MSigDB gene sets", "org.Hs.eg.db": "Gene identifier mapping", "AnnotationDbi": "Gene identifier mapping", "igraph": "Network centrality, communities and drawing",
       "drugfindR": "iLINCS signature preparation and queries", "enrichR": "Enrichr annotation of drug targets", "glmnet": "Elastic-net classifier",
       "ChIPseeker": "ATAC-seq peak annotation", "TxDb.Hsapiens.UCSC.hg38.knownGene": "Transcript annotation for peaks", "GEOquery": "GEO downloads",
       "ggplot2": "Figures", "ggrepel": "Figure labels", "pheatmap": "Heatmaps", "UpSetR": "UpSet plots", "GenomicRanges": "Genomic intervals", "IRanges": "Genomic intervals"}
rows = [["Environment", "R", "4.5.0", "All R analyses", ""], ["Environment", "Bioconductor", "3.22", "Bioconductor packages", ""]]
for r in pv.itertuples():
    if str(r.Version).strip().lower() == "not installed": continue
    rows.append(["R package", r.Package, r.Version, use.get(r.Package, "Support"), ""])
for p, u in [("Python 3", "Drug complementarity, figure assembly and tables"), ("RDKit", "Structure standardisation and Morgan fingerprints"), ("NetworkX", "Network distances for drug pairs"),
             ("pandas, numpy, matplotlib, seaborn, scikit-learn, openpyxl", "Tables, figures and clustering")]:
    rows.append(["Python", p, "not recorded in the repository", u, ""])
ext = rd("reproducibility/external_resource_versions.csv")
for r in ext.itertuples():
    rows.append(["Online resource or database", r.resource, r.version_evidence, f"Run: {r.run}", re.sub(r"\s*\(scripts/[^)]*\)", "", str(r.settings)) if isinstance(r.settings, str) and not r.settings.startswith("scripts/") else ""])
P = [("Differential expression", "Threshold", "adjusted p < 0.05 and |log2FC| ≥ 1 (Benjamini-Hochberg)"),
     ("RNA-seq filtering", "Threshold", "≥ 10 reads in at least as many samples as the smallest group"),
     ("Arrays", "Normalisation and filtering", "log2, quantile normalisation; probe detected (p < 0.05) in at least as many samples as the smallest group; most expressed probe per gene"),
     ("WGCNA", "Network", "unsigned network and TOM; minimum module size 30; merge eigengenes correlated > 0.75"),
     ("WGCNA", "Soft-thresholding power", "lowest power with scale-free fit ≥ 0.85, else 6: " + "; ".join(
         f"{lab} {int(pd.read_csv(ROOT / run / '10_wgcna_results' / f'{acc}_WGCNA_soft_threshold.csv').power_used.iloc[0])}"
         for lab, (run, acc) in [("TB", ("GSE114192", "GSE114192")), ("TB with diabetes", ("GSE114192/TB_DM", "GSE114192")), ("TB with HIV", ("GSE37250/TB_HIV", "GSE37250"))])),
     ("WGCNA", "Module selection", "three modules with the strongest |correlation| at FDR < 0.05"),
     ("STRING", "Confidence", "combined score ≥ 400 (0.4), Homo sapiens"),
     ("Random walk with restart", "Parameters", "restart probability 0.5; tolerance 1e-10; at most 200 iterations; seeds = disease genes"),
     ("Network diffusion", "Parameters", "α = 0.5; tolerance 1e-10; at most 200 iterations; seeds = disease genes"),
     ("Evidence families", "Scoring", "percentile per measure (−log10 p capped at 10); family = mean; support ≥ 0.80; ≥ 3 families required"),
     ("Evidence families", "Tiers", "High: ≥ 4 of 5 families and score ≥ 0.75; Medium: ≥ 2 families or score ≥ 0.60; otherwise Low"),
     ("GSEA", "Parameters", "ranking sign(log2FC) × −log10 p; gene sets of 15 to 500 genes; eps = 0; seed 42; BH"),
     ("Classification", "Parameters", "elastic net α = 0.5, lambda.min by AUC in 5-fold CV; nested CV 5 folds × 10 repeats, stratified outer folds; 1,000 bootstrap resamples; 1,000 random gene sets"),
     ("Drug signature reversal", "Parameters", "genes with padj < 0.05 in L1000 space; |log2FC| ≥ 0.5 per half; halves with < 5 genes not queried; reversal similarity ≤ −0.321"),
     ("Drug targets", "Filter", "CYP, UGT, SULT, SLC and ABC families and NFE2L2 excluded"),
     ("Network proximity", "Parameters", "closest measure within the combined module network; 1,000 degree-matched random target sets; disease genes fixed; BH"),
     ("Drug pairs", "Parameters", "Morgan fingerprints radius 2, 2,048 bits; 500 degree-matched random sets for nearest-target distance; separation only with ≥ 2 network targets per compound; silhouette ≥ 0.25"),
     ("Shortlist", "Rule", "strong reversal of both halves and up-half reversal in ≥ 60% of LINCS signatures; no cell-division target = recommended"),
     ("Random seed", "All steps", "42")]
for a, b, c in P: rows.append(["Parameter", a, "", b, c])
s11 = pd.DataFrame(rows, columns=["Category", "Software, resource or step", "Version", "Used for / setting", "Parameters or notes"])
add("S11_software", s11, "Software, versions, online resources and key parameters. Python package versions were not recorded in the repository. Online resources were queried at run time; STRING version for the networks reported here was not recorded (current scripts pin 12.5).", "reproducibility/package_versions.csv; reproducibility/external_resource_versions.csv; <run>/10_wgcna_results/*_WGCNA_soft_threshold.csv; step scripts")

# ---- extended main tables
EXT = [("Table3_top30_genes", "Table3_priority_genes", "Extended Table 3: the 30 best-ranked genes of each context with all evidence."),
       ("Table4_all_tests", "Table4_classifier", "Extended Table 4: every classification test, with cases, controls, AUC, 95% CI and comparison with random gene sets."),
       ("Table5_all_candidates", "Table5_candidate_drugs", "Extended Table 5: the 25 best-ranked compounds of each context with PubChem CID and targets in the disease network."),
       ("Table6_all_candidates", "Table6_drug_shortlist", "Extended Table 6: evidence group, rank range, proximity and targets of all 25 compounds per context.")]
EXTS = {}
for new, old, desc in EXT:
    EXTS[new] = main[old]; DESC.append((new, desc, len(main[old]), "scripts/manuscript_tables.py"))

# ---- write
order = ["Contents"] + [n for n, *_ in EXT] + list(S)
contents = pd.DataFrame([(n, d, r, s) for n, d, r, s in [(e[0], e[1], e[2], e[3]) for e in DESC]], columns=["Sheet", "Description", "Rows", "Source in the repository"])
contents = contents.set_index("Sheet").loc[[n for n in order if n != "Contents"]].reset_index()
with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    head = pd.DataFrame({"A": ["Additional file 1. Supplementary Tables S1–S11 and extended main tables",
                               "Contexts: TB = TB vs healthy controls (GSE114192); TB-DM = TB with type 2 diabetes vs diabetes alone (GSE114192); TB-HIV = active vs latent TB in people with HIV (GSE37250).",
                               "All values are pipeline outputs; nothing was typed by hand. Empty cells = not available.", ""]})
    head.to_excel(w, sheet_name="Contents", index=False, header=False)
    contents.to_excel(w, sheet_name="Contents", index=False, startrow=4)
    for n in order[1:]:
        df = EXTS.get(n, S.get(n)); df.to_excel(w, sheet_name=n[:31], index=False)
wb = load_workbook(OUT)
for ws in wb.worksheets:
    if ws.title == "Contents":
        ws["A1"].font = Font(bold=True, size=12)
        for c in ws[5]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DDDDDD")
        for col, wd in zip("ABCD", (24, 110, 9, 70)): ws.column_dimensions[col].width = wd
        for row in ws.iter_rows(min_row=6):
            for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
        continue
    ws.freeze_panes = "A2"
    for c in ws[1]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DDDDDD")
    for i, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
        L = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[col[0].column_letter].width = min(50, max(9, L + 2))
wb.save(OUT)
print("sheets:", len(wb.sheetnames)); print(contents[["Sheet", "Rows"]].to_string())
