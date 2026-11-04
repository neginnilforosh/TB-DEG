#!/usr/bin/env python3
"""Manuscript tables built directly from the pipeline outputs.

Main tables (manuscript_tables/Main_tables.xlsx), numbered as in the manuscript guide:
  Table 1  datasets            Table 2  disease-associated modules     Table 3  prioritised genes (top 30 per context)
  Table 4  classifier          Table 5  candidate drugs                Table 6  final drug shortlist
  Table 7  evidence across omics layers
Supplementary tables (manuscript_tables/Supplementary_tables.xlsx):
  S1 differential expression  S2 network nodes and edges  S3 network scores and communities  S4 evidence-family scores
  S5 GSEA (Hallmark, Reactome) S6 iLINCS screen  S7 drug targets and network proximity  S8 drug pairs
  S9 macrophage infection and chromatin  S10 single-cell and metabolomics  S11 software and versions
Every value is read from a file; nothing is typed by hand except the dataset descriptions in Table 1.
Run from the repository root:  python3 scripts/manuscript_tables.py
"""
import re, sys
from pathlib import Path
import pandas as pd

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
OUT = ROOT / "manuscript_tables"; OUT.mkdir(exist_ok=True)
CTX = {"TB": ("GSE114192", "GSE114192", "HealthyControl_vs_TBOnly"),
       "TB-DM": ("GSE114192/TB_DM", "GSE114192", "DM_only_vs_TB_DM"),
       "TB-HIV": ("GSE37250/TB_HIV", "GSE37250", "LTBI_HIVpos_vs_ATB_HIVpos")}
NOTES = {}

def rd(*p, **kw):
    f = ROOT.joinpath(*p)
    if not f.exists(): print(f"  (missing: {f.relative_to(ROOT)})"); return None
    return pd.read_csv(f, **kw)
def deg_file(run, acc, name):
    for base in (ROOT / run / "06_deg_results" / f"{acc}_{name}_DEG.csv", ROOT / acc / "06_deg_results" / f"{acc}_{name}_DEG.csv"):
        if base.exists(): return pd.read_csv(base)
    print(f"  (missing DEG table for {run})"); return None
nice = lambda d: d if re.search(r"\d", str(d)) else str(d).title()

# ---------------------------------------------------------------- Table 1
T1 = pd.DataFrame([
 ["GSE114192", "Adults with TB, type 2 diabetes, both, or neither; South Africa, Romania, Indonesia, Peru", "Whole blood", "Healthy 36; TB 46; diabetes 52; TB with diabetes 61", "RNA-seq", "Discovery: TB and TB with diabetes"],
 ["GSE37250", "Adults with and without HIV; South Africa and Malawi", "Whole blood", "HIV+: active TB 98, latent TB 84, other diseases 92; HIV-: 97, 83, 83", "Illumina HT-12 v4 array", "Discovery: TB with HIV; TB x HIV interaction"],
 ["GSE161829", "Adults without TB, with latent or active TB", "Memory CD4+ T cells", "No TB 20; latent 40; active 25", "RNA-seq", "Replication of the full pipeline; external validation"],
 ["GSE222001", "Healthy, latent and active TB", "NK cells", "5 per group", "RNA-seq", "External validation (exploratory)"],
 ["GSE99374", "Healthy and latent TB", "CD8+ T cells", "Healthy 19; latent 21", "RNA-seq", "External validation (latent TB)"],
 ["GSE107104", "Adults with advanced HIV, Uganda (ART-naive)", "Whole blood", "TB-HIV 15; HIV only 16", "RNA-seq", "Validation: TB with HIV"],
 ["GSE162164", "Adults with advanced HIV, India", "Whole blood", "TB-HIV 14; HIV only 11", "RNA-seq", "Validation: TB with HIV"],
 ["GSE165709", "Healthy donors and people with HIV", "Alveolar macrophages, infected or not in vitro", "RNA-seq 36 donors; ATAC-seq 31 donors", "RNA-seq, ATAC-seq", "Response to infection; chromatin"],
 ["GSE229020", "Healthy, latent, drug-susceptible and drug-resistant TB", "Monocyte miRNA", "4; 6; 6; 4", "NanoString", "miRNA (exploratory)"],
 ["ST000974-6", "GC6-74 household contacts: future TB cases and matched non-progressors", "Plasma, serum, diluted plasma", "366, 340 and 45 samples", "Untargeted metabolomics", "Metabolites before TB diagnosis"],
 ["SCP1749", "Macaque lung granulomas, 4 weeks after infection", "Single cells", "6 granulomas, 2 animals", "scRNA-seq", "Cell-type origin of genes"]],
 columns=["Dataset", "Population", "Biological material", "Groups and n", "Platform", "Purpose"])

# ---------------------------------------------------------------- per-context tables
T2, T3, T5, T6 = [], [], [], []
S = {}
for k, (run, acc, deg) in CTX.items():
    mt = rd(run, "10_wgcna_results", f"{acc}_ModuleTrait_Correlation.csv")
    ce = rd(run, "12_string_ppi", f"{acc}_COMBINED_centrality.csv")
    ev = rd(run, "20_evidence_families", f"{acc}_Gene_evidence_families.csv")
    gp = rd(run, "19_consensus", f"{acc}_GenePriority_all_modules.csv")
    cards = rd(run, "26_drug_evidence_cards", f"{acc}_drug_evidence_cards.csv")
    dfam = rd(run, "20_evidence_families", f"{acc}_Drug_evidence_families.csv")
    names = {}
    if ev is not None: names = dict(zip(ev.Module, ev.Module_Name))
    if mt is not None:
        sel = mt[mt.Selected == True].copy()
        for r in sel.itertuples():
            hubs = ", ".join(ce[ce.Module == r.Module].sort_values("Degree", ascending=False).Gene.head(5)) if ce is not None else ""
            T2.append([k, r.Module, str(names.get(r.Module, "")).replace("_", " "), r.N_genes, round(r.Correlation, 2), r.FDR, hubs])
    if ev is not None:
        top = ev.sort_values("Evidence_rank").head(30).copy()
        top["Gene"] = top.Gene.fillna(top.Ensembl)
        g = gp.set_index("Ensembl") if gp is not None else None
        for r in top.itertuples():
            x = g.loc[r.Ensembl] if g is not None and r.Ensembl in g.index else None
            if isinstance(x, pd.DataFrame): x = x.iloc[0]
            T3.append([k, int(r.Evidence_rank), r.Gene, r.Module, round(r.log2FC, 2), r.FDR,
                       None if x is None else round(x.MM, 2), None if x is None else round(x.Gene_Significance, 2),
                       None if x is None or pd.isna(x.Rank_PPI_Degree) else int(x.Rank_PPI_Degree), None if x is None or pd.isna(x.Rank_RWR) else int(x.Rank_RWR),
                       round(r.Evidence_score, 3), r.Evidence_tier, "yes" if r.Stable_topN else "no", r.Rank_best_across_variants, r.Rank_worst_across_variants])
    if cards is not None:
        lib = rd("drug_complementarity", "drug_library.csv")
        cid = dict(zip(lib.Drug, lib.CID)) if lib is not None else {}
        for r in cards.sort_values("Evidence_rank").itertuples():
            T5.append([k, int(r.Evidence_rank), nice(r.Drug), cid.get(r.Drug), r.Best_sim_UP, r.Best_sim_DOWN, r.Evidence_tier, r.Targets if isinstance(r.Targets, str) else ""])
            consistent = bool(r.Reverses_both) and (r.Share_reversing_UP if pd.notna(r.Share_reversing_UP) else 0) >= 0.6
            cd = isinstance(r.Cell_division_targets, str) and r.Cell_division_targets.strip() != ""
            grp = "Recommended" if consistent and not cd else ("Caution: possible antiproliferative" if consistent else "Lower support")
            T6.append([k, grp, nice(r.Drug), int(r.Evidence_rank), r.Rank_range, r.Proximity_Z, r.Proximity_FDR,
                       r.High_priority_targets if isinstance(r.High_priority_targets, str) else "", r.Cell_division_targets if isinstance(r.Cell_division_targets, str) else "", r.Reading])
    # supplementary per context
    d = deg_file(run, acc, deg)
    if d is not None and ev is not None: d = d.merge(ev[["Ensembl", "Gene"]].drop_duplicates("Ensembl"), left_on="gene", right_on="Ensembl", how="left").drop(columns="Ensembl")
    S[f"S1_{k}"] = d
    net = rd(run, "12_string_ppi", f"{acc}_COMBINED_STRING_network.csv")
    S[f"S2_{k}_nodes"] = ce; S[f"S2_{k}_edges"] = net[["preferredName_A", "preferredName_B", "score"]] if net is not None else None
    sc = rd(run, "15_network_algorithms", f"{acc}_GeneScores_RWR_Diffusion.csv"); co = rd(run, "15_network_algorithms", f"{acc}_Communities.csv")
    S[f"S3_{k}"] = sc.merge(co[["Gene", "Community"]], on="Gene", how="left") if sc is not None and co is not None else sc
    S[f"S4_{k}"] = ev
    hm = rd(acc, "17_gsea", f"{acc}_{deg}_GSEA_HALLMARK.csv"); rc = rd(acc, "17_gsea", f"{acc}_{deg}_GSEA_REACTOME.csv")
    gs = pd.concat([x.assign(Collection=c) for x, c in [(hm, "Hallmark"), (rc, "Reactome")] if x is not None]) if (hm is not None or rc is not None) else None
    if gs is not None: gs = gs[["Collection", "ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment"]]
    S[f"S5_{k}"] = gs
    S[f"S6_{k}"] = rd(run, "13_ilincs_signature", f"{acc}_iLINCS_candidate_compounds.csv")
    tg = rd(run, "14_drug_targets", f"{acc}_drug_targets_ALL.csv"); pxf = ROOT / run / "15_network_algorithms" / f"{acc}_NetworkProximity_withFDR.csv"
    px = pd.read_csv(pxf) if pxf.exists() else rd(run, "15_network_algorithms", f"{acc}_NetworkProximity.csv")
    S[f"S7_{k}_targets"] = tg; S[f"S7_{k}_proximity"] = px
    S[f"S8_{k}"] = rd("drug_complementarity", k.replace("-", "_"), f"{k.replace('-', '_')}_drug_pairs.csv")
    S[f"S10_{k}_celltype"] = rd(run, "18_singlecell", f"{acc}_module_celltype_profile.csv")

T2 = pd.DataFrame(T2, columns=["Context", "Module", "Annotation", "Genes", "Correlation with disease", "FDR", "Hub genes (top 5 by degree)"])
T3 = pd.DataFrame(T3, columns=["Context", "Rank", "Gene", "Module", "log2FC", "FDR", "Module membership", "Gene significance", "Degree rank", "RWR rank",
                               "Evidence score", "Tier", "Stable in top 50", "Best rank (variants)", "Worst rank (variants)"])
T5 = pd.DataFrame(T5, columns=["Context", "Rank", "Compound", "PubChem CID", "Best iLINCS similarity (up genes)", "Best iLINCS similarity (down genes)", "Evidence tier", "Targets in disease network"])
T6 = pd.DataFrame(T6, columns=["Context", "Group", "Compound", "Evidence rank", "Rank range", "Proximity Z", "Proximity FDR", "High-priority genes targeted", "Cell-division targets", "Evidence summary"])

# ---------------------------------------------------------------- Table 4
T4 = []
for k, (run, acc, _) in CTX.items():
    a = rd(run, "24_ml_validation", f"{acc}_ml_validation_AUC.csv")
    if a is not None:
        for r in a.itertuples(): T4.append([k, r.Test, r.Model, r.N_case, r.N_control, r.AUC, r.CI_low, r.CI_high, r.Random_gene_set_median_AUC, r.P_vs_random])
hv = rd("GSE37250", "TB_HIV", "23_hiv_validation", "TBHIV_validation_signature_AUC.csv")
if hv is not None:
    for r in hv.itertuples(): T4.append(["TB-HIV", f"{r.Dataset} TB-HIV vs HIV only (independent)", "Signature score", None, None, r.AUC, r.CI_low, r.CI_high, r.Random_median, r.P_vs_random])
T4 = pd.DataFrame(T4, columns=["Context", "Test", "Model", "Cases", "Controls", "AUC", "95% CI low", "95% CI high", "Random gene sets (median AUC)", "P vs random"])
NOTES["Table 4"] = "Sensitivity and specificity were not computed; discrimination is reported as AUC. Nested cross-validation rows give the spread over repeats, not a confidence interval."

# ---------------------------------------------------------------- Table 7 and S9/S10
T7 = []
for k, (run, acc, deg) in CTX.items():
    hm = rd(acc, "17_gsea", f"{acc}_{deg}_GSEA_HALLMARK.csv")
    if hm is not None:
        for s in ["HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"]:
            x = hm[hm.ID == s]
            if len(x): T7.append(["Blood transcriptome", k, s.replace("HALLMARK_", "").replace("_", " ").title(), f"NES {x.NES.iloc[0]:+.2f}, FDR {x['p.adjust'].iloc[0]:.1e}"])
mac = rd("external_results", "GSE165709", "integration", "GSE165709_gene_integration.csv")
ev1 = rd("GSE114192", "20_evidence_families", "GSE114192_Gene_evidence_families.csv")
gp1 = rd("GSE114192", "19_consensus", "GSE114192_GenePriority_all_modules.csv")
if mac is not None and gp1 is not None:
    hg = set(gp1.Gene[gp1.Final_Priority == "High"]); m = mac[mac.Gene.isin(hg)]
    T7.append(["Alveolar macrophages (GSE165709)", "TB genes", "Response to M. tuberculosis infection", f"{int(m.RNA_responds_HC.sum())} of {len(hg)} High-tier genes respond (healthy donors); {int(m.RNA_responds_PLWH.sum())} (people with HIV)"])
at = rd("external_results", "GSE165709", "atac_geneset", "GSE165709_ATAC_geneset_tests.csv")
if at is not None: S["S9_chromatin_genesets"] = at
S["S9_macrophage_genes"] = mac
aa = rd("GSE114192", "21_metabolomics", "progression", "GC6-74_amino_acids_and_kynurenine.csv")
if aa is not None:
    S["S10_metabolomics_key"] = aa
    for st, nm in [("ST000974", "plasma"), ("ST000975", "serum")]:
        x = aa[(aa.Study == st) & (aa.Metabolite == "kynurenine/tryptophan ratio") & (aa.Contrast == "near_vs_control")]
        if len(x): T7.append(["Metabolomics (GC6-74)", f"Future TB cases, {nm}", "Kynurenine/tryptophan ratio < 5 months before diagnosis", f"log2FC {x.log2FC.iloc[0]:+.2f}, FDR {x.padj.iloc[0]:.2f}"])
summ = rd("GSE114192", "21_metabolomics", "progression", "GC6-74_progression_summary.csv")
if summ is not None: S["S10_metabolomics_summary"] = summ
T7 = pd.DataFrame(T7, columns=["Layer", "Context or group", "Evidence", "Result"])
NOTES["Table 7"] = "Only completed layers are shown. Metabolites and genes come from different people and disease stages."

# ---------------------------------------------------------------- S11 software
pv = rd("reproducibility", "package_versions.csv")
S["S11_software"] = pv

# ---------------------------------------------------------------- write
def write(path, sheets, notes=None):
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        if notes:
            pd.DataFrame({"Notes": [f"{a}: {b}" for a, b in notes.items()]}).to_excel(w, sheet_name="Notes", index=False)
        for name, df in sheets.items():
            if df is None or len(df) == 0: continue
            df.to_excel(w, sheet_name=name[:31], index=False)
            ws = w.sheets[name[:31]]; ws.freeze_panes = "A2"
            for i, col in enumerate(df.columns, 1):
                q = df[col].dropna().astype(str).str.len().quantile(0.9) if df[col].notna().any() else 10
                width = min(60, max(10, int(q), len(str(col)) + 2))
                ws.column_dimensions[ws.cell(1, i).column_letter].width = width
    print(f">>> wrote {path.relative_to(ROOT)} ({sum(df is not None and len(df) > 0 for df in sheets.values())} sheets)")

write(OUT / "Main_tables.xlsx", {"Table1_datasets": T1, "Table2_modules": T2, "Table3_priority_genes": T3, "Table4_classifier": T4,
                                 "Table5_candidate_drugs": T5, "Table6_drug_shortlist": T6, "Table7_cross_omics": T7}, NOTES)
def order(name):
    ctx = next((c for c in sorted(CTX, key=len, reverse=True) if f"_{c}_" in name + "_"), None)
    return (int(re.match(r"S(\d+)", name).group(1)), list(CTX).index(ctx) if ctx else 9, name)
write(OUT / "Supplementary_tables.xlsx", {n: S[n] for n in sorted(S, key=order)})
print(f"  Table 2: {len(T2)} modules | Table 3: {len(T3)} genes | Table 4: {len(T4)} tests | Table 5: {len(T5)} drugs | Table 7: {len(T7)} rows")
