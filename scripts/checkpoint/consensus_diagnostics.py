#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: diagnostics for the existing consensus gene and drug rankings (step 15).

DIAGNOSTIC ONLY. Official tiers in 19_consensus/*.csv are NOT changed. All alternative tiers written here are
DRAFTS for open decisions D4 (gene scheme) and D5 (proximity in drug tiers).
Outputs: <run>/19_consensus/diagnostics_2026-10-14/
"""
import os, itertools
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUNS = {"reference": "GSE114192", "TB_DM_preliminary": "GSE114192/TB_DM"}
TOP = 0.20

def top_flag(v, frac=TOP, higher=True):
    v = pd.Series(v, dtype=float)
    if not higher: v = -v
    thr = v.quantile(1 - frac)          # numpy/pandas 'linear' == R quantile type 7
    return (v >= thr) & v.notna()

def pct(v):                              # within-context percentile rank (0-1], NA kept
    return pd.Series(v, dtype=float).rank(pct=True, method="average")

summary = []
for label, rel in RUNS.items():
    base = os.path.join(ROOT, rel, "19_consensus"); out = os.path.join(base, "diagnostics_2026-10-14")
    os.makedirs(out, exist_ok=True)
    g = pd.read_csv(os.path.join(base, "GSE114192_GenePriority_all_modules.csv"))
    # ---- 1. reproduce the official gene criteria exactly ----
    crit = {"MM": top_flag(g.MM.abs()), "GS": top_flag(g.Gene_Significance.abs()),
            "Degree": top_flag(g.PPI_Degree), "RWR": top_flag(g.RWR_score), "Diffusion": top_flag(g.Diffusion_score)}
    for k, c in [("MM", "Top_MM"), ("GS", "Top_GS"), ("Degree", "Top_PPI_Degree"), ("RWR", "Top_RWR"), ("Diffusion", "Top_Diffusion")]:
        assert (crit[k].values == g[c].astype(bool).values).all(), f"{label}: {k} flag not reproduced"
    votes = pd.DataFrame({"MM": crit["MM"], "GS": crit["GS"], "Degree": crit["Degree"],
                          "Propagation": crit["RWR"] | crit["Diffusion"]})
    n = votes.sum(1).astype(int)
    assert (n.values == g.N_criteria.values).all(), f"{label}: N_criteria not reproduced"
    avail = (g[["MM", "Gene_Significance"]].notna().sum(1) + g.PPI_Degree.notna().astype(int)
             + (g.RWR_score.notna() | g.Diffusion_score.notna()).astype(int))
    g2 = g[["Gene", "Ensembl", "Module", "Final_Priority", "DEG_WGCNA_overlap", "In_PPI_network"]].copy()
    g2["N_criteria"] = n; g2["max_achievable_criteria"] = avail
    g2["DRAFT_tier_D4"] = np.select([n >= 4, n >= 2, n == 1], ["High", "Medium", "Low"], "Unsupported (0 criteria)")
    g2["DRAFT_assessability"] = np.where(avail < 4, "partially assessable (no STRING node)", "fully assessable")
    # ---- 2. distribution ----
    dist = pd.crosstab([g2.N_criteria], [g2.In_PPI_network], margins=True)
    dist.to_csv(os.path.join(out, "gene_criteria_distribution.csv"))
    # ---- 3. sensitivity analyses (sets keyed by Ensembl ID; some genes lack symbols) ----
    base_high = set(g.Ensembl[g.Final_Priority == "High"])
    rows = []
    def record(name, high_set, n_med=None):
        j = len(base_high & high_set) / max(1, len(base_high | high_set))
        rows.append(dict(scenario=name, n_High=len(high_set), overlap_with_official_High=len(base_high & high_set),
                         jaccard_vs_official=round(j, 3)))
    record("official (top 20%, RWR|Diffusion merged, High>=4)", base_high)
    for frac in (0.10, 0.30):
        v = pd.DataFrame({"MM": top_flag(g.MM.abs(), frac), "GS": top_flag(g.Gene_Significance.abs(), frac),
                          "Degree": top_flag(g.PPI_Degree, frac),
                          "Prop": top_flag(g.RWR_score, frac) | top_flag(g.Diffusion_score, frac)}).sum(1)
        record(f"top {int(frac*100)}% instead of 20%", set(g.Ensembl[v >= 4]))
    for drop in votes.columns:
        v = votes.drop(columns=drop).sum(1)
        record(f"leave out {drop} (High = all 3 remaining)", set(g.Ensembl[v >= 3]))
    v5 = pd.DataFrame(crit).sum(1)
    record("RWR and Diffusion counted separately (5 criteria, High>=4)", set(g.Ensembl[v5 >= 4]))
    # family-weighted percentile draft (Technical doc sec.21): WGCNA, Topology, Propagation; equal weights;
    # missing family = unknown -> mean over available families, coverage reported
    fam = pd.DataFrame({"WGCNA": (pct(g.MM.abs()) + pct(g.Gene_Significance.abs())) / 2,
                        "Topology": pct(g.PPI_Degree),
                        "Propagation": pd.concat([pct(g.RWR_score), pct(g.Diffusion_score)], axis=1).mean(1)})
    g2["DRAFT_family_score"] = fam.mean(1, skipna=True).round(4); g2["family_coverage"] = fam.notna().sum(1)
    k = len(base_high)
    fam_top = set(g2.sort_values("DRAFT_family_score", ascending=False).Ensembl.head(k))
    record(f"family-weighted percentile, top {k} genes (same size as official High)", fam_top)
    rho = pd.Series(g2.DRAFT_family_score).rank().corr(pd.Series(n).rank())
    sens = pd.DataFrame(rows); sens.to_csv(os.path.join(out, "gene_consensus_sensitivity.csv"), index=False)
    g2.to_csv(os.path.join(out, "gene_consensus_DRAFT_tiers.csv"), index=False)
    # metric redundancy
    cor = g[["MM", "Gene_Significance", "PPI_Degree", "PPI_Betweenness", "RWR_score", "Diffusion_score"]].abs().corr(method="spearman")
    cor.round(3).to_csv(os.path.join(out, "gene_metric_spearman_correlations.csv"))

    # ---- 4. drug provenance and alternative proximity rules ----
    d = pd.read_csv(os.path.join(base, "GSE114192_DrugPriority.csv"))
    cc = ["Strong_reversal", "Negative_prox_Z", "Significant_prox_P", "Targets_in_TBnet", "Module_enriched", "Hits_high_gene", "Top_target_RWR"]
    B = d[cc].astype(str).str.upper().eq("TRUE") if False else d[cc].apply(lambda c: c.astype(str).str.upper().eq("TRUE"))
    assert (B.sum(1).values == d.N_criteria.values).all(), f"{label}: drug N_criteria not reproduced"
    prox_votes = B.Negative_prox_Z.astype(int) + B.Significant_prox_P.astype(int)
    tier = lambda s: np.select([s >= 4, s >= 2], ["High", "Medium"], "Low")
    alt_A = d.Final_Priority
    alt_B = tier(B.drop(columns=["Negative_prox_Z", "Significant_prox_P"]).sum(1) + (B.Negative_prox_Z | B.Significant_prox_P).astype(int))
    fdr = pd.read_csv(os.path.join(ROOT, rel, "15_network_algorithms", "GSE114192_NetworkProximity_withFDR.csv"))
    fdr_sig = d.Drug.map(dict(zip(fdr.Drug, fdr["FDR_lt_0.05"].astype(str).str.upper().eq("TRUE")))).fillna(False)
    alt_B2 = tier(B.drop(columns=["Negative_prox_Z", "Significant_prox_P"]).sum(1) + fdr_sig.astype(int))
    alt_C = tier(B.drop(columns=["Negative_prox_Z", "Significant_prox_P"]).sum(1))
    prov = d[["Drug", "Final_Priority", "N_criteria"]].copy()
    prov["criteria_met"] = B.apply(lambda r: "+".join(c for c in cc if r[c]), axis=1)
    prov["n_proximity_votes"] = prox_votes
    prov["proximity_FDR_significant"] = fdr_sig
    prov["N_criteria_without_proximity"] = B.drop(columns=["Negative_prox_Z", "Significant_prox_P"]).sum(1)
    prov["DRAFT_A_current"] = alt_A
    prov["DRAFT_B_one_proximity_vote_nominal"] = alt_B
    prov["DRAFT_B2_one_proximity_vote_FDR"] = alt_B2
    prov["DRAFT_C_proximity_tiebreak_only"] = alt_C
    prov["tie_break_proximity_Z"] = d.Network_Proximity_Z
    prov["note"] = "Computational candidate (hypothesis), not a validated treatment"
    prov.to_csv(os.path.join(out, "drug_criteria_provenance_and_DRAFT_rules.csv"), index=False)
    hi = d.Final_Priority == "High"
    summary.append(dict(run=label, genes=len(g), official_High=len(base_high),
        genes_0_criteria=int((n == 0).sum()), genes_1_criterion=int((n == 1).sum()),
        genes_not_fully_assessable=int((avail < 4).sum()), spearman_family_vs_count=round(rho, 3),
        drugs=len(d), drugs_High=int(hi.sum()), High_with_any_proximity_vote=int((hi & (prox_votes > 0)).sum()),
        High_dropping_below_High_without_proximity=int((hi & (prov.N_criteria_without_proximity < 4)).sum()),
        drugs_nominal_prox=int(B.Significant_prox_P.sum()), drugs_FDR_prox=int(fdr_sig.sum()),
        High_A=int((alt_A == "High").sum()), High_B=int((alt_B == "High").sum()),
        High_B2=int((alt_B2 == "High").sum()), High_C=int((alt_C == "High").sum())))
    print(label, "done:", sens.to_string(index=False))
pd.DataFrame(summary).to_csv(os.path.join(ROOT, "reproducibility", "checkpoint_2026-10-14", "consensus_diagnostics_summary.csv"), index=False)
print(pd.DataFrame(summary).T.to_string())
