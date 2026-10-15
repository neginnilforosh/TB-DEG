#!/usr/bin/env python3
"""Build the GSE114192 TB-DM design record.

NON-INFERENTIAL. Reads existing metadata/QC files only; fits no model; writes new files only.
Outputs (GSE114192/TB_DM/):
  00_sample_manifest.csv            one row per GEO sample (249)
  00_design/site_by_group.csv       site x GEO label counts
  00_design/design_identifiability.csv  rank / estimability checks for candidate designs
Inputs:
  00_design/GSE114192_geo_series_matrix_metadata.tsv  (GEO series matrix)
  00_design/GSE114192_htseq_library_qc.tsv            (HTSeq special counters per sample)
  00_design/GSE114192_supp_table1_alignment_match.csv (Eckold et al. Suppl. Table 1 matched by aligned reads)
  02_metadata/GSE114192_TB_DM_sample_metadata.csv, ../02_metadata/GSE114192_sample_metadata.csv
  figures_manuscript/GSE114192_A_sample_outliers.csv
"""
import csv, os, itertools
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBDM = os.path.join(ROOT, "GSE114192", "TB_DM")
D = os.path.join(TBDM, "00_design")

geo = pd.read_csv(os.path.join(D, "GSE114192_geo_series_matrix_metadata.tsv"), sep="\t")
qc = pd.read_csv(os.path.join(D, "GSE114192_htseq_library_qc.tsv"), sep="\t")
supp = pd.read_csv(os.path.join(D, "GSE114192_supp_table1_alignment_match.csv"))
tbdm = pd.read_csv(os.path.join(TBDM, "02_metadata", "GSE114192_TB_DM_sample_metadata.csv"))
ref = pd.read_csv(os.path.join(ROOT, "GSE114192", "02_metadata", "GSE114192_sample_metadata.csv"))
outl = pd.read_csv(os.path.join(TBDM, "figures_manuscript", "GSE114192_A_sample_outliers.csv"))

assert len(geo) == 249 and geo.gsm.is_unique
PIPE = {"Healthy_Control": "Healthy_Control", "TB_only": "TB_Only", "TB_DM": "TB_DM",
        "DM_only": "DM_only", "TB_IH": "Exclude", "IH": "Exclude"}
m = geo.merge(qc, on="gsm", how="left", validate="1:1")
m["analysis_group"] = m.disease_state_GEO.map(PIPE)
m["unique_reads_htseq"] = m.assigned + m.no_feature + m.ambiguous
m["multimapped_reads_htseq"] = m.alignment_not_unique
m["pct_unique_assigned"] = (100 * m.assigned / m.unique_reads_htseq).round(1)

# cross-checks against the existing pipeline metadata
t = dict(zip(tbdm.sample_id, tbdm.group)); r = dict(zip(ref.sample_id, ref.group))
m["in_TBDM_run_metadata"] = m.gsm.isin(t)
m["TBDM_run_group"] = m.gsm.map(t)
m["in_reference_run_metadata"] = m.gsm.isin(r)
assert (m.loc[m.in_TBDM_run_metadata, "TBDM_run_group"] == m.loc[m.in_TBDM_run_metadata, "analysis_group"]).all()
assert set(m.loc[m.analysis_group != "Exclude", "gsm"]) == set(t)

# Suppl. Table 1 match (by STAR-aligned unique reads ~ HTSeq assigned+no_feature+ambiguous)
best = supp.sort_values("absdiff").drop_duplicates("gsm")
best = best[best.absdiff <= 1000]
m = m.merge(best[["gsm", "supp_row", "supp_aligned", "absdiff"]].rename(columns={
    "supp_row": "suppT1_row", "supp_aligned": "suppT1_aligned_reads", "absdiff": "suppT1_absdiff"}), on="gsm", how="left")
m["matched_in_suppT1"] = m.suppT1_row.notna()

o = dict(zip(outl.Sample, outl.Outlier)); oc = dict(zip(outl.Sample, outl.Mean_correlation))
m["outlier_flag_meancorr"] = m.gsm.map(o)
m["mean_correlation"] = m.gsm.map(oc)

core = m.analysis_group != "Exclude"
m["included_core_195"] = core
m["included_primary_TBDM_vs_DM"] = m.analysis_group.isin(["TB_DM", "DM_only"])
m["included_secondary_TBDM_vs_TB"] = m.analysis_group.isin(["TB_DM", "TB_Only"])
m["included_interaction_4groups"] = core
m["site_has_HC_and_DMonly"] = m.field_site_GEO.isin(["South_Africa", "Romania"])

def reason(row):
    if row.disease_state_GEO == "IH":
        return "Excluded: IH (intermediate hyperglycaemia, no TB); not part of the specified TB-DM comparisons; not reported in Eckold et al."
    if row.disease_state_GEO == "TB_IH":
        return "Excluded: TB with intermediate hyperglycaemia (HbA1c 5.7-<6.5%); not part of the specified comparisons; candidate secondary analysis"
    return "Included: core group for the specified TB-DM comparisons"
m["inclusion_reason"] = m.apply(reason, axis=1)

flags = []
for _, x in m.iterrows():
    f = []
    if x.outlier_flag_meancorr is True: f.append("outlier_meancorr")
    if x.field_site_GEO in ("Indonesia", "Peru") and x.analysis_group in ("TB_DM", "TB_Only"):
        f.append("site_without_HC_or_DMonly")
    if not x.matched_in_suppT1: f.append("no_match_in_EckoldSuppT1")
    if str(x.rseq_id_local_file).endswith("_p1"): f.append("rseq_p1_resequenced_id")
    flags.append(";".join(f))
m["qc_and_design_flags"] = flags
m["clinical_covariates_available"] = "none in GEO (age, sex, BMI, HbA1c, glucose, treatment, HIV not deposited)"

cols = ["gsm", "patient_id", "rseq_id_local_file", "disease_state_GEO", "field_site_GEO", "analysis_group",
        "included_core_195", "included_primary_TBDM_vs_DM", "included_secondary_TBDM_vs_TB",
        "included_interaction_4groups", "in_TBDM_run_metadata", "in_reference_run_metadata",
        "inclusion_reason", "site_has_HC_and_DMonly", "assigned", "unique_reads_htseq",
        "multimapped_reads_htseq", "pct_unique_assigned", "matched_in_suppT1", "suppT1_row",
        "suppT1_aligned_reads", "suppT1_absdiff", "outlier_flag_meancorr", "mean_correlation",
        "qc_and_design_flags", "clinical_covariates_available"]
m[cols].to_csv(os.path.join(TBDM, "00_sample_manifest.csv"), index=False)

pd.crosstab(m.disease_state_GEO, m.field_site_GEO, margins=True).to_csv(os.path.join(D, "site_by_group.csv"))

# ---- design identifiability (rank + estimability of the contrast of interest) ----
def X_of(df, terms):
    cols = [np.ones(len(df))]; names = ["(Intercept)"]
    for t_ in terms:
        if t_ == "site":
            for lev in sorted(df.field_site_GEO.unique())[1:]:
                cols.append((df.field_site_GEO == lev).astype(float).values); names.append(f"site{lev}")
        elif t_ == "group":
            levs = sorted(df.analysis_group.unique())
            for lev in levs[1:]:
                cols.append((df.analysis_group == lev).astype(float).values); names.append(f"group{lev}")
        elif t_ in ("TB", "DM", "TB:DM"):
            tb = df.analysis_group.isin(["TB_Only", "TB_DM"]).astype(float).values
            dm = df.analysis_group.isin(["DM_only", "TB_DM"]).astype(float).values
            cols.append({"TB": tb, "DM": dm, "TB:DM": tb * dm}[t_]); names.append(t_)
        elif t_ == "site:group":
            levs = sorted(df.analysis_group.unique())
            for s in sorted(df.field_site_GEO.unique())[1:]:
                for g in levs[1:]:
                    cols.append(((df.field_site_GEO == s) & (df.analysis_group == g)).astype(float).values)
                    names.append(f"site{s}:group{g}")
    return np.column_stack(cols), names

def estimable(X, names, target):
    c = np.zeros(len(names)); c[names.index(target)] = 1
    return np.linalg.matrix_rank(np.vstack([X, c])) == np.linalg.matrix_rank(X)

cases = []
two = ["South_Africa", "Romania"]
specs = [
  ("Primary TB-DM vs DM", ["TB_DM", "DM_only"], None, ["group"], "groupTB_DM", "current (no site)"),
  ("Primary TB-DM vs DM", ["TB_DM", "DM_only"], None, ["site", "group"], "groupTB_DM", "option B: all sites + site"),
  ("Primary TB-DM vs DM", ["TB_DM", "DM_only"], two, ["site", "group"], "groupTB_DM", "option A: SA+RO + site"),
  ("Primary TB-DM vs DM", ["TB_DM", "DM_only"], None, ["site", "group", "site:group"], "groupTB_DM", "site-specific effects, all sites"),
  ("Secondary TB-DM vs TB", ["TB_DM", "TB_Only"], None, ["group"], "groupTB_Only", "current (no site)"),
  ("Secondary TB-DM vs TB", ["TB_DM", "TB_Only"], None, ["site", "group"], "groupTB_Only", "all sites + site"),
  ("Interaction TB x DM", ["Healthy_Control", "TB_Only", "DM_only", "TB_DM"], None, ["TB", "DM", "TB:DM"], "TB:DM", "current (no site)"),
  ("Interaction TB x DM", ["Healthy_Control", "TB_Only", "DM_only", "TB_DM"], None, ["site", "TB", "DM", "TB:DM"], "TB:DM", "all sites + site"),
  ("Interaction TB x DM", ["Healthy_Control", "TB_Only", "DM_only", "TB_DM"], two, ["site", "TB", "DM", "TB:DM"], "TB:DM", "option A: SA+RO + site"),
]
mc = m[core]
for name, groups, sites, terms, target, label in specs:
    df = mc[mc.analysis_group.isin(groups)]
    if sites: df = df[df.field_site_GEO.isin(sites)]
    X, nm = X_of(df, terms)
    rk = np.linalg.matrix_rank(X)
    est = estimable(X, nm, target) if target in nm else False
    n_by = df.groupby(["analysis_group", "field_site_GEO"]).size()
    # samples in sites where all compared groups are present = sites that inform the group effect once site is in the model
    sites_all = [s for s in df.field_site_GEO.unique() if df[df.field_site_GEO == s].analysis_group.nunique() == len(groups)]
    n_inf = int(df.field_site_GEO.isin(sites_all).sum())
    cases.append(dict(comparison=name, design="~ " + " + ".join(terms), variant=label, n_samples=len(df),
        n_columns=X.shape[1], rank=rk, full_rank=rk == X.shape[1], coefficient=target, coefficient_estimable=est,
        sites_with_all_groups=";".join(sorted(sites_all)), n_in_sites_with_all_groups=n_inf,
        counts_group_site="; ".join(f"{g}@{s}={v}" for (g, s), v in n_by.items())))
pd.DataFrame(cases).to_csv(os.path.join(D, "design_identifiability.csv"), index=False)

# Cramer's V between group and site in the 195 core samples
ct = pd.crosstab(mc.analysis_group, mc.field_site_GEO).values.astype(float)
n = ct.sum(); exp = np.outer(ct.sum(1), ct.sum(0)) / n
chi2 = ((ct - exp) ** 2 / exp).sum(); V = np.sqrt(chi2 / (n * (min(ct.shape) - 1)))
with open(os.path.join(D, "site_group_association.txt"), "w") as fh:
    fh.write(f"Core 195 samples: chi-square = {chi2:.1f} (df = {(ct.shape[0]-1)*(ct.shape[1]-1)}), Cramer's V = {V:.3f}\n")
    fh.write("Healthy_Control and DM_only occur only in South_Africa and Romania.\n")
print("manifest rows:", len(m), "| core:", int(core.sum()), "| V =", round(V, 3))
print(pd.DataFrame(cases)[["comparison", "variant", "n_samples", "rank", "full_rank", "coefficient_estimable", "n_in_sites_with_all_groups"]].to_string(index=False))
