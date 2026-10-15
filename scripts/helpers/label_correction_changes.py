#!/usr/bin/env python3
"""Rebuild the record of results before and after the sample-label correction (GSE161829, GSE222001).

"Before" values are read from the repository history (the commit just before the label fix);
"After" values from the files as they are now. Writes, in reproducibility/:
    label_correction_changes.csv   one row per result
    label_correction_changes.md    short explanation + the same table
Run from anywhere inside the repository:
    python3 scripts/helpers/label_correction_changes.py
    python3 scripts/helpers/label_correction_changes.py --before <commit>   # to compare with another commit
"""
import argparse, io, os, re, subprocess, sys
import pandas as pd

FIX_MESSAGE = "Name-based sample matching"   # start of the commit message that applied the label fix

def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout

root = git("rev-parse", "--show-toplevel").strip()
os.chdir(root)
ap = argparse.ArgumentParser()
ap.add_argument("--before", help="commit holding the results computed with the wrong labels (default: parent of the fix commit)")
args = ap.parse_args()
if args.before:
    BEFORE = args.before
else:
    fix = git("log", "--format=%h", f"--grep={FIX_MESSAGE}").split()
    if not fix:
        sys.exit(f"No commit message starting with '{FIX_MESSAGE}' found; pass --before <commit>.")
    BEFORE = git("rev-parse", "--short", fix[-1] + "^").strip()
print(f">>> 'Before' = commit {BEFORE}; 'After' = files on disk")

def old(path):
    try:
        return io.StringIO(git("show", f"{BEFORE}:{path}"))
    except subprocess.CalledProcessError:
        return None
def now(path):
    return path if os.path.exists(path) else None

rows = []
def add(area, result, before, after, note=""):
    rows.append(dict(Area=area, Result=result, Before=before, After=after, Note=note))
def both(path, fn):
    o, n = old(path), now(path)
    return (fn(pd.read_csv(o)) if o else "not available", fn(pd.read_csv(n)) if n else "not available")

# ---- differential expression ----
nsig = lambda d: int(((d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)).sum())
notes = {"TBneg_vs_ATB": "7 samples had the wrong group before the fix",
         "TBneg_vs_LTBI": "7 samples had the wrong group before the fix",
         "Healthy_vs_Active": "groups swapped before the fix: this was latent vs active TB",
         "Healthy_vs_Latent": "groups swapped before the fix: sign was reversed"}
for ds, comps in [("GSE161829", ["TBneg_vs_ATB", "TBneg_vs_LTBI"]), ("GSE222001", ["Healthy_vs_Active", "Healthy_vs_Latent"])]:
    for c in comps:
        b, a = both(f"{ds}/06_deg_results/{ds}_{c}_DEG.csv", lambda d: nsig(d.dropna(subset=["padj"])))
        add("Differential expression", f"{ds} {c}: significant genes (padj<0.05, |log2FC|>=1)", b, a, notes[c])

# ---- GSEA ----
for ds, c in [("GSE161829", "TBneg_vs_ATB"), ("GSE222001", "Healthy_vs_Active")]:
    p = f"{ds}/17_gsea/{ds}_{c}_GSEA_HALLMARK.csv"
    for s in ["HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"]:
        def f(d, s=s):
            r = d[d.ID == s]
            return f"NES {r.NES.iloc[0]:+.2f}, padj {r['p.adjust'].iloc[0]:.1e}" if len(r) else "not tested"
        b, a = both(p, f); add("GSEA", f"{ds} {c}: {s.replace('HALLMARK_', '')}", b, a)
    b, a = both(p, lambda d: int((d["p.adjust"] < 0.05).sum())); add("GSEA", f"{ds} {c}: Hallmark sets with padj<0.05", b, a)

# ---- replication summary (step 14) ----
sp = "comparisons/GSE114192_TB_vs_GSE161829_ATB/SUMMARY.txt"
so = old(sp).read() if old(sp) else ""; sn = open(sp).read() if now(sp) else ""
grab = lambda t, pat: (re.search(pat, t, re.S).group(1) if re.search(pat, t, re.S) else "n/a")
for lab, pat in [("log2FC correlation between cohorts", r"Spearman correlation of log2FC = ([\-\d.]+)"),
                 ("genes significant in both (|log2FC|>=1)", r"in both (\d+) \(expected"),
                 ("Hallmark sets replicated (FDR<0.05 in both, same sign)", r"GSEA HALLMARK:.*?same sign\) (\d+)"),
                 ("rank agreement of compound scores (Spearman)", r"compounds scored in both: Spearman of similarity ([\-\d.]+)"),
                 ("same LINCS signatures, Spearman of similarity", r"same LINCS signatures scored in both: \d+; Spearman of similarity ([\-\d.]+)")]:
    add("Replication (GSE161829 vs reference)", lab, grab(so, pat), grab(sn, pat))
mod = lambda t: (lambda m: f"{m.group(1)}: {m.group(2)} shared genes, {m.group(3)}% same direction" if m else "n/a")(
    re.search(r"green \(\d+ genes\) -> (\w+).*?(\d+) shared.*?same direction ([\d.]+)%", t, re.S))
add("Replication (GSE161829 vs reference)", "best match of the interferon module (reference green)", mod(so), mod(sn))

# ---- GSE161829 network run ----
r = "GSE161829/TBneg_vs_ATB"
b, a = both(f"{r}/10_wgcna_results/GSE161829_ModuleTrait_Correlation.csv",
            lambda d: int((d[[c for c in d.columns if "FDR" in c][0]] < 0.05).sum()))
add("GSE161829 network run", "modules passing FDR<0.05", b, a, "module colours are not comparable between runs")
b, a = both(f"{r}/19_consensus/GSE161829_GenePriority_all_modules.csv", lambda d: int((d.Final_Priority == "High").sum()))
add("GSE161829 network run", "High-priority genes (criteria count)", b, a)
hb, ha = both(f"{r}/19_consensus/GSE161829_DrugPriority.csv", lambda d: set(d.Drug[d.Final_Priority == "High"]))
if isinstance(hb, set) and isinstance(ha, set):
    add("GSE161829 network run", "High-priority drugs", len(hb), len(ha), f"{len(hb & ha)} in both")

# ---- cross-dataset overlap ----
up = "cross_dataset_comparison/cross_dataset_shared_unique_DEGs.csv"
def shared(d):
    cols = [c for c in d.columns if c.startswith("GSE")]
    n = d["N_comparisons"] if "N_comparisons" in d else d[cols].sum(axis=1)
    return sorted(d.gene[n >= 2])
b, a = both(up, shared)
add("Cross-dataset overlap", "genes significant in 2+ comparisons", len(b) if isinstance(b, list) else b, len(a) if isinstance(a, list) else a)
add("Cross-dataset overlap", "list", ", ".join(b) if isinstance(b, list) else b, ", ".join(a) if isinstance(a, list) else a)

# ---- external validation ----
for comp in ["GSE161829:TBneg_vs_ATB", "GSE222001:Healthy_vs_Active", "GSE222001:Healthy_vs_Latent"]:
    def f(d, comp=comp):
        x = d[(d.Comparison == comp) & (d.Tier == "High")]
        return f"{x.Pct_same_direction.iloc[0]}% same direction, {x.Pct_reproduced.iloc[0]}% reproduced" if len(x) else "n/a"
    b, a = both("GSE114192/23_external_validation/GSE114192_external_validation_summary.csv", f)
    add("External validation (High genes)", comp, b, a)

# ---- machine learning ----
for t in ["GSE161829 ATB vs TBneg", "GSE161829 ATB vs LTBI", "GSE161829 LTBI vs TBneg", "GSE222001 Active vs Healthy (n=5+5)"]:
    def f(d, t=t):
        x = d[(d.Test == t) & (d.Model == "Signature score")]
        return x.AUC.iloc[0] if len(x) else "n/a"
    b, a = both("GSE114192/24_ml_validation/GSE114192_ml_validation_AUC.csv", f)
    add("Machine learning (signature-score AUC)", t, b, a, "nested cross-validation introduced at the same time")

# ---- write ----
out = pd.DataFrame(rows)
os.makedirs("reproducibility", exist_ok=True)
out.to_csv("reproducibility/label_correction_changes.csv", index=False)
md = [f"# Results before and after the sample-label correction", "",
      "Two validation cohorts had sample labels assigned by position instead of by name. In GSE161829, from Sample74 "
      "onward each count column was shifted by one GEO sample (7 samples in the wrong group); in GSE222001 the healthy "
      "and latent groups were swapped. Both run scripts now match samples by name and check against the GEO sample "
      "descriptions, and every result that depended on these cohorts was regenerated.", "",
      f"\"Before\" values come from commit `{BEFORE}`, the last version computed with the wrong labels; \"After\" values "
      "from the current files. The machine-learning validation also switched to nested cross-validation at the same time.", "",
      "| Area | Result | Before | After | Note |", "| --- | --- | --- | --- | --- |"]
md += [f"| {r.Area} | {r.Result} | {r.Before} | {r.After} | {r.Note} |" for r in out.itertuples()]
open("reproducibility/label_correction_changes.md", "w").write("\n".join(md) + "\n")
print(out.to_string(index=False, max_colwidth=55))
print("\n>>> wrote reproducibility/label_correction_changes.csv and .md")
