#!/usr/bin/env python3
"""Compare the evidence-family ranking with a ranking by the expression family alone
(a differential-expression-only ranking) and count genes without a network score.
Reads <run>/20_evidence_families/<ACC>_Gene_evidence_families.csv of each context (written by evidence_family_consensus_step25.R).
Writes <run>/20_evidence_families/<ACC>_expression_only_comparison.csv:
  N_genes, N_without_network_score, Top50_overlap_with_expression_only, Top50_strong_expression
  (top-50 genes with expression family >= 0.80), Spearman_rank_vs_expression_family.
Ties in the expression family at the top-50 boundary are broken by the order in the file."""
import sys
from pathlib import Path
import pandas as pd

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
RUNS = {"TB": ("GSE114192", "GSE114192"), "TB-DM": ("GSE114192/TB_DM", "GSE114192"), "TB-HIV": ("GSE37250/TB_HIV", "GSE37250")}
TOP, STRONG = 50, 0.80
for ctx, (run, acc) in RUNS.items():
    f = ROOT / run / "20_evidence_families" / f"{acc}_Gene_evidence_families.csv"
    if not f.exists(): print(f"  (missing: {f.relative_to(ROOT)})"); continue
    d = pd.read_csv(f)
    top = d.sort_values("Evidence_rank", kind="stable").head(TOP)
    expr = d.sort_values("Fam_Expression", ascending=False, kind="stable").head(TOP)
    rho = d["Evidence_rank"].corr(d["Fam_Expression"].rank(ascending=False), method="spearman")
    out = pd.DataFrame([{"Context": ctx, "N_genes": len(d), "N_without_network_score": int(d.Fam_Network.isna().sum()),
                         "Top50_overlap_with_expression_only": len(set(top.Ensembl) & set(expr.Ensembl)),
                         "Top50_strong_expression": int((top.Fam_Expression >= STRONG).sum()),
                         "Spearman_rank_vs_expression_family": round(rho, 3)}])
    o = f.with_name(f"{acc}_expression_only_comparison.csv"); out.to_csv(o, index=False)
    print(out.to_string(index=False))
