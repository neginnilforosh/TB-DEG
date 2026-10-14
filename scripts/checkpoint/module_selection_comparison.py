#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: compare WGCNA module-selection rules on the EXISTING module-trait tables (no WGCNA rerun).
Rule in use (WGCNA_TBmodules_step1.R): BH-FDR < 0.05, then top 3 by |r| (MAX_SIG_MODULES = 3).
Alternatives are DRAFTS for decision D3. Grey (unassigned) module is never eligible.
Output: <run>/10_wgcna_results/checkpoint_2026-10-14/module_selection_comparison.csv (+ summary)."""
import os, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
summ = []
for label, rel in {"reference": "GSE114192", "TB_DM_preliminary": "GSE114192/TB_DM"}.items():
    d = os.path.join(ROOT, rel, "10_wgcna_results"); out = os.path.join(d, "checkpoint_2026-10-14"); os.makedirs(out, exist_ok=True)
    mt = pd.read_csv(os.path.join(d, "GSE114192_ModuleTrait_Correlation.csv"))
    mt = mt[mt.Module != "grey"].copy(); mt["abs_r"] = mt.Correlation.abs()
    sig = mt.FDR < 0.05
    rules = {"R0_current_FDR_top3_by_absr": mt.Selected.astype(str).str.upper().eq("TRUE"),
             "R1_all_FDR_lt_0.05": sig,
             "R2_FDR_and_absr_ge_0.3": sig & (mt.abs_r >= 0.3),
             "R3_FDR_and_absr_ge_0.5": sig & (mt.abs_r >= 0.5)}
    rank = mt[sig].abs_r.rank(ascending=False)
    mt["rank_absr_among_FDR_sig"] = rank
    for k, v in rules.items(): mt[k] = v
    mt.sort_values("abs_r", ascending=False).to_csv(os.path.join(out, "module_selection_comparison.csv"), index=False)
    for k, v in rules.items():
        summ.append(dict(run=label, rule=k, n_modules=int(v.sum()), n_genes=int(mt.N_genes[v].sum()),
                         modules=";".join(mt.Module[v].sort_values()), non_grey_modules=len(mt)))
    print(label, "non-grey:", len(mt), "FDR<0.05:", int(sig.sum()))
s = pd.DataFrame(summ); s.to_csv(os.path.join(ROOT, "reproducibility", "checkpoint_2026-10-14", "module_selection_comparison_summary.csv"), index=False)
print(s.to_string(index=False))
