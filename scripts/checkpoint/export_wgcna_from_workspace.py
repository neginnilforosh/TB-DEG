#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: export WGCNA objects already saved in *_WGCNA_workspace.rds (no WGCNA rerun).

Outputs per run (reference = GSE114192/10_wgcna_results, TB-DM = GSE114192/TB_DM/10_wgcna_results):
  GSE114192_ModuleEigengenes_from_workspace.csv   samples x module eigengenes (+ group)
  GSE114192_module_merging_summary.csv            module counts before / after merging
Validation: correlation of each eigengene with the case trait must reproduce ModuleTrait_Correlation.csv
(max abs difference < 1e-12), otherwise the script stops.
Soft-threshold power: NOT stored in the workspace. See ../../reproducibility/checkpoint_2026-10-14/
WGCNA_soft_threshold_reconstructed_from_figure.csv (reconstructed from the saved figure, approximate).
Uses rds_min.py (minimal R serialization reader) because R is not available in the checkpoint environment.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from rds_min import read_rds, as_dataframe

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUNS = [("GSE114192/10_wgcna_results", "TB_Only"), ("GSE114192/TB_DM/10_wgcna_results", "TB_DM")]
for rel, trait in RUNS:
    d_ = os.path.join(ROOT, rel)
    o = read_rds(os.path.join(d_, "GSE114192_WGCNA_workspace.rds"))
    d = dict(zip(o.attr("names").value, o.value))
    rn, me = as_dataframe(d["MEs"]); tr, tc = as_dataframe(d["datTraits"])
    MEs = pd.DataFrame(me, index=rn); T = pd.DataFrame(tc, index=tr)
    assert list(MEs.index) == list(T.index)
    mt = pd.read_csv(os.path.join(d_, "GSE114192_ModuleTrait_Correlation.csv"))
    r = MEs.apply(lambda c: np.corrcoef(c, T[trait])[0, 1])
    dev = (mt.Correlation - mt.Module.map(lambda m: r["ME" + m])).abs().max()
    assert dev < 1e-12, f"eigengene check failed ({dev})"
    out = MEs.copy(); out.insert(0, "group", T.idxmax(axis=1)); out.index.name = "sample_id"
    out.to_csv(os.path.join(d_, "GSE114192_ModuleEigengenes_from_workspace.csv"))
    net = dict(zip(d["net"].attr("names").value, d["net"].value))
    unm = pd.Series(net["unmergedColors"].value).astype(str); merged = pd.Series(d["moduleColors"].value)
    pd.DataFrame({"stage": ["before merging (net$unmergedColors)", "after merging (moduleColors)"],
                  "n_modules_incl_grey": [unm.nunique(), merged.nunique()],
                  "n_genes": [len(unm), len(merged)],
                  "n_grey_genes": [int(unm.isin(["grey", "0", "0.0"]).sum()), int((merged == "grey").sum())]}
                 ).to_csv(os.path.join(d_, "GSE114192_module_merging_summary.csv"), index=False)
    print(rel, "eigengenes", MEs.shape, "check max |dr| =", dev)
