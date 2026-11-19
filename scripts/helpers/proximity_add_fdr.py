#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: add Benjamini-Hochberg FDR to existing network-proximity tables (no recomputation).
Writes *_NetworkProximity_withFDR.csv next to the original; original columns verified unchanged."""
import csv, os, math
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
def bh(ps):
    idx = [i for i, p in enumerate(ps) if p is not None]; m = len(idx)
    order = sorted(idx, key=lambda i: ps[i]); q = [None] * len(ps); prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]; prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q, m
for rel in ["GSE114192/15_network_algorithms/GSE114192_NetworkProximity.csv",
            "GSE114192/TB_DM/15_network_algorithms/GSE114192_NetworkProximity.csv",
            "GSE37250/TB_HIV/15_network_algorithms/GSE37250_NetworkProximity.csv"]:
    src = os.path.join(ROOT, rel); rows = list(csv.DictReader(open(src))); hdr = list(rows[0])
    ps = [None if r["P_empirical"] in ("NA", "") else float(r["P_empirical"]) for r in rows]
    q, m = bh(ps)
    out = src[:-4] + "_withFDR.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr + ["FDR_BH", "n_tests_in_FDR", "nominal_p_lt_0.05", "FDR_lt_0.05", "interpretation"]); w.writeheader()
        for r, qq, p in zip(rows, q, ps):
            r = dict(r)
            r["FDR_BH"] = "NA" if qq is None else f"{qq:.6g}"; r["n_tests_in_FDR"] = m
            r["nominal_p_lt_0.05"] = p is not None and p < 0.05
            r["FDR_lt_0.05"] = qq is not None and qq < 0.05
            r["interpretation"] = ("not tested (no reachable targets)" if p is None else
                                   "FDR-significant proximity" if qq < 0.05 else
                                   "nominal p<0.05 only - NOT significant after FDR" if p < 0.05 else "not significant")
            w.writerow(r)
    back = list(csv.DictReader(open(out)))
    for a, b in zip(rows, back):
        for k in hdr: assert a[k] == b[k]
    print(rel, "tests:", m, "nominal<0.05:", sum(1 for p in ps if p is not None and p < 0.05),
          "FDR<0.05:", sum(1 for x in q if x is not None and x < 0.05),
          "min p:", min(p for p in ps if p is not None), "min FDR:", round(min(x for x in q if x is not None), 4))
