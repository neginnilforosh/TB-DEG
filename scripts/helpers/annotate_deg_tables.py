#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: add Symbol / Direction / Comparison columns to existing DESeq2 tables.

NO REANALYSIS. Reads existing CSVs, appends annotation columns, writes *_annotated.csv next to them.
The original numeric columns are copied as text and verified to be identical after writing.
Symbol source priority: (1) GSE114192 FINAL_GeneTable (pipeline's own org.Hs.eg.db 3.22 mapping),
(2) GSE114192 TB_DM FINAL_GeneTable, (3) GSE165709 RNA results table. Unmapped -> empty.
"""
import csv, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
G = os.path.join(ROOT, "GSE114192")

def load_map(path, ecol, scol):
    m = {}
    if not os.path.exists(path): return m
    for r in csv.DictReader(open(path)):
        e, s = r.get(ecol, ""), r.get(scol, "")
        if e.startswith("ENSG") and s and not s.startswith("ENSG"): m.setdefault(e, s)
    return m
sources = [
    ("GSE114192_FINAL_GeneTable", load_map(os.path.join(G, "16_final_tables", "GSE114192_FINAL_GeneTable.csv"), "Ensembl", "Gene")),
    ("GSE114192_TB_DM_FINAL_GeneTable", load_map(os.path.join(G, "TB_DM", "16_final_tables", "GSE114192_FINAL_GeneTable.csv"), "Ensembl", "Gene")),
]
# GSE165709 RNA table: find its Ensembl and symbol columns
p165 = os.path.join(ROOT, "external_results", "GSE165709", "RNA", "GSE165709_RNA_INF_vs_NEG_in_HC.csv")
if os.path.exists(p165):
    h = next(csv.reader(open(p165)))
    ecol = next((c for c in h if c.lower() in ("ensembl", "ensembl_id", "gene_id")), None)
    scol = next((c for c in h if c.lower() in ("gene", "symbol", "gene_symbol")), None)
    if ecol and scol: sources.append(("GSE165709_RNA_table", load_map(p165, ecol, scol)))

def symbol(e):
    for name, m in sources:
        if e in m: return m[e], name
    return "", "unmapped"

JOBS = [  # path, comparison label, case, reference, positive-LFC meaning, lfc type
  ("06_deg_results/GSE114192_HealthyControl_vs_TBOnly_DEG.csv", "TB_Only vs Healthy_Control", "TB_Only", "Healthy_Control", "apeglm-shrunk", "reference (completed)"),
  ("07_significant_degs/GSE114192_HealthyControl_vs_TBOnly_sigDEGs.csv", "TB_Only vs Healthy_Control", "TB_Only", "Healthy_Control", "apeglm-shrunk", "reference (completed)"),
  ("TB_DM/07_significant_degs/GSE114192_DM_only_vs_TB_DM_sigDEGs.csv", "TB_DM vs DM_only", "TB_DM", "DM_only", "apeglm-shrunk", "TB-DM PRELIMINARY (~ group, site-confounded)"),
  ("TB_DM/07_significant_degs/GSE114192_TB_Only_vs_TB_DM_sigDEGs.csv", "TB_DM vs TB_Only", "TB_DM", "TB_Only", "unshrunk MLE", "TB-DM PRELIMINARY (~ group)"),
  ("TB_DM/07_significant_degs/GSE114192_Healthy_Control_vs_TB_DM_sigDEGs.csv", "TB_DM vs Healthy_Control", "TB_DM", "Healthy_Control", "unshrunk MLE", "TB-DM PRELIMINARY (~ group, site-confounded)"),
  ("TB_DM/07_significant_degs/GSE114192_Healthy_Control_vs_DM_only_sigDEGs.csv", "DM_only vs Healthy_Control", "DM_only", "Healthy_Control", "unshrunk MLE", "TB-DM PRELIMINARY (~ group)"),
]
for rel, comp, case, ref, lfc_type, status in JOBS:
    src = os.path.join(G, rel)
    if not os.path.exists(src): print("missing", rel); continue
    rows = list(csv.DictReader(open(src))); hdr = list(rows[0].keys()) if rows else []
    out = src[:-4] + "_annotated.csv"
    extra = ["gene_symbol", "symbol_source", "direction", "padj_status", "comparison", "case_group",
             "reference_group", "lfc_type", "test", "status"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr + extra); w.writeheader()
        for r in rows:
            s, sname = symbol(r["gene"])
            lfc = r.get("log2FoldChange", "NA")
            try: d = f"Up in {case}" if float(lfc) > 0 else (f"Down in {case}" if float(lfc) < 0 else "zero")
            except ValueError: d = "NA"
            pa = r.get("padj", "NA")
            pst = "not tested (padj NA: independent filtering/outlier)" if pa in ("NA", "") else "tested"
            r.update(gene_symbol=s, symbol_source=sname, direction=d, padj_status=pst, comparison=comp,
                     case_group=case, reference_group=ref, lfc_type=lfc_type, test="DESeq2 Wald test, BH-adjusted",
                     status=status)
            w.writerow(r)
    # verify original columns unchanged
    back = list(csv.DictReader(open(out)))
    assert len(back) == len(rows)
    for a, b in zip(rows, back):
        for k in hdr: assert a[k] == b[k], (rel, k)
    nm = sum(1 for b in back if b["gene_symbol"])
    print(f"{rel}: {len(rows)} rows, original columns identical, symbols mapped {nm}/{len(rows)}")
