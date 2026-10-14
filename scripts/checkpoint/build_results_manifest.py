#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: results manifest. One row per result file, with a status label.
Status rules (path-based) are documented in the STATUS_RULES list. Raw inputs and git-ignored untarred copies are skipped.
Output: reproducibility/checkpoint_2026-10-14/RESULTS_MANIFEST.csv"""
import os, re, csv, subprocess, datetime
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CK = "checkpoint_2026-10-14"
base = set()
try:
    rows = open(os.path.join(ROOT, "reproducibility", CK, "baseline_sha256.tsv")).read().splitlines()
    hdr = rows[0].split("\t"); pi = hdr.index("path"); base = {r.split("\t")[pi] for r in rows[1:]}
except Exception: pass
STATUS_RULES = [
    (r"(^|/)(00_checkpoint_2026-10-14|diagnostics_2026-10-14|checkpoint_2026-10-14)/|_annotated\.csv$|_withFDR\.csv$|_from_workspace\.csv$|module_merging_summary|_checkpoint_2026-10-14\.md$|/STATUS\.md$|PRELIMINARY_STATUS|00_contrast_plan|00_dataset_decisions|00_sample_manifest|TB_DM/00_design/|^GSE50834/", "NEW at checkpoint (diagnostic/annotation; no inferential rerun)"),
    (r"^GSE222001/", "INVALID until rerun: Healthy/Latent labels swapped (D11)"),
    (r"^GSE161829/", "AFFECTED: 7/85 samples mislabelled by positional column mapping; rerun pending (D11)"),
    (r"^GSE114192/(23_external_validation|24_ml_validation)/", "AFFECTED: uses GSE222001/GSE161829 results; rerun pending (D11, D12)"),
    (r"^GSE114192/17_gsea/GSE114192_(DM_only_vs_TB_DM|Healthy_Control_vs_DM_only|Healthy_Control_vs_TB_DM|TB_Only_vs_TB_DM)", "PRELIMINARY TB-DM file stored in reference folder (copy in TB_DM/17_gsea)"),
    (r"^GSE114192/TB_DM/", "PRELIMINARY (TB-DM; pending D1-D4, D13)"),
    (r"^GSE114192/2[12]_", "SUPPORTING layer; ST000974 case definition unverified (D7)"),
    (r"^GSE114192/", "REFERENCE (TB vs HC; existing output, not rerun)"),
    (r"^external_results/GSE165709/", "SUPPORTING validation layer (GSE165709)"),
    (r"^GSE99374/", "PILOT / LTBI validation (CD8 T cells); labels verified"),
    (r"^GSE229020/", "miRNA layer (labels from sample-name prefixes); not audited at checkpoint"),
]
EXT = {".csv", ".tsv", ".txt", ".rds", ".png", ".pdf", ".tif", ".tiff", ".svg", ".md", ".xlsx"}
out = []
for top in sorted(os.listdir(ROOT)):
    if not (re.match(r"GSE\d+", top) or top == "external_results"): continue
    for dp, dn, fn in os.walk(os.path.join(ROOT, top)):
        rel_dp = os.path.relpath(dp, ROOT)
        if re.search(r"(^|/)(01_raw_counts|untarred)(/|$)", rel_dp): continue
        for f in sorted(fn):
            if f.startswith(".") or os.path.splitext(f)[1].lower() not in EXT: continue
            rel = os.path.join(rel_dp, f)
            st = next((s for p, s in STATUS_RULES if re.search(p, rel)), "UNCLASSIFIED")
            parts = rel.split("/")
            out.append(dict(path=rel, dataset=parts[0], step_folder=parts[-2] if len(parts) > 1 else "",
                            type=os.path.splitext(f)[1].lstrip(".").lower(), bytes=os.path.getsize(os.path.join(ROOT, rel)),
                            in_baseline=rel in base, status=st))
p = os.path.join(ROOT, "reproducibility", CK, "RESULTS_MANIFEST.csv")
with open(p, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
from collections import Counter
for k, v in Counter(r["status"] for r in out).most_common(): print(v, k)
print("total", len(out), "| new (not in baseline):", sum(not r["in_baseline"] for r in out))
