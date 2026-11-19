#!/usr/bin/env python3
"""Supplementary figures and Additional file 2.
  Fig. S1  soft-threshold diagnostics of the three co-expression networks (from the *_WGCNA_soft_threshold.csv tables)
  Fig. S2  drug evidence cards of the 25 best-ranked compounds per context (from drug_evidence_cards_step26.R)
  Additional_file_2.pdf = Fig. S1 + Fig. S2 + Fig. S3 (Fig. S3 is assembled by assemble_figures.py; run it first)
Outputs in manuscript_figures/assembled/ (PDF, and PNG at 600 dpi)."""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
OUT = ROOT / "manuscript_figures" / "assembled"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42})
NODATE = {"CreationDate": None, "ModDate": None}       # no timestamps in the PDFs, so re-runs give identical files
CTX = {"TB": ("GSE114192", "GSE114192", "TB vs healthy controls", "TB vs healthy"),
       "TB_DM": ("GSE114192/TB_DM", "GSE114192", "TB with diabetes vs diabetes alone", "TB-DM vs diabetes"),
       "TB_HIV": ("GSE37250/TB_HIV", "GSE37250", "Active vs latent TB in people with HIV", "TB-HIV, active vs latent")}

# ---- Fig. S1 ----
fig, ax = plt.subplots(3, 2, figsize=(7.2, 8.6))
for i, (k, (run, acc, name, short)) in enumerate(CTX.items()):
    f = ROOT / run / "10_wgcna_results" / f"{acc}_WGCNA_soft_threshold.csv"
    if not f.exists(): sys.exit(f"Missing {f.relative_to(ROOT)}: run wgcna_network/WGCNA_TBmodules_step1.R or WGCNA_soft_threshold_record.R first.")
    d = pd.read_csv(f); used = int(d.power_used.iloc[0]); fallback = pd.isna(d.power_estimate.iloc[0])
    y = d.signed_R2
    a = ax[i, 0]; a.plot(d.Power, y, color="0.6", lw=0.6)
    for p, v in zip(d.Power, y):
        a.text(p, v, str(int(p)), ha="center", va="center", fontsize=7, color="#B2182B" if p == used else "black", fontweight="bold" if p == used else "normal")
    a.axhline(0.85, color="#B2182B", lw=0.6, ls="--"); a.set_ylim(min(-0.4, y.min() - 0.05), 1.0)
    a.set_xlabel("Soft-thresholding power"); a.set_ylabel("Scale-free fit, signed R²")
    a.set_title(f"{'ACE'[i]}  {short}: power {used}" + (" (fallback)" if fallback else ""), loc="left", fontsize=7.5, fontweight="bold")
    b = ax[i, 1]; b.plot(d.Power, d["mean.k."], color="0.6", lw=0.6)
    for p, v in zip(d.Power, d["mean.k."]):
        b.text(p, v, str(int(p)), ha="center", va="center", fontsize=7, color="#B2182B" if p == used else "black", fontweight="bold" if p == used else "normal")
    b.set_yscale("log"); b.set_xlabel("Soft-thresholding power"); b.set_ylabel("Mean connectivity")
    b.set_title("BDF"[i], loc="left", fontsize=8, fontweight="bold")
    for z in (a, b): z.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(OUT / "FigS1.pdf", metadata=NODATE); fig.savefig(OUT / "FigS1.png", dpi=600); plt.close(fig)
print("  FigS1: soft-threshold diagnostics")

# ---- Fig. S2 ----
def nice(s): return s.title() if s.isupper() and not re.search(r"\d|-", s) else s
cols = [("Evidence\nscore", "Evidence_score", 1), ("Rank\nrange", "Rank_range", 1), ("Best\nsimilarity,\nup genes", "Best_sim_UP", 1),
        ("Best\nsimilarity,\ndown genes", "Best_sim_DOWN", 1), ("Signatures\nreversing\nup genes (%)", "Share_reversing_UP", 1),
        ("Cell lines,\nstrong\nreversal (up)", "Cell_lines_reversing_UP", 1), ("Targets\nin\nnetwork", "Targets_in_network", 1),
        ("High-priority\ntargets", "High_priority_targets", 1.8), ("Cell-division\ntargets", "Cell_division_targets", 3.2), ("Proximity\nZ", "Proximity_Z", 1)]
XS = [2.4]
for c in cols: XS.append(XS[-1] + c[2])
with PdfPages(OUT / "FigS2.pdf", metadata=NODATE) as pdf:
    for j, (k, (run, acc, name, short)) in enumerate(CTX.items()):
        c = pd.read_csv(ROOT / run / "26_drug_evidence_cards" / f"{acc}_drug_evidence_cards.csv").sort_values("Evidence_rank"); n = len(c)
        fig, ax = plt.subplots(figsize=(11, 0.28 * n + 1.9)); ax.set_xlim(0, XS[-1]); ax.set_ylim(n, -1.9); ax.axis("off")
        ax.text(0, -1.6, f"{'ABC'[j]}  Drug evidence cards: {name} ({n} best-ranked compounds)", fontweight="bold", fontsize=9)
        ax.text(0, -0.3, "Compound", fontweight="bold", va="center")
        for ci, (lab, _, w) in enumerate(cols): ax.text(XS[ci] + w / 2, -0.45, lab, ha="center", va="center", fontweight="bold", fontsize=5.6)
        cm = plt.get_cmap("RdBu")
        for ri, r in enumerate(c.itertuples()):
            y = ri + 0.5
            both = bool(r.Reverses_both) and pd.notna(r.Share_reversing_UP) and r.Share_reversing_UP >= 0.6
            ax.text(0, y, f"{int(r.Evidence_rank)}. {nice(str(r.Drug))}", va="center", fontsize=7, fontweight="bold" if both else "normal")
            for ci, (lab, col, w) in enumerate(cols):
                v = getattr(r, col); x = XS[ci]
                txt = "" if (v is None or (isinstance(v, float) and np.isnan(v))) else v; face = "white"
                if col in ("Best_sim_UP", "Best_sim_DOWN") and txt != "": face = cm(0.5 + 0.5 * max(-1, min(1, float(v) / 0.8))); txt = f"{float(v):.2f}"
                elif col == "Share_reversing_UP" and txt != "": txt = f"{100 * float(v):.0f}"; face = plt.get_cmap("Greens")(float(v) * 0.8)
                elif col == "Evidence_score": txt = f"{float(v):.3f}"
                elif col == "Proximity_Z" and txt != "": txt = f"{float(v):.2f}"
                elif col == "Cell_division_targets" and txt != "": face = "#F4C7C3"; txt = str(v).replace(";", ", ")
                elif col == "High_priority_targets" and txt != "": face = "#D9E7F5"; txt = str(v).replace(";", ", ")
                elif col in ("Cell_lines_reversing_UP", "Targets_in_network") and txt != "": txt = str(int(v))
                ax.add_patch(plt.Rectangle((x + 0.03, ri + 0.06), w - 0.06, 0.88, facecolor=face, edgecolor="0.85", lw=0.3))
                s = str(txt); ax.text(x + w / 2, y, s, ha="center", va="center", fontsize=5.6 if len(s) > 12 else 6.5)
        ax.text(0, n + 0.9, "Bold compound names: strongly reversed both halves of the signature and the up-regulated half in at least 60% of LINCS signatures. "
                "Similarity is the iLINCS concordance (negative = reversal; strong at ≤ −0.321). Proximity is shown for information only; no compound passed FDR.", fontsize=6, va="top", wrap=True)
        pdf.savefig(fig, bbox_inches="tight"); fig.savefig(OUT / f"FigS2{'ABC'[j]}.png", dpi=600, bbox_inches="tight"); plt.close(fig)
print("  FigS2: drug evidence cards (3 pages)")

# ---- Additional file 2 ----
try:
    import pymupdf
except ImportError:
    import fitz as pymupdf
s3 = OUT / "FigS3.pdf"
if not s3.exists(): sys.exit("Missing manuscript_figures/assembled/FigS3.pdf: run scripts/assemble_figures.py first.")
doc = pymupdf.open()
for p in (OUT / "FigS1.pdf", OUT / "FigS2.pdf", s3): doc.insert_pdf(pymupdf.open(p))
doc.set_metadata({}); doc.save(OUT / "Additional_file_2.pdf", garbage=3, no_new_id=True); print(f"  Additional_file_2.pdf: {len(doc)} pages (Fig. S1, Fig. S2 A-C, Fig. S3)")
