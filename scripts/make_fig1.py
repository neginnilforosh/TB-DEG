#!/usr/bin/env python3
"""Figure 1: study workflow. Run: python3 scripts/make_fig1.py  -> manuscript_figures/Fig1_workflow.png/.pdf"""
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), here)
OUT = ROOT / "manuscript_figures"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7})
INK, EDGE, BLUE, RED, GREY = "#1f1f1f", "#7a7a7a", "#2166AC", "#B2182B", "#f2f2f2"
fig = plt.figure(figsize=(180 / 25.4, 128 / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 180); ax.set_ylim(0, 128); ax.axis("off")
def box(x, y, w, h, title, sub="", fc="white", ec=EDGE, bold=True, ts=7, ss=6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.6", fc=fc, ec=ec, lw=0.8))
    if sub:
        ax.text(x + w / 2, y + h * 0.64, title, ha="center", va="center", fontsize=ts, fontweight="bold" if bold else "normal", color=INK)
        ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center", fontsize=ss, color="#444444", linespacing=1.15)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=ts, fontweight="bold" if bold else "normal", color=INK, linespacing=1.15)
def arrow(x1, y1, x2, y2, c=EDGE):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=c))
def band(y, h, label):
    ax.add_patch(FancyBboxPatch((2, y), 176, h, boxstyle="round,pad=0.3,rounding_size=2", fc=GREY, ec="none"))
    ax.text(4, y + h - 2.6, label, fontsize=6.5, fontweight="bold", color="#555555", va="top")
# 1 inputs
band(103, 23, "1  Clinical contexts, each analysed independently")
for i, (t, s) in enumerate([("TB vs healthy controls", "GSE114192, whole blood RNA-seq\n46 vs 36"),
                            ("TB with diabetes vs diabetes", "GSE114192, whole blood RNA-seq\n61 vs 52"),
                            ("Active vs latent TB in HIV", "GSE37250, whole blood arrays\n98 vs 84")]):
    box(14 + i * 56, 105.5, 48, 13, t, s, ec=BLUE)
# 2 pipeline
band(57, 43, "2  Same workflow in every context")
steps = [("Differential\nexpression", "DESeq2 / limma"), ("Co-expression\nmodules", "WGCNA"), ("Protein\nnetwork", "STRING, centrality,\nRWR, diffusion"),
         ("Evidence-family\nranking", "5 families,\nstability tests"), ("Pathways and\nclassifier", "GSEA, nested CV")]
for i, (t, s) in enumerate(steps):
    x = 6 + i * 34.5; box(x, 78, 29, 16, t, s)
    if i: arrow(x - 5.3, 86, x - 0.3, 86)
drugs = [("Signature\nreversal", "iLINCS / LINCS"), ("Targets and\nproximity", "DGIdb, random\ntarget sets"), ("Drug\nsimilarity", "chemistry, targets,\nnetwork, pathways"), ("Evidence\ncards", "drug shortlist")]
for i, (t, s) in enumerate(drugs):
    x = 23 + i * 34.5; box(x, 60, 29, 14.5, t, s)
    if i: arrow(x - 5.3, 67.2, x - 0.3, 67.2)
arrow(20.5, 78, 32, 74.8)
# 3 validation
band(29, 25, "3  Tested against data the ranking did not use")
val = [("Independent\ncohorts", "CD4+, NK, CD8+;\n2 HIV cohorts"), ("Infected alveolar\nmacrophages", "RNA-seq,\nATAC-seq"), ("Single-cell\nlung granulomas", "cell-type origin"),
       ("Metabolomics", "GC6-74, before\ndiagnosis"), ("Curated\nassociations", "CTD")]
for i, (t, s) in enumerate(val): box(6 + i * 34.5, 31.5, 29, 15, t, s)
# 4 comparison and outputs
band(2, 24, "4  Comparison across contexts")
box(8, 5, 52, 13, "Shared host response", "genes, modules and pathways\nin all three contexts", ec=RED)
box(64, 5, 52, 13, "Effect of the comorbidity", "TB x diabetes and TB x HIV\ninteraction models", ec=RED)
box(120, 5, 52, 13, "Drug hypotheses", "for experimental testing", ec=RED)
for x in (38, 94, 150): arrow(x, 103, x, 98.6, c=BLUE)
arrow(90, 57.2, 90, 52.5); arrow(90, 29.2, 90, 22.5)
for ext in ("png", "pdf"): fig.savefig(OUT / f"Fig1_workflow.{ext}", dpi=600)
print("saved", OUT / "Fig1_workflow.png")
