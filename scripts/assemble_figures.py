#!/usr/bin/env python3
"""Assemble the manuscript's multi-panel figures from existing panels.

Panels come from the pipeline's figure folders and from manuscript_figures/ (made by manuscript_figures.R).
Each figure is a stack of rows; panels in a row are scaled to a common height and the row to the page width
(180 mm, double column, 300 dpi). Panel letters are added in bold. A missing panel is replaced by a grey box
naming the file, so the layout can be checked before every panel exists.
Run from the repository root:  python3 scripts/assemble_figures.py
Output: manuscript_figures/assembled/Fig2.png ... Fig9.png, FigS1.png (and PDF versions).
"""
import glob, sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
MF = ROOT / "manuscript_figures"; OUT = MF / "assembled"; OUT.mkdir(parents=True, exist_ok=True)
DPI, WIDTH_MM, GAP = 300, 180, 30
MAX_ROW = {1: 760, 2: 1000}          # px: tallest allowed row for one panel / several panels (~64 and 85 mm)
W = round(WIDTH_MM / 25.4 * DPI)
try:
    FONT = ImageFont.truetype("DejaVuSans-Bold.ttf", 54)
except OSError:
    try: FONT = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 54)
    except OSError: FONT = ImageFont.load_default()

def find(pattern):
    hits = sorted(glob.glob(str(ROOT / pattern)))
    return Path(hits[0]) if hits else None

def load(pattern):
    p = find(pattern)
    if p: return Image.open(p).convert("RGB"), p.name
    im = Image.new("RGB", (900, 600), (235, 235, 235)); d = ImageDraw.Draw(im)
    d.text((30, 280), f"missing:\n{pattern}", fill=(90, 90, 90)); return im, None

def row(panels, letters):
    ims = [load(p)[0] for p in panels]
    h = 1000; ims = [im.resize((round(im.width * h / im.height), h), Image.LANCZOS) for im in ims]
    total = sum(im.width for im in ims) + GAP * (len(ims) - 1)
    scale = (W - 0) / total
    ims = [im.resize((max(1, round(im.width * scale)), max(1, round(h * scale))), Image.LANCZOS) for im in ims]
    limit = MAX_ROW[1 if len(ims) == 1 else 2]
    if max(im.height for im in ims) > limit:                       # keep rows short enough for one page
        f = limit / max(im.height for im in ims); ims = [im.resize((round(im.width * f), round(im.height * f)), Image.LANCZOS) for im in ims]
    width = sum(im.width for im in ims) + round(GAP * scale) * (len(ims) - 1)
    rh = max(im.height for im in ims); canvas = Image.new("RGB", (W, rh + 10), "white"); d = ImageDraw.Draw(canvas); x = (W - width) // 2
    for im, L in zip(ims, letters):
        canvas.paste(im, (x, 10))
        if L: d.text((x + 4, 0), L, fill="black", font=FONT)
        x += im.width + round(GAP * scale)
    return canvas

def figure(name, rows):
    parts = [row(p, l) for p, l in rows]
    H = sum(p.height for p in parts) + GAP * (len(parts) - 1)
    fig = Image.new("RGB", (W, H), "white"); y = 0
    for p in parts: fig.paste(p, (0, y)); y += p.height + GAP
    fig.save(OUT / f"{name}.png", dpi=(DPI, DPI)); fig.save(OUT / f"{name}.pdf", resolution=DPI)
    missing = [pp for p, _ in rows for pp in p if not find(pp)]
    print(f"  {name}: {fig.width}x{fig.height}px" + (f" | missing: {', '.join(missing)}" if missing else ""))

R, M = "GSE114192", "manuscript_figures"
figure("Fig2", [([f"{M}/Fig2A_volcano_TB.png", f"{R}/figures/{R}_HealthyControl_vs_TBOnly_DEG_heatmap.png"], ["A", "B"]),
                ([f"{M}/Fig2C_module_trait_TB.png"], ["C"]),
                ([f"{M}/Fig2D_dendrogram_TB.png"], ["D"])])
figure("Fig3", [([f"{R}/figures_manuscript/{R}_E_network_labelled.png", f"{M}/Fig3B_degree_distribution_TB.png"], ["A", "B"]),
                ([f"{R}/figures_manuscript/{R}_E_network_communities.png", f"{M}/Fig3D_score_agreement_TB.png"], ["C", "D"])])
figure("Fig4", [([f"{M}/Fig4A_upset_evidence_families_TB.png"], ["A"]),
                ([f"{M}/Fig4B_top_genes_evidence_TB.png", f"{R}/figures_manuscript/{R}_F_GSEA_dotplot.png"], ["B", "C"]),
                ([f"{M}/Fig4D_GSEA_INTERFERON_GAMMA_RESPONSE_TB.png", f"{M}/Fig4D_GSEA_INTERFERON_ALPHA_RESPONSE_TB.png"], ["D", ""])])
figure("Fig5", [([f"{R}/24_ml_validation/{R}_ml_validation_ROC.png"], ["A"]),
                ([f"{R}/23_external_validation/{R}_external_validation_heatmap.png"], ["B"])])
figure("Fig6", [([f"{R}/figures_manuscript/{R}_G_iLINCS_top_reversal.png", f"{R}/figures_manuscript/{R}_G_drug_module_heatmap.png"], ["A", "B"]),
                ([f"{R}/figures_manuscript/{R}_G_drug_target_network.png", f"{M}/Fig6D_reversal_vs_proximity_TB.png"], ["C", "D"])])
figure("Fig7", [([f"{M}/Fig7A_TBHIV_volcano.png", f"{M}/Fig7B_TBHIV_PCA.png"], ["A", "B"]),
                ([f"{M}/Fig7C_interferon_comorbidity_models.png"], ["C"])])
figure("Fig8", [([f"{M}/Fig8A_core_genes_across_contexts.png", f"{M}/Fig8B_hallmark_across_contexts.png"], ["A", "B"]),
                ([f"{M}/Fig8C_module_correspondence_TB_vs_TB_DM.png", f"{M}/Fig8C_module_correspondence_TB_vs_TB_HIV.png"], ["C", ""])])
figure("Fig9", [([f"{M}/Fig9A_macrophage_infection_response.png"], ["A"]),
                ([f"{M}/Fig9B_kynurenine_tryptophan.png"], ["B"]),
                ([f"{M}/Fig9C_celltype_origin_TB.png"], ["C"])])
pages = []
for ctx, tag in [("TB_DM", "TB with diabetes"), ("TB_HIV", "TB with HIV")]:
    letters = iter("ABCDEFGH")
    rows = [([f"{M}/Fig2C_module_trait_{ctx}.png", f"{M}/Fig3B_degree_distribution_{ctx}.png"], None),
            ([f"{M}/Fig2D_dendrogram_{ctx}.png"], None),
            ([f"{M}/Fig4A_upset_evidence_families_{ctx}.png", f"{M}/Fig4B_top_genes_evidence_{ctx}.png"], None),
            ([f"{M}/Fig9C_celltype_origin_{ctx}.png"], None)]
    rows = [(p, [next(letters) for _ in p]) for p, _ in rows]
    figure(f"FigS1_{ctx}", rows); pages.append(Image.open(OUT / f"FigS1_{ctx}.png").convert("RGB"))
pages[0].save(OUT / "FigS1.pdf", save_all=True, append_images=pages[1:], resolution=DPI)
print("  FigS1.pdf: 2 pages (TB with diabetes, TB with HIV)")
print(">>> DONE:", OUT)
