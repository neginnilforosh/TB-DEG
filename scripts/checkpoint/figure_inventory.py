#!/usr/bin/env python3
"""Pre-TB-HIV checkpoint: inventory of raster/vector figures with pixel size and embedded DPI (PNG pHYs / TIFF tags).
Checks that final raster figures are exported at 600 DPI. Output: reproducibility/FIGURE_INVENTORY.csv"""
import os, re, csv
from PIL import Image
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
rows = []
for top in sorted(os.listdir(ROOT)):
    if not (re.match(r"GSE\d+", top) or top in ("external_results", "figures")): continue
    for dp, _, fn in os.walk(os.path.join(ROOT, top)):
        for f in sorted(fn):
            ext = os.path.splitext(f)[1].lower()
            if ext not in (".png", ".tif", ".tiff", ".pdf", ".svg"): continue
            rel = os.path.relpath(os.path.join(dp, f), ROOT); stem = os.path.splitext(rel)[0]
            r = dict(path=rel, type=ext[1:], width_px="", height_px="", dpi_embedded="", width_in_at_dpi="",
                     meets_600dpi="", has_pdf_or_svg=os.path.exists(stem + ".pdf") or os.path.exists(stem + ".svg"),
                     has_png_or_tiff=any(os.path.exists(stem + e) for e in (".png", ".tif", ".tiff")))
            if ext in (".png", ".tif", ".tiff"):
                try:
                    im = Image.open(os.path.join(ROOT, rel)); r["width_px"], r["height_px"] = im.size
                    d = im.info.get("dpi")
                    if d:
                        r["dpi_embedded"] = round(float(d[0])); r["width_in_at_dpi"] = round(im.size[0] / float(d[0]), 2)
                        r["meets_600dpi"] = round(float(d[0])) >= 600
                    else: r["meets_600dpi"] = "unknown (no DPI tag)"
                except Exception as e: r["dpi_embedded"] = f"error: {e}"
            rows.append(r)
rows.sort(key=lambda r: r["path"])          # deterministic order across file systems
p = os.path.join(ROOT, "reproducibility", "FIGURE_INVENTORY.csv")
with open(p, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
import collections
ras = [r for r in rows if r["type"] in ("png", "tif", "tiff")]
print("figures:", len(rows), "raster:", len(ras), "vector:", len(rows) - len(ras))
print("DPI distribution:", collections.Counter(str(r["dpi_embedded"]) for r in ras).most_common())
print("raster without vector companion:", sum(not r["has_pdf_or_svg"] for r in ras))
ref = [r for r in ras if r["path"].startswith("GSE114192/") and "/TB_DM/" not in r["path"]]
print("reference raster:", len(ref), "meeting 600:", sum(r["meets_600dpi"] is True for r in ref))
