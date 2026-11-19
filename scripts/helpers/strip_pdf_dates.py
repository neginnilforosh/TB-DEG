#!/usr/bin/env python3
"""Remove creation and modification timestamps from every PDF in the repository.
R, matplotlib and Pillow write the time a PDF was made into its metadata; this step clears those two fields
(other metadata such as the title is kept) and also removes the XMP metadata stream, so the public files carry
no dates and a re-run that changes nothing gives an identical file. Page content is not touched.
Run after regenerating figures:   python3 scripts/helpers/strip_pdf_dates.py              (every PDF)
                                  python3 scripts/helpers/strip_pdf_dates.py a.pdf b.pdf  (only these)
Needs PyMuPDF (pip3 install pymupdf)."""
import sys, subprocess
from pathlib import Path
try:
    import pymupdf
except ImportError:
    import fitz as pymupdf

here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
if len(sys.argv) > 1:                                   # only the PDFs given on the command line
    files = [Path(a).resolve() for a in sys.argv[1:]]
else:
    files = sorted(p for p in ROOT.rglob("*.pdf") if ".git" not in p.parts and "external_data" not in p.parts)
changed = 0
for f in files:
    doc = pymupdf.open(f)
    md = dict(doc.metadata or {})
    if not (md.get("creationDate") or md.get("modDate")) and not doc.get_xml_metadata():
        doc.close(); continue
    md["creationDate"] = ""; md["modDate"] = ""
    doc.set_metadata({k: v for k, v in md.items() if k in ("title", "author", "subject", "keywords", "creator", "producer", "creationDate", "modDate")})
    doc.del_xml_metadata()
    tmp = f.with_suffix(".tmp.pdf"); doc.save(tmp, garbage=3, deflate=True, no_new_id=True); doc.close()
    tmp.replace(f); changed += 1
print(f"PDFs checked: {len(files)} | timestamps removed: {changed}")
