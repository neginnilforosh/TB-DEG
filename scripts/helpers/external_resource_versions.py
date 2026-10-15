#!/usr/bin/env python3
"""Record which external databases and web services each analysis used, with the best available version evidence.
Web services (STRING, iLINCS, DGIdb, Enrichr, KEGG REST) are not version-pinned by the scripts; for those the
record gives the date the outputs were first committed (the query happened on or before that date).
Reads existing files and `git log` only. Output: reproducibility/external_resource_versions.csv"""
import os, csv, gzip, subprocess, glob
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
def first_commit(path):
    try:
        out = subprocess.run(["git", "-C", ROOT, "log", "--diff-filter=A", "--format=%ad", "--date=short", "--", path],
                             capture_output=True, text=True, check=True).stdout.split()
        return out[-1] if out else "not committed"
    except Exception:
        return "git unavailable"
pkg = {}
pv = os.path.join(ROOT, "reproducibility", "package_versions.csv")
if os.path.exists(pv):
    pkg = {r["Package"]: r["Version"] for r in csv.DictReader(open(pv))}
def enrichr_libs(run):
    p = os.path.join(ROOT, run, "17_pathway_annotation", "GSE114192_DrugPathways_Enrichr_long.csv")
    if not os.path.exists(p): return ""
    return ";".join(sorted({r["Library"] for r in csv.DictReader(open(p))}))
def ctd_date(name):
    p = os.path.join(ROOT, "external_data", "CTD", name)
    if not os.path.exists(p): return "file not present locally"
    with gzip.open(p, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if "Report created" in line: return line.split(":", 1)[1].strip()
            if i > 60: break
    return "no header date"
rows = []
for run in ("GSE114192", "GSE114192/TB_DM"):
    rows += [
        dict(resource="STRING (REST API string-db.org)", run=run, version_evidence="version served on query date; not pinned",
             settings="species 9606; required_score 400 (scripts/wgcna_network/WGCNA_STRING_step4.R, step5)",
             outputs_first_committed=first_commit(f"{run}/12_string_ppi")),
        dict(resource="iLINCS (via drugfindR)", run=run, version_evidence=f"drugfindR {pkg.get('drugfindR', 'not recorded')}; iLINCS live service",
             settings="scripts/wgcna_network/WGCNA_iLINCS_signature_step6.R, step7", outputs_first_committed=first_commit(f"{run}/13_ilincs_signature")),
        dict(resource="DGIdb (GraphQL API dgidb.org)", run=run, version_evidence="live service; not pinned",
             settings="scripts/wgcna_network/WGCNA_drug_targets_step8.R", outputs_first_committed=first_commit(f"{run}/14_drug_targets")),
        dict(resource="Enrichr (per-drug target annotation)", run=run, version_evidence=f"libraries: {enrichr_libs(run)}; enrichR {pkg.get('enrichR', 'not recorded')}",
             settings="scripts/wgcna_network/WGCNA_Enrichr_drugs_step12.R", outputs_first_committed=first_commit(f"{run}/17_pathway_annotation")),
    ]
rows += [
    dict(resource="MSigDB (via msigdbr)", run="all GSEA", version_evidence=f"msigdbr {pkg.get('msigdbr', 'not recorded')}",
         settings="scripts/wgcna_network/WGCNA_GSEA_step11.R; seed 42", outputs_first_committed=first_commit("GSE114192/17_gsea")),
    dict(resource="GO / KEGG / Reactome over-representation", run="module enrichment",
         version_evidence=f"clusterProfiler {pkg.get('clusterProfiler', 'not recorded')}; ReactomePA {pkg.get('ReactomePA', 'not recorded')}; org.Hs.eg.db {pkg.get('org.Hs.eg.db', 'not recorded')}; KEGG via clusterProfiler online query",
         settings="scripts/wgcna_network/WGCNA_TBmodules_step2_enrichment.R", outputs_first_committed=first_commit("GSE114192/11_module_enrichment")),
    dict(resource="CTD chemical-gene interactions", run="CTD annotation", version_evidence=f"report created {ctd_date('CTD_chem_gene_ixns.csv.gz')}",
         settings="external_data/CTD/CTD_chem_gene_ixns.csv.gz (not tracked)", outputs_first_committed=first_commit("GSE114192/25_ctd_annotation")),
    dict(resource="CTD curated gene-disease associations", run="CTD annotation", version_evidence=f"report created {ctd_date('CTD_curated_genes_diseases.csv.gz')}",
         settings="external_data/CTD/CTD_curated_genes_diseases.csv.gz (not tracked)", outputs_first_committed=first_commit("GSE114192/25_ctd_annotation")),
    dict(resource="GEO series matrices and supplementary counts", run="all DEG runs", version_evidence=f"GEOquery {pkg.get('GEOquery', 'not recorded')}; downloaded at run time",
         settings="scripts/deg_pipeline/*_run.R", outputs_first_committed=first_commit("GSE114192/06_deg_results")),
    dict(resource="KEGG REST (metabolomics bridge)", run="ST000974", version_evidence="live service; not pinned",
         settings="scripts/wgcna_network/metabolomics_step19.R, step20", outputs_first_committed=first_commit("GSE114192/21_metabolomics")),
    dict(resource="SCP1749 macaque granuloma single-cell data", run="cell-type profiles", version_evidence="4-week matrix from the repository named in README",
         settings="scripts/wgcna_network/WGCNA_singlecell_step13.R", outputs_first_committed=first_commit("GSE114192/18_singlecell")),
]
p = os.path.join(ROOT, "reproducibility", "external_resource_versions.csv")
with open(p, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} rows -> {os.path.relpath(p, ROOT)}")
