#!/usr/bin/env python3
"""Drug similarity and complementarity among the candidate compounds of each clinical context.

For the 25 candidates of each context (the drug evidence cards), this script characterises how the compounds relate
to each other, using only outputs the pipeline has already produced plus chemical structures from PubChem:

  1. compound library  PubChem CID, SMILES and InChIKey per compound; structures standardised with RDKit (largest
                       fragment, neutralised); different salts or forms of one substance are flagged as duplicates;
                       compounds without a structure are kept and flagged, never given a similarity of zero
  2. target similarity Jaccard index of documented DGIdb targets (all targets, and targets in the disease network)
  3. network           shortest-path distances between the two target sets in the context's STRING network:
                       symmetric nearest-target distance (primary), mean pairwise distance and separation
                       s_AB = d_AB - (d_AA + d_BB) / 2 (when both sets have >= 2 targets), each pair compared with
                       degree-matched random target sets of the same sizes; disconnected targets are recorded;
                       overlap of the network communities and co-expression modules the targets fall in
  4. chemistry         Morgan fingerprints (radius 2, 2048 bits), Tanimoto similarity, hierarchical clustering with
                       silhouette for 2-4 clusters (no clusters reported when the structure is weak)
  5. pathways          Jaccard index of the significant Enrichr terms (GO BP, KEGG, Reactome, WikiPathways) of the
                       targets of each compound
  6. exploratory pairs every dimension reported separately, no combined score; each pair carries an evidence label
                       from the support of its two compounds (a pair of weakly supported compounds is never upgraded)
  7. sensitivity       agreement between distance summaries, excluding weakly supported compounds or compounds with
                       few network targets, and stability of chemical clusters across fingerprint settings
Nothing here predicts synergy or efficacy; distant targets are a possible sign of complementary mechanisms only.

Run from the repository root:  python3 scripts/drug_complementarity.py   (add --offline to skip PubChem)
Outputs in drug_complementarity/<context>/ and drug_complementarity/drug_library.csv.
"""
import argparse, json, re, ssl, sys, time, urllib.error, urllib.parse, urllib.request
from itertools import combinations
from pathlib import Path
import warnings; import numpy as np, pandas as pd
warnings.filterwarnings("ignore", category=RuntimeWarning)

CONTEXTS = {"TB": ("GSE114192", "GSE114192", "TB vs healthy"),
            "TB_DM": ("GSE114192/TB_DM", "GSE114192", "TB with diabetes"),
            "TB_HIV": ("GSE37250/TB_HIV", "GSE37250", "TB with HIV")}
N_NULL, SEED, STRONG_SIM = 500, 42, -0.321
CELL_DIVISION = re.compile(r"^(CDK[0-9]+|PLK[0-9]|AURK[ABC]|TOP1|TOP2A|TOP2B|TYMS|RRM[12]|DHFR|KIF11|WEE1|CHEK[12]|TTK|BUB1B?|CDC25[ABC]|HDAC[1-3])$")
ENRICHR_LIBS = ("GO_Biological_Process", "KEGG", "Reactome", "WikiPathway")
SALT = re.compile(r"\b(HYDROCHLORIDE|DIHYDROCHLORIDE|HCL|MESYLATE|MALEATE|SODIUM|POTASSIUM|CITRATE|TARTRATE|SULFATE|SULPHATE|ACETATE|FUMARATE|SUCCINATE|BESYLATE|TOSYLATE|PHOSPHATE|HYDRATE|MONOHYDRATE|DIHYDRATE|BROMIDE|CHLORIDE|NITRATE|LACTATE|CALCIUM)\b")

ap = argparse.ArgumentParser(); ap.add_argument("--offline", action="store_true", help="use the cached library only")
ap.add_argument("--refresh", action="store_true", help="look up every compound again")
args = ap.parse_args()
here = Path(__file__).resolve().parent
ROOT = next((p for p in [here, *here.parents] if (p / "GSE114192").is_dir()), None)
if ROOT is None: sys.exit("Run from inside the repository.")
OUT = ROOT / "drug_complementarity"; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(SEED)

def rd(*p):
    f = ROOT.joinpath(*p)
    if not f.exists(): print(f"  (missing: {f})"); return None
    return pd.read_csv(f)

def group_of(r):
    consistent = bool(r.Reverses_both) and (r.Share_reversing_UP if pd.notna(r.Share_reversing_UP) else 0) >= 0.6
    cd = isinstance(r.Cell_division_targets, str) and r.Cell_division_targets.strip() != ""
    return "1 Recommended" if consistent and not cd else ("2 Caution: possible antiproliferative" if consistent else "3 Lower support")

# ---------------------------------------------------------------- 1. compound library (shared by all contexts)
cards = {}
for k, (run, acc, _) in CONTEXTS.items():
    c = rd(run, "26_drug_evidence_cards", f"{acc}_drug_evidence_cards.csv")
    if c is not None: c["Group"] = c.apply(group_of, axis=1); cards[k] = c
names = sorted({d for c in cards.values() for d in c.Drug})
lib_f = OUT / "drug_library.csv"
lib = pd.read_csv(lib_f) if lib_f.exists() else pd.DataFrame(columns=["Drug"])
if args.refresh: lib = pd.DataFrame(columns=["Drug"])
retry = set(lib.Drug[~lib.get("Structure_status", pd.Series(dtype=str)).isin(["found", "not found in PubChem"])]) if len(lib) else set()
lib = lib[~lib.Drug.isin(retry)]                                  # failed lookups are retried
todo = [n for n in names if n not in set(lib.Drug)]
try:
    import certifi; SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CTX = ssl.create_default_context()
FIRST_ERROR = []

def pubchem(name):
    props = ["SMILES,InChIKey,IUPACName", "IsomericSMILES,InChIKey,IUPACName", "CanonicalSMILES,InChIKey,IUPACName"]
    for q in dict.fromkeys([name, SALT.sub("", name).strip()]):
        if not q: continue
        for pr in props:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(q)}/property/{pr}/JSON"
            try:
                with urllib.request.urlopen(url, timeout=30, context=SSL_CTX) as r: p = json.load(r)["PropertyTable"]["Properties"][0]
                smi = p.get("SMILES") or p.get("IsomericSMILES") or p.get("CanonicalSMILES")
                return dict(Drug=name, Query=q, CID=p.get("CID"), SMILES_PubChem=smi, InChIKey_PubChem=p.get("InChIKey"), Structure_status="found")
            except urllib.error.HTTPError as e:
                if e.code == 400: continue          # property name not supported: try the next form
                break                                # 404: name not found, try the stripped name
            except Exception as e:
                if not FIRST_ERROR: FIRST_ERROR.append(repr(e))
                return dict(Drug=name, Structure_status="lookup failed")
            finally: time.sleep(0.25)
    return dict(Drug=name, Structure_status="not found in PubChem")

if todo and not args.offline:
    print(f">>> PubChem lookup for {len(todo)} compounds")
    lib = pd.concat([lib, pd.DataFrame([pubchem(n) for n in todo])], ignore_index=True)
    if FIRST_ERROR:
        print("  !! PubChem lookups failed:", FIRST_ERROR[0])
        if "CERTIFICATE" in FIRST_ERROR[0].upper():
            print("  -> Python cannot verify HTTPS certificates. Run:  pip3 install certifi   (or, on macOS, open")
            print("     /Applications/Python 3.x/Install Certificates.command), then rerun this script.")
elif todo:
    lib = pd.concat([lib, pd.DataFrame({"Drug": todo, "Structure_status": "not looked up (offline)"})], ignore_index=True)

for col in ["Query", "CID", "SMILES_PubChem", "InChIKey_PubChem", "Structure_status"]:
    if col not in lib.columns: lib[col] = np.nan
print(">>> structure lookup:", ", ".join(f"{k} {v}" for k, v in lib.Structure_status.value_counts().items()))
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
chooser, uncharger = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()
def standardise(smi):
    if not isinstance(smi, str): return None
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return uncharger.uncharge(chooser.choose(m))
std = {}
for i, r in lib.iterrows():
    m = standardise(r.get("SMILES_PubChem"))
    std[r.Drug] = m
    lib.loc[i, "SMILES_standardised"] = Chem.MolToSmiles(m) if m else None
    lib.loc[i, "InChIKey_standardised"] = Chem.MolToInchiKey(m) if m else None
    if r.get("Structure_status") == "found" and m is None: lib.loc[i, "Structure_status"] = "invalid structure"
lib["Skeleton"] = lib.InChIKey_standardised.str[:14]
dup = lib.dropna(subset=["Skeleton"]).groupby("Skeleton").Drug.apply(lambda s: "; ".join(sorted(s)))
lib["Same_substance_as"] = lib.apply(lambda r: "; ".join(d for d in dup.get(r.Skeleton, "").split("; ") if d and d != r.Drug) if pd.notna(r.Skeleton) else "", axis=1)
lib.to_csv(lib_f, index=False)
print(f">>> library: {len(lib)} compounds | with structure: {lib.SMILES_standardised.notna().sum()} | flagged as the same substance: {(lib.Same_substance_as != '').sum()}")

def fps(radius=2, bits=2048):
    g = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    return {d: g.GetFingerprint(m) for d, m in std.items() if m is not None}

# ---------------------------------------------------------------- per context
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns, networkx as nx
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score, adjusted_rand_score

def jacc(a, b):
    if not a or not b: return np.nan
    return len(a & b) / len(a | b)

def save(fig, path):
    for ext in ("png", "pdf"): fig.savefig(f"{path}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)

summary = []
for k, (run, acc, label) in CONTEXTS.items():
    if k not in cards: continue
    print(f"\n========== {k} ({label}) ==========")
    o = OUT / k; o.mkdir(exist_ok=True)
    c = cards[k].copy(); drugs = list(c.Drug); grp = dict(zip(c.Drug, c.Group))
    # targets
    tall = rd(run, "14_drug_targets", f"{acc}_drug_targets_ALL.csv")
    T_all = tall.groupby("Drug").Target_Gene.apply(lambda s: set(s.dropna())).to_dict() if tall is not None else {}
    # network and distances
    net = rd(run, "12_string_ppi", f"{acc}_COMBINED_STRING_network.csv")
    G = nx.from_pandas_edgelist(net, "preferredName_A", "preferredName_B")
    nodes = list(G.nodes); idx = {n: i for i, n in enumerate(nodes)}
    D = np.full((len(nodes), len(nodes)), np.inf)
    for s, lens in nx.all_pairs_shortest_path_length(G):
        for t, l in lens.items(): D[idx[s], idx[t]] = l
    deg = np.array([G.degree(n) for n in nodes]); bins = np.digitize(deg, np.unique(np.quantile(deg, np.linspace(0, 1, 11)))[1:-1])
    by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    T_net = {d: {t for t in T_all.get(d, set()) if t in idx} for d in drugs}
    comm = rd(run, "15_network_algorithms", f"{acc}_Communities.csv")
    cmap = dict(zip(comm.Gene, comm.Community)) if comm is not None else {}
    mmap = dict(zip(comm.Gene, comm.Module)) if comm is not None else {}

    def dist_stats(A, B):
        a = np.array([idx[t] for t in A]); b = np.array([idx[t] for t in B]); M = D[np.ix_(a, b)]
        fin = np.isfinite(M)
        if not fin.any(): return np.nan, np.nan, 1.0
        near = np.mean(np.concatenate([np.where(np.isfinite(M).any(1), np.nanmin(np.where(fin, M, np.nan), 1), np.nan),
                                       np.where(np.isfinite(M).any(0), np.nanmin(np.where(fin, M, np.nan), 0), np.nan)]).astype(float))
        mean = np.nanmean(np.where(fin, M, np.nan)); unreach = 1 - fin.mean()
        return near, mean, unreach
    def self_d(A):
        if len(A) < 2: return np.nan
        a = np.array([idx[t] for t in A]); M = D[np.ix_(a, a)].astype(float); np.fill_diagonal(M, np.inf)
        m = np.where(np.isfinite(M), M, np.nan); return np.nanmean(np.nanmin(m, 1)) if np.isfinite(M).any() else np.nan
    def random_like(A):
        return [nodes[rng.choice(by_bin[bins[idx[t]]])] for t in A]

    rows = []
    for a, b in combinations(drugs, 2):
        A, B = T_net[a], T_net[b]; r = dict(Drug_A=a, Drug_B=b, Group_A=grp[a], Group_B=grp[b])
        ga, gb = grp[a][0], grp[b][0]
        r["Evidence_label"] = ("both well supported" if ga == gb == "1" else
                               "supported; antiproliferative caution" if {ga, gb} <= {"1", "2"} else "exploratory (weakly supported compound)")
        r["Targets_A"], r["Targets_B"] = len(T_all.get(a, set())), len(T_all.get(b, set()))
        r["Target_Jaccard_all"] = jacc(T_all.get(a, set()), T_all.get(b, set()))
        r["Shared_targets"] = ";".join(sorted(T_all.get(a, set()) & T_all.get(b, set())))
        r["Network_targets_A"], r["Network_targets_B"] = len(A), len(B)
        r["Target_Jaccard_network"] = jacc(A, B)
        if A and B:
            near, mean, unreach = dist_stats(A, B); r.update(Nearest_target_distance=near, Mean_pairwise_distance=mean, Unreachable_fraction=round(unreach, 3))
            dA, dB = self_d(A), self_d(B)
            r["Separation"] = mean - (dA + dB) / 2 if np.isfinite([dA, dB]).all() else np.nan
            null = np.array([dist_stats(random_like(A), random_like(B))[0] for _ in range(N_NULL)]); null = null[np.isfinite(null)]
            if len(null) > 10 and np.isfinite(near):
                r["Nearest_null_mean"] = null.mean(); r["Nearest_Z"] = (near - null.mean()) / (null.std() or np.nan)
                r["P_more_separated"] = (np.sum(null >= near) + 1) / (len(null) + 1)
            ca, cb = {cmap[t] for t in A if t in cmap}, {cmap[t] for t in B if t in cmap}
            r["Community_Jaccard"] = jacc(ca, cb); r["Modules_A"] = ";".join(sorted({mmap[t] for t in A if t in mmap})); r["Modules_B"] = ";".join(sorted({mmap[t] for t in B if t in mmap}))
        rows.append(r)
    P = pd.DataFrame(rows)
    # chemistry
    F = fps(); have = [d for d in drugs if d in F]
    tan = pd.DataFrame(np.nan, index=drugs, columns=drugs)
    for a in have:
        for b in have: tan.loc[a, b] = DataStructs.TanimotoSimilarity(F[a], F[b])
    P["Tanimoto"] = [tan.loc[a, b] for a, b in zip(P.Drug_A, P.Drug_B)]
    P["Structure_missing"] = [", ".join(d for d in (a, b) if d not in F) for a, b in zip(P.Drug_A, P.Drug_B)]
    sameS = dict(zip(lib.Drug, lib.Same_substance_as))
    P["Same_substance"] = [b in str(sameS.get(a, "")).split("; ") for a, b in zip(P.Drug_A, P.Drug_B)]
    # pathways
    enr = rd(run, "17_pathway_annotation", f"{acc}_DrugPathways_Enrichr_long.csv")
    paths = {}
    if enr is not None:
        e = enr[enr.Library.astype(str).str.contains("|".join(ENRICHR_LIBS)) & (enr["Adjusted.P.value"] < 0.05)]
        paths = e.groupby("Drug").Term.apply(set).to_dict()
    P["Pathway_Jaccard"] = [jacc(paths.get(a, set()), paths.get(b, set())) for a, b in zip(P.Drug_A, P.Drug_B)]
    P["Pathways_A"] = [len(paths.get(a, set())) for a in P.Drug_A]; P["Pathways_B"] = [len(paths.get(b, set())) for b in P.Drug_B]
    # candidate complementary pairs: both supported, no shared target, not the same substance, targets further apart than random, little pathway overlap
    P["Complementary_candidate"] = ((P.Evidence_label != "exploratory (weakly supported compound)") & (P.Target_Jaccard_all.fillna(1) == 0) & (~P.Same_substance)
                                    & (P.Nearest_Z.fillna(-9) > 0) & (P.Pathway_Jaccard.fillna(0) < 0.2))
    P = P.sort_values(["Complementary_candidate", "Evidence_label", "Nearest_Z"], ascending=[False, True, False])
    P.to_csv(o / f"{k}_drug_pairs.csv", index=False)
    # chemical clusters
    clus = {}
    if len(have) >= 5:
        dm = 1 - tan.loc[have, have].values; np.fill_diagonal(dm, 0); Z = linkage(squareform(dm, checks=False), "average")
        sil = {kk: silhouette_score(dm, fcluster(Z, kk, "maxclust"), metric="precomputed") for kk in range(2, min(5, len(have)))}
        kbest = max(sil, key=sil.get); lab = fcluster(Z, kbest, "maxclust")
        clus = dict(zip(have, lab)) if sil[kbest] >= 0.25 else {}
        ari = {}
        for rad, bits in [(3, 2048), (2, 1024)]:
            F2 = fps(rad, bits); t2 = np.array([[DataStructs.TanimotoSimilarity(F2[a], F2[b]) for b in have] for a in have])
            d2 = 1 - t2; np.fill_diagonal(d2, 0); ari[f"radius {rad}, {bits} bits"] = adjusted_rand_score(lab, fcluster(linkage(squareform(d2, checks=False), "average"), kbest, "maxclust"))
        chem_note = (f"best k={kbest}, silhouette {sil[kbest]:.2f}" + ("" if clus else " (weak: no clusters reported)") +
                     "; stability (ARI) " + ", ".join(f"{a}: {v:.2f}" for a, v in ari.items()))
    else: chem_note = "too few structures for clustering"
    # library view for this context
    lv = c[["Drug", "Group", "Evidence_rank", "Targets"]].merge(lib[["Drug", "CID", "SMILES_standardised", "InChIKey_standardised", "Structure_status", "Same_substance_as"]], on="Drug", how="left")
    lv["All_known_targets"] = [len(T_all.get(d, set())) for d in lv.Drug]; lv["Network_targets"] = [len(T_net[d]) for d in lv.Drug]
    lv["Significant_pathways"] = [len(paths.get(d, set())) for d in lv.Drug]; lv["Chemical_cluster"] = [clus.get(d, "") for d in lv.Drug]
    lv.to_csv(o / f"{k}_standardized_drug_library.csv", index=False)
    # sensitivity
    sens = {"Pairs": len(P), "Pairs with network distances": int(P.Nearest_target_distance.notna().sum()),
            "Spearman nearest vs mean distance": round(P[["Nearest_target_distance", "Mean_pairwise_distance"]].corr("spearman").iloc[0, 1], 2),
            "Spearman nearest vs separation": round(P[["Nearest_target_distance", "Separation"]].corr("spearman").iloc[0, 1], 2),
            "Candidate pairs": int(P.Complementary_candidate.sum()),
            "Candidate pairs, both compounds well supported": int((P.Complementary_candidate & (P.Evidence_label == "both well supported")).sum()),
            "Candidate pairs with >= 2 network targets each": int((P.Complementary_candidate & (P.Network_targets_A >= 2) & (P.Network_targets_B >= 2)).sum()),
            "Candidate pairs also more separated at p < 0.05": int((P.Complementary_candidate & (P.P_more_separated < 0.05)).sum()),
            "Chemical clustering": chem_note,
            "Not tested": "alternative STRING confidence thresholds (would require re-querying STRING)"}
    pd.Series(sens).to_csv(o / f"{k}_sensitivity.csv", header=["Value"])
    print("\n".join(f"  {a}: {v}" for a, v in sens.items()))
    # figures
    order = drugs if not have else (have[i] for i in leaves_list(linkage(squareform(np.clip(1 - tan.loc[have, have].fillna(0).values, 0, None) * (1 - np.eye(len(have))), checks=False), "average")))
    order = list(order) + [d for d in drugs if d not in have]
    nice = lambda d: d if re.search(r"\d", d) else d.title()
    def mat(col):
        M = pd.DataFrame(np.nan, index=order, columns=order)
        for r in P.itertuples(): M.loc[r.Drug_A, r.Drug_B] = M.loc[r.Drug_B, r.Drug_A] = getattr(r, col)
        M.index = M.columns = [nice(d) for d in order]; return M
    for col, title, cm in [("Target_Jaccard_all", "Target overlap (Jaccard)", "Reds"), ("Tanimoto", "Chemical similarity (Tanimoto, Morgan r=2)", "Blues"),
                           ("Nearest_Z", "Network separation of targets (Z vs random; higher = further apart)", "PuOr"), ("Pathway_Jaccard", "Pathway overlap (Jaccard)", "Greens")]:
        M = mat(col); fig, ax = plt.subplots(figsize=(7.2, 6.2))
        sns.heatmap(M, cmap=cm, center=0 if col == "Nearest_Z" else None, square=True, linewidths=0.3, linecolor="white", cbar_kws={"shrink": 0.6},
                    mask=M.isna(), ax=ax, xticklabels=True, yticklabels=True)
        ax.set_title(f"{label}: {title}", fontsize=9, loc="left"); ax.tick_params(labelsize=6); ax.set_facecolor("#eeeeee")
        save(fig, o / f"{k}_heatmap_{col}")
    cand = P[P.Complementary_candidate]
    if len(cand):
        H = nx.Graph(); gcol = {"1": "#B2182B", "2": "#F4A582", "3": "#BBBBBB"}
        for r in cand.itertuples(): H.add_edge(nice(r.Drug_A), nice(r.Drug_B), w=r.Nearest_Z)
        fig, ax = plt.subplots(figsize=(6.5, 5.5)); pos = nx.spring_layout(H, seed=SEED, k=0.9)
        nx.draw_networkx_edges(H, pos, width=[0.5 + H[u][v]["w"] / 2 for u, v in H.edges], alpha=0.5, ax=ax)
        nx.draw_networkx_nodes(H, pos, node_color=[gcol[grp[next(d for d in drugs if nice(d) == n)][0]] for n in H.nodes], node_size=260, ax=ax)
        nx.draw_networkx_labels(H, pos, font_size=6, ax=ax); ax.axis("off")
        ax.set_title(f"{label}: exploratory complementary pairs\n(no shared targets, targets further apart than random, low pathway overlap)", fontsize=8, loc="left")
        save(fig, o / f"{k}_candidate_pair_network")
    summary.append(dict(Context=k, Compounds=len(drugs), With_structure=len(have), Pairs=len(P), Candidate_pairs=int(P.Complementary_candidate.sum()),
                        Both_well_supported=sens["Candidate pairs, both compounds well supported"]))
    top = cand.head(10)
    if len(top):
        print("  top exploratory pairs:"); print(top[["Drug_A", "Drug_B", "Evidence_label", "Nearest_Z", "Tanimoto", "Pathway_Jaccard"]].to_string(index=False))

pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
print("\n>>> DONE. Outputs in", OUT)
