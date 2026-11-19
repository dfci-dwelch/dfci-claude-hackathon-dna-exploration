"""Investigate the S100A8 gene effect ~ 5S rDNA association (review item).
1. Do S100/1q21 guides resemble the 5S or 45S unit (<=3 mismatches, both strands)?
2. Does it hold per screen library (Avana / KY / Humagne)?
3. Is it explained by copy number at the S100 locus (1q21.3) or 1q42, or by raw vs CN-corrected 5S?
"""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
unit = lambda f: "".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
U = {"5S": unit("data/rdna/X12811.fa"), "45S": unit("data/rdna/KY962518.fa")}

def best_mm(g, seq):
    L = len(g); best = L
    for s in (seq, rc(seq)):
        a = np.frombuffer(s.encode(), np.uint8); q = np.frombuffer(g.encode(), np.uint8)
        w = np.lib.stride_tricks.sliding_window_view(a, L)
        best = min(best, int((w != q).sum(1).min()))
    return best

genes = ["S100A8", "S100A12", "S100A4", "S100A9", "CLIP4"]
print("1. guide similarity to rDNA units (min mismatches over the unit, both strands)")
for lib in ["Avana", "KY", "Humagne"]:
    gm = pd.read_csv(f"data/depmap/{lib}GuideMap.csv")
    gm = gm[gm.Gene.fillna("").str.split(" ").str[0].isin(genes)]
    for r in gm.itertuples():
        g = r.sgRNA.upper()
        print(f"  {lib:7s} {r.Gene.split(' ')[0]:8s} {g}  nAlign={r.nAlignments}  used={r.UsedByChronos}  "
              f"mm5S={best_mm(g, U['5S'])}  mm45S={best_mm(g, U['45S'])}  {r.GenomeAlignment}")

d = pd.read_csv("results/analysis/per_line.tsv", sep="\t", index_col=0)
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0, usecols=lambda c: c == "" or c.startswith("Unnamed") or c.split(" ")[0] in genes)
ge.columns = [c.split(" ")[0] for c in ge.columns]
x = d.join(ge, how="inner")
sm = pd.read_csv("data/depmap/CRISPRScreenMap.csv"); sm["lib"] = sm.ScreenID.str.extract(r"\.([A-Z]+)\d*$")[0]
libs = sm.groupby("ModelID").lib.apply(lambda v: "+".join(sorted(set(v))))
x["libs"] = libs.reindex(x.index)
print("\n2. per screen library: spearman(S100A8 effect, r5_1q)")
print("  library mix:", x.libs.value_counts().to_dict())
for lib in ["AV", "KY", "HU"]:
    s = x[x.libs.fillna("").str.contains(lib)].dropna(subset=["S100A8", "r5_1q"])
    if len(s) > 8: print(f"  lines screened with {lib}: n={len(s)} rho={stats.spearmanr(s.r5_1q, s.S100A8)[0]:+.3f} p={stats.spearmanr(s.r5_1q, s.S100A8)[1]:.2g}")

seg = pd.read_csv("data/depmap/OmicsAbsoluteCNSegmentsProfile.csv", dtype={"Chromosome": str})
prof = pd.read_csv("data/depmap/OmicsProfiles.csv"); prof = prof[prof.ProfileID.isin(seg.ProfileID)].copy()
prof["rank"] = prof.Datatype.map({"wgs": 0, "wes": 1}); pick = prof.sort_values("rank").drop_duplicates("ModelID").set_index("ModelID").ProfileID
def cn_at(m, chrom, a, b):
    if m not in pick.index: return np.nan
    s = seg[(seg.ProfileID == pick[m]) & (seg.Chromosome == chrom) & (seg.End > a) & (seg.Start < b)]
    if s.empty: return np.nan
    L = np.minimum(s.End, b) - np.maximum(s.Start, a); return float((s.SegmentAbsoluteCN * L).sum() / L.sum())
x["cn_1q21"] = [cn_at(m, "1", 153_300_000, 153_600_000) for m in x.index]   # S100 cluster, hg38
x["rel_1q21"] = x.cn_1q21 / x.P_ctrl
x["rel_1q42"] = x.cn_1q42 / x.P_ctrl
z = x.dropna(subset=["S100A8", "r5_1q", "cn_1q21"])
print(f"\n3. copy-number confounding (n={len(z)})")
for a, b in [("r5_1q", "S100A8"), ("r5_raw", "S100A8"), ("rel_1q21", "S100A8"), ("rel_1q42", "S100A8"), ("rel_1q21", "r5_1q"), ("rel_1q42", "r5_1q"), ("rel_1q42", "r5_raw")]:
    r, p = stats.spearmanr(z[a], z[b]); print(f"  spearman({a}, {b}) = {r:+.3f}  p={p:.2g}")
for f in ["S100A8 ~ np.log2(r5_1q)", "S100A8 ~ np.log2(r5_1q) + rel_1q21 + rel_1q42", "S100A8 ~ np.log2(r5_1q) + rel_1q21 + rel_1q42 + C(lineage) + C(batch)",
          "S100A8 ~ np.log2(r5_raw) + C(lineage) + C(batch)"]:
    m = smf.ols(f, data=z).fit(); k = [i for i in m.params.index if "r5" in i][0]
    print(f"  {f:70s} beta={m.params[k]:+.4f} p={m.pvalues[k]:.2g}")
print("\n  S100 cluster gene effects correlate with each other?")
print(z[["S100A8", "S100A12", "S100A4", "S100A9"]].corr(method="spearman").round(2).to_string())
print("\n  S100A8 effect by 1q21 relative CN (gain = rel > 1.25):")
print(z.assign(gain=np.where(z.rel_1q21 > 1.25, "gain", "no_gain")).groupby("gain")[["S100A8", "r5_1q"]].median().round(3).to_string())
