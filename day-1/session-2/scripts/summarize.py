"""Convert hit tallies to array size relative to CHM13v2 (haploid; chrY from HG002).

rel = hits / (sum_home_multiplicity * c_k_haploid), c_k_haploid = median_control_depth/2 * (R-k+1)/R.
units_median = median per-k-mer hits / c_k_haploid (rDNA only).\nOptional argv[1]: comma-separated modes to show.\nrel ~= 1 means "one haploid copy the size of the CHM13 array"; diploid normal autosomal ~2; male Y ~1.
"""
import sys, json, glob, numpy as np, pandas as pd
R, K = 151, 31
d = np.load("data/kmers/kmers_v1.npz"); names = list(d["group_names"])
keep = (d["n_out"] == 0) & (d["n_home"] > 0)
home = pd.Series(d["n_home"][keep]).groupby(d["group"][keep]).sum(); home.index = [names[i] for i in home.index]
nk = pd.Series(keep[keep]).groupby(d["group"][keep]).size(); nk.index = home.index  # kept k-mers per group
rows = []
for sj in sorted(glob.glob("results/*/*/stats.json")):
    st = json.load(open(sj)); p = sj.rsplit("/", 1)[0]
    h = pd.read_csv(f"{p}/hits.tsv", sep="\t"); c = pd.read_csv(f"{p}/coverage.tsv", sep="\t")
    depth = (c.aligned_bases / (c.end - c.start)).median(); ck = depth / 2 * (R - K + 1) / R
    hits = h.groupby("group").hits.sum().reindex(home.index, fill_value=0)
    unm = h[h.unmapped].groupby("group").hits.sum().reindex(home.index, fill_value=0)
    import os
    med = pd.read_csv(f"{p}/kmer_counts.tsv", sep="\t").groupby("group").hits.median() if os.path.exists(f"{p}/kmer_counts.tsv") else pd.Series(dtype=float)
    for g in home.index:
        rows.append(dict(line=st["line"], sex=st["sex"], mode=p.rsplit("/", 1)[1], depth=round(depth, 1), group=g,
                         hits=int(hits[g]), frac_unmapped=round(unm[g] / max(hits[g], 1), 4), rel_CHM13=hits[g] / (home[g] * ck),
                         units_per_genome=hits[g] / (nk[g] * ck),  # meaningful for rDNA (single-unit k-mer source)
                         units_median=med.get(g, np.nan) / ck))  # robust to polymorphic k-mers; needs kmer_counts.tsv
df = pd.DataFrame(rows)
if len(sys.argv) > 1: df = df[df["mode"].isin(sys.argv[1].split(","))]  # e.g. full,region_v2
df.to_csv("results/summary.tsv", sep="\t", index=False)
print(df.pivot_table(index=["line", "sex", "mode", "depth"], columns="group", values="rel_CHM13").round(4).to_string())
r = df[df.group.str.startswith("rDNA")].pivot_table(index=["line", "mode"], columns="group", values=["units_per_genome", "units_median"]).round(0)
print("\nrDNA units per genome (diploid-equivalent):\n" + r.to_string())
