"""Convert hit tallies to array size relative to CHM13v2 (haploid; chrY from HG002).

k-mer coverage per haploid copy:  c_k = D/2 * (R-K+1)/R
  D = median aligned depth over control windows (windows overlapping the v2 target regions excluded)
  R = this line's read length (control-window reads; older runs: 101 bp BAMs / 151 bp CRAMs, see read_length)
rel_CHM13         = hits / (sum of CHM13 home multiplicity * c_k)   (all groups; ~1 = one CHM13-sized haploid copy)
units_per_genome  = hits / (in-unit occurrences of kept k-mers * c_k)         (rDNA; mean-based)
units_median      = median per-k-mer hits / c_k                                 (rDNA; robust to polymorphic k-mers)
units_gc          = median per-k-mer hits / c_k over k-mers that occur once per unit and whose GC
                    (11-14 of 31 bases) matches the control windows (median 38.7%)  -> default rDNA metric;
                    removes GC-dependent coverage loss, which is batch-dependent (worse in the 101 bp BAMs)
units_best        = units_gc for 45S; units_median for 5S (5S k-mer counts show no GC dependence in either batch,
                    and only 73 5S k-mers fall in the GC band)
Y arrays have no per-k-mer counts, so they get no GC correction; analyze.py adds sequencing batch as a covariate.
Optional argv[1]: comma-separated modes to keep, e.g. full,region_v2.
"""
import sys, os, json, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kmerlib import canonical_kmers, encode

K = 31
GC_BAND = (11, 14)  # GC bases per 31-mer, matched to control windows
d = np.load("data/kmers/kmers_v1.npz"); names = list(d["group_names"])
keep = (d["n_out"] == 0) & (d["n_home"] > 0)
KM, GR = d["kmer"][keep], d["group"][keep]
home = pd.Series(d["n_home"][keep]).groupby(GR).sum(); home.index = [names[i] for i in home.index]

# per-k-mer properties for rDNA groups, in the same order as kmer_counts.tsv (ascending kept index)
UNITS = {"rDNA_45S": "data/rdna/KY962518.fa", "rDNA_5S": "data/rdna/X12811.fa"}
pk = np.nonzero([names[g].startswith("rDNA") for g in GR])[0]
gc = np.zeros(len(pk), int)
for j in range(K):
    b = (KM[pk] >> np.uint64(2 * j)) & np.uint64(3); gc += ((b == 1) | (b == 2))
mult = np.zeros(len(pk), int); unit_occ = {}
for g, f in UNITS.items():
    s = "".join(l.strip() for l in open(f) if not l.startswith(">"))
    km, _ = canonical_kmers(encode(s)); km = np.sort(km)
    sel = GR[pk] == names.index(g)
    mult[sel] = np.searchsorted(km, KM[pk][sel], side="right") - np.searchsorted(km, KM[pk][sel], side="left")
    unit_occ[g] = int(mult[sel].sum())
KPROP = pd.DataFrame({"gc": gc, "mult": mult})
band = (KPROP.gc >= GC_BAND[0]) & (KPROP.gc <= GC_BAND[1]) & (KPROP.mult == 1)

regs = pd.read_csv("data/kmers/hg38_regions_v2.tsv", sep="\t")
def overlaps_target(r):
    s = regs[regs.contig == r.contig]; return bool(((s.start < r.end) & (s.end > r.start)).any())

def read_length(st, cov):
    """Sequenced read length. Prefer control-window query bases (recorded by quantify.py >= 2026-11-19.3).
    Older runs: bases/reads over fetched reads is biased low by short unplaced-unmapped reads, so snap it to the
    two CCLE read lengths (BAMs 101 bp, CRAMs 151 bp), verified directly on control-window reads of 8 lines."""
    if "query_bases" in cov and cov.reads.sum() > 0: return cov.query_bases.sum() / cov.reads.sum()
    return 101.0 if st["bases"] / st["reads"] < 120 else 151.0


rows = []
for sj in sorted(glob.glob("results/*/*/stats.json")):
    st = json.load(open(sj)); p = sj.rsplit("/", 1)[0]; mode = p.rsplit("/", 1)[1]
    tag = (st.get("input_fingerprint") or [None, "v1"])[1]   # runs before 2026-11-19 have no fingerprint: all v1
    if tag != "v1": continue  # other k-mer sets (e.g. 'dj') are summarized by their own scripts
    h = pd.read_csv(f"{p}/hits.tsv", sep="\t"); c = pd.read_csv(f"{p}/coverage.tsv", sep="\t")
    c = c[~c.apply(overlaps_target, axis=1)]
    depth = (c.aligned_bases / (c.end - c.start)).median()
    R = read_length(st, c); ck = depth / 2 * (R - K + 1) / R
    hits = h.groupby("group").hits.sum().reindex(home.index, fill_value=0)
    unm = h[h.unmapped].groupby("group").hits.sum().reindex(home.index, fill_value=0)
    med = medgc = pd.Series(dtype=float)
    if os.path.exists(f"{p}/kmer_counts.tsv"):
        kc = pd.read_csv(f"{p}/kmer_counts.tsv", sep="\t"); assert len(kc) == len(KPROP)
        med = kc.groupby("group").hits.median(); medgc = kc[band.values].groupby("group").hits.median()
    for g in home.index:
        rows.append(dict(line=st["line"], model_id=p.split("/")[1], sex=st["sex"], mode=mode, read_len=round(R, 1),
                         depth=round(depth, 2), n_ctrl=len(c), group=g, hits=int(hits[g]),
                         frac_unmapped=round(unm[g] / max(hits[g], 1), 4), rel_CHM13=hits[g] / (home[g] * ck),
                         units_per_genome=hits[g] / (unit_occ[g] * ck) if g in unit_occ else np.nan,
                         units_median=med.get(g, np.nan) / ck, units_gc=medgc.get(g, np.nan) / ck,
                         units_best=(medgc if g == "rDNA_45S" else med).get(g, np.nan) / ck))
df = pd.DataFrame(rows)
if len(sys.argv) > 1: df = df[df["mode"].isin(sys.argv[1].split(","))]
df.to_csv("results/summary.tsv", sep="\t", index=False)
print(f"GC-matched rDNA k-mers: 45S {int((band & (GR[pk] == names.index('rDNA_45S'))).sum())}, "
      f"5S {int((band & (GR[pk] == names.index('rDNA_5S'))).sum())}; read lengths: {sorted(df.read_len.round().unique())}")
print(df.pivot_table(index=["line", "sex", "mode", "read_len", "depth"], columns="group", values="rel_CHM13").round(4).to_string())
r = df[df.group.str.startswith("rDNA")].pivot_table(index=["line", "mode"], columns="group", values=["units_median", "units_gc", "units_best"]).round(0)
print("\nrDNA units per genome (diploid-equivalent, before ploidy correction):\n" + r.to_string())
