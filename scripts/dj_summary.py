"""Distal-junction (DJ) dosage per line = direct count of rDNA-bearing acrocentric p-arms.

DJ_copies = median per-k-mer hits over DJ k-mers with n_home == 5 (once per DJ in CHM13), GC-matched to the control
windows (11-14 GC of 31), divided by c_k_haploid (same normalization as summarize.py; per-line read length).
Assumes median control depth = 2 copies, like the rDNA raw values; DJ_copies * P_ctrl/2 = absolute DJ copies.
Compared with DepMap acro_CN (sum of proximal q-arm CN of chr13/14/15/21/22; normal = 10).
Usage: dj_summary.py [outname=dj_region]
"""
import sys, json, glob, numpy as np, pandas as pd
K = 31
out = sys.argv[1] if len(sys.argv) > 1 else "dj_region"
d = np.load("data/kmers/kmers_dj.npz"); keep = (d["n_out"] == 0) & (d["n_home"] > 0)
km, nh = d["kmer"][keep], d["n_home"][keep]
gc = np.zeros(len(km), int)
for j in range(K):
    b = (km >> np.uint64(2 * j)) & np.uint64(3); gc += ((b == 1) | (b == 2))
once = nh == 5; band = once & (gc >= 11) & (gc <= 14)
regs = pd.read_csv("data/kmers/hg38_regions_v2.tsv", sep="\t")
cov_ = pd.read_csv("data/meta/covariates.tsv", sep="\t").set_index("model_id")
def read_length(st, cov):
    """Sequenced read length. Prefer control-window query bases (recorded by quantify.py >= 2026-11-19.3).
    Older runs: bases/reads over fetched reads is biased low by short unplaced-unmapped reads, so snap it to the
    two CCLE read lengths (BAMs 101 bp, CRAMs 151 bp), verified directly on control-window reads of 8 lines."""
    if "query_bases" in cov and cov.reads.sum() > 0: return cov.query_bases.sum() / cov.reads.sum()
    return 101.0 if st["bases"] / st["reads"] < 120 else 151.0


rows = []
for sj in sorted(glob.glob(f"results/*/{out}/stats.json")):
    st = json.load(open(sj)); p = sj.rsplit("/", 1)[0]; m = p.split("/")[1]
    c = pd.read_csv(f"{p}/coverage.tsv", sep="\t")
    c = c[~c.apply(lambda r: bool(((regs.contig == r.contig) & (regs.start < r.end) & (regs.end > r.start)).any()), axis=1)]
    depth = (c.aligned_bases / (c.end - c.start)).median(); R = read_length(st, c); ck = depth / 2 * (R - K + 1) / R
    kc = pd.read_csv(f"{p}/kmer_counts.tsv", sep="\t"); assert len(kc) == len(km)
    h = kc.hits.values
    P = cov_.P_ctrl.get(m, np.nan)
    rows.append(dict(model_id=m, line=st["line"], sex=st["sex"], read_len=round(R, 1), depth=round(depth, 2),
                     DJ_raw=np.median(h[once]) / ck, DJ_gc=np.median(h[band]) / ck, P_ctrl=P,
                     DJ_abs=np.median(h[band]) / ck * P / 2, acro_CN=cov_.acro_CN.get(m, np.nan)))
df = pd.DataFrame(rows); df.to_csv(f"results/dj_summary_{out}.tsv", sep="\t", index=False)
print(f"DJ k-mers once-per-DJ: {once.sum():,}; GC-matched: {band.sum():,}")
print(df.round(2).to_string(index=False))
