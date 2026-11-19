"""Compare full-stream vs region-mode hits for a line; show where region mode misses hits."""
import sys, pandas as pd
m = sys.argv[1]; RD = sys.argv[2] if len(sys.argv) > 2 else "region"; REGFILE = sys.argv[3] if len(sys.argv) > 3 else "data/kmers/hg38_regions_v1.tsv"
F = pd.read_csv(f"results/{m}/full/hits.tsv", sep="\t"); Rg = pd.read_csv(f"results/{m}/{RD}/hits.tsv", sep="\t")
cf = pd.read_csv(f"results/{m}/full/coverage.tsv", sep="\t"); cr = pd.read_csv(f"results/{m}/{RD}/coverage.tsv", sep="\t")
print("control coverage identical:", cf.aligned_bases.tolist() == cr.aligned_bases.tolist())
t = pd.DataFrame({"full": F.groupby("group").hits.sum(), "region": Rg.groupby("group").hits.sum()}).fillna(0)
t["region/full"] = (t.region / t.full).round(4); print(t.to_string())
reg = pd.read_csv(REGFILE, sep="\t")
def in_reg(r):
    if r.contig == "*": return True
    s = reg[reg.contig == r.contig]
    return ((s.start < (r.bin_Mb + 1) * 1e6) & (s.end > r.bin_Mb * 1e6)).any()  # approximate: any overlap credits the whole 1 Mb bin
F["in_regions"] = F.apply(in_reg, axis=1)
miss = F[~F.in_regions]
print("\nhits in 1Mb bins outside region list, top per group:")
for g, s in miss.groupby("group"):
    s = s.groupby(["contig", "bin_Mb", "unmapped"]).hits.sum().sort_values(ascending=False)
    print(f"{g}: {s.sum():,} ({s.sum()/t.loc[g,'full']:.2%} of full)"); print(s.head(6).to_string())
