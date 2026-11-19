"""Region list v1: per group, the 100 kb hg38 bins holding 99.9% of k-mer hits (scan_hg38.py output),
adjacent bins (gap <= 1 bin) merged, padded by one bin each side. Writes data/kmers/hg38_regions_<tag>.tsv
(argv[1], default v1)."""
import sys, pandas as pd
TAG = sys.argv[1] if len(sys.argv) > 1 else "v1"   # e.g. 'dj' -> hg38_hits_dj.tsv -> hg38_regions_dj.tsv
d = pd.read_csv(f"data/kmers/hg38_hits_{TAG}.tsv", sep="\t"); B = 100_000; rows = []
for g, s in d.groupby("group"):
    s = s.sort_values("n", ascending=False); s["cf"] = s.n.cumsum() / s.n.sum()
    top = s[s.cf.shift(fill_value=0) < 0.999]
    for c, t in top.groupby("contig"):
        bins = sorted(t["bin"]); start = prev = bins[0]
        for b in bins[1:] + [None]:
            if b is not None and b <= prev + 2: prev = b; continue
            rows.append((c, max(0, (start - 1) * B), (prev + 2) * B, g, int(t[(t["bin"] >= start) & (t["bin"] <= prev)].n.sum())))
            if b is not None: start = prev = b
r = pd.DataFrame(rows, columns=["contig", "start", "end", "group", "hits"])
r.to_csv(f"data/kmers/hg38_regions_{TAG}.tsv", sep="\t", index=False)
print(len(r), "regions,", round((r.end - r.start).sum() / 1e6, 1), "Mb")
