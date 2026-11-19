"""Region list v2 = v1 + 1 Mb bins (or whole small contigs) that held missed hits in full streams.

Y groups: bins from the male full stream (HCC2935, ACH-000150) only, so female background doesn't drive
the choice. rDNA: bins from both full streams. Add the smallest set of missed bins that brings
region-mode capture to >=99% (TARGET) of full-stream hits per group (counting unplaced unmapped '*', always fetched).

Caveat: full-stream hits are only binned at 1 Mb, so a bin counts as "covered" if any v1 region overlaps it.
That makes the printed capture an optimistic estimate (a partly covered bin is credited in full; a strict
"fully contained" bound is uninformative because v1 regions are mostly < 1 Mb). The exact, read-level capture comes
from compare_modes.py (region_v2 hits / full hits): HCC2935 DYZ2 98.95%, A427 DYZ2 99.85%, all others >= 99.4%.
Selection logic is unchanged so that data/kmers/hg38_regions_v2.tsv (used for all 89 lines) is reproduced exactly.
"""
import pandas as pd
v1 = pd.read_csv("data/kmers/hg38_regions_v1.tsv", sep="\t")
fai = pd.read_csv("data/ref/Homo_sapiens_assembly38.fasta.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "L"]).set_index("c").L
src = {"ACH-000150": None, "ACH-000861": ["rDNA_45S", "rDNA_5S"]}
H = []
for m, groups in src.items():
    h = pd.read_csv(f"results/{m}/full/hits.tsv", sep="\t")
    if groups: h = h[h.group.isin(groups)]
    H.append(h)
H = pd.concat(H).groupby(["group", "contig", "bin_Mb"]).hits.sum().reset_index()
def covered(r):
    s = v1[v1.contig == r.contig]; return ((s.start < (r.bin_Mb + 1) * 1e6) & (s.end > r.bin_Mb * 1e6)).any()
H["cov"] = H.apply(lambda r: r.contig == "*" or covered(r), axis=1)  # '*' (unplaced unmapped) is always fetched

TARGET = 0.99
add = []
for g, s in H.groupby("group"):
    tot = s.hits.sum(); have = s[s["cov"]].hits.sum(); miss = s[~s["cov"]].sort_values("hits", ascending=False)
    need = TARGET * tot - have; cum = miss.hits.cumsum()
    pick = miss[cum.shift(fill_value=0) < need]
    print(f"{g}: v1 capture {have/tot:.4f}; adding {len(pick)} bins -> {(have+pick.hits.sum())/tot:.4f} "
          f"(optimistic bin-level estimate; exact capture: compare_modes.py)")
    for r in pick.itertuples():
        add.append((r.contig, int(r.bin_Mb * 1e6), int(min((r.bin_Mb + 1) * 1e6, fai[r.contig])), g, int(r.hits)))
a = pd.DataFrame(add, columns=v1.columns)
v2 = pd.concat([v1, a]); v2.to_csv("data/kmers/hg38_regions_v2.tsv", sep="\t", index=False)
print(f"v1 {len(v1)} regions {((v1.end-v1.start).sum()/1e6):.1f} Mb -> v2 {len(v2)} regions {((v2.end-v2.start).sum()/1e6):.1f} Mb")
