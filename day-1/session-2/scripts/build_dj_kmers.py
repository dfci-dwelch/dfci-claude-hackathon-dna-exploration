"""Distal-junction (DJ) k-mers: a direct measure of rDNA-bearing acrocentric p-arm dosage.

In CHM13 every rDNA array is flanked distally by the same ~300 kb block (non-satellite 'ct' ~120 kb + bsat 12 kb +
'ct' ~170 kb). The two 'ct' segments next to each of the 5 arrays (10 intervals) are the DJ home. Keep 31-mers whose
every CHM13 occurrence lies in a DJ interval (shared between acrocentrics is fine: one group). A k-mer seen once per
DJ (n_home == 5) counts p-arms: median hits / c_k_haploid = DJ copies per genome (normal diploid = 10).
Output: data/kmers/kmers_dj.npz (same format as kmers_v1.npz) via build_kmers.main.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import build_kmers as bk

bed = bk.bed
home = []
for c in ["chr13", "chr14", "chr15", "chr21", "chr22"]:
    r_start = int(bed[(bed.chrom == c) & bed.name.str.startswith("rDNA")].start.iloc[0])
    flank = bed[(bed.chrom == c) & (bed.end <= r_start) & (bed.end > r_start - 400_000)].sort_values("start")
    ct = flank[flank.name.str.startswith("ct_")].tail(2)
    assert len(ct) == 2, (c, flank)
    home += ct[["chrom", "start", "end"]].values.tolist()
print("DJ home intervals:", [(c, s, e, e - s) for c, s, e in home])
bk.GROUPS = [("DJ", ("bed", home), home)]

if __name__ == "__main__":
    bk.main(tag="dj")
