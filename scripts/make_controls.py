"""Pick 200 x 50kb single-copy-ish control windows on hg38 chr1-22: no N within 3 Mb, outside v1 target regions.

Windows are chosen against region list v1 (kept so the file used for all 89 lines is reproduced exactly).
4 of the 200 overlap v2 regions; summarize.py drops windows that overlap v2 before taking the median depth.
"""
import numpy as np, pandas as pd, pysam
fa = pysam.FastaFile("data/ref/Homo_sapiens_assembly38.fasta"); rng = np.random.default_rng(42)
tgt = pd.read_csv("data/kmers/hg38_regions_v1.tsv", sep="\t")
W, PAD, rows = 50_000, 3_000_000, []
chroms = [f"chr{i}" for i in range(1, 23)]; L = np.array([fa.get_reference_length(c) for c in chroms], float)
while len(rows) < 200:
    c = rng.choice(chroms, p=L / L.sum()); n = fa.get_reference_length(c)
    s = int(rng.integers(PAD, n - PAD - W)) // 1000 * 1000
    ctx = fa.fetch(c, s - PAD, s + W + PAD).upper()
    if "N" in ctx: continue
    if ((tgt.contig == c) & (tgt.start < s + W + PAD) & (tgt.end > s - PAD)).any(): continue
    if any(r[0] == c and abs(r[1] - s) < W for r in rows): continue
    rows.append((c, s, s + W))
pd.DataFrame(rows, columns=["contig", "start", "end"]).sort_values(["contig", "start"]).to_csv("data/kmers/control_windows.tsv", sep="\t", index=False)
print(pd.DataFrame(rows)[0].value_counts().sort_index().to_dict())
