"""Build array-specific 31-mer sets and test their uniqueness against all of CHM13 v2.0.

For each target group, collect its canonical k-mers (from an external unit sequence or from
CenSat intervals in CHM13), then scan every CHM13 base and count, per k-mer, occurrences
inside the group's "home" intervals vs. anywhere else. K-mers found anywhere outside home
(or shared between groups) are discarded.

Outputs (data/kmers/):
  kmers_<tag>.npz      kmer, group, n_home, n_out, unit_pos (offset in source unit; -1 for CenSat)
  out_hits_<tag>.tsv   where non-home hits land (group, chrom, 1Mb bin, count) - for review
  summary_<tag>.tsv
"""
import sys, re, collections
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, pysam
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from kmerlib import encode, canonical_kmers, lookup

CHM13 = "data/ref/chm13v2.0.fa.gz"
BED = "data/ref/chm13v2.0_censat_v2.1.bed"
OUT = "data/kmers"
CHUNK = 10_000_000

bed = pd.read_csv(BED, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "name"])


def fasta_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def intervals(pattern):
    return bed[bed.name.str.contains(pattern, regex=True)][["chrom", "start", "end"]].values.tolist()


# group: (name, source, home_intervals)
GROUPS = [
    ("rDNA_45S", ("seq", "data/rdna/KY962518.fa"), intervals(r"^rDNA_")),
    ("rDNA_5S", ("seq", "data/rdna/X12811.fa"), [["chr1", 227_740_000, 228_030_000]]),  # 5S array 1q42, not in CenSat; bounds from X12811 k-mer density
    ("Y_DYZ1_HSat3", ("bed", intervals(r"^hsat3_Y_.*DYZ1")), intervals(r"^hsat3_Y_.*DYZ1")),
    ("Y_DYZ2_HSat1B", ("bed", intervals(r"^hsat1B_Y_.*DYZ2")), intervals(r"^hsat1B_Y_.*DYZ2")),
    ("Y_DYZ3_aSat", ("bed", intervals(r"^hor_Y_")), intervals(r"^hor_Y_")),
]


def collect_targets():
    fa = pysam.FastaFile(CHM13)
    ks, gs, ps = [], [], []
    for gi, (name, (kind, src), _) in enumerate(GROUPS):
        if kind == "seq":
            km, pos = canonical_kmers(encode(fasta_seq(src)))
        else:
            parts = [canonical_kmers(encode(fa.fetch(c, s, e)))[0] for c, s, e in src]
            km = np.concatenate(parts); pos = np.full(len(km), -1, np.int64)
        km, first = np.unique(km, return_index=True)
        ks.append(km); gs.append(np.full(len(km), gi, np.int16)); ps.append(pos[first])
        print(f"  {name}: {len(km):,} distinct k-mers", flush=True)
    k = np.concatenate(ks); g = np.concatenate(gs); p = np.concatenate(ps)
    order = np.argsort(k, kind="stable"); k, g, p = k[order], g[order], p[order]
    # k-mers in more than one group -> drop
    dup = np.zeros(len(k), bool); same = k[1:] == k[:-1]; dup[1:] |= same; dup[:-1] |= same
    print(f"  shared between groups (dropped): {dup.sum():,}")
    return k[~dup], g[~dup], p[~dup]


def scan_chunk(args):
    chrom, s, e, targets, groups = args
    fa = pysam.FastaFile(CHM13)
    km, pos = canonical_kmers(encode(fa.fetch(chrom, s, e)))
    pos += s
    hit, idx = lookup(targets, km)
    hpos = pos[hit]; hgrp = groups[idx]
    home = np.zeros(len(idx), bool)
    for gi, (_, _, iv) in enumerate(GROUPS):
        for c, a, b in iv:
            if c == chrom:
                home |= (hgrp == gi) & (hpos >= a) & (hpos + 31 <= b)
    n_home = np.bincount(idx[home], minlength=len(targets))
    n_out = np.bincount(idx[~home], minlength=len(targets))
    bins = collections.Counter(zip(hgrp[~home].tolist(), (hpos[~home] // 1_000_000).tolist()))
    return chrom, n_home, n_out, bins


def main(tag="v1", workers=8):
    import os; os.makedirs(OUT, exist_ok=True)
    print("collecting target k-mers", flush=True)
    k, g, p = collect_targets()
    fa = pysam.FastaFile(CHM13)
    jobs = [(c, s, min(s + CHUNK + 30, L), k, g) for c, L in zip(fa.references, fa.lengths)
            for s in range(0, L, CHUNK)]
    print(f"scanning CHM13 in {len(jobs)} chunks", flush=True)
    n_home = np.zeros(len(k), np.int64); n_out = np.zeros(len(k), np.int64); rows = []
    with ProcessPoolExecutor(workers) as ex:
        for i, (chrom, h, o, bins) in enumerate(ex.map(scan_chunk, jobs, chunksize=1)):
            n_home += h; n_out += o
            rows += [(GROUPS[gi][0], chrom, b, n) for (gi, b), n in bins.items()]
            if i % 50 == 0: print(f"  {i}/{len(jobs)}", flush=True)
    np.savez(f"{OUT}/kmers_{tag}.npz", kmer=k, group=g, n_home=n_home, n_out=n_out, unit_pos=p,
             group_names=np.array([x[0] for x in GROUPS]))
    oh = pd.DataFrame(rows, columns=["group", "chrom", "bin_Mb", "n"]).groupby(["group", "chrom", "bin_Mb"]).n.sum()
    oh.reset_index().sort_values(["group", "n"], ascending=[True, False]).to_csv(f"{OUT}/out_hits_{tag}.tsv", sep="\t", index=False)
    s = pd.DataFrame({"group": [GROUPS[i][0] for i in g], "home": n_home, "out": n_out})
    summ = s.groupby("group").agg(distinct=("out", "size"), absent_from_chm13=("out", lambda x: ((x == 0) & (s.loc[x.index, "home"] == 0)).sum()),
                                  unique_to_home=("out", lambda x: ((x == 0) & (s.loc[x.index, "home"] > 0)).sum()),
                                  found_elsewhere=("out", lambda x: (x > 0).sum()))
    summ.to_csv(f"{OUT}/summary_{tag}.tsv", sep="\t"); print(summ.to_string())


if __name__ == "__main__":
    main()
