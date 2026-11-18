"""Locate the kept (CHM13-home-unique) k-mers in the Broad hg38 reference used by the CCLE CRAMs.

Output data/kmers/hg38_hits_<tag>.tsv: group, contig, 100kb bin, hit count. Used to choose the
hg38 regions to fetch reads from (reads carrying these k-mers likely map there, or stay unmapped).
"""
import sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, pysam
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from kmerlib import encode, canonical_kmers, lookup

HG38 = "data/ref/Homo_sapiens_assembly38.fasta"
CHUNK, BIN = 10_000_000, 100_000


def load_kept(tag):
    d = np.load(f"data/kmers/kmers_{tag}.npz")
    keep = (d["n_out"] == 0) & (d["n_home"] > 0)
    return d["kmer"][keep], d["group"][keep], list(d["group_names"])


def scan(args):
    contig, s, e, k, g = args
    km, pos = canonical_kmers(encode(pysam.FastaFile(HG38).fetch(contig, s, e)))
    hit, idx = lookup(k, km)
    df = pd.DataFrame({"g": g[idx], "bin": (pos[hit] + s) // BIN})
    return contig, df.value_counts().reset_index(name="n"), np.bincount(idx, minlength=len(k)) > 0


def main(tag="v1"):
    k, g, names = load_kept(tag)
    fa = pysam.FastaFile(HG38)
    jobs = [(c, s, min(s + CHUNK + 30, L), k, g) for c, L in zip(fa.references, fa.lengths) for s in range(0, L, CHUNK)]
    out, seen = [], np.zeros(len(k), bool)
    with ProcessPoolExecutor(8) as ex:
        for contig, df, found in ex.map(scan, jobs, chunksize=4):
            seen |= found
            if len(df): df.insert(0, "contig", contig); out.append(df)
    df = pd.concat(out); df["group"] = [names[i] for i in df.g]
    df = df.drop(columns="g").sort_values(["group", "n"], ascending=[True, False])
    df.to_csv(f"data/kmers/hg38_hits_{tag}.tsv", sep="\t", index=False)
    pres = pd.DataFrame({"group": [names[i] for i in g], "in_hg38": seen}).groupby("group").in_hg38.agg(["size", "mean"])
    print("fraction of kept k-mers present anywhere in hg38:\n", pres.to_string())
    for grp, sub in df.groupby("group"):
        tot = sub.n.sum(); c = sub.groupby("contig").n.sum().sort_values(ascending=False)
        print(f"\n{grp}: {tot:,} hits in {len(c)} contigs; top contigs (cum. fraction):")
        print((c.head(12).cumsum() / tot).round(3).to_string())


if __name__ == "__main__":
    main()
