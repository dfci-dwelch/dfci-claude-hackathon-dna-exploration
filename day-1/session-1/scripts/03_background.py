"""Step 5a: background ingredients for the ncORF anomaly test (SNVs only, all donors kept).

Per ORF: flank = ORF genomic span +/-50 kb, minus canonical CDS, minus all lifted ORF blocks.
Counts the 32 pyrimidine-centred trinucleotide contexts (from hg19.2bit; N-containing positions
dropped) in each ORF's territory and flank, extracts deduplicated PCAWG SNVs in the flank union,
and writes per-ORF x project observed counts for territory and flank.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import twobitreader

ROOT = Path(sys.argv[1]).resolve()
RAW, DER, RES = ROOT / "data/raw", ROOT / "data/derived", ROOT / "results"
MAF = RAW / "pcawg/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
FLANK = 50_000
CHROMS = [str(i) for i in range(1, 23)] + ["X"]  # MAF has no Y
KEY = ["Donor_ID", "Project_Code", "Chromosome", "Start_position", "Reference_Allele",
       "Tumor_Seq_Allele2"]

# ---- 32 contexts: base-5 code (A0 C1 G2 T3 N4) of prev/mid/next -> 0..31 or -1 ----
BASES = "ACGT"
CTX_NAMES = [p + m + n for m in "CT" for p in BASES for n in BASES]  # mid C: 0-15, T: 16-31
CTX_IDX = {c: i for i, c in enumerate(CTX_NAMES)}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def collapse(tri):
    if tri[1] in "GA":
        tri = "".join(COMP[b] for b in reversed(tri))
    return tri


LUT = np.full(125, -1, dtype=np.int8)
for p in range(5):
    for m in range(5):
        for n in range(5):
            if 4 in (p, m, n):
                continue
            LUT[25 * p + 5 * m + n] = CTX_IDX[collapse(BASES[p] + BASES[m] + BASES[n])]
ENC = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(BASES):
    ENC[ord(b)] = ENC[ord(b.lower())] = i

tb = twobitreader.TwoBitFile(str(RAW / "ucsc/hg19.2bit"))
chrom_len = {c: tb.sequence_sizes()["chr" + c] for c in CHROMS}


def ctx_array(chrom):
    a = ENC[np.frombuffer(tb["chr" + chrom][:].encode(), dtype=np.uint8)]
    ctx = np.full(len(a), -1, dtype=np.int8)
    ctx[1:-1] = LUT[25 * a[:-2] + 5 * a[1:-1] + a[2:]]
    return ctx


# ---- regions ----
blocks = pd.read_csv(DER / "orf_blocks_grch37.tsv", sep="\t", dtype={"chrom37": str})
terr = pd.read_csv(DER / "ncorf_territory_grch37.tsv", sep="\t", dtype={"Chromosome": str})
terr = terr[terr.Chromosome.isin(CHROMS)]
names = set(terr["name"])

cds = []
with gzip.open(RAW / "gencode_annot/gencode.v35lift37.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if f[2] == "CDS" and 'transcript_type "protein_coding"' in f[8]:
            cds.append((f[0].replace("chr", ""), int(f[3]) - 1, int(f[4])))
cds_df = pd.DataFrame(cds, columns=["Chromosome", "Start", "End"])
all_orf = blocks.rename(columns={"chrom37": "Chromosome", "start37": "Start", "end37": "End"})
mask = pr.PyRanges(pd.concat([cds_df, all_orf[["Chromosome", "Start", "End"]]])).merge()

span = (all_orf[all_orf["name"].isin(names)].groupby(["name", "Chromosome"])
        .agg(Start=("Start", "min"), End=("End", "max")).reset_index())
span["Start"] = (span.Start - FLANK).clip(lower=0)
span["End"] = np.minimum(span.End + FLANK, span.Chromosome.map(chrom_len))
flanks = pr.PyRanges(span).subtract(mask).df
flanks.to_csv(DER / "flanks_grch37.tsv", sep="\t", index=False)
flank_union = pr.PyRanges(flanks).merge()

# ---- context counts ----
terr_ctx, flank_ctx, union_ctx = {}, {}, np.zeros(32, dtype=np.int64)


def add(store, df, ctx):
    for nm, s, e in zip(df["name"], df.Start, df.End):
        v = ctx[s:e]
        store[nm] = store.get(nm, 0) + np.bincount(v[v >= 0], minlength=32)


ctx_cache = {}
for c in CHROMS:
    ctx = ctx_array(c)
    add(terr_ctx, terr[terr.Chromosome == c], ctx)
    add(flank_ctx, flanks[flanks.Chromosome == c], ctx)
    for s, e in zip(*(flank_union.df.query("Chromosome == @c")[k] for k in ("Start", "End"))):
        v = ctx[s:e]
        union_ctx += np.bincount(v[v >= 0], minlength=32)
    ctx_cache[c] = ctx  # kept for SNV context lookup below (~3 GB int8 total)
    print(f"  contexts counted: chr{c}", flush=True)

pd.DataFrame(terr_ctx, index=CTX_NAMES).T.rename_axis("name").to_csv(DER / "ctx_bases_territory.tsv", sep="\t")
pd.DataFrame(flank_ctx, index=CTX_NAMES).T.rename_axis("name").to_csv(DER / "ctx_bases_flank.tsv", sep="\t")
pd.Series(union_ctx, index=CTX_NAMES, name="bases").to_csv(DER / "ctx_bases_flank_union.tsv", sep="\t")


def snv_ctx(df):
    out = np.full(len(df), -1, dtype=np.int8)
    for c, idx in df.groupby("Chromosome").groups.items():
        out[df.index.get_indexer(idx)] = ctx_cache[c][df.loc[idx, "Start_position"].values - 1]
    return out


# ---- flank SNVs (stream MAF again) ----
COLS = KEY + ["Variant_Type", "ref_context"]
parts, n = [], 0
for chunk in pd.read_csv(MAF, sep="\t", usecols=COLS, chunksize=2_000_000,
                         dtype={"Chromosome": str}, low_memory=False):
    n += len(chunk)
    chunk = chunk[(chunk.Variant_Type == "SNP") & chunk.Chromosome.isin(CHROMS)]
    g = pr.PyRanges(chunk.assign(Start=chunk.Start_position - 1, End=chunk.Start_position))
    ov = g.overlap(flank_union).df
    if len(ov):
        parts.append(ov.drop(columns=["Start", "End"]))
    print(f"  MAF rows read: {n:,}", flush=True)
fs = pd.concat(parts, ignore_index=True).drop_duplicates(KEY).reset_index(drop=True)
fs["ctx"] = snv_ctx(fs)

# territory SNVs from step 02 hits (re-key with Project_Code for multi-project donors)
hits = pd.read_csv(RES / "ncorf_hits.tsv", sep="\t", dtype={"Chromosome": str})
ts = hits[(hits.Variant_Type == "SNP") & hits.Chromosome.isin(CHROMS)]
ts = ts.drop_duplicates(KEY + ["name"]).reset_index(drop=True)
ts["ctx"] = snv_ctx(ts)


def check(df, label):
    rc = df.ref_context.str.upper().str[9:12]
    ok = rc.map(lambda t: CTX_IDX.get(collapse(t), -2) if set(t) <= set(BASES) else -2)
    agree = (ok == df.ctx).mean()
    print(f"{label}: {len(df):,} SNVs; ctx valid {(df.ctx >= 0).mean():.4%}; "
          f"agree with ref_context {agree:.4%}")


check(fs, "flank SNVs")
check(ts, "territory SNVs")
fs["ctx_name"] = fs.ctx.map(dict(enumerate(CTX_NAMES)))
ts["ctx_name"] = ts.ctx.map(dict(enumerate(CTX_NAMES)))
fs.to_csv(DER / "flank_snvs.tsv.gz", sep="\t", index=False)
ts.to_csv(DER / "territory_snvs.tsv", sep="\t", index=False)

# per-ORF flank observed counts per project (an SNV can sit in several ORFs' flanks)
fj = pr.PyRanges(fs[fs.ctx >= 0].assign(Start=fs.Start_position - 1, End=fs.Start_position)) \
    .join(pr.PyRanges(flanks[["Chromosome", "Start", "End", "name"]])).df
fj.groupby(["name", "Project_Code"]).size().rename("flank_obs").to_csv(
    DER / "flank_obs_by_orf_project.tsv", sep="\t")

print(f"\nORFs with territory (no chrY): {len(names):,}")
print(f"territory bases w/ valid ctx: {sum(v.sum() for v in terr_ctx.values()):,}")
print(f"flank per-ORF bases (sum): {sum(v.sum() for v in flank_ctx.values()):,}; "
      f"flank union bases: {union_ctx.sum():,}")
print(f"flank SNVs (dedup): {len(fs):,}; territory SNV-ORF pairs: {len(ts):,}")
print(f"flank SNVs per Mb union per project:\n"
      f"{(fs.groupby('Project_Code').size() / union_ctx.sum() * 1e6).round(1).to_string()}")
