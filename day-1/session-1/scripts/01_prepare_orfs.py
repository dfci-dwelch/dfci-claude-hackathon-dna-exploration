"""Expand GENCODE Phase 1 Ribo-seq ORFs (GRCh38 BED12) into exon blocks, lift to GRCh37,
and subtract canonical CDS (GENCODE v35lift37) to get ncORF-only territory."""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr
from pyliftover import LiftOver

ROOT = Path(sys.argv[1]).resolve()
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "derived"
OUT.mkdir(parents=True, exist_ok=True)

bed = pd.read_csv(RAW / "gencode" / "Ribo-seq_ORFs.all.bed", sep="\t", header=None,
                  names=["chrom", "start", "end", "name", "score", "strand", "thickStart",
                         "thickEnd", "rgb", "blockCount", "blockSizes", "blockStarts"])
bed["orf_id"] = bed["name"].str.split(":").str[0]
bed["gene"] = bed["name"].str.split(":").str[1]

biotype = {}
for line in open(RAW / "gencode" / "Ribo-seq_ORFs.all.faa"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        biotype[name] = (re.search(r"orf_biotype=(\S+)", line).group(1),
                         re.search(r"peptide_support=(\S+)", line).group(1))
bed["orf_biotype"] = bed["name"].map(lambda n: biotype[n][0])
bed["peptide_support"] = bed["name"].map(lambda n: biotype[n][1])

# Blocks in transcript (5'->3') order so orf_offset is the position within the spliced ORF.
rows = []
for r in bed.itertuples():
    sizes = [int(x) for x in r.blockSizes.strip(",").split(",")]
    starts = [int(x) for x in r.blockStarts.strip(",").split(",")]
    blocks = [(r.start + s, r.start + s + z) for s, z in zip(starts, sizes)]
    if r.strand == "-":
        blocks = blocks[::-1]
    offset = 0
    for i, (bs, be) in enumerate(blocks):
        rows.append(dict(name=r.name, orf_id=r.orf_id, gene=r.gene, orf_biotype=r.orf_biotype,
                         peptide_support=r.peptide_support, strand=r.strand, block_idx=i,
                         chrom38=r.chrom, start38=bs, end38=be, orf_offset=offset))
        offset += be - bs
blocks38 = pd.DataFrame(rows)

lo = LiftOver(str(RAW / "ucsc" / "hg38ToHg19.over.chain.gz"))


def lift_block(chrom, start, end, strand):
    a = lo.convert_coordinate(chrom, start)
    b = lo.convert_coordinate(chrom, end - 1)
    if not a or not b or len(a) != 1 or len(b) != 1:
        return None
    (ca, pa, sa, _), (cb, pb, sb, _) = a[0], b[0]
    if ca != cb or sa != sb:
        return None
    lo_, hi_ = sorted((pa, pb))
    if hi_ - lo_ + 1 != end - start:
        return None
    new_strand = strand if sa == "+" else {"+": "-", "-": "+"}[strand]
    return ca, lo_, hi_ + 1, new_strand


lifted = blocks38.apply(lambda r: lift_block(r.chrom38, r.start38, r.end38, r.strand), axis=1)
blocks38["ok"] = lifted.notna()
ok_by_orf = blocks38.groupby("name")["ok"].all()
report = bed[["name", "orf_id", "gene", "orf_biotype", "peptide_support"]].copy()
report["lift_ok"] = report["name"].map(ok_by_orf)
report.to_csv(OUT / "orf_lift_report.tsv", sep="\t", index=False)

good = blocks38[blocks38["name"].map(ok_by_orf)].copy()
good[["chrom37", "start37", "end37", "strand37"]] = pd.DataFrame(
    lifted[good.index].tolist(), index=good.index)
good["chrom37"] = good["chrom37"].str.replace("^chr", "", regex=True)
good = good.drop(columns=["ok"])
good.to_csv(OUT / "orf_blocks_grch37.tsv", sep="\t", index=False)

cds = []
with gzip.open(RAW / "gencode_annot" / "gencode.v35lift37.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if f[2] != "CDS" or 'transcript_type "protein_coding"' not in f[8]:
            continue
        cds.append((f[0].replace("chr", ""), int(f[3]) - 1, int(f[4])))
cds_gr = pr.PyRanges(pd.DataFrame(cds, columns=["Chromosome", "Start", "End"])).merge()

blk_gr = pr.PyRanges(good.rename(columns={"chrom37": "Chromosome", "start37": "Start",
                                          "end37": "End"})[
    ["Chromosome", "Start", "End", "name", "orf_biotype", "peptide_support"]])
terr = blk_gr.subtract(cds_gr).df
terr.to_csv(OUT / "ncorf_territory_grch37.tsv", sep="\t", index=False)

n_orf, n_ok = len(bed), int(ok_by_orf.sum())
print(f"ORFs: {n_orf}; lifted intact: {n_ok} ({n_ok / n_orf:.1%})")
print(f"Lifted ORF bp: {int((good.end37 - good.start37).sum()):,}; after CDS subtraction: "
      f"{int((terr.End - terr.Start).sum()):,} (merged: {int(pr.PyRanges(terr).merge().length):,})")
print(f"ORFs retaining territory: {terr['name'].nunique()}")
print((terr.assign(bp=terr.End - terr.Start).groupby("orf_biotype")["bp"].sum()
       .sort_values(ascending=False)).to_string())
lost = set(good["name"]) - set(terr["name"])
print(report[report["name"].isin(lost)]["orf_biotype"].value_counts().to_string())
