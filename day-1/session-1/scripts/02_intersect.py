"""Stream the PCAWG ICGC public MAF, compute per-donor burden, intersect with ncORF territory."""
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr

ROOT = Path(sys.argv[1]).resolve()
MAF = ROOT / "data/raw/pcawg/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
DER, RES = ROOT / "data/derived", ROOT / "results"
RES.mkdir(exist_ok=True)

COLS = ["Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification",
        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
        "Donor_ID", "Project_Code", "i_VAF", "ref_context"]
KEY = ["Donor_ID", "Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2"]

terr = pd.read_csv(DER / "ncorf_territory_grch37.tsv", sep="\t", dtype={"Chromosome": str})
terr_gr = pr.PyRanges(terr)
terr_bp = int(terr_gr.merge().length)

burden_parts, hit_parts, n_rows = [], [], 0
for chunk in pd.read_csv(MAF, sep="\t", usecols=COLS, chunksize=2_000_000,
                         dtype={"Chromosome": str}, low_memory=False):
    n_rows += len(chunk)
    chunk = chunk.drop_duplicates(KEY)  # collapse multi-sample donors
    burden_parts.append(chunk.groupby(["Donor_ID", "Project_Code", "Variant_Type"]).size())
    g = pr.PyRanges(chunk.assign(Start=chunk.Start_position - 1, End=chunk.End_position))
    j = g.join(terr_gr, suffix="_orf").df
    if len(j):
        hit_parts.append(j)
    print(f"  rows read: {n_rows:,}", flush=True)

burden = (pd.concat(burden_parts).groupby(level=[0, 1, 2]).sum()
          .unstack("Variant_Type", fill_value=0).reset_index())
burden["total"] = burden.drop(columns=["Donor_ID", "Project_Code"]).sum(axis=1)
burden.to_csv(RES / "donor_burden.tsv", sep="\t", index=False)

hits = pd.concat(hit_parts, ignore_index=True)
hits = hits.drop(columns=["Start", "End", "Start_orf", "End_orf"]).drop_duplicates(KEY + ["name"])
hits.to_csv(RES / "ncorf_hits.tsv", sep="\t", index=False)

uniq = hits.drop_duplicates(KEY)
GENOME_BP = 2.8e9  # crude; step 5 replaces this with a proper background model
proj = burden.groupby("Project_Code").agg(donors=("Donor_ID", "nunique"), total=("total", "sum"))
proj["ncorf_muts"] = uniq.groupby("Project_Code").size()
proj["donors_with_hit"] = uniq.groupby("Project_Code")["Donor_ID"].nunique()
proj = proj.fillna(0).astype(int)
proj["crude_rate_ratio"] = (proj.ncorf_muts / terr_bp) / (proj.total / GENOME_BP)
proj.sort_values("donors", ascending=False).to_csv(RES / "project_summary.tsv", sep="\t")

print(f"\nMAF rows: {n_rows:,}; donors: {burden.Donor_ID.nunique()}; projects: {len(proj)}")
print(f"territory: {terr_bp:,} bp; unique muts in territory: {len(uniq):,}; overlaps: {len(hits):,}")
print(f"donors with >=1 hit: {uniq.Donor_ID.nunique()}")
print(uniq["Variant_Type"].value_counts().to_string())
print(proj.to_string(float_format=lambda x: f"{x:.2f}"))
