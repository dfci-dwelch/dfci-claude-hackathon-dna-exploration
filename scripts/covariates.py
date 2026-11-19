"""Per-model copy-number and mutation covariates from DepMap 24Q4.

Absolute CN segments (prefer WGS profile, else WES):
  P_ctrl  : median absolute CN over our 200 control windows (the CN that our depth normalization assumes = 2)
  acro_CN : sum over chr13/14/15/21/22 of the length-weighted CN of the proximal 10 Mb of the q arm
            (closest assayed sequence to the p-arm rDNA arrays; p arms themselves are not segmented)
  cn_1q42 : CN at the 5S array, chr1:228.0-229.0 Mb
TP53:
  TP53_mut   : damaging (OmicsSomaticMutationsMatrixDamaging) or hotspot (MatrixHotspot) > 0
  TP53_cn, TP53_LOH : absolute CN and LOH status at the TP53 locus (chr17:7.66-7.69 Mb)
  p53_score  : CRISPR TP53 gene effect - MDM2 gene effect (high = p53 pathway active: TP53 KO gives a growth
               advantage and MDM2 is essential). Calibrated on all DepMap 24Q4 lines with CRISPR data.
  p53_functional : p53_score > 95th percentile of TP53-mutant lines (threshold 0.95; 72% of TP53-WT lines pass).
               NaN for lines without CRISPR data
  p53_loss   : TP53_mut OR TP53_cn < 0.5 (homozygous deletion) OR p53_functional == False
Also ploidy / WGD / aneuploidy from OmicsSignatures.
"""
import numpy as np, pandas as pd

CEN_END = {"13": 18_100_000, "14": 18_200_000, "15": 19_800_000, "21": 12_950_000, "22": 15_100_000}  # hg38 approx.


def wcn(seg, chrom, a, b):
    s = seg[(seg.Chromosome == chrom) & (seg.End > a) & (seg.Start < b)]
    if s.empty: return np.nan
    L = np.minimum(s.End, b) - np.maximum(s.Start, a)
    return float((s.SegmentAbsoluteCN * L).sum() / L.sum())


def main():
    seg = pd.read_csv("data/depmap/OmicsAbsoluteCNSegmentsProfile.csv", dtype={"Chromosome": str})
    prof = pd.read_csv("data/depmap/OmicsProfiles.csv")
    prof = prof[prof.ProfileID.isin(seg.ProfileID)].copy()
    prof["rank"] = prof.Datatype.map({"wgs": 0, "wes": 1})
    pick = prof.sort_values("rank").drop_duplicates("ModelID").set_index("ModelID")
    ctrl = pd.read_csv("data/kmers/control_windows.tsv", sep="\t")
    sel = pd.read_csv("data/meta/selected_lines.tsv", sep="\t")
    sig = pd.read_csv("data/depmap/OmicsSignatures.csv", index_col=0)
    dmg = pd.read_csv("data/depmap/OmicsSomaticMutationsMatrixDamaging.csv", index_col=0, usecols=lambda c: c == "" or c.startswith("Unnamed") or c.startswith("TP53 "))
    hot = pd.read_csv("data/depmap/OmicsSomaticMutationsMatrixHotspot.csv", index_col=0, usecols=lambda c: c == "" or c.startswith("Unnamed") or c.startswith("TP53 "))
    # functional p53 call, calibrated on every DepMap line with CRISPR data
    ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0,
                     usecols=lambda c: c == "" or c.startswith("Unnamed") or c.split(" ")[0] in ("TP53", "MDM2"))
    ge.columns = [c.split(" ")[0] for c in ge.columns]
    allmut = (dmg.iloc[:, 0].reindex(ge.index).fillna(0) + hot.iloc[:, 0].reindex(ge.index).fillna(0)) > 0
    score = ge.TP53 - ge.MDM2
    thr = float(score[allmut].quantile(0.95))
    print(f"p53_score threshold (95th pct of {int(allmut.sum())} TP53-mutant lines): {thr:.3f}; "
          f"WT lines above it: {(score[~allmut] > thr).mean():.2%}")
    rows = []
    for m in sel.model_id:
        r = {"model_id": m}
        if m in pick.index:
            s = seg[seg.ProfileID == pick.loc[m, "ProfileID"]]
            r["cn_source"] = pick.loc[m, "Datatype"]
            r["P_ctrl"] = np.nanmedian([wcn(s, c.replace("chr", ""), a, b) for c, a, b in ctrl[["contig", "start", "end"]].values])
            per = {c: wcn(s, c, e, e + 10_000_000) for c, e in CEN_END.items()}
            r.update({f"cn_{c}q": v for c, v in per.items()}); r["acro_CN"] = sum(per.values())
            r["cn_1q42"] = wcn(s, "1", 228_000_000, 229_000_000)
            r["TP53_cn"] = wcn(s, "17", 7_661_779, 7_687_538)
            t = s[(s.Chromosome == "17") & (s.Start < 7_687_538) & (s.End > 7_661_779)]
            r["TP53_LOH"] = bool((t.LoHStatus == "LOH").any()) if len(t) else np.nan
        if m in sig.index: r.update(sig.loc[m, ["Ploidy", "WGD", "Aneuploidy"]].to_dict())
        tp = [d.loc[m].iloc[0] for d in (dmg, hot) if m in d.index and d.shape[1]]
        r["TP53_mut"] = (np.nansum(tp) > 0) if tp else np.nan
        if m in score.index and pd.notna(score[m]):
            r["p53_score"] = float(score[m]); r["p53_functional"] = bool(score[m] > thr)
        r["p53_loss"] = bool(r.get("TP53_mut") is True or r.get("TP53_mut") == True or
                             (pd.notna(r.get("TP53_cn")) and r["TP53_cn"] < 0.5) or r.get("p53_functional") is False)
        rows.append(r)
    out = pd.DataFrame(rows); out.to_csv("data/meta/covariates.tsv", sep="\t", index=False)
    print(out.describe(include="all").T[["count", "mean", "min", "max"]].to_string())


if __name__ == "__main__":
    main()
