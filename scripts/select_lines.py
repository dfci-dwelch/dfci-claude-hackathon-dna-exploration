"""Select ~7 WGS lines per lineage (12 lineages) with 24Q4 CRISPR data, sex known, sex-balanced, seed 7.
Always include the 5 slice-1 lines. Writes data/meta/selected_lines.tsv."""
import pandas as pd, numpy as np
inv = pd.read_csv("data/meta/inventory.csv"); w = inv[inv.datatype == "wgs"].copy()
crispr = set(pd.read_csv("data/depmap/CRISPRGeneEffect.csv", usecols=[0]).iloc[:, 0])
sig = pd.read_csv("data/depmap/OmicsSignatures.csv", index_col=0)
eng = w.cell_line_name.str.contains(r"_|\bKO\b|\bKD\b", regex=True)
print("WGS:", len(w), "with CRISPR:", w.model_id.isin(crispr).sum(), "with signatures:", w.model_id.isin(sig.index).sum(), "engineered-looking names:", eng.sum(), w[eng].cell_line_name.tolist()[:8])
LIN = ["Lung", "Lymphoid", "CNS/Brain", "Skin", "Breast", "Esophagus/Stomach", "Myeloid", "Ovary/Fallopian Tube", "Bowel", "Pancreas", "Kidney", "Head and Neck"]
print("lineage names present:", [l for l in LIN if l not in set(w.oncotree_lineage)], "<- missing")
pool = w[w.model_id.isin(crispr) & ~eng & w.sex.isin(["male", "female"])]
rng = np.random.default_rng(7); base = ["HOP-62", "SW 1573", "HCC2935", "A427", "U-2 OS"]
sel = [w[w.cell_line_name == b].iloc[0] for b in base]
for lin in LIN:
    p = pool[(pool.oncotree_lineage == lin) & ~pool.cell_line_name.isin(base)]
    picks = []
    for sex in ["female", "male"]:
        q = p[p.sex == sex]; picks.append(q.iloc[rng.permutation(len(q))])
    # alternate sexes, fill from whichever has more
    f, m = picks; out = []
    for i in range(7):
        for d in ([f, m] if i % 2 == 0 else [m, f]):
            d2 = d[~d.model_id.isin([o.model_id for o in out])]
            if len(d2): out.append(d2.iloc[0]); break
    sel += out
S = pd.DataFrame(sel)[["cell_line_name", "model_id", "sex", "oncotree_lineage", "oncotree_primary_disease", "primary_or_metastasis", "cram_bam_size"]]
S["has_crispr"] = S.model_id.isin(crispr); S["has_sig"] = S.model_id.isin(sig.index)
S.to_csv("data/meta/selected_lines.tsv", sep="\t", index=False)
print(len(S), "lines;", S.sex.value_counts().to_dict(), "; total CRAM/BAM", round(S.cram_bam_size.sum() / 1e12, 1), "TB (not downloaded)")
print(S.groupby("oncotree_lineage").sex.value_counts().unstack(fill_value=0).to_string())
print("missing CRISPR among selected:", S[~S.has_crispr].cell_line_name.tolist(), "missing signatures:", S[~S.has_sig].cell_line_name.tolist())
