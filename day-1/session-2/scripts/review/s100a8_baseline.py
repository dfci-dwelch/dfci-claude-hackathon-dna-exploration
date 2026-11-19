"""Is S100A8 just tracking a per-line baseline shared by non-expressed / non-essential genes?
Baseline proxies: per-line median gene effect over (a) olfactory-receptor genes (never expressed), (b) all genes,
and the per-line spread (SD) of gene effects. If r5_1q tracks a baseline, S100A8 association is a screen-level artefact."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
d = pd.read_csv("results/analysis/per_line.tsv", sep="\t", index_col=0)
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0); ge.columns = [c.split(" ")[0] for c in ge.columns]
ge = ge.loc[ge.index.intersection(d.index)]
orx = [c for c in ge.columns if c.startswith("OR") and c[2:3].isdigit()]
base = pd.DataFrame({"OR_median": ge[orx].median(1), "all_median": ge.median(1), "all_sd": ge.std(1),
                     "S100A8": ge["S100A8"]})
x = d.join(base, how="inner").dropna(subset=["r5_1q"])
print(f"n={len(x)}; olfactory-receptor genes used: {len(orx)}")
for a, b in [("r5_1q", "OR_median"), ("r5_1q", "all_median"), ("r5_1q", "all_sd"), ("S100A8", "OR_median"), ("S100A8", "all_sd")]:
    r, p = stats.spearmanr(x[a], x[b], nan_policy="omit"); print(f"  spearman({a}, {b}) = {r:+.3f} p={p:.2g}")
m0 = smf.ols("S100A8 ~ np.log2(r5_1q)", data=x).fit()
m1 = smf.ols("S100A8 ~ np.log2(r5_1q) + OR_median + all_sd", data=x).fit()
print(f"  S100A8 ~ log2 r5_1q: beta={m0.params.iloc[1]:+.4f} p={m0.pvalues.iloc[1]:.2g};  + baseline terms: beta={m1.params.iloc[1]:+.4f} p={m1.pvalues.iloc[1]:.2g}")
# how many genes would 'baseline' produce? genome-wide: are top r5 hits enriched for non-expressed-type genes (OR / S100 / immune)?
scan = pd.read_csv("results/analysis/crispr_scan.tsv", sep="\t"); s = scan[scan.metric == "r5_1q"].copy()
s["is_OR"] = s.gene.isin(orx)
print("  genome-wide r5_1q scan: median rho for OR genes %.3f vs all genes %.3f; fraction OR with rho>0: %.2f" %
      (s[s.is_OR].rho.median(), s.rho.median(), (s[s.is_OR].rho > 0).mean()))
# S100A8 residual after baseline: which lines drive it?
x["resid"] = m1.resid
print(x.sort_values("S100A8", ascending=False)[["cell_line_name", "oncotree_lineage", "r5_1q", "S100A8", "OR_median"]].head(6).round(3).to_string())
