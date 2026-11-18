"""Slice-3 analysis on region_v2 results joined with DepMap 24Q4 covariates.

Derived metrics (per line):
  r45_raw, r5_raw   : rDNA units per genome assuming median control depth = 2 copies (units_median)
  r45_acro          : 45S units per acrocentric chromosome copy  = r45_raw * P_ctrl/2 / acro_CN
  r5_1q             : 5S units per 1q42 copy                      = r5_raw  * P_ctrl/2 / cn_1q42
  log2_5S_45S       : log2(r5_1q / r45_acro)   (Wang & Lemos: 5S up / 45S down, tied to TP53 loss)
  Y_*               : array size vs HG002 Y; P_ctrl-corrected (x P_ctrl/2) so 1 ~ one Y per modal set
Outputs results/analysis/*.tsv and prints the key tests. Exploratory: no multiple-testing claims beyond BH.
"""
import os, numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

os.makedirs("results/analysis", exist_ok=True)
S = pd.read_csv("results/summary.tsv", sep="\t"); S = S[S["mode"] == "region_v2"]
sel = pd.read_csv("data/meta/selected_lines.tsv", sep="\t")
cov = pd.read_csv("data/meta/covariates.tsv", sep="\t")
name2id = dict(zip(sel.cell_line_name, sel.model_id))
S["model_id"] = S.line.map(name2id)
w = S.pivot_table(index="model_id", columns="group", values="units_median")[["rDNA_45S", "rDNA_5S"]]
w.columns = ["r45_raw", "r5_raw"]
y = S[S.group.str.startswith("Y_")].pivot_table(index="model_id", columns="group", values="rel_CHM13")
d = sel.set_index("model_id").join(w, how="inner").join(y).join(cov.set_index("model_id"))
d["depth"] = S.groupby("model_id").depth.first()
d["r45_acro"] = d.r45_raw * d.P_ctrl / 2 / d.acro_CN
d["r5_1q"] = d.r5_raw * d.P_ctrl / 2 / d.cn_1q42
d["log2_5S_45S"] = np.log2(d.r5_1q / d.r45_acro)
for c in [c for c in d if c.startswith("Y_")]: d[c + "_adj"] = d[c] * d.P_ctrl / 2
d["lineage"] = d.oncotree_lineage.str.replace(r"[^A-Za-z]", "", regex=True)
d.to_csv("results/analysis/per_line.tsv", sep="\t")
print(f"n lines = {len(d)}  (with CN covariates: {d.acro_CN.notna().sum()})")

# --- sex / LOY control ---
print("\nY arrays (P_ctrl-adjusted) by reported sex:")
print(d.groupby("sex")[["Y_DYZ1_HSat3_adj", "Y_DYZ2_HSat1B_adj", "Y_DYZ3_aSat_adj"]].agg(["count", "min", "median", "max"]).T.round(4).to_string())
fem_max = d.loc[d.sex == "female", "Y_DYZ3_aSat_adj"].max()
d["Y_status"] = np.select([d.sex == "female", d.Y_DYZ3_aSat_adj.isna(), d.Y_DYZ3_aSat_adj < 0.2], ["female", "male_no_CN", "male_LOY(<0.2)"], "male_Y_present")
print("female max DYZ3:", round(fem_max, 4), "| male Y status:", d[d.sex == "male"].Y_status.value_counts().to_dict())
print("males with LOY:", d[d.Y_status == "male_LOY(<0.2)"].cell_line_name.tolist())

# --- rDNA distributions ---
print("\nrDNA summaries (median [IQR]):")
for c in ["r45_raw", "r45_acro", "r5_raw", "r5_1q", "log2_5S_45S"]:
    q = d[c].quantile([.25, .5, .75]); print(f"  {c:12s} {q[.5]:8.2f} [{q[.25]:.2f}, {q[.75]:.2f}]")
print("\nCN correction matters? spearman(r45_raw, acro_CN*2/P_ctrl):",
      np.round(stats.spearmanr(d.r45_raw, d.acro_CN * 2 / d.P_ctrl, nan_policy="omit"), 4))

# --- Wang & Lemos replication: 5S vs 45S, and TP53 ---
cc = d.dropna(subset=["r45_acro", "r5_1q"])
print("\n5S vs 45S spearman  raw:", np.round(stats.spearmanr(d.r5_raw, d.r45_raw), 4),
      " CN-corrected:", np.round(stats.spearmanr(cc.r5_1q, cc.r45_acro), 4))
rows = []
for c in ["r45_acro", "r5_1q", "log2_5S_45S", "r45_raw", "r5_raw"]:
    a = d.loc[d.TP53_mut == True, c].dropna(); b = d.loc[d.TP53_mut == False, c].dropna()
    u = stats.mannwhitneyu(a, b)
    m = smf.ols(f"np.log2({c}) ~ C(TP53_mut) + C(lineage) + C(sex)" if c != "log2_5S_45S" else f"{c} ~ C(TP53_mut) + C(lineage) + C(sex)",
                data=d.dropna(subset=[c, "TP53_mut"]).assign(TP53_mut=lambda x: x.TP53_mut.astype(bool))).fit()
    k = [i for i in m.params.index if i.startswith("C(TP53_mut)")][0]
    rows.append(dict(metric=c, n_mut=len(a), n_wt=len(b), median_mut=a.median(), median_wt=b.median(), MWU_p=u.pvalue,
                     ols_log2_effect=m.params[k], ols_p=m.pvalues[k]))
tp = pd.DataFrame(rows); tp.to_csv("results/analysis/tp53_tests.tsv", sep="\t", index=False)
print("\nTP53 mutant vs WT (OLS adjusts lineage + sex; log2 scale):"); print(tp.round(4).to_string(index=False))

# --- by lineage ---
lin = d.groupby("oncotree_lineage")[["r45_acro", "r5_1q", "log2_5S_45S"]].median().round(2)
lin["n"] = d.groupby("oncotree_lineage").size(); lin.to_csv("results/analysis/by_lineage.tsv", sep="\t")
print("\nBy lineage (medians):"); print(lin.sort_values("r45_acro").to_string())
print("Kruskal r45_acro across lineages:", np.round(stats.kruskal(*[g.r45_acro.dropna() for _, g in d.groupby("lineage") if g.r45_acro.notna().sum() >= 3]), 4))

# --- CRISPR: candidate genes + genome-wide scan ---
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
dd = d.join(ge, how="inner"); print(f"\nCRISPR: {len(dd)} lines with gene effect")
cand = ["POLR1A", "POLR1B", "POLR1C", "POLR1D", "RRN3", "UBTF", "TAF1B", "TCOF1", "NCL", "NOP56", "FBL", "POLR3A", "POLR3B", "GTF3C1", "TP53", "MDM2", "MYC"]
res = []
for metric in ["r45_acro", "r5_1q", "log2_5S_45S"]:
    x = dd[metric]; ok = x.notna()
    X = dd.loc[ok, ge.columns]; xr = x[ok].rank()
    # Spearman with per-gene pairwise-complete n (many genes are NaN in a subset of screens)
    n = X.notna().sum()
    Xr = X.rank(); xr_g = Xr.notna().mul(x[ok], axis=0).where(Xr.notna())  # metric values restricted to each gene's lines
    rho = Xr.corrwith(xr_g.rank())
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2)); p = 2 * stats.t.sf(np.abs(t), n - 2)
    r = pd.DataFrame({"metric": metric, "gene": rho.index, "n": n.values, "rho": rho.values, "p": p}).dropna()
    r = r[r.n >= 30]  # require >= 30 lines with data
    r["q_BH"] = multipletests(r.p, method="fdr_bh")[1]; res.append(r)
R = pd.concat(res); R.to_csv("results/analysis/crispr_scan.tsv", sep="\t", index=False)
print("candidate genes (spearman rho, p):")
print(R[R.gene.isin(cand)].pivot_table(index="gene", columns="metric", values=["n", "rho", "p"]).round(3).to_string())
print("\ngenome-wide top hits per metric (q_BH):")
for m_, g in R.groupby("metric"): print(m_, g.nsmallest(5, "p")[["gene", "n", "rho", "p", "q_BH"]].round(4).to_string(index=False), sep="\n")
