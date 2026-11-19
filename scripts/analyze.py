"""Slice-3 analysis on region_v2 results joined with DepMap 24Q4 covariates.

Derived metrics (per line):
  r45_raw, r5_raw   : rDNA units per genome assuming median control depth = 2 copies (summary units_best:
                      GC-matched median for 45S, plain median for 5S; per-line read length)
  batch             : 101 bp (CCLE BAMs) vs 151 bp (CRAMs); added to every adjusted model
  r45_acro          : 45S units per acrocentric chromosome copy  = r45_raw * P_ctrl/2 / acro_CN  (q-arm CN proxy)
  r45_dj            : 45S units per distal-junction copy = r45_raw / DJ_gc (DJ counted directly, same normalization,
                      so ploidy and the control-depth assumption cancel). Preferred 45S metric. Needs dj_summary.py.
  DJ_abs            : DJ copies per cell (~10 in a normal diploid); compared with acro_CN
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
w = S.pivot_table(index="model_id", columns="group", values="units_best")[["rDNA_45S", "rDNA_5S"]]
w.columns = ["r45_raw", "r5_raw"]
y = S[S.group.str.startswith("Y_")].pivot_table(index="model_id", columns="group", values="rel_CHM13")
d = sel.set_index("model_id").join(w, how="inner").join(y).join(cov.set_index("model_id"))
d["depth"] = S.groupby("model_id").depth.first()
d["read_len"] = S.groupby("model_id").read_len.first()
d["batch"] = np.where(d.read_len < 120, "bp101", "bp151")
d["r45_acro"] = d.r45_raw * d.P_ctrl / 2 / d.acro_CN
d["r5_1q"] = d.r5_raw * d.P_ctrl / 2 / d.cn_1q42
d["log2_5S_45S"] = np.log2(d.r5_1q / d.r45_acro)
for c in [c for c in d if c.startswith("Y_")]: d[c + "_adj"] = d[c] * d.P_ctrl / 2
if os.path.exists("results/dj_summary_dj_region.tsv"):
    dj = pd.read_csv("results/dj_summary_dj_region.tsv", sep="\t").set_index("model_id")[["DJ_gc", "DJ_abs"]]
    d = d.join(dj)
else:
    d["DJ_gc"] = d["DJ_abs"] = np.nan
d["r45_dj"] = d.r45_raw / d.DJ_gc
d["log2_5S_45S_dj"] = np.log2(d.r5_1q / d.r45_dj)
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
for c in ["r45_raw", "r45_acro", "r45_dj", "r5_raw", "r5_1q", "log2_5S_45S", "log2_5S_45S_dj", "DJ_abs"]:
    q = d[c].quantile([.25, .5, .75]); print(f"  {c:14s} {q[.5]:8.2f} [{q[.25]:.2f}, {q[.75]:.2f}]  n={d[c].notna().sum()}")
z = d.dropna(subset=["DJ_abs", "acro_CN"])
if len(z) > 5:
    print(f"\nDJ (direct p-arm count) vs q-arm proxy (n={len(z)}): spearman {stats.spearmanr(z.DJ_abs, z.acro_CN)[0]:+.3f}; "
          f"median DJ_abs/acro_CN {np.median(z.DJ_abs / z.acro_CN):.2f}; lines with DJ_abs < 0.75*acro_CN (p-arm loss): "
          f"{int((z.DJ_abs < 0.75 * z.acro_CN).sum())}")
    print(f"  r45_dj vs r45_acro spearman {stats.spearmanr(z.r45_dj, z.r45_acro)[0]:+.3f}; r45_dj vs r45_raw {stats.spearmanr(z.r45_dj, z.r45_raw)[0]:+.3f}")
print("\nraw vs CN-corrected spearman: 45S", np.round(stats.spearmanr(d.r45_raw, d.r45_acro, nan_policy="omit"), 4),
      " 5S", np.round(stats.spearmanr(d.r5_raw, d.r5_1q, nan_policy="omit"), 4))
print("\nsequencing batch (n):", d.batch.value_counts().to_dict())
for c in ["r45_dj", "r45_acro", "r5_1q", "DJ_abs", "Y_DYZ3_aSat_adj"]:
    sub = d.dropna(subset=[c])
    if len(sub) < 20: continue; sub = sub[sub[c] > 0] if c.startswith("Y_") else sub
    if c.startswith("Y_"): sub = sub[(sub.sex == "male") & (sub[c] >= 0.2)]
    m = smf.ols(f"np.log2({c}) ~ C(batch) + C(lineage)" + ("" if c.startswith("Y_") else " + C(sex)"), data=sub).fit()
    k = [i for i in m.params.index if i.startswith("C(batch)")][0]
    print(f"  residual batch effect on {c} (151 vs 101 bp, log2, lineage-adjusted): {m.params[k]:+.3f}  p={m.pvalues[k]:.3g}  n={int(m.nobs)}")

# --- Wang & Lemos replication: 5S vs 45S, and TP53 ---
cc = d.dropna(subset=["r45_acro", "r5_1q"])
print("\n5S vs 45S spearman  raw:", np.round(stats.spearmanr(d.r5_raw, d.r45_raw), 4),
      " CN-corrected:", np.round(stats.spearmanr(cc.r5_1q, cc.r45_acro), 4))
cj = d.dropna(subset=["r45_dj", "r5_1q"])
if len(cj) > 10: print("5S vs 45S spearman, 45S per DJ copy:", np.round(stats.spearmanr(cj.r5_1q, cj.r45_dj), 4), f"n={len(cj)}")
rows = []
# contrasts: mutation only; p53 loss (mutation, homozygous deletion, or no CRISPR-measured p53 activity)
for status in ["TP53_mut", "p53_loss"]:
    for c in ["r45_dj", "r45_acro", "r5_1q", "log2_5S_45S_dj", "log2_5S_45S", "r45_raw", "r5_raw"]:
        dd_ = d.dropna(subset=[c, status]).assign(st=lambda x: x[status].astype(bool))
        if len(dd_) < 10: continue
        a = dd_.loc[dd_.st, c]; b = dd_.loc[~dd_.st, c]
        u = stats.mannwhitneyu(a, b)
        lhs = c if c.startswith("log2_") else f"np.log2({c})"
        m = smf.ols(f"{lhs} ~ C(st) + C(lineage) + C(sex) + C(batch)", data=dd_).fit()
        k = [i for i in m.params.index if i.startswith("C(st)")][0]
        rows.append(dict(status=status, metric=c, n_loss=len(a), n_intact=len(b), median_loss=a.median(), median_intact=b.median(),
                         MWU_p=u.pvalue, ols_log2_effect=m.params[k], ols_p=m.pvalues[k]))
tp = pd.DataFrame(rows); tp.to_csv("results/analysis/tp53_tests.tsv", sep="\t", index=False)
print("\nTP53 loss vs intact (OLS adjusts lineage + sex + batch; log2 scale):"); print(tp.round(4).to_string(index=False))
print("continuous p53 activity (p53_score = TP53 effect - MDM2 effect), spearman:")
for c in ["r45_dj", "r45_acro", "r5_1q", "log2_5S_45S_dj"]:
    z = d.dropna(subset=[c, "p53_score"])
    if len(z) < 10: continue; r_, p_ = stats.spearmanr(z.p53_score, z[c]); print(f"  {c:12s} rho={r_:+.3f} p={p_:.2g} n={len(z)}")

# --- by lineage ---
lin = d.groupby("oncotree_lineage")[["r45_dj", "r45_acro", "r5_1q", "log2_5S_45S_dj"]].median().round(2)
lin["n"] = d.groupby("oncotree_lineage").size(); lin.to_csv("results/analysis/by_lineage.tsv", sep="\t")
print("\nBy lineage (medians):"); print(lin.sort_values("r45_dj").to_string())
for c in ["r45_dj", "r45_acro", "r5_1q"]:
    grp = [g[c].dropna() for _, g in d.groupby("lineage") if g[c].notna().sum() >= 3]
    if len(grp) > 2: print(f"Kruskal {c} across lineages:", np.round(stats.kruskal(*grp), 4))

# --- CRISPR: candidate genes + genome-wide scan ---
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
dd = d.join(ge, how="inner"); print(f"\nCRISPR: {len(dd)} lines with gene effect")
cand = ["POLR1A", "POLR1B", "POLR1C", "POLR1D", "RRN3", "UBTF", "TAF1B", "TCOF1", "NCL", "NOP56", "FBL", "POLR3A", "POLR3B", "GTF3C1", "TP53", "MDM2", "MYC"]
res = []
for metric in ["r45_dj", "r45_acro", "r5_1q", "log2_5S_45S_dj"]:
    x = dd[metric]; ok = x.notna()
    if ok.sum() < 30: continue
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
