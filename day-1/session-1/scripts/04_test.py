"""Step 5b: per-ORF SNV excess test vs. a trinucleotide-context + local-flank background.

expected(ORF, project) = scale * sum_ctx territory_bases[ctx] * rate[project, ctx]
  rate[project, ctx] = flank-union SNVs in ctx / flank-union bases in ctx  (all donors kept)
  scale = (flank_obs + K) / (flank_exp + K)   (local mutation-rate adjustment, shrunk to 1)
Test: one-sided Poisson P(X >= obs), BH FDR per project and pan-cancer (sum of per-project
expectations). Poisson ignores overdispersion between donors/loci -> anti-conservative.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests

ROOT = Path(sys.argv[1]).resolve()
DER, RES = ROOT / "data/derived", ROOT / "results"
K = 10

tctx = pd.read_csv(DER / "ctx_bases_territory.tsv", sep="\t", index_col="name")
fctx = pd.read_csv(DER / "ctx_bases_flank.tsv", sep="\t", index_col="name").reindex(tctx.index, fill_value=0)
uctx = pd.read_csv(DER / "ctx_bases_flank_union.tsv", sep="\t", index_col=0)["bases"]
fs = pd.read_csv(DER / "flank_snvs.tsv.gz", sep="\t", usecols=["Project_Code", "ctx_name"])
ts = pd.read_csv(DER / "territory_snvs.tsv", sep="\t", dtype={"Chromosome": str})
fobs = pd.read_csv(DER / "flank_obs_by_orf_project.tsv", sep="\t")

rate = (fs.groupby(["Project_Code", "ctx_name"]).size().unstack(fill_value=0)
        .reindex(columns=uctx.index, fill_value=0) / uctx)
rate.to_csv(RES / "context_rates.tsv", sep="\t")
projects = rate.index

terr_exp = pd.DataFrame(tctx.values @ rate.T.values, index=tctx.index, columns=projects)
flank_exp = pd.DataFrame(fctx.values @ rate.T.values, index=tctx.index, columns=projects)
flank_obs = (fobs.pivot(index="name", columns="Project_Code", values="flank_obs")
             .reindex(index=tctx.index, columns=projects).fillna(0))
scale = (flank_obs + K) / (flank_exp + K)
expected = terr_exp * scale

obs = ts.groupby(["name", "Project_Code"]).size().unstack(fill_value=0) \
    .reindex(index=tctx.index, columns=projects, fill_value=0)
don = ts.groupby(["name", "Project_Code"])["Donor_ID"].nunique().unstack(fill_value=0) \
    .reindex(index=tctx.index, columns=projects, fill_value=0)

meta = ts.drop_duplicates("name").set_index("name")[["orf_biotype", "peptide_support"]]
info = pd.DataFrame(index=tctx.index)
info["gene"] = info.index.str.split(":").str[1]
info = info.join(meta)
info["terr_bp"] = tctx.sum(axis=1)


def per_orf_extras(sub):
    g = sub.groupby("name")
    return pd.DataFrame({
        "positions": g.Start_position.nunique(),
        "max_per_donor": g.Donor_ID.agg(lambda s: s.value_counts().max()),
        "top_projects": g.Project_Code.agg(lambda s: ",".join(f"{p}:{n}" for p, n in s.value_counts().head(3).items())),
    })


# ---- per project ----
rows = []
for p in projects:
    d = pd.DataFrame({"project": p, "obs": obs[p], "donors": don[p], "expected": expected[p],
                      "local_scale": scale[p]})
    d["p"] = poisson.sf(d.obs - 1, d.expected)
    d["q"] = multipletests(d.p, method="fdr_bh")[1]
    rows.append(d)
byp = pd.concat(rows).join(info)
byp = byp.join(per_orf_extras(ts).rename(columns=lambda c: c), how="left")  # pan-wide extras
byp["obs_exp"] = byp.obs / byp.expected
byp = byp.sort_values("p")
byp.to_csv(RES / "anomaly_by_project.tsv", sep="\t")

# ---- pan-cancer ----
pan = info.copy()
pan["obs"] = obs.sum(axis=1)
pan["donors"] = ts.groupby("name")["Donor_ID"].nunique().reindex(pan.index, fill_value=0)
pan["expected"] = expected.sum(axis=1)
pan["obs_exp"] = pan.obs / pan.expected
pan["p"] = poisson.sf(pan.obs - 1, pan.expected)
pan["q"] = multipletests(pan.p, method="fdr_bh")[1]
pan = pan.join(per_orf_extras(ts)).sort_values("p")
pan.to_csv(RES / "anomaly_pancancer.tsv", sep="\t")

# ---- report ----
print(f"ORFs tested: {len(pan):,}; territory SNVs observed: {int(pan.obs.sum()):,}; "
      f"expected: {pan.expected.sum():,.1f} (global obs/exp {pan.obs.sum() / pan.expected.sum():.2f})")
print(f"local scale (pan, flank_obs/flank_exp weighted): median per ORF-project "
      f"{np.nanmedian(scale.values):.2f}")
cols = ["gene", "orf_biotype", "obs", "donors", "positions", "max_per_donor", "expected",
        "obs_exp", "q", "top_projects"]
for thr in (0.05, 0.1):
    print(f"pan-cancer q<{thr}: {(pan.q < thr).sum()}; per-project ORF-project pairs q<{thr}: "
          f"{(byp.q < thr).sum()}")
fmt = {"expected": "{:.2f}".format, "obs_exp": "{:.1f}".format, "q": "{:.2g}".format}
print("\nPan-cancer top 25:")
print(pan.head(25)[cols].to_string(formatters=fmt))
print("\nPer-project q<0.1:")
print(byp[byp.q < 0.1][["project"] + cols[:-1] + ["local_scale"]].to_string(
    formatters={**fmt, "local_scale": "{:.2f}".format}))
