import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import linregress, spearmanr
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Checking_Statistics_fig6.py
#
# Consolidated assumption checks for every MS-vs-cumulative-BMB regression
# behind Figure 6 (main text: Amundsen, Wilkes, no calving-group split) and
# Figure S ms_vs_bmb_other_sectors_bygroup (AIS, Ross, FilchnerRonne, Aurora,
# pooled and Group 2; Group 1 is excluded here since no fit line is shown for
# it and it is not used in this regression).
#
# For every case this computes:
#   (i)   linearity   - Pearson r vs Spearman rho on per-model-aggregated data
#                        (one point per model, cumulative BMB averaged across
#                        the 4 pooled experiments); a rho-r gap > 0.1 flags
#                        curvature the straight line is missing.
#   (ii)  homoscedasticity - Breusch-Pagan p-value on the same per-model fit;
#                        p<0.05 flags residual variance that changes with MS.
#   (iii) independence - cluster-robust p-value (clustered by model) on the
#                        full pooled (experiment-level) data, since each model
#                        contributes up to 4 points sharing one MS value; the
#                        naive OLS p-value is reported alongside for reference
#                        of how much it was overstated.
#
# Saves one combined table: reviewing/tables/Checking_Statistics_fig6.csv
# and one diagnostic residual-plot figure per case group:
#   reviewing/figures/Checking_Statistics_fig6_residuals.pdf
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT = os.path.join(PROJECT_ROOT, "ComputedScalars")
TABLES = os.path.join(PROJECT_ROOT, "analysis", "tables")
FIGDIR = os.path.join(PROJECT_ROOT, "reviewing", "figures")
OUTTABLES = os.path.join(PROJECT_ROOT, "reviewing", "tables")

EXPS = ['expAE02', 'expAE03', 'expAE04', 'expAE05']
MAIN_REGIONS = ['Amundsen', 'Wilkes']                              # Figure 6, no grouping
GROUPED_REGIONS = ['AIS', 'Ross', 'FilchnerRonne', 'Aurora']        # Figure S, by calving group
REGION_KEYS = {
    'Amundsen': ['_sector_4'], 'Wilkes': ['_sector_7'],
    'AIS': [''], 'Ross': ['_sector_2', '_sector_6'],
    'FilchnerRonne': ['_sector_5', '_sector_11'], 'Aurora': ['_sector_8'],
}
FACTOR = -1 * 365.25 * 24 * 60 * 60 / 10**12

# calving grouping (confirmed correct by Ronja), same list used throughout reviewing/
models_group1 = ['VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3', 'VUW_PISM1_s4', 'VUW_PISM2',
                 'VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2',
                 'LSCE_GRISLI', 'UCM_Yelmo', 'IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4']
GROUP = {m: 1 for m in models_group1}

df_best = pd.read_csv(os.path.join(TABLES, "best_ms_method_per_model.csv"))
models_in_table = set(df_best['Model'])

ALL_REGIONS = MAIN_REGIONS + GROUPED_REGIONS

# --- assemble all (model, experiment, region) points, matching fig6_pooled.py /
# fig6_SI_other_sectors_groups.py exactly (main regions keep VUW_PISM; grouped
# regions exclude it, as in the existing by-group robustness check) ---
all_recs = []
for exp in EXPS:
    ddir = os.path.join(DATA_ROOT, exp, "shelfmelt")
    files = sorted(f for f in os.listdir(ddir) if f.endswith(".nc"))
    for region in ALL_REGIONS:
        ms_col = f"melt_sensetivity_{region}"
        for f in files:
            model = f.replace("computed_shelfmelt_AIS_", "").replace(f"_{exp}.nc", "")
            if model not in models_in_table:
                continue
            if region in GROUPED_REGIONS and model.startswith("VUW_PISM"):
                continue  # excluded from the by-group analysis (see robustness check)
            ms = df_best.loc[df_best['Model'] == model, ms_col].values
            if len(ms) == 0 or np.isnan(ms[0]):
                continue
            d = xr.open_dataset(os.path.join(ddir, f), decode_times=False)
            y = sum(d[f"shelfmelt{k}"].values for k in REGION_KEYS[region]) * FACTOR
            time = d["time"].values
            d.close()
            cumulative_change = np.trapz(y - y[0], time)
            all_recs.append(dict(region=region, experiment=exp, model=model, ms=ms[0],
                                  bmb_change=cumulative_change, group=GROUP.get(model, 2)))

df_all = pd.DataFrame(all_recs)


def check_case(sub, label):
    """Run all three assumption checks on one (region, group) subset."""
    if sub.model.nunique() < 3:
        return None

    # (iii) independence: cluster-robust vs naive, on the full pooled (experiment-level) data
    X = sm.add_constant(sub['ms'])
    ols = sm.OLS(sub['bmb_change'], X).fit()
    p_naive = ols.pvalues['ms']
    clu = sm.OLS(sub['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['model']})
    p_cluster = clu.pvalues['ms']

    # (i) linearity and (ii) homoscedasticity: per-model-aggregated data (one point per
    # model; cumulative BMB averaged across the pooled experiments), so a single model's
    # repeated points don't distort the shape/variance diagnostics
    per_model = sub.groupby('model').agg(ms=('ms', 'first'), bmb_change=('bmb_change', 'mean')).reset_index()
    slope, intercept, r_value, p_lin, se = linregress(per_model.ms, per_model.bmb_change)
    rho, p_rho = spearmanr(per_model.ms, per_model.bmb_change)
    resid = per_model.bmb_change - (slope * per_model.ms + intercept)
    try:
        bp_stat, bp_p, _, _ = het_breuschpagan(resid, sm.add_constant(per_model['ms']))
    except Exception:
        bp_p = np.nan

    return dict(
        case=label, n_pts=len(sub), n_models=sub.model.nunique(),
        r=r_value, r2=r_value**2, rho=rho, gap=rho - r_value,
        bp_p=bp_p, p_naive=p_naive, p_cluster=p_cluster,
    ), per_model, slope, intercept, resid


results = []
plot_data = []  # (label, per_model, slope, intercept, resid) for the residual-plot grid

for region in MAIN_REGIONS:
    sub = df_all[df_all.region == region]
    out = check_case(sub, region)
    if out:
        res, per_model, slope, intercept, resid = out
        results.append(res)
        plot_data.append((region, per_model, slope, intercept, resid))

for region in GROUPED_REGIONS:
    sub_region = df_all[df_all.region == region]
    # Group 1 excluded here: it is not used in the linear regression / Fig. 6 analysis
    # (no fit line is shown for Group 1 in Fig. S ms_vs_bmb_other_sectors_bygroup), so
    # only the pooled (ungrouped) and Group 2 cases are of interest for this table.
    for gname, gsub in [('pooled', sub_region), ('group2', sub_region[sub_region.group == 2])]:
        out = check_case(gsub, f"{region}_{gname}")
        if out:
            res, per_model, slope, intercept, resid = out
            results.append(res)
            plot_data.append((f"{region} ({gname})", per_model, slope, intercept, resid))

results_df = pd.DataFrame(results)
pd.set_option('display.width', 160)
pd.set_option('display.float_format', lambda x: f'{x:.3g}')
print(results_df.to_string(index=False))

os.makedirs(OUTTABLES, exist_ok=True)
results_df.to_csv(os.path.join(OUTTABLES, "Checking_Statistics_fig6.csv"), index=False)
print("\nsaved", os.path.join(OUTTABLES, "Checking_Statistics_fig6.csv"))


# --- formatted LaTeX SI table, purely quantitative: gap, BP p, cluster p, with a
# pass/fail flag for each of the three assumptions (thresholds: |gap|<=0.1 for
# linearity, BP p>=0.05 for homoscedasticity, cluster p<0.05 for significance) ---
def fmt_p(p):
    """Format a p-value for LaTeX: plain decimal if >=1e-3, else proper
    scientific notation in math mode (e.g. $3.6\\times10^{-15}$) rather than
    a bare 'e-15' string, which LaTeX does not typeset well."""
    if p >= 1e-3:
        return f"{p:.2g}"
    mant, exp = f"{p:.1e}".split("e")
    return f"${float(mant):.1f}\\times10^{{{int(exp)}}}$"


def flag(ok):
    return "OK" if ok else "flag"


rows_tex = []
for _, r in results_df.iterrows():
    lin_ok = abs(r['gap']) <= 0.1
    homo_ok = r['bp_p'] >= 0.05
    sig = r['p_cluster'] < 0.05
    case_label = r['case'].replace('_', ' (') + ')' if '_' in r['case'] else r['case']
    rows_tex.append(
        f"{case_label} & {r['n_models']:.0f} & {r['r']:.2f} & {r['rho']:.2f} & "
        f"{r['gap']:.2f} ({flag(lin_ok)}) & {fmt_p(r['bp_p'])} ({flag(homo_ok)}) & "
        f"{fmt_p(r['p_naive'])} & {fmt_p(r['p_cluster'])} ({'sig.' if sig else 'n.s.'}) \\\\"
    )

tex = (
    "% ---- Auto-generated by Checking_Statistics_fig6.py. Do not edit by hand; re-run the script. ----\n"
    "\\begin{table}[h]\n\\centering\n\\small\n"
    "\\begin{tabular}{lccccccc}\n\\hline\n"
    "Case & $n$ models & $r$ & $\\rho$ & $\\rho-r$ (linearity) & BP $p$ (homoscedasticity) & "
    "naive $p$ & cluster-robust $p$ (independence) \\\\\n\\hline\n"
    + "\n".join(rows_tex) +
    "\n\\hline\n\\end{tabular}\n"
    "\\caption{Quantitative assumption checks for all 10 region/calving-group combinations "
    "underlying Fig.~\\ref{figure_6} (Amundsen, Wilkes) and Fig.~S\\ref{figure_SI_ms_vs_bmb_other_sectors_bygroup} "
    "(AIS, Ross, Filchner-Ronne, Aurora $\\times$ pooled/Group~2). Linearity: gap between Spearman "
    "$\\rho$ and Pearson $r$ on per-model-aggregated data; flagged if $|\\rho-r|>0.1$. Homoscedasticity: "
    "Breusch-Pagan test $p$-value; flagged if $p<0.05$. Independence: naive OLS $p$-value vs.\\ cluster-robust "
    "$p$-value (clustered by model), since each model contributes up to 4 pooled-experiment points; "
    "significance is reported using the cluster-robust value. Pooled (ungrouped) rows are shown for "
    "reference only; the substantive result in each sector is the Group~2 row (Group~1 is excluded from "
    "this regression throughout, see main text).}\n"
    "\\label{table_checking_statistics_fig6}\n\\end{table}\n"
)

tex_path = os.path.join(OUTTABLES, "Table_SI_checking_statistics_fig6.tex")
with open(tex_path, "w") as f:
    f.write(tex)
print("saved", tex_path)

# --- residual-plot grid, one panel per case ---
n = len(plot_data)
ncols = 4
nrows = -(-n // ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows), squeeze=False)
for idx, (label, per_model, slope, intercept, resid) in enumerate(plot_data):
    ax = axes[idx // ncols, idx % ncols]
    ax.axhline(0, color='gray', lw=1)
    ax.scatter(per_model.ms, resid, s=14, color='firebrick')
    ax.set_title(label, fontsize=9)
    ax.set_xlabel("MS", fontsize=8)
    ax.set_ylabel("residual", fontsize=8)
    ax.tick_params(labelsize=7)
for idx in range(n, nrows * ncols):
    axes[idx // ncols, idx % ncols].axis('off')
fig.suptitle("Figure 6 / Figure S ms_vs_bmb_bygroup: residuals vs MS, all cases", fontsize=11, y=1.01)
plt.tight_layout()
os.makedirs(FIGDIR, exist_ok=True)
outfig = os.path.join(FIGDIR, "Checking_Statistics_fig6_residuals.pdf")
plt.savefig(outfig, dpi=200, bbox_inches='tight')
print("saved", outfig)
