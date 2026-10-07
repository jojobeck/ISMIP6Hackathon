import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import statsmodels.api as sm

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT = os.path.join(PROJECT_ROOT, "ComputedScalars")
TABLES = os.path.join(PROJECT_ROOT, "analysis", "tables")
EXPS = ['expAE02', 'expAE03', 'expAE04', 'expAE05']
REGIONS = ['AIS', 'Ross', 'FilchnerRonne', 'Aurora']
REGION_KEYS = {
    'AIS': [''], 'Ross': ['_sector_2', '_sector_6'],
    'FilchnerRonne': ['_sector_5', '_sector_11'], 'Aurora': ['_sector_8'],
}
FACTOR = -1 * 365.25 * 24 * 60 * 60 / 10**12

df_best = pd.read_csv(os.path.join(TABLES, "best_ms_method_per_model.csv"))
models_in_table = set(df_best['Model'])

all_recs = []
for exp in EXPS:
    ddir = os.path.join(DATA_ROOT, exp, "shelfmelt")
    files = sorted(f for f in os.listdir(ddir) if f.endswith(".nc"))
    for region in REGIONS:
        ms_col = f"melt_sensetivity_{region}"
        for f in files:
            model = f.replace("computed_shelfmelt_AIS_", "").replace(f"_{exp}.nc", "")
            if model not in models_in_table:
                continue
            ms = df_best.loc[df_best['Model'] == model, ms_col].values
            if len(ms) == 0 or np.isnan(ms[0]):
                continue
            d = xr.open_dataset(os.path.join(ddir, f), decode_times=False)
            y = sum(d[f"shelfmelt{k}"].values for k in REGION_KEYS[region]) * FACTOR
            time = d["time"].values
            d.close()
            cumulative_change = np.trapz(y - y[0], time)
            all_recs.append(dict(region=region, experiment=exp, model=model, ms=ms[0], bmb_change=cumulative_change,
                                  is_vuw=model.startswith("VUW_PISM")))

df_all = pd.DataFrame(all_recs)
VUW_COLOR = '#888888'

# --- Table: per-experiment stats + a "pooled (all exps)" row per sector, with and without VUW_PISM ---
# p-value reported is cluster-robust (clustered by model) for the pooled (all 4 exps) rows, since each
# model contributes up to 4 non-independent points there; per-experiment rows are one point per model
# already, so the naive p-value is valid as-is.
def stat_row(sub, region, exp_label, vuw_status, cluster=False):
    slope, intercept, r, p_lin, se = linregress(sub.ms, sub.bmb_change)
    rho, p_spear = spearmanr(sub.ms, sub.bmb_change)
    p_report = p_lin
    if cluster and sub.model.nunique() >= 3:
        X = sm.add_constant(sub['ms'])
        p_report = sm.OLS(sub['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['model']}).pvalues['ms']
    return dict(region=region, experiment=exp_label, vuw_pism=vuw_status, n=len(sub),
                pearson_r2=r**2, pearson_p=p_report, spearman_rho=rho, spearman_p=p_spear)

rows = []
for region in REGIONS:
    for exp in EXPS:
        sub = df_all[(df_all.region == region) & (df_all.experiment == exp)]
        rows.append(stat_row(sub, region, exp, 'included'))
        rows.append(stat_row(sub[~sub.is_vuw], region, exp, 'excluded'))
    sub = df_all[df_all.region == region]
    rows.append(stat_row(sub, region, 'pooled (all 4 exps)', 'included', cluster=True))
    rows.append(stat_row(sub[~sub.is_vuw], region, 'pooled (all 4 exps)', 'excluded', cluster=True))

results = pd.DataFrame(rows)
results.to_csv("Table_SI_ms_vs_bmb_other_sectors.csv", index=False)
pd.set_option('display.width', 140)
print(results.to_string(index=False))

# --- Figure: pooled scatter, one panel per sector ---
fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5))
for ax, region in zip(axes.flat, REGIONS):
    sub = df_all[df_all.region == region]
    non_vuw = sub[~sub.is_vuw]
    vuw = sub[sub.is_vuw]

    # all points, colored to flag VUW_PISM, but the fit/stats use the FULL sample (all points)
    ax.scatter(non_vuw.ms, non_vuw.bmb_change, color='black', s=12, alpha=0.5, label='other models')
    ax.scatter(vuw.ms, vuw.bmb_change, color=VUW_COLOR, marker='s', s=16, alpha=0.7,
               edgecolors='none', label='VUW_PISM')

    # fit including all points (black dashed) vs. fit excluding VUW_PISM (blue dashed)
    slope, intercept, r, p_lin, se = linregress(sub.ms, sub.bmb_change)
    rho, p_spear = spearmanr(sub.ms, sub.bmb_change)
    slope_x, intercept_x, r_x, p_lin_x, se_x = linregress(non_vuw.ms, non_vuw.bmb_change)
    rho_x, p_spear_x = spearmanr(non_vuw.ms, non_vuw.bmb_change)

    # cluster-robust p (clustered by model): both fits pool all 4 experiments, so each
    # model contributes up to 4 non-independent points sharing one MS value
    X = sm.add_constant(sub['ms'])
    p_cluster = sm.OLS(sub['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['model']}).pvalues['ms']
    X_x = sm.add_constant(non_vuw['ms'])
    p_cluster_x = sm.OLS(non_vuw['bmb_change'], X_x).fit(cov_type='cluster', cov_kwds={'groups': non_vuw['model']}).pvalues['ms']

    xline = np.linspace(sub.ms.min(), sub.ms.max(), 50)
    ax.plot(xline, slope * xline + intercept, '--', color='gray', lw=1.5, label='fit, all points')
    xline_x = np.linspace(non_vuw.ms.min(), non_vuw.ms.max(), 50)
    ax.plot(xline_x, slope_x * xline_x + intercept_x, '--', color='tab:blue', lw=1.5, label='fit, VUW excl.')

    ymin, ymax = sub.bmb_change.min(), sub.bmb_change.max()
    pad = 0.45 * (ymax - ymin)
    ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + pad)
    ax.text(0.03, 0.97,
            f"All: R²={r**2:.2f}, p={p_cluster:.2g} (cluster-robust), ρ={rho:.2f}, p={p_spear:.2g}\n"
            f"VUW excl.: R²={r_x**2:.2f}, p={p_cluster_x:.2g} (cluster-robust), ρ={rho_x:.2f}, p={p_spear_x:.2g}",
            transform=ax.transAxes, fontsize=7.5, va='top')
    ax.set_title(region, fontsize=11)
    ax.set_xlabel("Melt sensitivity (m yr⁻¹ K⁻¹)")
    ax.set_ylabel("Cumulative ΔBMB\n2015-2300 (Gt)")

handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.02))
plt.tight_layout()
plt.savefig("Figure_SI_ms_vs_bmb_other_sectors.pdf", dpi=300, bbox_inches='tight')
print("saved Figure_SI_ms_vs_bmb_other_sectors.png")

caption = """Figure Sxx. Melt sensitivity factor versus cumulative basal mass balance (BMB) change (2015-2300) for AIS, Ross, Filchner-Ronne, and Aurora, the four sectors not shown in Figure 6 (all four climate-forcing experiments, expAE02-05, pooled per sector; n=172 points, 43 models per panel). Gray squares mark the VUW_PISM ensemble (9 members, one model family with multiple parameter perturbations); all other models are shown as black circles. Gray dashed lines show ordinary least-squares fits to the full sample (all points); blue dashed lines show fits with the VUW_PISM ensemble excluded, to test whether any apparent relationship is robust to the disproportionate weight of a single model family in the pooled sample. Annotated statistics give the Pearson R^2 and p-value and the Spearman rank correlation (rho) and p-value for both fits; reported p-values are cluster-robust (clustered by model), since each model contributes up to four points sharing a single MS value across the pooled experiments. In AIS and Ross, melt sensitivity is unrelated to cumulative basal melt regardless of whether VUW_PISM is included (R^2 <= 0.09 throughout), consistent with calving-driven shelf area change dominating cumulative basal melt in these sectors. Filchner-Ronne and Aurora show an apparent relationship when VUW_PISM is included (positive and negative, respectively), but this relationship disappears once VUW_PISM is excluded (Filchner-Ronne: R^2 drops from 0.20 to 0.05; Aurora: R^2 drops from 0.25 to 0.003), showing that it was driven by the disproportionate influence of this one model family -- which combines high melt sensitivity with severe shelf-area collapse that caps total melt volume -- rather than reflecting a genuine, ensemble-wide physical relationship. Accounting for this, no sector besides Amundsen and Wilkes (Figure 6, robust with or without VUW_PISM: R^2 = 0.31 and 0.26 respectively when VUW_PISM is excluded) shows a melt-sensitivity-to-cumulative-BMB relationship that holds across the full model ensemble. See Table Sxx for full per-experiment and pooled statistics, with and without VUW_PISM."""

with open("Figure_SI_ms_vs_bmb_other_sectors_caption.txt", "w") as f:
    f.write(caption)
print("saved caption")
