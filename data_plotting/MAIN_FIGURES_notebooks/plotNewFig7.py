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
REGIONS = ['Amundsen', 'Wilkes']
REGION_KEYS = {'Amundsen': ['_sector_4'], 'Wilkes': ['_sector_7']}
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
            all_recs.append(dict(region=region, experiment=exp, model=model, ms=ms[0], bmb_change=cumulative_change))

df_all = pd.DataFrame(all_recs)

fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))
for ax, region in zip(axes, REGIONS):
    sub = df_all[df_all.region == region]
    ax.scatter(sub.ms, sub.bmb_change, color='black', s=12, alpha=0.5)
    slope, intercept, r, p_lin, se = linregress(sub.ms, sub.bmb_change)
    rho, p_spear = spearmanr(sub.ms, sub.bmb_change)

    # cluster-robust p-value (clustered by model): this regression pools all 4 experiments,
    # so each model contributes up to 4 non-independent points sharing the same MS value;
    # naive OLS p-values above overstate significance (see Methods, regression assumptions)
    X = sm.add_constant(sub['ms'])
    clu = sm.OLS(sub['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['model']})
    p_cluster = clu.pvalues['ms']

    xline = np.linspace(sub.ms.min(), sub.ms.max(), 50)
    ax.plot(xline, slope * xline + intercept, '--', color='gray', lw=1.5)
    # extra headroom above the data so the annotation text doesn't sit on top of any points
    ymin, ymax = sub.bmb_change.min(), sub.bmb_change.max()
    pad = 0.4 * (ymax - ymin)
    ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + pad)
    ax.text(0.03, 0.97,
            f"lin $R^2$={r**2:.2f}, $p$={p_cluster:.2g} (cluster-robust)\n"
            f"Spearman ρ={rho:.2f}, $p$={p_spear:.2g}",
            transform=ax.transAxes, fontsize=8, va='top')
    ax.set_title(region, fontsize=11)
    ax.set_xlabel("Melt sensitivity (m yr⁻¹ K⁻¹)")
    ax.set_ylabel("Cumulative ΔBMB\n2015-2300 (Gt)")

plt.tight_layout()
plt.savefig("fig6_pooled.pdf", dpi=150)
print("saved fig6_pooled.png")
