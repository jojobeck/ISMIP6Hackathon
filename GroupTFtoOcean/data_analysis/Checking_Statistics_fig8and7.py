import os
import numpy as np
import pandas as pd
import xarray as xr
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
            all_recs.append(dict(region=region, experiment=exp, model=model, ms=ms[0],
                                  bmb_change=cumulative_change, is_vuw=model.startswith("VUW_PISM")))

df_all = pd.DataFrame(all_recs)

rows = []
for region in REGIONS:
    sub_full = df_all[df_all.region == region]
    for label, sub in [('all models', sub_full), ('VUW_PISM excluded', sub_full[~sub_full.is_vuw])]:
        slope, intercept, r, p_lin, se = linregress(sub.ms, sub.bmb_change)
        rho, p_spear = spearmanr(sub.ms, sub.bmb_change)
        X = sm.add_constant(sub['ms'])
        clu = sm.OLS(sub['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['model']})
        p_cluster = clu.pvalues['ms']
        rows.append(dict(region=region, subset=label, n_pts=len(sub), n_models=sub.model.nunique(),
                          slope=slope, r2=r**2, pearson_p=p_cluster,
                          spearman_rho=rho, spearman_p=p_spear))

out = pd.DataFrame(rows)
pd.set_option('display.width', 150)
print(out.to_string(index=False))

out_csv = os.path.join(PROJECT_ROOT, "reviewing", "tables", "Table_fig6_pooled_vuw_robustness.csv")
out.to_csv(out_csv, index=False)
print("saved", out_csv)
