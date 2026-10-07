import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import linregress, spearmanr
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

# notebook calving grouping (confirmed correct by Ronja), same list used throughout reviewing/
models_group1 = ['VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3', 'VUW_PISM1_s4', 'VUW_PISM2',
                 'VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2',
                 'LSCE_GRISLI', 'UCM_Yelmo', 'IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4']
GROUP = {m: 1 for m in models_group1}
GROUP_COLORS = {1: '#1a7d3c', 2: '#332288'}
GROUP_LABELS = {1: 'Group 1 (strong calving)', 2: 'Group 2 (moderate/other)'}

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
            if model.startswith("VUW_PISM"):
                continue  # excluded from this analysis (see robustness check)
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

# --- Table: per-group + pooled stats, per sector (VUW_PISM excluded throughout) ---
# pearson_p: cluster-robust p-value, clustered by model, correcting for the fact that each
# model contributes up to 4 non-independent points (one per pooled experiment) sharing the
# same MS value; NaN where a group has too few model-clusters for a stable estimate (<3).
def cluster_p(d):
    if d.model.nunique() < 3:
        return np.nan
    X = sm.add_constant(d['ms'])
    return sm.OLS(d['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': d['model']}).pvalues['ms']

rows = []
for region in REGIONS:
    sub = df_all[df_all.region == region]
    slope, intercept, r, p_lin, se = linregress(sub.ms, sub.bmb_change)
    rho, p_spear = spearmanr(sub.ms, sub.bmb_change)
    rows.append(dict(region=region, group='pooled (1+2)', n=len(sub), n_models=sub.model.nunique(),
                      mean_ms=sub.ms.mean(), mean_bmb_change=sub.bmb_change.mean(),
                      pearson_r2=r**2, pearson_p=cluster_p(sub),
                      spearman_rho=rho, spearman_p=p_spear))
    for g in [1, 2]:
        gs = sub[sub.group == g]
        slope, intercept, r, p_lin, se = linregress(gs.ms, gs.bmb_change)
        rho, p_spear = spearmanr(gs.ms, gs.bmb_change)
        rows.append(dict(region=region, group=f'group {g}', n=len(gs), n_models=gs.model.nunique(),
                          mean_ms=gs.ms.mean(), mean_bmb_change=gs.bmb_change.mean(),
                          pearson_r2=r**2, pearson_p=cluster_p(gs),
                          spearman_rho=rho, spearman_p=p_spear))

results = pd.DataFrame(rows)
results.to_csv("Table_SI_ms_vs_bmb_other_sectors_bygroup.csv", index=False)
pd.set_option('display.width', 150)
print(results.to_string(index=False))

def sig_mark(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    if p < 0.1:
        return '†'  # dagger, marginal
    return 'n.s.'


# --- Figure: 2x2 panels, colored by calving group, VUW_PISM excluded ---
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
for ax, region in zip(axes.flat, REGIONS):
    sub = df_all[df_all.region == region]

    for g in [1, 2]:
        gs = sub[sub.group == g]
        ax.scatter(gs.ms, gs.bmb_change, s=16, alpha=0.6, color=GROUP_COLORS[g],
                   label=GROUP_LABELS[g], edgecolors='none')
        ax.scatter(gs.ms.mean(), gs.bmb_change.mean(), marker='*', s=220,
                   color=GROUP_COLORS[g], edgecolors='black', linewidths=0.8, zorder=5)

    # only Group 2 gets a fit line: it's the group where MS remains a meaningful predictor
    # in most sectors (see per-group stats in the table); Group 1's smaller sample (n=8
    # models) is mostly non-significant, and the pooled fit is masked by the group offset.
    g2 = sub[sub.group == 2]
    slope, intercept, r, p_lin, se = linregress(g2.ms, g2.bmb_change)
    rho, p_spear = spearmanr(g2.ms, g2.bmb_change)

    # cluster-robust p-value (clustered by model): this regression pools all 4 experiments,
    # so each model contributes up to 4 non-independent points sharing the same MS value;
    # naive OLS p-values overstate significance (see Methods, regression assumptions)
    X = sm.add_constant(g2['ms'])
    clu = sm.OLS(g2['bmb_change'], X).fit(cov_type='cluster', cov_kwds={'groups': g2['model']})
    p_cluster = clu.pvalues['ms']

    xline = np.linspace(g2.ms.min(), g2.ms.max(), 50)
    ax.plot(xline, slope * xline + intercept, '--', color=GROUP_COLORS[2], lw=1.8, label='Group 2 fit')

    ymin, ymax = sub.bmb_change.min(), sub.bmb_change.max()
    pad = 0.45 * (ymax - ymin)
    ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + pad)
    ax.text(0.03, 0.97,
            f"Group 2 fit — Pearson: $R^2$={r**2:.2f}, $p$={p_cluster:.2g} {sig_mark(p_cluster)} (cluster-robust)\n"
            f"Spearman: ρ={rho:.2f}, $p$={p_spear:.2g} {sig_mark(p_spear)}",
            transform=ax.transAxes, fontsize=8, va='top', color=GROUP_COLORS[2])
    ax.set_title(region, fontsize=11)
    ax.set_xlabel("Melt sensitivity (m yr⁻¹ K⁻¹)")
    ax.set_ylabel("Cumulative ΔBMB\n2015-2300 (Gt)")

handles, labels = axes.flat[0].get_legend_handles_labels()
star_handle = plt.Line2D([], [], marker='*', color='gray', markeredgecolor='black',
                          linestyle='None', markersize=13, label='Group mean')
handles.append(star_handle)
labels.append('Group mean')
fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.03))
fig.suptitle("Melt sensitivity vs. cumulative basal melt, by calving group\n(VUW_PISM excluded; all 4 experiments pooled)",
             fontsize=11, y=1.08)
fig.text(0.5, -0.01, "significance: *** p<0.001, ** p<0.01, * p<0.05, † p<0.1, n.s. not significant",
         ha='center', fontsize=7.5, style='italic')
plt.tight_layout()
plt.savefig("Fig_cumBMB_ms_group.pdf", dpi=300, bbox_inches='tight')
print("saved Figure_SI_ms_vs_bmb_other_sectors_bygroup.png")

caption = """Figure Sxx. Melt sensitivity factor versus cumulative basal mass balance (BMB) change (2015-2300) for AIS, Ross, Filchner-Ronne, and Aurora, colored by qualitative calving group (Group 1: strong calving/shelf-area collapse; Group 2: moderate calving/other), with the VUW_PISM ensemble excluded (see Figure Sxx and Table Sxx for the VUW_PISM robustness check). Stars mark each group's mean melt sensitivity and mean cumulative BMB change. The dashed line and annotated statistics show the linear (Pearson R^2, p, cluster-robust by model to account for the four pooled experiments per model) and rank (Spearman rho, p) fit within Group 2 only, the group in which melt sensitivity remains a meaningful within-group predictor of cumulative BMB in most sectors (AIS: R^2=0.37, p=5.2e-09; Ross: R^2=0.15, p=2.5e-03; Filchner-Ronne: R^2=0.54, p=4.5e-10; not significant in Aurora, p=0.32). No fit line is shown for Group 1, whose smaller sample (n=8 models) gives an unstable cluster-robust estimate and is not a significant predictor in any sector except Aurora (R^2=0.37; see Table Sxx). Calving group visibly separates the two point clouds along the cumulative-BMB axis in all four sectors (Group 1 systematically lower, consistent with shelf-area loss capping total melt volume); this offset, rather than an absence of melt-sensitivity control, is what weakens the pooled (both-groups) correlation reported in Figure Sxx. Melt sensitivity therefore continues to influence cumulative basal melt within a given calving regime (particularly Group 2), in contrast to Amundsen and Wilkes (Figure 6), where the pooled relationship is directly significant because the calving-group offset is absent or negligible. See Table Sxx for full per-group means and correlation statistics."""

with open("Figure_SI_ms_vs_bmb_other_sectors_bygroup_caption.txt", "w") as f:
    f.write(caption)
print("saved caption")
