import os
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan

# ---------------------------------------------------------------------------
# Checking_Statistics_fig_MS_DSLR.py
#
# Quantitative assumption checks for the MS factor vs. dSLR (2300) regressions
# behind the SI figure "Figs SI Analyse dyn SLR  cum BMB-region of interest
# -latest.ipynb" (regression (2) in the consolidated Methods paragraph: DSLR
# directly against the MS factor across models, as an alternative predictor
# to cumulative BMB). Same design as Fig. 8 and the grounding-line-corridor
# regression, but with x = MS (melt_sensetivity_{region} from
# best_ms_method_per_model.csv, one value per model per region, identical
# across the four experiments) rather than a cumulative BMB integral.
#
# One regression per (region, experiment) pair, not pooled across
# experiments (the four per-experiment slopes are averaged afterward to get
# one representative slope per region, per Methods), so independence across
# experiments holds by construction.
#
# For every case this computes:
#   (i)   linearity        - Spearman rho vs Pearson r; a gap > 0.1 flags
#                             curvature the straight-line fit misses.
#   (ii)  homoscedasticity  - Breusch-Pagan test; p<0.05 flags residual
#                             variance that changes with MS.
#   (iii) independence      - satisfied by construction across experiments;
#                             supplementary check clustering by model family
#                             (perturbed-parameter variants of the same
#                             underlying code), as for Fig. 8 and the
#                             grounding-line-corridor regression.
#
# Saves reviewing/tables/Checking_Statistics_fig_MS_DSLR.csv
# and reviewing/tables/Table_SI_checking_statistics_fig_MS_DSLR.tex
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TABLES = os.path.join(PROJECT_ROOT, "analysis", "tables")
OUTTABLES = os.path.join(PROJECT_ROOT, "reviewing", "tables")

EXPS = ['expAE02', 'expAE03', 'expAE04', 'expAE05']
REGIONS = ['AIS', 'Amundsen', 'Ross', 'FilchnerRonne', 'Aurora', 'Wilkes']
REGION_SECTORS = {
    'AIS': None,  # use the precomputed AIS column directly
    'Amundsen': ['sector_4'],
    'Ross': ['sector_2', 'sector_6'],
    'FilchnerRonne': ['sector_5', 'sector_11'],
    'Aurora': ['sector_8'],
    'Wilkes': ['sector_7'],
}

FAMILY = {
    'DC_ISSM': 'DC_ISSM',
    'DOE_MALI_4km': 'DOE_MALI', 'DOE_MALI_8km_Ant95': 'DOE_MALI', 'DOE_MALI_8km_AntMean': 'DOE_MALI',
    'IGE_ElmerIce': 'IGE_ElmerIce',
    'ILTS_SICOPOLIS': 'ILTS_SICOPOLIS',
    'IMAU_UFEMISM1': 'IMAU_UFEMISM', 'IMAU_UFEMISM2': 'IMAU_UFEMISM',
    'IMAU_UFEMISM3': 'IMAU_UFEMISM', 'IMAU_UFEMISM4': 'IMAU_UFEMISM',
    'LSCE_GRISLI': 'LSCE_GRISLI', 'LSCE_GRISLI2': 'LSCE_GRISLI',
    'NCAR_CISM1': 'NCAR_CISM', 'NCAR_CISM2': 'NCAR_CISM',
    'NORCE_CISM2-MAR364-ERA-t1': 'NORCE_CISM', 'NORCE_CISM3-MAR364-ERA-t1': 'NORCE_CISM',
    'NORCE_CISM3-MAR364-ERA-t1-local': 'NORCE_CISM', 'NORCE_CISM3-MAR364-ERA-t1-nonlocal': 'NORCE_CISM',
    'NORCE_CISM4-MAR364-ERA-t1': 'NORCE_CISM', 'NORCE_CISM4-MAR364-ERA-t1-local': 'NORCE_CISM',
    'NORCE_CISM4-MAR364-ERA-t1-nonlocal': 'NORCE_CISM', 'NORCE_CISM4-MAR364-JRA-t1': 'NORCE_CISM',
    'NORCE_CISM5-MAR364-ERA-t1': 'NORCE_CISM', 'NORCE_CISM5-MAR364-ERA-t1-local': 'NORCE_CISM',
    'NORCE_CISM5-MAR364-ERA-t1-nonlocal': 'NORCE_CISM',
    'PIK_PISM': 'PIK_PISM',
    'UCM_Yelmo': 'UCM_Yelmo',
    'UCSD_ISSM': 'UCSD_ISSM',
    'ULB_fETISh-KoriBU1': 'ULB_fETISh-KoriBU', 'ULB_fETISh-KoriBU2': 'ULB_fETISh-KoriBU',
    'UNN_Ua': 'UNN_Ua',
    'UTAS_ElmerIce': 'UTAS_ElmerIce',
    'VUB_AISMPALEO': 'VUB_AISMPALEO',
    'VUW_PISM1': 'VUW_PISM', 'VUW_PISM1_s1': 'VUW_PISM', 'VUW_PISM1_s2': 'VUW_PISM',
    'VUW_PISM1_s3': 'VUW_PISM', 'VUW_PISM1_s4': 'VUW_PISM',
    'VUW_PISM2': 'VUW_PISM', 'VUW_PISM2_s1': 'VUW_PISM', 'VUW_PISM2_s2': 'VUW_PISM',
    'VUW_PISM2_s3': 'VUW_PISM', 'VUW_PISM2_s4': 'VUW_PISM',
}

df_best = pd.read_csv(os.path.join(TABLES, "best_ms_method_per_model.csv"))
models_in_table = set(df_best['Model'])

df_dslr = pd.read_csv(os.path.join(TABLES, "dslc_anom_2300.csv"))
df_bmb = pd.read_csv(os.path.join(TABLES, "cumulative_shelfmelt_2300.csv"))  # for model list per experiment only


def region_value_dslr(df, region):
    if region == 'AIS':
        return df['AIS']
    cols = REGION_SECTORS[region]
    return sum(df[c] for c in cols)


def check_case(sub, label):
    if len(sub) < 5:
        return None
    slope, intercept, r, p_lin, se = linregress(sub.x, sub.y)
    rho, p_rho = spearmanr(sub.x, sub.y)
    resid = sub.y - (slope * sub.x + intercept)
    try:
        bp_stat, bp_p, _, _ = het_breuschpagan(resid, sm.add_constant(sub['x']))
    except Exception:
        bp_p = np.nan

    X = sm.add_constant(sub['x'])
    clu = sm.OLS(sub['y'], X).fit(cov_type='cluster', cov_kwds={'groups': sub['family']})
    p_family_cluster = clu.pvalues['x']

    return dict(case=label, n=len(sub), n_families=sub.family.nunique(), r=r, r2=r**2, p_lin=p_lin,
                rho=rho, gap=rho - r, bp_p=bp_p, slope=slope, p_family_cluster=p_family_cluster)


results = []
for region in REGIONS:
    ms_col = f"melt_sensetivity_{region}"
    for exp in EXPS:
        dslr_exp = df_dslr[df_dslr.Experiment == exp].reset_index(drop=True)
        bmb_exp = df_bmb[df_bmb.Experiment == exp].reset_index(drop=True)  # gives the model list for this experiment
        x = df_best[['Model', ms_col]].rename(columns={ms_col: 'x'})
        y = pd.DataFrame({'Model': dslr_exp.Model, 'y': region_value_dslr(dslr_exp, region)})
        # restrict to models actually present in this experiment (matches source notebook)
        models_this_exp = set(bmb_exp.Model)
        merged = pd.merge(x, y, on='Model')
        merged = merged[merged.Model.isin(models_in_table) & merged.Model.isin(models_this_exp)]
        if exp == 'expAE02':
            merged = merged[merged.Model != 'VUW_PISM1_s1']
        merged = merged.dropna()
        merged['family'] = merged.Model.map(FAMILY)
        out = check_case(merged, f"{region}_{exp}")
        if out:
            results.append(out)

df_res = pd.DataFrame(results)
pd.set_option('display.width', 160)
pd.set_option('display.float_format', lambda x: f'{x:.3g}')
print(df_res.to_string(index=False))

os.makedirs(OUTTABLES, exist_ok=True)
df_res.to_csv(os.path.join(OUTTABLES, "Checking_Statistics_fig_MS_DSLR.csv"), index=False)
print("\nsaved", os.path.join(OUTTABLES, "Checking_Statistics_fig_MS_DSLR.csv"))

print("\n--- summary ---")
print("n cases:", len(df_res))
print("linearity flags (|gap|>0.1):", (df_res.gap.abs() > 0.1).sum())
print("homoscedasticity flags (bp_p<0.05):", (df_res.bp_p < 0.05).sum())
print("significant, naive p (p_lin<0.05):", (df_res.p_lin < 0.05).sum())
print("significant, family-cluster-robust p (<0.05):", (df_res.p_family_cluster < 0.05).sum())
print("max family-cluster-robust p:", df_res.p_family_cluster.max())


# --- formatted LaTeX SI table ---
def fmt_p(p):
    if p >= 1e-3:
        return f"{p:.2g}"
    mant, exp = f"{p:.1e}".split("e")
    return f"${float(mant):.1f}\\times10^{{{int(exp)}}}$"


def flag(ok):
    return "OK" if ok else "flag"


rows_tex = []
for _, r in df_res.iterrows():
    region, exp = r['case'].rsplit('_', 1)
    lin_ok = abs(r['gap']) <= 0.1
    homo_ok = r['bp_p'] >= 0.05
    rows_tex.append(
        f"{region} & {exp} & {r['n']:.0f} & {r['r']:.2f} & {r['rho']:.2f} & "
        f"{r['gap']:.2f} ({flag(lin_ok)}) & {fmt_p(r['bp_p'])} ({flag(homo_ok)}) & "
        f"{fmt_p(r['p_lin'])} & {fmt_p(r['p_family_cluster'])} \\\\"
    )

tex = (
    "% ---- Auto-generated by Checking_Statistics_fig_MS_DSLR.py. Do not edit by hand; re-run the script. ----\n"
    "\\begin{table}[h]\n\\centering\n\\small\n"
    "\\begin{tabular}{llccccccc}\n\\hline\n"
    "Region & Experiment & $n$ models & $r$ & $\\rho$ & $\\rho-r$ (linearity) & "
    "BP $p$ (homoscedasticity) & naive $p$ & family-cluster-robust $p$ \\\\\n\\hline\n"
    + "\n".join(rows_tex) +
    "\n\\hline\n\\end{tabular}\n"
    "\\caption{Quantitative assumption checks for all 24 region/experiment combinations "
    "underlying Fig.~S\\ref{figure_SI_dynslr_MS} (melt sensitivity (MS) vs.\\ dynamic "
    "sea-level-rise contribution, 2300). Each regression uses one point per model within a "
    "single experiment (not pooled across experiments; the four per-experiment slopes are "
    "averaged afterward to obtain one representative slope per region), so residual "
    "independence across experiments holds by construction. As a supplementary check on "
    "independence across models, since several models are perturbed-parameter variants of "
    "the same underlying code (e.g.\\ 10 VUW\\_PISM and 11 NORCE\\_CISM variants), we also "
    "report a $p$-value clustering the 43 models into 16 model families. Linearity: gap "
    "between Spearman $\\rho$ and Pearson $r$; flagged if $|\\rho-r|>0.1$. Homoscedasticity: "
    "Breusch-Pagan test $p$-value; flagged if $p<0.05$.}\n"
    "\\label{table_checking_statistics_fig_MS_DSLR}\n\\end{table}\n"
)

tex_path = os.path.join(OUTTABLES, "Table_SI_checking_statistics_fig_MS_DSLR.tex")
with open(tex_path, "w") as f:
    f.write(tex)
print("saved", tex_path)
