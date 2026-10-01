import os, time
import numpy as np
import xarray as xr
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT = os.path.join(PROJECT_ROOT, "ComputedScalars")
EXPS = ['expAE02', 'expAE03', 'expAE04', 'expAE05']

# region -> list of sector suffixes to sum (per Johanna's mapping / ANOVA_bmr_regional.m)
REGION_KEYS = {
    'AIS':           [''],
    'Amundsen':      ['_sector_4'],
    'Ross':          ['_sector_2', '_sector_6'],
    'FilchnerRonne': ['_sector_5', '_sector_11'],
    'Aurora':        ['_sector_8'],
    'Wilkes':        ['_sector_7'],
}

# calving groups per Johanna's plotting notebook (models_group1/2, group4 folded into group2)
models_group1 = ['VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3', 'VUW_PISM1_s4', 'VUW_PISM2',
                 'VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2',
                 'LSCE_GRISLI', 'UCM_Yelmo', 'IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4']
models_group2 = ['DC_ISSM', 'ILTS_SICOPOLIS', 'VUB_AISMPALEO', 'NCAR_CISM1', 'NORCE_CISM3-MAR364-ERA-t1-nonlocal',
                 'NORCE_CISM4-MAR364-ERA-t1-nonlocal', 'NORCE_CISM5-MAR364-ERA-t1-nonlocal', 'NCAR_CISM2',
                 'NORCE_CISM3-MAR364-ERA-t1-local', 'NORCE_CISM4-MAR364-ERA-t1-local',
                 'NORCE_CISM5-MAR364-ERA-t1-local', 'UNN_Ua', 'NORCE_CISM2-MAR364-ERA-t1',
                 'NORCE_CISM3-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-JRA-t1',
                 'NORCE_CISM5-MAR364-ERA-t1', 'UTAS_ElmerIce', 'ULB_fETISh-KoriBU2', 'IGE_ElmerIce',
                 'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean']
models_group4 = ['ULB_fETISh-KoriBU1', 'UCSD_ISSM']  # folded into group 2

GROUP = {m: 1 for m in models_group1}
GROUP.update({m: 2 for m in models_group2 + models_group4})

BMB_FACTOR = -1 * 365.25 * 24 * 60 * 60 / 10**12  # matches Johanna's notebook: -> Gt/yr


def load_curves(variable, exp):
    """variable: 'iareafl' or 'shelfmelt'. Returns dict region -> (curves[n_models,n_time], groups[n_models], time)"""
    ddir = os.path.join(DATA_ROOT, exp, variable)
    files = sorted(f for f in os.listdir(ddir) if f.endswith(".nc"))
    per_model = {}
    time_ref = None
    for f in files:
        fname_model = f.replace(f"computed_{variable}_AIS_", "").replace(f"_{exp}.nc", "")
        if fname_model not in GROUP:
            continue
        d = xr.open_dataset(os.path.join(ddir, f), decode_times=False)
        t = d["time"].values
        region_vals = {}
        for region, keys in REGION_KEYS.items():
            total = sum(d[f"{variable}{k}"].values.astype(float) for k in keys)
            region_vals[region] = total
        d.close()
        per_model[fname_model] = (region_vals, t)
        if time_ref is None or len(t) < len(time_ref):
            time_ref = t

    n_time = len(time_ref)
    out = {}
    for region in REGION_KEYS:
        curves, groups, names = [], [], []
        for model, (region_vals, t) in per_model.items():
            curves.append(region_vals[region][:n_time])
            groups.append(GROUP[model])
            names.append(model)
        out[region] = (np.array(curves), np.array(groups), time_ref[:n_time], names)
    return out


def fanova_stat(curves, groups):
    a = curves[groups == 1]
    b = curves[groups == 2]
    mean_a, mean_b = a.mean(0), b.mean(0) #mean along time 
    na, nb = a.shape[0], b.shape[0] #number of models 
    #ddof devides by e.g. na-1 not na for standaed sample variance , e degree of freedom already used for mean calculation
    #then multiply ueach goup variance by number of smaples fividing by total sample size _> weighting it 
    pooled_var = (a.var(0, ddof=1) * (na - 1) + b.var(0, ddof=1) * (nb - 1)) / (na + nb - 2)
    pooled_var = np.maximum(pooled_var, 1e-12) #Clamping it to a small positive floor 
    stat_curve = (mean_a - mean_b) ** 2 / pooled_var #normalized squared distance 
    return np.trapz(stat_curve, dx=1.0) #sum (integral) one number


def permutation_test(curves, groups, n_perm=5000, seed=42):
    obs = fanova_stat(curves, groups)
    rng = np.random.default_rng(seed)
    perm_stats = np.empty(n_perm)
    for p in range(n_perm):
        perm_stats[p] = fanova_stat(curves, rng.permutation(groups))#permutation of group labeling 1 &2 
    pval = (np.sum(perm_stats >= obs) + 1) / (n_perm + 1)
    return obs, pval, perm_stats


def run(variable, label):
    rows = []
    for exp in EXPS:
        t0 = time.time()
        data = load_curves(variable, exp)
        for region in REGION_KEYS:
            curves, groups, t, names = data[region]
            # normalize as % change from first time step for iareafl; keep physical units (scaled) for shelfmelt
            if variable == "iareafl":
                base = curves[:, [0]]
                base = np.where(base == 0, np.nan, base)
                proc = 100.0 * (curves - base) / base
            else:
                proc = curves * BMB_FACTOR
            valid = ~np.isnan(proc).any(axis=1)
            proc_v, groups_v = proc[valid], groups[valid]
            n1 = int((groups_v == 1).sum())
            n2 = int((groups_v == 2).sum())
            obs, pval, _ = permutation_test(proc_v, groups_v, n_perm=5000)
            rows.append(dict(variable=label, region=region, experiment=exp,
                              n_group1=n1, n_group2=n2, stat=obs, p_value=pval))
        print(f"{variable} {exp} done in {time.time()-t0:.1f}s")
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df_area = run("iareafl", "shelf area (% chg)")
    df_bmb = run("shelfmelt", "BMB (Gt/yr)")
    df = pd.concat([df_area, df_bmb], ignore_index=True)
    df.to_csv(os.path.join(PROJECT_ROOT, "reviewing", "tables", "fanova_all_regions_results.csv"), index=False)
    print(df.to_string(index=False))
