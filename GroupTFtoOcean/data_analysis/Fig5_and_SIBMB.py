import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT = os.path.join(PROJECT_ROOT, "ComputedScalars")
TABLES = os.path.join(PROJECT_ROOT, "analysis", "tables")
ALL_EXPS = ["expAE02", "expAE03", "expAE04", "expAE05"]
MAIN_EXP = "expAE04"
SI_EXPS = [e for e in ALL_EXPS if e != MAIN_EXP]
REGIONS = ['AIS', 'Amundsen', 'Ross', 'FilchnerRonne', 'Aurora', 'Wilkes']
REGION_KEYS = {
    'AIS': [''], 'Amundsen': ['_sector_4'], 'Ross': ['_sector_2', '_sector_6'],
    'FilchnerRonne': ['_sector_5', '_sector_11'], 'Aurora': ['_sector_8'], 'Wilkes': ['_sector_7'],
}
BMB_FACTOR = -1 * 365.25 * 24 * 60 * 60 / 10**12  # -> Gt/yr

models_group1 = ['VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3', 'VUW_PISM1_s4', 'VUW_PISM2',
                 'VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2',
                 'LSCE_GRISLI', 'UCM_Yelmo', 'IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4']
GROUP = {m: 1 for m in models_group1}
GROUP_COLORS = {1: '#117733', 2: '#88CCEE'}
GROUP_LABELS = {1: 'Group 1 (strong calving)', 2: 'Group 2 (moderate/other)'}

df_best = pd.read_csv(os.path.join(TABLES, "best_ms_method_per_model.csv"))
models_in_table = set(df_best['Model'])

fanova = pd.read_csv(os.path.join(PROJECT_ROOT, "reviewing", "tables", "fanova_all_regions_results.csv"))
fanova_bmb = fanova[fanova.variable == "BMB (Gt/yr)"].set_index(["region", "experiment"])
fanova_area = fanova[fanova.variable == "shelf area (% chg)"].set_index(["region", "experiment"])

HALO = [pe.Stroke(linewidth=4.5, foreground='white'), pe.Normal()]


def load_curve(exp, region, variable):
    """variable: 'area' or 'bmb'. Returns dict[model] -> (t, y), plus group assignment."""
    ddir = os.path.join(DATA_ROOT, exp, "shelfmelt" if variable == "bmb" else "iareafl")
    prefix = "computed_shelfmelt_AIS_" if variable == "bmb" else "computed_iareafl_AIS_"
    files = sorted(f for f in os.listdir(ddir) if f.endswith(".nc"))
    out = {}
    for f in files:
        model = f.replace(prefix, "").replace(f"_{exp}.nc", "")
        if model not in models_in_table:
            continue
        d = xr.open_dataset(os.path.join(ddir, f), decode_times=False)
        key = "shelfmelt" if variable == "bmb" else "iareafl"
        y = sum(d[f"{key}{k}"].values for k in REGION_KEYS[region])
        if variable == "bmb":
            y = y * BMB_FACTOR
        else:
            y = y / y[0] * 100.0
        t = d["time"].values
        d.close()
        out[model] = (t, y)
    return out


def plot_panel(ax, exp, region, variable, annotate=True, highlight=None):
    highlight = highlight or {}
    curves = load_curve(exp, region, variable)
    group_curves = {1: [], 2: []}
    group_times = {1: [], 2: []}
    for model, (t, y) in curves.items():
        g = GROUP.get(model, 2)
        group_curves[g].append(y)
        group_times[g].append(t)
        if model in highlight:
            continue  # drawn separately, on top of the spaghetti + group means
        col = GROUP_COLORS[g]
        ax.plot(t, y, color=col, lw=0.6, alpha=0.35)

    for g in [1, 2]:
        if not group_curves[g]:
            continue
        min_len = min(len(c) for c in group_curves[g])
        mean_y = np.mean([c[:min_len] for c in group_curves[g]], axis=0)
        t_mean = group_times[g][0][:min_len]
        ax.plot(t_mean, mean_y, color=GROUP_COLORS[g], lw=2.6,
                 solid_capstyle='round', zorder=5, path_effects=HALO)

    for model, col in highlight.items():
        if model in curves:
            t, y = curves[model]
            ax.plot(t, y, color=col, lw=1.8, alpha=0.95, zorder=6)

    if annotate:
        fdf = fanova_area if variable == "area" else fanova_bmb
        row = fdf.loc[(region, exp)]
        p_val = row["p_value"]
        is_sig = row["sig"] != "ns"
        label = f"fANOVA p={p_val:.3f} {'*' if is_sig else '(n.s.)'}"
        box_color = 'black' if is_sig else '#B00020'
        ax.text(0.97, 0.97, label, transform=ax.transAxes, fontsize=7,
                fontweight='bold' if is_sig else 'normal',
                ha='right', va='top', color=box_color, zorder=10,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=box_color, alpha=0.95, linewidth=0.7))


def legend_elements():
    els = [plt.Line2D([0], [0], color=GROUP_COLORS[g], lw=2, label=GROUP_LABELS[g]) for g in [1, 2]]
    els.append(plt.Line2D([0], [0], color='dimgray', lw=2.8, label='Group mean'))
    return els


# ---------------------------------------------------------------
# 1) Figure 4: shelf area, all 4 experiment columns (same layout as Figure 5)
# ---------------------------------------------------------------
fig, axes = plt.subplots(len(REGIONS), len(ALL_EXPS), figsize=(13, 11), sharex=True)
for i, region in enumerate(REGIONS):
    for j, exp in enumerate(ALL_EXPS):
        ax = axes[i, j]
        plot_panel(ax, exp, region, "area")
        ax.tick_params(axis='both', labelsize=11)
        if i == 0:
            ax.set_title(exp, fontsize=13)
        if j == 0:
            ax.set_ylabel(f"{region}\ndA (%)", fontsize=12)
        if i == len(REGIONS) - 1:
            ax.set_xlabel("Year", fontsize=12)
fig.legend(handles=legend_elements(), loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.03), fontsize=12, frameon=False)
plt.tight_layout()
plt.savefig("Figure4_shelfarea_allexps.pdf", dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close(fig)
print("saved Figure4_shelfarea_allexps.png")

# ---------------------------------------------------------------
# 3) Figure 5: BMB, all 4 experiment columns
# ---------------------------------------------------------------
fig, axes = plt.subplots(len(REGIONS), len(ALL_EXPS), figsize=(12, 10), sharex=True)
for i, region in enumerate(REGIONS):
    for j, exp in enumerate(ALL_EXPS):
        ax = axes[i, j]
        plot_panel(ax, exp, region, "bmb")
        if i == 0:
            ax.set_title(exp, fontsize=10)
        if j == 0:
            ax.set_ylabel(f"{region}\ndBMB (Gt yr⁻¹)", fontsize=9)
        if i == len(REGIONS) - 1:
            ax.set_xlabel("Year", fontsize=9)
fig.legend(handles=legend_elements(), loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.04), fontsize=9, frameon=False)
fig.suptitle("ΔBMB (Gt yr⁻¹), all experiments", fontsize=11, y=1.07)
plt.tight_layout()
plt.savefig("Figure5_BMB_allexps.pdf", dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close(fig)
print("saved Figure5_BMB_allexps.png")

# ---------------------------------------------------------------
# 4) SI Figure: shelf area, ALL 4 experiments (rows) x all 6 regions (cols),
#    same layout logic as Figure 5 (BMB, all experiments), but for shelf area,
#    with three individual models highlighted on top of the group spaghetti/means.
#    Only ULB_fETISh-KoriBU1 is highlighted (not KoriBU2), since it's the one
#    discussed in the text as showing (inconsistent) area growth; KoriBU2 is
#    left in the ordinary Group 2 spaghetti.
# ---------------------------------------------------------------
HIGHLIGHT_MODELS = {
    'UCSD_ISSM': '#FF69B4',            # pink
    'ULB_fETISh-KoriBU1': '#FF8C00',   # orange
    'UTAS_ElmerIce': '#E41A1C',        # red
}

fig, axes = plt.subplots(len(ALL_EXPS), len(REGIONS), figsize=(16, 9), sharex=True)
for i, exp in enumerate(ALL_EXPS):
    for j, region in enumerate(REGIONS):
        ax = axes[i, j]
        plot_panel(ax, exp, region, "area", annotate=False, highlight=HIGHLIGHT_MODELS)
        if i == 0:
            ax.set_title(region, fontsize=10)
        if j == 0:
            ax.set_ylabel(f"{exp}\ndA (%)", fontsize=9)
        if i == len(ALL_EXPS) - 1:
            ax.set_xlabel("Year", fontsize=9)

highlight_handles = []
seen_labels = set()
for model, col in HIGHLIGHT_MODELS.items():
    label = 'ULB_Kori1' if model == 'ULB_fETISh-KoriBU1' else model
    if label in seen_labels:
        continue
    seen_labels.add(label)
    highlight_handles.append(plt.Line2D([0], [0], color=col, lw=2, label=label))

fig.legend(handles=legend_elements() + highlight_handles, loc='upper center', ncol=3,
           bbox_to_anchor=(0.5, 1.06), fontsize=9, frameon=False)
fig.suptitle("Shelf area (% of 2015), all experiments, with individual models highlighted",
             fontsize=12, y=1.1)
plt.tight_layout()
plt.savefig("Figure_SI_shelfarea_allexps_highlight.pdf", dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close(fig)
print("saved Figure_SI_shelfarea_allexps_highlight.png")

# ---------------------------------------------------------------
# 5) All models, shelf area, all 4 experiments (rows) x all 6 regions (cols),
#    same layout/style as section 4, but highlighting the four IMAU_UFEMISM
#    models on top of the ordinary Group 1/Group 2 spaghetti + means, to
#    inspect their behaviour ahead of the Group 1 vs Group 2 classification
#    question for this model family
# ---------------------------------------------------------------
IMAU_MODELS = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4']
IMAU_HIGHLIGHT = {
    'IMAU_UFEMISM1': '#1b9e77',
    'IMAU_UFEMISM2': '#d95f02',
    'IMAU_UFEMISM3': '#7570b3',
    'IMAU_UFEMISM4': '#e7298a',
}

fig, axes = plt.subplots(len(ALL_EXPS), len(REGIONS), figsize=(13, 11), sharex=True)
for i, exp in enumerate(ALL_EXPS):
    for j, region in enumerate(REGIONS):
        ax = axes[i, j]
        plot_panel(ax, exp, region, "area", annotate=False, highlight=IMAU_HIGHLIGHT)
        ax.tick_params(axis='both', labelsize=11)
        if i == 0:
            ax.set_title(region, fontsize=13)
        if j == 0:
            ax.set_ylabel(f"{exp}\ndA (%)", fontsize=12)
        if i == len(ALL_EXPS) - 1:
            ax.set_xlabel("Year", fontsize=12)

imau_handles = [plt.Line2D([0], [0], color=IMAU_HIGHLIGHT[m], lw=2, label=m) for m in IMAU_MODELS]
fig.legend(handles=legend_elements() + imau_handles, loc='upper center', ncol=4,
           bbox_to_anchor=(0.5, 1.05), fontsize=11, frameon=False)
plt.tight_layout()
plt.savefig("Figure_SI_shelfarea_IMAU_UFEMISM.pdf", dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close(fig)
print("saved Figure_SI_shelfarea_IMAU_UFEMISM.png")

# ---------------------------------------------------------------
# 6) Overview figure: all models, shelf area, all 4 experiments (rows) x all 6
#    regions (cols), highlighting ULB_Kori1, UCSD_ISSM, UTAS_ElmerIce (the G2
#    fixed/near-fixed-front examples) together with all four IMAU_UFEMISM
#    models (the G1 retreat-only counterexample) on top of the ordinary
#    Group 1/Group 2 spaghetti + means. This is the combined figure for the
#    Methods "Calving group classification" illustrative examples.
# ---------------------------------------------------------------
OVERVIEW_HIGHLIGHT = {
    'UCSD_ISSM': '#FF69B4',            # pink
    'ULB_fETISh-KoriBU1': '#FF8C00',   # orange
    'UTAS_ElmerIce': '#E41A1C',        # red
    'IMAU_UFEMISM1': '#08519c',        # dark blue
    'IMAU_UFEMISM2': '#6a3d9a',        # purple
    'IMAU_UFEMISM3': '#000000',        # black
    'IMAU_UFEMISM4': '#b15928',        # brown
}
OVERVIEW_LABELS = {
    'UCSD_ISSM': 'UCSD_ISSM',
    'ULB_fETISh-KoriBU1': 'ULB_Kori1',
    'UTAS_ElmerIce': 'UTAS_ElmerIce',
    'IMAU_UFEMISM1': 'IMAU_UFEMISM1',
    'IMAU_UFEMISM2': 'IMAU_UFEMISM2',
    'IMAU_UFEMISM3': 'IMAU_UFEMISM3',
    'IMAU_UFEMISM4': 'IMAU_UFEMISM4',
}

fig, axes = plt.subplots(len(ALL_EXPS), len(REGIONS), figsize=(13, 11), sharex=True)
for i, exp in enumerate(ALL_EXPS):
    for j, region in enumerate(REGIONS):
        ax = axes[i, j]
        plot_panel(ax, exp, region, "area", annotate=False, highlight=OVERVIEW_HIGHLIGHT)
        ax.tick_params(axis='both', labelsize=11)
        if i == 0:
            ax.set_title(region, fontsize=13)
        if j == 0:
            ax.set_ylabel(f"{exp}\ndA (%)", fontsize=12)
        if i == len(ALL_EXPS) - 1:
            ax.set_xlabel("Year", fontsize=12)

overview_handles = [plt.Line2D([0], [0], color=col, lw=2, label=OVERVIEW_LABELS[m])
                     for m, col in OVERVIEW_HIGHLIGHT.items()]
fig.legend(handles=legend_elements() + overview_handles, loc='upper center', ncol=5,
           bbox_to_anchor=(0.5, 1.07), fontsize=11, frameon=False)
plt.tight_layout()
plt.savefig("Figure_SI_shelfarea_overview_highlight.pdf", dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close(fig)
print("saved Figure_SI_shelfarea_overview_highlight.png")
