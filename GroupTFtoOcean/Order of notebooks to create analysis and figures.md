### Melt sensitivity and BMB  Figures
1. Beware I deleted in ComputedScalarsHelene/exp/shelfmelt the `computed_shelfmelt_AIS_IMAU_UFEMISM1_expAE{}.nc `and `ctrlAE` data because I had new calculations for IMAU. You have to do that too in order to create the correct tables!
   Then you first need to find the time in which BMB is maximum (2), fit the linear regression of the avg BMB to TF (3) analysis the PRMSE of each method (4) and evaluate the best method of each model (5)
2. Create necessary tables with notebooks:
		`Timing of maximum BMB` for `tables/time_maximum_BMB.csv`
			`get cumulative BMB 2100 ,2200 2300`  for `./tables/cumulative_shelfmelt_2200.csv`
1. Pre-analisys for melt sensetivity with three different methods:
	1. `determine_the_linear_melt_sensitivity_factor-until_maximum_BMB` creates `tables/linear_meltsensetivity_from_BMB_max.csv` and `tables/linear_meltsensetivity_from_BMB_max_intercept.csv`
	2. `determine_the_linear_melt_sensitivity_factor-until_maximum_BMB-until 2100` for `'tables/linear_meltsensetivity_from_BMB_max_time2100.csv` and `tables/linear_meltsensetivity_from_BMB_max_intercept_time2100.csv
	3. `determine_the_linear_melt_sensitivity_factor-until_maximum_BMB-max_2K` for `tables/linear_meltsensetivity_from_BMB_max_TFuntil2K.csv` and ` tables/linear_meltsensetivity_from_BMB_max_intercept_TFuntil2K.csv`
	4. 
2. Predict BMB evolution and get PRMSE for each method (and model)
	1. notebook: `plot_predict variable BMB, recalculate MS times shelf area times TF`  creates `tables/prmse_recalc_linear_meltsensetivity_from_BMB_max.csv`
	2. notebook: `plot_predict variable BMB, recalculate MS times shelf area times TF-until 2100` creates `'tables/prmse_recalc_linear_meltsensetivity_from_BMB_max_time2100.csv`
	3. `plot_predict variable BMB, recalculate MS times shelf area times TF_TFuntil2K` creates `tables/prmse_recalc_linear_meltsensetivity_from_BMB_max_TFuntil2K.csv`
	
3. Evaluate best method;
	1. notebook `Evaluate Best method for MS by PRMSE for Paperfig 3 and SI` creates `figs/SI/Fig1.pdf` and ` 'figs/SI/predictive_mean_error_and_sensitivity.pdf'` `figs/paper/melt_sensetivity_region_m_K_a.pdf'`
		creates  table: `tables/best_ms_method_per_model.csv`
					`tables/ms_best_min_max.csv`
1. Plots 
	1. Fig. 3: `Fig 3 - MS all regions and part for Table 2`
	2. all SI MS figs for single regions : `Evaluate Best method for MS by PRMSE and SI-one region`
### Dyn slr and cumulative BMB Figure
1. Create necessary tables and plot
	1. notebook `get cumulative BMB 2100 ,2200 2300 ` and `get dyn slr 2100 ,2200 2300
	create tables `dslc_anom_2300.csv `(2100 and 2200) and `cumulative_shelfmelt_2300.csv` (2200 and 2300)
	2. notebook `Analyse dyn SLR  cum BMB-region of interest Fig. 4 .ipynb` created Figure `figs/paper/dynslr_cum_BMB.pdf`
	
### BMB evolution and grouping

1. Plot Figures for BMB evolution also with ms for SI
	1. `Fig 5 BMB evolution Figures SI `
## Shelf area evolution 
1 . Analayze evolution/grouping and ms SI
	`Fig SI shelf area evolution and grouping`




# MAIN Figure notebook (plots only)
in folder data_plotting/MAIN_FIGURES_notebooks

Fig. 1 -`Fig.1_Intropric_Recreating_IMSIP6_thermalforcing_BMB_melting_Figure-dSLR`
 Fig. 3: `Fig 3 - MS all regions and part for Table 2`
 Fig. 4: `Fig 4 shelf area evolution Figures`
 
 Fig 5, 6: `Fig 5 and 6 Main BMB evolution and scatter with MS`
 
 