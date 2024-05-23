
import numpy as np
import xarray as xr
import os
import pandas as pd
import time
import function_sampling_Thermal_forcing_files as fn

res =16
mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
forcing_data_path = '/home/565/jb1863/ismip6_2300/'
# Example experiment and corresponding ocean forcing names this can be drawn from the INFOS Dataframe

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
expnames_plot = ['CCSM4','HadGEM2','CESM2','UKESM','ctrl']
expnames_path = ['/1995-2300/CCSM4_RCP85','/1995-2299/HadGEM2-ES_RCP85','/1995-2299/CESM2-WACCM_SSP585','/1995-2300/UKESM1-0-LL_SSP585','climatology_from_obs_1995-2017/']


removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

metadata_file = 'Metadata.txt'
# Generate loop_info DataFrame
loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# do 8km for now
loop_info=loop_info[loop_info.Grid ==res].reset_index(drop=True)


start_time = time.time()
# for i_e,expo in enumerate(expnames):
for i_e in range(4,5): 
    experiment =expnames[i_e]
    name = expnames_path[i_e]
    results_exp =loop_info[loop_info['Experiment'] == experiment].reset_index(drop =True)
    print('running all for experiment',experiment)
    #return correct foricng for exp kind
    print('thermal forcing ' ,name)
    tf,zbnd = fn.tf_z_for_exp(forcing_data_path,experiment,name,res)
    for i in range(0,results_exp.shape[0]):
        print(i)
        path_save,printm = fn.process_model_tf(results_exp,i,experiment,tf,zbnd)
        print('save data to', path_save);
        print('save data to', printm);
end_time = time.time()
# Calculate elapsed time
elapsed_time = end_time - start_time

print("Elapsed time:", elapsed_time, "seconds")
