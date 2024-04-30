import numpy as np
import xarray as xr
import pandas as pd
import os
import time
import function_sampling_Thermal_forcing_files as fn
#################################### 
#  0. Path, expereiment names etc
############################################
#{{{ oath experiment etc 
forcing_data_path = '/home/565/jb1863/ismip6_2300/'
# Example experiment and corresponding ocean forcing names this can be drawn from the INFOS Dataframe

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05']
expnames_plot = ['CCSM4','HadGEM2','CESM2','UKESM']
expnames_path = ['/1995-2300/CCSM4_RCP85','/1995-2299/HadGEM2-ES_RCP85','/1995-2299/CESM2-WACCM_SSP585','/1995-2300/UKESM1-0-LL_SSP585']

i_e= 0 #choose experiment 2 -4 ,corresponds to the right climate modelforcign 
experiment =expnames[i_e]
name = expnames_path[i_e]

mnt_pth = '/home/565/jb1863/' #mount path,actual path if script is run on cluster
pth_calc_output = mnt_pth + 'ismip6_hackathon/' #write folder
dirPath =mnt_pth + 'ismip6_2300' #read folder
# }}}
###########################
# 1. Get forcing data
#############################
#{{{ forcing data 
tf = xr.open_dataset(
    forcing_data_path+name+"_thermal_forcing_8km_x_60m.nc")


# Get the time from the forcing

time_forcing = np.empty(len(tf['time'].values))
for i in range(len(tf['time'].values)):
    time_forcing[i] = (int(tf['time'].values[i].strftime().split('-')[0]))
zbnd = tf.get('z_bnds').values[0,:]
# }}}
#############################
#2. Create Table
######################
#{{{ table 
# create table to bu used the loop of neccesary model output to be used stroe in results
metadata = pd.read_csv('Metadata.txt')






expFilter = [experiment] # specify an experiment
# gridFilter = ['04', '4', '8', '_08']
removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

models = metadata.loc[( metadata['Experiment'].isin(expFilter)) &
        ~(metadata['fileID'].isin(removeFileID)),]


results_data = []
for i in range(models.shape[0]):

    modelPath = models.iloc[i].fileID
    exp = models.iloc[i].Experiment
    modelFiles = os.listdir(dirPath + '/' + modelPath)
    expDirPos = [i for i, dirName in enumerate(modelFiles) if dirName.startswith(exp)]



    # For each relevant directory, loop
    for expdir_subgrid in expDirPos:
        expDir = modelFiles[expdir_subgrid]
        expFiles = os.listdir(dirPath + '/' + modelPath + '/' + expDir)

        # Get grid size from expDir: format expAE01_04 where 04 is the grid size and convert to int
        grid = int(expDir.split('_')[1])

        # Find position of relevant files
        draftPos = [i for i, fileName in enumerate(expFiles) if fileName.startswith('base')]
        meltPos = [i for i, fileName in enumerate(expFiles) if fileName.startswith('libmassbffl')]
        maskPos = [i for i, fileName in enumerate(expFiles) if fileName.startswith('sftflf')]

        # Construct dictionary and append to list
        results_data.append({'base_file': expFiles[draftPos[0]],
                             'mask_file': expFiles[maskPos[0]],
                             'path': os.path.join(dirPath, modelPath, expDir),
                             'Experiment': exp,
                             'Model': modelPath,
                             'Grid': grid})


pre_results = pd.DataFrame(results_data)
#HAack only use 8km file ,todo subsample for 4 and 16
results = pre_results[pre_results['Grid'] == 8]
# }}}
#######################
# 3. Run the loop preprocessing
###############################
#{{{ runing script 
start_time = time.time()
for j in range(results.shape[0]):
    i=j+1 #already made the first one
    maskData = xr.open_dataset(results.iloc[i].path + '/' + results.iloc[i].mask_file)

    draftData = xr.open_dataset(results.iloc[i].path + '/' + results.iloc[i].base_file)

    d_calc,d_bin = fn.create_empty_dummies_nc(maskData)
    time_model = fn.get_time_model(maskData)

    #get path save and names for the new nc data
    path_save, file_end=fn.get_outpath_tf_filend(results,i)
    #process data
    d_calcnew,d_binnew = fn.process_data(maskData, draftData, zbnd, tf, time_model,d_calc,d_bin)
    print('save data to', path_save);
    d_calcnew.to_netcdf(path_save +'/shelf_thermalforcingBin'+file_end)
    d_binnew.to_netcdf(path_save +'/shelf_numberBin'+file_end)
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time

print("Elapsed time:", elapsed_time, "seconds")
#}}}


