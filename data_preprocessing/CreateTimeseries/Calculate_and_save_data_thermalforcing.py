#!/usr/bin/env python
# coding: utf-8

# # Calculate average thermal forcing per region 

# In[1]:


# Calculate average thermal forcing, melting and draft depth

import numpy as np
import xarray as xr
import os
from matplotlib import pylab as plt
import pandas as pd


# In[2]:


mnt_pth = '/mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars/' #write folder

pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05']


# In[3]:


dirPath = pth_ismip6

metadata = pd.read_csv('Metadata.txt')

removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

models = metadata.loc[( metadata['Experiment'].isin(expnames)) &
        ~(metadata['fileID'].isin(removeFileID)),]

# initialize empty dataframe to store processed_data, containing file, pathtofile, experiment, model, grid
loop_info = pd.DataFrame(columns=['basefile', 'libmassbfflfile', 'maskfile' ,
                                 'path', 'Experiment', 'Model', 'Grid'])

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
    
        #save into processed_data
        loop_info = loop_info.append({#'basefile': expFiles[draftPos[0]],
                                    #'thermalforcingfile': expFiles[meltPos[0]],
                                    'maskfile': expFiles[maskPos[0]],
                                    'path': dirPath + '/' + modelPath + '/' + expDir,
                                    'Experiment': exp,
                                    'Model': modelPath,
                                    'Grid': grid},
                                    ignore_index=True)


# In[4]:


loop_info


# In[ ]:


# Iterate over loop info

variable_name = 'thermalforcing' #please don't change


for fi in range(len(loop_info.index)):
    
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['maskfile'][fi]
    
    print(model)
    
    # FIXME REMOVE THIS PART ONCE WE HAVE MORE DATA!
    if 'DC' in model:
        print('Calc')
    else:
        continue
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #
    mask = xr.open_dataset(pth+'/'+maskfile)
    
    tf = xr.open_dataset(mnt_pth+'/ismip6_hackathon/dataset_2HD/'+
                         model+'/'+exp+'_0'+str(grid)+'/'+
                         maskfile.split('/')[-1].replace('sftflf','shelf_thermalforcingBin')
                              )
    
 

    # create the new file     
    filename =  pth_helene+'/'+exp+'/shelfmelt/computed_shelfmelt_AIS_'+model+'_'+exp+'.nc'
    d_example = xr.open_dataset(filename,decode_times=False )
    d_calc =d_example.copy()
    d_calc.attrs ={}
    regions =['region_1','region_2','region_3']
    sectors =[]
    for i in range(18):
        sectors.append('sector_'+str(i+1))
    d_calc =d_calc.rename({'shelfmelt':variable_name })
    d_calc[variable_name].attrs["units"] = "K"
    for region in regions:
        varname_old = 'shelfmelt'+'_'+region
        varname_new = variable_name+'_'+region    
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "K"
    for sector in sectors:
        varname_old = 'shelfmelt'+'_'+sector
        varname_new = variable_name+'_'+sector 
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "K"
   

    # Loop through each variable in the dataset
    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'rhoi') and (var_name != 'rhow'):
            print(var_name)
            d_calc[var_name].values = np.zeros(len(d_example.time.values)) # assign you calc output not 0 :)

            time = range(len(mask.time)) # FIXME is this a good idea?
            thermalforcing = np.zeros([len(time)])

            if 'sector' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.sectors.values==si
            elif 'region' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.regions.values==si
            else:
                smask = d_region.sectors.values>0

            for ti in time:
                tf_slice = tf.isel(time=ti).thermalforcing_bin.values[:]
                thermalforcing[ti] = np.nanmean(tf_slice[smask])
              
            d_calc[var_name].values = thermalforcing[:len(d_calc[var_name])]
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"

    if not os.path.exists(f"{pth_calc_output}{exp}" ):
        os.mkdir(f"{pth_calc_output}{exp}" )
    if not os.path.exists(f"{pth_calc_output}{exp}/{variable_name}" ):
        os.mkdir(f"{pth_calc_output}{exp}/{variable_name}" )

    d_calc.to_netcdf(path_save+save_name)

    #break
    
    


# In[ ]:


print('Done :-)')


# In[ ]:





# In[ ]:




