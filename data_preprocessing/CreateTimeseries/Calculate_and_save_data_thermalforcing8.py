#!/usr/bin/env python
# coding: utf-8

# # Calculate average thermal forcing per region 


# Calculate average thermal forcing, melting and draft depth

import numpy as np
import xarray as xr
import os
from matplotlib import pylab as plt
import pandas as pd
import function_sampling_Thermal_forcing_files as fn
grid = 8

mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars/' #write folder
pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
removeFileID = []

metadata_file = 'Metadata.txt'
# Generate loop_info DataFrame
loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# do 8km for now
loop_info=loop_info[loop_info.Grid ==grid].reset_index(drop=True)

# CHeck why VUW PISM has weired thermal forcing
# loop_info = loop_info[loop_info['Model'].isin(['DOE_MALI_2','DOE_MALI_3' ])].reset_index(drop = True)

scalefac_pth= f'/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_{grid}km.nc'
#loadscaling mask
scalefac_model= xr.open_dataset(scalefac_pth)

variable_name = 'thermalforcing' #please don't change


for fi in range(len(loop_info.index)):
    
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['mask_file'][fi]
    maskfile_ice= maskfile.replace('sftflf','sftgif')
    
    print(model)
    
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #
    sector_masks = {si: d_region.sectors.values == si for si in range(1, 19)}
    region_masks = {si: d_region.regions.values == si for si in range(1, 4)}

    mask = xr.open_dataset(pth+'/'+maskfile)
    mask_ice = xr.open_dataset(pth+'/'+maskfile_ice)
    
    path_save, file_end=fn.get_outpath_tf_filend(loop_info,fi)
    tf = xr.open_dataset(path_save +'/shelf_thermalforcingInterp'+file_end)
# mn = tf.thermalforcing_interp.values <=-1e6
# tf.thermalforcing_interp.values[mn] =np.nan
    
 
    #ensure same time length and truncate dataset if neccessary
    mask_trun,tf_trun,mask_ice_trun=fn.assure_minimum_same_timelength(mask,tf,mask_ice)

    d_calc=fn.create_empty_dummy_t(variable_name,mask_trun.time)
   

    # Loop through each variable in the dataset
    mask_floating = mask_trun.sftflf.values
    m_ice = mask_ice_trun.sftgif.values
    #mask handling as Helene
    m_ice[m_ice<0]=0
    m_ice[m_ice>1]=1
    mask_floating[mask_floating<0]=0
    mask_floating[mask_floating>1]=1
    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'z') :
            print(var_name)
            smask = sector_masks[int(var_name.split('_')[-1])] if 'sector' in var_name else region_masks[int(var_name.split('_')[-1])] if 'region' in var_name else d_region.sectors.values > 0

            TF_total=np.nansum(tf_trun.thermalforcing_interp.values*mask_floating *m_ice*scalefac_model.af2.values*smask,axis=(1,2)) 
            shelf_total=np.nansum(mask_floating *m_ice*scalefac_model.af2.values*smask,axis=(1,2)) #in kg/s
            mean_values = TF_total/shelf_total
            d_calc[var_name].values = mean_values
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"


    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    d_calc.to_netcdf(path_save+save_name)

    
    
print('Done :-)')


