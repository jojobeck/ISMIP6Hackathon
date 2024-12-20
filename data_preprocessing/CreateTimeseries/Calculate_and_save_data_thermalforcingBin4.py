#!/usr/bin/env python
# coding: utf-8

# # Calculate average thermal forcing per region 


# Calculate average thermal forcing, melting and draft depth

import numpy as np
import xarray as xr
import os
import pandas as pd
import function_sampling_Thermal_forcing_files as fn


mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars/' #write folder
pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
# removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
# 'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

metadata_file = 'Metadata.txt'
# Generate loop_info DataFrame
loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# do 8km for now
# do 8km for now
loop_info=loop_info[loop_info.Grid ==8].reset_index(drop=True)

# CHeck why VUW PISM has weired thermal forcing
# loop_info = loop_info.iloc[41:]
# loop_info = loop_info.reset_index(drop = True)


variable_name = 'thermalforcingBin' #please don't change


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
    tf = xr.open_dataset(path_save +'/shelf_thermalforcingBin'+file_end)
    #attention, VUW PISM fill nan values was not detected as np.nan from xarray,
    #there replace with nan 
    mn = tf.thermalforcing_bin.values <=-1e6 
    tf.thermalforcing_bin.values[mn] =np.nan
    
 
    #ensure same time length and truncate dataset if neccessary
    mask_trun,tf_trun,tff=fn.assure_minimum_same_timelength(mask,tf,tf)

    d_calc=fn.create_empty_dummy_t(variable_name,mask_trun.time)
   

    # Loop through each variable in the dataset
    mask_floating = mask_trun.sftflf.values > 0 
    # adjust ice mask
    m_ice = mask_ice.sftgif.values >0; #because in the end Helen set all values bigger 1 to 1 and allows values between 0 and 1
    mask_floating =mask_floating & m_ice
    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'z') :
            print(var_name)
            smask = sector_masks[int(var_name.split('_')[-1])] if 'sector' in var_name else region_masks[int(var_name.split('_')[-1])] if 'region' in var_name else d_region.sectors.values > 0

            m_flotregion = mask_floating & smask  # mask region and shelf extent over time
            masked_values = np.where(m_flotregion, tf_trun.thermalforcing_bin.values, np.nan)
                                                                                                                                                                                                                                              
            # Compute the mean along the time axis, ignoring NaNs
            mean_values = np.nanmean(masked_values, axis=(1, 2))

            d_calc[var_name].values = mean_values
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"


    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    d_calc.to_netcdf(path_save+save_name)

    #break
    
    
print('Done :-)')


