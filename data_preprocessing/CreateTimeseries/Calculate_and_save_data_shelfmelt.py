#!/usr/bin/env python
# coding: utf-8

# # Calculate shelf area per region 

import numpy as np
import xarray as xr
import os
from matplotlib import pylab as plt
import pandas as pd

import function_sampling_Thermal_forcing_files as fn


mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars/' #write folder

pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']



dirPath = pth_ismip6

metadata_file ='Metadata.txt'

removeFileID = []


loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# comment line below if script is run for all models
loop_info= loop_info[loop_info['Model'].isin(['IMAU_UFEMISM1','IMAU_UFEMISM2','IMAU_UFEMISM3', 'IMAU_UFEMISM4'])].reset_index(drop = True)
variable_name = 'shelfmelt' #please don't change
variable = 'shelfmelt'


for fi in range(len(loop_info.index)):
    
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['mask_file'][fi]
    maskfile_ice= maskfile.replace('sftflf','sftgif')
    melt = xr.open_dataset(pth+'/'+maskfile.replace('sftflf','libmassbffl'))
    maskfile = maskfile.replace('sftflf_AIS_','new_sftflf_AIS_')
    scalefac_pth= f'/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_{grid}km.nc'
    #loadscaling mask
    scalefac_model= xr.open_dataset(scalefac_pth)
    
    print(model)
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #
    mask = xr.open_dataset(pth+'/'+maskfile)
    mask_ice = xr.open_dataset(pth+'/'+maskfile_ice)

    # create the new file     
    filename =  pth_helene+'/'+exp+f'/{variable}/computed_{variable}_AIS_'+model+'_'+exp+'.nc'
    d_example = xr.open_dataset(filename,decode_times=False )
    filename_orig_save= pth_helene+'/'+exp+f'/{variable}/orig_computed_{variable}_AIS_'+model+'_'+exp+'.nc'
    #save copy original data
    d_example.to_netcdf(filename_orig_save)

    d_calc =d_example.copy()
    d_calc.attrs ={}
    regions =['region_1','region_2','region_3']
    sectors =[]
    for i in range(18):
        sectors.append('sector_'+str(i+1))
    d_calc =d_calc.rename({f'{variable}':variable_name })
    d_calc[variable_name].attrs["units"] = "m2"
    for region in regions:
        varname_old = f'{variable}'+'_'+region
        varname_new = variable_name+'_'+region    
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "m2"
    for sector in sectors:
        varname_old = f'{variable}'+'_'+sector
        varname_new = variable_name+'_'+sector 
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "m2"
   

    # Loop through each variable in the dataset
    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'rhoi') and (var_name != 'rhow'):
            print(var_name)
            d_calc[var_name].values = np.zeros(len(d_example.time.values)) # assign you calc output not 0 :)

            time = range(len(mask.time)) # FIXME is this a good idea?
            # shelfarea = np.zeros([len(time)])
            if len(time)<len(d_example.time.values):
                print(f'Warning ice mask has shorter time then {variable}')
                print(len(time))
                print(len(d_example.time.values))

            shelfarea = np.zeros(len(d_example.time.values)) # assign you calc output not 0 :)


            if 'sector' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.sectors.values==si
            elif 'region' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.regions.values==si
            else:
                smask = d_region.sectors.values>0

            shelf_total_new=np.sum(melt['libmassbffl'].values* mask['sftflf'].values *mask_ice['sftgif'].values*scalefac_model.af2.values*smask,axis=(1,2))*((grid)*1000.0)**2 #in kg/s

# for ti in time:
# like Helene
# mask_slice = mask.isel(time=ti).sftflf.values[:] > 0.5
# mask_slice_ice = mask_ice.isel(time=ti).sftgif.values[:]  ==1
# msk = np.logical_and(smask,mask_slice,mask_slice_ice)

# if np.sum(msk) == 0:
# shelfarea[ti] = np.nan
# else:

# shelfarea[ti] = np.nansum(msk)*grid*grid*1e3*1e3
            d_calc[var_name].values = shelf_total_new[:len(d_calc[var_name])]
            
    # save data 
# path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    path_save =  pth_helene+'/'+exp+f'/{variable}/'
    save_name = f"jb_computed_{variable_name}_AIS_{model}_{exp}.nc"

    if not os.path.exists(f"{pth_calc_output}{exp}" ):
        os.mkdir(f"{pth_calc_output}{exp}" )
    if not os.path.exists(f"{pth_calc_output}{exp}/{variable_name}" ):
        os.mkdir(f"{pth_calc_output}{exp}/{variable_name}" )

    d_calc.to_netcdf(path_save+save_name)

    #break
    


print('Done :-D')

