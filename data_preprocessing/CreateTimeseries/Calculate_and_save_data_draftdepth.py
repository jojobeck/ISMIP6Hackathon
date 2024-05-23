
#!/usr/bin/env python
# coding: utf-8

# # Calculate shelf area per region 

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


dirPath = pth_ismip6

metadata_file ='Metadata.txt'

removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# Iterate over loop info

variable_name = 'averageDraftDepth' #please don't change


for fi in range(len(loop_info.index)):
    
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['mask_file'][fi]
    basefile = loop_info['base_file'][fi]
    
    print(model)
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #
    mask = xr.open_dataset(pth+'/'+maskfile)
    base = xr.open_dataset(pth+'/'+basefile)

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
    d_calc[variable_name].attrs["units"] = "m"
    for region in regions:
        varname_old = 'shelfmelt'+'_'+region
        varname_new = variable_name+'_'+region    
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "m"
    for sector in sectors:
        varname_old = 'shelfmelt'+'_'+sector
        varname_new = variable_name+'_'+sector 
        d_calc =d_calc.rename({varname_old:varname_new })
        d_calc[varname_new].attrs["units"] = "m"
   

    # Loop through each variable in the dataset
    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'rhoi') and (var_name != 'rhow'):
            print(var_name)
            d_calc[var_name].values = np.zeros(len(d_example.time.values)) # assign you calc output not 0 :)

            time_i= range(len(mask.time)) # FIXME is this a good idea?
            shelfdraft = np.zeros([len(time_i)])

            if 'sector' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.sectors.values==si
            elif 'region' in var_name:
                si = int(var_name.split('_')[-1]) # get sector number
                smask = d_region.regions.values==si
            else:
                #entire AIS
                smask = d_region.sectors.values>0

            for ti in time_i:
                mask_slice = mask.isel(time=ti).sftflf.values[:] > 0.5
                #mask_slice[mask_slice<0.5] = 0 # MAKE sure to only include cells that are more floating than grounded
                msk = np.logical_and(smask,mask_slice) 
                base_slice = base.isel(time=ti).base.values[:]
                # Check if the slice is empty
                if np.sum(msk) == 0:
                    shelfdraft[ti] = np.nan
                else:
                    shelfdraft[ti] = np.nanmean(base_slice[msk])

                
            d_calc[var_name].values = shelfdraft[:len(d_calc[var_name].values)]
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"

    if not os.path.exists(f"{pth_calc_output}{exp}" ):
        os.mkdir(f"{pth_calc_output}{exp}" )
    if not os.path.exists(f"{pth_calc_output}{exp}/{variable_name}" ):
        os.mkdir(f"{pth_calc_output}{exp}/{variable_name}" )

    d_calc.to_netcdf(path_save+save_name)

    #break
    
    


print('Done :-D')

#d_calc.averageDraftDepth_region_1.plot()





