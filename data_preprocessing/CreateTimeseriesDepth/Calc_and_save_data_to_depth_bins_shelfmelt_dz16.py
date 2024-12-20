
import concurrent.futures
import numpy as np
import xarray as xr
import os
import pandas as pd
import function_sampling_Thermal_forcing_files as fn

grid =16
mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars_bin/' #write folder
pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
# expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05']
# removeFileID = ['IMAU_UFEMISM1','IMAU_UFEMISM2','IMAU_UFEMISM3', 'IMAU_UFEMISM4']
removeFileID = []


metadata_file = 'Metadata.txt'
# Generate loop_info DataFrame
loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# do 8km for now
loop_info=loop_info[loop_info.Grid ==grid].reset_index(drop=True)
# loop_info = fn.correct_runs_missing_exp(loop_info)
# 2nd time runninf correcting the errors in the metadata.txt
# loop_info = fn.correct_runs_missing_exp(loop_info)
#3rdtime runnign chechking data in Yelmo
# loop_info=loop_info[loop_info.Model == 'UCM_Yelmo'].reset_index(drop = True)


# loop_info = loop_info[loop_info['Model'].isin(['DOE_MALI_2','DOE_MALI_3' ])].reset_index(drop = True)
# get depth bins which will be our new dimension for the dataset
d_z =  xr.open_dataset(pth_ismip6 +'/climatology_1995-2014/UKESM1-0-LL_SSP585_temperature_8km_x_60m.nc')
depth =d_z.z.values
max_bin_number = depth.size -1

scalefac_pth= f'/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_{grid}km.nc'
#loadscaling mask
scalefac_model= xr.open_dataset(scalefac_pth)


variable_name = 'shelfmelt' #please don't change

dA = (grid *1000 *grid *1000)


# def process_file(fi):
    
for fi in range(len(loop_info.index)):
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['mask_file'][fi]
    if model.startswith('IMAU'):
        maskfile_ice= maskfile.replace('new_sftflf','sftgif')
    else:
        maskfile_ice= maskfile.replace('sftflf','sftgif')
    basefile = loop_info['base_file'][fi]
    
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #

    # Precompute sector and region masks to avoid recomputation in the loop
    sector_masks = {si: d_region.sectors.values == si for si in range(1, 19)}
    region_masks = {si: d_region.regions.values == si for si in range(1, 4)}

    mask = xr.open_dataset(pth+'/'+maskfile)
    mask_ice = xr.open_dataset(pth+'/'+maskfile_ice)
    
    tf = xr.open_dataset(pth+'/'+basefile.replace('base','libmassbffl'))
    
    #get path save and names for the new nc data
    path_save, file_end=fn.get_outpath_tf_filend(loop_info,fi)
    bins = xr.open_dataset(path_save +'/shelf_numberBin'+file_end)
    filename =  pth_helene+'/'+exp+'/shelfmelt/computed_shelfmelt_AIS_'+model+'_'+exp+'.nc'

                           

    #ensure same time length and truncate dataset if neccessary
    mask_trun,bins_trun,tf_trun,mask_ice_trun=fn.assure_minimum_same_timelength4(mask,bins,tf,mask_ice)
# mask_trun,bins_trun,mask_ice_trun=fn.assure_minimum_same_timelength(mask,bins,mask_ice)
    # Create the mask for floating cells over time
    # create emptty dataset
    d_calc=fn.create_empty_dummy_zt(variable_name,depth,tf_trun.time)

    # Create a dictionary to store masks for each bin
    dic_mask_bin = {bini: bins_trun.bin_number.values == bini for bini in range(max_bin_number + 1)}

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

            m_floatregion = mask_floating * smask  *m_ice # mask region and shelf extent over time

            for bini in range(max_bin_number + 1):
                m = dic_mask_bin[bini]  # mask in time(x,y) for each bin
                # Apply the mask to base.base.values
                # masked_values = np.where(m, tf_trun.libmassbffl.values, np.nan) *dA

                # Compute the sum  along the time axis, ignoring NaNs
                shelfmelttotal=np.nansum(m_floatregion*m* scalefac_model.af2.values*tf_trun.libmassbffl.values*dA,axis=(1,2))#in kg/s
                d_calc[var_name][bini, :] = shelfmelttotal

                              
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    d_calc.to_netcdf(path_save+save_name)
    #close dataset to save memory 
    d_region.close()
    mask.close()
    bins.close()
    tf.close()
    mask_trun.close()
    bins_trun.close()
    tf_trun.close()
    d_calc.close()

# return f"Done processing {model}"

# Main execution
# with concurrent.futures.ProcessPoolExecutor() as executor:
# futures = [executor.submit(process_file, fi) for fi in range(len(loop_info.index))]
# for future in concurrent.futures.as_completed(futures):
# print(future.result())


print('All Done :-D')    



