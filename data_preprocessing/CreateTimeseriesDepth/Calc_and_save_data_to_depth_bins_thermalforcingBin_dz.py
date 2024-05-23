
import concurrent.futures
import numpy as np
import xarray as xr
import os
import pandas as pd
import function_sampling_Thermal_forcing_files as fn


mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
pth_calc_output = mnt_pth + 'ismip6_hackathon/ComputedScalars_bin/' #write folder
pth_ismip6 =mnt_pth + 'ismip6_2300/' #read folder
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'

expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

metadata_file = 'Metadata.txt'
# Generate loop_info DataFrame
loop_info = fn.create_loop_info(mnt_pth, expnames, removeFileID, metadata_file)
# do 8km for now
# loop_info=loop_info[loop_info.Grid ==16].reset_index(drop=True)

# get depth bins which will be our new dimension for the dataset
d_z =  xr.open_dataset(pth_ismip6 +'/climatology_1995-2014/UKESM1-0-LL_SSP585_temperature_8km_x_60m.nc')
depth =d_z.z.values
max_bin_number = depth.size -1
variable_name = 'thermalforcingBin' #please don't change






def process_file(fi):
    
    pth = loop_info['path'][fi]
    grid = loop_info['Grid'][fi]
    exp = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    maskfile = loop_info['mask_file'][fi]
    basefile = loop_info['base_file'][fi]
    
    
    # Load data to calc file
    d_region = xr.open_dataset(pth_ismip6+'masks/sectors_'+str(grid)+'km.nc') #

    # Precompute sector and region masks to avoid recomputation in the loop
    sector_masks = {si: d_region.sectors.values == si for si in range(1, 19)}
    region_masks = {si: d_region.regions.values == si for si in range(1, 4)}

    mask = xr.open_dataset(pth+'/'+maskfile)
    
    
    #get path save and names for the new nc data
    path_save, file_end=fn.get_outpath_tf_filend(loop_info,fi)
    bins = xr.open_dataset(path_save +'/shelf_numberBin'+file_end)
    tf = xr.open_dataset(path_save +'/shelf_thermalforcingBin'+file_end)

                           

    #ensure same time length and truncate dataset if neccessary
    mask_trun,bins_trun,tf_trun=fn.assure_minimum_same_timelength(mask,bins,tf)
    # Create the mask for floating cells over time
    mask_floating = mask_trun.sftflf.values > 0.5
    # create emptty dataset
    d_calc=fn.create_empty_dummy_zt(variable_name,depth,tf_trun.time)

    # Create a dictionary to store masks for each bin
    dic_mask_bin = {bini: bins_trun.bin_number.values == bini for bini in range(max_bin_number + 1)}


    for var_name in d_calc.variables:
        if (var_name != 'time') and (var_name != 'z') :
            print(var_name)
            smask = sector_masks[int(var_name.split('_')[-1])] if 'sector' in var_name else region_masks[int(var_name.split('_')[-1])] if 'region' in var_name else d_region.sectors.values > 0

            m_flotregion = mask_floating & smask  # mask region and shelf extent over time

            for bini in range(max_bin_number + 1):
                m = dic_mask_bin[bini] & m_flotregion  # mask in time(x,y) for each bin
                # Apply the mask to base.base.values
                masked_values = np.where(m, tf_trun.thermalforcing_bin.values, np.nan)

                # Compute the mean along the time axis, ignoring NaNs
                mean_values = np.nansum(masked_values, axis=(1, 2))
                d_calc[var_name][bini, :] = mean_values

                              
            
    # save data 
    path_save = f"{pth_calc_output}{exp}/{variable_name}/"
    save_name = f"computed_{variable_name}_AIS_{model}_{exp}.nc"
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    d_calc.to_netcdf(path_save+save_name)
    return f"Done processing {model}"

# Main execution
if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_file, fi) for fi in range(len(loop_info.index))]
        for future in concurrent.futures.as_completed(futures):
            print(future.result())

print('All Done :-D')    



