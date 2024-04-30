import numpy as np
import xarray as xr
import pandas as pd
import os
def create_empty_dummies_nc(d):
    if not isinstance(d, xr.Dataset):
        raise ValueError("Input 'd' must be an xarray Dataset.")

    # Check if 'sftflf' variable is present
    if 'sftflf' not in d:
        raise ValueError("Variable 'sftflf' not found in the Dataset.")

   
    d_calc =d.copy()
    d_calc.attrs ={}
    variable_name = 'thermalforcing_bin'
    varname_old = 'sftflf'
    d_calc =d_calc.rename({varname_old:variable_name })

    d_calc.thermalforcing_bin.values = np.zeros((d_calc.thermalforcing_bin.shape))
    d_bin =d.copy()
    d_bin.attrs={}
    d_bin =d_bin.rename({varname_old:'bin_number' })
    d_bin.bin_number.values = np.zeros((d_calc.thermalforcing_bin.shape))



    return d_calc, d_bin

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created successfully.")
    else:
        print(f"Directory '{directory}' already exists.")

def get_outpath_tf_filend(results,i):
    pend = results.iloc[i].path.split('ismip6_2300')[-1] 
    makedir=results.iloc[i].path.split('ismip6_2300')[0] +'ismip6_hackathon/dataset_2HD' +pend
    create_directory_if_not_exists(makedir)    
    file_end = results.iloc[i].mask_file.split('sftflf')[-1]
    return makedir, file_end 


def process_data(d, draft, zbnd, tf, time_model,d_calc,d_bin):
    for timei, time in enumerate(time_model[time_model <= 2300]):
        print(time)
        m_in = d.isel(time=timei).sftflf.values
        m = m_in > 0.5  # takes cells that are mostly floating
        shelf_draft = draft.isel(time=timei).base.values
        shelf_draft[~m] = 100
        shelf_draft_bin = np.empty(shelf_draft.shape)
        shelf_draft_bin[:] = np.nan

        # This will be the relevant thermal forcing
        shelf_draft_bin_temp = np.empty(shelf_draft.shape)
        shelf_draft_bin_temp[:] = np.nan

        # create a mask of the depth indices of the shelf draft
        for i, zbnd_i in enumerate(zbnd):
            zmax = zbnd_i[0]
            zmin = zbnd_i[-1]
            ind_mask = np.nonzero(np.logical_and(shelf_draft < zmax, shelf_draft >= zmin))
            shelf_draft_bin[ind_mask] = i
        d_bin.isel(time=timei).bin_number.values[:] = shelf_draft_bin.copy()

        yr = int(time_model[timei])
        tf_time = tf.loc[{'time': slice(str(yr) + '-01', str(yr) + '-12')}].thermal_forcing.values[0,]
        for i, z in enumerate(tf.thermal_forcing.z.values):
            maski = shelf_draft_bin == i
            shelf_draft_bin_temp[maski] = tf_time[i, maski]
        d_calc.isel(time=timei).thermalforcing_bin.values[:] = shelf_draft_bin_temp.copy()

    return d_calc ,d_bin
    
def get_time_model(d):
    time_model = np.empty(len(d['time'].values))
    for i in range(len(d['time'].values)):
        time_model[i] = (int(d.sftflf['time'].values[i].strftime().split('-')[0]))
    
    return(time_model)
    
