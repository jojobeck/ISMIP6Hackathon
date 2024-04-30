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

def create_empty_dummy_interp_nc(d):
    if not isinstance(d, xr.Dataset):
        raise ValueError("Input 'd' must be an xarray Dataset.")

    # Check if 'sftflf' variable is present
    if 'sftflf' not in d:
        raise ValueError("Variable 'sftflf' not found in the Dataset.")

   
    d_calc =d.copy()
    d_calc.attrs ={}
    variable_name = 'thermalforcing_interp'
    varname_old = 'sftflf'
    d_calc =d_calc.rename({varname_old:variable_name })

    d_calc.thermalforcing_interp.values = np.zeros((d_calc.thermalforcing_interp.shape))



    return d_calc
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
# check time end of forcing file is either 2299 or 2300
    if np.size(tf.time.values) > 305:
        time_end = 2300
    else:
        time_end=2299

    for timei, time in enumerate(time_model[time_model <= time_end]):
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
    
def process_data_withInterp(d, draft, zbnd, tf, time_model,d_calc):
# check time end of forcing file is either 2299 or 2300
    if np.size(tf.time.values) > 305:
        time_end = 2300
    else:
        time_end=2299
    nt,nz,ny,nx = tf.thermal_forcing.shape
    for timei, time in enumerate(time_model[time_model <= time_end]):
        print(time)
        m_in = d.isel(time=timei).sftflf.values
        m = m_in > 0.5  # takes cells that are mostly floating
        shelf_draft = draft.isel(time=timei).base.values
        shelf_draft[~m] = 100
        shelf_draft_bin = np.empty(shelf_draft.shape)
        shelf_draft_bin[:] = np.nan
        shelf_draft_upper = np.empty(shelf_draft.shape)
        shelf_draft_upper[:] = np.nan
        shelf_draft_lower = np.empty(shelf_draft.shape)
        shelf_draft_lower[:] = np.nan

        # This will be the relevant thermal forcing
        shelf_draft_bin_temp = np.empty(shelf_draft.shape)
        shelf_draft_tf_upper = np.empty(shelf_draft.shape)
        shelf_draft_tf_lower = np.empty(shelf_draft.shape)
        
        # find corretc time of thermal forcing
        yr = int(time_model[timei])
        tf_time = tf.loc[{'time': slice(str(yr) + '-01', str(yr) + '-12')}].thermal_forcing.values[0,]
        tf_bnds = np.empty([nz+1,ny,nx])
        tf_bnds[1:-1,:,:] = tf_time[0:-1,:,:] + 0.5*(tf_time[1:,:,:]-tf_time[0:-1,:,:])
        tf_bnds[0,:,:] = tf_time[0,:,:]
        tf_bnds[-1,:,:] = tf_time[-1,:,:]

        # create a mask of the depth indices of the shelf draft
        for i, zbnd_i in enumerate(zbnd):
            zmax = zbnd_i[0]
            zmin = zbnd_i[-1]
            ind_mask = np.nonzero(np.logical_and(shelf_draft < zmax, shelf_draft >= zmin))
            shelf_draft_bin[ind_mask] = i
            shelf_draft_upper[ind_mask] = zmax
            shelf_draft_lower[ind_mask] = zmin

        # find the thermal forcing at the correct deoth
        for i, z in enumerate(tf.thermal_forcing.z.values):
            maski = shelf_draft_bin == i
            shelf_draft_bin_temp[maski] = tf_time[i, maski]
            shelf_draft_tf_upper[maski] = tf_bnds[i,maski]
            shelf_draft_tf_lower[maski] = tf_bnds[i+1,maski]
        # now interpolate to find the weighted thermal forcing
        thermalforcing_interp_in = weightedInterpTemp(shelf_draft_upper,shelf_draft_lower,shelf_draft,shelf_draft_tf_upper,shelf_draft_tf_lower)
        d_calc.isel(time=timei).thermalforcing_interp.values[:] = thermalforcing_interp_in.copy()

    return d_calc 

def weightedInterpTemp(zmax,zmin,depth,tmax,tmin):
    bin_thickness = zmax-zmin # zmax is top cell, zmin is bottom
    depth_below_topcell = zmax - depth #calculate how far down the depth is from the top cell
    frac_thickness = depth_below_topcell/bin_thickness #fractional thickness through bin

    t_interp = (tmin-tmax)*frac_thickness + tmax   # calculate the deltaT between top and bottom,
                                        # scale by the fraction of deltaT equivalent to the fractional thickness,
                                        # add to the top cell temp
    return t_interp

def get_time_model(d):
    time_model = np.empty(len(d['time'].values))
    for i in range(len(d['time'].values)):
        time_model[i] = (int(d.sftflf['time'].values[i].strftime().split('-')[0]))
    
    return(time_model)
    
