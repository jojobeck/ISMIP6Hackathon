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
    if not np.issubdtype(d_calc[variable_name].values.dtype, np.floating):
        d_calc[variable_name] =d_calc[variable_name].astype(float)

    d_calc.thermalforcing_bin.values = np.full(d['sftflf'].shape, np.nan)
    d_bin =d.copy()
    d_bin.attrs={}
    d_bin =d_bin.rename({'sftflf':'bin_number' })
    if not np.issubdtype(d_bin.bin_number.values.dtype, np.floating):
        d_bin['bin_number'] =d_bin.bin_number.astype(float)
    d_bin.bin_number.values = np.full(d['sftflf'].shape, np.nan)



    return d_calc, d_bin
def process_model_tf(results_exp,i,experiment,tf,zbnd):

    maskData = xr.open_dataset(results_exp.iloc[i].path + '/' + results_exp.iloc[i].mask_file)
    maskfile = results_exp.iloc[i].mask_file
    maskfile_ice= maskfile.replace('sftflf','sftgif')
    mask_ice = xr.open_dataset(results_exp.iloc[i].path +'/'+maskfile_ice)
    no_ice=mask_ice.sftgif.values ==0

    # maskData.sftflf.values[no_ice]=np.nan


    draftData = xr.open_dataset(results_exp.iloc[i].path + '/' + results_exp.iloc[i].base_file)
    draftData.base.values[no_ice]=np.nan

    d_calc,d_bin = create_empty_dummies_nc(maskData)

    maskData_trun,draftData_trun,d_calctrun=assure_minimum_same_timelength(maskData,draftData,d_calc)
    maskData_trun,draftData_trun,d_bintrun=assure_minimum_same_timelength(maskData,draftData,d_bin)
    time_model = get_time_model(maskData_trun)
    #get path save and names for the new nc data
    path_save, file_end=get_outpath_tf_filend(results_exp,i)
    #process data
    if experiment!='ctrlAE':
        d_calcnew,d_binnew,printm = process_data_mask(maskData_trun, draftData_trun, zbnd, tf, time_model,d_calctrun,d_bintrun)
    else:
        d_calcnew,d_binnew = process_data_ctrl_mask(maskData_trun, draftData_trun, zbnd, tf, time_model,d_calctrun,d_bintrun)
        printm =''
    d_calcnew.to_netcdf(path_save +'/shelf_thermalforcingBin'+file_end)
    d_binnew.to_netcdf(path_save +'/shelf_numberBin'+file_end)
    return(path_save,printm)
def tf_z_for_exp(forcing_data_path,experiment,name,res):
    if experiment!='ctrlAE':
        tf = xr.open_dataset(
            forcing_data_path+name+"_thermal_forcing_"+str(int(res))+"km_x_60m.nc")



        zbnd = tf.get('z_bnds').values[0,:]
    else:
        tf = xr.open_dataset(
            forcing_data_path+name +"obs_thermal_forcing_1995-2017_"+str(int(res))+"km_x_60m.nc")
        zbnd = tf.get('z_bnds').values
    return(tf,zbnd)
def assure_minimum_same_timelength(mask,bins,base):
    #ensure same time length
    min_length = min(mask.time.size, bins.time.size, base.time.size)
    # Truncate datasets along the time axis if their sizes are different
    if mask.time.size > min_length:
        mask_trun = mask.isel(time=slice(0, min_length))
    else:
        mask_trun = mask

    if bins.time.size > min_length:
        bins_trun = bins.isel(time=slice(0, min_length))
    else:
        bins_trun = bins

    if base.time.size > min_length:
        base_trun = base.isel(time=slice(0, min_length))
    else:
        base_trun = base
    return(mask_trun,bins_trun,base_trun)

def create_empty_dummy_t(variable_name, time_coord):
    # Create an empty dataset with time coordinates
    d_calc = xr.Dataset(coords={'time': time_coord})

    # Get the length of the time coordinate
    time_length = len(time_coord)
    
    # Add the variable without sectors or regions
    d_calc[variable_name] = xr.DataArray(
        data=np.full((time_length,), np.nan),  # Initialize with NaNs
        dims=('time',),
        coords={'time': time_coord}
    )
    
    # Add variables for each region and sector
    regions = range(1, 4)
    sectors = range(1, 19)
    
    for region in regions:
        d_calc[f'{variable_name}_region_{region}'] = xr.DataArray(
            data=np.full((time_length,), np.nan),  # Initialize with NaNs
            dims=('time',),
            coords={'time': time_coord}
        )
    
    for sector in sectors:
        d_calc[f'{variable_name}_sector_{sector}'] = xr.DataArray(
            data=np.full((time_length,), np.nan),  # Initialize with NaNs
            dims=('time',),
            coords={'time': time_coord}
        )
    
    return d_calc
def create_empty_dummy_zt(variable_name,z_values,time_coord):
    # Create z coordinate
    z_coord = xr.DataArray(z_values, dims=('z',), coords={'z': z_values})
# time_coord = xr.DataArray(time_values, dims=('time',), coords={'time': time_values})

    # Create time coordinate
    # time_coord = xr.IndexVariable('time', range(min_length))

    # Create an empty dataset with z and time coordinates
    d_calc = xr.Dataset(coords={'z': z_coord, 'time': time_coord})

    # Add the variable without sectors or regions
    d_calc[variable_name] = xr.DataArray(
        data=None,
        dims=('z', 'time'),
        coords={'z': z_coord, 'time': time_coord}
        )

    # Add variables for each region and sector
    regions = range(1, 4)
    sectors = range(1, 19)

    for region in regions:
        d_calc[f'{variable_name}_region_{region}'] = xr.DataArray(
            data=None,
            dims=('z', 'time'),
            coords={'z': z_coord, 'time': time_coord}
            )

    for sector in sectors:
        d_calc[f'{variable_name}_sector_{sector}'] = xr.DataArray(
            data=None,
            dims=('z', 'time'),
            coords={'z': z_coord, 'time': time_coord}
            )
    return(d_calc)


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
    if not np.issubdtype(d_calc.thermalforcing_interp.values.dtype, np.floating):
        d_calc['thermalforcing_interp'] =d_calc.thermalforcing_interp.astype(float)

    d_calc.thermalforcing_interp.values =  np.full(d['sftflf'].shape, np.nan)   



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

def process_data_mask(d_in, draft_in, zbnd, tf_in, time_model,d_calc,d_bin):
# check time end of forcing file is either 2299 or 2300
    if np.size(tf_in.time.values) > 305:
        time_end = 2300
    else:
        time_end=2299
    if time_model[-1] <time_end:
        time_end = int(time_model[-1])
        printm = 'ohoh the mask is too short'+str(time_end)
    else:
        printm ='all good'
    tf = tf_in.sel(time=slice(str(int(time_model[0])),str(time_end)))
    d = d_in.sel(time=slice(str(int(time_model[0])),str(time_end)))
    draft = draft_in.sel(time=slice(str(int(time_model[0])),str(time_end)))
    m =  d.sftflf.values > 0.5 #take mainly floating cells
    draft.base.values[~m] = 100
    
    shelf_draft_bin = np.empty(draft.base.shape)
    shelf_draft_bin[:] = np.nan

    # This will be the relevant thermal forcing
    shelf_draft_bin_temp = np.empty(draft.base.shape)
    shelf_draft_bin_temp[:] = np.nan
    # create a mask of the depth indices of the shelf draft
    for i, zbnd_i in enumerate(zbnd):
        zmax = zbnd_i[0]
        zmin = zbnd_i[-1]
        ind_mask = (draft.base.values < zmax) & (draft.base.values >= zmin)
        shelf_draft_bin[ind_mask] = i
        masked_values = np.where(ind_mask, tf.thermal_forcing[:,i,].values, np.nan)
        shelf_draft_bin_temp[ind_mask] = masked_values[ind_mask]
    d_bin.bin_number.values[:d.time.size,] = shelf_draft_bin.copy()
    d_calc.thermalforcing_bin.values[:d.time.size,] = shelf_draft_bin_temp.copy()

    return d_calc ,d_bin,printm

def process_data_ctrl_mask(d, draft, zbnd, tf_in, time_model,d_calc,d_bin):
# check time end of forcing file is either 2299 or 2300
    tf_values = tf_in.thermal_forcing.values
# d = d_in.sel(time=slice(str(int(time_model[0])),str(time_end)))
# draft = draft_in.sel(time=slice(str(int(time_model[0])),str(time_end)))
    m =  d.sftflf.values > 0.5 #take mainly floating cells
    draft.base.values[~m] = 100
    
    shelf_draft_bin = np.empty(draft.base.shape)
    shelf_draft_bin[:] = np.nan

    # This will be the relevant thermal forcing
    shelf_draft_bin_temp = np.empty(draft.base.shape)
    shelf_draft_bin_temp[:] = np.nan
    # create a mask of the depth indices of the shelf draft
    for i, zbnd_i in enumerate(zbnd):
        zmax = zbnd_i[0]
        zmin = zbnd_i[-1]
        ind_mask = (draft.base.values < zmax) & (draft.base.values >= zmin)
        shelf_draft_bin[ind_mask] = i
        masked_values = np.where(ind_mask, tf_values[i,], np.nan)
        shelf_draft_bin_temp[ind_mask] = masked_values[ind_mask]
    d_bin.bin_number.values = shelf_draft_bin.copy()
    d_calc.thermalforcing_bin.values = shelf_draft_bin_temp.copy()

    return d_calc ,d_bin
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
    

def process_data_ctrlrun(d, draft, zbnd, tf, time_model, d_calc, d_bin):
    # Precompute the mask for floating cells for all times
    m_floating = d.sftflf.values > 0.5

    # Precompute the thermal forcing values
    tf_values = tf.thermal_forcing.values

    # Iterate through each time step within the model time range
    for timei, time in enumerate(time_model):
        print(time)

        # Get the mask and draft for the current time step
        m = m_floating[timei]
        shelf_draft = np.where(m, draft.isel(time=timei).base.values, 100)

        # Initialize arrays for the shelf draft bin and thermal forcing temperature
        shelf_draft_bin = np.full(shelf_draft.shape, np.nan)
        shelf_draft_bin_temp = np.full(shelf_draft.shape, np.nan)

        # Create a mask of the depth indices of the shelf draft
        for i, (zmax, zmin) in enumerate(zbnd):
            ind_mask = np.logical_and(shelf_draft < zmax, shelf_draft >= zmin)
            shelf_draft_bin[ind_mask] = i
        d_bin.isel(time=timei).bin_number.values[:] = shelf_draft_bin

        # Apply the thermal forcing to the appropriate depth bins
        for i in range(tf_values.shape[0]):
            maski = shelf_draft_bin == i
            shelf_draft_bin_temp[maski] = tf_values[i]

        d_calc.isel(time=timei).thermalforcing_bin.values[:] = shelf_draft_bin_temp

    return d_calc, d_bin


def process_data_withInterp_ctrlrun(d, draft, zbnd, tf, time_model,d_calc):
# check time end of forcing file is either 2299 or 2300
    nz,ny,nx = tf.thermal_forcing.shape
    #no time for tf here!
    for timei, time in enumerate(time_model):
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
        
        # find correct time of thermal forcing
        yr = int(time_model[timei])
        tf_time = tf.thermal_forcing.values
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

        # find the thermal forcing at the correct depth
        for i, z in enumerate(tf.thermal_forcing.z.values):
            maski = shelf_draft_bin == i
            shelf_draft_bin_temp[maski] = tf_time[i, maski]
            shelf_draft_tf_upper[maski] = tf_bnds[i,maski]
            shelf_draft_tf_lower[maski] = tf_bnds[i+1,maski]
        # now interpolate to find the weighted thermal forcing
        thermalforcing_interp_in = weightedInterpTemp(shelf_draft_upper,shelf_draft_lower,shelf_draft,shelf_draft_tf_upper,shelf_draft_tf_lower)
        d_calc.isel(time=timei).thermalforcing_interp.values[:] = thermalforcing_interp_in.copy()

    return d_calc 
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
    #no diving by zeor0
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
    
# function_sample_thermal_forcing.py

def create_loop_info(mnt_pth, expnames, removeFileID, metadata_file):
    pth_ismip6 = os.path.join(mnt_pth, 'ismip6_2300/')
    dirPath = pth_ismip6

    metadata = pd.read_csv(metadata_file)

    models = metadata.loc[(metadata['Experiment'].isin(expnames)) &
                          ~(metadata['fileID'].isin(removeFileID)),]

    # Initialize empty list to store processed data
    loop_info_in = []

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

            # Save into processed_data
            loop_info_in.append({
                'base_file': expFiles[draftPos[0]],
                'mask_file': expFiles[maskPos[0]],
                'path': os.path.join(dirPath, modelPath, expDir),
                'Experiment': exp,
                'Model':  modelPath,
                'Grid': grid
            })

    # Convert to DataFrame
    loop_info = pd.DataFrame(loop_info_in)
    ind= loop_info.index[loop_info.path =='/home/565/jb1863/ismip6_2300/NORCE_CISM5-MAR364-ERA-t1-local/expAE05_16_old'].tolist()
    loop_info = loop_info[loop_info.index != ind[0]].reset_index(drop=True)
    return loop_info
def correct_runs_missing_exp(df):
    df = df[(df.path == '/home/565/jb1863/ismip6_2300/NORCE_CISM5-MAR364-ERA-t1-local/expAE05_16') | (df.Model =='ULB_fETISh-KoriBU2')].reset_index(drop=True)
    return df
