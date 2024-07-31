'''
This script samples ISMIP6 data from close to the grounding line. It assumes a fixed corridor width (default = 8 km), sampling all shelf cells within this corridor.

Local relative paths assume that data exist at the following location: XXXX

Major and minor section breaks are denoted using the following line breaks, respectively:
############################################################
## ---------------------------------------------------------

REQUIREMENTS:
pandas (v xxx)
numpy (v xxx)
....

'''

import numpy as np
import xarray as xr
from matplotlib import pylab as plt
import pandas as pd
import os
import function_sampling_Thermal_forcing_files as fn

############################################################
## DEFINE NECESSARY FUNCTIONS
############################################################
def getBoundaryCells(mask):
    '''
    Identify border cells around a given ice mask, based on the first floating cell or open ocean.

    :param mask: xarray binary mask of ice/no-ice
    '''

    borderMaskFloat = np.zeros(mask.shape)
    difM = mask[1:, :] - mask[:-1, :]
    m1 = difM == 1
    mMinus1 = difM == -1
    borderMaskFloat[:-1, ][m1] = 1
    borderMaskFloat[1:, ][mMinus1] = 1

    difU = mask[:, 1:] - mask[:, :-1]
    m1u = difU == 1
    m_Minus1u = difU == -1
    borderMaskFloat[:, :-1][m1u] = 1
    borderMaskFloat[:, 1:][m_Minus1u] = 1

    return borderMaskFloat


def getGLMask(mask,mask_fl, nCells = 1):
    '''
    Get grounding line cells, using the getBoundaryCells function above.
    flaoting and grouned mask are needed at all the models
    don't have consitentgly handling with 0 and nans.

    :param mask: xarray binary mask of grounded ice
    :param mask_fl: xarray binary mask of floating ice
    :param nCells: number of cells to buffer the grounding line by (default = 1)
    '''

    m = np.zeros(mask.values.shape)
    m[mask.values == 1] = 1
    newM = m.copy()

    for i in range(nCells):
        boundaryCells = getBoundaryCells(newM)
        newM[boundaryCells == 1] = 1

    shelfGL = (newM == 1) & (mask_fl.values > 0)

    return shelfGL

############################################################
## 0. LOAD MODEL METADATA & FILTER MODELS
############################################################

# Define corridor size
glCorridor = 32
glCorridor_str=str(glCorridor)
expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']

# Define directories
mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
dataDir = f'{mnt_pth}/ismip6_2300'
dataDirTf = f'{mnt_pth}/ismip6_hackathon/dataset_2HD'
outputDir = f'{mnt_pth}/ismip6_hackathon/ComputedScalars/'
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'
# dataDir = './training/ISMIP6Hackathon'

# Load metadata table
metadata_file ='Metadata.txt'

removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove


# Only get local models
# TODO: Delete when running on cluster
# models = models.loc[models['fileID'].isin(['PIK_PISM'])]


############################################################
## 1. GET MODEL INFORMATION
############################################################
loop_info = fn.create_loop_info_gl(mnt_pth, expnames, removeFileID, metadata_file)
processedData = loop_info[loop_info['Grid'] == 4]

############################################################
## 2. LOAD SECTOR INFORMATION
############################################################
print('Loading sector information')
sectors = {}
for grid in processedData.Grid.unique():
    sectors[grid] = xr.open_dataset(dataDir + '/masks/sectors_' + str(grid) + 'km.nc')

# Get sector and region numbers
regNums = np.unique(sectors[list(sectors.keys())[0]].regions).astype(int)
secNums = np.unique(sectors[list(sectors.keys())[0]].sectors).astype(int)

############################################################
## 3. PROCESS ALL VARIABLES
############################################################
print('Processing variables...')
varNames = ['groundingline_depth', 'groundling_bmb', 'groundingline_thermal_forcing','groundingline_bmb_avg','groudningline_zone_area']

## Loop through processed_data and load each mask in
for i in range(processedData.shape[0]):
    print(i)
    dsBmb = xr.Dataset()
    dsAvgBmb = xr.Dataset()
    dsDraft = xr.Dataset()
    dsTf = xr.Dataset()
    dsZone = xr.Dataset()

    model = processedData.iloc[i].Model
    exp = processedData.iloc[i].Experiment
    grid = processedData.iloc[i].Grid
    print(model)

    writeLoc = outputDir + exp
    fnameBmb = 'computed_groundingline_bmb_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'
    fnameBmb_avg = 'computed_groundingline_bmb_avg_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'
    fnameDraft = 'computed_groundingline_depth_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'
    fnameZone = 'computed_groundingline_zone_area_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'
    fnameTf = 'computed_groundingline_thermal_forcing_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'

    maskData_float = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].maskfile)
    maskData_grounded = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].maskfile.replace('sftflf','sftgrf'))
    BmbData = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].libmassbfflfile)
    draftData = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].basefile)
    tfData = xr.open_dataset(processedData.iloc[i].pathTf + '/' + processedData.iloc[i].thermalforcingfile)
    filename =  pth_helene+'/'+exp+'/shelfmelt/computed_shelfmelt_AIS_'+model+'_'+exp+'.nc'
    d_example = xr.open_dataset(filename,decode_times=False )

    print('opening '+processedData.iloc[i].thermalforcingfile)

    # remove any values after 2300-01-01
    # TODO: Add tfData into this list when it's processed.

    glBmb = BmbData.libmassbffl.copy()  # Pre-assign output size
    glDraft = draftData.base.copy()  # Pre-assign output size
    glTf = tfData.thermalforcing_interp.copy()  # Pre-assign output size
    minimum_dataset=fn.min_time_size_dataset(maskData_grounded,BmbData,draftData,tfData)
    nT = minimum_dataset.time.size

    # Save time
    # TODO: Add dsTf into this list when it's processed
    for ds in (dsBmb, dsAvgBmb,dsZone,dsDraft, dsTf):
        ds.coords['time'] =d_example.time[:nT]

    # How many cells should be used for the GL corridor
    nCells = glCorridor // grid

    # Loop over each timestep, isolate the mask, sample parameters
    glDraftValues = draftData.base.values
    bmbTs = np.zeros(nT)
    avg_bmbTs = np.zeros(nT)
    draftTs = np.zeros(nT)
    zoneTs = np.zeros(nT)
    tfTs = np.zeros(nT)

    # Antarctica wide
    print('-    Compute ANT-wide...')
    for j in range(nT):
        mask = maskData_grounded.sftgrf.isel(time=j)
        mask_fl = maskData_float.sftflf.isel(time=j)
        maskT = getGLMask(mask,mask_fl, nCells)
        glBmb[j, :, :] = xr.where(maskT == 1, BmbData.libmassbffl.isel(time=j), np.nan)
        glDraft[j, :, :] = xr.where(maskT == 1, glDraftValues[j], np.nan)
        glTf[j, :, :] = xr.where(maskT == 1, tfData.thermalforcing_interp.isel(time=j), np.nan)

        bmbTs[j] = float(np.nansum(glBmb.isel(time=j)))
        avg_bmbTs[j] = float(np.nanmean(glBmb.isel(time=j)))
        nonans = ~np.isnan(glBmb.isel(time=j) )
        zoneTs[j] =nonans.sum()*grid

        draftTs[j] = float(np.nanmean(glDraft.isel(time=j)))
        tfTs[j] = float(np.nanmean(glTf.isel(time=j)))

    # Save Antarctica wide
    print('-    Save ANT-wide...')
    dsBmb['groundingline_bmb'] = xr.DataArray(bmbTs, dims=['time'])
    dsAvgBmb['groundingline_bmb_avg'] = xr.DataArray(bmbTs, dims=['time'])
    dsZone['groundingline_zone_area'] = xr.DataArray(zoneTs, dims=['time'])
    dsDraft['groundingline_depth'] = xr.DataArray(draftTs, dims=['time'])
    dsTf['groundingline_thermal_forcing'] = xr.DataArray(tfTs, dims=['time'])
    
    # Process Sectors
    print('-    Compute Sectors...')
    sectFrame = sectors[grid]
    for sec in secNums:
        secMask = sectFrame.sectors.values == sec
        bmb = np.zeros(nT)
        draft = np.zeros(nT)
        tf = np.zeros(nT)
        avg_bmb = np.zeros(nT)
        zone = np.zeros(nT)
        secBmb = np.where(secMask, glBmb.values, np.nan)
        secDraft = np.where(secMask, glDraft.values, np.nan)
        secTf = np.where(secMask, glTf.values, np.nan)
        for j in range(nT):


            nonans = ~np.isnan(secBmb[j,] )

            bmb[j] = float(np.nansum(secBmb[j,]))
            avg_bmb[j] = float(np.nanmean(secBmb[j,]))
            draft[j] = float(np.nanmean(secDraft[j,]))
            tf[j] = float(np.nanmean(secTf[j,]))
            zone[j] = float(np.sum(nonans)*grid)
        print('-    Save Sector ' + str(sec))
        dsBmb[f'groundingline_bmb_sector_{sec}'] = xr.DataArray(bmb, dims=['time'])
        dsAvgBmb[f'groundingline_bmb_avg_sector_{sec}'] = xr.DataArray(avg_bmb, dims=['time'])
        dsDraft[f'groundingline_depth_sector_{sec}'] = xr.DataArray(draft, dims=['time'])
        dsTf[f'groundingline_thermal_forcing_sector_{sec}'] = xr.DataArray(tf, dims=['time'])
        dsZone[f'groundingline_zone_area_sector_{sec}'] = xr.DataArray(zone, dims=['time'])


    # Process Regions
    print('-    Compute Regions...')
    for reg in regNums:
        regMask = sectFrame.regions.values  == reg
        regBmb = np.where(regMask, glBmb.values, np.nan)
        regDraft = np.where(regMask, glDraft.values, np.nan)
        regTf = np.where(regMask, glTf.values, np.nan)

        bmb = np.zeros(nT)
        draft = np.zeros(nT)
        tf = np.zeros(nT)
        avg_bmb = np.zeros(nT)
        zone = np.zeros(nT)
        for j in range(nT):
            nonans = ~np.isnan(regBmb[j,] )

            bmb[j] = float(np.nansum(regBmb[j,]))
            avg_bmb[j] = float(np.nanmean(regBmb[j,]))
            draft[j] = float(np.nanmean(regDraft[j,]))
            tf[j] = float(np.nanmean(regTf[j,]))
            zone[j] = float(np.sum(nonans)*grid)

        print('-    Save Region ' + str(reg))
        dsBmb[f'groundingline_bmb_region_{reg}'] = xr.DataArray(bmb, dims=['time'])
        dsDraft[f'groundingline_depth_region_{reg}'] = xr.DataArray(draft, dims=['time'])
        dsTf[f'groundingline_thermal_forcing_region_{reg}'] = xr.DataArray(tf, dims=['time'])
        dsAvgBmb[f'groundingline_bmb_avg_region_{reg}'] = xr.DataArray(avg_bmb, dims=['time'])
        dsZone[f'groundingline_zone_area_region_{reg}'] = xr.DataArray(zone, dims=['time'])

    # Add attributes
    # TODO: Add dsTf into this list when it's processed
    print('Add attributes')
    for ds in (dsBmb, dsDraft, dsTf,dsAvgBmb,dsZone):
        ds.attrs['model'] = model
        ds.attrs['experiment'] = exp
        ds.attrs['grid'] = grid

    # Save to file
    print('Save to file')
    saveDir = outputDir + exp + '/groundingline_corridor'+glCorridor_str
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dsBmb.to_netcdf(saveDir + '/' + fnameBmb)
    dsAvgBmb.to_netcdf(saveDir + '/' + fnameBmb_avg)
    dsZone.to_netcdf(saveDir + '/' + fnameZone)
    dsDraft.to_netcdf(saveDir + '/' + fnameDraft)
    dsTf.to_netcdf(saveDir + '/' + fnameTf)
    print(model + ' Complete...')

