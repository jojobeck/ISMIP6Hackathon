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
renamer = {}
renamer['DOE_MALI_1'] = 'DOE_MALI_4km'
renamer['DOE_MALI_2'] = 'DOE_MALI_8km_Ant95'
renamer['DOE_MALI_3'] = 'DOE_MALI_8km_AntMean'
# Define corridor size
grid =4
glCorridor = 32
glCorridor_str=str(glCorridor)
expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05', 'ctrlAE']
# expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05']

# Define directories
mnt_pth = '/home/565/jb1863/' #mnt/ronja/nci/' #mount path
dataDir = f'{mnt_pth}/ismip6_2300'
dataDirTf = f'{mnt_pth}/ismip6_hackathon/dataset_2HD'
outputDir = f'{mnt_pth}/ismip6_hackathon/ComputedScalars/'
pth_helene =mnt_pth + 'ismip6_hackathon/ComputedScalarsHelene/'
# dataDir = './training/ISMIP6Hackathon'

# Load metadata table
metadata_file ='Metadata.txt'

removeFileID = []

scalefac_pth= f'/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_{grid}km.nc'
#loadscaling mask
scalefac_model= xr.open_dataset(scalefac_pth)
dA= grid*grid*1000*1000

# Only get local models
# TODO: Delete when running on cluster
# models = models.loc[models['fileID'].isin(['PIK_PISM'])]


############################################################
## 1. GET MODEL INFORMATION
############################################################
loop_info = fn.create_loop_info_gl(mnt_pth, expnames, removeFileID, metadata_file)
processedData = loop_info[loop_info['Grid'] == grid]

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
# varNames = ['groundingline_depth', 'groundling_bmb', 'groundingline_thermal_forcing','groundingline_bmb_avg','groudningline_zone_area']
varNames = [ 'groundling_bmb']

## Loop through processed_data and load each mask in
for i in range(processedData.shape[0]):
# for i in range(45,processedData.shape[0]):
    
    dsBmb = xr.Dataset()

    modelname = processedData.iloc[i].Model
    if modelname in ['DOE_MALI_1', 'DOE_MALI_2', 'DOE_MALI_3']:
        model = renamer[modelname]
    else:
        model= modelname
    exp = processedData.iloc[i].Experiment
    grid = processedData.iloc[i].Grid

    writeLoc = outputDir + exp
    fnameBmb = 'computed_groundingline_bmb_corr_'+glCorridor_str+'_AIS_' + model + '_' + exp + '.nc'

    maskData_float = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].maskfile)
    if model.startswith('IMAU'):
        maskfile_ice=  processedData.iloc[i].maskfile.replace('new_sftflf','sftgif')
        maskData_grounded = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].maskfile.replace('new_sftflf','sftgrf'))
    else:
        maskfile_ice=  processedData.iloc[i].maskfile.replace('sftflf','sftgif')
        maskData_grounded = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].maskfile.replace('sftflf','sftgrf'))
    mask_ice = xr.open_dataset(processedData.iloc[i].path+'/'+maskfile_ice)
    BmbData = xr.open_dataset(processedData.iloc[i].path + '/' + processedData.iloc[i].libmassbfflfile)
    filename =  pth_helene+'/'+exp+'/shelfmelt/computed_shelfmelt_AIS_'+model+'_'+exp+'.nc'
    d_example = xr.open_dataset(filename,decode_times=False )

    print('opening '+processedData.iloc[i].thermalforcingfile)

    # remove any values after 2300-01-01
    # TODO: Add tfData into this list when it's processed.

    glBmb = BmbData.libmassbffl.copy()  # Pre-assign output size
    nT = BmbData.time.size

    # Save time
    # TODO: Add dsTf into this list when it's processed
    for ds in (dsBmb):
        ds.coords['time'] =d_example.time[:nT]

    # How many cells should be used for the GL corridor
    nCells = glCorridor // grid

    # Loop over each timestep, isolate the mask, sample parameters
    bmbTs = np.zeros(nT)

    # Antarctica wide
    print('-    Compute ANT-wide...')
    for j in range(nT):
        mask = maskData_grounded.sftgrf.isel(time=j)
        mask_fl = maskData_float.sftflf.isel(time=j)
        m_ice = mask_ice.sftgif.isel(time=j)
        area_dis=scalefac_model.af2.values
        #adjust similar to Helene
        #mask handling as Helene
        m_ice.values[m_ice.values<0]=0
        m_ice.values[m_ice.values>1]=1
        mask.values[mask.values<0]=0
        mask.values[mask.values>1]=1
        mask_fl.values[mask_fl.values<0]=0
        mask_fl.values[mask_fl.values>1]=1
        maskT = getGLMask(mask,mask_fl, nCells)
        glBmb[j, :, :] = xr.where(maskT == 1, BmbData.libmassbffl.isel(time=j), np.nan)

        bmbTs[j] = float(np.nansum(glBmb.isel(time=j)*mask_fl*m_ice*dA*area_dis))

    # Save Antarctica wide
    print('-    Save ANT-wide...')
    dsBmb['groundingline_bmb'] = xr.DataArray(bmbTs, dims=['time'])
    
    # Process Sectors
    print('-    Compute Sectors...')
    sectFrame = sectors[grid]
    for sec in secNums:
        secMask = sectFrame.sectors.values == sec
        secBmb = np.where(secMask, glBmb.values, np.nan)

        bmb = np.zeros(nT)
        for j in range(nT):
            mask_fl = maskData_float.sftflf.isel(time=j).values
            m_ice = mask_ice.sftgif.isel(time=j).values
            area_dis=scalefac_model.af2.values
            #adjust similar to Helene
            #mask handling as Helene
            m_ice[m_ice<0]=0
            m_ice[m_ice>1]=1
            mask_fl[mask_fl<0]=0
            mask_fl[mask_fl>1]=1

            bmb[j] = float(np.nansum(secBmb[j,]*mask_fl*m_ice*dA*area_dis))
        print('-    Save Sector ' + str(sec))
        dsBmb[f'groundingline_bmb_sector_{sec}'] = xr.DataArray(bmb, dims=['time'])

    # Process Regions
    print('-    Compute Regions...')
    for reg in regNums:
        regMask = sectFrame.regions.values == reg
        regBmb = np.where(regMask, glBmb.values, np.nan)

        bmb = np.zeros(nT)
        for j in range(nT):
            mask_fl = maskData_float.sftflf.isel(time=j).values
            m_ice = mask_ice.sftgif.isel(time=j).values
            area_dis=scalefac_model.af2.values
            #adjust similar to Helene
            #mask handling as Helene
            m_ice[m_ice<0]=0
            m_ice[m_ice>1]=1
            mask_fl[mask_fl<0]=0
            mask_fl[mask_fl>1]=1

            bmb[j] = float(np.nansum(secBmb[j,]*mask_fl*m_ice*dA*area_dis))

        print('-    Save Region ' + str(reg))
        dsBmb[f'groundingline_bmb_region_{reg}'] = xr.DataArray(bmb, dims=['time'])
    # Add attributes
    # TODO: Add dsTf into this list when it's processed
    print('Add attributes')
# for ds in (dsBmb):
    dsBmb.attrs['model'] = model
    dsBmb.attrs['experiment'] = exp
    dsBmb.attrs['grid'] = grid

    # Save to file
    print('Save to file')
    saveDir = outputDir + exp + '/groundingline_corridor'+glCorridor_str
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dsBmb.to_netcdf(saveDir + '/' + fnameBmb)
    print(model + ' Complete...')

