#!/usr/bin/env python
# coding: utf-8

# # Calculate the grounded area wrt control from Helene's computed scalars

# In[1]:


import os
import re
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as n

import xarray as xr
import seaborn as sns
import matplotlib.lines as mlines
import pandas as pd

sns.set_context("paper")
sns.set_style("whitegrid")

computed_name = 'ComputedScalarsHelene'


# In[2]:


import function_sampling_Thermal_forcing_files as fn # FIXME needed an update, unfortunately


# In[3]:


mnt_path = '/mnt/ronja/nci/'

path = '/mnt/ronja/nci/ismip6_hackathon/'

pth_calc_output = '/mnt/ronja/nci/ismip6_hackathon/ComputedScalars/'


# In[4]:


df = pd.read_csv('Metadata.txt') 


# In[5]:


expnames = ['expAE02', 'expAE03', 'expAE04', 'expAE05',]

removeFileID = ['IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
               'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean'] # specify models to remove

metadata_file = 'Metadata.txt'

loop_info = fn.create_loop_info(mnt_path, expnames, removeFileID, metadata_file)


# In[6]:


loop_info


# In[37]:


#loop_info = correct_runs_missing_exp(loop_info, mnt_path)  # selects only Kori runs


# In[9]:


yearlen = 360*24*60*60 # FIXME: is this the same for all models? if not, does it matter?


# this is the name of the variable to be generated
variable_name = 'iareagr_anom'


regions =['region_1','region_2','region_3']
sectors =[]
for i in range(18):
    sectors.append('sector_'+str(i+1))


for fi in range(len(loop_info.index[:1])): # FIXME remove [:1]
    print (fi)
    
    pth = loop_info['path'][fi]
    #grid = loop_info['Grid'][fi]
    expname = loop_info['Experiment'][fi]
    model = loop_info['Model'][fi]
    #maskfile = loop_info['mask_file'][fi]
    #basefile = loop_info['base_file'][fi]
    
    print(pth)
    filename = os.path.join(path,computed_name,expname,'iareagr','computed_iareagr_AIS_'+model+'_'+expname+'.nc')
    d_iareagr = xr.open_dataset(filename,decode_times=False )
    d_iareagr_ctrl = xr.open_dataset(filename.replace(expname, 'ctrlAE'),decode_times=False)
    
    
    # create a new xarray
    d_new =d_iareagr.copy()
    d_new.attrs ={}
    #rename
    d_new =d_new.rename({'iareagr':variable_name })
    for region in regions:
        varname_old = 'iareagr'+'_'+region
        varname_new = variable_name+'_'+region    
        d_new =d_new.rename({varname_old:varname_new })
    for sector in sectors:
        varname_old = 'iareagr'+'_'+sector
        varname_new = variable_name+'_'+sector 
        d_new =d_new.rename({varname_old:varname_new })
   
    # Loop through each variable in the dataset and add values
    for var_name in d_new.variables:
        if (var_name != 'time') and (var_name != 'rhoi') and (var_name != 'rhow'):
                
            var_name_iareagr = var_name.replace(variable_name, 'iareagr')
            data_iareagr = d_iareagr[var_name_iareagr].values # m2?
            data_iareagr_ctrl = d_iareagr_ctrl[var_name_iareagr].values # m2
                
            dt = 1 # FIXME we assume a yearly time step
            d_new[var_name].values =  np.squeeze(data_iareagr - data_iareagr_ctrl) 
        
        
    # save dataset 
    save_name = filename.split('/')[-1].split('/')[-1].replace('iareagr', variable_name)
    path_save = os.path.join(pth_calc_output,expname,variable_name)
        
    if not os.path.exists(os.path.join(pth_calc_output,expname)):
        os.mkdir( os.path.join(pth_calc_output,expname) )
    if not os.path.exists(path_save):
        os.mkdir(path_save)
            
        
    d_new.to_netcdf(path_save+'/'+save_name)    


# In[ ]:




