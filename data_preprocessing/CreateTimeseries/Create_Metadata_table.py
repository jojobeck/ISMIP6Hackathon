#!/usr/bin/env python
# coding: utf-8

# # Create pandas dataframe containing relevant ISMIP6 2300 infos
# 
# This notebook creates a pandas Dataframe containing the relevant metadata for the ISMIP6 AIS 2300 simulations

# In[1]:


import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


pd.set_option('display.max_columns', None)  
pd.set_option('display.max_rows', None)  
pd.set_option('display.max_colwidth', None)  


# In[3]:


# We here collect all information

# List of the experiments we will analyse
experiments = ['ctrlAE','expAE02', 'expAE03', 'expAE04', 'expAE05']
forcing_models = ['ctrl','CCSM4', 'HadGEM2', 'CESM2', 'UKESM']

## TODO: add here a list of the correspondong forcing names
#experiment_models = ['ukesm1-0-ll_ssp585']
#experiment_model_names = ['UKESM1-0-LL_SSP585'] 

# Groups that submitted runs and the models they used
groups =['DC', 'DOE', 'IGE', 'ILTS', 'IMAU', 'LSCE', 'NCAR', 'NORCE', 'PIK', 'UGM', 'UCSD', 'ULB',
        'UNN', 'UTAS', 'VUB', 'VUW']
models = ['ISSM', 'MALI', 'Elmer/Ice', 'SICOPOLIS', 'UFEMISM', 'Grisli', 'CISM', 'CISM', 'PISM', 'Yelmo',
         'ISSM', 'Kori', 'Úa', 'Elmer/Ice', 'AISMPALEO', 'PISM']

# number of submissions per group
n_submissions = [1, 3, 1, 1, 4, 2, 2, 11,1, 1, 1, 2, 1, 1, 1, 10 ]

# names of all submissions for plotting purposes
names_of_submissions = ['DC\_ISSM', 
                        'DOE\_MALI\_4km', 'DOE\_MALI\_8km\_Ant95', 'DOE\_MALI\_8km\_AntMean',
                        'IGE\_ElmerIce',
                        'ILTS\_SICOPOLIS',
                        'IMAU\_UFEMISM1','IMAU\_UFEMISM2', 'IMAU\_UFEMISM3', 'IMAU\_UFEMISM4',
                        'LSCE\_GRISLI','LSCE\_GRISLI2',
                        'NCAR\_CISM1', 'NCAR\_CISM2*',
                        'NORCE\_CISM2', 'NORCE\_CISM3','NORCE\_CISM3\_nonlocal*','NORCE\_CISM3\_local',
                        'NORCE\_CISM4', 'NORCE\_CISM4\_nonlocal', 'NORCE\_CISM4\_local', 'NORCE\_CISM4\_JRA',
                        'NORCE\_CISM5', 'NORCE\_CISM5\_nonlocal', 'NORCE\_CISM5\_local',
                        'PIK_PISM',
                        'UCM\_Yelmo',
                        'UCSD\_ISSM',
                        'ULB\_Kori1', 'ULB\_Kori2',
                        'UNN\_\'Ua',
                        'UTAS\_ElmerIce',
                        'VUB\_AISMPALEO',
                        'VUW\_PISM1','VUW\_PISM1\_s1', 'VUW\_PISM1\_s2','VUW\_PISM1\_s3', 'VUW\_PISM1\_s4',
                        'VUW\_PISM2', 'VUW\_PISM2\_s1', 'VUW\_PISM2\_s2', 'VUW\_PISM2\_s3', 'VUW\_PISM2\_s4'                      
                       ]

# names of all submissions for accessing purposes
names_of_files = ['DC_ISSM', 
                        'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean',
                        'IGE_ElmerIce',
                        'ILTS_SICOPOLIS',
                        'IMAU_UFEMISM1','IMAU_UFEMISM2', 'IMAU_UFEMISM3', 'IMAU_UFEMISM4',
                        'LSCE_GRISLI','LSCE_GRISLI2',
                        'NCAR_CISM1', 'NCAR_CISM2',
                        'NORCE_CISM2-MAR364-ERA-t1', 'NORCE_CISM3-MAR364-ERA-t1','NORCE_CISM3-MAR364-ERA-t1-nonlocal', 'NORCE_CISM3-MAR364-ERA-t1-local', 
                        'NORCE_CISM4-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-ERA-t1-nonlocal', 'NORCE_CISM4-MAR364-ERA-t1-local', 'NORCE_CISM4-MAR364-JRA-t1',
                        'NORCE_CISM5-MAR364-ERA-t1', 'NORCE_CISM5-MAR364-ERA-t1-nonlocal', 'NORCE_CISM5-MAR364-ERA-t1-local',
                        'PIK_PISM',
                        'UCM_Yelmo',
                        'UCSD_ISSM',
                        'ULB_fETISh-KoriBU1', 'ULB_fETISh-KoriBU1',
                        'UNN_Ua',
                        'UTAS_ElmerIce',
                        'VUB_AISMPALEO',
                        'VUW_PISM1','VUW_PISM1_s1', 'VUW_PISM1_s2','VUW_PISM1_s3', 'VUW_PISM1_s4',
                        'VUW_PISM2', 'VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4'                      
                       ]


# Groups with several submissions did indicate their main submission
main_submission = [
    True,
    True, False, False, 
    True,
    True, 
    True, False, False, False, 
    True, False, 
    False, True, 
    False, False, True, False,
    False, False, False, False,
    False, False, False,
    True,
    True,
    True,
    True, False, 
    True, 
    True, 
    True,
    True, False, False, False, False,
    False, False, False, False, False, 
    
]


# The overall grid resolution, if unstructured I chose the lowest one
resolution = ['1km', '2km','4km', '8km', '16km', '32km']
resolution_per_submission = [
                        2,
                        2,3,3,
                        0, 
                        3,
                        5,5,4,4,
                        4,4,
                        2,2,
                        2, 3, 3, 3, 
                        4, 4, 4, 4, 
                        5, 5, 5,
                        3, 
                        4,
                        0,
                        4,
                        4,
                        0,
                        0,
                        4,
                        3, 3, 3, 3, 3,
                        3, 3, 3, 3, 3            
]


# The grounding line grid resolution, do not account for sub-grid friction treatmeat here 
resolution_gl = ['0.5km', '1km','2km','4km', '8km', '16km']
resolution_gl_per_submission = [
                        3,
                        2,2,2,
                        0, 
                        4,
                        5,5,4,4,
                        5,5,
                        2,2,
                        2,3, 3, 3, 
                        4, 4, 4, 4, 
                        5, 5, 5,
                        3, 
                        4,
                        0,
                        5,
                        5,
                        0,
                        0,
                        5,
                        3, 3, 3, 3, 3,
                        3, 3, 3, 3, 3         
]


# This is a list of melt parameterisations the groups did use
melt_parameterisations = ['Quad. Non-local', 'PICO', 'Quad. Local', 
                         'Quad. Non-local Slope','PICOP', 'Lin']
melt_parameterisation_per_submission = [
                        0,
                        0,0,0,
                        1, 
                        0,
                        2,2,2,2,
                        0,0,
                        0,2,
                        3,3, 0, 2, 
                        3, 0, 2, 3, 
                        3, 0, 2,
                        1, 
                        0,
                        4,
                        1,
                        0,
                        2,
                        2,
                        0,
                        5, 5, 5, 5, 5,
                        5, 5, 5, 5, 5                   
                 ]


# Parameters used in the melt parameterisation
melt_parameterisation_pamas = [ 'AntMean median', 'AntMean 5th', 'AntMean 95th',
                                'PIGL median', 'PIGL 5th', 'PIGL 95th',
                                'custom' ]

melt_parameterisation_params_per_submission = [
                            0,
                            0,2,0,
                            -1,
                            0,
                            -1,-1,-1,-1,
                            0,0,
                            0,0,
                            0,0,0,0,
                            0,0,0,0,
                            0,0,0,
                            -1,
                            0,
                            -1,
                            -1,0,
                            0,
                            0,
                            -1,
                            -1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,               
                 ]


# Treatment of the grounding line - subgrid interpolation etc, probabaly not needed
grounding_line_treatment = ['Sub-Grid', 'Floating condition', 'No', 'N/A', 'Yes']
grounding_line_treatment_per_submission = [
    0,
    1,1,1,
    2,
    1,
    1,1,1,1,
    1,2,
    0,0,
    1,1,1,1,
    1,1,1,1,
    1,1,1,
    1,
    0,
    0,
    2,2,
    2,
    0,
    3,
    2,2,2,2,2,
    4,4,4,4,4,
]


# Treatment of melting across the grounding line
grounding_line_melting = ['Sub-Grid', 'No melt', 'N/A']
grounding_line_melting_per_submission = [
    0,
    1,1,1,
    1,
    1,
    1,1,1,1,
    0,1,
    0,0,
    0,0,0,0,
    0,0,0,0,
    0,0,0,
    1,
    0,
    1,
    1,1,
    1,
    1,
    1,
    1,1,1,1,1,
    0,0,0,0,0,    
]


# In[4]:


# We put all information into one Dataframe

MD = pd.DataFrame({})

tot_count = 0

for i,group in enumerate(groups):
    print(group)
    model = models[i]
    for j in range(n_submissions[i]):
        #print(j)
        for k,exp in enumerate(experiments):
            #print(exp)
            temp_dict ={'Group': group, 
                        'Model': model, 
                        'Experiment': exp,
                        'Climate Forcing': forcing_models[k],
                        'Model setup':names_of_submissions[tot_count],
                        'fileID':names_of_files[tot_count],
                        'Main submission': main_submission[tot_count],
                        'Melt parameterisation': melt_parameterisations[melt_parameterisation_per_submission[tot_count]],
                        'Melt parameters':melt_parameterisation_pamas[melt_parameterisation_params_per_submission[tot_count]],
                        'Resolution':resolution[resolution_per_submission[tot_count]],
                        'Resolution GL': resolution_gl[resolution_gl_per_submission[tot_count]],
                        'GL treatment': grounding_line_treatment[grounding_line_treatment_per_submission[tot_count]],
                        'GL melt': grounding_line_melting[grounding_line_melting_per_submission[tot_count]],                       
                       }
            MD = MD.append(temp_dict, ignore_index=True)
        tot_count = tot_count+1


# In[5]:


MD


# In[6]:


4*43 # This is how long it should be


# In[7]:


# Save the dataframe

MD.to_csv('Metadata.txt')


# # Examples on how to use the dataframe

# In[8]:


MD = pd.read_csv('Metadata.txt', index_col=0)


# In[9]:


display(MD)


# In[12]:


MD.loc[MD['Experiment']==exp,'fileID']


# In[13]:


MD['fileID']


# In[14]:


# test loading data using the names 

yearlen = 360*24*60*60

exp='expAE02'
var = 'shelfmelt'

MD_exp = MD.loc[MD['Experiment']==exp,:]

path_to_files = '/home/ronja/projects/MeltSensitivity/10528582/ComputedScalars/ComputedScalars/'+exp+'/'+var+'/'

plt.figure(figsize=[20,10])

for i in MD_exp.index:
    filename = path_to_files+'computed_'+var+'_AIS_'+MD['fileID'][i]+'_'+exp+'.nc'
    #print(filename)
    ncf = nc.Dataset(filename)
    melt = np.squeeze(ncf['shelfmelt']) *yearlen/1e12*-1 # kg/s -> Gt/a, make it positive
    
    plt.plot(melt)


# In[ ]:




