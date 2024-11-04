import numpy as np 
import xarray as xr
imaus = ['IMAU_UFEMISM1','IMAU_UFEMISM2','IMAU_UFEMISM3','IMAU_UFEMISM4']
resolutions = ['32','32','16','16']
exps = ['expAE02','expAE03','expAE04','expAE05','ctrlAE']
# imaus = ['IMAU_UFEMISM1']
# resolutions = ['32']
# exps = ['expAE02']
for j,IMAU in enumerate(imaus):
    res = resolutions[j]
    for exp in exps:
        pth = f'/home/565/jb1863/ismip6_2300/{IMAU}/{exp}_{res}/'
        rho_ice = 910
        rho_ocean = 1028
        dic={}
        dic['floating_ice_mask']=f'sftflf_AIS_{IMAU}_{exp}.nc'
        dic['ice_mask']=f'sftgif_AIS_{IMAU}_{exp}.nc'
        dic['thk']=f'lithk_AIS_{IMAU}_{exp}.nc'
        dic['bed']=f'topg_AIS_{IMAU}_{exp}.nc'
        #create thick calc 
        d_icemask= xr.open_dataset(pth +dic['ice_mask'])
        d_thk = xr.open_dataset(pth +dic['thk'])
        d_bed = xr.open_dataset(pth + dic['bed'])
        d_maskfloat  =xr.open_dataset(pth+ dic['floating_ice_mask'])

        #make copy of original false netcdf data for safety
        d_maskfloat.to_netcdf(pth + 'orig_'+dic['floating_ice_mask'])
        # calc floating base from thk
        base_float = d_thk['lithk'].values *-1*rho_ice/rho_ocean
        m_float = base_float.copy() *0
        m_float[base_float>d_bed['topg'].values] =1
        m_float[d_icemask['sftgif'].values !=1]=0
        d_maskfloat['sftflf'].values = m_float
        d_maskfloat.to_netcdf(pth +'new_'+dic['floating_ice_mask'])

