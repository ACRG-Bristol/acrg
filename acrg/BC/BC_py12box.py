import pandas as pd
import numpy as np
import glob
import xarray as xr
from datetime import datetime as dt
import getpass

def make_BCs_from_py12box(species, pathtobm, years, domain, outdir, speciestocopy='cfc-12', var_name='Semihemispheric_modelled_mole_fractions'):
    """
    Function to create boundary conditions from py12box data.
    
    Args:
    species (str): Chemical species name or formula. This should match the py12box naming convention, e.g., 'HFC-23', 'HCFC-141b', 'CCl4'.
    pathtobm (str): Path to the py12box_agage repository. This can be cloned from https://github.com/mrghg/py12box_agage
    years (np.ndarray): Array of years to process.
    domain (str): Domain for the boundary conditions (e.g. EUROPE). This is suggested to be a domain from openghg_defs. 
    outdir (str): Output directory where the boundary condition files will be saved. Should already have a {domain} subdirectory.
    speciestocopy (str): Species to copy data for. Defaults to 'cfc-12'. 
    var_name (str): Variable name to copy from py12box_agage outputs. Defaults to 'Semihemispheric_modelled_mole_fractions'.
    
    Returns:
    None
    """
    # Load the CSV file from the box model outputs
    csv = pd.read_csv(pathtobm+f"{species}/outputs/{species}_{var_name}.csv", comment="#")

    # Load in a template file 
    # 
    fn = sorted(glob.glob('/group/chem/acrg/LPDM/bc/'+domain+'/'+speciestocopy+'_'+domain+'_'+str(201701)+'*'))
    with xr.open_dataset(fn[0]) as load:
        bcds = load.load()
    lons = bcds.lon.values
    lats = bcds.lat.values
    minlat =  np.min(lats)
    maxlat = np.max(lats) 

    for year in years:
        yind = np.floor(csv["Year"].values) == float(year)

        for yx in range(12):

            if minlat <= -30:
                bcds.vmr_s.values[:] = csv[f"{var_name}_box3"][yind].values[yx]*1e-12
            elif minlat > -30 and minlat <= 0:
                bcds.vmr_s.values[:] = csv[f"{var_name}_box2"][yind].values[yx]*1e-12
            elif minlat > 0 and minlat <= 30:
                bcds.vmr_s.values[:] = csv[f"{var_name}_box1"][yind].values[yx]*1e-12
            elif minlat >= 30:
                bcds.vmr_s.values[:] = csv[f"{var_name}_box0"][yind].values[yx]*1e-12
            
            # bcds.vmr_s.values[:] = vals
            
            if maxlat <= -30:
                bcds.vmr_n.values[:] = csv[f"{var_name}_box3"][yind].values[yx]*1e-12
            elif maxlat > -30 and maxlat <= 0:
                bcds.vmr_n.values[:] = csv[f"{var_name}_box2"][yind].values[yx]*1e-12
            elif maxlat > 0 and maxlat <= 30:
                bcds.vmr_n.values[:] = csv[f"{var_name}_box1"][yind].values[yx]*1e-12
            elif maxlat >= 30:
                bcds.vmr_n.values[:] = csv[f"{var_name}_box0"][yind].values[yx]*1e-12
        
            bcds.vmr_e.values[:,np.where(lats <= -30),:] = csv[f"{var_name}_box3"][yind].values[yx]*1e-12
            bcds.vmr_e.values[:,np.where(np.logical_and(lats > -30, lats <= 0)),:] = csv[f"{var_name}_box2"][yind].values[yx]*1e-12
            bcds.vmr_e.values[:,np.where(np.logical_and(lats > 0, lats <= 30)),:] = csv[f"{var_name}_box1"][yind].values[yx]*1e-12
            bcds.vmr_e.values[:,np.where(lats >= 30),:] = csv[f"{var_name}_box0"][yind].values[yx]*1e-12

            bcds.vmr_w.values[:,np.where(lats <= -30),:] = csv[f"{var_name}_box3"][yind].values[yx]*1e-12
            bcds.vmr_w.values[:,np.where(np.logical_and(lats > -30, lats <= 0)),:] = csv[f"{var_name}_box2"][yind].values[yx]*1e-12
            bcds.vmr_w.values[:,np.where(np.logical_and(lats > 0, lats <= 30)),:] = csv[f"{var_name}_box1"][yind].values[yx]*1e-12
            bcds.vmr_w.values[:,np.where(lats >= 30),:] = csv[f"{var_name}_box0"][yind].values[yx]*1e-12
        
            bcds = bcds.assign_coords(time=np.atleast_1d((pd.to_datetime(f"{int(year)}-{int(yx+1)}-01"))))

            bcds.attrs.update({'author': f"{getpass.getuser()}@bristol.ac.uk", 
                               'date_created' : dt.today().strftime("%D"), 
                               'source' : 'AGAGE 12-box model',
                               'species':speciestocopy, 
                               'title':f'{speciestocopy} mixing ratio at domain edges',
                               'time period':'monthly' })
            bcds.to_netcdf(path = f"{outdir}/{domain}/{speciestocopy}_{domain}_{str(year)}{str(yx+1).zfill(2)}.nc", mode="w") 