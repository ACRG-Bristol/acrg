import iris
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
from os.path import join
import sys

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH")
    data_path = os.getenv("DATA_PATH")
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg
    data_path = paths.data


folder = join(data_path,'LPDM/Model/NAMEIII_v7_2_particlelocation_satellite/Resources/Topog/')

# topog function loads UM file in Iris and saves as NetCDF in the same folder you are working in

def global_topog(filename):
    '''where filename is the name of the UM file as a str'''
    cube = iris.load(folder+filename)[0]
    netcdf_file = filename.replace('.pp','.nc')

    iris.save(cube, netcdf_file)

# corrects the longitude and saves new netcdf file in topog_NAME in the shared area. The initial netcdf in your folder can be deleted (optional)

    ds = xr.open_dataset(netcdf_file)
    ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
    new_netcdf = netcdf_file.replace('.nc','_global.nc')

    ds.to_netcdf(path = data_path + '/LPDM/topog_NAME/' + new_netcdf)
