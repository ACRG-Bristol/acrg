import iris
import xarray as xr
import os

from acrg.config.paths import Paths


data_path = Paths.data

folder = data_path / 'LPDM/Model/NAMEIII_v7_2_particlelocation_satellite/Resources/Topog/'

# topog function loads UM file in Iris and saves as NetCDF in the same folder you are working in

def global_topog(filename, output_directory=None):
    '''where filename is the name of the UM file as a str'''
    cube = iris.load(folder / filename)[0]
    netcdf_file = filename.replace('.pp','.nc')

    iris.save(cube, netcdf_file)
    
# corrects the longitude and saves new netcdf file in topog_NAME in the shared area. The initial netcdf in your folder can be deleted (optional)

    ds = xr.open_dataset(netcdf_file)
    ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
    new_netcdf = netcdf_file.replace('.nc','_global.nc') 
    
    ds.to_netcdf(path = output_directory / new_netcdf)
    
