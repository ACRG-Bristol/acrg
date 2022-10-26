
import numpy as np
import acrg.name.name as name
import xarray as xr
from acrg.config.paths import Paths
import os

data_path = Paths.data

def regrid_country(ds, country_filename, domain):
    '''
    Regrid a dataset or dataarray to the lat/lons of a desired country file. 

    Args:
        ds (xarray dataset or dataarray)
            Dataset to be changed to new lat/lons. Must have lat/lon dimensions
        country_filename (str)
            Name of country file to be regridded to 
        domain (str)
            All caps name of domain
    
    '''

    country_file = os.path.join(data_path, 'LPDM/countries/country_'+country_filename+ '_'+ domain +'.nc')

    c_object = name.get_country('ARCTIC', country_file=country_file)
    cntryds = xr.Dataset({'country': (['lat','lon'], c_object.country), 
                        'name' : (['ncountries'],c_object.name) },
                                        coords = {'lat': (c_object.lat),
                                        'lon': (c_object.lon)})

    lonmin_cds = np.min(cntryds.lon.values)
    lonmax_cds = np.max(cntryds.lon.values)
    latmin_cds = np.min(cntryds.lat.values)
    latmax_cds = np.max(cntryds.lat.values)

    lon_ds = ds.lon
    lat_ds = ds.lat

    lonmin_ds = lon_ds[np.where(np.isclose(lon_ds, lonmin_cds, atol = 0.01, rtol=0))[0][0]]
    lonmax_ds = lon_ds[np.where(np.isclose(lon_ds, lonmax_cds, atol = 0.01, rtol=0))[0][0]]
    latmin_ds = lat_ds[np.where(np.isclose(lat_ds, latmin_cds, atol = 0.01, rtol=0))[0][0]]
    latmax_ds = lat_ds[np.where(np.isclose(lat_ds, latmax_cds, atol = 0.01, rtol=0))[0][0]]    

    ds_country = ds.sel(lon=slice(lonmin_ds,lonmax_ds),lat=slice(latmin_ds,latmax_ds))

    return ds_country


