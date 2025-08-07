import xarray as xr
import numpy as np
import getpass
from datetime import datetime as dt

def create_land_sea_file(domain):

    """
    Creates a land-sea mask file for a specified domain.
    
    Args:
    domain (str): The domain for which the land-sea mask is created (e.g., 'EUROPE'). There must be a 
                  corresponding file in the LPDM/countries directory named 'country_{domain}.nc'.

    Returns:
    None: The function saves the land-sea mask as a NetCDF file named 'country-land-sea_{domain}.nc' in the LPDM/countries directory.
    """
    

    ds = xr.open_dataset(f'/group/chem/acrg/LPDM/countries/country_{domain}.nc')
    mask = np.where(ds.country.values == 0, 0, 1)
    ds_out = ds.copy()
    ds_out.country.values = mask
    ds_out.attrs['title'] = f'grid of land-sea split across {domain} domain'
    ds_out.attrs['created from'] = f'/group/chem/acrg/LPDM/countries/country_{domain}.nc'
    ds_out.attrs['created_by'] = f'{getpass.getuser()}@bristol.ac.uk'
    ds_out.attrs['created'] = dt.today().strftime("%D")
    ds_out.to_netcdf(f'/group/chem/acrg/LPDM/countries/country-land-sea_{domain}.nc')