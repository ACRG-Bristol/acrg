import xarray as xr
import pandas as pd
from datetime import datetime as dt
import getpass

def outer_regions_from_UKMO(filepath, domain, outdir='/user/home/qq23840/openghg_inversions/openghg_inversions/basis'):
    """
    Reads outer region definitions from a UK Met Office file and converts it to an dataset for use in OpenGHG Inversions.
    
    Args:
        filepath (str): Path to the UK Met Office outer region definition file.
        domain (str): The domain for which the outer regions are defined (e.g., 'EASTASIA').
        outdir (str): The output directory for the NetCDF file. Defaults to openghg_inversions.

    """

    ds_to_copy = xr.open_dataset(f'/group/chem/acrg/LPDM/countries/countries_{domain}.nc')


    intem = pd.read_csv(filepath, delim_whitespace=True, skiprows=10).rename(columns={'Lat':'lat', 'Lon':'lon', 'Code':'region'}).set_index(['lat','lon']).drop(columns=['i', 'j'])
    ds = intem.to_xarray().assign_coords({'lat':ds_to_copy.lat.values, 'lon':ds_to_copy.lon.values})
    ds.attrs["info"] = f"outer region definition for {domain} domain, aligned with UK Met Office definition"
    ds.attrs['created_by'] = f'{getpass.getuser()}@bristol.ac.uk, using UK Met Office data'
    ds.attrs['created'] = dt.today().strftime("%D")
    ds.attrs["created_from"] = filepath
    ds.attrs["domain"] = domain
    ds.to_netcdf(f'{outdir}/outer_region_definition_{domain}.nc')