""""
Script / how-to for adding subregions (e.g. specific states/provinces/counties) to an existing country file.

Requires a shapefile for the subregions of interest, and the country file must already exist.
Shapefiles can be found online - mine came from https://data.humdata.org/. When you download them make sure that
you also download the accompanying files (.shx, .dbf etc) and keep them all in the same folder. This allows the shapefile
to be read properly by the geopandas package.

"""

import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import numpy as np
import getpass
from datetime import datetime as dt

# Load your NetCDF file
ds = xr.open_dataset('/group/chem/acrg/LPDM/countries/country_eastchina_EASTASIA.nc')

# Load the shapefile (replace with your actual file path)
gdf = gpd.read_file('/user/home/qq23840/inversions_intercomparisons/haklim_h2402/rus_admbnda_adm2_gadm_2022_v02.shp')

# You then need to identify which polygon(s) in the shapefile correspond to your subregion of interest. The `gdf` variable
# is a GeoDataFrame, which is like a pandas DataFrame but with a geometry column. You can inspect it by printing the first few rows, and 
# searching through the columns to find one that contains the names of the subregions you're interested in adding. You can also plot the 
# geometries to help identify the right one:
# gdf.plot(column='ADM2_EN', legend=True)  # Replace 'ADM2_EN' with the appropriate column name for your shapefile. 

# In my case, it was just row 37, which corresponds to Primorsky Krai in Russia. 

primorsky_poly = gdf.geometry.iloc[37]  # Primorsky Krai polygon 

# Prepare grid points
lat_grid, lon_grid = np.meshgrid(ds.lat.values, ds.lon.values, indexing='ij')
points = [Point(lon, lat) for lat, lon in zip(lat_grid.ravel(), lon_grid.ravel())]

# Mask: True if point is inside Primorsky Krai
mask = np.array([primorsky_poly.contains(pt) for pt in points]).reshape(lat_grid.shape)

# Assign new country code
new_country_code = ds.country.values.max() + 1
ds.country.values[mask] = new_country_code

# Convert to list, append the name(s) of the new territorie(s), and assign back
names = [str(n) for n in ds.name.values] 
names.append("PRIMORSKY KRAI")
ds["name"] = ("ncountries", np.array(names, dtype=object))

# add some relevant attributes and save

attrs = {'description':'EASTASIA country map with East China and Primorsky Krai (Russia) added, for Haklim Choi H-2402 inversion',
         'source': 'Existing EASTASIA country map and shapefile for Primorsky Krai https://data.humdata.org/dataset/cod-ab-rus',
         'author': getpass.getuser(),
         'created': dt.now().strftime("%Y-%m-%d %H:%M:%S")
         }
ds.to_netcdf('/group/chem/acrg/LPDM/countries/country_eastchina_primorsky_EASTASIA.nc')