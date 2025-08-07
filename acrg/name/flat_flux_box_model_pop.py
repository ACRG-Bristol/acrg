import numpy as np
import xarray as xr
from openghg_inversions.utils import areagrid
import pandas as pd
from openghg_defs import species_info_file, domain_info_file
import json
from datetime import datetime
import getpass

def population_scaling(domain):
    """
    Calculates the fraction of the global population in a specified domain, according to the UN WPP-Adjusted Population Density dataset for 2020.

    Args:
        domain (str): The domain for which the population fraction is calculated (e.g., 'EASTASIA'). Must correspond to a domain in openghg_defs. 
        
    Returns:
        float: The fraction of the global population in the specified domain.
    """
    # open the domain info file

    with domain_info_file.open() as f:
        domain_info = json.load(f)
    
    domain_lat_min = domain_info[domain]['latitude_range'][0]
    domain_lat_max = domain_info[domain]['latitude_range'][1]
    domain_lon_min = domain_info[domain]['longitude_range'][0]
    domain_lon_max = domain_info[domain]['longitude_range'][1]


    # find the fraction of the global population in the EASTASIA domain
    population_filepath = '/group/chem/acrg/LPDM/world_population_density.nc'  # this is a gridded file of population density. Raster layer 5 has values for 2020.
    ds_population = xr.open_dataset(population_filepath).sel(raster=5)
    total_population = ds_population.sum(dim=['longitude', 'latitude'])['UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 30 arc-minutes'].values

    domain_population = ds_population.sel(latitude=slice(domain_lat_max, domain_lat_min), longitude=slice(domain_lon_min, domain_lon_max)).sum(dim=['longitude', 'latitude'])['UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 30 arc-minutes'].values
    domain_fraction = domain_population / total_population
    
    return domain_fraction


def extract_box_model_total(species, pathtobm, year_to_scale=2023):
    """
    Extracts and box model emissions for a given domain and species based on population fraction.
    
    Args:
        pathtobm (str): Path to the py12box_agage repository. This can be cloned from https://github.com/mrghg/py12box_agage
        species (str): The species for which the emissions are extracted. Use IUPAC names (e.g. 'CFC-12', 'HFC-134a', 'CCl4').
        year_to_scale (int): The year for which the emissions are extracted. Defaults to 2023.
        
    Returns:
        total (float): The total emissions for the specified domain and species in Gg/yr.
    """

    csv = pd.read_csv(pathtobm+f"{species}/outputs/{species}_Global_annual_emissions.csv", comment="#", index_col='Year')
    year_to_scale = 2023
    total = csv.loc[year_to_scale]['Global_annual_emissions']
    return total

def create_flat_flux_box_model_pop(domain, species, pathtobm, year_to_scale, year_range, outdir, outname):
    """
    Creates a flat flux box model for a specified domain and species, scaling emissions based on population fraction.
    This scales a single year of emissions (year_to_scale) to the population fraction in the specified domain, and creates a 
    prior for the specified year range.

    Args:
        domain (str): The domain for which the box model is created (e.g., 'EASTASIA').
        species (str): The species for which the box model is created. Should use IUPAC names (e.g. 'CFC-12', 'HFC-134a', 'CCl4').
        pathtobm (str): Path to the py12box_agage repository. This can be cloned from https://github.com/mrghg/py12box_agage
        year_to_scale (int): The year for which the emissions are scaled.
        year_range (tuple): A tuple containing the start and end years for the prior to be created (e.g., (2008, 2025)).
        outdir (str): Output directory where the NetCDF file will be saved. Should already have a {domain} subdirectory.
        outname (str): Name of the output file. File will start with species and domain (e.g., if outname=='flat_population', the output filename 
                       will be '{species}_{domain}_flat_population.nc').
    """

    # open a sample dataset to copy

    pathtocopy = f'/group/chem/acrg/LPDM/countries/country_{domain}.nc' # just a sample dataset to copy the grid from
    ds = xr.open_dataset(pathtocopy).copy()
    ds = ds.rename({'country': 'flux'})
    ds['flux'].values = 1 # set all fluxes to 1, we will scale them later  

    # find the fraction of the global population in the EASTASIA domain

    domain_fraction = population_scaling(domain)

    # find the total emissions for the species in the domain

    global_total = extract_box_model_total(species, pathtobm, year_to_scale)
    total_target = global_total * domain_fraction
    
    with species_info_file.open() as f:
        species_info = json.load(f)


    areas = areagrid(ds.lat.values, ds.lon.values)
    ds_total = np.sum(ds.flux.values[:,:,0] * areas) * 365.25 * 24 * 60 * 60 * np.float64(species_info[species.replace('-','')]['mol_mass']) * 1e-9 # what values currently are in Gg/yr
    ds_out = ds.copy()
    ds_out['flux'] = ds.flux * (total_target / ds_total)  # scale the fluxes to match the target total emissions


    # do all years in the specified range

    years = pd.date_range(f'{year_range[0]}-01-01', f'{year_range[1]}-01-01', freq='YS')

    repeated = []
    for t in years:
        ds_new = ds_out.copy(deep=True)
        ds_new = ds_new.assign_coords(time=[t])
        repeated.append(ds_new)

    ds_out = xr.concat(repeated, dim='time')

    # Add attributes to the dataset

    attrs = {'title': f'Flat prior for {species} emissions in EASTASIA domain',
            'comment': f'Box model emissions for {year_to_scale} scaled by population fraction in EASTASIA domain ({domain_fraction.round(4)})',
            'emissions_total': f'{total_target.round(4)} Gg/yr',
            'created_by': f'{getpass.getuser()}@bristol.ac.uk', 
            'date_created': datetime.today().strftime('%Y-%m-%d'),
            'regridder_used': 'openghg_inversions.utils.areagrid',
            'species': species,
            'population_file_used': 'UN WPP-Adjusted Population Count, v4.11 for 2020',
            'population_file_reference':'https://doi.org/10.7927/H4F47M65'}
    ds_out.attrs = attrs

    # Write

    ds_out.to_netcdf(f'{outdir}/{domain}/{species.lower()}_{domain}_{outname}.nc', mode='w')