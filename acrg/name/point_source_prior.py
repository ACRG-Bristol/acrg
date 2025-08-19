import xarray as xr
import numpy as np
from openghg_inversions.utils import areagrid
import getpass
from datetime import datetime
from openghg_defs import species_info_file
import json
import os
import glob
from acrg.config.paths import Paths

LPDM_path = Paths.lpdm

def point_source_prior(species, 
                       locs, 
                       domain, 
                       emissions_total_target, 
                       outdir, 
                       outname, 
                       expand = 1, 
                       add_sources = True):
    """
    Create a prior for point source emissions.

    Args:
        species (str) : The species to create the prior for.
        locs (list of tuples) : The (longitude, latitude) locations of the point sources.
        domain (str) : The model domain to use.
        emissions_total_target (float) : The total target emissions for the point sources.
        expand (int) : The number of grid cells to expand the point sources.
        add_sources (bool) : Whether to add to existing sources or replace them, in the case that sources overlap.
        outdir (str) : The output directory for the NetCDF file.
        outname (str) : The name of the output NetCDF file.
    """

    with species_info_file.open() as f:
        species_info = json.load(f)
    molar_mass = float(species_info[species.upper()]['mol_mass'])

    template_path = os.path.join(LPDM_path, 'emissions', domain, '*')

    fn = sorted(glob.glob(template_path))
    with xr.open_dataset(fn[0]) as load:
        ds_template = load.load()


    ds = ds_template.copy()
    ds.flux.values = np.zeros_like(ds.flux.values)

    areas = areagrid(ds.lat.values, ds.lon.values)

    emissions_per_single_source = emissions_total_target / len(locs)

    n_cells_per_source = (expand * 2 + 1) ** 2  # Total number of cells covered by each source (central cell + surrounding cells)
    index_range = list(np.arange(-expand, expand + 1))


    # Set flux values to the calculated emissions per source at specified locations and surrounding 8 cells
    for lon_target, lat_target in locs:
        # Find the closest grid cell indices
        lon_idx = np.argmin(np.abs(ds.lon.values - lon_target))
        lat_idx = np.argmin(np.abs(ds.lat.values - lat_target))
        area = areas[lat_idx, lon_idx] # set this to area of the central cell. 
        
        # Set the central cell and 8 surrounding cells to 1
        # Flux dimensions are (lat, lon, time)
        for di in index_range:  # longitude offsets
            for dj in index_range:  # latitude offsets
                new_lon_idx = lon_idx + di
                new_lat_idx = lat_idx + dj
                
                
                if (0 <= new_lon_idx < len(ds.lon) and 0 <= new_lat_idx < len(ds.lat)):
                    if add_sources: # add to existing values
                        ds.flux.values[new_lat_idx, new_lon_idx, :] += emissions_per_single_source / (area * n_cells_per_source * 365.25 * 24 * 60 * 60 * molar_mass * 1e-9)
                    else: # replace existing values
                        ds.flux.values[new_lat_idx, new_lon_idx, :] = emissions_per_single_source / (area * n_cells_per_source * 365.25 * 24 * 60 * 60 * molar_mass * 1e-9)


    locs_string = ""
    for loc in locs:
        locs_string += f"({loc[0]:.2f}, {loc[1]:.2f}), "

    emissions_total_actual = np.sum(ds.flux.values[:,:,0]*areas) * 365.25* 24 * 60 * 60 * molar_mass * 1e-9

    attrs = {'title': f'Point source prior for {species} in {domain} domain',
            'sources': locs_string,
            'emissions_total_target': f'{round(emissions_total_target, 4)} Gg/yr',
            'emissions_total_actual': f'{round(emissions_total_actual, 4)} Gg/yr',
            'comment': 'target and actual totals may disagree if overlapping sources and add_sources=False',
            'created_by': f'{getpass.getuser()}@bristol.ac.uk',
            'date_created': datetime.today().strftime('%Y-%m-%d'),
            'species': species,
            }

    ds.attrs = attrs
    outpath = os.path.join(outdir, outname)
    ds.to_netcdf(outpath, mode='w')