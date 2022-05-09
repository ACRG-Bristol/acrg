#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:00:20 2019

Use makeCAMSBC

@author: rt17603

---------------------------------------------------------------------------
Example:

make_climatology = False

# start and end dates to for estimating climatology
clim_start = None
clim_end   = None

# start and end dates of the output BC files
start      = "2015-01-01"
end        = "2016-01-01"

# geographical domain to create BC for
domain  = "EUROPE"

# directory in which to save the BC files
outdir  = os.path.join(data_path)

# use makeCAMSBC to create the CAMS BC files with the inputs defined above
cams_bc = makeCAMSBC(domain           = domain,
                     start            = start,
                     end              = end,
                     outdir           = outdir,
                     cams_directory   = cams_directory,
                     clim_start       = clim_start,
                     clim_end         = clim_end,
                     make_climatology = make_climatology)
---------------------------------------------------------------------------

"""
from datetime import datetime as dt
from acrg.countrymask import domain_volume
from acrg.satellite.gosat import extract_files
from acrg.convert import concentration
import json
from collections import OrderedDict
from acrg.config.paths import Paths
import os
import numpy as np
import xarray as xr
from acrg.name.name import open_ds
from acrg.config.paths import Paths
from . import climatology

acrg_path = Paths.acrg
data_path = Paths.data
cams_directory_default = os.path.join(data_path, 'ECMWF_CAMS', 'CAMS_inversion/')

# latest version is currently v19, update as necessary
latest_version = 'v19'

# default variable names in the CAMS files for each species and CAMS versions
default_inputs = {'ch4': {'v19': {'altitude': 'altitude',
                                  'file_start_str': 'cams73',
                                  'height': 'height',
                                  'hlevel': 'hlevel',
                                  'lat': 'latitude',
                                  'level': 'level',
                                  'lon': 'longitude',
                                  'time': 'time',
                                  'z': 'z'}
                         }
                   }
default_inputs['ch4']['latest'] = default_inputs['ch4'][latest_version]

with open(os.path.join(acrg_path, "data/species_info.json")) as f:
    species_info=json.load(f,object_pairs_hook=OrderedDict)

def readCAMSInversion(start, end, species='ch4', cams_directory=cams_directory_default,
                      version='latest', variables=None):
    '''
    Args:
        start, end (str)
            Start and end dates
        species (str, optional)
            defaults to 'ch4'
        cams_directory (str, optional)
            Path to CAMS files
        version (str, int, optional))
            version number of inversion, e.g. 'v19'
            if an int is entered, it will be converted to a str starting with 'v'
            defaults to 'latest'
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
            
    Returns:
        xarray dataset of combined CAMS inversions
    '''
    # set default inputs if not given
    cams_directory = cams_directory_default if cams_directory is None else cams_directory
    version = f'v{version}' if not isinstance(version, str) else version

    variables = default_inputs[species][version] if variables is None else variables
    
    search_str = f"{variables['file_start_str']}*{version}*{species}*.nc"
    
    files = extract_files(cams_directory, search_str, start=start, end=end, day=False)
    ds_list = [open_ds(file) for file in files]
    
    ds = ds_list[0] if len(ds_list)==1 else xr.concat(ds_list, dim=variables['time'])
    if species.upper() in ds:
        ds  = ds.rename({species.upper(): species})
    
    return ds

def convertCAMSaltitude(ds, species='ch4', version='latest'):
    '''
    Convert altitude coordinate to level
    
    Args
        ds (xarray.Dataset)
            CAMS data
        species, version (str, optional)
            used to infer the variable names from the 'default_inputs' dictionary
            default to species='ch4' and version='latest'
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
            
    Returns
        xarray.Dataset
    '''
    version = f'v{version}' if not isinstance(version, str) else version
    # get default variable names
    variables = default_inputs[species][version]

    hlevels = ds[variables['hlevel']].values
    h_min = hlevels[0]
    h_max = hlevels[-1]

    low_alt = ds[variables['altitude']].sel(**{variables['hlevel']:slice(h_min,h_max-1)})
    high_alt = ds[variables['altitude']].sel(**{variables['hlevel']:slice(h_min+1,h_max)})
    
    altitude_diff = high_alt.values - low_alt.values
    z = low_alt + altitude_diff/2.
    
    z_dims = tuple([dim if dim!=variables['hlevel'] else "level" for dim in ds[variables['altitude']].dims])
    
    ds = ds.assign(**{"z":(z_dims,z.data)})
    ds["z"] = ds["z"].transpose(*(variables['time'],"level",variables['lat'],variables['lon'])) if variables['time'] in ds[variables['altitude']].dims else \
              ds["z"].transpose(*("level",variables['lat'],variables['lon']))
    
    return ds

def interpheight(nesw, fp_height, lonorlat, species='ch4', version='latest', reverse=None, variables=None):
    """
    Interpolates the CAMS data to the NAME heights
    
    Args:
        nesw (dataset)
            The N, E, S or W BC CAMS boundary data
        fp_height (array)
            NAME footprint heights
        lonorlat (string)
            Whether you're interpolating along the longitude (N or S) or the latitude (E or W).
        species (string)
            The gas species of interest, defaults to 'ch4'
        version (string)
            The version to use, used to infer the variable names from the 'default_inputs' dictionary
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
        reverse (bool/None)
            Whether height values within is nesw input are in reverse order 
            (i.e. nesw["z"] level 1 values > nesw["z"] level 2 values).
            Default = None. If this is set to None this will be automatically determined.
            
    Returns:
        dataset : CAMS BC data interpolated to NAME heights
        
    """
    version = f'v{version}' if not isinstance(version, str) else version
    # get default variable names
    variables = default_inputs[species][version] if variables is None else variables

    if lonorlat.lower() in ['longitude', 'lon']:     
        interp = np.zeros((len(fp_height),len(nesw.longitude) ))
    elif lonorlat.lower() in ['latitude', 'lat']:
        interp = np.zeros((len(fp_height),len(nesw.latitude) ))
    else:
        print("Please specify either lonorlat='longitude' (or 'lon') or 'latitude' (or 'lat')")
        return None
    
    if reverse is None:
        z_coords = nesw[variables['level']].values
        z_0 = nesw[variables['z']].sel(**{variables['level']:z_coords[0]}).values[0]
        z_1 = nesw[variables['z']].sel(**{variables['level']:z_coords[1]}).values[0]

        if z_1 >= z_0:
            reverse=False
        elif z_0 > z_1:
            reverse=True        
    
    for jj in range(len(nesw[variables['z']][0,:])):
        if reverse == True:
            interp[:,jj] = np.interp(fp_height, nesw[variables['z']][:,jj][::-1], nesw[species][:,jj][::-1]).astype(np.float)
        elif reverse == False:
            interp[:,jj] = np.interp(fp_height, nesw[variables['z']][:,jj], nesw[species][:,jj]).astype(np.float)
    
    ds2 = xr.DataArray(interp, coords=[fp_height, nesw[lonorlat].values], dims=['height', lonorlat])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def interplonlat(nesw, fp_lonorlat, lonorlat, species='ch4', reverse=None, version='latest', variables=None, verbose=False):
    """
    Interpolates the CAMS data to the NAME longitudes and latitudes
    
    Args:
        nesw (dataset)
            The N, E, S or W BC CAMS boundary data
        fp_lonorlat (array)
            NAME footprint longitudes or latitudes
        species (string)
            The gas species of interest
        lonorlat (string)
            Whether you're interpolating along the 'longitude' 
            (N or S) or the 'latitude' (E or W).
        reverse (bool/None)
            Whether lon or lat values within nesw input are in reverse order
            (i.e. nesw[lonorlat][0] > nesw[lonorlat][1]).
            Default = None. If this is set to None this will be automatically determined.
        version (string)
            The version to use, used to infer the variable names from the 'default_inputs' dictionary
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
            
    Returns:
        dataset : CAMS BC data interpolated to NAME longitudes or latitudes
        
    """
    version = f'v{version}' if not isinstance(version, str) else version
    # get default variable names
    variables = default_inputs[species][version] if variables is None else variables

    if reverse is None:
        if nesw[lonorlat].values[1] >= nesw[lonorlat].values[0]:
            reverse=False
        elif nesw[lonorlat].values[0] > nesw[lonorlat].values[1]:
            reverse=True
    
    if verbose: print("Reversing heights")
    
    interp = np.zeros(( len(nesw[variables['height']]),len(fp_lonorlat) ))
    for jj in range(len(nesw[variables['height']])):
        if reverse == True:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values[::-1], nesw[species][jj,:][::-1]).astype(np.float)
        elif reverse == False:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values, nesw[species][jj,:]).astype(np.float)
            
    ds2 = xr.DataArray(interp, coords=[nesw[variables['height']].values, fp_lonorlat], dims=[variables['height'], lonorlat[0:3]])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def bc_filename(domain, start_date, species='ch4', from_climatology=False):
    '''
    Create a standardised filename for CAMS boundary conditions files
    
    Args:
        domain (str)
            Georgaphical domain, e.g. 'EUROPE', 'PACIFIC', 'USA'
        start_date (str)
            Start date of BC, e.g. '2015-01-01'
        species (str, optional)
            Gas species, e.g. 'ch4', 'co2'
            defaults to 'ch4'
        from_climatology (bool)
            Whether the BC has been produced from climatologies
    
    Returns:
        A standardised filename (str)
    '''
    clim_str = "_climatology" if from_climatology else ""
    date_str = dt.strptime(start_date, '%Y-%m-%d').strftime('%Y%m')
    
    return f"{species.lower()}_{domain}_{date_str}{clim_str}.nc"


def write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, st_date, domain, outdir, gridsize,
                           species='ch4', version='latest', variables=None, from_climatology=False):
    """
    Writes the CAMS BC data to a netcdf file.
    This uses bc_filename to produce a standardised file name for the boundary conditions file.
    
    Args:
        vmr_n (array)
            Molar ratio at northern boundary
        vmr_e (array)
            Molar ratio at eastern boundary
        vmr_s (array)
            Molar ratio at western boundary
        vmr_w (array)
            Molar ratio at southern boundary
        st_date (str)
            Start date of form "YYYY-MM-dd"
        domain (str)
            The domain which you want the boundary conditions for.
        outdir (str)
            Directory for writing output file.
        gridsize (int/float)
            Resolution of CAMS output in degrees.
            Possible are: 0.125, 0.25, 0.4, 0.5, 0.75, 1, 1.125, 1.5, 2, 2.5, 3
        species (str, optional)
            The gas species e.g. 'ch4', 'co2'
            defaults to 'ch4'
        version (str, optional)
            CAMS version to use, used to infer the variable names from the 'default_inputs' dictionary
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
        from_climatology (bool)
            Whether BCs were produced from climatologies
    
    Returns
        netcdf file: Boundary conditions at domain boundaries
    """
    version = f'v{version}' if not isinstance(version, str) else version
    # get default variable names
    variables = default_inputs[species][version] if variables is None else variables

    BC_edges = vmr_n.merge(vmr_e).merge(vmr_s).merge(vmr_w)
    BC_edges.expand_dims(variables['time'], 2)
    BC_edges.coords[variables['time']] = (dt.strptime(st_date, '%Y-%m-%d'))
    
    BC_edges.attrs['title']           = f"ECMWF CAMS {species} volume mixing ratios at domain edges"
    BC_edges.attrs['CAMS_resolution'] = gridsize
    BC_edges.attrs['author']          = os.getenv('USER')
    BC_edges.attrs['date_created']    = np.str(dt.today())
    
    if not os.path.isdir(outdir): os.makedirs(outdir)
    
    BC_filename = bc_filename(domain = domain,
                              species = species,
                              start_date = st_date,
                              from_climatology = from_climatology)
    BC_edges.to_netcdf(path = os.path.join(outdir, BC_filename), mode = 'w')

def create_CAMS_BC(ds, fp_lat, fp_lon, fp_height, date, domain, species=None, version='latest', variables=None,
                   outdir=None, from_climatology=False, verbose=False, test=False):
    '''
    Create CAMS boundary conditions and write to a netcdf file
    
    This uses interplonlat and interpheight to regrid the boundary conditions to match the NAME grid.
    Then writes to a netcdf file using write_CAMS_BC_tonetcdf.
    
    Args
        ds (xarray.Dataset)
            CAMS dataset
        fp_lat, fp_lon (numpy.ndarray)
            NAME footprint latitudes and longitudes
        date (str)
            Start date of form "YYYY-MM-dd"
        domain (str)
            Geographical region of interest, e.g. 'EUROPE', 'PACIFIC', 'USA'
        species (str)
            Gas species of interest, e.g. 'ch4', 'co2'
            defaults to 'ch4'
        version (str, optional)
            CAMS version to use, used to infer the variable names from the 'default_inputs' dictionary
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
        outdir (str, optional)
            Directory in which to save outputs
            Outputs will be saved to a subdirectory : outdir/LPDM/bc/<domain>
        from_climatology (bool, optional)
            Whether BCs are being produced from climatologies
        verbose (bool, optional)
            Whether to print any updates
        test (bool, optional)
            If creating boundary conditions in a test units will not be converted
            Should never be True outside of testing
    
    Returns
        netcdf file
            Writes CAMS BC xarray.Datasets to a netcdf file in outdir/LPDM/bc/<domain>
            with a standardised filename:
            <species>_<domain>_<start_year><start_month>.nc'
    
    '''
    version = f'v{version}' if not isinstance(version, str) else version
    # get default variable names
    variables = default_inputs[species][version] if variables is None else variables

    lat_grid = np.mean(ds[variables['lat']][1:] - ds[variables['lat']][:-1])
    lon_grid = np.mean(ds[variables['lon']][1:] - ds[variables['lon']][:-1])
    gridsize = f"{lat_grid} x {lon_grid}"
    
    ds = convertCAMSaltitude(ds)
    ds = ds.mean(variables['time'])
    
    # find the correct unit conversion between mol/mol and species specific parts-per- units
    if not test: 
        conversion = concentration(species_info[species.upper()]['units'])
        ds[species].values *= conversion
    else:
        if verbose: print('Warning: testing, units will not be converted')
    
    #Select the gridcells closest to the edges of the  domain and make sure outside of fp
    lat_n = (np.abs(ds.coords[variables['lat']].values - max(fp_lat))).argmin()
    if ds.coords[variables['lat']].values[lat_n] < np.max(fp_lat) and lat_n != 0:
        lat_n -= 1
    lat_s = (np.abs(ds.coords[variables['lat']].values - min(fp_lat))).argmin()
    if ds.coords[variables['lat']].values[lat_s] > np.min(fp_lat) and lat_s != (len(ds.coords[variables['lat']].values)-1):
        lat_s += 1
    lon_e = (np.abs(ds.coords[variables['lon']].values - max(fp_lon))).argmin()
    if ds.coords[variables['lon']].values[lon_e] < max(fp_lon) and lon_e != (len(ds.coords[variables['lon']].values)-1):
        lon_e += 1
    lon_w = (np.abs(ds.coords[variables['lon']].values - min(fp_lon))).argmin()
    if ds.coords[variables['lon']].values[lon_w] > min(fp_lon) and lon_w != 0:
        lon_e -= 1
    #Cut to these and then interpolate
    north = ds.sel(latitude  = ds.coords[variables['lat']][lat_n],
                   longitude = slice(ds.coords[variables['lon']][lon_w],ds.coords[variables['lon']][lon_e])).drop_vars([variables['lat']])
    south = ds.sel(latitude  = ds.coords[variables['lat']][lat_s],
                   longitude = slice(ds.coords[variables['lon']][lon_w],ds.coords[variables['lon']][lon_e])).drop_vars([variables['lat']])
    east  = ds.sel(longitude = ds.coords[variables['lon']][lon_e],
                   latitude  = slice(ds.coords[variables['lat']][lat_s],ds.coords[variables['lat']][lat_n])).drop_vars([variables['lon']])
    west  = ds.sel(longitude = ds.coords[variables['lon']][lon_w],
                   latitude  = slice(ds.coords[variables['lat']][lat_s],ds.coords[variables['lat']][lat_n])).drop_vars([variables['lon']])
    
    vmr_n = interplonlat(interpheight(north, fp_height, species=species, lonorlat=variables['lon']),
                         fp_lon, species=species, lonorlat=variables['lon'], verbose=False).rename({species : 'vmr_n'})   
    vmr_s = interplonlat(interpheight(south, fp_height, species=species, lonorlat=variables['lon']),
                         fp_lon, species=species, lonorlat=variables['lon'], verbose=False).rename({species : 'vmr_s'}) 
    vmr_e = interplonlat(interpheight(east, fp_height, species=species, lonorlat=variables['lat']),
                         fp_lat, species=species, lonorlat=variables['lat'], verbose=False).rename({species : 'vmr_e'}) 
    vmr_w = interplonlat(interpheight(west, fp_height, species=species, lonorlat=variables['lat']),
                         fp_lat, species=species, lonorlat=variables['lat'], verbose=False).rename({species : 'vmr_w'})      

    write_CAMS_BC_tonetcdf(vmr_n = vmr_n,
                           vmr_e = vmr_e,
                           vmr_s = vmr_s,
                           vmr_w = vmr_w,
                           st_date = date,
                           domain = domain,
                           outdir = outdir,
                           gridsize = gridsize,
                           species = species,
                           from_climatology = from_climatology) 

def makeCAMSBC(domain, start, end,
               species = 'ch4',
               cams_version = 'v19',
               outdir = None,
               cams_directory = cams_directory_default,
               clim_start = None, clim_end = None,
               make_climatology = False,
               fp_directory = None,
               variables=None,
               verbose = False, overwrite = False, test = False):
    '''
    Make boundary condition files from the CAMS inversion product
    
    This uses:
        - readCAMSInversion to find and open the CAMS inversion files for the given date range
        - climatology.monthly_cycle and climatology.add_full_time to produce climatologies (if required)
        - create_CAMS_BC to produce boundary conditions for multiple date-time intervals within the date range
    
    Args:
        domain (str): 
            Domain name, e.g., "EUROPE"
        start (str): 
            Start date of output, e.g., "2015-01-01"
        end (str): 
            End date of output, e.g., "2016-01-01"
        species (str, optional)
            Gas species, defaults to 'ch4'
        cams_version (str, int)
            version number of CAMS inversion, e.g. 'v19'
            if an int is entered, it will be converted to a str starting with 'v'
        outdir (str, optional)
            Output directory to save output. Output will automatically be written to outdir/LPDM/bc/DOMAIN
        cams_directory (str, optional)
            Location of CAMS inversion output
        clim_start (str, optional): 
            Start date to average fields into a climatology, e.g., "2010-01-01"
            Default of none will assume start
        clim_end (str, optional): 
            End date to average fields into a climatology, e.g., "2010-01-01"    
            Default of none will assume end
        make_climatology (bool, optional)
            If True climatologies will be created from the dates clim_start & clim_end, or start & end if clim_start and clim_end are None
        fp_directory (str/None, optional):
            Footprint folder to use when extracting domain volume. By default, this will use your data_path and expected file structure by default.
            If specified, sub-folder of domain name will be assumed.
        variables (dict, optional)
            dictionary containing variable names including:
                'altitude', 'file_start_str', 'height', 'hlevel', 'lat', 'lon', 'time', & 'z'
            if not given variables will be assumed to be the same as for CAMS CH4 v19
        verbose (bool, optional)
            If False do not print updates
        overwrite (bool, optional)
            If True, files will automatically be overwritten
            If False and if the file already exists, the user will be asked to confirm that the file should be overwritten
        test (bool, optional)
            If creating boundary conditions in a test units will not be converted
            Should never be True outside of testing
 
    Returns:
        netcdf files of monthly boundary condition curtain files corresponding to the edges of the corresponding NAME domain
        
    Example:
        makeCAMSBC(domain, start, end, clim_start, clim_end)
    
    Todo: 
        This code currently works as such:  if you want real CAMS inversion for a month (not a climatology), then you can only specify one year
        at a time. If multiple years are given in start and end, it will average them together. This code needs to be split up with a 
        climatology keywork so that it does not automatically do averaging.
    '''
    cams_version = f'v{cams_version}' if not isinstance(cams_version, str) else cams_version
    # get default variable names
    variables = default_inputs[species][cams_version] if variables is None else variables
    cams_directory = cams_directory_default if cams_directory is None else cams_directory

    # rename clim_start and clim_end if None
    if (clim_start is None and clim_end is None) or not make_climatology:
        clim_start = start
        clim_end   = end
    
    # Use fp_directory if specified.
    kwargs = {'domain': domain} if fp_directory is None else {'domain': domain, 'fp_directory': fp_directory}
    fp_lat,fp_lon,fp_height = domain_volume(**kwargs)
    
    outdir  = os.path.join(data_path, 'LPDM', 'bc', domain) if outdir is None else outdir
    
    cams_ds = readCAMSInversion(clim_start, clim_end, cams_directory = cams_directory,
                                version=cams_version)
    
    # create climatology if required
    if make_climatology:
        cams_seasonal = climatology.monthly_cycle(cams_ds)
        cams_ds       = climatology.add_full_time(cams_seasonal, start = start, end = end)
    
    date_range = np.arange(start, end, dtype="datetime64[M]")
    date_range = [np.datetime_as_string(date)+"-01" for date in date_range]
    date_range += [end]

    for start, end in zip(date_range[:-1], date_range[1:]):

        # Check if file exists so a current file doesn't get overwritten
        out_filename = bc_filename(domain=domain, start_date=start, species=species, from_climatology=make_climatology)
        
        if overwrite and verbose:
            print(f'Boundary condition file {os.path.join(outdir, out_filename)} already exists and is being overwritten.')
        if os.path.isfile(os.path.join(outdir, out_filename)) and not overwrite:
            print(f'Boundary condition file {os.path.join(outdir, out_filename)} already exists.')
            answer = input("You are about to overwrite an existing file, do you want to continue? Y/N ")
            if answer.upper() == 'N':
                continue
            elif answer.upper() == 'Y':
                pass
        
        # select the data for the correct date range
        ds = cams_ds.sel(**{variables['time'] : slice(start, end)}) 

        create_CAMS_BC(ds               = ds,
                       fp_lat           = fp_lat,
                       fp_lon           = fp_lon,
                       fp_height        = fp_height,
                       date             = start,
                       species          = species,
                       version          = cams_version,
                       domain           = domain,
                       outdir           = outdir,
                       verbose          = verbose,
                       from_climatology = make_climatology,
                       test             = test)
   