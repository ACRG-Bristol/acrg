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

with open(os.path.join(acrg_path, "data/species_info.json")) as f:
    species_info=json.load(f,object_pairs_hook=OrderedDict)

def readCAMSInversion(start, end, species=None, cams_directory=None,
                      version=None, file_start_str=None, time_var=None):
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
            defaults to 'v19'
        file_start (str, optional)
            string at the start of the filename, used to search for files
            defaults to 'cams73'
        time_var (str, optional)
            name of the time variable, defaults to 'time'
            
    Returns:
        xarray dataset of combined CAMS inversions
    '''
    # set default inputs if not given
    species = 'ch4' if species is None else species
    file_start_str = 'cams73' if file_start_str is None else file_start_str
    cams_directory = os.path.join(data_path, 'ECMWF_CAMS', 'CAMS_inversion/') if cams_directory is None else cams_directory
    version = 'v19' if version is None else version
    version = f'v{version}' if not isinstance(version, str) else version
    time_var = 'time' if time_var is None else time_var
    
    search_str = f"{file_start_str}*{version}*{species}*.nc"
    
    files = extract_files(cams_directory, search_str, start=start, end=end, day=False)
    ds_list = [open_ds(file) for file in files]
    
    ds = ds_list[0] if len(ds_list)==1 else xr.concat(ds_list, dim=time_var)
    if species.upper() in ds:
        ds  = ds.rename({species.upper(): species})
    
    return ds

def convertCAMSaltitude(ds, altitude_var=None, hlevel_var=None,
                        lat_var=None, lon_var=None, time_var=None):
    '''
    Convert altitude coordinate to level
    
    Args
        ds (xarray.Dataset)
            CAMS data
        altitude_var, hlevel_var, lat_var, lon_var, time_var (str, optional)
            name of the altitude, hlevel, latitude, longitude, and time variables
            default to: 'altitude', 'hlevel', 'latitude', 'longitude', and 'time'
            
    Returns
        xarray.Dataset
    '''
    # set defaults if not given
    altitude_var = 'altitude' if altitude_var is None else altitude_var
    hlevel_var = 'hlevel' if hlevel_var is None else hlevel_var
    lat_var = 'latitude' if lat_var is None else lat_var
    lon_var = 'longitude' if lon_var is None else lon_var
    time_var = 'time' if time_var is None else time_var

    hlevels = ds[hlevel_var].values
    h_min = hlevels[0]
    h_max = hlevels[-1]

    low_alt = ds[altitude_var].sel(**{hlevel_var:slice(h_min,h_max-1)})
    high_alt = ds[altitude_var].sel(**{hlevel_var:slice(h_min+1,h_max)})
    
    altitude_diff = high_alt.values - low_alt.values
    z = low_alt + altitude_diff/2.
    
    z_dims = tuple([dim if dim!=hlevel_var else "level" for dim in ds[altitude_var].dims])
    
    ds = ds.assign(**{"z":(z_dims,z)})
    ds["z"] = ds["z"].transpose(*(time_var,"level",lat_var,lon_var)) if time_var in ds[altitude_var].dims else \
              ds["z"].transpose(*("level",lat_var,lon_var))
    
    return ds

def interpheight(nesw, fp_height, species=None, lonorlat=None, reverse=None, z_var=None, level_var=None):
    """
    Interpolates the CAMS data to the NAME heights
    
    Args:
        nesw (dataset)
            The N, E, S or W BC CAMS boundary data
        fp_height (array)
            NAME footprint heights
        species (string)
            The gas species of interest
        lonorlat (string)
            Whether you're interpolating along the longitude (N or S) or the latitude (E or W).
        reverse (bool/None)
            Whether height values within is nesw input are in reverse order 
            (i.e. nesw["z"] level 1 values > nesw["z"] level 2 values).
            Default = None. If this is set to None this will be automatically determined.
        z_var: (str, optional)
            name of z (altitude) variable with level coordinate, default to 'z'
        level_var (str, optional) : Used if reverse is not defined to extract appropriate
            height values to compare.
            deaults to 'level'
            
    Returns:
        dataset : CAMS BC data interpolated to NAME heights
        
    """
    # set defaults if not given
    species = 'ch4' if species is None else species
    z_var = "z" if z_var is None else z_var
    level_var = "level" if level_var is None else level_var

    if lonorlat.lower() in ['longitude', 'lon']:     
        interp = np.zeros((len(fp_height),len(nesw.longitude) ))
    elif lonorlat.lower() in ['latitude', 'lat']:
        interp = np.zeros((len(fp_height),len(nesw.latitude) ))
    else:
        print("Please specify either lonorlat='longitude' or 'latitude'")
        return None
    
    if reverse is None:
        z_coords = nesw[level_var].values
        z_0 = nesw[z_var].sel(**{level_var:z_coords[0]}).values[0]
        z_1 = nesw[z_var].sel(**{level_var:z_coords[1]}).values[0]

        if z_1 >= z_0:
            reverse=False
        elif z_0 > z_1:
            reverse=True        
    
    for jj in range(len(nesw[z_var][0,:])):
        if reverse == True:
            interp[:,jj] = np.interp(fp_height, nesw[z_var][:,jj][::-1], nesw[species][:,jj][::-1]).astype(np.float)
        elif reverse == False:
            interp[:,jj] = np.interp(fp_height, nesw[z_var][:,jj], nesw[species][:,jj]).astype(np.float)
    
    ds2 = xr.DataArray(interp, coords=[fp_height, nesw[lonorlat].values], dims=['height', lonorlat])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def interplonlat(nesw, fp_lonorlat, species=None, lonorlat=None, reverse=None, height_var=None, verbose=False):
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
        height_var (str, optional)
            name of the height variable, defaults to 'height'
            
    Returns:
        dataset : CAMS BC data interpolated to NAME longitudes or latitudes
        
    """
    # set defaults if not given
    species = 'ch4' if species is None else species
    height_var = 'height' if height_var is None else height_var

    if reverse is None:
        if nesw[lonorlat].values[1] >= nesw[lonorlat].values[0]:
            reverse=False
        elif nesw[lonorlat].values[0] > nesw[lonorlat].values[1]:
            reverse=True
    
    if verbose: print("Reversing heights")
    
    interp = np.zeros(( len(nesw[height_var]),len(fp_lonorlat) ))
    for jj in range(len(nesw[height_var])):
        if reverse == True:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values[::-1], nesw[species][jj,:][::-1]).astype(np.float)
        elif reverse == False:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values, nesw[species][jj,:]).astype(np.float)
            
    ds2 = xr.DataArray(interp, coords=[nesw[height_var].values, fp_lonorlat], dims=[height_var, lonorlat[0:3]])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def bc_filename(domain, start_date, species=None, from_climatology=False):
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
    # set defaults if not given
    species = 'ch4' if species is None else species

    clim_str = "_climatology" if from_climatology else ""
    date_str = dt.strptime(start_date, '%Y-%m-%d').strftime('%Y%m')
    
    return f"{species.lower()}_{domain}_{date_str}{clim_str}.nc"


def write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, st_date, domain, outdir, gridsize,
                           species=None, from_climatology=False, time_var=None):
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
        from_climatology (bool)
            Whether BCs were produced from climatologies
    
    Returns
        netcdf file: Boundary conditions at domain boundaries
    """
    # set defaults if not given
    time_var = 'time' if time_var is None else time_var
    species = 'ch4' if species is None else species

    BC_edges = vmr_n.merge(vmr_e).merge(vmr_s).merge(vmr_w)
    BC_edges.expand_dims(time_var, 2)
    BC_edges.coords[time_var] = (dt.strptime(st_date, '%Y-%m-%d'))
    
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

def create_CAMS_BC(ds, fp_lat, fp_lon, fp_height, date, domain, species=None,
                   outdir=None, from_climatology=False, verbose=False, test=False,
                   time_var=None, lat_var=None, lon_var=None):
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
        time_var, lat_var, lon_var (str, optional)
            names of the time, latitude and longitude variables
            default to 'time', 'latitude', 'longitude'
    
    Returns
        netcdf file
            Writes CAMS BC xarray.Datasets to a netcdf file in outdir/LPDM/bc/<domain>
            with a standardised filename:
            <species>_<domain>_<start_year><start_month>.nc'
    
    '''
    # set defaults if not given
    time_var = 'time' if time_var is None else time_var
    lat_var = 'latitude' if lat_var is None else lat_var
    lon_var = 'longitude' if lon_var is None else lon_var
    species = 'ch4' if species is None else species

    lat_grid = np.mean(ds[lat_var][1:] -  ds[lat_var][:-1])
    lon_grid = np.mean(ds[lon_var][1:] - ds[lon_var][:-1])
    gridsize = f"{lat_grid} x {lon_grid}"
    
    ds = convertCAMSaltitude(ds)
    ds = ds.mean(time_var)
    
    # find the correct unit conversion between mol/mol and species specific parts-per- units
    if not test: 
        conversion = concentration(species_info[species.upper()]['units'])
        ds[species].values *= conversion
    else:
        if verbose: print('Warning: testing, units will not be converted')
    
    #Select the gridcells closest to the edges of the  domain and make sure outside of fp
    lat_n = (np.abs(ds.coords[lat_var].values - max(fp_lat))).argmin()
    if ds.coords[lat_var].values[lat_n] < np.max(fp_lat) and lat_n != 0:
        lat_n -= 1
    lat_s = (np.abs(ds.coords[lat_var].values - min(fp_lat))).argmin()
    if ds.coords[lat_var].values[lat_s] > np.min(fp_lat) and lat_s != (len(ds.coords[lat_var].values)-1):
        lat_s += 1
    lon_e = (np.abs(ds.coords[lon_var].values - max(fp_lon))).argmin()
    if ds.coords[lon_var].values[lon_e] < max(fp_lon) and lon_e != (len(ds.coords[lon_var].values)-1):
        lon_e += 1
    lon_w = (np.abs(ds.coords[lon_var].values - min(fp_lon))).argmin()
    if ds.coords[lon_var].values[lon_w] > min(fp_lon) and lon_w != 0:
        lon_e -= 1
    #Cut to these and then interpolate
    north = ds.sel(latitude  = ds.coords[lat_var][lat_n],
                   longitude = slice(ds.coords[lon_var][lon_w],ds.coords[lon_var][lon_e])).drop_vars([lat_var])
    south = ds.sel(latitude  = ds.coords[lat_var][lat_s],
                   longitude = slice(ds.coords[lon_var][lon_w],ds.coords[lon_var][lon_e])).drop_vars([lat_var])
    east  = ds.sel(longitude = ds.coords[lon_var][lon_e],
                   latitude  = slice(ds.coords[lat_var][lat_s],ds.coords[lat_var][lat_n])).drop_vars([lon_var])
    west  = ds.sel(longitude = ds.coords[lon_var][lon_w],
                   latitude  = slice(ds.coords[lat_var][lat_s],ds.coords[lat_var][lat_n])).drop_vars([lon_var])
    
    vmr_n = interplonlat(interpheight(north, fp_height, species=species, lonorlat=lon_var),
                         fp_lon, species=species, lonorlat=lon_var, verbose=False).rename({species : 'vmr_n'})   
    vmr_s = interplonlat(interpheight(south, fp_height, species=species, lonorlat=lon_var),
                         fp_lon, species=species, lonorlat=lon_var, verbose=False).rename({species : 'vmr_s'}) 
    vmr_e = interplonlat(interpheight(east, fp_height, species=species, lonorlat=lat_var),
                         fp_lat, species=species, lonorlat=lat_var, verbose=False).rename({species : 'vmr_e'}) 
    vmr_w = interplonlat(interpheight(west, fp_height, species=species, lonorlat=lat_var),
                         fp_lat, species=species, lonorlat=lat_var, verbose=False).rename({species : 'vmr_w'})      

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
               species = None, 
               outdir = None,
               cams_directory = None,
               cams_version = None,
               clim_start = None, clim_end = None,
               make_climatology = False,
               fp_directory = None,
               time_var = None,
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
        outdir (str, optional)
            Output directory to save output. Output will automatically be written to outdir/LPDM/bc/DOMAIN
        cams_directory (str, optional)
            Location of CAMS inversion output
        cams_version (str, int)
            version number of CAMS inversion, e.g. 'v19'
            if an int is entered, it will be converted to a str starting with 'v'
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
        time_var (str, optional)
            name of the time variable, defaults to 'time'
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
    # set defaults if not given
    time_var = 'time' if time_var is None else time_var
    species = 'ch4' if species is None else species
    cams_version = 'v19' if cams_version is None else cams_version
    cams_directory = os.path.join(data_path, 'ECMWF_CAMS', 'CAMS_inversion/') if cams_directory is None else cams_directory

    # rename clim_start and clim_end if None
    if (clim_start is None and clim_end is None) or not make_climatology:
        clim_start = start
        clim_end   = end
    
    # Use fp_directory if specified.
    if fp_directory is not None:
        fp_lat,fp_lon,fp_height = domain_volume(domain, fp_directory=fp_directory)
    else:
        fp_lat,fp_lon,fp_height = domain_volume(domain)
    
    outdir  = os.path.join(data_path, 'LPDM', 'bc', domain) if outdir is None else outdir
    
    cams_ds = readCAMSInversion(clim_start, clim_end, cams_directory = cams_directory, version=cams_version)
    
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
        ds = cams_ds.sel(**{time_var : slice(start, end)}) 

        create_CAMS_BC(ds               = ds,
                       fp_lat           = fp_lat,
                       fp_lon           = fp_lon,
                       fp_height        = fp_height,
                       date             = start,
                       species          = species,
                       domain           = domain,
                       outdir           = outdir,
                       verbose          = verbose,
                       from_climatology = make_climatology,
                       test             = test)
   