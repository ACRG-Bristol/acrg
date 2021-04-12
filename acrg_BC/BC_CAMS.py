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
from acrg_countrymask import domain_volume
from acrg_satellite.gosat import extract_files
from acrg_convert import concentration
import json
from collections import OrderedDict
from .__init__ import *

with open(os.path.join(acrg_path, "acrg_species_info.json")) as f:
    species_info=json.load(f,object_pairs_hook=OrderedDict)

# data_path      = paths.data
# cams_directory = os.path.join(data_path, 'ECMWF_CAMS', 'CAMS_inversion/')

def readCAMSInversion(start, end, species='ch4', cams_directory=cams_directory,
                      version='v19'):
    '''
    Args:
        start, end (str)
            Start and end dates
        species (str)
        cams_directory (str)
            Path to CAMS files
        version (str, int)
            version number of inversion, e.g. 'v19'
            if an int is entered, it will be converted to a str starting with 'v'
            
    Returns:
        xarray dataset of combined CAMS inversions
    '''
    version = f'v{version}' if not isinstance(version, str) else version
    
    search_str = f"cams73*{version}*{species}*.nc"
    
    files   = extract_files(cams_directory, search_str, start=start, end=end, day=False)
    ds_list = [open_ds(file) for file in files]
    
    ds      = ds_list[0] if len(ds_list)==1 else xr.concat(ds_list, dim="time")
    if species.upper() in ds:
        ds  = ds.rename({species.upper(): species})
    
    return ds

def convertCAMSaltitude(ds):
    '''
    Convert altitude coordinate to level
    
    Args
        ds (xarray.Dataset)
            CAMS data
            
    Returns
        xarray.Dataset
    '''
    
    hlevels  = ds.hlevel.values
    h_min    = hlevels[0]
    h_max    = hlevels[-1]

    low_alt  = ds["altitude"].sel(**{"hlevel":slice(h_min,h_max-1)})
    high_alt = ds["altitude"].sel(**{"hlevel":slice(h_min+1,h_max)})
    
    altitude_diff = high_alt.values - low_alt.values
    z = low_alt + altitude_diff/2.
    
    z_dims = tuple([dim if dim!="hlevel" else "level" for dim in ds["altitude"].dims])
    
    ds = ds.assign(**{"z":(z_dims,z)})
    ds["z"] = ds["z"].transpose(*("time","level","latitude","longitude"))
    
    return ds

def interpheight(nesw, fp_height, species, lonorlat=None, reverse=None, z_dim="level"):
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
            Whether you're interpolating along the 'longitude' 
            (N or S) or the 'latitude' (E or W).
        reverse (bool/None)
            Whether height values within is nesw input are in reverse order 
            (i.e. nesw["z"] level 1 values > nesw["z"] level 2 values).
            Default = None. If this is set to None this will be automatically determined.
        z_dim (str, optional) : Used if reverse is not defined to extract appropriate
            height values to compare.
            
    Returns:
        dataset : CAMS BC data interpolated to NAME heights
        
    """
    if lonorlat == 'longitude':     
        interp = np.zeros((len(fp_height),len(nesw.longitude) ))
    elif lonorlat == 'latitude':
        interp = np.zeros((len(fp_height),len(nesw.latitude) ))
    else:
        print("Please specify either lonorlat='longitude' or 'latitude'")
        return None
    
    if reverse is None:
        z_coords = nesw[z_dim].values
        z_0 = nesw['z'].sel(**{z_dim:z_coords[0]}).values[0]
        z_1 = nesw['z'].sel(**{z_dim:z_coords[1]}).values[0]

        if z_1 >= z_0:
            reverse=False
        elif z_0 > z_1:
            reverse=True        
    
    for jj in range(len(nesw['z'][0,:])):
        if reverse == True:
            interp[:,jj] = np.interp(fp_height, nesw['z'][:,jj][::-1], nesw[species][:,jj][::-1]).astype(np.float)
        elif reverse == False:
            interp[:,jj] = np.interp(fp_height, nesw['z'][:,jj], nesw[species][:,jj]).astype(np.float)
    
    ds2 = xr.DataArray(interp, coords=[fp_height, nesw[lonorlat].values], dims=['height', lonorlat])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def interplonlat(nesw, fp_lonorlat, species, lonorlat=None, reverse=None, verbose=False):
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
            
    Returns:
        dataset : CAMS BC data interpolated to NAME longitudes or latitudes
        
    """

    if reverse is None:
        if nesw[lonorlat].values[1] >= nesw[lonorlat].values[0]:
            reverse=False
        elif nesw[lonorlat].values[0] > nesw[lonorlat].values[1]:
            reverse=True
    
    if verbose: print("Reversing heights")
    
    interp = np.zeros(( len(nesw.height),len(fp_lonorlat) ))
    for jj in range(len(nesw.height)):
        if reverse == True:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values[::-1], nesw[species][jj,:][::-1]).astype(np.float)
        elif reverse == False:
            interp[jj, :] = np.interp(fp_lonorlat, nesw[lonorlat].values, nesw[species][jj,:]).astype(np.float)
            
    ds2 = xr.DataArray(interp, coords=[nesw.height.values, fp_lonorlat], dims=['height', lonorlat[0:3]])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def bc_filename(domain, species, start_date, from_climatology=False):
    '''
    Create a standardised filename for CAMS boundary conditions files
    
    Args:
        domain (str)
            Georgaphical domain, e.g. 'EUROPE', 'PACIFIC', 'USA'
        species (str)
            Gas species, e.g. 'ch4', 'co2'
        start_date (str)
            Start date of BC, e.g. '2015-01-01'
        from_climatology (bool)
            Whether the BC has been produced from climatologies
    
    Returns:
        A standardised filename (str)
    '''
    clim_str = "_climatology" if from_climatology else ""
    date_str = dt.strptime(start_date, '%Y-%m-%d').strftime('%Y%m')
    
    return f"{species.lower()}_{domain}_{date_str}{clim_str}.nc"


def write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, st_date, species, domain, outdir, gridsize,
                           from_climatology=False):
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
        species (str)
            The gas species e.g. 'ch4', 'co2'
        domain (str)
            The domain which you want the boundary conditions for.
        outdir (str)
            Directory for writing output file.
        gridsize (int/float)
            Resolution of CAMS output in degrees.
            Possible are: 0.125, 0.25, 0.4, 0.5, 0.75, 1, 1.125, 1.5, 2, 2.5, 3
        from_climatology (bool)
            Whether BCs were produced from climatologies
    
    Returns
        netcdf file: Boundary conditions at domain boundaries
    """
    BC_edges = vmr_n.merge(vmr_e).merge(vmr_s).merge(vmr_w)
    BC_edges.expand_dims('time', 2)
    BC_edges.coords['time'] = (dt.strptime(st_date, '%Y-%m-%d'))
    
    BC_edges.attrs['title']           = f"ECMWF CAMS {species} volume mixing ratios at domain edges"
    BC_edges.attrs['CAMS_resolution'] = gridsize
    BC_edges.attrs['author']          = os.getenv('USER')
    BC_edges.attrs['date_created']    = np.str(dt.today())
    
    if not os.path.isdir(outdir): os.makedirs(outdir)
    
    BC_filename = bc_filename(domain           = domain,
                              species          = species,
                              start_date       = st_date,
                              from_climatology = from_climatology)
    BC_edges.to_netcdf(path = os.path.join(outdir, BC_filename), mode = 'w')

def create_CAMS_BC(ds, fp_lat, fp_lon, fp_height, date, domain, species="ch4",
                   outdir=None, cams_directory=cams_directory, 
                   from_climatology=False, verbose=False, test=False):
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
        outdir (str, optional)
            Directory in which to save outputs
            Outputs will be saved to a subdirectory : outdir/LPDM/bc/<domain>
        cams_directory (str, optional)
            Directory in which the CAMS data is stored
            Defaults to <data_path>/ECMWF_CAMS/CAMS_inversion/
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
    lat_grid = np.mean(ds["latitude"][1:] -  ds["latitude"][:-1])
    lon_grid = np.mean(ds["longitude"][1:] - ds["longitude"][:-1])
    gridsize = f"{lat_grid} x {lon_grid}"
    
    ds = convertCAMSaltitude(ds)
    ds = ds.mean('time')
    
    # find the correct unit conversion between mol/mol and species specific parts-per- units
    if not test: 
        conversion = concentration(species_info[species.upper()]['units'])
        ds[species].values *= conversion
    else:
        if verbose: print('Warning: testing, units will not be converted')
    
    #Select the gridcells closest to the edges of the  domain and make sure outside of fp
    lat_n = (np.abs(ds.coords['latitude'].values - max(fp_lat))).argmin()
    if ds.coords['latitude'].values[lat_n] < np.max(fp_lat) and lat_n != 0:
        lat_n -= 1
    lat_s = (np.abs(ds.coords['latitude'].values - min(fp_lat))).argmin()
    if ds.coords['latitude'].values[lat_s] > np.min(fp_lat) and lat_s != (len(ds.coords['latitude'].values)-1):
        lat_s += 1
    lon_e = (np.abs(ds.coords['longitude'].values - max(fp_lon))).argmin()
    if ds.coords['longitude'].values[lon_e] < max(fp_lon) and lon_e != (len(ds.coords['longitude'].values)-1):
        lon_e += 1
    lon_w = (np.abs(ds.coords['longitude'].values - min(fp_lon))).argmin()
    if ds.coords['longitude'].values[lon_w] > min(fp_lon) and lon_w != 0:
        lon_e -= 1
    #Cut to these and then interpolate
    north = ds.sel(latitude  = ds.coords['latitude'][lat_n],
                   longitude = slice(ds.coords['longitude'][lon_w],ds.coords['longitude'][lon_e])).drop_vars(['latitude'])
    south = ds.sel(latitude  = ds.coords['latitude'][lat_s],
                   longitude = slice(ds.coords['longitude'][lon_w],ds.coords['longitude'][lon_e])).drop_vars(['latitude'])
    east  = ds.sel(longitude = ds.coords['longitude'][lon_e],
                   latitude  = slice(ds.coords['latitude'][lat_s],ds.coords['latitude'][lat_n])).drop_vars(['longitude'])
    west  = ds.sel(longitude = ds.coords['longitude'][lon_w],
                   latitude  = slice(ds.coords['latitude'][lat_s],ds.coords['latitude'][lat_n])).drop_vars(['longitude'])
    
    vmr_n = interplonlat(interpheight(north, fp_height, species, lonorlat='longitude'),
                         fp_lon, species, lonorlat='longitude', verbose=False).rename({species : 'vmr_n'})   
    vmr_s = interplonlat(interpheight(south, fp_height, species, lonorlat='longitude'),
                         fp_lon, species, lonorlat='longitude', verbose=False).rename({species : 'vmr_s'}) 
    vmr_e = interplonlat(interpheight(east, fp_height, species, lonorlat='latitude'),
                         fp_lat, species, lonorlat='latitude', verbose=False).rename({species : 'vmr_e'}) 
    vmr_w = interplonlat(interpheight(west, fp_height, species, lonorlat='latitude'),
                         fp_lat, species, lonorlat='latitude', verbose=False).rename({species : 'vmr_w'})      

    write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, date, species, domain, outdir, gridsize, from_climatology=from_climatology) 

def makeCAMSBC(domain, start, end,
               species = 'ch4',
               outdir = None, cams_directory = cams_directory,
               cams_version = 'v19',
               clim_start = None, clim_end = None,
               make_climatology = False,
               verbose = False, overwrite = False, test=False):
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
    
    # rename clim_start and clim_end if None
    if (clim_start is None and clim_end is None) or not make_climatology:
        clim_start = start
        clim_end   = end
    
    # find a file with the NAME grid to match the input file to
    # is testing, this will be in the test files
    if test:
        fp_directory = os.path.join(acrg_path, 'tests', 'files', 'LPDM', 'fp_NAME')
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
        out_filename = bc_filename(domain, species, start, from_climatology=make_climatology)
        
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
        ds = cams_ds.sel(**{"time" : slice(start, end)}) 

        create_CAMS_BC(ds, fp_lat, fp_lon, fp_height, start,
                       species          = species,
                       domain           = domain,
                       outdir           = outdir,
                       cams_directory   = cams_directory,
                       verbose          = verbose,
                       from_climatology = make_climatology,
                       test             = test)
   