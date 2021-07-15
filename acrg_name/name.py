# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 10:45:51 2014

"""
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import netCDF4 as nc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import datetime as dt
import os
import sys
import glob
from matplotlib import ticker
import pandas as pd
import bisect
import subprocess
import json
from os.path import join
import xarray as xr
from acrg_time import convert
import calendar
import pickle
from scipy.interpolate import interp1d
import dateutil.relativedelta
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits import mplot3d
from collections import OrderedDict
import acrg_obs as obs
from acrg_config.paths import paths
from acrg_utils import is_number
from tqdm import tqdm
import dask.array as da
from numba import jit

import dask
dask.config.set({"array.slicing.split_large_chunks": True})

acrg_path = paths.acrg
data_path = paths.data

# Get acrg_site_info file
with open(acrg_path / "acrg_site_info.json") as f:
    site_info=json.load(f,object_pairs_hook=OrderedDict)

    
def open_ds(path, chunks=None):
    """
    Function efficiently opens xray datasets.

    Args:
        path (str)
        chunks (dict)
            size of chunks for each dimension
            e.g. {'lat': 50, 'lon': 50}
            opens dataset with dask, such that it is opened 'lazily'
            and all of the data is not loaded into memory
            defaults to None - dataset is opened with out dask
    """
    # use a context manager, to ensure the file gets closed after use
    if chunks is not None:
        ds = xr.open_mfdataset(path, chunks=chunks)
    else:
        with xr.open_dataset(path) as ds:
            ds.load()
    return ds 


def filenames(site, domain, start, end, height, fp_directory, met_model = None, network=None, species=None):
    """
    The filenames function outputs a list of available footprint file names,
    for given site, domain, directory and date range.
    
    Expect filenames of the form:
        [fp_directory]/domain/site*-height-species*domain*ym*.nc or [fp_directory]/domain/site*-height_domain*ym*.nc
        e.g. /data/shared/LPDM/fp_NAME/EUROPE/HFD-UKV-100magl-rn_EUROPE_202012.nc or /data/shared/LPDM/fp_NAME/EUROPE/MHD-10magl_EUROPE_201401.nc 
    
    Args:
        site (str) : 
            Site name. Full list of site names should be defined within acrg_site_info.json
        domain (str) : 
            Domain name. The footprint files should be sub-categorised by the NAME domain name.
        start (str) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
        end (str) : 
            End date in same format as start for range of files to find.
        height (str) : 
            Height related to input data. 
        fp_directory (str) :
            fp_directory can be specified if files are not in the default directory must point to a directory 
            which contains subfolders organized by domain.
        met_model (str):
            Met model used to run NAME
            Default is None and implies the standard global met model
            Alternates include 'UKV' and this met model must be in the outputfolder NAME        
        network (str, optional) : 
            Network for site. 
            If not specified, first entry in acrg_site_info.json file will be used (if there are multiple).
        species (str, optional)
            If specified, will search for species specific footprint files.
    Returns:
        list (str): matched filenames
    """
        
    # Read site info for heights
    if height is None:
        if not site in list(site_info.keys()):
            print("Site code not found in arcg_site_info.json to get height information. " + \
                  "Check that site code is as intended. "+ \
                  "If so, either add new site to file or input height manually.")
            return None
        if network is None:
            network = list(site_info[site].keys())[0]
        height = site_info[site][network]["height_name"][0]
    
    with open(os.path.join(acrg_path,"acrg_species_info.json")) as f:
        species_info=json.load(f)
 
    if species:
        species_obs = obs.read.synonyms(species, species_info)
        
        if 'lifetime' in species_info[species_obs].keys():
            lifetime = species_info[species_obs]["lifetime"]
            lifetime_hrs = convert.convert_to_hours(lifetime)
            # if a monthly lifetime is a list, use the minimum lifetime 
            # in the list to determine whether a species specific footprint is needed
            if type(lifetime) == list:
                lifetime_hrs = min(lifetime_hrs)
        else:
            lifetime_hrs = None
    else:
        lifetime_hrs = None
    
    # Convert into time format
    months = pd.date_range(start = start, end = end, freq = "M").to_pydatetime()
    yearmonth = [str(d.year) + str(d.month).zfill(2) for d in months]

    # first search for species specific footprint files.
    # if does not exist, use integrated files if lifetime of species is over 2 months
    # if lifetime under 2 months and no species specific file exists, fail

    files = []
    for ym in yearmonth:

        if species:

            if met_model:
                f=glob.glob(os.path.join(fp_directory,domain,f"{site}-{height}_{met_model}_{species}_{domain}_{ym}*.nc"))
            else:
                f=glob.glob(os.path.join(fp_directory,domain,f"{site}-{height}_{species}_{domain}_{ym}*.nc"))


        else:
            #manually create empty list if no species specified
            f = []
            
        if len(f) == 0:

            if met_model:
                glob_path = os.path.join(fp_directory,domain,f"{site}-{height}_{met_model}_{domain}_{ym}*.nc")
            else:
                glob_path = os.path.join(fp_directory,domain,f"{site}-{height}_{domain}_{ym}*.nc")
            
            if lifetime_hrs is None:
                print("No lifetime defined in species_info.json or species not defined. WARNING: 30-day integrated footprint used without chemical loss.")
                f=glob.glob(glob_path)
            elif lifetime_hrs <= 1440:
                print("This is a short-lived species. Footprints must be species specific. Re-process in process.py with lifetime")
                return []
            else:
                print("Treating species as long-lived.")
                f=glob.glob(glob_path)        
            
        if len(f) > 0:
            files += f

    files.sort()

    if len(files) == 0:
        print(f"Can't find footprints file: {glob_path}")
    return files


def read_netcdfs(files, dim = "time", chunks=None, verbose=True):
    """
    The read_netcdfs function uses xarray to open sequential netCDF files and 
    and concatenates them along the specified dimension.
    Note: this function makes sure that file is closed after open_dataset call.
    
    Args:
        files (list) : 
            List of netCDF filenames.
        dim (str, optional) : 
            Dimension of netCDF to use for concatenating the files.
            Default = "time".
        chunks (dict)
            size of chunks for each dimension
            e.g. {'lat': 50, 'lon': 50}
            opens dataset with dask, such that it is opened 'lazily'
            and all of the data is not loaded into memory
            defaults to None - dataset is opened with out dask
    
    Returns:
        xarray.Dataset : 
            All files open as one concatenated xarray.Dataset object    
    """
    if verbose:
        print("Reading and concatenating files: ")
        for fname in files:
            print(fname)
    
    datasets = [open_ds(p, chunks=chunks) for p in sorted(files)]
    
    # reindex all of the lat-lon values to a common one to prevent floating point error differences
    with xr.open_dataset(files[0]) as temp:
        fields_ds = temp.load()
    fp_lat = fields_ds["lat"].values
    fp_lon = fields_ds["lon"].values

    datasets = [ds.reindex(indexers={"lat":fp_lat, "lon":fp_lon}, method="nearest", tolerance=1e-5) for ds in datasets]

    combined = xr.concat(datasets, dim)
    return combined   


def footprints(sitecode_or_filename, met_model = None, fp_directory = None, 
               start = None, end = None, domain = None, height = None, network = None,
               species = None, HiTRes = False, chunks = None, verbose=True):

    """
    The footprints function loads a NAME footprint netCDF files into an xarray Dataset.
    Flux and boundary conditions are also loaded if species or emissions_name is specified.
    
    Either specify:
        a) A file name:
            fp = footprints(filename)
        b) A site code, domain, and date range:
            fp = footprints("MHD", 
                    start = "2014-01-01", end = "2014-01-01",
                    domain = "EUROPE")
    
    See filenames() function for expected format of files.

    Args:
        sitecode_or_filename (str) : 
            Site (e.g. 'MHD') or a netCDF filename (*.nc) (str)
        met_model (dict or str)     : Met model used to run NAME for each site in data
                               Default is None and implies the standard global met model used for all sites
                               met_model = 'UKV' will work if all sites are the same 'UKV'
                               {"MHD":'UKV', "JFJ":None} or {"MHD":'UKV', "TAC":"UKV"}
                               {"MHD":None, "TAC":None} is equivalent to met_model = None
                               {"MHD":'UKV', "TAC":'UKV'} is equivalent to met_model = 'UKV'
        fp_directory (str/dict) : 
            If the high time resolution footprints are used (HiTRes = True) fp_directory must be a dictionary
            of the form :
                fp_directory = {"integrated":PATH_TO_INTEGRATED_FP,"HiTRes":PATH_TO_HIGHTRES_FP}
            otherwise can be a single string if only integrated FPs are used and non-default.
        start (str) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
        end (str) : 
            End date in same format as start for range of files to find.
        domain (str) : 
            Domain name. The footprint files should be sub-categorised by the domain.
        height (str, optional) : 
            Height related to NAME. If the HEIGHT keyword is not specified, the default height from the 
            acrg_site_info.json file is assumed.
        network (str, optional) :
            Network associated with the site. This is only needed if there are multiple network options.
            This is used to extract site information e.g. height from acrg_site_info.json. 
        species (str) : 
            Species name. All species names are defined acrg_species_info.json.
        HiTRes (bool, optional) : 
            Whether to include high time resolution footprints.
            Default = False.
        chunks (dict)
            size of chunks for each dimension
            e.g. {'lat': 50, 'lon': 50}
            opens dataset with dask, such that it is opened 'lazily'
            and all of the data is not loaded into memory
            defaults to None - dataset is opened with out dask
        
    Returns:
        xarray.Dataset : 
            Combined footprint files
    """
    # Chose whether we've input a site code or a file name
    # If it's a three-letter site code, assume it's been processed
    # into an annual footprint file in (mol/mol) / (mol/m2/s)
    # using acrg_name_process
    
    if fp_directory is None:
        fp_integrated_directory = join(data_path, 'LPDM/fp_NAME/')
        fp_HiTRes_directory = join(data_path,'LPDM/fp_NAME_high_time_res/')
        fp_directory = {'integrated': fp_integrated_directory,
                        'HiTRes': fp_HiTRes_directory}

    #   Error message if looking for HiTRes files and fp_directory is not a dictionary        
        if HiTRes is True:
            if type(fp_directory) is not dict:
                print("fp_directory needs to be a dictionary containing paths \
                       to integrated and HiTRes footprints \
                       {integrated:path1, HiTRes:path2}")
                return None
    
    if '.nc' in sitecode_or_filename:
        if not '/' in sitecode_or_filename:
            files = [os.path.join(fp_directory, sitecode_or_filename)]
        else:
            files=[sitecode_or_filename]
    else:
        site=sitecode_or_filename[:]

        # Finds integrated footprints if specified as a dictionary with multiple entries (HiTRes = True) 
        # or a string with one entry        
        if type(fp_directory) is dict:
            files = filenames(site, domain, start, end, height, fp_directory["integrated"], met_model = met_model, network=network, species=species)
        else:
            files = filenames(site, domain, start, end, height, fp_directory, met_model = met_model, network=network, species=species)

    if len(files) == 0:
        print(f"\nWarning: Can't find footprint files for {sitecode_or_filename}. \
              This site will not be included in the output dictionary.\n")
        return None

    else:
        fp=read_netcdfs(files, chunks=chunks, verbose=verbose)  

        # if HiTRes == True:
        #     HiTRes_files = filenames(site, domain, start, end, height, fp_directory["HiTRes"])
        #     HiTRes_ds = read_netcdfs(HiTRes_files)
        #     fp = combine_datasets(fp, HiTRes_ds, method='ffill')

        return fp


def flux(domain, species, start = None, end = None, flux_directory=None, chunk=True, verbose=True):
    """
    The flux function reads in all flux files for the domain and species as an xarray Dataset.
    Note that at present ALL flux data is read in per species per domain or by emissions name.
    To be consistent with the footprints, fluxes should be in mol/m2/s.
    
    Expect filenames of the form:
        [flux_directory]/domain/species.lower()_*.nc
        e.g. [/data/shared/LPDM/emissions]/EUROPE/ch4_EUROPE_2013.nc
    
    TODO: This may get slow for very large flux datasets, and we may want to subset.
    
    Args:
        domain (str) : 
            Domain name. The flux files should be sub-categorised by the domain.
        species (str) : 
            Species name. All species names are defined acrg_species_info.json.
            The source can also be specified as part of the species name e.g. "co2-ff"
            will find all co2 emissions file for the fossil fuel sector.
        start (str, optional) : 
            Start date in format "YYYY-MM-DD" to output only a time slice of all the flux files.
            The start date used will be the first of the input month. I.e. if "2014-01-06" is input,
            "2014-01-01" will be used.  This is to mirror the time slice functionality of the filenames 
            function.
        end (str, optional) : 
            End date in same format as start to output only a time slice of all the flux files.
            The end date used will be the first of the input month and the timeslice will go up
            to, but not include, this time. I.e. if "2014-02-25' is input, "2014-02-01" will be used.
            This is to mirror the time slice functionality of the filenames function.
        flux_directory (str, optional) : 
            flux_directory can be specified if files are not in the default directory.
            Must point to a directory which contains subfolders organized by domain.
    Returns:
        xarray.Dataset : combined dataset of all matching flux files
    """
    
    if flux_directory is None:
        flux_directory = join(data_path, 'LPDM/emissions/')
    
    if "-" not in species:
        species_search_list = [(species+"-total").lower(), species.lower()]
    else:
        species_search_list = [species.lower()]
    
    for species_search in species_search_list:
        filename = os.path.join(flux_directory,domain,f"{species_search}_*.nc")
        if verbose:
            print("\nSearching for flux files: {}".format(filename))
    
        files = sorted(glob.glob(filename))
        
        if len(files) > 0:
            break
    else:
        raise IOError("\nError: Can't find flux files for domain '{0}' and species '{1}' (using search: {2})".format(domain,species,species_search_list))
    
    chunks = None if start == None or end == None or not chunk else {'lat': 50, 'lon': 50}
    flux_ds = read_netcdfs(files, chunks=chunks)
    # Check that time coordinate is present
    if not "time" in list(flux_ds.coords.keys()):
        raise KeyError(f"ERROR: No 'time' coordinate in flux dataset for domain '{domain}' species '{species}'")

    # Check for level coordinate. If one level, assume surface and drop
    if "lev" in list(flux_ds.coords.keys()):
        print(f"WARNING: Can't support multi-level fluxes. Trying to remove 'lev' coordinate \
              from flux dataset for {domain}, {species}")
        if len(flux_ds.lev) > 1:
            print("ERROR: More than one flux level")
        else:
            return flux_ds.drop("lev")
        
    if start == None or end == None:
        print("To get fluxes for a certain time period you must specify a start or end date.")
        return flux_ds
    else:

        #Change timeslice to be the beginning and end of months in the dates specified.
        start = pd.to_datetime(start)
        month_start = dt.datetime(start.year, start.month, 1, 0, 0)
        
        end = pd.to_datetime(end)
        month_end = dt.datetime(end.year, end.month, 1, 0, 0) - \
                        dt.timedelta(seconds = 1)
           
        if 'climatology' in species:
            ndate = pd.to_datetime(flux_ds.time.values)
            if len(ndate) == 1:  #If it's a single climatology value
                dateadj = ndate - month_start  #Adjust climatology to start in same year as obs  
            else: #Else if a monthly climatology
                dateadj = ndate[month_start.month-1] - month_start  #Adjust climatology to start in same year as obs  
            ndate = ndate - dateadj
            flux_ds = flux_ds.update({'time' : ndate})  
            flux_tmp = flux_ds.copy()
            while month_end > ndate[-1]:
                ndate = ndate + pd.DateOffset(years=1)      
                flux_ds = xr.merge([flux_ds, flux_tmp.update({'time' : ndate})])
                    
        flux_timeslice = flux_ds.sel(time=slice(month_start, month_end))
        if np.logical_and(month_start.year != month_end.year, len(flux_timeslice.time) != dateutil.relativedelta.relativedelta(end, start).months):
            month_start = dt.datetime(start.year, 1, 1, 0, 0)
            flux_timeslice = flux_ds.sel(time=slice(month_start, month_end))
        if len(flux_timeslice.time)==0:
            flux_timeslice = flux_ds.sel(time=start, method = 'ffill')
            flux_timeslice = flux_timeslice.expand_dims('time',axis=-1)
            print("Warning: No fluxes available during the time period specified so outputting\
                          flux from %s" %flux_timeslice.time.values[0])
        else:
            if verbose:
                print("Slicing time to range {} - {}".format(month_start,month_end))
        
        flux_timeslice = flux_timeslice.compute()
        return flux_timeslice


def flux_for_HiTRes(domain, emissions_dict, start=None, end=None, flux_directory=None, verbose=True):
    """
    Creates a dictionary of high and low frequency fluxes for use with HiTRes footprints.
    
    Args:
        domain (str) : 
            Domain name. The flux files should be sub-categorised by the domain.
        emissions_dict (dict) : 
            This should be a dictionary of the form:
                {'high_freq':high_freq_emissions_name, 'low_freq':low_freq_emissions_name}
                e.g. {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}.
        start (str, optional) : 
            Start date in format "YYYY-MM-DD" to output only a time slice of all the flux files.
            The start date used will be the first of the input month. I.e. if "2014-01-06" is input,
            "2014-01-01" will be used.  This is to mirror the time slice functionality of the filenames 
            function.
        end (str, optional): 
            End date in same format as start to output only a time slice of all the flux files.
            The end date used will be the first of the input month and the timeslice will go up
            to, but not include, this time. I.e. if "2014-02-25' is input, "2014-02-01" will be used.
            This is to mirror the time slice functionality of the filenames function.
        flux_directory (str, optional) : 
            flux_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by domain.
    Returns:
        dictionary {'high_freq': flux_xray_dataset, 'low_freq': flux_xray_dataset} :
            Dictionary containing xray datasets for both high and low frequency fluxes.
    """
    
    if 'low_freq' not in list(emissions_dict.keys()):
        print("low_freq must be a key in the emissions_dict in order to combine with HiTRes footprints.")
        return None
    elif 'high_freq' not in list(emissions_dict.keys()):
        print("high_freq must be a key in the emissions_dict in order to use HiTRes footprints.")
        return None
    
    flux_dict = {}
    
    if start:
        #Get the month before the one one requested because this will be needed for the first few
        #days in timeseries_HiTRes to calculate the modleed molefractions for times when the footprints
        #are in the previous month.
        start = str(pd.to_datetime(start) - dateutil.relativedelta.relativedelta(months=1))
    
    fluxes_highfreq = flux(domain, emissions_dict['high_freq'], start = start, end = end,flux_directory=flux_directory,
                           verbose=verbose)
    fluxes_lowfreq = flux(domain, emissions_dict['low_freq'], start = start, end = end,flux_directory=flux_directory,
                          verbose=verbose)
    
    flux_dict['low_freq'] = fluxes_lowfreq
    flux_dict['high_freq'] = fluxes_highfreq
    
    return flux_dict


def boundary_conditions(domain, species, start = None, end = None, bc_directory=None):
    """
    The boundary_conditions function reads in the files with the global model vmrs at the domain edges 
    to give the boundary conditions as an xarray Dataset.

    Expect filenames of the form:
        [bc_directory]/domain/species.lower()_*.nc
        e.g. [/data/shared/LPDM/bc]/EUROPE/ch4_EUROPE_201301.nc

    Args:
        domain (str) : 
            Domain name. The boundary condition files should be sub-categorised by the domain.
        species (str) : 
            Species name. All species names are defined acrg_species_info.json.
        start (str, optional) : 
            Start date in format "YYYY-MM-DD" to output only a time slice of all the flux files.
            The start date used will be the first of the input month. I.e. if "2014-01-06" is input,
            "2014-01-01" will be used.  This is to mirror the time slice functionality of the filenames 
            function.
        end (str, optional): 
            End date in same format as start to output only a time slice of all the flux files.
            The end date used will be the first of the input month and the timeslice will go up
            to, but not include, this time. I.e. if "2014-02-25' is input, "2014-02-01" will be used.
            This is to mirror the time slice functionality of the filenames function.
        bc_directory (str, optional) : 
            bc_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by domain.

    Returns:
        xarray.Dataset : 
            Combined dataset of matching boundary conditions files
    """
    
    if bc_directory is None:
        bc_directory = join(data_path, 'LPDM/bc/')
    
    filenames = os.path.join(bc_directory,domain,f"{species.lower()}_*.nc")
    
    files       = sorted(glob.glob(filenames))
    file_no_acc = [ff for ff in files if not os.access(ff, os.R_OK)]
    files       = [ff for ff in files if os.access(ff, os.R_OK)]
    if len(file_no_acc)>0:
        print('Warning: unable to read all boundary conditions files which match this criteria:')
        [print(ff) for ff in file_no_acc]
    
    if len(files) == 0:
        print("Cannot find boundary condition files in {}".format(filenames))
        raise IOError(f"\nError: Cannot find boundary condition files for domain '{domain}' and species '{species}' ")

    bc_ds = read_netcdfs(files)

    if start == None or end == None:
        print("To get boundary conditions for a certain time period you must specify an end date.")
        return bc_ds
    else:
        #Change timeslice to be the beginning and end of months in the dates specified.
        start = pd.to_datetime(start)
        month_start = dt.datetime(start.year, start.month, 1, 0, 0)
        
        end = pd.to_datetime(end)
        month_end = dt.datetime(end.year, end.month, 1, 0, 0) - \
                    dt.timedelta(seconds = 1)
            
        bc_timeslice = bc_ds.sel(time=slice(month_start, month_end))
        if len(bc_timeslice.time)==0:
            bc_timeslice = bc_ds.sel(time=start, method = 'ffill')
            bc_timeslice = bc_timeslice.expand_dims('time',axis=-1)
            print(f"No boundary conditions available during the time period specified so outputting\
                    boundary conditions from {bc_timeslice.time.values[0]}")
        return bc_timeslice


def basis(domain, basis_case, basis_directory = None):
    """
    The basis function reads in the all matching files for the basis case and domain as an xarray Dataset.
    
    Expect filenames of the form:
        [basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/LPDM/basis_functions]/EUROPE/sub_transd_EUROPE_2014.nc

    TODO: More info on options for basis functions.

    Args:
        domain (str) : 
            Domain name. The basis files should be sub-categorised by the domain.
        basis_case (str) : 
            Basis case to read in. Examples of basis cases are "voroni","sub-transd",
            "sub-country_mask","INTEM".
        basis_directory (str, optional) : 
            basis_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by domain.
    
    Returns:
        xarray.Dataset : combined dataset of matching basis functions
    """
    if basis_directory is None:
        basis_directory = join(data_path, 'LPDM/basis_functions/')
        
    file_path = os.path.join(basis_directory,domain,f"{basis_case}_{domain}*.nc")
        
    files = sorted(glob.glob(file_path))
    
    if len(files) == 0:
        raise IOError(f"\nError: Can't find basis function files for domain '{domain}' \
                          and basis_case '{basis_case}' ")

    basis_ds = read_netcdfs(files)

    return basis_ds

def basis_boundary_conditions(domain, basis_case, bc_basis_directory = None):
    """
    The basis_boundary_conditions function reads in all matching files for the boundary conditions 
    basis case and domain as an xarray Dataset.
    
    Expect filesnames of the form:
        [bc_basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/LPDM/bc_basis_directory]/EUROPE/NESW_EUROPE_2013.nc

    TODO: More info on options for basis functions.
    
    Args:
        domain (str) : 
            Domain name. The basis files should be sub-categorised by the domain.
        basis_case (str) : 
            Basis case to read in. Examples of basis cases are "NESW","stratgrad".
        bc_basis_directory (str, optional) : 
            bc_basis_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by domain.
    
    Returns:
        xarray.Datset : combined dataset of matching basis functions
    """
    
    if bc_basis_directory is None:
        bc_basis_directory = join(data_path,'LPDM/bc_basis_functions/')
    
    file_path = os.path.join(bc_basis_directory,domain,f"{basis_case}_{domain}*.nc")
    
    files       = sorted(glob.glob(file_path))
    file_no_acc = [ff for ff in files if not os.access(ff, os.R_OK)]
    files       = [ff for ff in files if os.access(ff, os.R_OK)]
    if len(file_no_acc)>0:
        print('Warning: unable to read all boundary conditions basis function files which match this criteria:')
        [print(ff) for ff in file_no_acc]

    if len(files) == 0:
        raise IOError("\nError: Can't find boundary condition basis function files for domain '{0}' "
              "and basis_case '{1}' ".format(domain,basis_case))

    basis_ds = read_netcdfs(files)

    return basis_ds

def indexesMatch(dsa, dsb):
    """
    Check if two datasets need to be reindexed_like for combine_datasets
    
    Args:
        dsa (xarray.Dataset) : 
            First dataset to check
        dsb (xarray.Dataset) : 
            Second dataset to check
            
    Returns:
        boolean:
            True if indexes match, False if datasets must be reindexed
    """
    
    commonIndicies  = [key for key in dsa.indexes.keys() if key in dsb.indexes.keys()]
    
    #test if each comon index is the same
    for index in commonIndicies:
        #first check lengths are the same to avoid error in second check
        if not len(dsa.indexes[index])==len(dsb.indexes[index]):
            return False
        
        #check number of values that are not close (testing for equality with floating point)
        if index == "time":
            #for time iverride the default to have ~ second precision
            rtol = 1e-10
        else:
            rtol = 1e-5
        if not np.sum(~np.isclose(dsa.indexes[index].values.astype(float),dsb.indexes[index].values.astype(float), rtol=rtol ))==0:
            return False
        
    return True
        

def combine_datasets(dsa, dsb, method = "ffill", tolerance = None):
    """
    The combine_datasets function merges two datasets and re-indexes to the FIRST dataset.
    If "fp" variable is found within the combined dataset, the "time" values where the "lat","lon"
    dimensions didn't match are removed.
    
    Example:
        ds = combine_datasets(dsa, dsb)

    Args:
        dsa (xarray.Dataset) : 
            First dataset to merge
        dsb (xarray.Dataset) : 
            Second dataset to merge
        method (str, optional) : 
            One of {None, ‘nearest’, ‘pad’/’ffill’, ‘backfill’/’bfill’}
            See xarray.DataArray.reindex_like for list of options and meaning.
            Default = "ffill" (forward fill)
        tolerance (int/float??) : 
            Maximum allowed tolerance between matches.

    Returns:
        xarray.Dataset: 
            Combined dataset indexed to dsa
    """
    # merge the two datasets within a tolerance and remove times that are NaN (i.e. when FPs don't exist)
    
    if not indexesMatch(dsa, dsb):
        dsb_temp = dsb.reindex_like(dsa, method, tolerance = tolerance)
    else:
        dsb_temp = dsb
    
    ds_temp = dsa.merge(dsb_temp)
    if 'fp' in list(ds_temp.keys()):
        flag = np.where(np.isfinite(ds_temp.fp.mean(dim=["lat","lon"]).values))
        ds_temp = ds_temp[dict(time = flag[0])]
    return ds_temp


def align_datasets(ds1, ds2, platform=None, resample_to_ds1=False):
    """
    Slice and resample two datasets to align along time
    
    Args:
        ds1, ds2 (xarray.Dataset) :
            Datasets with time dimension. It is assumed that ds1 is obs data and ds2 is footprint data
            
        platform (str) :
            obs platform used to decide whether to resample
            
        resample_to_ds1 (boolean) :
            Override resampling to coarser resolution and resample to ds1 regardless
    
    Returns:
        2 xarray.dataset with aligned time dimensions
    """
    platform_skip_resample = ("satellite","flask")
    if platform in platform_skip_resample:
        return ds1, ds2
    #lw13938: 12/04/2018 - This should slice the date to the smallest time frame
    # spanned by both the footprint and obs, then resamples the data 
    #using the mean to the one with coarsest median resolution 
    #starting from the sliced start date. 
    ds1_timeperiod = np.nanmedian((ds1.time.data[1:] - ds1.time.data[0:-1]).astype('int64')) 
    ds2_timeperiod = np.nanmedian((ds2.time.data[1:] - ds2.time.data[0:-1]).astype('int64')) 
    ds1_st = ds1.time[0]
    ds1_et = ds1.time[-1]
    ds2_st = ds2.time[0]
    ds2_et = ds2.time[-1]
    if int(ds1_st.data) > int(ds2_st.data):
        start_date = ds1_st
    else:  
        start_date = ds2_st
    if int(ds1_et.data) < int(ds2_et.data):
        end_date = ds1_et
    else:
        end_date = ds2_et
    

    start_s = str(np.round(start_date.data.astype(np.int64)-5e8,-9).astype('datetime64[ns]')) # subtract half a second to ensure lower range covered
    end_s = str(np.round(end_date.data.astype(np.int64)+5e8,-9).astype('datetime64[ns]')) # add half a second to ensure upper range covered
    
    ds1 = ds1.sel(time=slice(start_s,end_s))
    ds2 = ds2.sel(time=slice(start_s,end_s))
  
    
    #only non satellite datasets with different periods need to be resampled
    if not np.isclose(ds1_timeperiod, ds2_timeperiod):
        base = start_date.dt.hour.data + start_date.dt.minute.data/60. + start_date.dt.second.data/3600.
        if (ds1_timeperiod >= ds2_timeperiod) or (resample_to_ds1 == True):
            resample_period = str(round(ds1_timeperiod/3600e9,4))+'H' # rt17603: Added 24/07/2018 - stops pandas frequency error for too many dp.
            ds2 = ds2.resample(indexer={'time':resample_period}, base=base).mean()
        elif ds1_timeperiod < ds2_timeperiod or (resample_to_ds1 == False):
            resample_period = str(round(ds2_timeperiod/3600e9,4))+'H' # rt17603: Added 24/07/2018 - stops pandas frequency error for too many dp.
            ds1 = ds1.resample(indexer={'time':resample_period}, base=base).mean()
    
    return ds1, ds2


def footprints_data_merge(data, domain, met_model = None, load_flux = True, load_bc = True,
                          calc_timeseries = True, calc_bc = True, HiTRes = False,
                          site_modifier = {}, height = None, network = None,
                          emissions_name = None,
                          fp_directory = None,
                          flux_directory = None,
                          bc_directory = None,
                          resample_to_data = False,
                          species_footprint = None,
                          chunks = False,
                          verbose = True):

    """
    Output a dictionary of xarray footprint datasets, that correspond to a given
    dictionary of Pandas dataframes.
    
    Args:
        data (dict)          : Input dictionary of dataframes with the sites as the keys.
                               Should match output from acrg_agage.get_obs() function. For example:
                               data = {"MHD": MHD_dataframe, "TAC": TAC_dataframe}
        domain (str)         : Domain name. The footprint files should be sub-categorised by the domain.
        met_model (dict or str)     : Met model used to run NAME for each site in data
                               Default is None and implies the standard global met model used for all sites
                               met_model = 'UKV' will work if all sites are the same 'UKV'
                               {"MHD":'UKV', "JFJ":None} or {"MHD":'UKV', "TAC":"UKV"}
                               {"MHD":None, "TAC":None} is equivalent to met_model = None
                               {"MHD":'UKV', "TAC":'UKV'} is equivalent to met_model = 'UKV'
        load_flux (bool)     : True includes fluxes in output, False does not. Default True.
        load_bc (bool)       : True includes boundary conditions in output, False does not. Default True.
        calc_timeseries (bool) : True calculates modelled mole fractions for each site using fluxes, False does not. Default True.
        calc_bc (bool)       : True calculates modelled baseline for each site using boundary conditions, False does not. Default True.
        HiTRes (bool)        : Set to True to include HiTRes footprints in output. Default False.
        site_modifier        : An optional site modifier dictionary is used that maps the site name in the
                               obs file to the site name in the footprint file, if they are different. This
                               is useful for example if the same site FPs are run with a different met and 
                               they are named slightly differently from the obs file. E.g.
                               site_modifier = {"DJI":"DJI-SAM"} - station called DJI, FPs called DJI-SAM
        height (dict)        : Height related to input data. Should be a dictionary with
                               {site: height} key:value pairs. Can be found from acrg_sites_info.json
        network (str/list)   : Network associated with the site. This is only needed if there are multiple 
                               network options. This is used to extract site information 
                               e.g. height from acrg_site_info.json. 
                               Can be specified as a string for one network for all sites or as a list for
                               multiple networks.
        emissions_name (dict): Allows emissions files with filenames that are longer than just the species name
                               to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
                               with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
                               multiple sources can be read in simultaneously if they are added as separate entries to
                               the emissions_name dictionary.
                               If using HiTRes footprints, both the high and low frequency emissions files must be specified
                               in a second dictionary like so: {'anth': {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}}.
                               It is not a problem to have a mixture of sources, with some that use HiTRes footprints and some
                               that don't.
        fp_directory         : fp_directory must be a dictionary of the form 
                               fp_directory = {"integrated":PATH_TO_INTEGRATED_FP, 
                                               "HiTRes":PATH_TO_HIGHTRES_FP}
                               if the high time resolution footprints are used (HiTRes = True)
                               otherwise can be a single string if only integrated FPs are used and 
                               non-default.
        flux_directory (str) : flux_directory can be specified if files are not in the default directory. 
                               Must point to a directory which contains subfolders organized by domain.
                               (optional)
        bc_directory (str)   : Same sytax as flux_directory (optional)
        resample_to_data (bool) : If set to True, the footprints are resampled to the data time series.
                                  If set to False, the data is resampled to the data
                                  If set to None (default), then the footprints and data are sampled to the 
                                  coarset resolution of the two.
                                  (optional)
        chunks (dict)
            size of chunks for each dimension of fp
            e.g. {'lat': 50, 'lon': 50}
            opens dataset with dask, such that it is opened 'lazily'
            and all of the data is not loaded into memory
            defaults to None - dataset is opened with out dask
                                    
    
    Returns:
        Dictionary of the form {"MHD": MHD_xarray_dataset, "TAC": TAC_xarray_dataset, ".flux": dictionary_of_flux_datasets, ".bc": boundary_conditions_dataset}:
            combined dataset for each site
    """

    sites = [key for key in data]
    
    species_list = []
    for site in sites:
        species_list += [item.species for item in data[site]]
    if not all(s==species_list[0] for s in species_list):
        raise Exception("Species do not match in for all measurements")
    else:
        species = species_list[0]
        
    species_footprint = species if species_footprint is None else species_footprint
    if species_footprint!=species:
        print(f'Finding footprints files for {species_footprint}')

    if load_flux:
        if emissions_name is not None:
            if type(emissions_name) != dict:
                print("emissions_name should be a dictionary: {source_name: emissions_file_identifier}.\
                      Setting load_flux to False.")
                load_flux=False   
        else:
            emissions_name = {'all':species}    


    # Output array
    fp_and_data = {}

    #empty variables to fill with earliest start and latest end dates
    flux_bc_start = None
    flux_bc_end = None    

    

    for i, site in enumerate(sites):

        site_ds_list = []

        for site_ds in data[site]:

            if network is None:
                network_site = list(site_info[site].keys())[0]
            elif not isinstance(network,str):
                network_site = network[i]
            else:
                network_site = network

            # Dataframe for this site            
            # Get time range
            df_start = pd.to_datetime(min(site_ds.time.values)) #.to_pydatetime()
            start = dt.datetime(df_start.year, df_start.month, 1, 0, 0)

            df_end = pd.to_datetime(max(site_ds.time.values)) #max(site_df.index).to_pydatetime()
            month_days = calendar.monthrange(df_end.year, df_end.month)[1]
            end = dt.datetime(df_end.year, df_end.month, 1, 0, 0) + \
                    dt.timedelta(days = month_days)

            #Get earliest start and latest end dates from all sites for loading fluxes and bcs
            if flux_bc_start == None or flux_bc_start > start:
                flux_bc_start = start
            if flux_bc_end == None or flux_bc_end < end:
                flux_bc_end = end


            ## Convert to dataset
            #site_ds = xr.Dataset.from_dataframe(site_df)
            if site in list(site_modifier.keys()):
                site_modifier_fp = site_modifier[site]
            else:    
                site_modifier_fp = site

            if "platform" in site_info[site][network_site]:
                platform = site_info[site][network_site]["platform"]
            else:
                platform = None
            
            if height is not None:
                if type(height) is not dict:
                    print("Height input needs to be a dictionary with {sitename:height}")
                    return None 
                height_site = height[site] 
            else:
                if platform == "satellite":
                    height_site = site_info[site][network_site]["height_name"][0]
                else:
                    #Find height closest to inlet
                    siteheights = [int(sh[:-4]) for sh in site_info[site][network_site]["height_name"]]
                    wh_height = np.where(abs(np.array(siteheights) - int(site_ds.inlet[:-1])) == np.min(abs(np.array(siteheights) - int(site_ds.inlet[:-1]))))
                    height_site = site_info[site][network_site]["height_name"][wh_height[0][0]] #NB often different to inlet
            
            if type(met_model) is dict:
                met_model_site = met_model[site]
            else:
                met_model_site = met_model
            
            # Get footprints
            site_fp = footprints(site_modifier_fp, met_model = met_model_site, fp_directory = fp_directory, 
                                 start = start, end = end,
                                 domain = domain,
                                 species = species_footprint.lower(),
                                 height = height_site,
                                 network = network_site,
                                 HiTRes = HiTRes,
                                 chunks = chunks,
                                 verbose = verbose)

            mfattrs = [key for key in site_ds.mf.attrs]
            if "units" in mfattrs:
                if is_number(site_ds.mf.attrs["units"]):
                    units = float(site_ds.mf.attrs["units"])
                else: units = None
            else: units = None
        
            if site_fp is not None:
                # If satellite data, check that the max_level in the obs and the max level in the processed FPs are the same
                # Set tolerance tin time to merge footprints and data   
                # This needs to be made more general to 'satellite', 'aircraft' or 'ship'                

                if platform == "satellite":
                #if "GOSAT" in site.upper():
                    ml_obs = site_ds.max_level
                    ml_fp = site_fp.max_level
                    tolerance = 60e9 # footprints must match data with this tolerance in [ns]
                    if ml_obs != ml_fp:
                        print("ERROR: MAX LEVEL OF SAT OBS DOES NOT EQUAL MAX LEVEL IN FP")
                        print("max_level_fp =",ml_fp)
                        print("max_level_obs =",ml_obs)
                        #return None
                elif "GAUGE-FERRY" in site.upper():
                    tolerance = '5min'
                elif "GAUGE-FAAM" in site.upper():
                    tolerance = '1min'    
                else:
                    tolerance = None

                #gets number of unsorted times in time dimensions, sorting is expensive this is cheap
                if np.sum(np.diff(site_fp.time.values.astype(float))<0) > 0:
                    site_fp = site_fp.sortby("time")

                site_ds, site_fp = align_datasets(site_ds, site_fp, platform=platform, resample_to_ds1=resample_to_data)
                
                site_ds = combine_datasets(site_ds, site_fp,
                                           method = "ffill",
                                           tolerance = tolerance)

                #transpose to keep time in the last dimension position in case it has been moved in resample
                expected_dim_order = ['height','lat','lon','lev','time','H_back']
                for d in expected_dim_order[:]:
                    if d not in list(site_ds.dims.keys()):
                        expected_dim_order.remove(d)
                site_ds = site_ds.transpose(*expected_dim_order)

                # If units are specified, multiply by scaling factor
                if units:
                    site_ds.update({'fp' : (site_ds.fp.dims, old_div(site_ds.fp, units))})
                    if HiTRes:
                        site_ds.update({'fp_HiTRes' : (site_ds.fp_HiTRes.dims, 
                                                       old_div(site_ds.fp_HiTRes, units))})

                site_ds_list += [site_ds]
    
        fp_and_data[site] = xr.merge(site_ds_list)

    if load_flux:

        flux_dict = {} 
        basestring = (str, bytes)    
        for source in list(emissions_name.keys()):

            if isinstance(emissions_name[source], basestring):
                flux_dict[source] = flux(domain, emissions_name[source], start=flux_bc_start, end=flux_bc_end,
                                         flux_directory=flux_directory, verbose=verbose)

            elif isinstance(emissions_name[source], dict):
                if HiTRes == False:
                    print("HiTRes is set to False and a dictionary has been found as the emissions_name dictionary value\
                          for source %s. Either enter your emissions names as separate entries in the emissions_name\
                          dictionary or turn HiTRes to True to use the two emissions files together with HiTRes footprints." %source)
                    #return None
                else:
                    flux_dict[source] = flux_for_HiTRes(domain, emissions_name[source], start=flux_bc_start,
                                                        end=flux_bc_end, flux_directory=flux_directory, verbose=verbose)

        fp_and_data['.flux'] = flux_dict
            
        
    if load_bc:       
        bc = boundary_conditions(domain, species, start=flux_bc_start, end=flux_bc_end, bc_directory=bc_directory)
        if units:
            fp_and_data['.bc'] = old_div(bc, units)               
        else:
            fp_and_data['.bc'] = bc
            
    # Calculate model time series, if required
    if calc_timeseries:
        fp_and_data = add_timeseries(fp_and_data, load_flux, verbose=verbose)
        
    # Calculate boundary conditions, if required         
    if calc_bc:
        fp_and_data = add_bc(fp_and_data, load_bc, species)
    
    #Add "." attributes manually to be back-compatible
    fp_and_data[".species"] = species
    if platform != "satellite":
        scales = {}
        for site in sites:
            scales_list = []
            for site_ds in data[site]:
                scales_list += [site_ds.scale]
                if not all(s==scales_list[0] for s in scales_list):
                    rt = []
                    for i in scales_list:
                        if isinstance(i,list): rt.extend(flatten(i))
                    else: 
                        rt.append(i)
                    scales[site] = rt
                else:
                    scales[site] = scales_list[0]
        fp_and_data[".scales"] = scales
    if units:
        fp_and_data[".units"] = units
  
    return fp_and_data


def add_timeseries(fp_and_data, load_flux, verbose=True):
    """
    Add timeseries mole fraction values in footprint_data_merge
    
    Args:
        fp_and_data [dict]:
            output created during footprint_data_merge
        load_flux [boolean]:
            whether the flux was loaded
    """
    if load_flux == False:
        print("Can't get modelled mole fraction timeseries because load_flux is set to False.")
    else:
        sites = [key for key in list(fp_and_data.keys()) if key[0] != '.']
        sources = list(fp_and_data['.flux'].keys())
        for site in sites:                    
            for source in sources:
                if type(fp_and_data['.flux'][source]) == dict:
                    # work out the lowest resolution of the flux and the fp, to use estimating the timeseries
                    res = [(fp_and_data['.flux'][source].time[1] - fp_and_data['.flux'][source].time[0]).values.astype('timedelta64[h]').astype(int),
                           (fp_and_data[site].time[1] - fp_and_data[site].time[0]).values.astype('timedelta64[h]').astype(int)]
                    res = np.nanmin(res)
                    # estimate the timeseries
                    fp_and_data[site]['mf_mod_'+source] = timeseries_HiTRes(fp_HiTRes_ds = fp_and_data[site],
                                                                            flux_dict = fp_and_data['.flux'][source],
                                                                            output_fpXflux = False,
                                                                            verbose = verbose,
                                                                            time_resolution = res,
                                                                            output_type = 'DataArray')
                else:
                    flux_reindex = fp_and_data['.flux'][source].reindex_like(fp_and_data[site], 'ffill')
                    if source == 'all':
                        fp_and_data[site]['mf_mod'] = xr.DataArray((fp_and_data[site].fp*flux_reindex.flux).sum(["lat", "lon"]), coords = {'time':fp_and_data[site].time})
                    else:
                        fp_and_data[site]['mf_mod_'+source] = xr.DataArray((fp_and_data[site].fp*flux_reindex.flux).sum(["lat", "lon"]), coords = {'time':fp_and_data[site].time})
    return fp_and_data


def add_bc(fp_and_data, load_bc, species):
    """
    Add boundary condition mole fraction values in footprint_data_merge
    Boundary conditions are multipled by any loss (exp(-t/lifetime)) for the species
    
    Args:
        fp_and_data [dict]:
            output created during footprint_data_merge
        load_bc [boolean]:
            whether the boundary conditon was loaded
        species [str]
            name of species in the dataset
    """
    if load_bc == False:
        print("Can't get modelled baseline timeseries because load_bc is set to False.")
    else:
        
        with open(os.path.join(acrg_path,"acrg_species_info.json")) as f:
            species_info=json.load(f)
 
        species_obs = obs.read.synonyms(species, species_info)
    

        sites = [key for key in list(fp_and_data.keys()) if key[0] != '.']
        for site in sites:
            bc_reindex = fp_and_data['.bc'].reindex_like(fp_and_data[site], 'ffill')
            
            if 'lifetime' in species_info[species_obs].keys():
                lifetime = species_info[species_obs]["lifetime"]
                lifetime_hrs_list_or_float = convert.convert_to_hours(lifetime)

                # calculate the lifetime_hrs associated with each time point in fp_and_data
                # this is because lifetime can be a list of monthly values
                
                time_month = fp_and_data[site].time.dt.month
                if type(lifetime_hrs_list_or_float) is list:
                    lifetime_hrs = [lifetime_hrs_list_or_float[item-1] for item in time_month.values]
                else:
                    lifetime_hrs = lifetime_hrs_list_or_float
                                
                loss_n = np.exp(-1*fp_and_data[site].mean_age_particles_n/lifetime_hrs).rename('loss_n')
                loss_e = np.exp(-1*fp_and_data[site].mean_age_particles_e/lifetime_hrs).rename('loss_e')
                loss_s = np.exp(-1*fp_and_data[site].mean_age_particles_s/lifetime_hrs).rename('loss_s')
                loss_w = np.exp(-1*fp_and_data[site].mean_age_particles_w/lifetime_hrs).rename('loss_w')
            else:
                loss_n = 1
                loss_e = 1
                loss_s = 1
                loss_w = 1

            fp_and_data[site]['bc'] = (fp_and_data[site].particle_locations_n*bc_reindex.vmr_n*loss_n).sum(["height", "lon"]) + \
                                        (fp_and_data[site].particle_locations_e*bc_reindex.vmr_e*loss_e).sum(["height", "lat"]) + \
                                        (fp_and_data[site].particle_locations_s*bc_reindex.vmr_s*loss_s).sum(["height", "lon"]) + \
                                        (fp_and_data[site].particle_locations_w*bc_reindex.vmr_w*loss_w).sum(["height", "lat"])
    return fp_and_data


def fp_sensitivity(fp_and_data, domain, basis_case,
                   basis_directory = None, verbose=True):
    """
    The fp_sensitivity function adds a sensitivity matrix, H, to each site xarray dataframe in fp_and_data.

    Basis function data in an array: lat, lon, no. regions. In each 'region'
    element of array there is a lt lon grid with 1 in region and 0 outside region.
    
    Region numbering must start from 1
    
    Args:
        fp_and_data (dict)    : Output from footprints_data_merge() function. Dictionary of datasets.
        domain (str)          : Domain name. The footprint files should be sub-categorised by the domain.
        basis_case            : Basis case to read in. Examples of basis cases are "NESW","stratgrad".
                                String if only one basis case is required. Dict if there are multiple
                                sources that require separate basis cases. In which case, keys in dict should
                                reflect keys in emissions_name dict used in fp_data_merge.
        basis_directory (str) : basis_directory can be specified if files are not in the default 
                                directory. Must point to a directory which contains subfolders organized 
                                by domain. (optional)
    
    Returns:
        dict (xarray.Dataset) : Same format as fp_and_data with sensitivity matrix and basis function grid added.
    """    
    
    sites = [key for key in list(fp_and_data.keys()) if key[0] != '.']
    
    flux_sources = list(fp_and_data['.flux'].keys())
    
    if type(basis_case) is not dict:
        if len(flux_sources) == 1:
            basis_case = {flux_sources[0]:basis_case}
        else:
            basis_case = {'all':basis_case}
    
    if len(list(basis_case.keys())) != len(flux_sources):
        if len(list(basis_case.keys())) == 1:
            print("Using %s as the basis case for all sources" %basis_case[list(basis_case.keys())[0]])
        else:
            print("There should either only be one basis_case, or it should be a dictionary the same length\
                  as the number of sources.")
            return None
    
    
    for site in sites:
        
        for si, source in enumerate(flux_sources):
        
            if source in list(basis_case.keys()):
                basis_func = basis(domain = domain, basis_case = basis_case[source], basis_directory = basis_directory)
            else:
                basis_func = basis(domain = domain, basis_case = basis_case['all'], basis_directory = basis_directory)
            
            if type(fp_and_data['.flux'][source]) == dict:
                if 'fp_HiTRes' in list(fp_and_data[site].keys()):
                    site_bf = xr.Dataset({"fp_HiTRes":fp_and_data[site]["fp_HiTRes"],
                                          "fp":fp_and_data[site]["fp"]})
                    res = [(fp_and_data['.flux'][source].time[1] - fp_and_data['.flux'][source].time[0]).values.astype('timedelta64[h]').astype(int),
                           (site_bf.time[1] - site_bf.time[0]).values.astype('timedelta64[h]').astype(int)]
                    res = np.nanmin(res)
                    H_all = timeseries_HiTRes(fp_HiTRes_ds = site_bf, flux_dict = fp_and_data['.flux'][source], output_TS = False,
                                              output_fpXflux = True, output_type = 'DataArray', time_resolution=f'{res}H', verbose = verbose)
                else:
                    print("fp_and_data needs the variable fp_HiTRes to use the emissions dictionary with high_freq and low_freq emissions.")
        
            else:
                site_bf = combine_datasets(fp_and_data[site]["fp"].to_dataset(), fp_and_data['.flux'][source])
                H_all=site_bf.fp*site_bf.flux
            
            H_all_v=H_all.values.reshape((len(site_bf.lat)*len(site_bf.lon),len(site_bf.time)))        
        
        
            if 'region' in list(basis_func.dims.keys()):
            
                if 'time' in basis_func.basis.dims:
                    basis_func = basis_func.isel(time=0)
            
                site_bf = xr.merge([site_bf, basis_func])
                
                H = np.zeros((len(site_bf.region),len(site_bf.time)))
            
                base_v = site_bf.basis.values.reshape((len(site_bf.lat)*len(site_bf.lon), len(site_bf.region)))
            
                for i in range(len(site_bf.region)):
                    H[i,:] = np.sum(H_all_v*base_v[:,i,np.newaxis], axis = 0)
                
                if source == all:
                    if (sys.version_info < (3,0)):
                        region_name = site_bf.region
                    else:
                        region_name = site_bf.region.decode('ascii')
                else:
                    if (sys.version_info < (3,0)):
                        region_name = [source+'-'+reg for reg in site_bf.region.values]
                    else:
                        region_name = [source+'-'+reg.decode('ascii') for reg in site_bf.region.values]

                sensitivity = xr.DataArray(H, 
                                             coords=[('region', region_name), 
                                                     ('time', fp_and_data[site].coords['time'])])
        
            else:
                print("Warning: Using basis functions without a region dimension may be deprecated shortly.")
        
                site_bf = combine_datasets(site_bf,basis_func, method='ffill')
 
                H = np.zeros((int(np.max(site_bf.basis)),len(site_bf.time)))

                basis_scale = xr.Dataset({'basis_scale': (['lat','lon','time'],

                                                    np.zeros(np.shape(site_bf.basis)))},
                                       coords = site_bf.coords)
                site_bf = site_bf.merge(basis_scale)

                base_v = np.ravel(site_bf.basis.values[:,:,0])
                for i in range(int(np.max(site_bf.basis))):
                    wh_ri = np.where(base_v == i+1)
                    H[i,:]=np.sum(H_all_v[wh_ri[0],:], axis = 0)      
                  
                if source == all:
                    region_name = list(range(1,np.max(site_bf.basis.values)+1))
                else:
                    region_name = [source+'-'+str(reg) for reg in range(1,int(np.max(site_bf.basis.values)+1))]

                sensitivity = xr.DataArray(H, 
                                             coords=[('region', region_name), 
                                                     ('time', fp_and_data[site].coords['time'])])
                                     
            if si == 0:
                concat_sensitivity = sensitivity
            else:
                concat_sensitivity = xr.concat((concat_sensitivity,sensitivity), dim='region')
            
            sub_basis_cases = 0
            if basis_case[source].startswith('sub'):
                """
                To genrate sub_lon and sub_lat grids basis case must start with 'sub'
                e.g.
                'sub-transd', 'sub_transd', sub-intem' will work
                'transd' or 'transd-sub' won't work
                """
                sub_basis_cases += 1
                if sub_basis_cases > 1:
                    print("Can currently only use a sub basis case for one source. Skipping...")
                else:
                    sub_fp_temp = site_bf.fp.sel(lon=site_bf.sub_lon, lat=site_bf.sub_lat,
                                                 method="nearest") 
                    sub_fp = xr.Dataset({'sub_fp': (['sub_lat','sub_lon','time'], sub_fp_temp)},
                                           coords = {'sub_lat': (site_bf.coords['sub_lat']),
                                                     'sub_lon': (site_bf.coords['sub_lon']),
                                                     'time' : (fp_and_data[site].coords['time'])})
                                
                    sub_H_temp = H_all.sel(lon=site_bf.sub_lon, lat=site_bf.sub_lat,
                                           method="nearest")                             
                    sub_H = xr.Dataset({'sub_H': (['sub_lat','sub_lon','time'], sub_H_temp)},
                                          coords = {'sub_lat': (site_bf.coords['sub_lat']),
                                                    'sub_lon': (site_bf.coords['sub_lon']),
                                                    'time' : (fp_and_data[site].coords['time'])},
                                          attrs = {'flux_source_used_to_create_sub_H':source})
       
                    fp_and_data[site] = fp_and_data[site].merge(sub_fp)
                    fp_and_data[site] = fp_and_data[site].merge(sub_H)
            
        fp_and_data[site]['H'] = concat_sensitivity                             
        fp_and_data['.basis'] = site_bf.basis[:,:,0]
                    
    return fp_and_data


def bc_sensitivity(fp_and_data, domain, basis_case, bc_basis_directory = None):

    """
    The bc_sensitivity adds H_bc to the sensitivity matrix, to each site xarray dataframe in fp_and_data.
    
    Args:
        fp_and_data (dict)       : Output from footprints_data_merge() function. Dictionary of datasets.
        domain (str)             : Domain name. The footprint files should be sub-categorised by the domain.
        basis_case (str)         : Basis case to read in. Examples of basis cases are "NESW","stratgrad".
        bc_basis_directory (str) : bc_basis_directory can be specified if files are not in the default 
                                   directory. Must point to a directory which contains subfolders organized 
                                   by domain. (optional)
    
    Returns:
        dict (xarray.Dataset) : Same format as fp_and_data with sensitivity matrix added.
    """    
    
    
    sites = [key for key in list(fp_and_data.keys()) if key[0] != '.']

    basis_func = basis_boundary_conditions(domain = domain,
                                           basis_case = basis_case, bc_basis_directory=bc_basis_directory)
# sort basis_func into time order    
    ind = basis_func.time.argsort()                                        
    timenew = basis_func.time[ind]
    basis_func = basis_func.reindex({"time":timenew})
    
    species = fp_and_data[".species"]
    with open(os.path.join(acrg_path,"acrg_species_info.json")) as f:
        species_info=json.load(f)
    species = obs.read.synonyms(species, species_info)
            
    for site in sites:
        # compute any chemical loss to the BCs, use lifetime or else set loss to 1 (no loss)
        if 'lifetime' in species_info[species].keys():
            lifetime = species_info[species]["lifetime"]
            lifetime_hrs_list_or_float = convert.convert_to_hours(lifetime)

            # calculate the lifetime_hrs associated with each time point in fp_and_data
            # this is because lifetime can be a list of monthly values

            time_month = fp_and_data[site].time.dt.month
            if type(lifetime_hrs_list_or_float) is list:
                lifetime_hrs = [lifetime_hrs_list_or_float[item-1] for item in time_month.values]
            else:
                lifetime_hrs = lifetime_hrs_list_or_float

            loss_n = np.exp(-1*fp_and_data[site].mean_age_particles_n/lifetime_hrs).rename('loss_n')
            loss_e = np.exp(-1*fp_and_data[site].mean_age_particles_e/lifetime_hrs).rename('loss_e')
            loss_s = np.exp(-1*fp_and_data[site].mean_age_particles_s/lifetime_hrs).rename('loss_s')
            loss_w = np.exp(-1*fp_and_data[site].mean_age_particles_w/lifetime_hrs).rename('loss_w')
        else:
            loss_n = fp_and_data[site].particle_locations_n.copy()
            loss_e = fp_and_data[site].particle_locations_e.copy()
            loss_s = fp_and_data[site].particle_locations_s.copy()
            loss_w = fp_and_data[site].particle_locations_w.copy()
            loss_n[:]=1
            loss_e[:]=1
            loss_s[:]=1
            loss_w[:]=1

        DS_particle_loc = xr.Dataset({"particle_locations_n":fp_and_data[site]["particle_locations_n"],
                                "particle_locations_e":fp_and_data[site]["particle_locations_e"],
                                "particle_locations_s":fp_and_data[site]["particle_locations_s"],
                                "particle_locations_w":fp_and_data[site]["particle_locations_w"],
                                "loss_n":loss_n,
                                "loss_e":loss_e,
                                "loss_s":loss_s,
                                "loss_w":loss_w})                                     
#                                 "bc":fp_and_data[site]["bc"]})

        DS_temp = combine_datasets(DS_particle_loc, fp_and_data[".bc"], method='ffill')
                        
        DS = combine_datasets(DS_temp, basis_func, method='ffill')                                    
                               
        DS = DS.transpose('height','lat','lon','region','time')

        part_loc = np.hstack([DS.particle_locations_n,
                                DS.particle_locations_e, 
                                DS.particle_locations_s,
                                DS.particle_locations_w])
        
        loss = np.hstack([DS.loss_n,
                                DS.loss_e,
                                DS.loss_s,
                                DS.loss_w])
        
        vmr_ed = np.hstack([DS.vmr_n,
                           DS.vmr_e,
                           DS.vmr_s,
                           DS.vmr_w])
        
        bf = np.hstack([DS.bc_basis_n,
                        DS.bc_basis_e,
                        DS.bc_basis_s,
                        DS.bc_basis_w])
        
        H_bc = np.zeros((len(DS.coords['region']),len(DS["particle_locations_n"]["time"])))
        
        for i in range(len(DS.coords['region'])):
            reg = bf[:,:,i,:]
            H_bc[i,:] = np.sum((part_loc*loss*vmr_ed*reg), axis=(0,1))
        
        sensitivity = xr.Dataset({'H_bc': (['region_bc','time'], H_bc)},
                                    coords = {'region_bc': (DS.coords['region'].values),
                                              'time' : (DS.coords['time'])})

        fp_and_data[site] = fp_and_data[site].merge(sensitivity)
    
    return fp_and_data


def filtering(datasets_in, filters, keep_missing=False):
    """
    Applies filtering (in time dimension) to entire dataset.
    Filters supplied in a list and then applied in order. For example if you wanted a daily, daytime 
    average, you could do this:
    
        datasets_dictionary = filtering(datasets_dictionary, 
                                    ["daytime", "daily_median"])
    
    The order of the filters reflects the order they are applied, so for 
    instance when applying the "daily_median" filter if you only wanted
    to look at daytime values the filters list should be 
    ["daytime","daily_median"]                

    Args:
        datasets_in         : Output from footprints_data_merge(). Dictionary of datasets.
        filters (list)      : Which filters to apply to the datasets. 
                              All options are:
                                 "daytime"           : selects data between 1100 and 1500 local solar time
                                 "daytime9to5"       : selects data between 0900 and 1700 local solar time
                                 "nighttime"         : Only b/w 23:00 - 03:00 inclusive
                                 "noon"              : Only 12:00 fp and obs used
                                 "daily_median"      : calculates the daily median
                                 "pblh_gt_threshold" : 
                                 "local_influence"   : Only keep times when localness is low
                                 "six_hr_mean"       :
                                 "local_lapse"       :
        keep_missing (bool) : Whether to reindex to retain missing data.
    
    Returns:
       Same format as datasets_in : Datasets with filters applied. 
    """

    if type(filters) is not list:
        filters = [filters]

    datasets = datasets_in.copy()

    def local_solar_time(dataset):
        """
        Returns hour of day as a function of local solar time
        relative to the Greenwich Meridian. 
        """
        sitelon = dataset.release_lon.values[0]
        # convert lon to [-180,180], so time offset is negative west of 0 degrees
        if sitelon > 180:
            sitelon = sitelon - 360.
        dataset["time"] = dataset.time + pd.Timedelta(minutes=float(24*60*sitelon/360.))
        hours = dataset.time.to_pandas().index.hour
        return hours
    
    def local_ratio(dataset):
        """
        Calculates the local ratio in the surrounding grid cells
        """
        release_lons = dataset.release_lon[0].values
        release_lats = dataset.release_lat[0].values
        dlon = dataset.lon[1].values - dataset.lon[0].values
        dlat = dataset.lat[1].values-dataset.lat[0].values
        local_sum=np.zeros((len(dataset.mf)))

        for ti in range(len(dataset.mf)):
            release_lon=dataset.release_lon[ti].values
            release_lat=dataset.release_lat[ti].values
            wh_rlon = np.where(abs(dataset.lon.values-release_lon) < dlon/2.)
            wh_rlat = np.where(abs(dataset.lat.values-release_lat) < dlat/2.)
            if np.any(wh_rlon[0]) and np.any(wh_rlat[0]):
                local_sum[ti] = old_div(np.sum(dataset.fp[
                        wh_rlat[0][0]-2:wh_rlat[0][0]+3,wh_rlon[0][0]-2:wh_rlon[0][0]+3,ti].values),np.sum(
                        dataset.fp[:,:,ti].values))  
            else:
                local_sum[ti] = 0.0 
        
        return local_sum
    
    
    
    # Filter functions
    def daily_median(dataset, keep_missing=False):
        """ Calculate daily median """
        if keep_missing:
            return dataset.resample(indexer={'time':"1D"}).median()
        else:
            return dataset.resample(indexer={'time':"1D"}).median().dropna(dim="time")
        
    def six_hr_mean(dataset, keep_missing=False):
        """ Calculate six-hour median """
        if keep_missing:
            return dataset.resample(indexer={'time':"6H"}).mean()
        else:
            return dataset.resample(indexer={'time':"6H"}).mean().dropna(dim="time")
    

    def daytime(dataset, site,keep_missing=False):
        """ Subset during daytime hours (11:00-15:00) """
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 11 and h <= 15]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
        
    def daytime9to5(dataset, site,keep_missing=False):
        """ Subset during daytime hours (9:00-17:00) """
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 9 and h <= 17]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def nighttime(dataset, site,keep_missing=False):
        """ Subset during nighttime hours (23:00 - 03:00) """
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 23 or h <= 3]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def noon(dataset, site,keep_missing=False):
        """ Select only 12pm data """
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h == 12]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)] 
        
            
    def local_influence(dataset,site, keep_missing=False):
        """
        Subset for times when local influence is below threshold.       
        Local influence expressed as a fraction of the sum of entire footprint domain.
        """        
        if not dataset.filter_by_attrs(standard_name="local_ratio"):
            lr = local_ratio(dataset)
        else:
            lr = dataset.local_ratio
            
        pc = 0.1
        ti = [i for i, local_ratio in enumerate(lr) if local_ratio <= pc]
        if keep_missing is True: 
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')

            dataarray_temp = mf_data_array[dict(time = ti)]   

            mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})

            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
                     
        
    filtering_functions={"daily_median":daily_median,
                         "daytime":daytime,
                         "daytime9to5":daytime9to5,
                         "nighttime":nighttime,
                         "noon":noon,
                         "local_influence":local_influence,
                         "six_hr_mean":six_hr_mean}

    # Get list of sites
    sites = [key for key in list(datasets.keys()) if key[0] != '.']
    
    # Do filtering
    for site in sites:
    
            for filt in filters:
                if filt == "daily_median" or filt == "six_hr_mean":
                    datasets[site] = filtering_functions[filt](datasets[site], keep_missing=keep_missing)
                else:
                    datasets[site] = filtering_functions[filt](datasets[site], site, keep_missing=keep_missing)

    return datasets


def plot(fp_data, date, out_filename=None, out_format = 'pdf',
         lon_range=None, lat_range=None, log_range = [5., 9.], plot_borders = False,
         zoom = False, colormap = 'YlGnBu', tolerance = None, interpolate = False, dpi = 300,
         figsize=None, nlevels=256):
    """
    Plot footprint for a given timestamp.
    
    Args:
        fp_data (dict): 
            Output of footprints_data_merge(). Dictionary of xarray datasets containing   
            footprints and other variables
        date (str): 
            Almost any time format should work (datetime object, string, etc). An example
            time format is '2014-01-01 00:00'. Footprints from all sites in dictionary that 
            have time indices nearest to the specified time will be plotted.
        out_filename (str, optional):
            Full path to filname to save figure
        out_format (str, optional):
            Format to save figure (e.g., png or pdf)
        lon_range (list, optional): 
            list of min and max longitudes [min, max] to plot
        lat_range (list, optional): 
            list of min and max latitudes [min, max] to plot
        log_range (list, optional): 
            list of min and max LOG10(footprints) for color scale       
        plot_borders (bool, optional) :
            Plot country borders. Default = False.
        zoom (bool, optional): 
            True will plot a zoomed map (+/- 10 degrees around all site in fp_data)
        colormap (str, optional): 
            Color map to use for contour plot 
            (https://matplotlib.org/examples/color/colormaps_reference.html)
        tolerance (float or str, optional): 
            Maximum distance between date specified and closest available footprint
            Default is in nanosec if a float or in any units (e.g., '1H') if a string
        interpolate (bool, optional): 
            If True, interpolates footprint between the nearest two footprints to the date specified
            If False, uses the nearest footprint avaialble
        dpi (int, optional):
            Dots per square inch resolution to save image format such as png
        figsize (tuple, optional):
            Specify figure size as width, height in inches. e.g. (12,9). Default = None.
        nlevels (int):
            Number of levels in contour plot.
            
    Returns
        None
        
        If out_filename == None:
            produces as interactive plot
        Else:
            saves plot as specified by the full path in out_filename and format in out_format
    """
    
    def fp_nearest(fp, tolerance = None):
        return fp.reindex_like( \
                            xr.Dataset(coords = {"time": [date]}),
                            method = "nearest",
                            tolerance = tolerance)

    
    # Looks for nearest time point aviable in footprint   
    date = np.datetime64(convert.reftime(date))

    # Get sites
    sites = [key for key in list(fp_data.keys()) if key[0] != '.']
    
    # Find lat and lon range of the footprints
    if lon_range is None:
        lons = fp_data[sites[0]].lon.values
    else:
        lons_all = fp_data[sites[0]].lon.values
        indx = np.where((lons_all >= lon_range[0]) & (lons_all <= lon_range[-1]))[0]
        lons = lons_all[indx]
        indx = np.where((lons_all < lon_range[0]) | (lons_all > lon_range[-1]))[0]
                
    if lat_range is None:
        lats = fp_data[sites[0]].lat.values
    else:
        lats_all = fp_data[sites[0]].lat.values
        indy = np.where((lats_all >= lat_range[0]) & (lats_all <= lat_range[-1]))[0]
        lats = lats_all[indy]
        indy = np.where((lats_all < lat_range[0]) | (lats_all > lat_range[-1]))[0]       
    
#    Zoom in. Get min and max release lat lons to zoom the map around the data (+/- 10 degrees)
    if zoom:
        release_lons = [fp_data[key].release_lon for key in list(fp_data.keys()) if key[0] != '.']     
        release_lats = [fp_data[key].release_lat for key in list(fp_data.keys()) if key[0] != '.']
        
        release_lons = [y for x in release_lons for y in x]
        release_lats = [y for x in release_lats for y in x]
        lon_range = [np.min(release_lons)-10, np.max(release_lons)+10]
        lat_range = [np.min(release_lats)-10, np.max(release_lats)+10]

        lons_all = fp_data[sites[0]].lon.values
        lats_all = fp_data[sites[0]].lat.values
        indx = np.where((lons_all >= lon_range[0]) & (lons_all <= lon_range[-1]))[0]
        lons = lons_all[indx]
        indx = np.where((lons_all < lon_range[0]) | (lons_all > lon_range[-1]))[0]
        indy = np.where((lats_all >= lat_range[0]) & (lats_all <= lat_range[-1]))[0]
        lats = lats_all[indy]
        indy = np.where((lats_all < lat_range[0]) | (lats_all > lat_range[-1]))[0]  
    
    if figsize:
        if isinstance(figsize,tuple) and len(figsize) == 2:
            plt.figure(figsize=figsize)
        elif len(figsize) == 2:
            plt.figure(figsize=(figsize[0],figsize[1]))
        else:
            print("Could not apply figure size: {}. Expected two item tuple.".format(figsize))
    
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=np.median(lons)))
    ax.set_extent([lons[0], lons[-1], lats[0], lats[-1]], crs=ccrs.PlateCarree())
    ax.coastlines()
    
    if plot_borders:
        ax.add_feature(cartopy.feature.BORDERS,linewidth=0.5)
    
    #Calculate color levels
    rp_color = {"SURFACE": "black",
                "SHIP": "purple",
                "AIRCRAFT": "red",
                "SATELLITE": "green"}

            
    levels = MaxNLocator(nbins=nlevels).tick_values(log_range[0], log_range[1])

    # Create dictionaries and arrays    
    release_lon = {}
    release_lat = {}

    data = {}

    data = np.zeros(np.shape(
            fp_data[sites[0]].fp[dict(time = [0])].values.squeeze()))

    
    time = None

    # Generate plot data
    for site in sites:
    
        time = fp_data[site].time.values.squeeze().copy()
        
        if tolerance is None:
            tol = np.median(time[1:] - time[:-1])
        
        if interpolate:

            ti = bisect.bisect(time, date)
            
            if (ti > 0) and (ti < len(time)):

                if np.float(date - time[ti-1]) < np.float(tol):

                    dt = old_div(np.float(date - time[ti-1]),np.float(time[ti] - time[ti-1]))
                    fp_ti_0 = fp_data[site][dict(time = ti-1)].fp.values.squeeze()
                    fp_ti_1 = fp_data[site][dict(time = ti)].fp.values.squeeze()
                    fp_ti = fp_ti_0 + (fp_ti_1 - fp_ti_0)*dt
                    
                    data += np.nan_to_num(fp_ti)
    
                    # Store release location to overplot later
                    if "release_lat" in list(fp_data[site].keys()):
                        lon_0 = fp_data[site][dict(time = ti-1)].release_lon.values.squeeze()
                        lon_1 = fp_data[site][dict(time = ti)].release_lon.values.squeeze()                
                        lat_0 = fp_data[site][dict(time = ti-1)].release_lat.values.squeeze()
                        lat_1 = fp_data[site][dict(time = ti)].release_lat.values.squeeze()
                        release_lon[site] = lon_0 + (lon_1 - lon_0)*dt
                        release_lat[site] = lat_0 + (lat_1 - lat_0)*dt
                    
            else:

                fp_data_ti = fp_nearest(fp_data[site], tolerance = tolerance)
                data += np.nan_to_num(fp_data_ti.fp.values.squeeze())
                # Store release location to overplot later
                if "release_lat" in dir(fp_data_ti):
                    release_lon[site] = fp_data_ti.release_lon.values
                    release_lat[site] = fp_data_ti.release_lat.values

        else:
            fp_data_ti = fp_nearest(fp_data[site], tolerance = tolerance)
            data += np.nan_to_num(fp_data_ti.fp.values.squeeze())

            # Store release location to overplot later
            if "release_lat" in dir(fp_data_ti):
                release_lon[site] = fp_data_ti.release_lon.values
                release_lat[site] = fp_data_ti.release_lat.values

    #Set very small elements to zero
    data = np.log10(data)
    data[np.where(data <  log_range[0])]=np.nan

    # If lat range, lon range or zoom input, then subset data    
    if lon_range or zoom:
        data = np.delete(data, indx, axis=1)
    if lat_range or zoom:
        data = np.delete(data, indy, axis=0)
    
    #Plot footprint
    plt.contourf(lons, lats, data, transform=ccrs.PlateCarree(), cmap = colormap, levels=levels)

    # over-plot release location
    if len(release_lon) > 0:
        for site in sites:
            network = list(site_info[site].keys())[0]
            if site in release_lon:
                if "platform" in site_info[site][network]:
                    color = rp_color[site_info[site][network]["platform"].upper()]
                else:
                    color = rp_color["SURFACE"]
                plt.plot(release_lon[site], release_lat[site], color = color, marker = 'o', markersize=4,  transform=ccrs.PlateCarree())

    plt.title(str(pd.to_datetime(str(date))), fontsize=12)

    cb = plt.colorbar(orientation='horizontal', pad=0.05)
   
    tick_locator = ticker.MaxNLocator(nbins=np.ceil(log_range[1] - log_range[0]).astype(int))
    cb.locator = tick_locator
    cb.update_ticks()
 
    cb.set_label('log$_{10}$( (nmol/mol) / (mol/m$^2$/s) )', 
                fontsize=12)
    cb.ax.tick_params(labelsize=12) 
    
    if out_filename is not None:
        plt.savefig(out_filename, dpi = dpi, format = out_format)
        plt.close()
    else:
        plt.show()


def plot_particle_location(fp_data, date, particle_direction = 'nw', out_filename=None,
                           out_format = 'pdf', tolerance = None, log_range = [5., 9.], 
                           colormap_fp = 'inferno_r', colormap_part = 'GnBu',
                           particle_clevs = [0., 0.009, 0.001], dpi = 300,
                           figsize = None):

    """
    3D plot showing the footprint and particle exit locations for a given timestamp.
    
    Args:
        fp_data (dict): 
            Dictionary of xarray datasets containing footprints and other variables. 
            Currently only one site will be plotted.
        date (str): 
            Almost any time format should work (datetime object, string, etc). An example
            time format is '2014-01-01 00:00'. Footprints from all sites in dictionary that 
            have time indices nearest to the specified time will be plotted.
        particle_direction (str):
            One of two options 'nw' = north-west or 'se' = south-west to determine which
            edges of the domain are plotted for particle exit locations
        out_filename (str, optional):
            Full path to filname to save figure
        out_format (str, optional):
            Format to save figure (e.g., png or pdf)
        tolerance (float or str, optional): 
            Maximum distance between date specified and closest available footprint
            Default is in nanosec if a float or in any units (e.g., '1H') if a string
        log_range (list, optional): 
            list of min and max LOG10(footprints) for color scale       
        colormap_fp (str, optional): 
            Color map to use for contour plot of footprint
            (https://matplotlib.org/examples/color/colormaps_reference.html)
        colormap_part (str, optional): 
            Color map to use for contour plot of particle location
            (https://matplotlib.org/examples/color/colormaps_reference.html)                    
        particle_clevs (list, optional):
            List of [min, max, interval] to plot the particle color levels
        dpi (int, optional):
            Resolution to save png
            
            
    Returns
        None
        
        If out_filename == None:
            produces as interactive plot
        Else:
            saves plot as specified by the full path in out_filename and format in out_format
    """
    def fp_nearest(fp, tolerance = None):
        return fp.reindex_like( \
                        xr.Dataset(coords = {"time": [date]}),
                        method = "nearest",
                        tolerance = tolerance)
    
    sites = [key for key in list(fp_data.keys()) if key[0] != '.']
    fp_data = fp_data[sites[0]]

    date = np.datetime64(convert.reftime(date))
        
    fp_data = fp_nearest(fp_data, tolerance = tolerance)
    fp_data = fp_data.squeeze(dim='time')

    lon_range = (fp_data.lon.min().values, fp_data.lon.max().values)
    lat_range = (fp_data.lat.min().values, fp_data.lat.max().values)

    #Set very small elements to zero
    fp_data.where(np.log10(fp_data["fp"]) < log_range[0])
    
    if figsize:
        if len(figsize) == 2:
            if not isinstance(figsize,tuple):
                figsize=(figsize[0],figsize[1])
        else:
            print("Could not apply figure size: {}. Using default: {}".format(figsize,(8,6)))
            figsize = (8,6)
    else:
        figsize = (8,6)
        print("Using default figsize: {}".format(figsize))

    
    figure = plt.figure(figsize=figsize, facecolor='w')
    ax = figure.gca(projection='3d')
    
    ax.set_ylim(lat_range)
    ax.set_xlim(lon_range)
    ax.set_zlim((min(fp_data.height), max(fp_data.height)))

    fpX, fpY = np.meshgrid(fp_data.lon, fp_data.lat)
    levels = np.arange(log_range[0], log_range[1], 0.05)

    plnX, plnY = np.meshgrid(fp_data.lon.values.squeeze(), fp_data.height.values.squeeze())
    plwX, plwY = np.meshgrid(fp_data.lat.values.squeeze(), fp_data.height.values.squeeze())
    pllevs = np.arange(particle_clevs[0], particle_clevs[1], particle_clevs[2])

    if particle_direction == 'nw':
        plfp = ax.contourf(fpX, fpY, np.log10(fp_data["fp"]), levels, offset = 0., cmap = colormap_fp)
        plnvals = fp_data.particle_locations_n.values
        plnvals[np.where(plnvals == 0.)]=np.nan
        plwvals = fp_data.particle_locations_w.values
        plwvals[np.where(plwvals == 0.)]=np.nan
        plpln = ax.contourf(plnvals,plnX, plnY,
                            zdir = 'y', offset = max(fp_data.lat.values), levels = pllevs, cmap = colormap_part)
        ax.text(150,0,23500,"North", color='red', zdir=None)
        plplw = ax.contourf(plwvals, plwX, plwY,
                            zdir = 'x', offset = min(fp_data.lon.values), levels = pllevs, cmap = colormap_part)    
        ax.text(0,0,10000,"West", color='red', zdir=None)
    elif particle_direction == 'se':
        plfp = ax.contourf(fpX, fpY, np.fliplr(np.flipud(np.log10(fp_data["fp"]))), levels, offset = 0., cmap = colormap_fp)
        plnvals = fp_data.particle_locations_s.values
        plnvals[np.where(plnvals == 0.)]=np.nan
        plwvals = fp_data.particle_locations_e.values
        plwvals[np.where(plwvals == 0.)]=np.nan
        plpln = ax.contourf(plnvals,plnX, plnY,
                            zdir = 'y', offset = max(fp_data.lat.values), levels = pllevs, cmap = colormap_part)
        ax.text(150,0,23500,"South", color='red', zdir=None)
        plplw = ax.contourf(plwvals, plwX, plwY,
                            zdir = 'x', offset = min(fp_data.lon.values), levels = pllevs, cmap = colormap_part)    
        ax.text(0,0,10000,"East", color='red', zdir=None)
    ax.view_init(50)

    plt.title(str(pd.to_datetime(str(date))), fontsize = 12, y=1.08)
    plt.xlabel('longitude', color='blue')
    plt.ylabel('latitude', color='blue')
    ax.set_zlabel('altitude (m)', color='blue')
    cb = plt.colorbar(plfp, orientation='horizontal', shrink=0.4)
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label('log$_{10}$( (mol/mol) / (mol/m$^2$/s))', 
             fontsize=12)
    cb.ax.tick_params(labelsize=12) 
    
    if out_filename is not None:
        plt.savefig(out_filename, dpi = dpi, format = out_format)
        plt.close()
    else:
        plt.show()


def animate(fp_data, output_directory, plot_function = "plot", file_label = 'fp', 
            video_os="mac", time_regular = False,        
            lon_range = None, lat_range = None, log_range = [5., 9.],plot_borders=False,
            colormap_fp = 'inferno_r', colormap_part = 'GnBu', zoom = False,
            particle_clevs = [0., 0.009, 0.001], overwrite = True, 
            framerate=10, delete_png=False, ffmpeg_only = False,
            frame_max = None, dpi = 300,figsize=None):

    """
    Animate footprints into a movie.
    
    Args:
        fp_data (dict): 
            Dictionary of xarray datasets containing footprints and other variables
        output_directory (str):
            Full path to directory in which movie and images will be saved
        plot_function (str):
            Either 'plot' or 'plot_particle_location'. If plots are being generated 
            (ffpmeg_only = False), call to this plotting function.
        file_label (str):
            File lavel of animation and and images of the frames 
        video_os (str):
            Operating system to play video - must be either mac or pc
        time_regular (str): 
            Frequency between minumum and maximum to set the time values within fp_data. (e.g., '1H') 
            Set at False to not apply this step (Default = False).
        lon_range (list, optional): 
            list of min and max longitudes [min, max] to plot
        lat_range (list, optional): 
            list of min and max latitudes [min, max] to plot
        log_range (list, optional): 
            list of min and max LOG10(footprints) for color scale       
        plot_borders (bool, optional) :
            Plot country borders. Only applicable to plot_function="plot" Default = False.
        colormap_fp (str, optional): 
            Color map to use for contour plot of footprint
            (https://matplotlib.org/examples/color/colormaps_reference.html)
        colormap_part (str, optional): 
            Color map to use for contour plot of particle location
            (https://matplotlib.org/examples/color/colormaps_reference.html)
        zoom (bool, optional): 
            True will plot a zoomed map (+/- 10 degrees around all site in fp_data)
            if plot_function is 'plot'
        particle_clevs (list, optional):
            List of [min, max, interval] to plot the particle color levels if 
            plot_function is 'plot_particle_location'
        overwrite (bool, optional):
            True will overwrite any existing files
        framerate (int, optional):
            Framerate of animation in frames per second
        delete_png (bool. optional):
            True will delete the images of each frame
        ffmpeg_only (bool, optional):
            True will generate movie from pre-saved images in the output directory 
            named with file_label
        frame_max (int, optional):
            Set the maximum number of frames in the datset to animate. Animation will plot first n 
            frames up to the frame_max. Useful for testing.
        dpi (int, optional):
            Dots per square inch resolution for each image generate as for example png 
        figsize (tuple, optional):
            Specify figure size in width, height in inches. e.g. (12,9)
            
    Returns
        None
        
        If out_filename == None:
            produces as interactive plot
        Else:
            saves plot as specified by the full path in out_filename and format in out_format
    """
    
    def time_unique(fp_data, time_regular = False):
    
        sites = [key for key in list(fp_data.keys()) if key[0] != '.']
        
        time_array = fp_data[sites[0]].time
        time_array.name = "times"
        time = time_array.to_dataset()
        if len(sites) > 1:
            for site in sites[1:]:
                time.merge(fp_data[site].time.to_dataset(), inplace = True)
        
        time = time.time.to_pandas().index.to_pydatetime()
        if time_regular is not False:
            time = pd.date_range(min(time), max(time), freq = time_regular)
        
        return time
  
        """
        The time_unique function creates one set of time entries within fp_data (dictionary of datasets)
        merging across multiple sites if necessary.
        Time can also be reindexed to be evenly spaced between the minimum and maximum values if time_regular
        if specified.
        
        Args:
            fp_data (dict)     : 
                Output from footprints_data_merge(). Dictionary of datasets.
            time_regular (str) : 
                Frequency between minumum and maximum to set the time values within fp_data. 
                Set at False to not apply this step (Default = False).
        Returns:
            xarray.Dataset : Time values extracted for this first site within fp_data
        """    
    
    if ffmpeg_only is False:
        
        # Find unique times
        times = time_unique(fp_data,
                            time_regular = time_regular)
        
        if frame_max:
            times = times[:min([frame_max, len(times)])]

        # Start progress bar
        #pbar = ProgressBar(maxval=len(times)).start()

        # Plot each timestep
        for ti, t in enumerate(times):
            
            fname=os.path.join(output_directory, 
                               file_label + '_' + str(ti).zfill(5) + '.png')
            
            if len(glob.glob(fname)) == 0 or overwrite == True:
                if plot_function == "plot":
                    plot(fp_data, t, out_filename = fname, out_format = 'png',
                         lon_range = lon_range, lat_range = lat_range,
                         log_range = log_range, plot_borders = plot_borders, 
                         zoom = zoom, colormap = colormap_fp,
                         dpi = dpi, figsize = figsize)
                elif plot_function == "plot_particle_location":
                    plot_particle_location(fp_data, t, out_filename = fname, out_format = 'png',
                                           log_range = log_range, 
                                           particle_direction = 'nw', colormap_fp = colormap_fp,
                                           colormap_part = colormap_part,
                                           particle_clevs = particle_clevs, dpi = dpi,
                                           figsize = figsize)
                     
            #pbar.update(ti)
            #print("")
        #pbar.finish()
    
    print("")
    print("... running ffmpeg")

    if video_os.lower() == "mac":
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            "-f mp4 -vcodec libx264 " + \
            "-pix_fmt yuv420p -intra -qscale 0 -y " + \
            os.path.join(output_directory, file_label) + ".mp4", shell=True)
    elif video_os.lower() == "pc":
#        os.remove(os.path.join(output_directory, file_label) + ".wmv")
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            " -b 5000k -f asf -vcodec wmv2 -acodec wmav2 " + \
            os.path.join(output_directory, file_label) + ".wmv", shell=True)
    else:
        print("ERROR: video_os must be mac or pc")
        return None
    
    print("... done with status " + str(ffmpeg_status))

    if delete_png:
        filelist = glob.glob(os.path.join(output_directory, "*.png"))
        for f in filelist:
            os.remove(f)


class get_country(object):
  def __init__(self, domain, country_file=None):
        
        if country_file is None:
            filename=glob.glob(join(data_path,'LPDM/countries/',f"country_{domain}.nc"))
            f = xr.open_dataset(filename[0])
        else:
            filename = country_file
            f = xr.open_dataset(filename)
    
        lon = f.variables['lon'][:].values
        lat = f.variables['lat'][:].values
    
        #Get country indices and names
        if "country" in f.variables:
            country = f.variables['country'][:, :]
        elif "region" in f.variables:
            country = f.variables['region'][:, :]
        
#         if (ukmo is True) or (uk_split is True):
#             name_temp = f.variables['name'][:]  
#             f.close()
#             name=np.asarray(name_temp)
        
#         else:
        name_temp = f.variables['name'].values
        f.close()

        name_temp = np.ma.filled(name_temp,fill_value=None)

        name=[]
        for ii in range(len(name_temp)):
            if type(name_temp[ii]) is not str:
                name.append(''.join(name_temp[ii].decode("utf-8")))
            else:
                name.append(''.join(name_temp[ii]))
        name=np.asarray(name)
    
    
        self.lon = lon
        self.lat = lat
        self.lonmax = np.max(lon)
        self.lonmin = np.min(lon)
        self.latmax = np.max(lat)
        self.latmin = np.min(lat)
        self.country = np.asarray(country)
        self.name = name 

def timeseries_HiTRes(flux_dict, fp_HiTRes_ds=None, fp_file=None, output_TS = True, output_fpXflux = True,
                      output_type='Dataset', output_file=None, verbose=False, chunks=None,
                      time_resolution='2H'):
    """
    The timeseries_HiTRes function computes flux * HiTRes footprints.
    
    HiTRes footprints record the footprint at each 2 hour period back in time for the first 24 hours.
    Need a high time resolution flux to multiply the first 24 hours back of footprints.
    Need a residual flux to multiply the residual integrated footprint for the remainder of the 30 
    day period.
    
    Args:
        fp_HiTRes_ds (xarray.Dataset)
            Dataset of High Time resolution footprint. HiTRes footprints record the footprint at 
            each 2 hour period back in time for the first 24 hours.
        domain (str)
            Domain name. The footprint files should be sub-categorised by the domain.
        flux_dict (dict)
            This should be a dictionary of the form output in the format
            flux_dict: {'high_freq': flux_dataset, 'low_freq': flux_dataset}.
            This is because this function needs two time resolutions of fluxes as
            explained in the header.
            
            If there are multiple sectors, the format should be:
            flux_dict: {sector1 : {'high_freq' : flux_dataset, 'low_freq'  : flux_dataset},
                        sector2 : {'high_freq' : flux_dataset, 'low_freq'  : flux_dataset}}
        output_TS (bool)
            Output the timeseries. Default is True.
        output_fpXflux (bool)
            Output the sensitivity map. Default is True.
        verbose (bool)
            show progress bar throughout loop
        chunks (dict)
            size of chunks for each dimension
            e.g. {'lat': 50, 'lon': 50}
            opens dataset with dask, such that it is opened 'lazily'
            and all of the data is not loaded into memory
            defaults to None - dataset is opened with out dask
    
    Returns:
        xarray.Dataset or dict
            Same format as flux_dict['high_freq']:
                If flux_dict['high_freq'] is an xarray.Dataset then an xarray.Dataset is returned
                If flux_dict['high_freq'] is a dict of xarray.Datasets then a dict of xarray.Datasets
                is returned (an xarray.Dataset for each sector)
        
            If output_TS is True:
                Outputs the timeseries
            If output_fpXflux is True:
                Outputs the sensitivity map   
    """
    if verbose:
        print(f'\nCalculating timeseries with {time_resolution} resolution, this might take a few minutes')
    ### get the high time res footprint
    if fp_HiTRes_ds is None and fp_file is None:
        print('Must provide either a footprint Dataset or footprint filename')
        return None
    elif fp_HiTRes_ds is None:
        fp_HiTRes_ds = read_netcdfs(fp_file, chunks=chunks)
        fp_HiTRes    = fp_HiTRes_ds.fp_HiTRes
    else:
        fp_HiTRes = fp_HiTRes_ds if type(fp_HiTRes_ds)==xr.core.dataarray.DataArray else \
                    fp_HiTRes_ds.fp_HiTRes.chunk(chunks) \
                    if fp_HiTRes_ds.chunks is None and chunks is not None else fp_HiTRes_ds.fp_HiTRes
    
    # resample fp and extract array to make the loop quicker
    fp_HiTRes = fp_HiTRes.resample(time=time_resolution).ffill()
    fp_HiTRes = da.array(fp_HiTRes)

    # convert fluxes to dictionaries with sectors as keys
    # flux = {freq: flux_freq if isinstance(flux_freq, dict) else {'total': flux_freq}
    #         for freq, flux_freq in flux_dict.items()}
    flux = {'total': flux_dict} if any([ff in list(flux_dict.keys()) for ff in ['high_freq', 'low_freq']]) else flux_dict
    flux = {freq: {sector: None if flux_sector is None else flux_sector
                   for sector, flux_sector in flux_freq.items()}
            for freq, flux_freq in flux.items()}
    # extract the required time data
    flux = {sector: {freq: flux_freq.sel(time=slice(fp_HiTRes_ds.time[0] - np.timedelta64(24, 'h'),
                                                    fp_HiTRes_ds.time[-1])).flux
                   for freq, flux_freq in flux_sector.items()}
            for sector, flux_sector in flux.items()}
    # convert to array to use in numba loop
    flux = {sector: {freq: flux_freq.values if flux_freq.chunks is None else da.array(flux_freq)
                   for freq, flux_freq in flux_sector.items()}
            for sector, flux_sector in flux.items()}
    
    # create time array to loop through, with the required resolution
    # fp_HiTRes.time is the release time of particles into the model
    time_array = np.arange(fp_HiTRes_ds.time[0].values, fp_HiTRes_ds.time[-1].values,
                           int(time_resolution[0]),
                           dtype=f'datetime64[{time_resolution[1].lower()}]')
    
    # Set up a numpy array to calculate the product of the footprint (H matrix) with the fluxes
    if output_fpXflux:
        fpXflux   = {sector: da.zeros((len(fp_HiTRes_ds.lat),
                                       len(fp_HiTRes_ds.lon),
                                       len(time_array)))
                     for sector in flux.keys()}
    
    elif output_TS:
        timeseries = {sector: da.zeros(len(time_array)) for sector in flux.keys()}
    
    # month and year of the start of the data - used to index the low res data
    start = {dd: getattr(time_array[0].astype(object), dd) for dd in ['month', 'year']}
    
    # get the data timesteps - used for resampling the hourly footprints
    num = int(time_resolution[0])
    
    # put the time array into tqdm if we want a progress bar to show throughout the loop
    iters = tqdm(time_array) if verbose else time_array
    ### iterate through the time coord to get the total mf at each time step using the H back coord
    # at each release time we disaggregate the particles backwards over the previous 24hrs
    for tt, time in enumerate(iters):
        # get 4 dimensional chunk of high time res footprint for this timestep
        # units : mol/mol/mol/m2/s
        # reverse the H_back coordinate to be chronological, and resample
        fp_time   = fp_HiTRes[:,:,tt,::-num]
                
        # get the correct index for the low res data
        # estimated using the difference between the current and start month and year
        current = {dd: getattr(time.astype(object), dd) for dd in ['month', 'year']}
        tt_low = current['month'] - start['month'] + 12*(current['year']-start['year'])

        # select the high res emissions for the corresponding 24 hours
        # if there aren't any high frequency data it will select from the low frequency data
        # this is so that we can compare emissions data with different resolutions e.g. ocean species
        emissions     = {sector: flux_sector['high_freq'][:,:,tt+1:tt+fp_time.shape[2]]
                                 if flux_sector['high_freq'] is not None else
                                 flux_sector['low_freq'][:,:,tt_low]
                         for sector, flux_sector in flux.items()}
        # add an axis if the emissions is array is 2D so that it can be multiplied by the fp
        emissions     = {sector: em_sec[:,:,np.newaxis] if len(em_sec.shape)==2 else em_sec
                        for sector, em_sec in emissions.items()}
        # select average monthly emissions for the start of the month
        emissions_end = {sector: flux_sector['low_freq'][:,:,tt_low]
                         for sector, flux_sector in flux.items()}
        
        # Multiply the HiTRes footprint with the HiTRes emissions to give mf
        # we take all but the slice for H_back==24 as these are the hourly disaggregated fps
        # flux units : mol/m2/s;       fp units : mol/mol/mol/m2/s
        # --> mol/mol/mol/m2/s * mol/m2/s === mol / mol
        fpXflux_time  = {sector: em_sec * fp_time[:,:,1:] for sector, em_sec in emissions.items()}
        # multiply the monthly flux by the residual fp, at H_back==24
        fpXflux_end   = {sector: em_end * fp_time[:,:,0] for sector, em_end in emissions_end.items()}
        # append the residual emissions
        fpXflux_time  = {sector: np.dstack((fp_fl, fpXflux_end[sector]))
                         for sector, fp_fl in fpXflux_time.items()}
        
        for sector, fp_fl in fpXflux_time.items():
            if output_fpXflux:
                # Sum over time (H back) to give the total mf at this timestep
                fpXflux[sector][:,:,tt] = fp_fl.sum(axis=2)
            
            elif output_TS:
                # work out timeseries by summing over lat, lon, & time (24 hrs)
                timeseries[sector][tt] = fp_fl.sum()
    
    if output_fpXflux and not output_TS:
        fpXflux = xr.Dataset(fpXflux,
                             coords={'lat'  : fp_HiTRes.lat.values,
                                     'lon'  : fp_HiTRes.lon.values,
                                     'time' : time_array}) \
                  if output_type.lower()=='dataset' else fpXflux
        return fpXflux.compute()
    
    else:
        if output_fpXflux:
            # if not already done then calculate the timeseries
            timeseries = {sector: fp_fl.sum(axis=(0,1)) for sector, fp_fl in fpXflux.items()}
            fpXflux = xr.Dataset(fpXflux,
                                 coords={'lat'  : fp_HiTRes.lat.values,
                                         'lon'  : fp_HiTRes.lon.values,
                                         'time' : time_array}) \
                      if output_type.lower()=='dataset' else \
                      {sector: xr.DataArray(data = fpXflux_sector,
                                   dims = ['lat', 'lon', 'time'],
                                   coords = {'lat' : fp_HiTRes.lat.values,
                                             'lon'  : fp_HiTRes.lon.values,
                                             'time' : time_array})
                       for sector, fpXflux_sector in fpXflux.items()} \
                      if output_type.lower()=='dataarray' else fpXflux
        # for each sector create a tuple of ['time'] (coord name) and the timeseries
        # if the output required is a dataset
        timeseries = {sec : (['time'], ts_sec) for sec, ts_sec in timeseries.items()} \
                     if output_type.lower()=='dataset' else timeseries
        timeseries = xr.Dataset(timeseries,
                                coords={'time' : time_array}) \
                     if output_type.lower()=='dataset' else \
                     {sector: xr.DataArray(data = ts_sector,
                                           dims = ['time'],
                                           coords = {'time' : time_array})
                      for sector, ts_sector in timeseries.items()} \
                     if output_type.lower()=='dataarray' else timeseries
        if output_type.lower()=='dataarray' and list(flux.keys())==['total']:
            if output_fpXflux: fpXflux = fpXflux['total']
            timeseries = timeseries['total']
        
        if output_file is not None and output_type.lower()=='dataset':
            print(f'Saving to {output_file}')
            timeseries.compute().to_netcdf(output_file)
        elif output_file is not None:
            print(f'output type must be dataset to save to file')
        
        if output_fpXflux:
            return timeseries.compute(), fpXflux.compute()
        elif output_TS:
            return timeseries.compute()