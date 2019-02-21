# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 10:45:51 2014

"""
from __future__ import print_function

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import datetime as dt
import os
import glob
from matplotlib import ticker
import pandas as pd
import bisect
import subprocess
#from progressbar import ProgressBar
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

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

if acrg_path is None:
    acrg_path = os.getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")
if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")


# Get acrg_site_info file
with open(join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f)

def open_ds(path):
    
    """
    Function efficiently opens xray datasets.
    """
    # use a context manager, to ensure the file gets closed after use
    with xr.open_dataset(path) as ds:
        ds.load()
    return ds 

def filenames(site, domain, start, end, height, fp_directory):
    """
    The filenames function outputs a list of available footprint file names,
    for given site, domain, directory and date range.
    
    Expect filenames of the form:
        [fp_directory]/domain/site*-height*domain*yearmonth*.nc
        e.g. [/data/shared/NAME/fp]/EUROPE/MHD-10magl_EUROPE_201401.nc
    
    Args:
        site (str)         : Site name. Full list of site names should be 
                             defined within acrg_site_info.json
        domain (str)       : Domain name. The footprint files should be 
                             sub-categorised by the NAME domain name.
        start (str)        : Start date in format "YYYY-MM-DD" for range of files
                             to find.
        end (str)          : End date in same format as start for range of files
                             to find.
        height (str)       : Height related to input data. 
        fp_directory (str) : fp_directory can be specified if files are not in 
                             the default directory must point to a directory 
                             which contains subfolders organized by domain.
    Returns:
        list (str): matched filenames
    """

    baseDirectory = fp_directory
        
    # Read site info for heights
    if height is None:
        if not site in site_info.keys():
            print("Site code not found in arcg_site_info.json to get height information. " + \
                  "Check that site code is as intended. "+ \
                  "If so, either add new site to file or input height manually.")
            return None
        height = site_info[site]["height_name"][0]
    
    # Convert into time format
    months = pd.DatetimeIndex(start = start, end = end, freq = "M").to_pydatetime()
    yearmonth = [str(d.year) + str(d.month).zfill(2) for d in months]

    files = []
    for ym in yearmonth:
        f=glob.glob(baseDirectory + \
            domain + "/" + \
            site + "*" + "-" + height + "*" + domain + "*" + ym + "*.nc")

        if len(f) > 0:
            files += f

    files.sort()

    if len(files) == 0:
        print("Can't find file: " + baseDirectory + \
            domain + "/" + \
            site + "*" + height + "*" + domain + "*" + "*.nc")
    return files

def read_netcdfs(files, dim = "time"):
    """
    The read_netcdfs function uses xarray to open sequential netCDF files and 
    and concatenates them along the specified dimension.
    Note: this function makes sure that file is closed after open_dataset call.
    
    Args:
        files (list) : List of netCDF filenames.
        dim (str)    : Dimension of netCDF to use for concatenating the files.
                       Default= "time".
    
    Returns:
        xarray.Dataset : all files open as one concatenated xarray.Dataset object    
    """
    
    #def process_one_path(path):
    #    with xr.open_dataset(path) as ds:
    #        ds.load()
    #    return ds
    
    print("Reading and concatenating files: ")
    for fname in files:
        print(fname)
    
    datasets = [open_ds(p) for p in sorted(files)]
    combined = xr.concat(datasets, dim)
    return combined   

def interp_time(bc_ds,vmr_var_names, new_times):
    """
    The interp_time function interpolates the times of the VMR variable 
    'vmr_var_name' in the xarray.Dataset 'bc_ds' to the times specified in 
    'interp_times'. The variable must have dimensions (height, lat_or_lon, time) 
    in that order. 
    Note: This function was created to convert MOZART monthly averages into 
    same frequency as NAME footprints.

    TODO: Add details for vmr_var_names and new_times

    Args:
        bc_ds (xarray.Dataset)   : Output from boundary_conditions() function
        vmr_var_names (iterable) : ???
        new_times                : ???
    
    Returns:
        xarray.Dataset : new dataset with the VMRs recalculated at interpolated times.

    """

    vmr_dict={}

    for vi,vmr_var_name in enumerate(vmr_var_names):

        x_id= np.arange(len(bc_ds.time))
        new_times_id = np.linspace(0.,np.max(x_id), num=len(new_times)) 
        vmr_new = np.zeros((len(bc_ds.height),len(bc_ds[vmr_var_name][0,:,0]),len(new_times)))
        for j in range(len(bc_ds.height)):
            for i in range(len(bc_ds[vmr_var_name][0,:,0])):
                y = bc_ds[vmr_var_name][j,i,:]
                f = interp1d(x_id,y, bounds_error = False,kind='linear', 
                                         fill_value = np.max(y))
                vmr_new[j,i,:] = f(new_times_id)

        vmr_dict[vmr_var_name]=vmr_new
        
    ds2 = xr.Dataset({"vmr_n": (["height", "lon", "time"],vmr_dict["vmr_n"]),
                        "vmr_e": (["height", "lat", "time"],vmr_dict["vmr_e"]),
                        "vmr_s": (["height", "lon", "time"],vmr_dict["vmr_s"]),
                        "vmr_w": (["height", "lat", "time"],vmr_dict["vmr_w"])},
                        coords={"lon":bc_ds.lon, "lat": bc_ds.lat, "time": new_times,
                                "height":bc_ds.height})

    return ds2


def footprints(sitecode_or_filename, fp_directory = None, 
               flux_directory = None, bc_directory = None,
               start = None, end = None, domain = None, height = None,
               species = None, emissions_name = None, HiTRes = False,interp_vmr_freq=None):

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
    Note: fp_directory, flux_directory and bc_directory can point to specified directories
    but if not specified, default directories will be used (set at the top of file).

    Args:
        sitecode_or_filename : Site (e.g. 'MHD') or a netCDF filename (*.nc) (str)
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
        start (str)          : Start date in format "YYYY-MM-DD" for range of files to find.
        end (str)            : End date in same format as start for range of files to find.
        domain (str)         : Domain name. The footprint files should be sub-categorised by the domain.
        height (str)         : Height related to NAME. If the HEIGHT keyword is not specified, the default 
                               height from the acrg_site_info.json file is assumed.
        species (str)        : Species name. All species names are defined acrg_species_info.json.
        emissions_name (str) : Allows emissions files such as co2nee_EUROPE_2012.nc to be read in. 
                               In this case EMISSIONS_NAME would be 'co2nee'
        HiTRes (bool)        : Whether to include high time resolution footprints.
        interp_vmr_freq      : Frequency to interpolate vmr time. (float/int)
        
    Returns:
        xarray.Dataset : combined footprint files
    """
    # Chose whether we've input a site code or a file name
    # If it's a three-letter site code, assume it's been processed
    # into an annual footprint file in (mol/mol) / (mol/m2/s)
    # using acrg_name_process
    
    if fp_directory is None:
        fp_integrated_directory = join(data_path, 'NAME/fp/')
        fp_HiTRes_directory = join(data_path,'NAME/fp_high_time_res/')
        fp_directory = {'integrated': fp_integrated_directory,
                        'HiTRes': fp_HiTRes_directory}

    #   Error message if looking for HiTRes files and fp_directory is not a dictionary        
        if HiTRes is True:
            if type(fp_directory) is not dict:
                print("fp_directory needs to be a dictionary containing paths \
                       to integrated and HiTRes footprints \
                       {integrated:path1, HiTRes:path2}")
                return None
                
                print("As HiTRes is set to True, make sure that the high and low time resolution emissions name\
                      pairs are input correctly for the emissions sources where HiTRes applies. They should look like:\
                      emissions_name = {source_name:{'high_res':emissions_file_identifier,\
                      'low_res':emissions_file_identifier}.")
    
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
            files = filenames(site, domain, start, end, height, fp_directory["integrated"])
        else:
            files = filenames(site, domain, start, end, height, fp_directory)

    if len(files) == 0:
        print("Can't find files, " + sitecode_or_filename)
        return None

    else:
        fp=read_netcdfs(files)  
        
        # If a species is specified, also get flux and vmr at domain edges
        if emissions_name is not None:
            flux_ds = flux(domain, emissions_name, flux_directory=flux_directory)
            if flux_ds is not None:
                fp = combine_datasets(fp, flux_ds, method='ffill')
        elif species is not None:
            flux_ds = flux(domain, species, flux_directory=flux_directory)
            if flux_ds is not None:
                fp = combine_datasets(fp, flux_ds, method='ffill')
        
        if species is not None:
            bc_ds = boundary_conditions(domain, species, bc_directory=bc_directory)
            if bc_ds is not None:                   
                if interp_vmr_freq is not None:
                    # Interpolate bc_ds between months to same timescale as footprints
                    dum_ds = bc_ds.resample(interp_vmr_freq, "time")
                    new_times=dum_ds.time            
                    vmr_var_names=["vmr_n", "vmr_e", "vmr_s", "vmr_w"]
                    bc_ds = interp_time(bc_ds,vmr_var_names, new_times)  
                fp = combine_datasets(fp, bc_ds, method='ffill')

        if HiTRes == True:
            HiTRes_files = filenames(site, domain, start, end, height, fp_directory["HiTRes"])
            HiTRes_ds = read_netcdfs(HiTRes_files)
            fp = combine_datasets(fp, HiTRes_ds, method='ffill')

        return fp


def flux(domain, species, start = None, end = None, flux_directory=None):
    """
    The flux function reads in all flux files for the domain and species as an xarray Dataset.
    Note that at present ALL flux data is read in per species per domain or by emissions name.
    To be consistent with the footprints, fluxes should be in mol/m2/s.
    
    Expect filenames of the form:
        [flux_directory]/domain/species.lower()_*.nc
        e.g. [/data/shared/NAME/emissions]/EUROPE/ch4_EUROPE_2013.nc
    
    TODO: This may get slow for very large flux datasets, and we may want to subset.
    
    Args:
        domain (str)         : Domain name. The flux files should be sub-categorised by the domain.
        species (str)        : Species name. All species names are defined acrg_species_info.json.
        start (str)          : Start date in format "YYYY-MM-DD" to output only a time slice of all the flux files.
                               The start date used will be the first of the input month. I.e. if "2014-01-06" is input,
                               "2014-01-01" will be used.  This is to mirror the time slice functionality of the filenames function.
        end (str)            : End date in same format as start to output only a time slice of all the flux files.
                               The end date used will be the first of the input month and the timeslice will go up
                               to, but not include, this time. I.e. if "2014-02-25' is input, "2014-02-01" will be used.
                               This is to mirror the time slice functionality of the filenames function.
        flux_directory (str) : flux_directory can be specified if files are not in 
                               the default directory. Must point to a directory 
                               which contains subfolders organized by domain. (optional)
    Returns:
        xarray.Dataset : combined dataset of all matching flux files
    """
    
    if flux_directory is None:
        flux_directory = join(data_path, 'NAME/emissions/')

    print(("filename",flux_directory + domain + "/" + species.lower() + "_" + "*.nc"))
    files = sorted(glob.glob(flux_directory + domain + "/" + 
                   species.lower() + "_" + "*.nc"))
    if len(files) == 0:
        print("Can't find flux: " + domain + " " + species)
        return None
    
    flux_ds = read_netcdfs(files)
    # Check that time coordinate is present
    if not "time" in flux_ds.coords.keys():
        print("ERROR: No 'time' coordinate " + \
              "in flux dataset for " + domain + ", " + species)
        return None

    # Check for level coordinate. If one level, assume surface and drop
    if "lev" in flux_ds.coords.keys():
        print("WARNING: Can't support multi-level fluxes. Trying to remove 'lev' coordinate " + \
              "from flux dataset for " + domain + ", " + species)
        if len(flux_ds.lev) > 1:
            print("ERROR: More than one flux level")
        else:
            return flux_ds.drop("lev")
        
    if start == None:
        return flux_ds
    else:
        if end == None:
            print("To get fluxes for a certain time period you must specify an end date.")
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
                print("Slicing time to range {} - {}".format(month_start,month_end))
            
            return flux_timeslice


def flux_for_HiTRes(domain, emissions_dict, start=None, end=None, flux_directory=None):
    """
    Creates a dictionary of high and low frequency fluxes for use with HiTRes footprints.
    
    Args:
        domain (str)          : Domain name. The flux files should be sub-categorised by the 
                                domain.
        emissions_dict (dict) : This should be a dictionary of the form:
                                {'high_freq':high_freq_emissions_name, 'low_freq':low_freq_emissions_name}
                                e.g. {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}.
        start (str)           : Start date in format "YYYY-MM-DD" to output only a time slice of all the flux files.
                                The start date used will be the first of the input month. I.e. if "2014-01-06" is input,
                                "2014-01-01" will be used.  This is to mirror the time slice functionality of the filenames function.
        end (str)             : End date in same format as start to output only a time slice of all the flux files.
                                The end date used will be the first of the input month and the timeslice will go up
                                to, but not include, this time. I.e. if "2014-02-25' is input, "2014-02-01" will be used.
                                This is to mirror the time slice functionality of the filenames function.
        flux_directory (str)  : flux_directory can be specified if files are not in 
                                the default directory. Must point to a directory 
                                which contains subfolders organized by domain. (optional)
    Returns:
        dictionary {'high_freq': flux_xray_dataset, 'low_freq': flux_xray_dataset} :
            dictionary containing xray datasets for both high and low frequency fluxes.
    """
    
    if 'low_freq' not in emissions_dict.keys():
        print("low_freq must be a key in the emissions_dict in order to combine with HiTRes footprints.")
        return None
    elif 'high_freq' not in emissions_dict.keys():
        print("high_freq must be a key in the emissions_dict in order to use HiTRes footprints.")
        return None
    
    flux_dict = {}
    
    if start:
        #Get the month before the one one requested because this will be needed for the first few
        #days in timeseries_HiTRes to calculate the modleed molefractions for times when the footprints
        #are in the previous month.
        start = str(pd.to_datetime(start) - dateutil.relativedelta.relativedelta(months=1))
    
    fluxes_highfreq = flux(domain, emissions_dict['high_freq'], start = start, end = end,flux_directory=flux_directory)
    fluxes_lowfreq = flux(domain, emissions_dict['low_freq'], start = start, end = end,flux_directory=flux_directory)
    
    flux_dict['low_freq'] = fluxes_lowfreq
    flux_dict['high_freq'] = fluxes_highfreq
    
    return flux_dict


def boundary_conditions(domain, species, start = None, end = None, bc_directory=None):
    """
    The boundary_conditions function reads in the files with the global model vmrs at the domain edges 
    to give the boundary conditions as an xarray Dataset.

    Expect filenames of the form:
        [bc_directory]/domain/species.lower()_*.nc
        e.g. [/data/shared/NAME/bc]/EUROPE/ch4_EUROPE_201301.nc

    Args:
        domain (str)       : Domain name. The boundary condition files should be sub-categorised by the 
                             domain.
        species (str)      : Species name. All species names are defined acrg_species_info.json.
        bc_directory (str) : bc_directory can be specified if files are not in 
                             the default directory. Must point to a directory 
                             which contains subfolders organized by domain. (optional)
    Returns:
        xarray.Dataset : combined dataset of matching boundary conditions files
    """
    
    if bc_directory is None:
        bc_directory = join(data_path, 'NAME/bc/')
    
    files = sorted(glob.glob(bc_directory + domain + "/" + 
                   species.lower() + "_" + "*.nc"))
    if len(files) == 0:
        print("Can't find boundary condition files: " + domain + " " + species)
        return None

    bc_ds = read_netcdfs(files)

    if start == None:
        return bc_ds
    else:
        if end == None:
            print("To get boundary conditions for a certain time period you must specify an end date.")
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
                print("No boundary conditions available during the time period specified so outputting\
                          boundary conditions from %s" %bc_timeslice.time.values[0])
            return bc_timeslice


def basis(domain, basis_case, basis_directory = None):
    """
    The basis function reads in the all matching files for the basis case and domain as an xarray Dataset.
    
    Expect filenames of the form:
        [basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/NAME/basis_functions]/EUROPE/sub_transd_EUROPE_2014.nc

    TODO: More info on options for basis functions.

    Args:
        domain (str)       : Domain name. The basis files should be sub-categorised by the domain.
        basis_case (str)   : Basis case to read in. Examples of basis cases are "voroni","sub-transd",
                             "sub-country_mask","INTEM".
        basis_directory (str) : basis_directory can be specified if files are not in 
                             the default directory. Must point to a directory 
                             which contains subfolders organized by domain. (optional)
    Returns:
        xarray.Dataset : combined dataset of matching basis functions
    """
    if basis_directory is None:
        basis_directory = join(data_path, 'NAME/basis_functions/')
        
    files = sorted(glob.glob(basis_directory + domain + "/" +
                    basis_case + "_" + domain + "*.nc"))
    if len(files) == 0:
        print("Can't find basis functions: " + domain + " " + basis_case)
        return None

    basis_ds = read_netcdfs(files)

    return basis_ds

def basis_boundary_conditions(domain, basis_case, bc_basis_directory = None):
    """
    The basis_boundary_conditions function reads in all matching files for the boundary conditions 
    basis case and domain as an xarray Dataset.
    
    Expect filesnames of the form:
        [bc_basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/NAME/bc_basis_directory]/EUROPE/NESW_EUROPE_2013.nc

    TODO: More info on options for basis functions.
    
    Args:
        domain (str)             : Domain name. The basis files should be sub-categorised by the domain.
        basis_case (str)         : Basis case to read in. Examples of basis cases are "NESW","stratgrad".
        bc_basis_directory (str) : bc_basis_directory can be specified if files are not in 
                                   the default directory. Must point to a directory 
                                   which contains subfolders organized by domain. (optional)
    Returns:
        xarray.Datset : combined dataset of matching basis functions
    """
    
    if bc_basis_directory is None:
        bc_basis_directory = join(data_path,'NAME/bc_basis_functions/')
    
    files = sorted(glob.glob(bc_basis_directory + domain + "/" +
                    basis_case + '_' + domain + "*.nc"))

    if len(files) == 0:
        print("Can't find boundary condition basis functions: " + domain + " " + basis_case)
        return None
    
    basis_ds = read_netcdfs(files)

    return basis_ds


def combine_datasets(dsa, dsb, method = "ffill", tolerance = None):
    """
    The combine_datasets function merges two datasets and re-indexes to the FIRST dataset.
    If "fp" variable is found within the combined dataset, the "time" values where the "lat","lon"
    dimensions didn't match are removed.
    
    Example:
        ds = combine_datasets(dsa, dsb)

    Args:
        dsa (xarray.Dataset) : First dataset to merge
        dsb (xarray.Dataset) : Second dataset to merge
        method (str)         : One of {None, ‘nearest’, ‘pad’/’ffill’, ‘backfill’/’bfill’}
                               See xarray.DataArray.reindex_like for list of options and meaning.
        tolerance (??)      : Maximum allowed tolerance between matches.

    Returns:
        xarray.Dataset: combined dataset indexed to dsa
    """
    # merge the two datasets within a tolerance and remove times that are NaN (i.e. when FPs don't exist)
    
    ds_temp = dsa.merge(dsb.reindex_like(dsa, method, tolerance = tolerance))
    if 'fp' in ds_temp.keys():
        flag = np.where(np.isfinite(ds_temp.fp.mean(dim=["lat","lon"]).values))
        ds_temp = ds_temp[dict(time = flag[0])]
    return ds_temp

def timeseries(ds):
    """
    The timeseries function compute flux * footprint time series.

    Example:
        ts = timeseries(dataset)

    TODO: There are almost certainly much more efficient ways of doing this.

    Args:
        ds (xarray.Dataset) : Dataset with both the flux and footprint fields present
    Returns:
        xarray.Dataset        
    """

    if "flux" in ds.keys():
        return (ds.fp*ds.flux).sum(["lat", "lon"])
    else:
        print("Can't calculate time series " + \
              "no fluxes. Check flux file.")
        return None

def timeseries_HiTRes(fp_HiTRes_ds, flux_dict, output_TS = True, output_fpXflux = True):
    """
    The timeseries_HiTRes function computes flux * HiTRes footprints.
    
    HiTRes footprints record the footprint at each 2 hour period back in time for the first 24 hours.
    Need a high time resolution flux to multiply the first 24 hours back of footprints.
    Need a residual flux to multiply the residual integrated footprint for the remainder of the 30 
    day period.
    
    Args:
        fp_HiTRes_ds (xarray.Dataset) : Dataset of High Time resolution footprint. HiTRes footprints 
                                        record the footprint at each 2 hour period back in time for the 
                                        first 24 hours.
        domain (str)                  : Domain name. The footprint files should be sub-categorised by the 
                                        domain.
        flux_dict (dict)              : This should be a dictionary of the form output in the function
                                        flux_for_HiTRes: {'high_freq': flux_dataset, 'low_freq': flux_dataset}.
                                        This is because this function needs two time resolutions of fluxes as
                                        explained in the header.
        output_TS (bool)              : Output the timeseries. Default is True.
        output_fpXflux (bool)         : Output the sensitivity map. Default is True.

    
    Returns:
        xarray.Dataset / tuple(xarray.Dataset,xarray.Dataset)
        
        If output_TS is True:
            Outputs the timeseries
        If output_fpXflux is True:
            Outputs the sensitivity map   
    """
    
#    # Get time range
#    fp_start = fp_HiTRes_ds.time[0]
#    start = dt.datetime(int(fp_start['time.year']), int(fp_start['time.month']), 1, 0, 0)
#        
#    fp_end = fp_HiTRes_ds.time[-1]
#    month_days = calendar.monthrange(int(fp_end['time.year']), int(fp_end['time.month']))[1]
#    end = dt.datetime(int(fp_end['time.year']), int(fp_end['time.month']), 1, 0, 0) + \
#            dt.timedelta(days = month_days)
    
    flux_HiTRes = flux_dict['high_freq']
    flux_resid = flux_dict['low_freq']
    
    """
    Probably need a check in here to make sure dates of fluxes correspond to dates of footprints
    """
    
    fp_HiTRes = fp_HiTRes_ds.fp_HiTRes.to_dataset()
    fpXflux = np.zeros((len(fp_HiTRes.lat), len(fp_HiTRes.lon), len(fp_HiTRes.time)))
    
    for ti, time in enumerate(fp_HiTRes.time):
        fp = fp_HiTRes.sel(time=time).fp_HiTRes.to_dataset()
        time_back = [fp.time.values - np.timedelta64(i,'h') for i in fp.H_back.values]
        fp = fp.update({'H_back':time_back})
        fp = fp.drop('time')
        fp = fp.rename({'H_back':'time'})

        new_fp = fp.fp_HiTRes[:,:,::-1]
        new_time = fp.time[::-1]
        new_ds = xr.Dataset({'fp_HiTRes':(['lat','lon','time'], new_fp)},
                               coords={'lat':fp.lat,
                                       'lon':fp.lon,
                                       'time':new_time})

        em = flux_HiTRes.reindex_like(new_ds, method='ffill')
        
        #Use end of hours back as closest point for finding the emissions file
        print("flux_resid",flux_resid)
        print("new_ds.time[0]",new_ds.time[0])
        import pdb
        pdb.set_trace()
        emend = flux_resid.sel(time = new_ds.time[0], method = 'ffill')
        em.flux[:,:,0] = emend.flux
        fpXflux[:,:,ti] = (new_ds.fp_HiTRes*em.flux).sum(["time"])
        
    timeseries= np.sum(fpXflux, axis = (0,1))
    
    if output_fpXflux == True and output_TS ==True:
        return timeseries, fpXflux
    
    elif output_fpXflux == False and output_TS ==True:
        return timeseries
        
    elif output_fpXflux == True and output_TS ==False:
        return fpXflux       

def timeseries_boundary_conditions(ds):
    """
    The timeseries_boundary_conditions function compute particle location * global model edges time series.
    
    Args:
        ds (xarray.Dataset) : Dataset with both the particle locations and vmr at domain edge fields 
                              present.
    
    Returns:
        xarray.Dataset
    """ 

    return (ds.particle_locations_n*ds.vmr_n).sum(["height", "lon"]) + \
           (ds.particle_locations_e*ds.vmr_e).sum(["height", "lat"]) + \
           (ds.particle_locations_s*ds.vmr_s).sum(["height", "lon"]) + \
           (ds.particle_locations_w*ds.vmr_w).sum(["height", "lat"])

    
def footprints_data_merge(data, domain, load_flux = True, load_bc = True,
                          calc_timeseries = True, calc_bc = True, HiTRes = False,
                          average = None, site_modifier = {}, height = None,
                          emissions_name = None, interp_vmr_freq = None,
                          fp_directory = None,
                          flux_directory = None,
                          bc_directory = None,
                          resample_to_data = False):
#                          perturbed=False, fp_dir_pert=None, pert_year=None, pert_month=None):

    """
    Output a dictionary of xarray footprint datasets, that correspond to a given
    dictionary of Pandas dataframes.
    
    Args:
        data (dict)          : Input dictionary of dataframes with the sites as the keys.
                               Should match output from acrg_agage.get_obs() function. For example:
                               data = {"MHD": MHD_dataframe, "TAC": TAC_dataframe}
        domain (str)         : Domain name. The footprint files should be sub-categorised by the domain.
        load_flux (bool)     : True includes fluxes in output, False does not. Default True.
        load_bc (bool)       : True includes boundary conditions in output, False does not. Default True.
        calc_timeseries (bool) : True calculates modelled mole fractions for each site using fluxes, False does not. Default True.
        calc_bc (bool)       : True calculates modelled baseline for each site using boundary conditions, False does not. Default True.
        HiTRes (bool)        : Set to True to include HiTRes footprints in output. Default False.
        average (dict)       : Averaging period for each dataset (for each site). Should be a dictionary with
                               {site: averaging_period} key:value pairs.
                               Each value should be a string of the form e.g. "2H", "30min" (should match
                               pandas offset aliases format).
        site_modifier        : An optional site modifier dictionary is used that maps the site name in the
                               obs file to the site name in the footprint file, if they are different. This
                               is useful for example if the same site FPs are run with a different met and 
                               they are named slightly differently from the obs file. E.g.
                               site_modifier = {"DJI":"DJI-SAM"} - station called DJI, FPs called DJI-SAM
        height (dict)        : Height related to input data. Should be a dictionary with
                               {site: height} key:value pairs. Can be found from acrg_sites_info.json
        emissions_name (dict): Allows emissions files with filenames that are longer than just the species name
                               to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
                               with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
                               multiple sources can be read in simultaneously if they are added as separate entries to
                               the emissions_name dictionary.
                               If using HiTRes footprints, both the high and low frequency emissions files must be specified
                               in a second dictionary like so: {'anth': {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}}.
                               It is not a problem to have a mixture of sources, with some that use HiTRes footprints and some
                               that don't.
        interp_vmr_freq      : Frequency to interpolate vmr time. (float/int)
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
                                    
    
    Returns:
        Dictionary of the form {"MHD": MHD_xarray_dataset, "TAC": TAC_xarray_dataset, ".flux": dictionary_of_flux_datasets, ".bc": boundary_conditions_dataset}:
            combined dataset for each site
    """
    
    sites = [key for key in data.keys() if key[0] != '.']
    attributes = [key for key in data.keys() if key[0] == '.']
    
    if average is not None:
        if type(average) is not dict:
            print("WARNING: average list must be a dictionary with {site: averaging_period}\
                  key value pairs. Ignoring. Output dataset will not be resampled.")
            average = {x:None for x in sites}
    else:
        average = {x:None for x in sites}

    # If not given, check if species is defined in data dictionary:
#    if species is None:
    if ".species" in data.keys():
        species = data[".species"]
    else:
        print("Species can't be found in data dictionary.")

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

        # Dataframe for this site            
        site_df = data[site] 
        # Get time range
        df_start = min(site_df.index).to_pydatetime()
        start = dt.datetime(df_start.year, df_start.month, 1, 0, 0)
        
        df_end = max(site_df.index).to_pydatetime()
        month_days = calendar.monthrange(df_end.year, df_end.month)[1]
        end = dt.datetime(df_end.year, df_end.month, 1, 0, 0) + \
                dt.timedelta(days = month_days)
      
        #Get earliest start and latest end dates from all sites for loading fluxes and bcs
        if flux_bc_start == None or flux_bc_start > start:
            flux_bc_start = start
        if flux_bc_end == None or flux_bc_end < end:
            flux_bc_end = end
            
        # Convert to dataset
        site_ds = xr.Dataset.from_dataframe(site_df)
        
        if site in site_modifier.keys():
            site_modifier_fp = site_modifier[site]
        else:    
            site_modifier_fp = site
         
        if height is not None:
            
            if type(height) is not dict:
                print("Height input needs to be a dictionary with {sitename:height}")
                return None
                
            height_site = height[site] 
        else:
            height_site = height
        
        # Get footprints

        site_fp = footprints(site_modifier_fp, fp_directory = fp_directory, 
                             flux_directory = flux_directory, 
                             bc_directory = bc_directory,
                             start = start, end = end,
                             domain = domain,
                             species = None,
                             height = height_site,
                             emissions_name = None,
                             HiTRes = HiTRes,
                             interp_vmr_freq=interp_vmr_freq)                         
                        
        if site_fp is not None:                        
            # If satellite data, check that the max_level in the obs and the max level in the processed FPs are the same
            # Set tolerance tin time to merge footprints and data   
            # This needs to be made more general to 'satellite', 'aircraft' or 'ship'                
            if "GOSAT" in site.upper():
                ml_obs = site_df.max_level
                ml_fp = site_fp.max_level
                tolerance = 60e9 # footprints must match data with this tolerance in [ns]
                if ml_obs != ml_fp:
                    print("ERROR: MAX LEVEL OF SAT OBS DOES NOT EQUAL MAX LEVEL IN FP")
                    print("max_level_fp =",ml_fp)
                    print("max_level_obs =",ml_obs)
                    return None
            elif "GAUGE-FERRY" in site.upper():
                tolerance = '5min'
            elif "GAUGE-FAAM" in site.upper():
                tolerance = '1min'    
            else:
                tolerance = None
            
            # rt17603: 06/04/2018 - Added sort as some footprints weren't sorted by time for satellite data.
            site_fp = site_fp.sortby("time")
            
            #lw13938: 12/04/2018 - This should slice the date to the smallest time frame
            # spanned by both the footprint and obs, then resamples the data 
            #using the mean to the one with coarsest median resolution 
            #starting from the sliced start date. 
            ds_timeperiod = np.nanmedian((site_ds.time.data[1:] - site_ds.time.data[0:-1]).astype('int64')) 
            fp_timeperiod = np.nanmedian((site_fp.time.data[1:] - site_fp.time.data[0:-1]).astype('int64')) 
            ds_st = site_ds.time[0]
            ds_et = site_ds.time[-1]
            fp_st = site_fp.time[0]
            fp_et = site_fp.time[-1]
            if int(ds_st.data) > int(fp_st.data):
                start_date = ds_st
            else:  
                start_date = fp_st
            if int(ds_et.data) < int(fp_et.data):
                end_date = ds_et
            else:
                end_date = fp_et
            
            # rt17603: 24/07/2018 - Rounding to the nearest second(+/-1). Needed for sub-second dates otherwise sel was giving a KeyError
            start_s = str(np.round(start_date.data.astype(np.int64)-5e8,-9).astype('datetime64[ns]')) # subtract half a second to ensure lower range covered
            end_s = str(np.round(end_date.data.astype(np.int64)+5e8,-9).astype('datetime64[ns]')) # add half a second to ensure upper range covered
            
            site_ds = site_ds.sel(time=slice(start_s,end_s))
            site_fp = site_fp.sel(time=slice(start_s,end_s))
            
            #site_ds = site_ds.sel(time=slice(str(start_date.data),str(end_date.data)))
            #site_fp = site_fp.sel(time=slice(str(start_date.data),str(end_date.data)))
            
            base = start_date.dt.hour.data + start_date.dt.minute.data/60. + start_date.dt.second.data/3600.
            if (ds_timeperiod >= fp_timeperiod) or (resample_to_data == True):
                resample_period = str(round(fp_timeperiod/3600e9,5))+'H' # rt17603: Added 24/07/2018 - stops pandas frequency error for too many dp.
                site_fp = site_fp.resample(resample_period, dim='time', how='mean', base=base)
            elif ds_timeperiod < fp_timeperiod or (resample_to_data == False):
                resample_period = str(round(fp_timeperiod/3600e9,5))+'H' # rt17603: Added 24/07/2018 - stops pandas frequency error for too many dp.
                site_ds = site_ds.resample(resample_period, dim='time', how='mean', base=base)
                        
            site_ds = combine_datasets(site_ds, site_fp,
                                       method = "ffill",
                                       tolerance = tolerance)
            
            #transpose to keep time in the last dimension position in case it has been moved in resample
            if 'H_back' in site_ds.dims.keys():
                site_ds = site_ds.transpose('height','lat','lon','lev','time', 'H_back')
            else:
                site_ds = site_ds.transpose('height','lat','lon','lev','time')
                
            # If units are specified, multiply by scaling factor
            if ".units" in attributes:
                site_ds.update({'fp' : (site_ds.fp.dims, site_ds.fp / data[".units"])})
#                if calc_bc:
#                    for key in site_ds.keys():
#                        if "vmr" in key:
#                            site_ds.update({key :
#                                            (site_ds[key].dims, site_ds[key] / \
#                                            data[".units"])})
                if HiTRes:
                    site_ds.update({'fp_HiTRes' : (site_ds.fp_HiTRes.dims, 
                                                   site_ds.fp_HiTRes / data[".units"])})
        
            # Resample, if required
            if average[site] is not None:
                site_ds = site_ds.resample(average[site], dim = "time")
            
            fp_and_data[site] = site_ds
            
        
    if load_flux:
            
        flux_dict = {} 
            
        for source in emissions_name.keys():
            if type(emissions_name[source]) == str:
                flux_dict[source] = flux(domain, emissions_name[source], start=flux_bc_start, end=flux_bc_end, flux_directory=flux_directory)
            elif type(emissions_name[source]) == dict:
                if HiTRes == False:
                    print("HiTRes is set to False and a dictionary has been found as the emissions_name dictionary value\
                          for source %s. Either enter your emissions names as separate entries in the emissions_name\
                          dictionary or turn HiTRes to True to use the two emissions files together with HiTRes footprints." %source)
                    return None
                else:
                    flux_dict[source] = flux_for_HiTRes(domain, emissions_name[source], start=flux_bc_start, end=flux_bc_end, flux_directory=flux_directory)
                        
        fp_and_data['.flux'] = flux_dict
            
        
    if load_bc:
            
        bc = boundary_conditions(domain, species, start=flux_bc_start, end=flux_bc_end, bc_directory=bc_directory)

        if  ".units" in attributes:
            fp_and_data['.bc'] = bc / data[".units"]               
        else:
            fp_and_data['.bc'] = bc
            
        
                                           
    # Calculate model time series, if required
    if calc_timeseries:
        if load_flux == False:
            print("Can't get modelled mole fraction timeseries because load_flux is set to False.")
        else:
            sites = [key for key in fp_and_data.keys() if key[0] != '.']
            sources = fp_and_data['.flux'].keys()
            for site in sites:                    
                for source in sources:
                    if type(fp_and_data['.flux'][source]) == dict:
                        fp_and_data[site]['mf_mod_'+source] = xr.DataArray(timeseries_HiTRes(fp_and_data[site],fp_and_data['.flux'][source], output_fpXflux=False), coords = {'time': fp_and_data[site].time})
                    else:
                        flux_reindex = fp_and_data['.flux'][source].reindex_like(fp_and_data[site], 'ffill')
                        if source == 'all':
                            fp_and_data[site]['mf_mod'] = xr.DataArray((fp_and_data[site].fp*flux_reindex.flux).sum(["lat", "lon"]), coords = {'time':fp_and_data[site].time})
                        else:
                            fp_and_data[site]['mf_mod_'+source] = xr.DataArray((fp_and_data[site].fp*flux_reindex.flux).sum(["lat", "lon"]), coords = {'time':fp_and_data[site].time})
        
    # Calculate boundary conditions, if required         
    if calc_bc:
        if load_bc == False:
            print("Can't get modelled baseline timeseries because load_bc is set to False.")
        else:
            sites = [key for key in fp_and_data.keys() if key[0] != '.']
            for site in sites:
                bc_reindex = fp_and_data['.bc'].reindex_like(fp_and_data[site], 'ffill')
                fp_and_data[site]['bc'] = (fp_and_data[site].particle_locations_n*bc_reindex.vmr_n).sum(["height", "lon"]) + \
                                            (fp_and_data[site].particle_locations_e*bc_reindex.vmr_e).sum(["height", "lat"]) + \
                                            (fp_and_data[site].particle_locations_s*bc_reindex.vmr_s).sum(["height", "lon"]) + \
                                            (fp_and_data[site].particle_locations_w*bc_reindex.vmr_w).sum(["height", "lat"])
        
    for a in attributes:
        fp_and_data[a] = data[a]
  
    return fp_and_data


def fp_sensitivity(fp_and_data, domain, basis_case,
                   basis_directory = None):
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
        dict (xarray.Dataset) : Same format as fp_and_data with sensitivity matrix added.
    """    
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']
    
    flux_sources = fp_and_data['.flux'].keys()
    
    if type(basis_case) is not dict:
        if len(flux_sources) == 1:
            basis_case = {flux_sources[0]:basis_case}
        else:
            basis_case = {'all':basis_case}
    
    if len(basis_case.keys()) != len(flux_sources):
        if len(basis_case.keys()) == 1:
            print("Using %s as the basis case for all sources" %basis_case[basis_case.keys()[0]])
        else:
            print("There should either only be one basis_case, or it should be a dictionary the same length\
                  as the number of sources.")
            return None
    
    
    for site in sites:
        
        for si, source in enumerate(flux_sources):
        
            if source in basis_case.keys():
                basis_func = basis(domain = domain, basis_case = basis_case[source], basis_directory = basis_directory)
            else:
                basis_func = basis(domain = domain, basis_case = basis_case['all'], basis_directory = basis_directory)
            
            if type(fp_and_data['.flux'][source]) == dict:
                if 'fp_HiTRes' in fp_and_data[site].keys():
                    site_bf = xr.Dataset({"fp_HiTRes":fp_and_data[site]["fp_HiTRes"],
                                         "fp":fp_and_data[site]["fp"]})
                    H_all_arr=timeseries_HiTRes(site_bf, fp_and_data['.flux'][source], output_TS = False, output_fpXflux = True)
                    H_all = xr.DataArray(H_all_arr, coords=[site_bf.lat, site_bf.lon, site_bf.time], dims = ['lat','lon','time'])
                else:
                    print("fp_and_data needs the variable fp_HiTRes to use the emissions dictionary with high_freq and low_freq emissions.")
        
            else:
                site_bf = combine_datasets(fp_and_data[site]["fp"].to_dataset(), fp_and_data['.flux'][source])
                H_all=site_bf.fp*site_bf.flux 
                H_all_arr = H_all.values
            
            H_all_v=H_all.values.reshape((len(site_bf.lat)*len(site_bf.lon),len(site_bf.time)))        
        
        
            if 'region' in basis_func.dims.keys():
            
                if 'time' in basis_func.basis.dims:
                    basis_func = basis_func.isel(time=0)
            
                site_bf = xr.merge([site_bf, basis_func])
                
                H = np.zeros((len(site_bf.region),len(site_bf.time)))
            
                base_v = site_bf.basis.values.reshape((len(site_bf.lat)*len(site_bf.lon), len(site_bf.region)))
            
                for i in range(len(site_bf.region)):
                    H[i,:] = np.sum(H_all_v*base_v[:,i,np.newaxis], axis = 0)
                
                if source == all:
                    region_name = site_bf.region
                else:
                    region_name = [source+'-'+reg for reg in site_bf.region.values]

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
                    region_name = range(1,np.max(site_bf.basis)+1)
                else:
                    region_name = [source+'-'+str(reg) for reg in range(1,np.max(site_bf.basis)+1)]

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
    
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']

    basis_func = basis_boundary_conditions(domain = domain,
                                           basis_case = basis_case, bc_basis_directory=bc_basis_directory)
# sort basis_func into time order    
    ind = basis_func.time.argsort()                                        
    timenew = basis_func.time[ind]
    basis_func = basis_func.reindex({"time":timenew})
    
    for site in sites:

        DS_particle_loc = xr.Dataset({"particle_locations_n":fp_and_data[site]["particle_locations_n"],
                                "particle_locations_e":fp_and_data[site]["particle_locations_e"],
                                "particle_locations_s":fp_and_data[site]["particle_locations_s"],
                                "particle_locations_w":fp_and_data[site]["particle_locations_w"],
                                "bc":fp_and_data[site]["bc"]})
        
        DS_temp = combine_datasets(DS_particle_loc, fp_and_data[".bc"], method='ffill')
                        
        DS = combine_datasets(DS_temp, basis_func, method='ffill')                                    

        DS = DS.transpose('height','lat','lon','region','time')

        part_loc = np.hstack([DS.particle_locations_n,
                                DS.particle_locations_e,
                                DS.particle_locations_s,
                                DS.particle_locations_w])
       
        vmr_ed = np.hstack([DS.vmr_n,
                           DS.vmr_e,
                           DS.vmr_s,
                           DS.vmr_w])
        
        bf = np.hstack([DS.bc_basis_n,
                        DS.bc_basis_e,
                        DS.bc_basis_s,
                        DS.bc_basis_w])
        
        H_bc = np.zeros((len(DS.coords['region']),len(DS.bc)))
        
        for i in range(len(DS.coords['region'])):
            reg = bf[:,:,i,:]
            H_bc[i,:] = np.sum((part_loc*vmr_ed*reg), axis=(0,1))
        
        sensitivity = xr.Dataset({'H_bc': (['region_bc','time'], H_bc)},
                                    coords = {'region_bc': (DS.coords['region'].values),
                                              'time' : (DS.coords['time'])})

        fp_and_data[site] = fp_and_data[site].merge(sensitivity)
    
    return fp_and_data


def merge_sensitivity(fp_data_H,
                      out_filename = None,
                      remove_nan = True):
    """
    The merge_sensitivity function outputs y, y_site, y_time in a single array for all sites
    (as opposed to a dictionary) and H and H_bc if present in dataset.
    
    Args:
        fp_data_H (dict)   : Output from footprints_data_merge() function. Dictionary of datasets.
        out_filename (str) : If specified the output will be writen to this filename. Otherwise values are
                             returned.
        remove_nan (bool)  : Whether to remove NaN values in y. Default = True.
        
    Returns:
        tuple : each variable as an array (y, y_error, y_site, y_time [...])
        
        H, H_bc are also returned if present otherwise None is returned in their place.
        e.g.
            (y, y_error, y_site, y_time, H, H_bc)
            (y, y_error, y_site, y_time, None, H_bc)
    """

    y = []
    y_error = []
    y_site = []
    y_time = []
    H = []
    H_bc = []
    
    sites = [key for key in fp_data_H.keys() if key[0] != '.']
    
    for si, site in enumerate(sites):
        
        if remove_nan:
            fp_data_H[site] = fp_data_H[site].dropna("time", how="all")
        
        y_site.append([site for i in range(len(fp_data_H[site].coords['time']))])
        y_time.append(fp_data_H[site].coords['time'].values)        
        
        if 'mf' in fp_data_H[site].data_vars:
            y.append(fp_data_H[site].mf.values)
        
            # Approximate y_error
            if "vmf" in fp_data_H[site].keys():
                y_error.append(fp_data_H[site].vmf.values)
            elif "dmf" in fp_data_H[site].keys():
                y_error.append(fp_data_H[site].dmf.values)
            else:
                print("Measurement error not found in dataset for site %s" %site)

        if 'H' in fp_data_H[site].data_vars:        
            # Make sure H matrices are aligned in the correct dimensions
            if fp_data_H[site].H.dims[0] == "time":
                H.append(fp_data_H[site].H.values)
            else:
                H.append(fp_data_H[site].H.values.T)
                
        if 'H_bc' in fp_data_H[site].data_vars:         
            if fp_data_H[site].H_bc.dims[0] == "time":
                H_bc.append(fp_data_H[site].H_bc.values)
            else:
                H_bc.append(fp_data_H[site].H_bc.values.T)


    out_variables = ()

    y_site = np.hstack(y_site)
    y_time = np.hstack(y_time)

    if len(y) > 0:
        y = np.hstack(y)
        out_variables += (y,)
    else:
        out_variables += (None,)
    
    if len(y_error) > 0:
        y_error = np.hstack(y_error)
        out_variables += (y_error,)
    else:
        out_variables += (None,)
        
    out_variables += (y_site, y_time)
    
    if len(H_bc) > 0:
        H_bc = np.vstack(H_bc)
        out_variables += (H_bc,)
    else:
        out_variables += (None,) 
    
    if len(H) > 0:
        H = np.vstack(H)
        out_variables += (H,)
    else:
        out_variables += (None,)

    # Save or return y, y_error, y_site, y_time, H_bc, H
    if out_filename is None:
        return out_variables
    else:
        with open(out_filename, "w") as outfile:
            pickle.dump(out_variables, outfile)
        print("Written " + out_filename)
        return out_variables


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
                                 "daytime"           : selects data between 1100 and 1500 UTC
                                 "daytime9to5"       : selects data between 0900 and 1700 UTC
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

    # Filter functions
    def daily_median(dataset, keep_missing=False):
        """ Calculate daily median """
        return dataset.resample("1D", "time", how = "median")
        
    def six_hr_mean(dataset, keep_missing=False):
        """ Calculate daily median """
        return dataset.resample("6H", "time", how = "mean")
    

    def daytime(dataset, site,keep_missing=False):
        """ Subset during daytime hours (11:00-15:00) """
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 11 and h <= 15]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
        
    def daytime9to5(dataset, site,keep_missing=False):
        """ Subset during daytime hours (9:00-17:00) """
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 9 and h <= 17]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def nighttime(dataset, site,keep_missing=False):
        """ Subset during nighttime hours (23:00 - 03:00) """
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 23 or h <= 3]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def noon(dataset, site,keep_missing=False):
        """ Select only 12pm data """
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h == 12]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)] 
        

    def pblh_gt_threshold(dataset,site, keep_missing=False):
        """
        Subset for times when boundary layer height > threshold
        Threshold needs to be set in dataset as pblh_threshold
        """
        threshold = dataset.pblh_threshold.values
        ti = [i for i, pblh in enumerate(dataset.PBLH) if pblh > threshold]
        
        if keep_missing:
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')
            
            dataarray_temp = mf_data_array[dict(time = ti)]   
            
            mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def local_lapse(dataset,site, keep_missing=False):
        """
        Subset for times when linear combination of lapse rate and local influence
        is below threshold.
        
        Both normalized by 500 m. Thus, for low inlets the local influence is more dominant.
        For higher inlets the vertical profile of potential temperature (stability) is more dominant. 
        
        This combination correlates to variation in mole fraction difference between inlets.
        """
        in_height = dataset.inlet
        lapse_norm = dataset.theta_slope*in_height/500.
        lr_norm = dataset.local_ratio*500./in_height
        comb_norm = lr_norm + lapse_norm
        cutoff=0.5
        ti = [i for i, lr in enumerate(comb_norm) if lr < cutoff]
        
        if len(ti) > 0:
            if keep_missing:
                mf_data_array = dataset.mf            
                dataset_temp = dataset.drop('mf')
                
                dataarray_temp = mf_data_array[dict(time = ti)]   
                
                mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                      coords = {'time' : (dataarray_temp.coords['time'])})
                
                dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
                return dataset_out
            else:
                return dataset[dict(time = ti)]   
        else:
            return None       

    def local_lapse_045(dataset,site, keep_missing=False):
        """
        Subset for times when linear combination of lapse rate and local influence
        is below threshold.
        
        Both normalized by 500 m. Thus, for low inlets the local influence is more dominant.
        For higher inlets the vertical profile of potential temperature (stability) is more dominant. 
        
        This combination correlates to variation in mole fraction difference between inlets.
        """
        in_height = dataset.inlet
        lapse_norm = dataset.theta_slope*in_height/500.
        lr_norm = dataset.local_ratio*500./in_height
        comb_norm = lr_norm + lapse_norm
        cutoff=0.45
        ti = [i for i, lr in enumerate(comb_norm) if lr < cutoff]
        
        if len(ti) > 0:
            if keep_missing:
                mf_data_array = dataset.mf            
                dataset_temp = dataset.drop('mf')
                
                dataarray_temp = mf_data_array[dict(time = ti)]   
                
                mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                      coords = {'time' : (dataarray_temp.coords['time'])})
                
                dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
                return dataset_out
            else:
                return dataset[dict(time = ti)]   
        else:
            return None                       
            
    def local_influence(dataset,site, keep_missing=False):
        """
        Subset for times when local influence is below threshold.       
        Local influence expressed as a fraction of the sum of entire footprint domain.
        """
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
                         "pblh_gt_threshold": pblh_gt_threshold,
                         "local_influence":local_influence,
                         "six_hr_mean":six_hr_mean,
                         "local_lapse":local_lapse,
                         "local_lapse_045":local_lapse_045}

    # Get list of sites
    sites = [key for key in datasets.keys() if key[0] != '.']
    
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
            Dictionary of xarray datasets containing footprints and other variables
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
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
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
        release_lons = [fp_data[key].release_lon for key in fp_data.keys() if key[0] != '.']     
        release_lats = [fp_data[key].release_lat for key in fp_data.keys() if key[0] != '.']
                      
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

                    dt = np.float(date - time[ti-1])/np.float(time[ti] - time[ti-1])
                    fp_ti_0 = fp_data[site][dict(time = ti-1)].fp.values.squeeze()
                    fp_ti_1 = fp_data[site][dict(time = ti)].fp.values.squeeze()
                    fp_ti = fp_ti_0 + (fp_ti_1 - fp_ti_0)*dt
                    
                    data += np.nan_to_num(fp_ti)
    
                    # Store release location to overplot later
                    if "release_lat" in fp_data[site].keys():
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
            if site in release_lon:
                if "platform" in site_info[site]:
                    color = rp_color[site_info[site]["platform"].upper()]
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
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
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
    
        sites = [key for key in fp_data.keys() if key[0] != '.']
        
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
        pbar = ProgressBar(maxval=len(times)).start()

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
                     
            pbar.update(ti)
            print("")
        pbar.finish()
    
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


class get_country:
  def __init__(self, domain, ocean=False, ukmo=False, uk_split=False, country_dir = None):
      
        if country_dir is None:
            countryDirectory=data_path +'NAME/countries/' 
        else:
            countryDirectory = country_dir
            
        if ocean is False:
            filename=glob.glob(countryDirectory + 
                 "/" + "country_" 
                 + domain + ".nc")
             
        else:
            if uk_split is True:
                filename=glob.glob(countryDirectory + 
                     "/" + "country-ukmo-split_"
                     + domain + ".nc")
            else:
                if ukmo is False:
                    filename=glob.glob(countryDirectory + 
                         "/" + "country_ocean_"
                         + domain + ".nc")
                else:
                    filename=glob.glob(countryDirectory + 
                         "/" + "country-ukmo_"
                         + domain + ".nc")
        
        f = nc.Dataset(filename[0], 'r')
    
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get country indices and names
        country = f.variables['country'][:, :]
        if (ukmo is True) or (uk_split is True):
            name_temp = f.variables['name'][:]  
            f.close()
            name=np.asarray(name_temp)
        
        else:
            name_temp = f.variables['name'][:]
            f.close()
    
            name=[]
            for ii in range(len(name_temp)):
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
