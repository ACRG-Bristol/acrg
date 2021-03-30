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
import os
import pandas as pd
import numpy as np
import xarray as xr
from acrg_countrymask import domain_volume
from acrg_config.paths import paths
import cams_inversion
import climatology
import CAMS_BC

data_path      = paths.data
cams_directory = os.path.join(data_path, "ECMWF_CAMS/CAMS_inversion/")

def add_full_time(ds_seasonal, start, end, verbose=False):
    '''
    Expand out climatology dataset to include full date range from start to end.
    Expect input to be a climatology dataset with dimension of "month". This dimension will be expanded
    and replaced with the "time" dimension.
    
    Args:
        ds_seasonal (xarray.Dataset)
            Monthly averages with 'month' coordinate
        start, end (str)
            Start and end date range required to be added to ds_seasonal
    
    Returns:
        xarray.Dataset:
            dataset with a time coordinate in the range start:end, and the month coordinate removed
    '''
    
    date_range = []
    start_year = int(start.split("-")[0])
    end_year   = int(end.split("-")[0])
    
    ss   = start
    year = start_year
    while year < end_year+1:
        end  = end if year==end_year else str(year+1)+"-01-01"
        date = np.arange(ss, end, dtype="datetime64[M]").astype("datetime64[ns]")
        
        if len(date) > 0:
            date_range.append(date)
        ss = end
        year+=1
    
    if verbose: print("date_range",date_range)
    
    for i, date in enumerate(date_range):
        
        ds_copy = ds_seasonal.copy()

        if verbose: print(f'date {date}')

        months = pd.DatetimeIndex(date).month
        if verbose: print("months",months)
        
        ds_copy = ds_copy.sel(**{"month":months})

        ds_copy = ds_copy.assign(**{"time":("month", date)})
        ds_copy = ds_copy.swap_dims({"month":"time"})
        ds_copy = ds_copy.drop("month")
        
        ds = ds_copy if i==0 else xr.concat([ds, ds_copy],dim="time")
        
#         if i==0:
#             ds = ds_copy
#         else:
#             ds = xr.concat([ds, ds_copy],dim="time")
        
    return ds

def makeCAMSBC(domain, start, end,
               species = 'ch4',
               outdir = None, cams_directory = cams_directory,
               clim_start = None, clim_end = None,
               make_climatology = False,
               verbose = False):
    '''
    Make boundary condition files from the CAMS inversion product
    
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
            CHECK ON THIS
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
        verbose (bool, optional)
            If False do not print updates
 
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
    
    fp_lat,fp_lon,fp_height = domain_volume(domain)
    
    outdir = os.path.join(data_path, f"LPDM/bc/{domain}/") if outdir is None else outdir
    
    cams_ds = cams_inversion.readCAMSInversion(clim_start, clim_end)
    
    # create climatology if required
    if make_climatology:
        cams_seasonal = climatology.seasonal_cycle(cams_ds)
        cams_ds       = add_full_time(cams_seasonal, start = start, end = end)
    
    date_range = np.arange(start, end, dtype="datetime64[M]")
    date_range = [np.datetime_as_string(date)+"-01" for date in date_range]
    date_range += [end]

    for start, end in zip(date_range[:-1], date_range[1:]):

        # Check if file exists so a current file doesn't get overwritten
        out_filename = cams_inversion.bc_filename(domain, species, start)

        if os.path.isfile(os.path.join(outdir, out_filename)):
            print(f'Boundary condition file {out_filename} already exists.')
            print('Delete old one first to replace it.')
            continue
        
        # select the data for the correct date range
        ds = cams_ds.sel(**{"time" : slice(start, end)}) 

        cams_inversion.create_CAMS_BC(ds, fp_lat, fp_lon, fp_height, start,
                                      species          = species,
                                      domain           = domain,
                                      outdir           = outdir,
                                      cams_directory   = cams_directory,
                                      verbose          = verbose,
                                      from_climatology = make_climatology)
   