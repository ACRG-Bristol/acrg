#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:00:20 2019

Use makeCAMSBC

@author: rt17603
"""

import cams_inversion
import climatology
import os
import pandas as pd
import numpy as np
import xarray as xr
from acrg_countrymask import domain_volume

data_path = os.getenv("DATA_PATH")
cams_directory = os.path.join(data_path,"ECMWF_CAMS/CAMS_inversion/")

def add_full_time(ds_seasonal,start,end):
    '''
    Expand out climatology dataset to include full date range from start to end.
    Expect input to be a climatology dataset with dimension of "month". This dimension will be expanded
    and replaced with the "time" dimension.
    '''
    
    date_range = []
    start_year = int(start.split("-")[0])
    end_year = int(end.split("-")[0])
    
    s = start
    year = start_year
    while year < end_year+1:
        if year == end_year:
            e = end
        else:
            e = str(year+1)+"-01-01"         
        dr = np.arange(s,e,dtype="datetime64[M]").astype("datetime64[ns]")
        if len(dr) > 0:
            date_range.append(dr)
        s = e
        year+=1
    
    print("date_range",date_range)
    
    for i,dr in enumerate(date_range):
        
        ds_copy = ds_seasonal.copy()

        print("dr",dr)

        months = pd.DatetimeIndex(dr).month
        print("months",months)
        
        ds_copy = ds_copy.sel(**{"month":months})

        ds_copy = ds_copy.assign(**{"time":("month",dr)})
        ds_copy = ds_copy.swap_dims({"month":"time"})
        ds_copy = ds_copy.drop("month")
        
        if i == 0:
            ds = ds_copy
        else:
            ds = xr.concat([ds,ds_copy],dim="time")

    return ds

def makeCAMSBC(domain,st_date,st_end,
                          clim_start=None,clim_end=None,
                          outdir=None,cams_directory=cams_directory):
    ''' Make boundary condition files from the CAMS inversion product"
    
    Args:
        domain (str): 
            Domain name, e.g., "EUROPE"
        st_date (str): 
            Start date of output, e.g., "2015-01-01"
        st_end (str): 
            End date of output, e.g., "2016-01-01"
        clim_start (str, optional): 
            Start date to average fields into a climatology, e.g., "2010-01-01"
            Default of none will assume st_date
        clim_end (str, optional): 
            End date to average fields into a climatology, e.g., "2010-01-01"    
            Default of none will assume st_end
        outdir (str, optional)
            Output directory to save output. Output will automatically be written to outdir/LPDM/bc/DOMAIN
            CHECK ON THIS
        cams_directory (str, optional)
            Location of CAMS inversion output
 
    Returns:
        netcdf files of monthly boundary condition curtain files corresponding to the edges of the corresponding NAME domain
        
    Example:
        makeCAMSBC(domain,output_start,output_end,clim_start,clim_end)
    
    Todo: 
        Add some error checking (e.g. check that domain is correct)
    '''
    if clim_start is None and clim_end is None:
        clim_start = st_date
        clim_end = st_end
    
    species="ch4"

    fp_lat,fp_lon,fp_height = domain_volume(domain)
    
    if outdir is None:
        outdir = os.path.join(data_path,"LPDM/bc/{}/".format(domain))
    
    cams_ds = cams_inversion.readCAMSInversion(clim_start,clim_end)
    cams_seasonal = climatology.seasonal_cycle(cams_ds)    
    cams_bc = add_full_time(cams_seasonal,start=output_start,end=output_end)

    date_range = np.arange(st_date,st_end,dtype="datetime64[M]")
    date_range = [np.datetime_as_string(date)+"-01" for date in date_range]
    date_range += [st_end]

    for start,end in zip(date_range[:-1],date_range[1:]):
        
        # Check if file exists so a current file doesn't get overwritten
        out_filename = cams_inversion.bc_filename(domain,species,start)

        if os.path.isfile(os.path.join(outdir,out_filename)):
            print('Boundary condition file {} already exists.'.format(out_filename))
            print('Delete old one first to replace it.')
            continue
        
        ds = cams_bc.sel(**{"time":slice(start,end)})

        cams_inversion.create_CAMS_BC(ds,fp_lat,fp_lon,fp_height,start,species=species,domain=domain,
                                      outdir=outdir,cams_directory=cams_directory)

if __name__=="__main__":

#     clim_start = "2010-01-01"
#     clim_end = "2019-01-01"

    clim_start = None
    clim_end = None
    
    output_start = "2015-01-01"
    output_end = "2016-01-01"

    domain = "SOUTHAFRICA"

    outdir = os.path.join(data_path)

    cams_bc = makeCAMSBC(domain,output_start,output_end,clim_start,clim_end,outdir=outdir)
    