#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:54:43 2020

@author: lw13938
"""
import numpy as np
import acrg_obs as getobs
import pandas as pd
import xarray as xr

def opends(fn):
    '''
    Open a netcdf dataset with xarray
    Args:
        fn (string)
            netcdf file to be opened
    Returns:
        xarray.Dataset: 
            netcdf file as  dataset 
    '''
    with xr.open_dataset(fn) as load:
        ds = load.load()
        return ds

def addaveragingerror(fp_all, sites, species, start_date, end_date, meas_period,  
                      inlet=None, instrument=None, obs_directory=None):
    """
    Adds the variablility within the averaging period to the mole fraction error.
    
    Args:
        fp_all (dict):
            Output from footprint_data_merge
        sites (list):
            List of site names
        species (str):
            Species of interest
        start_date (str):
            Start time of inversion "YYYY-mm-dd"
        end_date (str):
            End time of inversion "YYYY-mm-dd"
        meas_period (list):
            Averaging period of measurements
        inlet (str/list, optional):
            Specific inlet height for the site (must match number of sites).
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites).
        obs_directory (str, optional):
            Directory containing the obs data (with site codes as subdirectories)
            if not default.    
            
    Returns:
        fp_all (dict):
            fp_all from input with averaging error added to the mole fraction
            error
    """
    #Add variability in measurement averaging period to repeatability 
    dataerr = getobs.get_obs(sites, species, start_date = start_date, end_date = end_date,  
                          keep_missing=False,inlet=inlet, instrument=instrument,
                          data_directory=obs_directory)
    for si, site in enumerate(sites):
        if min(dataerr[site].index) > pd.to_datetime(start_date):
            dataerr[site].loc[pd.to_datetime(start_date)] = \
                [np.nan for col in dataerr[site].columns]           
        # Pad with an empty entry at the end date
        if max(dataerr[site].index) < pd.to_datetime(end_date):
            dataerr[site].loc[pd.to_datetime(end_date)] = \
                [np.nan for col in dataerr[site].columns]
        # Now sort to get everything in the right order
        dataerr[site] = dataerr[site].sort_index()
        fp_all[site].dmf.values = fp_all[site].dmf.values + dataerr[site].mf.resample(meas_period[si]).std(ddof=0).dropna().values
    return fp_all

def monthly_bcs(start_date, end_date, site, fp_data):
    """
    Creates a sensitivity matrix (H-matrix) for the boundary conditions, which
    will map monthly boundary condition scalings to the observations. This is 
    for a single site.
    
    Args:
        start_date (str):
            Start time of inversion "YYYY-mm-dd"
        end_date (str):
            End time of inversion "YYYY-mm-dd"
        site (str):
            Site that you're creating it for
        fp_data (dict):
            Output from acrg_name.bc_sensitivity
            
    Returns:
        Hmbc (array):
            Sensitivity matrix by month for observations
            
    """
    allmonth= pd.date_range(start_date, end_date, freq="MS")[:-1] 
    nmonth = len(allmonth)
    curtime = pd.to_datetime(fp_data[site].time.values).to_period("M")
    pmonth = pd.to_datetime(fp_data[site].resample(time="MS").mean().time.values)
    Hmbc = np.zeros((4*nmonth, len(fp_data[site].time.values)))
    cnt=0
    for cord in range(4):
        for m in range(0,nmonth):
            if allmonth[m] not in pmonth:
                cnt += 1
                continue
            mnth = allmonth[m].month
            yr = allmonth[m].year
            mnthloc = np.where(np.logical_and(curtime.month == mnth,curtime.year == yr))[0]
            Hmbc[cnt,mnthloc] = fp_data[site].H_bc.values[cord,mnthloc] 
            cnt += 1
    return Hmbc