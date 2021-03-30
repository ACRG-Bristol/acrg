#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:28:40 2019

@author: rt17603
"""

import numpy as np
import pandas as pd
import xarray as xr
from collections import OrderedDict
from acrg_satellite.gosat import coord_order
import pdb
    
def mean_month(ds,time_col="time",):
    '''
    Calculate mean across month axis only for each data variable. Do not collapse other axes.
    Used on Grouper object created by ds.groupby("month").
    '''
    
    if isinstance(ds,xr.core.dataset.Dataset):
        
        data_vars = list(ds.data_vars)
        
        for dv in data_vars:
            ds[dv] = ds[dv].mean(dim="month",keep_attrs=True)
        ds[time_col] = ds[time_col].swap_dims({"month":time_col})
            
        ds["month"] = ds["month"].mean(dim="month",dtype=np.int)

    else:
        if "month" not in ds.dims:
            # If only one value is grouped, no "month" dims will be present but may still exist within the 
            # DataArray.
            # In this case don't need to average but do want to make sure "month" was a part of this 
            # DataArray before grouping. Expect KeyError to be raised otherwise.
            # Seems like an xarray bug as should really be labelled as a dim.
            check = ds["month"]
        else:
            ds = ds.mean(dim="month",keep_attrs=True)

    return ds


def seasonal_cycle(ds,time_col="time"):
    '''
    Calculate the monthly cycle across time period of input object.
    
    Args:
        ds (xarray.Dataset) :
            Dataset object to look at seasonal cycle. Must contain a time coordinate.
        time_col (str) :
            Name of time co-ordinate in dataset.
            Default = "time".
    
    Output:
        xarray.Dataset :
            Dataset with each variable averaged for each repeated month.
            "time" dimension is replaced by "month".
            "time" coordinate still contained within dataset but no longer attached to any data variables.
    '''
    time = ds[time_col].values
    months = pd.DatetimeIndex(ds.time.values).month

    if isinstance(ds,xr.core.dataset.Dataset):
        dims_list = []
        for dv in ds.data_vars:
            dv_dims = list(ds[dv].dims)
            if time_col in dv_dims:
                i = dv_dims.index(time_col)
                dv_dims[i] = "month"
                
            dims_list.append(dv_dims)
    else:
        dims_list = [list(ds.dims)]
        if time_col in ds.dims:
            i = ds.dims.index(time_col)
            dims_list[0][i] = "month"

    ds_new = ds.copy(deep=True)
    
    ds_new = ds_new.assign_coords(**{"month":(time_col,months)})
    ds_new = ds_new.swap_dims({time_col:"month"})
    
    group = ds_new.groupby("month")
    
    ds_mean = group.apply(mean_month,**{"time_col":time_col})
    
    if isinstance(ds,xr.core.dataset.Dataset):
        for i,dv in enumerate(ds.data_vars):
            ds_mean[dv] = ds_mean[dv].transpose(*dims_list[i])
    else:
        ds_mean = ds_mean.transpose(*dims_list[0])
    
    ds_mean.attrs["timeframe"] = "Climatology over time range: {} - {}".format(time[0].astype(str),time[-1].astype(str))
    
    return ds_mean


def calc_sector_seasonal(sector_flux):
    '''
    Calculate seasonal cycle for each Dataset in sector_flux dictionary
    
    Args:
        sector_flux (dict/collection.OrderdDict) :
            Dictionary containing datasets of values. Expect this to be split by sector names.
            Should match format of output from sector_prior() and sector_posterior() functions 
            (or from parse_sector_file() if reading from output sector total file).
    
    Returns:
        dict (xarray.Dataset) :
            Same format as passed to function. See seasonal_cycle() for format of individual Datasets.
    '''
    season_sector_dict = OrderedDict([])

    for sector,ds in sector_flux.items():
        ds_season = seasonal_cycle(ds)
        season_sector_dict[sector] = ds_season
    
    return season_sector_dict


def reformat_time(time):
    '''
    Reformat time axis which is in the new NoLeap format to an array containing np.datetime64 values.
    
    Args:
        time (xarray.DataArray):
            Time axis as a DataArray object.
            TODO: Can update this to allow arrays if needs be as well.
    
    Returns:
        numpy.array:
            Numpy array with time values in np.datetime64 format.
            TODO: May want to update this to allow xr.DataArray to be returned?
    '''
    
    time = np.array(time.values,dtype=str)
    time = np.array(time,dtype=np.datetime64)
    
    return time

