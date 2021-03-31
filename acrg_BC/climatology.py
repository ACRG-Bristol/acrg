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
        ds_copy = ds_copy.drop_vars("month")
        
        ds = ds_copy if i==0 else xr.concat([ds, ds_copy],dim="time")
        
    return ds


def mean_month(ds,time_col="time",):
    '''
    Calculate mean across month axis only for each data variable. Do not collapse other axes.
    Used on Grouper object created by ds.groupby("month").
    
    Args:
        ds (xarray.Dataset)
            Emissions data 
        time_col (str)
            Name of the time coordinate
            Defaults to 'time'
    
    Returns:
        xarray.Dataset
            Average emissions per month
            "time" dimension is replaced by integer "month" coordinates.
    '''
    
    if isinstance(ds, xr.core.dataset.Dataset):
        
        data_vars = list(ds.data_vars)
        
        for dv in data_vars:
            if 'month' in ds[dv].dims:
                    ds[dv] = ds[dv].mean(dim="month",keep_attrs=True)
        ds[time_col] = ds[time_col].swap_dims({"month":time_col})
            
        ds["month"]  = ds["month"].mean(dim="month", dtype=np.int)

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

def monthly_cycle(ds, time_col="time"):
    '''
    Calculate the monthly cycle across time period of input object.
    
    Args:
        ds (xarray.Dataset)
            Dataset object to look at seasonal cycle. Must contain a time coordinate.
        time_col (str) 
            Name of time co-ordinate in dataset.
            Default = "time".
    
    Output:
        xarray.Dataset
            Dataset with each variable averaged for each repeated month.
            "time" dimension is replaced by "month".
            "time" coordinate still contained within dataset but no longer attached to any data variables.
    '''
    time = ds[time_col].values
    months = pd.DatetimeIndex(ds.time.values).month

    if isinstance(ds, xr.core.dataset.Dataset):
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

    ds_new  = ds.copy(deep=True)
    
    ds_new  = ds_new.assign_coords(**{"month":(time_col,months)})
    ds_new  = ds_new.swap_dims({time_col:"month"})
    
    group   = ds_new.groupby("month")
    ds_mean = group.map(mean_month,**{"time_col":time_col})
    
    if isinstance(ds,xr.core.dataset.Dataset):
        for i,dv in enumerate(ds.data_vars):
            if 'month' in ds_mean[dv].dims and 'month' not in dims_list[i]:
                dims_list[i].insert(1, 'month')
            ds_mean[dv] = ds_mean[dv].transpose(*dims_list[i])
    else:
        ds_mean = ds_mean.transpose(*dims_list[0])
    
    ds_mean.attrs["timeframe"] = f"Climatology over time range: {time[0].astype(str)} - {time[-1].astype(str)}"
    
    return ds_mean
