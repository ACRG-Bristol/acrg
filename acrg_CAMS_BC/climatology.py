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
from errors import semi_corr


def da_semi_corr_error(da,p_bounds=[2.5,97.5],correlation_coeff=0.5,reduce_col="time"):
    '''
    Calculate semi-correalted errors across a DataArray object based on input percentiles.
    Must contain at least the 2 percentiles specified by p_bounds and the median (percentile=50.0).
    Expect dimension called "percentile".
    '''
    p_dim = "percentile"
    
    if p_dim in da.dims:
        for p in p_bounds:
            if p not in da[p_dim]:
                raise Exception("Percentile bound {} not found. Unable to calculate semi correlated error".format(p))
        
        med_p = 50.0
        if med_p not in da[p_dim]:
            raise Exception("Cannot calculate semi-correlated error without median ({} percentile) defined".format(med_p))
    else:
        raise Exception("Cannot calculate semi-correlated error without percentile dimension at the moment")

    ## Inherant assumptions being made here about the shape but not sure how to make better
    nreduce = len(da[reduce_col])

    diff_low = da.sel(**{"percentile":med_p}) - da.sel(**{"percentile":p_bounds[0]})
    diff_upp = da.sel(**{"percentile":p_bounds[1]}) - da.sel(**{"percentile":med_p})
    
    std = np.zeros((nreduce,2))
    std[:,0] = diff_low
    std[:,1] = diff_upp
    std = np.mean(std,axis=1) # 1D array of length ntime
    
    smoothed_uncertainty = semi_corr(std,correlation_coeff=correlation_coeff)
    
    da.values[:] = smoothed_uncertainty
    
    da = da.mean(dim=reduce_col)
    
    da.attrs["error_type"] = "semi-corr"
    da.attrs["error_type_description"] = "Combined error has been calculated as semi-correlated (correlation coefficient = {})".format(correlation_coeff)
    
    return da
    

def mean_month(ds,time_col="time",error_col=None,error_type="corr",correlation_coeff=0.5):
    '''
    Calculate mean across month axis only for each data variable. Do not collapse other axes.
    Used on Grouper object created by ds.groupby("month").
    
    error_col (str/bool/None) :
        Name of error data variable in input xarray.Dataset objects (str) 
        OR whether the input xarray.DatasArray objects are error columns (bool)
        Only needed if you want to calculated something other than a correlated error.
        Default = None.
    error_type (str) :
        Type of error to calculate if an error column is specified.
        At present, this can be one of: "corr" or "semi-corr".
        Correlated - just takes the mean across all the percentiles/
        Semi-correlated - creates a covariance matrix where the off-diagonal elements are a fraction
        of the variance based on the correlation coefficient..
        Default = "corr"
    correlation_coeff (float) :
        Correlation coeffient to use when calculating semi-correlated errors. (i.e. error_type="semi-corr").
        Should be between 0.0 (uncorrelated) and 1.0 (fully correlated).
        Default = 0.5
    
    '''
    
    if isinstance(ds,xr.core.dataset.Dataset):
        
        data_vars = list(ds.data_vars)
        if error_col is not None and error_type == "semi-corr":
            ds[error_col] = da_semi_corr_error(ds[error_col],correlation_coeff=correlation_coeff,time_col=time_col)
            data_vars.remove(error_col)
        
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
            if error_col is True and error_type == "semi-corr":
                ds = da_semi_corr_error(ds,correlation_coeff=correlation_coeff,reduce_col="month")
                if "percentile" in ds.dims:
                    ds = ds.mean(dim="percentile",keep_attrs=True)
            else:
                ds = ds.mean(dim="month",keep_attrs=True)

    return ds

def seasonal_cycle(ds,time_col="time",error_col=None,error_type="corr"):
    '''
    Calculate the monthly cycle across time period of input object.
    
    Args:
        ds (xarray.Dataset) :
            Dataset object to look at seasonal cycle. Must contain a time coordinate.
        time_col (str) :
            Name of time co-ordinate in dataset.
            Default = "time".
        error_col (str/bool/None) :
            Name of error data variable in input xarray.Dataset objects (str) 
            OR whether the input xarray.DatasArray objects are error columns (bool)
            Only needed if you want to calculated something other than a correlated error.
            Default = None.
        error_type (str) :
            Type of error to calculate if an error column is specified.
            At present, this can be one of: "corr" or "semi-corr".
            Correlated - just takes the mean across all the percentiles/
            Semi-correlated - creates a covariance matrix where the off-diagonal elements are some fraction
            of the variance (default is 50%).
            Default = "corr"
    
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
            if error_col is not None and error_type == "semi-corr":
                if "percentile" in dv_dims:
                    dv_dims.remove("percentile")
                
            dims_list.append(dv_dims)
    else:
        dims_list = [list(ds.dims)]
        if time_col in ds.dims:
            i = ds.dims.index(time_col)
            dims_list[0][i] = "month"
        if error_col is not None and error_type == "semi-corr":
            if "percentile" in dims_list[0]:
                dims_list[0].remove("percentile")

    ds_new = ds.copy(deep=True)
    
    ds_new = ds_new.assign_coords(**{"month":(time_col,months)})
    ds_new = ds_new.swap_dims({time_col:"month"})
    
    group = ds_new.groupby("month")
    
    ds_mean = group.apply(mean_month,**{"time_col":time_col,"error_col":error_col,"error_type":error_type})
    
    if isinstance(ds,xr.core.dataset.Dataset):
        for i,dv in enumerate(ds.data_vars):
            ds_mean[dv] = ds_mean[dv].transpose(*dims_list[i])
    else:
        ds_mean = ds_mean.transpose(*dims_list[0])
    
    ds_mean.attrs["timeframe"] = "Climatology over time range: {} - {}".format(time[0].astype(str),time[-1].astype(str))
    
    return ds_mean

def calc_sector_seasonal(sector_flux,error_col=None,error_type="corr"):
    '''
    Calculate seasonal cycle for each Dataset in sector_flux dictionary
    
    Args:
        sector_flux (dict/collection.OrderdDict) :
            Dictionary containing datasets of values. Expect this to be split by sector names.
            Should match format of output from sector_prior() and sector_posterior() functions 
            (or from parse_sector_file() if reading from output sector total file).
        error_col (str/bool/None) :
            Name of error data variable in input xarray.Dataset objects (str) 
            OR whether the input xarray.DatasArray objects are error columns (bool)
            Only needed if you want to calculated something other than a correlated error.
            Default = None.
        error_type (str) :
            Type of error to calculate if an error column is specified.
            At present, this can be one of: "corr" or "semi-corr".
            Correlated - just takes the mean across all the percentiles/
            Semi-correlated - creates a covariance matrix where the off-diagonal elements are some fraction
            of the variance (default is 50%).
            Default = "corr"
    
    Returns:
        dict (xarray.Dataset) :
            Same format as passed to function. See seasonal_cycle() for format of individual Datasets.
    '''
    season_sector_dict = OrderedDict([])

    for sector,ds in sector_flux.items():
        ds_season = seasonal_cycle(ds,error_col=error_col,error_type=error_type)
        season_sector_dict[sector] = ds_season
    
    return season_sector_dict



def mean_season(da):
    '''
    Function which is applied to xarray.DataArray object grouped by season to extract mean of 
    one dimension rather than across the entire data array.
    Expect the input to have been grouped by month and to contain a dimension called "season".
    '''

    if "season" not in da.dims:
        # If only one value is grouped, no "month" dims will be present but may still exist within the 
        # DataArray.
        # In this case don't need to average but do want to make sure "month" was a part of this 
        # DataArray before grouping. Expect KeyError to be raised otherwise.
        # Seems like an xarray bug as should really be labelled as a dim.
        check = da["season"]
    else:
        da = da.mean(dim="season",keep_attrs=True)

    return da


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


def define_season_months(season_area="Amazon"):
    
    if season_area == "Amazon":
        months = OrderedDict([("wet",[2,3,4]),
                              ("dry_early",[5,6,7]),
                              ("dry",[8,9,10]),
                              ("wet_early",[11,12,1])
                             ])        
    else:
        months = OrderedDict([("winter",[12,1,2]),
                               ("spring",[3,4,5]),
                               ("summer",[6,7,8]),
                               ("autumn",[9,10,11])
                              ])
    return months


def assign_season(ds,name="season",season_months=None,season_area="Amazon",
                  set_as_coord=True,time_col="time",format_time=False):
    '''
    Assign season values to dataset based on the time column. Seasons are defined based on the month values
    using season_months input (dictionary of season to month numbers) or using the season_area keyword.
    
    Note:
        If some dates cannot be assigned to a season, these will be labelled as "unknown".
    
    Args:
        ds (xarray.Dataset) :
            Input Dataset (or dataarray?). Should contain a time axis (name defined by time_col input).
        name (str, optional) :
            Name to use for new season co-ordinate.
            Default = "season"
        season_months (dict, optional):
            Dictionary relating season names to months. Season names can be anything but
            months should be represented as numbers 1-12 representing Jan-Dec.
            Default = None. (area input will be used to define this)
        season_area (str, optional) :
            Area name to use when extracted associated months and season names.
            See define_season_months() function for options.
            Either an area or an explicit season_months dictionary should be specified where season_months
            will take precedence.
            Default = "Amazon"
        set_as_coord (bool, optional) :
            Whether to set new season column as the coordinate axis instead of time.
            Default = True
        time_col (str, optional) :
            Name of time dimension within ds variable.
            Default = "time"
        format_time (bool, optional):
            If time is in new NoLeap format, reformat this to match accepted time values to be recognised by
            pandas.
            Default = False
    
    Returns:
        (xarray.Dataset,numpy.array):
            Dataset with new co-ordinate for the season.
            Season names contained within the dataset including "unknown" if some dates could not be assigned.
    '''

    if format_time:
        time = reformat_time(ds[time_col])
    else:
        time = ds[time_col].values

    if season_months is None:
        season_months = define_season_months(season_area)
    
    all_season_names = list(season_months.keys())
    
    month_values = pd.DatetimeIndex(time).month
    
    max_s_len = max([len(sn) for sn in all_season_names])
    season = np.zeros(len(month_values),dtype=f"U{max_s_len}")
    for s,months in season_months.items():
        for m in months:
            index = np.where(month_values == m)[0]
            if len(index) > 0:
                season[index] = s

    season_names = np.unique(season)
    season_names = season_names[season_names!='']
    if '' in season:
         unknown = "unknown"
         season = np.where(season=='',unknown,season)
    
    ds = ds.assign_coords(**{name:(time_col,season)})
    if set_as_coord:
        ds = ds.swap_dims({time_col:name})
    
    return ds,season_names

    
def apply_seasonal_mean(ds,season_months=None,season_area="Amazon",select_season=None,
                        drop_unknown=False,time_col="time",format_time=False):
    '''
    Apply mean across grouped season values. Season is decided based on input months. Will be labelled based
    on season_name. At the moment this will only group into two seasons: the name defined by season_name 
    (e.g. "wet") and "other".
    
    Note: any months which cannot be assigned to a season will be labelled as "unknown". These entries will 
    be removed if drop_unknown=True.
    
    Args:
        ds (xarray.Dataset) :
            Dataset (or dataarray?) to average across. Should contain a time axis.
        season_months (dict, optional):
            Dictionary relating season names to months. Season names can be anything but
            months should be represented as numbers 1-12 representing Jan-Dec.
            Default = None. (area input will be used to define this)
        season_area (str, optional) :
            Area name to use when extracted associated months and season names.
            See define_season_months() function for options.
            Either an area or an explicit season_months dictionary should be specified where season_months
            will take precedence.
            Default = "Amazon"
        select_season (str/list/None, optional):
            Filter to include details for a subset of seasons.
            Default = None
        drop_unknown (bool, optional) :
            If any dates cannot be matched with a season these will be dropped.
            Default = False.
        time_col (str, optional) :
            Name of time dimension within ds variable.
            Default = "time"
        format_time (bool, optional):
            If time is in new NoLeap format, reformat this to match accepted time values to be recognised by
            pandas.
            Default = False

    Returns:
        xarray.Dataset:
            Dataset with new "season" coordinate. Values will be averaged across this dimension.
    '''
    
    name = "season"
    ds,season_names = assign_season(ds,name=name,season_months=season_months,season_area=season_area,
                                    time_col=time_col,format_time=format_time)
    
    ds_group = ds.groupby(name)
    ds_season_mean = ds_group.apply(mean_season)

    if drop_unknown:
        ds_season_mean = ds_season_mean.sel({name:season_names},drop=True)
    
    if select_season is not None:
        print("select_season",select_season)
        ds_season_mean = ds_season_mean.sel({name:select_season},drop=True)

    return ds_season_mean


def apply_seasonal_mean_simple(ds,months=[1,2,3],season_name="wet",time_col="time",format_time=True):
    '''
    ** CAN BE LARGELY DEPRECEATED IN FAVOUR OF apply_seasonal_mean() FUNCTION. **
    
    Apply mean across grouped season values. Season is decided based on input months. Will be labelled based
    on season_name. At the moment this will only group into two seasons: the name defined by season_name 
    (e.g. "wet") and "other".
    
    Args:
        ds (xarray.Dataset) :
            Dataset (or dataarray?) to average across. Should contain a time axis.
        months (iterable):
            Months to include in the main season to keep. Other months with be labelled "other".
            Default = [1,2,3]
        season_name (str) :
            Name of the season covered by months.
            Default = "wet"
        time_col (str, optional) :
            Name of time dimension within ds variable.
            Default = "time"
        format_time (bool, optional):
            If time is in new NoLeap format, reformat this to match accepted time values to be recognised by
            pandas.
            Default = False

    Returns:
        xarray.Dataset:
            Dataset with new "season" coordinate split into season_name (e.g. "wet") and "other". Values will
            be averaged across these dimensions.
    '''
    if format_time:
        time = reformat_time(ds[time_col])
    else:
        time = ds[time_col].values

    month_values = pd.DatetimeIndex(time).month
    season = np.array([season_name if month in months else "other" for month in month_values])
    
    print("month_values",month_values)
    
    ds = ds.assign_coords(**{"season":(time_col,season)})
    ds = ds.swap_dims({time_col:"season"})

    print("ds",ds)
    print("ds.season",ds.season)
    
    ds_group = ds.groupby("season")
    ds_season_mean = ds_group.apply(mean_season)
    ds_season_mean = ds_season_mean.sel(**{"season":season_name})

    return ds_season_mean
