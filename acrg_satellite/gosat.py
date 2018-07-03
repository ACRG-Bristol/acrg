# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 10:25:39 2017

@author: rt17603

This module provides a set of functions for processing GOSAT data.

Within this module there are two main summary functions used for processing GOSAT CH4 Level 2 Data Product files:
    * gosat_process(directory, ...) - process a directory of files
    * gosat_process_file(filename, ...) - process an individual file

These two functions allow the following processes to be applied to the data:
 - data can be extracted for a certain date range
 - data can be cut down based on latitude and/or longitude upper and lower bounds
   OR
   data can be cut down based on a specified domain (latitude and longitude bounds will be extracted for 
   this domain if a footprint for this domain already exists)
 - data can be binned into specified latitude and longitude bin sizes
 - data can be filtered based on the quality flag present within the input file (performed by default)
 - data can be filtered to remove any data points containing a pressure_levels value of -9999.99 
 (performed by default)
 - data can be filtered based on "land" or "glint" mode
 - data can be filtered based on surface pressure from NAME (cutoff percent, layer order or extent of 
 surface layer).
 - data can be written out to netCDF files - one for each binned data point or one file per day.
 - data can be written out to csv files for input into NAME - one for each binned data point or one file per day.
 
Example 1:
    Want to extract data for India (site="GOSAT-INDIA")
    Want to extract data between dates 2012-01-01 - 2012-03-01
    Want to cut down the global data to:
        Latitude range of 6.0 - 36.5 degrees
        Longitude range of 68.0 - 93.0 degrees
    Want to bin the data into:
        Latitude bins of 0.234 degrees
        Longitude bins of 0.352 degrees
    Want to apply filters based on:
        - quality filter flag (only keep xch4_quality_filter = 0)
        - bad pressure filter (remove any data points where any pressure_levels value = -9999.99)
        - NAME surface pressure (filter out any values where GOSAT surface pressure is >5% different) (3 associated parameters)
    Want to use SOUTHASIA domain to access pressure files.
    Want to use surface pressure values from NAME rather than GOSAT values for level 0 (use_name_pressure )
    Want to write out one file per day
    Want to overwrite any files which may already be present

This could be run as:
    gosat_process("GOSAT-INDIA",
                  species="ch4",
                  start="2012-01-01",
                  end="2012-03-01",
                  lat_bounds=[6.0,36.5],
                  lon_bounds=[68.0,93.0],
                  coord_bin=[0.234,0.352],
                  quality_filt=True,
                  bad_pressure_filt=True,
                  pressure_domain="SOUTHASIA",
                  name_sp_filt=True,
                  filters=["cutoff"],
                  cutoff=5.0,
                  use_name_pressure=True
                  write_nc=True,
                  output_directory="/data/rt17603/obs/",
                  write_name=True,
                  name_directory="/data/rt17603/NAME_files/",
                  file_per_day=True,
                  overwrite=True)

Example 2:
    Want to extract data for EUROPE domain (GOSAT-EUROPE site) from one file into an xarray dataset 
    (no filtering)

This could be run as:
    ds = gosat_process_file("/data/shared/obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/2012/ESACCI-GHG-L2-CH4-GOSAT-OCPR-20121231-fv7.2.nc"
                            quality_filt=False,
                            bad_pressure_filt=False,
                            site="GOSAT-EUROPE",
                            write_nc=False,
                            write_name=False)
"""

import glob
import os
import re
import random
import math
import numpy as np
import pandas as pd
import xarray as xray
import datetime as dt
from collections import OrderedDict
import itertools
import acrg_agage as agage
from barometric import pressure_at_height
from acrg_countrymask import domain_volume

data_path = os.getenv("DATA_PATH")
home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/")
fp_directory = os.path.join(data_path,'NAME/fp/')
obs_directory = os.path.join(data_path,'obs/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"NAME/surface_pressure/")

def apply_filter(ds,filter_array,dim_apply='time'):
    '''
    The apply_filter function applies a filter array of indices to a dataset.
    All data variables with the dim_apply value as their first dimension will be filtered.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset to filter. 
        filter_array (np.array) : 
            Array of indices to keep within the dataset. All other indices will be removed.
        dim_apply (str, optional) : 
            First order dimension to apply the filter to.
            This was almost always be "time" (default) but is included here for generality.
            
    Returns:
        xarray.Dataset : 
            Copy of ds with data variables (data_vars) filtered.    
    '''
    
    data_vars = ds.data_vars
    ds = ds.copy(deep=True)

    # Apply filter to data variables
    if len(filter_array) > 0:
        for name in data_vars:
            try:
                dim_order = coord_order(ds,data_vars=[name])
                if dim_apply == dim_order["1"][0]:
                    ds[name] = ds[name][filter_array]
            except IndexError:
                raise IndexError("Indices from filter_array are out of bounds for data variable {0}.\n Filter array: {1}".format(name,filter_array))
            except KeyError:
                raise KeyError("Input for filter_array ({0}) should be an array of indices to keep. Current type: {1}".format(filter_array,type(filter_array)))
                
    else:
        # If filter is an empty array, we assume all elements should be removed.
        print('WARNING: Applying this filter has removed all elements.')
        for name in data_vars:
            dim_order = coord_order(ds,data_vars=[name])
            if dim_apply == dim_order["1"][0]:
                ds[name].values.fill(np.nan)

    ds = ds.dropna(dim_apply)
    
    return ds

def add_history_attr(ds,mod_attr):
    '''
    The add_history_attr adds or modifies the "history" attributes to include details of the modifications
    made to the original data.
    When history attribute is first created a timestamp is added of the current time.
    Note: Expect an additional attribute based on mod_attr has been added which contains details of the
    modification.
    
    Args:
        ds (xarray.Dataset) :
            Input Dataset.
        mod_attr (str) :
            Name of attribute which includes details of modification. This string will be added to the
            history attribute.
    Returns:
        xarray.Dataset:
            Input dataset with history attribute added or modified.
    '''
    if "history" in ds.attrs:
        ds.attrs["history"] += "{}, ".format(mod_attr)
    else:
        now = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%d %H:%M:%S")
        ds.attrs["history"] = "File modified on {} by University of Bristol ACRG group. Modification details listed within global attributes: {}, ".format(now,mod_attr)
    
    return ds

def gosat_quality_filter(ds):
    '''
    The gosat_quality_filter function filters all data variables within a dataset by the "xch4_quality_flag" variable.
    The xch4_quality_flag equates to:
        - 0 indicates good_quality
        - 1 indicates potentially_bad_quality
    Note: the Dataset is filtered assuming all data variables have the same first dimension (will be true for GOSAT data
    but would need to check for other data)

    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
                                
    Returns:
        ds (xarray.Dataset) : 
            Filtered Dataset with indices associated with xch4_quality_flag=1 removed.
    '''

    dim_order,consistent = coord_order(ds,check_consistent=True)
    if not consistent:
        print('WARNING: Order of dimensions is not consistent. Unable to apply gosat_quality_filter function')
        return None
    
    species = 'ch4'
    filter_name = "x" + species.lower() + "_quality_flag"
    flag = 0
    #dim_apply = "time"

    # Filter based on filter_name and flag
    try:
        # Note: adds dimension associated with where condition if not already
        # present. Should be fine for GOSAT data (all variables should have
        # a first order of "time") but may need to changed for other data.
        ds_new = ds.where(ds[filter_name] == flag,drop=True)
        #filt = np.where(ds[filter_name] == flag)
    except KeyError:
        raise KeyError("Unable to apply filter based on quality flag. Input dataset does not contain data variable {0}".format(filter_name)) 
    
    # Add attribute describing modification made to original data
    mod_attr = "quality_filter"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Original data has been filtered using {0} variable to include only flag = {1} (indicates good data).".format(filter_name,flag)

        
    return ds_new

def gosat_mode_filter(ds,mode='land'):
    '''
    The gosat_mode_filter function filters all data variables within a dataset by the mode as determined by the "retr_flag" variable.
    Mode is either "land" or "glint" which corresponds to 0 and 1 within  "retr_flag" data variable.
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
        mode (str, optional) : 
            Mode to include. Should be one of "land" or "glint". Default="land"
    
    Returns:
        ds (xarray.Dataset) : 
            Filtered Dataset with indices relevant to specified mode remaining.
    '''
   
    filter_name = "retr_flag"
    
    if mode == 'land':
        flag = 0
    elif mode == 'glint':
        flag = 1
    else:
        print('WARNING: Did not recognise input for mode. Should be one of: {0} or {1}'.format('land','glint'))
        return None
    
    # Filter based on filter_name and land or glint flag
    try:
        # Note: adds dimension associated with where condition if not already
        # present. Should be fine for GOSAT data (all variables should have
        # a first order of "time") but may need to changed for other data.
        ds_new = ds.where(ds[filter_name] == flag,drop=True)
    except KeyError:
        raise Exception("Unable to apply filter based on mode. Input dataset does not contain data variable {0}".format(filter_name)) 

    # Add attribute describing modification made to original data
    mod_attr = "retrieval_flag"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Original data has been filtered to only include retrievals from mode={0}. ({1} = {2})".format(mode,filter_name,flag)
    
    return ds_new

def gosat_pressure_filter(ds):
    '''
    The gosat_pressure_filter function removes all points which have at least one pressure_level undefined (-9999.99).
    Unable to use these points to calculate relevant pressure levels and difference between levels within the atmosphere.
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
            Should contain "pressure_levels" as a data variable.
            
    Returns:
       ds (xarray.Dataset) : 
           Filtered Dataset with data points with undefined pressure levels removed. 
    '''
    
    filter_name = "pressure_levels"
    
    exclude = -9999.99
    
    # Filter based on points where any pressure_levels value is -9999.99
    try:
        # Note: adds dimension associated with where condition if not already
        # present. Should be fine for GOSAT data (all variables should have
        # a first order of "time") but may need to changed for other data.
        ds_new = ds.where(~(ds[filter_name] == exclude).any(axis=1),drop=True) 
    except KeyError:
        raise Exception("Unable to apply filter to {0} to exclude values of {1}. Input dataset does not contain data variable {0}".format(filter_name,exclude)) 

    # Add attribute describing modification made to original data
    mod_attr = "bad_pressure_filter"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Data points with values of {1} for any of the {0} values have been removed.".format(filter_name,exclude)

    return ds_new
        
 
def latlon_filter(ds,lat_bounds,lon_bounds,columns=["latitude","longitude"]):
    '''
    The latlon_filter function filters a Dataset by latitude and longitude bounds.
    The columns related to latitude and longitude should be data variables within the input Dataset.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing latitude and longitude data variables.
        lat_bounds (iterable) : 
            Upper and lower bounds of latitude in degrees to include (two item tuple e.g. (6.0,33.5))
        lon_bounds (iterable) : 
            Upper and lower bounds of longitude in degrees to include (two item tuple e.g. (23.0,71,4))
        columns (iterable, optional) : 
            Data variable names for latitude and longitude (two item object ["latitude","longitude"] by default)
    
    Returns:
        ds (xarray.Dataset) : 
            Filtered Dataset only including data between lat_bounds and lon_bounds        
    '''
    
    if lat_bounds[0] > lat_bounds[1]:
        print('WARNING: Latitude bounds have been specified as lower bound = {0}, above upper bound = {1}. Switching around.'.format(*lat_bounds))
        lat_bounds = [lat_bounds[1],lat_bounds[0]]
    if lon_bounds[0] > lon_bounds[1]:
        print('WARNING: Longitude bounds have been specified as lower bound = {0}, above upper bound = {1}. Switching around.'.format(*lon_bounds))
        lon_bounds = [lon_bounds[1],lon_bounds[0]]
    
    try:
        # Note: adds dimension associated with where condition if not already
        # present. Should be fine for GOSAT data (all variables should have
        # a first order of "time") but but may need to changed for other data.
        ds_new = ds.where((ds[columns[0]] >= lat_bounds[0]) & (ds[columns[0]] < lat_bounds[1]) & (ds[columns[1]] >= lon_bounds[0]) & (ds[columns[1]] < lon_bounds[1]),drop=True)
    except KeyError:
        if columns[0] not in ds.data_vars:
            raise KeyError("Unable to apply lat-lon filter. Input dataset does not contain latitude data variable {0}.".format(columns[0]))
        elif columns[1] not in ds.data_vars:
            raise KeyError("Unable to apply lat-lon filter. Input dataset does not contain longitude data variable {0}.".format(columns[1]))

    # Add attribute describing modification made to original data
    mod_attr = "latlon_filter"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Original data has been filtered to include latitude range {0} - {1} and longitude range {2} - {3} degrees".format(lat_bounds[0],lat_bounds[1],lon_bounds[0],lon_bounds[1])
    
    return ds_new
    
 
def coord_order(ds,data_vars=[],check_consistent=False):
    '''
    The coord_order function determines the order of the coordinates for the data variables.
    
    For example for the xarray.Dataset:
        <xarray.Dataset>
            Dimensions:                (level: 3, time: 10)
            Coordinates:
              * time                   (time) datetime64[ns] 2010-01-01T01:30:00 ...
              * level                  (level) int64 0 1 2
            Data variables:
                lat                    (time) float64 nan nan 30.67 31.0 31.33 31.67 ...
                xch4_averaging_kernel  (time, level) float64 nan nan nan nan nan nan 0.0 ...
    
    The two coordinate axes are "time" and "level". This function would determine that "time" is only ever a first
    order coordinate (only used for axis=0) and that "level" is a second order coordinate (only used for axis=1).
    
    Args:
        ds (xarray.Dataset) : 
            Dataset
        data_vars (list, optional) : 
            Which data variables within ds to check.
            If data_vars if left blank ([] by default) all data variables will be examined.
        check_consistent (bool, optional) : 
            If True: checks that each coordinate is only ever used at one order (e.g. time is only
            ever used as a first order coordinate, level as a second order coordinate etc.)
    
    Returns:
        If check_consistent == False:
            dict : orders as a keys each with a list of relevant coordinate names
                    e.g. {1: ["time"], 2: ["level"]}
        If check_consistent == True:
            dict, bool : dictionary as described above, boolean flag for whether order is consistent
    '''
    if not data_vars:
        data_vars = ds.data_vars # Extract data varaible names from Dataset
    dim_order = {} 
    consistent_dims = True
    for name in data_vars:
        dims = ds[name].dims
        for i,d in enumerate(dims):
            order = str(i+1) # Set order as 1 for first order, 2 for second order...
            if order in dim_order.keys(): # Check if key already exists
                if d not in dim_order[order]: # Check that dimension is not already present for that key
                    dim_order[order].append(d)
            else:
                dim_order[order] = [d] # Create key if it doesn't exist
            
            if check_consistent:
                for key in dim_order.keys():
                    if key != order: # Check all other keys except where we just added the dimension name
                        if d in dim_order[key]:
                            # If the dimension is being used for different orders by different data variables
                            # the dimensions are not consistent.
                            consistent_dims = False
    
    if check_consistent:
        return dim_order,consistent_dims
    else:
        return dim_order

def midpoint(datetimes,weight=False):
    '''
    The midpoint function calculates the midpoint between a set of dates. This can either be done between
    the earliest and latest dates within the list or taken recursively between each pair of dates.
    
    Args:
        datetimes (list/np.array) : 
            Datetime array. Each element should be a numpy.datetime64 object.
        weight (bool) : 
            If True: calculate the midpoint based on all elements in the list
            If False: calculate the midpoint based on the earliest and latest datetimes only.
    
    Return:
        numpy.datetime64 : midpoint
    '''
    
    datetimes = np.sort(datetimes) # Ensure datetimes are in ascending order
    if weight:
        dt_list = list(datetimes) # Create datetimes list (easier to manipulate size)
        while len(dt_list) > 1:
                # Repeatedly calculate midpoint between first and last element to create weighted midpoint
                dt_list[0] = dt_list[0] + (dt_list[-1] - dt_list[0])/2.
                dt_list.pop(-1) # Remove last element from list
        mid_point = dt_list[0]
    else:
        # Take the middle datetime between the earliest and latest dates in datetimes list
        mid_point = datetimes[-1] - (datetimes[-1] - datetimes[0])/2.
    
    # Add very small random value to avoid two midpoint values being the same (has happened in rare circumstances).
    #mid_point += np.timedelta64(random.randrange(-1000,1000,1),'us')

    return mid_point    

def calc_mid(ds,name,message=False):
    '''
    The calc_mid function returns either the midpoint for a set of datetimes or the mean on the first order dimension
    (axis=0) of an array.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing name data variable .
        name (str) : 
            Name of data variable within ds
        message (bool, optional) :
            If True, print a warning when the function is unable to find a mean based on the object type.
            Default = False
            
    Returns:
        mean or midpoint of ds[name] (same type as objects within ds[name] array)    
    '''
    try:
        ds[name]
    except KeyError:
        raise KeyError("Unable to calculate mean for data variable {0}. Not found within dataset".format(name))
        
    if isinstance(ds[name][0].values,np.datetime64):
        return midpoint(ds[name].values)
    else:
        try:
            return np.mean(ds[name],axis=0)
        except:
            if message:
                print('WARNING: Unable to find a mean for {0}, dtype = {1}'.format(name,ds[name][0].dtype))
            return None

def associated_error_name(ds,error_name,error_ident='uncertainty'):
    '''
    The associated_error_name function finds the name of a variable associated with an error from the name.
    Assumes name of variable and error is e.g. "xch4" and "xch4_uncertainty" where error_ident is the extra
    string added to the variable name.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing both error_name and associated variable
        error_name  (str) : 
            Name of error data variable
        error_ident (str, optional) : 
            string to be removed from error_name to extract data_variable.
            NOTE: Function assumes the original quantity with which the error is associated can be 
            extracted by stripping off the error_ident (e.g. "uncertainty"). E.g. for "xch4_uncertainty" 
            the name for the related quantity would be "xch4".
    Returns:
        str : 
            Variable name
        
        If extracted variable is not found within the dataset or new variable name matches error_name:
            returns None    
    '''
    variable = error_name.rstrip(error_ident)
    variable = variable.rstrip('_')

    if variable == error_name:
        print("WARNING: Unable to find variable related to {0} by removing error identifier = {1}.".format(error_name,error_ident))
        return None
    
    if variable in ds.data_vars:
        return variable
    else:
        print("WARNING: Unable to find variable related to {0}. Variable {1} does not exist.".format(error_name,variable))
        return None

#def calc_mid_err(ds,error_name,error_ident='uncertainty'):
def calc_mid_err(ds,name,error_name):
    '''
    The calc_mid_err function returns an error value appropriately scaled for multiple points.
    This uses the equation:
        (mean(error)**2 + std(variable)**2)**1/2
    
    Args:
        ds (xarray.Dataset) : 
            Dataset. Note that dimensions must be used consistently for each data variable
            E.g. "time" should always be used as the first dimension, "level" as the second (if present)
        name (str) : 
            Name of data variable for quantity related to error e.g. 'xch4'
        error_name (str) : 
            Name of the data variable for the error quantity e.g. 'xch4_uncertainty'
    
    Returns:
        xarray.DataArray (1 entry): 
            mean error of ds[name]    
    '''
    
    return np.sqrt(np.mean(ds[error_name])**2. + np.std(ds[name])**2.)

def concat_str(ds,name,separator=','):
    '''
    The concat_str function concatenates the passed data variable as a string (arguments separated by "separator" e.g. ',').
    e.g. if name="exposure_id" then if the dataset contains:
    array(['2011101600560430371003', '2011101600560430372000',
       '2011101600560430372001'], dtype='|S22')
    then a string would be returned of the form
        '2011101600560430371003,2011101600560430372000,2011101600560430372001'
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing data variable specified by name parameter.
        name (str) : 
            data variable name within ds
        separator (str, optional) : 
            Joining str for the string values when concatenated.
            Default = ','
    
    Returns:
        str : 
            Values from name data variable as a string
    '''
    #attrs = ds[name].attrs
    #out = ','.join(ds[name].astype(dtype=str).values)
    return separator.join(ds[name].astype(dtype=str).values)

radius = 6371 #km
#radius = 6367.5 #km
def distance_lat(distance,radius=radius):
    '''
    Calculate equivalent difference in latitude for a distance in km.
    
    Args:
        distance (int/float) : 
            Distance in km.
        radius (int/float, optional)   : 
            Radius of the Earth (default to 6371km)
        
    Returns:
        float : 
            Latitude difference in degrees
    '''
    # Haversine distance formula is given as distance = 2*asin(sqrt(a))*radius
    a = math.sin(distance/(2.*radius))*math.sin(distance/(2.*radius))
    
    # Full formula for a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    # When dlon=0, lon1=lon2, a reduces to a = sin(dlat/2)**2
    dlat = 2*math.asin(math.sqrt(a))
    dlat = np.abs(dlat)
    
    dlat = math.degrees(dlat)
    return dlat

def distance_lon(distance,lat,radius=radius):
    '''
    Calculate equivalent difference in longitude for a distance in km at a given latitude.
    
    Args:
        distance (int/float) : 
            Distance in km.
        lat (int/float) : 
            Latitude in degrees (between but not including -90 and +90. Note: max value = 89.995)
        radius (int/float) : 
            Radius of the Earth (default to 6371km)
        
    Returns:
        float : 
            Longitude difference in degrees
    '''
    lat = math.radians(lat)
    # Haversine distance formula is given as distance = 2*asin(sqrt(a))*radius
    a = math.sin(distance/(2.*radius))*math.sin(distance/(2.*radius))
    
    # Full formula for a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    # When dlat=0, lat1=lat2, a reduces to a = cos^2(lat)*sin^2(dlat/2)
    b = math.sqrt(a)/math.cos(lat)
    dlon = 2.*math.asin(b)
    dlon = np.abs(dlon)
    
    dlon = math.degrees(dlon)
    return dlon

def calc_dlat(ds,sounding=10.5,lat_column="latitude"):
    '''
    For a set of measurements to be averaged calculate the spread of latitude based on the maximum latitude difference
    between the points and the resolution (sounding) on each individual measurement.
        dlat = (sounding in degrees)*2 + (max(lat) - min(lat))
    
    Args:
        ds (xarray.Dataset) : 
            Dataset which contains a data variable of latitudes (specified by lat_column)
        sounding (float, optional) : 
            Footprint diameter for each individual measurement in km.
            Default = 10.5 km (this is the sounding for the GOSAT satellite)
        lat_column (str, optional) : 
            Name of data variable containing latitude values within ds.
            Default = "latitude"
    
    Returns:
        float : 
            dlat value
    '''
    dlat = distance_lat(sounding) # Full size of box around point
    lat_spread = max(ds[lat_column].values) - min(ds[lat_column].values)
    dlat += lat_spread
    
    return dlat

def calc_dlon(ds,sounding=10.5,columns=["latitude","longitude"]):
    '''
    For a set of measurements to be averaged calculate the spread of longitude based on the maximum longitude difference
    between the points and the resolution (sounding) on each individual measurement.
        dlon = (sounding in degrees at mean latitude)*2 + (max(lon) - min(lon))
    
    Args:
        ds (xarray.Dataset) : 
            Dataset which contains data variables of latitude and longitude values (specified by columns)
        sounding (float, optional) : 
            Footprint diameter for each individual measurement in km.
            Default = 10.5 km (this is the sounding for the GOSAT satellite)
        columns (list, optional) : 
            Names of the latitude and longitude data variables (two item list).
            Default = ["latitude","longitude"]
    
    Returns:
        float : 
            dlon value    
    '''
    lat_mid = np.mean(ds[columns[0]])
    dlon = distance_lon(sounding,lat_mid) # Full size of box around point

    lon_spread = max(ds[columns[1]].values) - min(ds[columns[1]].values)
    dlon += lon_spread
    
    return dlon
     

def mean_ds(ds,error_ident=["uncertainty"],ident='exposure_id',dlat="dlat",dlon="dlon",coord_columns=["latitude","longitude"]):
    '''
    The mean_ds function reduces down a Dataset to the mean or midpoint of the contained data variables.
    Note: Function expects all data variables within the dataset to contain the same first dimension (axis=0).
    Note: The mean is taken over the first dimension only and any higher orders will be retained.
    Note: For datetime objects the midpoint between the earliest and latest date will be calculated.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset with consistent dimensions for each data variable. Only one first order dimension should be present.
            E.g. "time" should always be used as the first dimension, "level" as the second (if present)
        error_ident (list, optional) : 
            Keywords to be used for calculting errors appropriately.
            NOTE: Function assumes the original quantity with which the error is associated can be 
            extracted by stripping off the error_ident (e.g. "uncertainty"). E.g. for "xch4_uncertainty" 
                        the name for the related quantity would be "xch4".
        ident (str, optional) : 
            Data variable to keep as an identifier for each point within original dataset. This will be 
            returned as a string concatention of the values.
                        Default = "exposure_id". If ident value is not present in dataset, no value will be used as identifier.
        dlat (str, optional) : 
            Data variable to treat as latitude spread value. Default="dlat". If value is not present in dataset, 
            no dlat value will be calcuated.
        dlon (str, optional) : 
            Data variable to treat as longitude spread value. Default="dlon". If value is not present in dataset, 
            no dlon value will be calcuated.
        coord_columns (list, optional) : 
            Data variable names associated with latitude and longitude values. Used to calculate dlat and dlon.
            Default = ["latitude","longitude"]. 
    
    Returns:
        xarray.Dataset : 
            Mean of all dataset values across first order dimension (axis=0)
    '''
    
    dim_order,consistent = coord_order(ds,check_consistent=True)
    if not consistent:
        print('WARNING: Order of dimensions is not consistent. Unable to apply mean_ds function')
        return None
    if len(dim_order["1"]) > 1:
        print('WARNING: More than one first order dimension within dataset ({0}). Unable to apply mean_ds function consistently.'.format(dim_order["1"]))
        return None
    
    #dim_apply = "time"
    
    data_vars = ds.data_vars
    
    # Error quantities and dlon, dlat values need to be calculated before other values are averaged.
    calculated = []
    attrs = {}
    for name in data_vars:
        attrs[name] = ds[name].attrs
        for error in error_ident:
            if name.find(error) != -1:
                variable = associated_error_name(ds,name,error)
                ds[name] = calc_mid_err(ds,variable,name)
                calculated.append(name)
        if name == dlat:
            ds[name] = calc_dlat(ds,lat_column=coord_columns[0])
            calculated.append(name)
        elif name == dlon:
            ds[name] = calc_dlon(ds,columns=coord_columns)
            calculated.append(name)
            
    for name in data_vars:
        #if coord_order(ds,[name])["1"][0] == dim_apply:
        if name not in calculated and name != ident:
            ds[name] = calc_mid(ds,name)
        elif name == ident:
            ds[name] = concat_str(ds,ident)
        ds[name].attrs = attrs[name]


    coords = ds.coords
    attrs = {}
    for name in coords:
        attrs[name] = ds[name].attrs
        #if name == dim_apply:
        if name in dim_order["1"]: # Only calculate mean over first order co-ordinates
            ds[name] = calc_mid(ds,name)
        ds[name].attrs = attrs[name]
    
    return ds


def gosat_add_coords(ds,data_vars=[]):
    '''
    The gosat_add_coords function adds co-ordinates to a GOSAT dataset.
    This is valid for (at least) v6 of GOSAT CH4 Proxy Level 2 Data Product downloaded from the CCI Open Data Portal
    (e.g. filename of the form: ESACCI-GHG-L2-CH4-GOSAT-OCPR-20101231-fv6.nc)
    
    This assumes the input Dataset will be of the form:
        <xarray.Dataset>
        Dimensions:                           (m: 20, n: 1962)
        Dimensions without coordinates: m, n
        Data variables:
            xch4_quality_flag                 (n) int8 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            ch4_profile_apriori               (n, m) float32 1751.95 1751.95 1751.94 ...
            xch4                              (n) float32 1644.78 1657.72 1648.04 ...
            ...
        Attributes:
            title:                     ESA CCI GOSAT OCPR CH4
            institution:               University of Leicester (UoL), UK
            ...
    
    i.e. with no coordinates assigned for the two dimensions within the Dataset.
    This function redefines "m" as "time" and "n" as "level" and adds the associated values as coords.

    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
        data_vars (list) : 
            Data variables which will have the time and level (where appropriate) dataset coords added
            If data_vars if left blank ([] by default) coords will be added to all data variables within ds.
                              
    Returns:
        ds (xarray.Dataset) : 
            Dataset with time and level coordinates assigned
            (only includes specified data variables unless left blank)
    '''
    
    dims_current = ['n','m']
    dims = ['time','lev']
    
    if ds.dims.keys() != dims_current and ds.dims.keys() != dims_current[::-1]: # Check dimensions match what we expect (check forward and reverse list)
        print('WARNING: Do not recognise dimensions of input gosat dataset. Unable to add new dimensions.')
        return None
        
    
    # Define co-ordinate values. Extract from Dataset if already present (e.g. time), construct if not (e.g. lev)
    coords = {}
    for i,d in enumerate(dims):
        try:
            coords[d] = ds[d].values # Extract values for each coord from Dataset if present (e.g. time)
        except (AttributeError,KeyError):
            current_dim = dims_current[i]
            dim_length = ds.dims[current_dim]
            coords[d] = np.arange(dim_length) # If no values defined, create integer list based on dimension size
    
    if not data_vars:
        data_vars = ds.data_vars
    
    data = OrderedDict()
    for name in data_vars:
        if name not in dims:
            if len(ds[name].dims) == 1:
                data[name] = (dims[0],ds[name].values)
            elif len(ds[name].dims) == 2:
                data[name] = ([dims[0],dims[1]],ds[name].values)
    
    #coord_dict = {dims[0]:dim1,dims[1]:dim2}
    
    ds_new = xray.Dataset(data,coords=coords,attrs=ds.attrs)
    for dv in ds_new.data_vars:
        ds_new[dv].attrs = ds[dv].attrs
    for coord in ds_new.coords:
        try:
            ds_new[coord].attrs = ds[coord].attrs
        except (AttributeError,KeyError):
            if coord == "lev" or coord == "level":
                ds_new[coord].attrs["short_description"] = "Number for each level within the vertically resolved data."
                ds_new[coord].attrs["long_name"] = "level"
    
    # Add attribute describing modification made to original data
    mod_attr = "nc_restructure"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Dataset rearranged to re-assign {} dimensions as {} coordinates, extracting values from dataset where present.".format(dims_current,dims).replace("'","")
    
    return ds_new

def create_zbin(ds,lat_bounds=[],lon_bounds=[],domain=None,columns=["latitude","longitude"],coord_bin=0.5):
    '''
    The create_zbin function adds an additional column to a dataset containing the midpoints of a specified lat-lon
    bin as a string. Each data point is assigned to a bin.
    The bin values will be of the form: "23.25,32.25" to indicate a position is within lat bin 23 - 23.5 and lon bin
    32 - 32.5. If a position is outside the latitude or longitude bounds the bin will contain "above" or "below" 
    e.g. "above,32.25" to indicate this position is above the upper latitude bound but within lon bin 32 - 32.5.
    
    Adding this data variable allows the data to be easily grouped within each lat-lon bin. 
    
    Note: if using imported GOSAT CH4 Proxy Level 2 Data Product (from the CCI Open Data Portal) you will need to
    have added the time and level coords using the gosat_add_coords() function before using this function.
    
    E.g. ds should be of the form:
        <xarray.Dataset>
        Dimensions:                           (level: 20, time: 1962)
        Coordinates:
          * level                             (level) int64 0 1 2 3 4 5 6 7 8 9 10 ...
          * time                              (time) datetime64[ns] 2010-12-31T00:31:51.000015 ...
        Data variables:
            xch4_quality_flag                 (time) int8 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            ch4_profile_apriori               (time, level) float32 1751.95 1751.95 ...
            xch4                              (time) float32 1644.78 1657.72 1648.04 ...  
            longitude                         (time) float32 152.298 152.301 152.305 ...
            latitude                          (time) float32 -70.4105 -70.4111 ...
            ...
   i.e. with coords assigned (in this case time and level)
    
    Args:
        ds (xarray.Dataset) : 
             Dataset with coordinates and dimensions assigned - see above.
        lat_bounds (list, optional) : 
            Upper and lower bounds for latitude bins. Upper bound is included in range.
        lon_bounds (list, optional) : 
            Upper and lower bounds for longitude bins. Upper bound is included in range.
        domain (str,optional) : 
            Domain to take lat_bounds and lon_bounds from. This will supercede any lat_bounds and lon_bounds 
            specified.
        columns (list, optional) : 
            column names within ds containing the latitude and longitude positions as a two item list
            Default = ["latitude","longitude"].
        coord_bin (float/list, optional) : 
            Size of each bin for latitude and longitude in degrees.
            To specify different bin sizes for latitude and longitude include a two item list 
            (e.g. [0.234,0.356])
            Default = 0.5.
    
    Returns:
        ds (xarray.Dataset) :
            Dataset with additional "zbin" data variable
    '''
    
    if domain:
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if not lat_bounds:
        lat_bounds = [min(ds[columns[0]]),max(ds[columns[0]])]
    if not lon_bounds:
        lon_bounds = [min(ds[columns[1]]),max(ds[columns[1]])]
    
    if isinstance(coord_bin,float) or isinstance(coord_bin,int):
        coord_bin = [coord_bin,coord_bin]
    elif len(coord_bin) == 1:
        coord_bin = [coord_bin[0],coord_bin[0]]
    
    if len(coord_bin) != 2:
        print("WARNING: Did not recognise input for coord_bin into create_zbin function: {0}. Should be a two item object".format(coord_bin))
        return None
    
        
    lat_bins = np.arange(lat_bounds[0],lat_bounds[1]+coord_bin[0],coord_bin[0])
    lon_bins = np.arange(lon_bounds[0],lon_bounds[1]+coord_bin[1],coord_bin[1])
    
    order = coord_order(ds)
    coord_1storder = order["1"][0]
    
    latitude = ds[columns[0]].values
    longitude = ds[columns[1]].values
    
    lat_bin_assign = np.digitize(latitude, lat_bins)
    lon_bin_assign = np.digitize(longitude, lon_bins)
    
    bins = []
    
    for lat_bin,lon_bin in zip(lat_bin_assign,lon_bin_assign):
        if lat_bin == 0:
            bin1 = 'below'
        elif lat_bin == len(lat_bins):
            bin1 = 'above'
        else:
            bin1 = lat_bins[lat_bin-1] + (lat_bins[lat_bin-1] - lat_bins[lat_bin-2])/2.
            #bin1 = lat_bin
        if lon_bin == 0:
            bin2 = 'below'
        elif lon_bin == len(lon_bins):
            bin2 = 'above'
        else:
            bin2 = lon_bins[lon_bin-1] + (lon_bins[lon_bin-1] - lon_bins[lon_bin-2])/2.
            #bin2=lon_bin

        bins.append("%s,%s"%(bin1,bin2))
    
    #bins = ["%s,%s"%(lat,lon) for lat, lon in zip(np.digitize(latitude, lat_bins),np.digitize(longitude, lon_bins))]
    bins = np.array(bins)
    
    #np.digitize - bins[i-1] <= x < bins[i]
    
    ds = ds.assign(zbin=(coord_1storder, bins))
    
    attributes = OrderedDict({})
    attributes["bin_size"] = "Binned to {0} degrees in latitude, {1} degrees in longitude".format(*coord_bin)
    attributes["latitude_bounds"] = "{0} - {1} degrees".format(lat_bounds[0],lat_bounds[1])
    attributes["longitude bounds"] = "{0} - {1} degrees".format(lon_bounds[0],lon_bounds[1])
    attributes["comment"] = "Bin values are the midpoints within each latitude and longitude bin." 
    
    ds["zbin"] = ds["zbin"].assign_attrs(attributes) # Add attributes to zbin coordinate
    
    return ds
 
def zbin_filter(ds,dim_apply="time"):
    '''
    The zbin_filter function uses the "zbin" column (added using create_zbin function) and removes any values which 
    fall outside the latitude or longitude range set when creating the column (for the specified data variables or 
    all if none are specified).
    
    See the create_zbin function for more information but a value is outside the range if the "zbin" column contains
    the word "above" or "below".
    
    Args:
        ds (xarray.Dataset) :
            Dataset which has "zbin" data variable created using create_zbin() function.
        dim_apply (str, time) : 
            Dimension to apply filter to. Default = "time"
                    
    Returns:
        ds (xarray.Dataset) : 
            Dataset with data variables outside range removed
    '''

    remove = ['above','below']
    filter_name = "zbin"
    #dim_apply = "time"
    
    attrs = ds[filter_name].attrs
    
    for string in remove:
        ds[filter_name] = ds[filter_name].astype(dtype=str) # Have to recast as explicit string type
        filt = np.where(np.char.find(ds[filter_name].values,string)==-1)
        ds = apply_filter(ds,filt,dim_apply=dim_apply)    
        #ds = ds.where(np.char.find(ds[filter_name].values,string)==-1,drop=True)
    
    ds[filter_name].attrs = attrs
    
    return ds    

def binned_mean(ds,lat_bounds=[],lon_bounds=[],domain=None,coord_bin=0.5,columns=["latitude","longitude"],add_spread=True):
    '''
    The binned_mean function finds the mean within bins specified by outer latitude and longitude bounds and 
    a given bin size for the input Dataset.
    
    Either lat_bounds AND lon_bounds OR domain must be specified.
    
    Args:
        ds (xarray.Dataset) :
            Dataset which contains latitude and longitude data variables (with names specified by columns).
        lat_bounds (list) : 
            Upper and lower bounds for latitude bins.
        lon_bounds (list) : 
            Upper and lower bounds for longitude bins.
            Either lat_bounds AND lon_bounds OR domain must be specified.
        domain (str) : 
            Domain to take lat_bounds and lon_bounds from. This will supercede any lat_bounds and lon_bounds 
            specified.
        coord_bin (float/list) : 
            Size of each bin in degrees for latitude and longitude
            To specify different bin sizes for latitude and longitude include a two item list 
            (e.g. [0.234,0.356])
            Default = 0.5
        columns (list, optional) : 
            Column names within ds containing the latitude and longitude positions.
            Default = ["latitude","longitude"]
        add_spread (bool, optional) : 
            Whether to add spread to dataset in the form of data variables called "dlat" and "dlon". 
            Default=True

    Returns:
        xarray.Dataset : 
            Dataset with mean values within each latitude and longitude bin            
    '''
    
    filter_name = "zbin"
    swap_dim = coord_order(ds)["1"][0] # Dimension to re-assign after grouping e.g. "time"
    
    ds = create_zbin(ds,lat_bounds,lon_bounds,domain,columns=columns,coord_bin=coord_bin)
    ds = zbin_filter(ds)

    if add_spread:
        dlat = "dlat"
        dlon = "dlon"
        
        t_vals = ds[swap_dim]
        dummy_vals = np.zeros(len(t_vals))
        dummy_vals.fill(np.nan)
        da = xray.DataArray(dummy_vals,coords={swap_dim:t_vals},dims={swap_dim:len(t_vals)})

        ds = ds.assign(**{dlat:da,dlon:da})
    
    if ds[swap_dim].size != 0: # xarray bug, if no elements for swap_dim (e.g. time) groupby returns a StopIteration
        ds = ds.groupby(filter_name).apply(mean_ds) # Creates xarray.Grouper object and then applies mean_ds function on those grouped datasets.
        
        if len(ds[swap_dim].dims) == 0: # If swap_dim has size 1 it loses dimensionality and filter_name dimension must be explicitly reassigned
            t = np.array([ds[swap_dim].values])
            ds = ds.drop(swap_dim)
            da = xray.DataArray(t,coords={filter_name:ds[filter_name].values},dims={filter_name:ds[filter_name].size})
            ds = ds.assign(**{swap_dim:da})
        
        ds = ds.swap_dims({filter_name:swap_dim}) # Reassign dimension e.g. "time" as the first dimension
        ds = ds.drop(filter_name) # Drop added filter_name (e.g. "zbin") data variable 
        
        if not lat_bounds:
            lat_bounds = [min(ds[columns[0]]),max(ds[columns[0]])]
        if not lon_bounds:
            lon_bounds = [min(ds[columns[1]]),max(ds[columns[1]])]
        
        # Add attribute describing modification made to original data
        mod_attr = "bins"
        ds = add_history_attr(ds,mod_attr)
        ds.attrs[mod_attr] = "Original data has been binned in latitude range: {0} - {1}, longitude range {2} - {3} in bins of {4} degrees.".format(lat_bounds[0],lat_bounds[1],lon_bounds[0],lon_bounds[1],coord_bin)
        
    else:
        print('MESSAGE: No values within lat and lon bounds to bin')
    
    return ds

    
def bin_check(filename,lat_bounds=[],lon_bounds=[],domain=None,coord_bin=0.5):
    '''
    The bin_check function uses a slower and simpler alternative method for binning (compared to xarray.Dataset.groupby function).
    This was built as a quick comparison method to check the output of method used in functions above.
    
    Args:
        filename (str) : 
            Filename of GOSAT CH4 Level 2 Data Product file (str)
        lat_bounds (list, optional) : 
            Upper and lower bounds for latitude. (tuple of values e.g. (6.0,36.5))
        lon_bounds (list, optional) : 
            Upper and lower bounds for longitude. (tuple of values e.g. (70.0,90.5))
        domain (str, optional) : 
            If lat_bounds and lon_bounds are not specified a domain can be specified instead.
        coord_bin (float/list) : 
            Size of each bin in degrees (default = 0.5).
    
    Returns:
        np.array (3) : binned latitude, binned longitude, binned xch4 column 
    '''
    
    gosat = gosat_process_file(filename,lat_bounds=[],lon_bounds=[],domain=None,coord_bin=None)
    latitude = gosat.latitude.values
    longitude = gosat.longitude.values
    xch4 = gosat['xch4'].values
    
    if domain:
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if not lat_bounds:
        lat_bounds = [min(gosat.latitude.values),max(gosat.latitude.values)]
    if not lon_bounds:
        lon_bounds = [min(gosat.longitude.values),max(gosat.longitude.values)]
    
    lat_range = np.arange(lat_bounds[0],lat_bounds[1]+coord_bin,coord_bin)
    lon_range = np.arange(lon_bounds[0],lon_bounds[1]+coord_bin,coord_bin)
    
    lat_bin = []
    lon_bin = []
    xch4_mean = []
    
    for i in range(len(lat_range)-1):
        lat1 = lat_range[i]
        lat2 = lat_range[i+1]
        for j in range(len(lon_range)-1):
            lon1 = lon_range[j]
            lon2 = lon_range[j+1]
            filt = np.where((latitude >= lat1) & (latitude < lat2) & (longitude >= lon1) & (longitude < lon2))
            extract_lat = latitude[filt]
            extract_lon = longitude[filt]
            extract_xch4 = xch4[filt]
            
            if extract_lat.any():
                mean_lat = np.mean(extract_lat)
                mean_lon = np.mean(extract_lon)
                mean_xch4 = np.mean(extract_xch4)
                lat_bin.append(mean_lat)
                lon_bin.append(mean_lon)
                xch4_mean.append(mean_xch4)
    
    lat_bin = np.array(lat_bin)
    lon_bin = np.array(lon_bin)
    xch4_mean = np.array(xch4_mean)
    
    return lat_bin,lon_bin,xch4_mean  

def units(name):
    '''
    The units function defines units for output variables.
    
    Data variable units defined within this function are:
        ch4_profile_apriori: ppb
        xch4: ppb
        xch4_uncertainty: ppb
        lat: degrees
        lon: degrees
        pressure_levels: hPa
        time: seconds since 1970-01-01 00:00:00
        pressure_weights: unitless
        xch4_averaging_kernel: unitless
    
    Args:
        name (str) : 
            Name of output variable
    
    Returns:
        str : 
            Unit of name (if present), None otherwise
    '''
    
    unit_data = {}
    
    unit_data["ch4_profile_apriori"] = "1e-9"
    unit_data["xch4"] = "1e-9"
    unit_data["xch4_uncertainty"] = "1e-9"
    
    unit_data["lat"] = "degrees"
    unit_data["lon"] = "degrees"
    
    unit_data["pressure_levels"] = "hPa"
    
    unit_data["time"] = "seconds since 1970-01-01 00:00:00"
    
    unit_data["pressure_weights"] = "unitless"
    unit_data["xch4_averaging_kernel"] = "unitless"
    
    try:
        unit = unit_data[name]
    except KeyError:
        return None
    else:
        return unit

def extract_dates(ds,dim="time",dtype='M8[D]'):
    '''
    The extract_dates function converts datetime objects into date strings for all values along a given dimension.
    Note: dimension should contain an array of np.datetime64 objects.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset which contains dimension specified by dim
        dim (str, optional): 
            String of dimension containing datetime objects. Default = "time"
        dtype (str, optional):
            Specific dtype value to use to cast datetime object. Default = 'M8[D]'
        
    Returns:
        np.array (str) : dates as an array of strings
    '''
    return ds[dim].values.astype(dtype).astype(str) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)

def extract_files(directory,search_str=None,start=None,end=None,date_separator=''):
    '''
    The extract_files function extracts filenames from a directory based on a search_str and/or based on a start
    and end dates (if filename contains details of a datestamp).
    
    Note:
        By default if using start and end dates files of the form "*YYYYMMDD*" are expected e.g. filename_20111001.csv.
        A date_separator can also be specified between YYYY/MM and MM/DD. 
        For example e.g. date_separator='-', will search for files of the form "*YYYY-MM-DD*".
    
    Args:
        directory (str) : 
            Directory to search (str)
        search_str (str/None, optional) : 
            String to use to search directory (using glob). If full filename is not specified this should 
            contain at least one wildcard character ("*").
            If no search_str specified, all files from directory will be returned and filtered by any start 
            and end criteria.
        start (str/None, optional) : 
            Start date of the form "YYYY-MM-DD".
        end (str/None, optional) : 
            End date of the form "YYYY-MM-DD". Date range is up to but not including this date.
            Either both or neither of start and end should be specified.
        date_separator (str) : 
            Date separator string between year, month and day in filename.
            By default this is "", meaning dates of the form "YYYYMMDD" will be searched for.
    
    Returns:
        list : 
            Filenames as string s(full path information).  
    '''
    
    if start:
        if start.find('-') == -1:
            print('WARNING: Start date to extract files is not in correct format: should be YYYY-MM-DD. Unable to set date range.')
            #start=None
            #end=None
            return None
    if end:
        if end.find('-') == -1:
            print('WARNING: End date to extract files is not in correct format: should be YYYY-MM-DD.  Unable to set date range.')
            #start=None
            #end=None
            return None
    
    if (start and not end) or (end and not start):
        print('WARNING: Start and end must both be specified to extract date range of files. Unable to set date range.')
        #start=None
        #end=None
        return None
    
    search_str_short = search_str
    if search_str:
        search_str = os.path.join(directory,search_str)
    else:
        search_str = os.path.join(directory ,"*")

    filenames = glob.glob(search_str)
    filenames.sort()
 
    if start and end:
        date_range = np.arange(start,end,dtype="datetime64[D]").astype(str)
        date_range = [date.replace('-',date_separator) for date in date_range]
        print('Finding files in range: {0} - {1} in directory {2} using search string {3}'.format(start,end,directory,search_str_short))
        
        files = []
        for filename in filenames:
            try:
                if date_separator:
                    d_sep = "[{0}]".format(date_separator)
                    re_str = "\d{4}"+d_sep+"\d{2}"+d_sep+"\d{2}"  # Creating a regular expression to find 8 digits separated by date_separator e.g. 2012-01-01
                else:
                    re_str = "\d{8}" # Creating a regular expression to find an 8 digit number (should be the date) e.g. 20120101
                d = re.search(re_str,filename) 
                d = d.group() # Extract value from regular expression compiler
            except AttributeError:
                pass
            else:
                if d in date_range:
                    files.append(filename)
    else:
        print('Finding files in directory {0} with search str {1}'.format(directory,search_str_short))
        files = filenames
    
    return files

def extract_files_dir_split(directory,search_str=None,start=None,end=None,date_separator=''):
    '''
    The extract_files_dir_split function looks for files within sub-directories split by year.
    E.g. if directory="/data/shared/obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/" then files for 2011 would be within
    "/data/shared/obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/2011/"
    
    Args:
        directory (str) :
            Top level directory.
        search_str (str/None, optional) :
            String to use to search directory (using glob). If full filename is not specified this should 
            contain at least one wildcard character ("*").
            If no search_str specified, all files from directory will be returned and filtered by any start 
            and end criteria.
        start (str/None, optional) :
            Start date of the form "YYYY-MM-DD".
        end (str/None, optional) :
            End date of the form "YYYY-MM-DD". Date range is up to but not including this date.
            Either both or neither of start and end should be specified.
        date_separator (str, optional) :
            Date separator string between year, month and day in filename.
            By default this is "", meaning dates of the form "YYYYMMDD" will be searched for.
    
    Returns:
        list :
            Filenames as string s(full path information).
    '''
    
    input_date_separator = '-'
    
    if start and end:
            year_start = start.split(input_date_separator)[0]
            year_end = end.split(input_date_separator)[0]
    elif not start and not end:
        all_dir = os.listdir(directory)
        all_dir_years = [d for d in all_dir if len(d) == 4 and re.match("\d{4}",d)]
        year_start = min(all_dir_years)
        year_end = max(all_dir_years)
        print "No start and end date specified, so processing all files from all date labelled sub-directories within input directory {} for year range: {}-{}".format(directory,year_start,year_end)
    year_range = range(int(year_start),int(year_end)+1)
    
    # Extracting files from directory based on start and end dates (if present). Otherwise extract all files.    
    files = []
    for year in year_range:
        directory_year = os.path.join(directory,str(year))
        if str(year) == year_start and start:
            start_1 = start
        elif year > int(year_start) or not start:
            start_1 = '{}-01-01'.format(year)
        
        if str(year) == year_end and end:
            end_1 = end
        elif int(year) < int(year_end) or not end:
            end_1 = '{}-01-01'.format(year+1)
        if start_1 != end_1:
            files.extend(extract_files(directory_year,search_str,start=start_1,end=end_1,date_separator=''))
    
    return files
    

def ds_check_unique(ds1,ds2,dim_apply="time"):
    '''
    The ds_check_unique function checks for any repeats between two datasets on the dim_apply axis.
    Any repeats are removed from the first dataset (ds1) and returned.
    
    Args:
        ds1 (xarray.Dataset) : 
            First Dataset with dim_apply dimension
        ds2 (xarray.Dataset) : 
            Second Dataset with dim_apply dimension to compare to.
        dim_apply (str, optional) : 
            Dimension to check for repeats (str). Default = "time"
    
    Returns:
        xarray.Dataset : 
            ds1 with any matching values removed.    
    '''
    
    t1 = ds1[dim_apply]
    t2 = ds2[dim_apply]
    
    filt = ~np.in1d(t1,t2) # Compare values in ds2 to ds1 and check for repeats. Want to keep values that aren't repeats.
    ds = apply_filter(ds1,filter_array=filt,dim_apply=dim_apply)
    
    return ds
    

def name_pressure_file(filename,name='surface_pressure',column_names=["latitude","longitude","time"],
                       set_columns=["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"]):
    '''
    The name_pressure_file function extracts the pressure values from one file of a special NAME run to find surface pressure 
    values across a given domain.
    
    NOTE: This function can also be used to extract Topography information from the same NAME run
    To do this set the filename to a topography file (e.g. "Topog_C1_T1_201201010000.txt") and the name to e.g. "topography".
    The set_columns should be the same between pressure and topography files.
    
    Output from:
    - BackRuns_SurfaceField.scr (with SurfaceFieldOutputRequests.txt)
    
    Args:
        filename (str) : 
            Pressure file produced by the appropriate BackRuns script including path information.
            Current script calls this: "Pressure_C1_'YYYYMMDD'_1d.txt" e.g. "Pressure_C1_20120101_15d.txt"
        name (str) : 
            Name to use when creating new pressure (or topography) data variable within output Dataset.
            Default = "surface_pressure".
        column_names (list) : 
            Names for co-ordinate axes within dataset. Should be a 3-item list for latitude, longitude and time.
            Default = ["latitude","longitude","time"]
        set_columns (list) : 
            Column names within input file which contain longitude, latitude and height (in that order).
            Function assumes these columns will be the first in the file and all other populated columns
            contain the surface pressure (or topography) values.
            Default = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"] based on current BackRuns script.
    Returns:
        xarray.Dataset : 
            Dataset containing extracted information with dimensions of latitude,longitude,time
    '''
    
    ## Extract data from input file
    lineskip = 35 # 35 lines to skip for a standard NAME file before reaching final header row
    set_columns = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"]
    num_set_columns = len(set_columns)
    df = pd.read_csv(filename,skipinitialspace=True,skiprows=lineskip)
        
    ## Remove empty column from data, assume this has been dropeped from the end and compare to original length    
    initial_len = len(df.columns)
    df = df.dropna(axis=1) # Normally contains empty column at the end
    if len(df.columns) != initial_len:
        drop_last = True
    else:
        drop_last = False
    
    ## Extract longitdue and latitude from dataframe and then drop all label columns
    lon = df[set_columns[0]].values
    lat = df[set_columns[1]].values
    df = df.drop(set_columns,axis=1)

    ## Extract pressure values from dataframe and rearrange from row-based to column-based
    pressures = df[:].values
    pressures = np.swapaxes(pressures,1,0) # Swap from row-based to column-based
    
    ## Extract dates lines from file (skip all rows expect this line).
    date_row = 33 # Date row should be on the 33rd (34th?) line
    dates_str = pd.read_csv(filename,skipinitialspace=True,skiprows=date_row-1,nrows=1)
    if drop_last:
        dates_str = dates_str[:].values[0][num_set_columns:-1]
    else:
        dates_str = dates_str[:].values[0][num_set_columns:]
    
    ## Convert dates to datetime object via dt.datetime module; expect date format of: "02/01/2012 00:00 UTC"
    dates = [np.datetime64(dt.datetime.strptime(date.rstrip("UTC"),'%d/%m/%Y %H:%M ')) for date in dates_str]

#   ### Previous method for converting to xarray Dataset but was slow    
#    ## Define axes and elements to create new Multi-Index dataframe - can be converted simply to an xarray dataset
#    num_elements = len(df.index)*len(df.columns)
#    
#    lon_index = np.zeros(num_elements)
#    lat_index = np.zeros(num_elements)
#    dates_index = np.array([np.datetime64("1970-01-01 00:00")]*num_elements,dtype=np.datetime64) # Enter dummy datetime value to initialise array.
#
#    num_rows = len(df.index) # Number of rows for each surface pressure column
#    for i,column in enumerate(pressures):
#        dates_index[i*num_rows:(i+1)*num_rows] = np.array([dates[i]]*num_rows)
#        lon_index[i*num_rows:(i+1)*num_rows] = lon
#        lat_index[i*num_rows:(i+1)*num_rows] = lat
#    
#    pressures_stack = np.ravel(pressures) # Stack all pressure values into one column
#    
#    df_split = pd.DataFrame(pressures_stack,index=[dates_index,lon_index,lat_index]) # Create MultiIndex DataFrame
#
#    ## Create xarray dataset and convert into expected format
#    ds = df_split.to_xarray() # Cast to xarray.dataset
#    ds = ds.rename({"level_0":column_names[2],"level_1":column_names[1],"level_2":column_names[0],0:name})
#    ds = ds.transpose(*column_names) # Rearrange dimensions to match typical footprint order
#   ######
    
    # Lon and Lat columns are presented as a grid of numbers, with one lon value for a set of lat values, so should be safe to do this.
    lon_index = np.unique(lon)
    lat_index = np.unique(lat)
    pressures = np.reshape(pressures,(pressures.shape[0],len(lon_index),len(lat_index)))

    p_col_names = [column_names[2],column_names[1],column_names[0]]
    coords = OrderedDict([(p_col_names[0],np.array(dates)),(p_col_names[1],lon_index),(p_col_names[2],lat_index)])
    
    ds = xray.Dataset(data_vars={name:(p_col_names,pressures)},coords=coords)
    ds = ds.transpose(*column_names) # Rearrange dimensions to match typical footprint order

    if len(np.where(ds["time"] == "1970-01-01 00:00")[0]) != 0:
        print('WARNING: Some time values have not been propagated correctly. Default entries of "1970-01-01 00:00" entries are still present')
        return None
    
    return ds

def name_pressure(directory,start_date=None,end_date=None,name="surface_pressure",
                  column_names=["latitude","longitude","time"],
                  set_columns=["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"],max_days=31):
    '''
    The name_pressure function extracts the pressure values from a special NAME run to find surface pressure 
    values across a given domain.
    These values will often be spread over multiple files and this function can take several files from a 
    directory and create one dataset.
    
    Args:
        directory (str) : 
            Directory containing the NAME output files for the SurfacePressure run.
        start_date (str, optional) : 
            Start date range of interest for the pressure values.
        end_date (str, optional) : 
            End date range of interest for the pressure values (optional)
        name (str, optional) : 
            Name to use when creating new pressure data variables within output Dataset.
            Default = "surface_pressure"
        column_names (list, optional) : 
            Names for co-ordinate axes within dataset. Should be a 3-item list for latitude, longitude and time.
            Default = ["latitude","longitude","time"]
        set_columns (list, optional) : 
            Column names within input file which contain longitude, latitude and height (in that order).
            Function assumes these columns will be the first in the file and all other populated columns
            contain the surface pressure (or topography) values.
            Default = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"] based on current BackRuns script.
        max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
    
    Returns:
        xarray.Dataset : 
            Surface pressure data covering at least start_date - end_date range (if specified)    
    '''
    
    #search_str = "Pressure_C1_*_{0}d.txt".format(max_days)
    search_str = "Pressure_C1_*.txt"
    
    if start_date and end_date:
        print('Setting tolerance of {0} days back when searching for NAME pressure files.'.format(max_days))
        start_date = (np.datetime64(start_date) - np.timedelta64(max_days-1,'D')).astype(str)
        end_date = (np.datetime64(end_date) + np.timedelta64(1,'D')).astype(str)
    
    files = extract_files(directory,search_str=search_str,start=start_date,end=end_date)
    
    if len(files) == 0:
        raise Exception('No NAME pressure files found within directory: {0} for date range {1} - {2} (using search_str {3})'.format(directory,start_date,end_date,search_str))
        #return None
    
    for i,filename in enumerate(files):
        if i == 0:
            name_ds = name_pressure_file(filename,name=name,column_names=column_names,set_columns=set_columns)
        else:
            ds = name_pressure_file(filename,name=name,column_names=column_names,set_columns=set_columns)
            ds = ds_check_unique(ds,name_ds,dim_apply="time")
            name_ds = name_ds.merge(ds)
    
    return name_ds   

def name_pressure_match(ds,pressure_domain,columns=["latitude","longitude","time"],p_column="surface_pressure",
                        pressure_base_dir=name_pressure_directory,pressure_convert=1/100.,
                        max_days=31,day_template=True):
    '''
    The name_pressure_match function finds the corresponding NAME surface pressure values based on a set of
    lat,lon and time values extracted from the input dataset.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing pressure values to match the NAME surface pressure values to. Within the 
            dataset the pressure values should be related to latitude, longitude and time values.
        pressure_domain (str) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Check $DATA_PATH/NAME/surface_pressure folder to see which domains currently exist.
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude, longitude and time values.
            Default = ["latitude","longitude","time"]
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_convert (float, optional) : 
            If pressure values extracted from NAME are not in the required units, pressure_convert
            should be set to the scaling factor to convert these units.
            By default, we assume we want to convert to hPa from NAME input in Pa (pressure_covert=1/100.)
        max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Note: Assumes NAME run for e.g. 2012-01-01 ends with datetime 2012-02-01 00:00:00.
            Default = True.
        
    Returns:
        numpy.array : 
            Array of pressure values. One for each lat,lon,time value.    
    '''
    
    #if len(lat) != len(lon) or len(lat) != len(time):
    #    print 'To extract matched name pressure values the same number of latitude, longitude and time values must be provided.'
    #    print 'Lat: {0}, Lon: {1}, Time: {2}'.format(len(lat),len(lon),len(time))
    #    return None
    
    lat = ds[columns[0]].values
    lon = ds[columns[1]].values
    time = ds[columns[2]].values
    
    # Can't applu a list tolerance with current version of xarray. May be able to add with newer versions.
    #lat_tolerance = 5.0
    #lon_tolerance = 5.0
    #time_tolerance = np.timedelta64(max_days,'D')
    #tolerance = [lat_tolerance,lon_tolerance,time_tolerance]
    
    start_date = min(time).astype('M8[D]').astype(str)
    end_date = max(time).astype('M8[D]').astype(str)
    
    pressure_dir = os.path.join(pressure_base_dir,pressure_domain)
    
    print('Extracting pressure_NAME from name_pressure function')
    pressure_NAME = name_pressure(pressure_dir,start_date=start_date,end_date=end_date,name=p_column,column_names=columns,max_days=max_days)
    #matched_pressure_NAME = pressure_NAME.sel(method="nearest",**{columns[0]:lat,columns[1]:lon,columns[2]:time}) # Returns grid of values at timexlatxlon
    #matched_pressure_NAME = matched_pressure_NAME[p_column].values
    #matched_pressure_NAME = np.diagonal(np.diagonal(matched_pressure_NAME)) # Only want values in this grid on the diagonal for our purposes (i.e. time1,lat1,lon1; time2,lat2,lon2, not time1,lat2,lon1 etc.).
    matched_pressure_NAME = np.zeros(len(lat))
    if day_template:
        end_datetime = np.max(pressure_NAME[columns[2]].values) # Extract maximum datetime value
    for i,la,lo,t in zip(range(len(lat)),lat,lon,time):
        if not day_template or t <= end_datetime:
            match = pressure_NAME.sel(method="nearest",**{columns[0]:la,columns[1]:lo,columns[2]:t})
        else:
            # Assumes offset is +1 day because we expect NAME run for e.g. 2012-01-01 to end with datetime 2012-02-01 00:00:00.
            max_offset_days = np.timedelta64((t - end_datetime),'D').astype(int)+1
            pressure_NAME_offset = pressure_NAME.copy(deep=True)
            pressure_NAME_offset[columns[2]] += np.timedelta64(max_offset_days,'D') # Create copy of dataset with date offset to date we're looking for.
            #day_tolerance = tolerance[:1] + [np.timedelta64(1,'D')]
            match = pressure_NAME_offset.sel(method="nearest",**{columns[0]:la,columns[1]:lo,columns[2]:t})
        matched_pressure_NAME[i] = match[p_column].values
    
    matched_pressure_NAME = matched_pressure_NAME*pressure_convert
    
    return matched_pressure_NAME

def name_topog_match(ds,pressure_domain,columns=["latitude","longitude"],name="topography",
                     topog_base_dir=name_pressure_directory):
    '''
    The name_topog_match function finds the corresponding topography values based on a set of lat,lon coordinates.
    Note: expect lat,lon arrays are the same length and each set of axis corresponds to a pressure value at 
    that lat,lon coord.
    Note: Uses first topography file found within relevant folder and does not check the date.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing latitude and longitude values values to match the NAME topography values to.
        pressure_domain (str) :
            Domain over which topography (and surface pressure) values have been extracted (can be 
            distinct from domain if pressure_domain contains area of domain).
            Check $DATA_PATH/NAME/surface_pressure folder to see which domains currently exist.
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude and longitude values.
            Default = ["latitude","longitude"]
        name (str, optional) :
            Name to use when creating new pressure data variables within output Dataset.
            Default = "topography"
        topog_base_dir (str, optional) :
            Base directory containing the NAME output files for the SurfacePressure run also
            containing topography information.
            Filename is assumed to be of the form "Topog_C1_T1_*.txt"
            See name_pressure_file() function for more details.
    
    Returns:
        numpy.array : 
            Array of topography values. One for each lat,lon pair.
    '''
    
    lat = ds[columns[0]]
    lon = ds[columns[1]]

    if len(columns) != 3:
        column_names = [columns[0],columns[1],"time"]
    else:
        column_names = columns

    topog_dir = os.path.join(topog_base_dir,pressure_domain)
    topog_search_str = "Topog_C1_T1_*.txt"
    topog_file = glob.glob(os.path.join(topog_dir,topog_search_str))[0]

    topog_NAME = name_pressure_file(topog_file,column_names=column_names,name=name)
    matched_topog_NAME = topog_NAME.sel(method="nearest",**{columns[0]:lat,columns[1]:lon}) # Returns grid of values at latxlon
    matched_topog_NAME = matched_topog_NAME[name].values
    matched_topog_NAME = np.diagonal(matched_topog_NAME) # Only want values in this grid on the diagonal for our purposes (i.e. lat1,lon1; lat2,lon2, not lat2,lon1 etc.).
    matched_topog_NAME = matched_topog_NAME[0] # Extract single time point.
    
    return matched_topog_NAME
    
def name_pressure_filter(ds,filters,pressure_NAME=None,columns=["latitude","longitude","time","pressure_levels"],
                         cutoff=5.0,layer_range=[50.,500.],
                         pressure_base_dir=name_pressure_directory,pressure_domain=None,max_days=31,
                         day_template=True,pressure_convert=1/100.):
    '''
    The name_pressure_filter function removes any data points from a dataset which does not match the criteria set out
    by the filters parameter (and associated values).
    Filter options include:
        - cutoff based on a percentage difference between NAME and first pressure level
        - check NAME surface pressure is greater than (i.e. below in height) the second pressure level
        - check rough extent of the surface layer based on NAME pressure and second pressure level is in a certain range
    
    Note:
        Surface pressure level is (currently) defined as p0, where p0 is the first pressure level.
        This is correct for the GOSAT Level 2 Product based on the CCI GHG Documentation for level-based extraction: 
            http://www.esa-ghg-cci.org/index.php?q=webfm_send/160
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing pressure values with multiple pressure_levels defined for each data point.
        filters (list) : 
            Which filters to apply based on NAME surface pressure
             - "cutoff": remove all points where surface pressure is outside a cutoff value compared to NAME
             - "level_order": remove all points where NAME surface pressure is less than pressure level 2
             - "dpressure_range": remove all points where NAME surface layer is outside a range of sizes.
        pressure_NAME (np.array, optional) : 
            If pressure from NAME run has already been extracted, this can be specified explicitly to save 
            computing time.
            If not specified, columns from ds and pressure_dir will be used to extract matching pressure values.
        columns (list, optional) : 
            Names of the latitude, longitude, time and pressure variables to extract from input Dataset.
            This should be a 4-item list. Default = ["latitude","longitude","time","pressure_levels"]
            Note: Latitude, Longitude and Time are used to match to NAME pressure values.
        cutoff (float, optional) : 
            Only used when "cutoff" is within filters. Percentage cutoff to apply from comparison between 
            input pressure data and NAME pressure. Default = 5.0
        layer_range (list, optional) : 
            Only used when "dpressure_range" is within filters. Range in metres the surface layer should have 
            (will be converted to pressure units using barometric equation). (two-item list). 
            Default = [50.,500.]
        pressure_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_domain (str/None,optional) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            * Must be specified if pressure_NAME has not been specified *
            Check $DATA_PATH/NAME/surface_pressure folder to see which domains currently exist.
        max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        pressure_convert (float, optional) : 
            If pressure values extracted from NAME are not in the required units, pressure_convert
            should be set to the scaling factor to convert these units.
            By default, we assume we want to convert to hPa from NAME input in Pa (pressure_covert=1/100.)
        
    Returns:
        xarray.Dataset : 
            Filtered Dataset with data points with pressure levels too different from the NAME values (based on 
            the filters specified) removed. 
    '''
    
    dim_apply = "time"
    
    if pressure_NAME is None:
        pressure_NAME = name_pressure_match(ds,columns=columns[:-1],pressure_domain=pressure_domain,
                                            pressure_base_dir=pressure_base_dir,
                                            max_days=max_days,day_template=day_template,
                                            pressure_convert=pressure_convert)
    
    dpressure_w_NAME = define_pressure_levels(ds,p_column=columns[-1],use_name_pressure=True,columns=columns,
                                              pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                                              max_days=max_days,pressure_NAME=pressure_NAME)[1]
    
    #pressure_levels = ds[columns[3]].values
    pressure_levels = ds[columns[3]]
    
    #filt = np.arange(0,pressure_levels.shape[0],1)
    attr = "Original data has been filtered to exclude any data points based on the following conditions: "
    ds_new = ds.copy(deep=True)
    
    
    if "cutoff" in filters:
        model_diff = pressure_NAME - pressure_levels[:,0]
        percent_diff = np.abs(model_diff*100./pressure_NAME)
        ds_new = ds_new.where(percent_diff <= cutoff, drop=True)
        #filt_1 = np.where(percent_diff <= cutoff)[0]
        #filt = filt_1
        attr += "the NAME surface pressure differs by more than {0}% from the first satellite level, ".format(cutoff)
    
    if "level_order"in filters:
        model_diff_level_2 = pressure_NAME - pressure_levels[:,1]
        ds_new = ds_new.where(model_diff_level_2 > 0.0,drop=True)
        #filt_2 = np.where(model_diff_level_2 > 0.0)[0]
        #filt = np.intersect1d(filt,filt_2)
        attr += "the NAME surface pressure is less than the second extracted satellite level, "
    
    if "dpressure_range" in filters:
        layer_size = dpressure_w_NAME[:,0]
        p_at_ground = pressure_at_height(0)
        pressure_range = [(p_at_ground - pressure_at_height(layer_range[0]))*pressure_convert,(p_at_ground - pressure_at_height(layer_range[1]))*pressure_convert]
        ds_new = ds_new.where((layer_size > pressure_range[0]) & (layer_size <= pressure_range[1]),drop=True)
        filt_3 = np.where((layer_size > pressure_range[0]) & (layer_size <= pressure_range[1]))[0]
        #filt = np.intersect1d(filt,filt_3)
        ds_new = apply_filter(ds,filter_array=filt_3,dim_apply=dim_apply)
        attr += "the inferred surface pressure layer based on NAME data would be outside the range {0} - {1}m (taken as {2}-{3}hPa).".format(layer_range[0],layer_range[1],pressure_range[0],pressure_range[1])
    
    #ds_new = apply_filter(ds,filter_array=filt,dim_apply=dim_apply)
    
    # Add attribute describing modification made to original data
    mod_attr = "name_pressure_filter"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = attr

    return ds_new

def midpoint_bounds(pressure_levels,delta_pressure=None,pressure_NAME=None,set_as_edges=False,
                    calc_boundaries=True):
    '''
    The midpoint_bounds function uses midpoint between pressure points to define the boundaries. The pressure levels are
    then scaled to be at the centre of that boundary.
    
    Note:
        The assumption is made that first level represents the surface pressure (unless pressure_NAME is specified, where NAME 
        surface pressure value replaces p_0).
    
    Args:
        pressure_levels (numpy.array) : 
            Array containing pressure levels values for each data point (2D array)
        delta_pressure (numpy.array, optional) : 
            Difference between the pressure levels.
        pressure_NAME (numpy.array, optional) : 
            NAME surface pressure values which correspond to pressure_levels latitude, longitude and time coordinates.
            (1D numpy.array). These values can be extracted using the name_pressure_match() function.
            If this is specified, these NAME values will be used as the surface pressure values rather than the 
            first level for each data point within the pressure_levels array.
        set_as_edges (bool, optional) : 
            Whether to treat the pressure levels as points at the edge of the layers rather than within the layers.
            Default = False.
        calc_boundaries (bool, optional) : 
            Whether to apply the necessary special conditions to calculate the bottom and top layers. 
            Note: if this value is set to False, the pressure_levels and dpressure output arrays will not be the 
            same length as the input pressure_levels array.
            Default = True.
    
    Returns:
        (np.array,np.array) : 
            new_pressure_levels,dpressure
    '''

    min_pressure = 0.0
    
    if delta_pressure is not None:
        delta_pressure = pressure_levels[:,:-1] - pressure_levels[:,1:]

    if pressure_NAME is not None:
        use_name_pressure = True
    else:
        use_name_pressure = False
    
    if calc_boundaries:
        dpressure = np.zeros(pressure_levels.shape)
        new_pressure_levels = np.zeros(pressure_levels.shape)
        
        if set_as_edges:
            if use_name_pressure:
                delta_pressure[:,0] = pressure_NAME - pressure_levels[:,1]
                pressure_levels[:,0] = pressure_NAME
            
            dpressure[:,0] = delta_pressure[:,0]
            new_pressure_levels[:,0] = pressure_levels[:,0] - dpressure[:,0]/2. # Define new pressure level from p0 (surface)
            
            dpressure[:,1:-1] = delta_pressure[:,1:]

            dpressure[:,-1] = delta_pressure[:,-1] # Set last dpressure as matching previous
            bound_below_min = np.where(pressure_levels[:,-1] - dpressure[:,-1] < min_pressure)[0]
            dpressure[bound_below_min,-1] = pressure_levels[:,-1] - min_pressure # If pressure would be < min_pressure, set dpressure as distance from min_pressure (e.g. 0.0)

            new_pressure_levels[:,1:] = pressure_levels[:,1:] - dpressure[:,1:]/2.
        else:
            if use_name_pressure:
                delta_pressure[:,0] = pressure_NAME - pressure_levels[:,1]
                pressure_levels[:,0] = pressure_NAME
            
            dpressure[:,0] = delta_pressure[:,0]/2.
            new_pressure_levels[:,0] = pressure_levels[:,0] - dpressure[:,0]/2.
            
            dpressure[:,1:-1] = delta_pressure[:,:-1]/2. + delta_pressure[:,1:]/2.
            
            dpressure[:,-1] = delta_pressure[:,-1] # Set last dpressure as matching previous
            bound_below_min = np.where(pressure_levels[:,-1] - dpressure[:,-1]/2. < min_pressure)[0]
            dpressure[bound_below_min,-1] = pressure_levels[:,-1] + delta_pressure[:,-1]/2. - min_pressure # If pressure would be < min pressure (e.g. 0), set dpressure as distance from min_pressure
            
            new_pressure_levels[:,1:] = (pressure_levels[:,1:] + delta_pressure/2.) - dpressure[:,1:]/2.
    else:
        dpressure = delta_pressure[:,:-1]/2. + delta_pressure[:,1:]/2.
        new_pressure_levels = pressure_levels[:,1:-1] + delta_pressure[:,:-1]/2. - dpressure/2.

    return new_pressure_levels,dpressure
    

def define_pressure_levels(ds,pressure_domain=None,p_column="pressure_levels",set_as_edges=False,
                           include_bounds=False,use_name_pressure=False,pressure_NAME=None,
                           columns=["latitude","longitude","time"],pressure_base_dir=name_pressure_directory,
                           max_days=31,day_template=True):
    '''
    The define_pressure_levels function scales the pressure_levels and dpressure values for input into NAME.
    Note: values to be output to obs/ files *should NOT* be run through this function.
    
    The pressure_levels provided with the GOSAT data (in particular) define instantaneous pressure points and not the boundaries or 
    midpoints of each of the layers.
    The first pressure level (pressure_0) is treated as the surface pressure value.
 
    To create layers from these pressure levels the boundaries are defined as the midpoint between each level. If the levels are 
    unevenly spaced the pressure levels themselves are then re-defined to be at the midpoint within each layer 
    (or as the edge if set_as_edges=True). The top and bottom layer are treated as follows:
         - Bottom layer - Either whole or 1/2 difference between NAME or satellite surface pressure is used to define layer extent and
                          boundaries depending on if set_as_edges is True or False.
         - Top layer - set layer extent and as matching previous level value OR set layer extent and boundaries based on a minumum 
                       pressure (set as 0 hPa), whichever dpressure is smaller.
     
    Args:
        ds (xarray.Dataset) : 
            Dataset containg multiple data points with a set of pressure levels for each.
        pressure_domain (str/None, optional) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Must be specified if pressure_NAME is not included explicitly.
            Check $DATA_PATH/NAME/surface_pressure folder to see which domains currently exist.
        p_column (str, optional) : 
            Name for the pressure_levels data (str). Default = "pressure_levels"
        set_as_edges (bool, optional) : 
            Whether to treat the pressure levels as points at the edge of the layers rather than within the layers .
            Default = False.
            Note: level 0 will ALWAYS be treated as being the surface pressure (i.e. at the boundary) regardless 
            of this parameter.
        include_bounds (bool, optional)   : 
            Whether to output the boundary values for the layers as well as the pressure_levels and dpressure.
            Default = False.
        use_name_pressure (bool, optional) : 
            If to use the NAME surface pressure rather than the GOSAT value included within the data. 
            Default = False.
        pressure_NAME (numpy.array, optional) : 
            This and the following parameters are only used when use_name_pressure=True.
            If pressure from NAME run has already been extracted, this can be specified with this parameter to 
            save computing time.
            If not specified, columns from ds and pressure_base_dir and pressure_domain will be used to extract matching 
            pressure values.
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude, longitude and time values 
            (3 item list). Default = ["latitude","longitude","time"]
        pressure_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        max_days (int, optional) : 
            Maximum number of days from time within ds to search for the relevant pressure data. 
            Default = 31 (days).
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
                            
    
    Returns:
        if include_bounds == True:
            (np.array,np.array,np.array):
                pressure_levels,dpressure,pressure_bounds
        else:
            (np.array,np.array): 
                pressure_levels,dpressure

    '''

    pressure_levels = ds[p_column].copy(deep=True).values
    delta_pressure = pressure_levels[:,:-1] - pressure_levels[:,1:]
    
    if use_name_pressure:
        if pressure_NAME is None:
            if pressure_domain is not None:
                pressure_NAME = name_pressure_match(ds,pressure_domain=pressure_domain,
                                                pressure_base_dir=pressure_base_dir,columns=columns,
                                                day_template=day_template,max_days=max_days)
            else:
                raise Exception("Pressure_domain must be specified if pressure values need to be \
                                extracted to use define_pressure_levels function.")
    else:
        pressure_NAME = None
    
    # Use midpoint between points to define the boundary
    # Scale the pressure points to be at the centre of that boundary
    pressure_levels,dpressure = midpoint_bounds(pressure_levels,delta_pressure,pressure_NAME,set_as_edges=set_as_edges)

    if include_bounds:
        pressure_bound = np.zeros((len(pressure_levels),len(pressure_levels[0])+1))
        pressure_bound[:,:-1] = pressure_levels + np.abs(dpressure/2.)
        pressure_bound[:,-1] = pressure_levels[:,-1] - np.abs(dpressure[:,-1]/2.)
        
        return pressure_levels,dpressure,pressure_bound
    else:
        return pressure_levels,dpressure    

def ds_check_internal_unique(ds,axis="time"):
    '''
    The ds_check_internal_unique function checks whether an xarray.Dataset has repeat values on the "time" axis.
    Any time values which are repeated are modified by a small time increment to make the values unique. 
    
    Args:
        ds (xarray.Dataset):
            Dataset with time axis.
        axis (str, optional):
            Name of time axis.
            Default = "time"
    
    Returns:
        xarray.Dataset:
            If no repeats are present:
                Original dataset is returned
            If repeats are present:
                The repeats of any time values are modified by applying a small random increment so the time values
                are not identical.
    '''
    
    if len(ds[axis].values) > len(np.unique(ds[axis].values)):
        axis_copy = ds[axis].values
        unique,inverse,counts = np.unique(axis_copy,return_inverse=True,return_counts=True)
        repeat_index_0 = np.where(counts > 1)[0] # Find first indices of all elements with repeats
        repeat_all_index = [np.where(inverse == i)[0] for i in repeat_index_0] # Extract all indices for these repeats
        for indices in repeat_all_index:
            for i in indices[1:]:
                # Add very small random value to repeat of a time value to avoid two times being exactly the same
                axis_copy[i] += np.timedelta64(random.randrange(-1000,1000,1),'us')
        ds[axis].values = axis_copy
        
        # Add attribute describing modification made to original data
        mod_attr = "repeat_time_modified"
        if mod_attr not in ds.attrs.keys():
            ds = add_history_attr(ds,mod_attr)
            ds.attrs[mod_attr] = "Repeated times were found at original indices {} (may not be the same if data has been binned). Small random increments were added to any repeat values to allow xarry to distinguish them.".format(repeat_all_index)
        else:
            ds.attrs[mod_attr] += " Also found at indices {}.".format(repeat_all_index)
    
    return ds
    
def gosat_output_filename(output_directory,network,instrument,date,species,inlet=None,num=None):
    '''
    The gosat_output_filename function creates an output filename of the correct format for gosat based on the 
    inputs.
    
    Filenames are of the form: "instrument"_gosat_"date (reformated)"-"num"_"species"-"inlet".nc
    e.g. /shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/gosat-fts_gosat_20120920-09_ch4-column.nc
    
    Args:
        output_directory (str) : 
            Top level for output directory e.g. "/shared_data/air/shared/obs/"
        network (str) : 
            Which network is being considered e.g. "GOSAT/GOSAT-INDIA"
            This will be used to create the output path e.g. "/shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/"
        instrument (str) : 
            Instrument on satellite being used. e.g. "gosat-fts"
        date (str) : 
            Date the measurements are relevant to e.g. "2010-01-01"
        species (str) : 
            Species being considered e.g. "ch4". Should be defined within "acrg_species_info.json" file
        inlet (str/None, optional) : 
            Additional information to add to ch4 measurement e.g. column
        num (str/None, optional) : 
            Number to append to the filename if gosat data is being written out as multiple files over one day.
            
     Returns:
         str: 
             GOSAT observation filename with path information    
    '''
    #Example /shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/gosat-fts_gosat_20120920-09_ch4-column.nc
    
    output_directory = os.path.join(output_directory,network)
    
    satellite = 'gosat'
    
    date = date.replace('-','') # Turn date from e.g. 2012-09-20 to 20120920
    
    if num:
        date = '-'.join([date,num]) # Create composite string of datetime and num (if present) e.g. "20120920-01"
    if inlet:
        species = '-'.join([species,inlet]) # Create composite string of species and inlet (if present) e.g. "ch4-column"
  
    filename = '_'.join([instrument,satellite,date,species]) # Create filename string joined by "_" e.g. "gosat-fts_gosat_20120920-09_ch4-column"
    filename += '.nc' # Add file extension
    
    filename = os.path.join(output_directory,filename)
       
#   From cf_data_process.output_filename   
#        return join(output_directory,
#                network + "/" + \
#                network + "-" + \
#                instrument + "_" + \
#                site + "_" + \
#                year + "0101_" + \
#                species + "-" + \
#                inlet + ".nc")
#    
    return filename
  
def gosat_split_output(ds,index,mapping=None,data_vars=[],split_dim="time",ident=None,ident_sep=','):
    '''
    The gosat_split_output function creates a new dataset based on input index values along the 
    split dimension (e.g. time).
    This is used to separate data points into subsets (usually to write out to file).
    
    If an ident column is specified the values within this column are also split out into a new "id" dimension.
    
    Note: units for each variable are set by the units() function (and so to include a unit must be
    defined there).
    
    Args:
        ds (xarray.Dataset) : 
            Dataset with consistent dimension to split along.
        index (int/list) : 
            Index value or values to use for the split.
        mapping (dict/OrderedDict, optional) : 
            Mapping between the data variable names within ds and the output data variable names. 
            If no mapping is specified this will use the same data variable names in input dataset 
            for the output dataset.
        data_vars (list, optional) : 
            Data variables within ds to include in output dataset. (list)
            If both mapping and data_vars are specified, mapping will supercede data_vars.
            If no mapping is included and data_vars if left blank ([] by default) all data variables will 
            be copied.
        split_dim (str, optional) : 
            Dimension to split along (i.e. that index is relevant to). Default = 'time' (str)
        ident (str, optional) : 
            Data variable which includes identifier values for all data points within a bin e.g. "exposure_id".
            If ident is set then this data variable will be split out to contain an extra dimension ("id")
            The ident data variable is expected to contain some identifier for each data point separated 
            by the ident_sep value (e.g. commas)
            e.g. ['2010123106380100271013,2010123106380100271014,2010123106380100271015']
        ident_sep (str, optional) : 
            Separator value to use for identifier column.
            Default is a comma (','). Only used if ident is set.
    
    Returns:
        xarray.Dataset : 
            Dataset containing values specified by indices along split_dim (e.g. time)
            Identifier column will be split into an extra dimension if ident is specified
    '''
    
    if not data_vars and not mapping:
        data_vars = ds.data_vars()
    
    if not mapping:
        mapping = OrderedDict([(name,name) for name in data_vars]) # If no mapping specified, assume data variable names in created ds match input names
   
    if data_vars:
        if ident not in data_vars:
            print 'WARNING: Identifier column {0} is not within input data_vars: {1}. No ident column will be included'.format(ident,data_vars)
    
    if isinstance(index,int):
        indices = [index]
    else:
        indices = index
    
    coords = OrderedDict([])
    data = OrderedDict([])
    
    # Create data and coords to be included within dataset
    for name,new_name in mapping.items():
        
        dims = ds[name].dims # Extract dimensions for variable from current dataset
        data_var = ds[name][indices]
        
        # Format the identifier data variable (e.g. exposure_id) to contain an extra "id" dimension to allow for multiple values
        if name == ident:
            identifiers = [value.split(ident_sep) for value in data_var.values] # Split identifier value by the ident_sep value (e.g. ',')
            identifiers = np.array(list(itertools.izip_longest(*identifiers,fillvalue=np.nan))).T # Create array with consistent dimensions for "id" and fill in any gaps with np.nan values
            
            id_dim_name = "id" # Define new dimension name
            split_dim_dim,id_dim_dim = identifiers.shape # Define dimensionality of new dimension and dimension we're splitting on
            id_dim_coord = np.arange(1,id_dim_dim+1)
            split_dim_coord = ds[split_dim][indices] # Extract associated split dimension coordinate e.g. time coord
            
            id_coords = OrderedDict([(split_dim,split_dim_coord),
                                     (id_dim_name,id_dim_coord)])
            id_dims = OrderedDict([(split_dim,split_dim_dim),
                                   (id_dim_name,id_dim_dim)])
            
            data_var = xray.DataArray(identifiers,coords=id_coords,dims=id_dims) # Reset data_var as new variable with extra dimension
            
            # Add suitable attributes for the new ident data variable
            if ident == "exposure_id":
                data_var.attrs["short_description"] = "Exposure identification number of the sounding for each data point within the bin."
                data_var.attrs["long_name"] = "exposure_id_mult"
            else:
                data_var.attrs["short_description"] = "{} label for each data point within the bin.".format(name)
                data_var.attrs["long_name"] = ident
            
            data_var[id_dim_name].attrs["short_description"] = "Number for each data point which has been combined within the bin."
            data_var[id_dim_name].attrs["long_name"] = "id"
                  
        # Define units for each data variable, if not already specified
        unit = units(new_name)
        if unit and ("units" not in data_var.attrs):
            data_var.attrs["units"] = unit
        data[new_name] = data_var
        
        # Extract relevant coords from input dataset
        for dim in dims:
            if dim not in coords.keys(): 
                if dim == split_dim:
                    coords[dim] = ds[dim][indices] # If dim is split_dim dimension (e.g. time), extract only the relevant values
                else:
                    coords[dim] = ds[dim] # Extract all coords for other dimensions
        
        # Add extra "id" dimension if identifier column is included.
        if name == ident:
            coords[id_dim_name] = id_coords[id_dim_name]
    
    ds_output = xray.Dataset(data,coords=coords) # Create new dataset
    
    # Add units for coordinates
    for dim in ds_output.coords:
        unit = units(dim)
        if unit and ("units" not in ds_output[dim].attrs):
            ds_output[dim].encoding["units"] = unit # For e.g. "time" dimension, "units" is a special label set via an initial encoding and so must be written this way (rather than as an attribute) or this will produce an error.
    ds_output = ds_output.assign_attrs(ds.attrs) # Copy across global attributes from input dataset
    
    
    return ds_output

def write_netcdf(ds,filename,overwrite=False):
    '''
    The write_netcdf function writes an xarray.Dataset to a netCDF file.
    The function will also checks whether the file exists and only overwrite if this option is set to True.
    
    Args:
        ds (xarray.Dataset) :
            Dataset to write to file
        filename (str) :
            Filename (including path information) for the output file
        overwrite (bool) :
            Whether to overwrite an existing file if present.
            Default = False
    
    Returns:
        None
        
        Writes dataset to file (filename)
    '''    
    if not overwrite:
        if os.path.isfile(filename):
            print '{0} already exists. Data has not been written to file because overwrite=False.'.format(filename)
        else:
            print 'Writing to filename:',filename    
            ds.to_netcdf(filename)
    else:
        print 'Writing to filename:',filename    
        ds.to_netcdf(filename,mode="w")

def gosat_output(ds,site,species="ch4",file_per_day=False,output_directory=obs_directory,
                 overwrite=False):
    '''
    The gosat_output function creates output netCDF files each containing data related to one data point per file.
    The data variables within each netCDF file are:
        "ch4_profile_apriori","lat","lon","pressure_levels","pressure_weights","xch4","xch4_averaging_kernel",
        "xch4_uncertainty"
    Units for these values are set by the units() function.
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file as an xarray Dataset.
            Should have had co-ordinates and dimensions assigned with gosat_add_coords() function.
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA (should be defined within acrg_sites_info.json)
        species (str, optional) : 
            Species of interest e.g. "ch4" (should be defined within acrg_species_info.json).
            Default = "ch4"
        file_per_day (bool, optional) : 
            Output all results to one file per day rather than splitting out per time point. Default = False.
        output_directory (str, optional) : 
            Top level directory to write output. Full path will be based on the "network" related to "site" .
            Default = obs_directory (defined at top of file).
        overwrite (bool, optional) : 
            Whether to overwrite any files already present in the output_directory + network folder.
            Default = False.
    
    Returns:
        None
        
        Writes output to multiple .nc files (split on time axis).
    
    '''
    # Define data variables to be written to output
    out_data_vars = ["xch4","xch4_uncertainty","lat","lon","pressure_levels","pressure_weights",
                     "xch4_averaging_kernel","ch4_profile_apriori","exposure_id"]
    
    # Map to input dataset (from GOSAT data)
    data_vars = ["latitude" if item=="lat" else item for item in out_data_vars]
    data_vars = ["longitude" if item=="lon" else item for item in data_vars]
    data_vars = ["pressure_weight" if item=="pressure_weights" else item for item in data_vars]
    
    data_var_mapping = OrderedDict([(name,new_name) for name,new_name in zip(data_vars,out_data_vars)])
    
    #if site == None:
    #    site="global"
    
    # Set dimension on which data will be split into separate files
    split_dim="time"
    
    # Set name of data variable which includes the data point identifiers
    ident = "exposure_id"    
    
    try:
        network = agage.site_info[site]["network"]
    except KeyError:
        network = "unknown"
    instrument = 'gosat-fts'
    inlet = 'column'
    species = species.lower()
    
    #all_dates = ds["time"].values.astype('M8[D]').astype(str) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)
    all_dates =  extract_dates(ds,dim=split_dim)
    dates = np.unique(all_dates)
    
    full_output_directory = os.path.join(output_directory,network)
    if not os.path.isdir(full_output_directory):
        os.makedirs(full_output_directory)
    
    for date in dates:
        wh_date = np.where(all_dates == date)[0] # Find indices for each date
        if file_per_day:
            ds_output = gosat_split_output(ds,index=wh_date,mapping=data_var_mapping,split_dim=split_dim,ident=ident)
            
            # Create filename and write dataset to file
            filename = gosat_output_filename(output_directory,network,instrument,date,species,inlet=inlet)
            ds_output.attrs["id"] = os.path.split(filename)[1]
            write_netcdf(ds_output,filename,overwrite=overwrite)
        else:
            for ID,index in enumerate(wh_date):
                ds_output = gosat_split_output(ds,index=index,mapping=data_var_mapping,split_dim=split_dim,ident=ident)

                # Create filename and write dataset to file
                ID_str = str(ID+1).zfill(3) # Number to add to filename - three digit with leading zeros
                filename = gosat_output_filename(output_directory,network,instrument,date,species,num=ID_str,inlet=inlet)
                ds_output.attrs["id"] = os.path.split(filename)[1]
                write_netcdf(ds_output,filename,overwrite=overwrite)
    

def gosat_output_name(ds,site,max_level=17,use_name_pressure=False,pressure_domain=None,
                      pressure_base_dir=name_pressure_directory,
                      pressure_max_days=31,pressure_day_template=True,name_directory=name_csv_directory,
                      file_per_day=False,overwrite=False):
    '''
    The gosat_output_name function creates a csv file in the correct format for input into NAME-III.
    Note: it is assumed input pressure values are in hPa values and so are converted to Pa for NAME (*100.)
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file as an xarray Dataset.
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within acrg_sites_info.json.
            Note: at the moment this is just used in the filename and not in the folder structure
        max_level (int, optional) : 
            Maximum level to include from input GOSAT data (up to 20).
            At the moment NAME footprints go up to 19km which corresponds to level 17 in GOSAT.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the GOSAT surface pressure value for each data 
            point.
        pressure_domain (str, optional) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Must be included if use_name_pressure=True.
            Check $DATA_PATH/NAME/surface_pressure folder to see which domains currently exist.
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
        pressure_day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        name_directory (str, optional) : 
            Top level directory to write files for NAME. Full path will be based on the "network" related to "site" 
            Default defined at the top of this module.
        file_per_day (bool, optional) : 
            Output all results to one file per day rather than splitting out per time point. 
            Default=False.
        overwrite (bool, optional) : 
            Allow any files already present within the full path to name_directory to be overwritten. 
            Default = False
    
    Returns:
        None
        
        Writes output to multiple .csv files (split on time axis).
    '''
    
    # Note column mapping here is out:in values rather than in:out because some columns in output do not map to input data variables.
    col_mapping = OrderedDict([('ID_Level','lev'),
                               ('Time','time'),
                               ('x','longitude'),
                               ('dx','dlon'),
                               ('y','latitude'),
                               ('dy','dlat'),
                               ('z','pressure_levels'),
                               ('dz',None)])
    
    try:
        ds[col_mapping['dx']]
        ds[col_mapping['dy']]
    except KeyError:
        sounding = 10.5
        dlat = np.zeros(len(ds[col_mapping['y']]))+distance_lat(sounding)
        dlon = np.array([distance_lon(sounding,lat) for lat in ds[col_mapping['y']]])
        axis1 = col_mapping['Time']
        dlat_da = xray.DataArray(dlat,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
        dlon_da = xray.DataArray(dlon,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
        ds = ds.assign(**{col_mapping['dx']:dlon_da,col_mapping['dy']:dlat_da})
    
    columns=col_mapping.keys()
    
    if site is None:
        site = 'global'
    
    split_dim = "time"
    out_index = "ID_Level"
    
    name_pressure_convert = 100. # Input from GOSAT is in hPa and need to convert to Pa for NAME
    
    all_dates = extract_dates(ds,dim=split_dim)
    dates = np.unique(all_dates)
    
    pressure_levels,dpressure = define_pressure_levels(ds,pressure_domain=pressure_domain,
                                                       use_name_pressure=use_name_pressure,
                                                       pressure_base_dir=pressure_base_dir,
                                                       max_days=pressure_max_days,
                                                       day_template=pressure_day_template)
    try:
        network = agage.site_info[site]["network"]
    except KeyError:
        network = "unknown"
    
    name_directory = os.path.join(name_directory,network)
    if not os.path.isdir(name_directory):
        os.makedirs(name_directory)
    
    for date in dates:
        wh_date = np.where(all_dates == date)[0] # Find indices for each different date
        for ID,index in enumerate(wh_date):
            data = OrderedDict()
            ID_str = '{num:03d}'.format(num=ID+1) # Number for each point in a day
            
            for col in columns:
                name = col_mapping[col]
                if name:
                    ds_column = ds[name]
                    if name == "pressure_levels":
                        data[col] = pressure_levels[index][:max_level]*name_pressure_convert
                    elif name == "lev":
                        
                        if file_per_day:
                            # Write out ID as a str combination of the point number and level each with leading zeros
                            data[col] = np.array([ID_str+value.zfill(2) for value in (ds_column[:max_level].values+1).astype(str)])
                        else:
                            # Write out ID (levels) as two digit numbers with a leading zero e.g. 01, 02, ... 21, 22
                            data[col] = np.array([value.zfill(2) for value in (ds_column[:max_level].values+1).astype(str)])
                    else:
                        data[col] = np.array([ds_column[index].values]*max_level) # Populate each level with the same values e.g. lat,lon,dlat,dlon
                else:
                    if col == "dz":
                        data[col] = dpressure[index][:max_level]*name_pressure_convert

            df_output = pd.DataFrame(data,columns=columns)
            df_output.set_index(out_index,inplace=True)
            
            if file_per_day:
                filename = site+"_"+date.replace('-','')+".csv"
                filename = os.path.join(name_directory,filename)
                
                if not overwrite:
                    if ID == 0 and os.path.isfile(filename):
                        print('{0} already exists. Data has not been written to file because overwrite=False.'.format(filename))
                        break
                    elif not os.path.isfile(filename):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename)
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,mode='a',header=False)    
                else:
                    if (not os.path.isfile(filename)) or (os.path.isfile(filename) and ID == 0):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename)
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,mode='a',header=False)
            else:
                filename = "{site}_{date}-{ID}.csv".format(site=site,date=date.replace('-',''),ID=ID_str)
                filename = os.path.join(name_directory,filename)
                if not overwrite:
                    if os.path.isfile(filename):
                        print('{0} already exists. Data has not been written to file because overwrite=False.'.format(filename))
                    else:
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename)
                else:
                        print('Writing to filename: {}'.format(filename))    
                        df_output.to_csv(filename)
                
    
def gosat_process_file(filename,site,species="ch4",lat_bounds=[],lon_bounds=[],domain=None,
                       coord_bin=None,quality_filt=True,bad_pressure_filt=True,name_sp_filt=False,
                       name_filters=[],cutoff=5.,layer_range=[50.,500.],                       
                       mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31.,pressure_day_template=True,
                       write_nc=False,output_directory=obs_directory,
                       write_name=False,name_directory=name_csv_directory,
                       file_per_day=False,overwrite=False,verbose=True):
    '''
    The gosat_process_file function processes input gosat data for one file, applying designated filters, 
    binning and writing output as required.
    Binning is done based on latitude and longitude values in degrees.
    
    Filter criteria that can be applied include:
        - Latitude and longitude range - based on input bounds or specified domain
        - Start and end date range
        - Quality flag (xch4_quality_flag)
        - Bad pressure flag
        - Mode (land or glint)
        - Comparison from NAME surface pressure including a cutoff, comparing surface to the first level and the extent of the layer.
    
    Output can be written to:
        - netCDF files (.nc) based on suitable format for acrg repository. See gosat_output() and gosat_output_filename() functions.
        - text files  (.csv) suitable for input into NAME. See gosat_output_name() function.
    
    Args:
        filename (str) : 
            Filename of GOSAT CH4 Level 2 Data Product file (str)
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within sites_info.json
        species (str, optional) : 
            Species of interest. Should be defined within species_info.json. Default = "ch4"
        lat_bounds (list, optional) : 
            Upper and lower bounds for latitude. (two-item list e.g. [6.0,36.5])
        lon_bounds (list, optional) : 
            Upper and lower bounds for longitude. (two-item list e.g. [70.0,90.5])
        domain (str, optional) :
            If lat_bounds and lon_bounds are not specified a domain can be specified instead.
            This must be a predetermined NAME domain and a footprint file must exist to extract the latitude
            and longitude bounds.
            If both domain and pair of lat_bounds and lon_bounds are specified, the explictly specified
            bounds will take precedence.
            This parameter will also be used to find folder containing surface pressure files (if 
            applicable) if pressure_domain is not specified.
        coord_bin (float/list, optional) : 
            If applying binning, this should be the size of each bin in degrees.
            To specify the same bin for latitude and longitude a single values (float) can be used.
            To specify different bins include a two item list (e.g. [0.234,0.356])
            Set to None if no binning is required.
        quality_filt (bool, optional) : 
            Whether to remove data points using the quality filter flag (xch4_quality_flag) 
            indicating possibly bad data (=1). 
            Default = True.
        bad_pressure_filt (bool, optional) : 
            Whether to remove data points where any pressure value = -9999.99 (seems to be a flag).
            Default = True.
        name_sp_filt (bool, optional) :
            Whether to remove data points based on the comparison between the GOSAT surface pressure value
            and the NAME surface pressure value. 
            Exact conditions are specified by name_filters with associated values specified by cutoff
            and layer_range parameters if applicable.
            Default = False.
        mode (str, optional) : 
            Which retrieval mode (retr_flag) to filter by (if any). Options: 'land' or 'glint'.
            Set to None if no filtering by mode is required.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the input GOSAT surface pressure
            (first pressure_levels value for each time point) when creating the NAME csv file.
            Note:
                 - This is independent of the name_sp_filt and both or either can be applied.
                 - This does not affect the pressure_levels in the netCDF output.
            Default = False.
        name_filters (list) : 
            If name_sp_filt=True, this is a list of which filters will be used.
            Options are: "cutoff", "level_order","dpressure_range". 
            This value must be specified if name_sp_filt=True.
            See name_pressure_filter() function for details of these filters.
        cutoff (float, optional) : 
            When applying the "cutoff" NAME filter, this is the percentage cutoff to apply from comparison 
            between input pressure data and NAME pressure.
            Default = 5.0
        layer_range (list, optional) : 
            When applying the "dpressure_range" NAME filter, this is the range in metres the surface layer 
            should have (will be converted to pressure units using barometric equation). (two-item list). 
            Default = [50.,500.]
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            If pressure_domain (or domain) is specified this will be appended to the pressure_base_dir
            Default = "$DATA_PATH/NAME/surface_pressure/"
        pressure_domain (str, optional) :
            Domain over which surface pressure has been created. Only needs to be specified if this is 
            different (i.e. wider domain) than domain, or domain has not been explictly specified.
            Will be used to find full pressure_dir path to Pressure_*.txt files.
        pressure_max_days (int, optional) : 
            Maximum number of days from time within ds to search for the relevant pressure data. 
            Default = 31 (days).
        pressure_day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        write_nc (bool, optional) : 
            Write output .nc files (one per bin along time axis by default)
            Default = False.
        output_directory (str, optional) : 
            Top level directory to write output. Full path will be based on the "network" related to "site" 
        write_name (bool, optional) : 
            Write output .csv file for input into NAME (one per bin along time axis by default)
        name_directory (str, optional) : 
            Top level directory to write files for NAME. Full path will be based on the "network" related to "site" 
            Default defined at the top of this module.
        file_per_day (bool, optional) : 
            Group together all points into create one output file per day rather than one per 
            data point (bool).
            Default = False.
        overwrite (bool, optional) : 
            Allow any files already present within the output_directory and name_directory to be overwritten.
            Default = False.
        verbose (bool, optional) : 
            Print details of file processing to screen.
            Default = True.
        
    Returns:
        xarray.Dataset: 
            processed GOSAT Dataset
        
        If write_nc is set to True:
            Set of .nc files will be written within path started by output_directory
        If write_name is set to True:
            Set of .csv files will be written to name_directory
    '''
    
    axis="time"
    
    if use_name_pressure or name_sp_filt:
        if domain and not pressure_domain:
            pressure_domain = domain
        elif not domain and not pressure_domain:
            raise Exception("To access NAME surface pressure files, pressure_domain (or domain) must be \
                            specified. Current pressure_domain={}.".format(pressure_domain))
    
    if lat_bounds:
        if len(lat_bounds) == 1 or len(lat_bounds) > 2:
            raise ValueError('Lat bounds must be specified as a two item iterable (e.g. list). Current value: {0}'.format(lat_bounds))
    
    if lon_bounds:
        if len(lon_bounds) == 1 or len(lon_bounds) > 2:
            raise ValueError('Lon bounds must be specified as a two item iterable (e.g. list). Current value: {0}'.format(lon_bounds))
    
    if coord_bin:
        if len(coord_bin) > 2:
            raise ValueError('Coordinate bin should be a one or two item iterable (e.g. list). Current value: {0}'.format(coord_bin))
    
    if verbose:
        print("========================")
        print('\nProcessing file: {0}\n'.format(filename))
    gosat = xray.open_dataset(filename)
    
    gosat = gosat_add_coords(gosat)
    gosat = gosat.sortby(axis)
    gosat = ds_check_internal_unique(gosat,axis) # Check time values are unique and slightly modify if necessary
    
    if quality_filt:
        if verbose:
            print('Applying filter based on quality flag')
        gosat = gosat_quality_filter(gosat)
    if bad_pressure_filt:
        if verbose:
            print('Applying filter based on bad pressure flag')
        gosat = gosat_pressure_filter(gosat) # Removes points with any pressure values of -9999.99
    if mode:
        if mode == 'land' or mode == 'glint':
            if verbose:
                print('Applying filter based on mode: {0}'.format(mode))
            gosat = gosat_mode_filter(gosat,mode=mode)
        else:
            print('WARNING: Did not recognise input for mode filtering: {0}. Should be one of "land" or "glint". No mode filtering applied.'.format(mode))
    if name_sp_filt:
        if not name_filters:
            raise Exception("If name_sp_filt=True, name_filters must be specified.")
        if verbose:
            print('Applying filters based on NAME surface pressure.')
        gosat = name_pressure_filter(gosat,filters=name_filters,pressure_domain=pressure_domain,
                                     cutoff=cutoff,layer_range=layer_range,pressure_base_dir=pressure_base_dir,
                                     max_days=pressure_max_days,day_template=pressure_day_template)
        
    columns=["latitude","longitude"]
    
    if domain and not (lat_bounds and lon_bounds):
        if verbose:
            print "Extracting latitude and longitude bounds from footprints associated with domain: {}".format(domain)
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if not lat_bounds:
        lat_bounds = [min(gosat.latitude.values),max(gosat.latitude.values)]
    if not lon_bounds:
        lon_bounds = [min(gosat.longitude.values),max(gosat.longitude.values)]
    
    if coord_bin:
        if verbose:
            print('Binning data based on {0} degree bins'.format(coord_bin))
            print('Looking at area within latitude and longitude bounds: {0},{1}'.format(lat_bounds,lon_bounds))
        gosat = binned_mean(gosat,lat_bounds,lon_bounds,domain,columns=columns,coord_bin=coord_bin)
        gosat = gosat.sortby(axis)
        gosat = ds_check_internal_unique(gosat,axis=axis) # After binning check time values are unique again
    else:
        if verbose:
            print('Looking at area within latitude and longitude bounds: {0},{1}'.format(lat_bounds,lon_bounds))
        gosat = latlon_filter(gosat,lat_bounds,lon_bounds,columns=columns)

    gosat.attrs["input_file"] = os.path.split(filename)[1]

    if len(gosat[axis].values) > 0:
        if write_nc:
            gosat_output(gosat,site=site,species=species,output_directory=output_directory,
                         file_per_day=file_per_day,overwrite=overwrite)
        
        if write_name:
            gosat_output_name(gosat,site=site,use_name_pressure=use_name_pressure,
                              pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                              pressure_max_days=pressure_max_days,pressure_day_template=pressure_day_template,
                              name_directory=name_directory,file_per_day=file_per_day,overwrite=overwrite)
    else:
        print('No points extracted from {} for specified parameters.'.format(filename))

    return gosat

def gosat_process(site,species="ch4",input_directory=input_directory,start=None,end=None,
                  lat_bounds=[],lon_bounds=[],domain=None,coord_bin=None,
                  quality_filt=True,bad_pressure_filt=True,
                  name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],
                  mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                  pressure_domain=None,pressure_max_days=31,pressure_day_template=True,
                  write_nc=False,output_directory=obs_directory,
                  write_name=False,name_directory=name_csv_directory,file_per_day=False,
                  overwrite=False):
    '''
    The gosat_process function processes GOSAT input files from a directory.
    As part of processing the GOSAT input data there are multiple options including filtering based on
    various criteria and binning. See gosat_process_file() function for more details.
    
    Args:
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within acrg_sites_info.json
        species (str, optional) : 
            Species of interest. Should be defined within acrg_species_info.json
            Default = "ch4"
        input_directory (str, optional) : 
            Top level directory containing GOSAT CH4 Level 2 Data Product files only.
            Assumes sub-directories are labelled by year and will find sub-directory based on start 
            and end dates, if specified.
            If no dates are specified all files within all year labelled sub-directories will be processed.
            See top of file for default.
        start (str, optional) : 
            Start date for range of files to be processed. Should be in format "YYYY-MM-DD"
        end (str, optional) : 
            End date for range of files to be processed (optional but must be specified if start specified).
            Should match format of start.
        lat_bounds (list, optional) : 
            Upper and lower bounds for latitude. (two-item list e.g. [6.0,36.5])
        lon_bounds (list, optional) : 
            Upper and lower bounds for longitude. (two-item list e.g. [70.0,90.5])
        domain (str, optional) :
            If lat_bounds and lon_bounds are not specified a domain can be specified instead.
            This must be a predetermined NAME domain and a footprint file must exist to extract the latitude
            and longitude bounds.
            If both domain and pair of lat_bounds and lon_bounds are specified, the explictly specified
            bounds will take precedence.
            This parameter will also be used to find folder containing surface pressure files (if 
            applicable) if pressure_domain is not specified.
        coord_bin (float/list, optional) : 
            If applying binning, this should be the size of each bin in degrees.
            To specify the same bin for latitude and longitude a single values (float) can be used.
            To specify different bins include a two item list (e.g. [0.234,0.356]).
            Set to None if no binning is required.
        quality_filt (bool, optional) : 
            Whether to remove data points using the quality filter flag (xch4_quality_flag) 
            indicating possibly bad data (=1). 
            Default = True.
        bad_pressure_filt (bool, optional) : 
            Whether to remove data points where any pressure value = -9999.99 (seems to be a flag).
            Default = True.
        name_sp_filt (bool, optional) : 
            Whether to remove data points based on the comparison between the GOSAT surface pressure value
            and the NAME surface pressure value. 
            Exact conditions are specified by name_filters with associated values specified by cutoff
            and layer_range parameters if applicable.
            Default = False.
        mode (str, optional) : 
            Which retrieval mode (retr_flag) to filter by (if any). Options: 'land' or 'glint'.
            Set to None if no filtering by mode is required. 
            Default = None.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the input GOSAT surface pressure
            (first pressure_levels value for each time point) when creating the NAME csv file.
            Note:
                 - This is independent of the name_sp_filt and both or either can be applied.
                 - This does not affect the pressure_levels variable in any netCDF (.nc) output.
        name_filters (list, optional) : 
            If name_sp_filt=True, this is a list of which filters which will be used. Options are: "cutoff",
            "level_order","dpressure_range".
            This value must be specified if name_sp_filt=True.
            See name_pressure_filter() function for details of these filters.
        cutoff (float, optional) : 
            When applying the "cutoff" NAME filter, this is the percentage cutoff to apply from comparison 
            between input pressure data and NAME pressure. 
            Default = 5.0
        layer_range (list, optional) : 
            When applying the "dpressure_range" NAME filter, this is the range in metres the surface layer 
            should have (will be converted to pressure units using barometric equation). (two-item list). 
            Default = [50.,500.]
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            If pressure_domain (or domain) is specified this will be appended to the pressure_base_dir.
            Default = "$DATA_PATH/NAME/surface_pressure/"
        pressure_domain (str, optional) :
            Domain over which surface pressure has been created. Only needs to be specified if this is 
            different (i.e. wider domain) than domain, or domain has not been explictly specified.
            Will be used to find full pressure_dir path to Pressure_*.txt files.
        pressure_max_days (int, optional) : 
            Maximum number of days from time within ds to search for the relevant pressure data. 
            Default = 31 (days).
        pressure_day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        write_nc (bool, optional) : 
            Write output .nc files (one per bin along time axis by default)
            Default = False.
        output_directory (str, optional) : 
            Top level directory to write output. Full path will be based on the "network" related to "site" 
        write_name (bool, optional) : 
            Write output .csv file for input into NAME (one per bin along time axis by default)
        name_directory (str, optional) : 
            Top level directory to write files for NAME. Full path will be based on the "network" related to "site" 
            Default defined at the top of this module.
        file_per_day (bool, optional) : 
            Group together all points into create one output file per day rather than one per data point.
            Default = False.
        overwrite (bool, optional) : 
            Allow any files already present within the full path to output_directory and name_directory 
            to be overwritten.
            Default = False.
    
    Returns:
        xarray.Dataset : merged dataset for all processed GOSAT files
        
        If write_nc is set to True:
            Set of .nc files will be written within path started by output_directory
        If write_name is set to True:
            Set of .csv files will be written within path started by name_directory
    '''

    
    if input_directory.find("$DATA_PATH"):
        input_directory = input_directory.replace("$DATA_PATH",data_path)
    
    input_directory_format = "year_split" # "year_split" (subdirectories are split into years) or None
    search_str = "*.nc"
    
    if input_directory_format == "year_split":
        files = extract_files_dir_split(input_directory,search_str=search_str,start=start,end=end,
                                        date_separator='')
    else:
        files = extract_files(input_directory,search_str,start=start,end=end,date_separator='')
    
    if domain and not (lat_bounds and lon_bounds):
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if (use_name_pressure or name_sp_filt) and domain and not pressure_domain:
        pressure_domain = domain
    
    if not write_nc and not write_name:
        print "\n**** WARNING: GOSAT_PROCESS IS SET TO NOT OUTPUT ANYTHING TO FILE. ****"
        print "Please cancel this run and restart with write_nc or write_name "
        print "parameters set to True if you want to write the processed output to disk."
        print "*************************************************************************"
    
    if len(files) > 0:
        for i,filename in enumerate(files):
            if i == 0:
                verbose = True
            else:
                verbose = False

            ds = gosat_process_file(filename,
                                    site=site,
                                    species=species,
                                    lat_bounds=lat_bounds,
                                    lon_bounds=lon_bounds,
                                    coord_bin=coord_bin,
                                    quality_filt=quality_filt,
                                    bad_pressure_filt=bad_pressure_filt,
                                    name_sp_filt=name_sp_filt,
                                    name_filters=name_filters,
                                    cutoff=cutoff,
                                    layer_range=layer_range,
                                    mode=mode,
                                    use_name_pressure=use_name_pressure,
                                    pressure_base_dir=pressure_base_dir,
                                    pressure_domain=pressure_domain,
                                    pressure_max_days=pressure_max_days,
                                    pressure_day_template=pressure_day_template,
                                    write_nc=write_nc,
                                    output_directory=output_directory,
                                    write_name=write_name,
                                    name_directory=name_directory,
                                    file_per_day=file_per_day,
                                    overwrite=overwrite,
                                    verbose=verbose)
            
            if i == 0:
                gosat = ds
            else:
                gosat = gosat.merge(ds)
                
    else:
        print "No gosat files found within directory: {0} using search_str: {1} for dates {2} - {3}".format(input_directory,search_str,start,end)
        gosat = None
    
    return gosat       
    

                 
        
