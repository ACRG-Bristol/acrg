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
import numpy as np
import xarray as xray
from collections import OrderedDict

from .barometric import pressure_at_height
from acrg.countrymask import domain_volume
from acrg.config.paths import Paths

from acrg.satellite.common import extract_files, name_pressure_file, name_pressure_match,apply_filter,coord_order,define_pressure_levels,distance_lat,distance_lon,output_name,output,add_coords,ds_check_internal_unique,add_history_attr

data_path = Paths.data

home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/")
fp_directory = os.path.join(data_path,'LPDM/fp_NAME/')
obs_directory = os.path.join(data_path,'obs/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure/")

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

    #import pdb
    #pdb.set_trace()

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


# def gosat_add_coords(ds,data_vars=[]):
#     '''
#     The gosat_add_coords function adds co-ordinates to a GOSAT dataset.
#     This is valid for (at least) v6 of GOSAT CH4 Proxy Level 2 Data Product downloaded from the CCI Open Data Portal
#     (e.g. filename of the form: ESACCI-GHG-L2-CH4-GOSAT-OCPR-20101231-fv6.nc)
    
#     This assumes the input Dataset will be of the form:
#         <xarray.Dataset>
#         Dimensions:                           (m: 20, n: 1962)
#         Dimensions without coordinates: m, n
#         Data variables:
#             xch4_quality_flag                 (n) int8 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
#             ch4_profile_apriori               (n, m) float32 1751.95 1751.95 1751.94 ...
#             xch4                              (n) float32 1644.78 1657.72 1648.04 ...
#             ...
#         Attributes:
#             title:                     ESA CCI GOSAT OCPR CH4
#             institution:               University of Leicester (UoL), UK
#             ...
    
#     i.e. with no coordinates assigned for the two dimensions within the Dataset.
#     This function redefines "m" as "time" and "n" as "level" and adds the associated values as coords.

#     Args:
#         ds (xarray.Dataset) : 
#             GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
#         data_vars (list) : 
#             Data variables which will have the time and level (where appropriate) dataset coords added
#             If data_vars if left blank ([] by default) coords will be added to all data variables within ds.
                              
#     Returns:
#         ds (xarray.Dataset) : 
#             Dataset with time and level coordinates assigned
#             (only includes specified data variables unless left blank)
#     '''
    
#     dims_current = ['n','m']
#     dims = ['time','lev']
    
#     if list(ds.dims.keys()) != dims_current and list(ds.dims.keys()) != dims_current[::-1]: # Check dimensions match what we expect (check forward and reverse list)
#         print('WARNING: Do not recognise dimensions of input gosat dataset. Unable to add new dimensions.')
#         return None
    
#     # Define co-ordinate values. Extract from Dataset if already present (e.g. time), construct if not (e.g. lev)
#     coords = {}
#     for i,d in enumerate(dims):
#         try:
#             coords[d] = ds[d].values # Extract values for each coord from Dataset if present (e.g. time)
#         except (AttributeError,KeyError):
#             current_dim = dims_current[i]
#             dim_length = ds.dims[current_dim]
#             coords[d] = np.arange(dim_length) # If no values defined, create integer list based on dimension size
    
#     if not data_vars:
#         data_vars = ds.data_vars
    
#     data = OrderedDict()
#     for name in data_vars:
#         if name not in dims:
#             if len(ds[name].dims) == 1:
#                 data[name] = (dims[0],ds[name].values)
#             elif len(ds[name].dims) == 2:
#                 data[name] = ([dims[0],dims[1]],ds[name].values)
    
#     #coord_dict = {dims[0]:dim1,dims[1]:dim2}
    
#     ds_new = xray.Dataset(data,coords=coords,attrs=ds.attrs)
#     for dv in ds_new.data_vars:
#         ds_new[dv].attrs = ds[dv].attrs
#     for coord in ds_new.coords:
#         try:
#             ds_new[coord].attrs = ds[coord].attrs
#         except (AttributeError,KeyError):
#             if coord == "lev" or coord == "level":
#                 ds_new[coord].attrs["short_description"] = "Number for each level within the vertically resolved data."
#                 ds_new[coord].attrs["long_name"] = "level"
    
#     # Add attribute describing modification made to original data
#     mod_attr = "nc_restructure"
#     ds_new = add_history_attr(ds_new,mod_attr)
#     ds_new.attrs[mod_attr] = "Dataset rearranged to re-assign {} dimensions as {} coordinates, extracting values from dataset where present.".format(dims_current,dims).replace("'","")
    
#     return ds_new

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
        ds[filter_name] = ds[filter_name].astype(dtype="str") # Update for new xarray version? - Have to recast as explicit string type
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
        all_dir_years = [d for d in all_dir if len(d) == 4 and re.match(r"\d{4}",d)]
        year_start = min(all_dir_years)
        year_end = max(all_dir_years)
        print("No start and end date specified, so processing all files from all date labelled sub-directories within input directory {} for year range: {}-{}".format(directory,year_start,year_end))
    year_range = list(range(int(year_start),int(year_end)+1))
    
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
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
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
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
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

def gosat_process_file(filename,site,species="ch4",lat_bounds=[],lon_bounds=[],domain=None,
                       coord_bin=None,quality_filt=True,bad_pressure_filt=True,
                       name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],                       
                       mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31.,pressure_day_template=True,
                       write_nc=False,output_directory=obs_directory,
                       write_name=False,name_directory=name_csv_directory,
                       file_per_day=False,max_name_level=17,max_name_points=None,overwrite=False,
                       verbose=True):
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
        - text files  (.csv) suitable for input into NAME. See output_name() function.
    
    Args:
        filename (str) : 
            Filename of GOSAT CH4 Level 2 Data Product file (str)
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within site_info.json
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
            Default = "$DATA_PATH/LPDM/surface_pressure/"
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
        max_name_level (int, optional) :
        	Maximum level to use when writing out NAME csv files.
        	Default = 17
        max_name_points (int/None, optional) :
            Only applicable if file_per_day is True.
            Maximum number of points to write to NAME csv file if writing file per day. 
            If number of points within the day exceeds max_points multiple files will be written out with 
            a character after the date of "-A","-B",...
            NOTE: This will not be applied to netCDF files even if number of points is greater than
            max_name_points, only to the NAME csv files. All points per day will be written to one file.
            A sensible number for this, if specified, would be 60.
            Default = None.
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
    
    gosat = add_coords(gosat,network='GOSAT')
    gosat = gosat.sortby(axis)
    #gosat = ds_check_internal_unique(gosat,axis) # Check time values are unique and slightly modify if necessary
    
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
        if len(gosat[axis].values) > 0:
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
            print("Extracting latitude and longitude bounds from footprints associated with domain: {}".format(domain))
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
            output(gosat,site=site,species=species,output_directory=output_directory,
                   file_per_day=file_per_day,overwrite=overwrite)
        
        if write_name:
            output_name(gosat,network='GOSAT',site=site,max_level=max_name_level,use_name_pressure=use_name_pressure,
                        pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                        pressure_max_days=pressure_max_days,pressure_day_template=pressure_day_template,
                        name_directory=name_directory,file_per_day=file_per_day,
                        max_points=max_name_points,overwrite=overwrite)
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
                  max_name_level=17,max_name_points=None,overwrite=False):
    '''
    The gosat_process function processes GOSAT input files from a directory.
    As part of processing the GOSAT input data there are multiple options including filtering based on
    various criteria and binning. See gosat_process_file() function for more details.
    
    Args:
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within site_info.json
        species (str, optional) : 
            Species of interest. Should be defined within species_info.json
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
            Default = "$DATA_PATH/LPDM/surface_pressure/"
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
        max_name_level (int, optional) :
            	Maximum level to use when writing out NAME csv files.
            Default = 17
        max_name_points (int/None, optional) :
            Only applicable if file_per_day is True.
            Maximum number of points to write to NAME csv file if writing file per day. 
            If number of points within the day exceeds max_points multiple files will be written out with 
            a character after the date of "-A","-B",...
            NOTE: This will not be applied to netCDF files even if number of points is greater than
            max_name_points, only to the NAME csv files. All points per day will be written to one file.
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
        input_directory = input_directory.replace("$DATA_PATH",str(data_path))
    
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
        print("\n**** WARNING: GOSAT_PROCESS IS SET TO NOT OUTPUT ANYTHING TO FILE. ****")
        print("Please cancel this run and restart with write_nc or write_name ")
        print("parameters set to True if you want to write the processed output to disk.")
        print("*************************************************************************")
    
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
                                    max_name_level=max_name_level,
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
                                    max_name_points=max_name_points,
                                    overwrite=overwrite,
                                    verbose=verbose)
            
            if i == 0:
                gosat = ds
            else:
                gosat = gosat.merge(ds)
                
    else:
        print("No gosat files found within directory: {0} using search_str: {1} for dates {2} - {3}".format(input_directory,search_str,start,end))
        gosat = None
    
    return gosat       
    

                 
        
