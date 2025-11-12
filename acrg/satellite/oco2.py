# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 10:25:39 2017

@author: rt17603

This module provides a set of functions for processing OCO2 data.

Within this module there are two main summary functions used for processing OCO2 CO2 Level 2 Data Product files:
    * oco2_process(directory, ...) - process a directory of files (Require more work)
    * oco2_process_file(filename, ...) - process an individual file

Currently the file processinf can be done using oco2_process_file function only.
The processing includes options for filtering based on quality flags, pressure levels, comparison with NAME surface pressure,
binning based on latitude and longitude and writing output files in netCDF format and/or text format for input into NAME.

To process multiple files the oco2_process_file function can be called within a loop over files in a directory.

Below is the example showing processing of one file:

STEP 1: 
    from acrg.satellite.oco2 import oco2_process_file
    
STEP 2: 

    NOTE: output_directory must be changed to a suitable path for writing output files.

    ds = oco2_process_file("/group/chem/acrg/obs_raw/OCO2/oco2_LtCO2_150131_B11210Ar_240820000231s.nc4",
                            quality_filt=False,
                            bad_pressure_filt=False,
                            site="oco2-CHINA",
                            lat_bounds=[3.97,53.55],
                            lon_bounds=[73.5,134.77],
                            coord_bin=[0.234,0.352],
                            pressure_domain="SOUTHASIA",
                            write_nc=True,
                            write_name=True,
                            output_directory="/group/chem/acrg/prasad/acrg_oco2_processed") 

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

from acrg.satellite.common import extract_files_oco2, name_pressure_file, name_pressure_match,apply_filter,coord_order,define_pressure_levels,distance_lat,distance_lon,output_name,output,add_coords,ds_check_internal_unique,add_history_attr

from acrg.satellite.gosat import latlon_filter, binned_mean, extract_files_dir_split

data_path = Paths.data

home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/OCO2/")
fp_directory = os.path.join(data_path,'LPDM/fp_NAME/')
obs_directory = os.path.join(data_path,'obs/OCO2/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files/OCO2") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure/")

data_path = Paths.data


def oco2_quality_filter(ds):
    '''
    The oco2_quality_filter function filters all data variables within a dataset by the "xco2_quality_flag" variable.
    The xco2_quality_flag equates to:
        - 0 indicates good_quality
        - 1 indicates potentially_bad_quality
    Note: the Dataset is filtered assuming all data variables have the same first dimension (will be true for oco2 data
    but would need to check for other data)

    Args:
        ds (xarray.Dataset) : 
            OCO2 CO2 Level 2 Data Product file opened with xarray as a Dataset.
                                
    Returns:
        ds (xarray.Dataset) : 
            Filtered Dataset with indices associated with xch4_quality_flag=1 removed.
    '''

    dim_order,consistent = coord_order(ds,check_consistent=True)
    if not consistent:
        print('WARNING: Order of dimensions is not consistent. Unable to apply oco2_quality_filter function')
        return None
    
    species = 'co2'
    filter_name = "x" + species.lower() + "_quality_flag"
    flag = 0
    #dim_apply = "time"

    #import pdb
    #pdb.set_trace()

    # Filter based on filter_name and flag
    try:
        # Note: adds dimension associated with where condition if not already
        # present. Should be fine for oco2 data (all variables should have
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

def oco2_pressure_filter(ds):
    '''
    The oco2_pressure_filter function removes all points which have at least one pressure_level undefined (-9999.99).
    Unable to use these points to calculate relevant pressure levels and difference between levels within the atmosphere.
    
    Args:
        ds (xarray.Dataset) : 
            OCO2 CO2 Level 2 Data Product file opened with xarray as a Dataset.
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
        # present. Should be fine for oco2 data (all variables should have
        # a first order of "time") but may need to changed for other data.
        ds_new = ds.where(~(ds[filter_name] == exclude).any(axis=1),drop=True) 
    except KeyError:
        raise Exception("Unable to apply filter to {0} to exclude values of {1}. Input dataset does not contain data variable {0}".format(filter_name,exclude)) 

    # Add attribute describing modification made to original data
    mod_attr = "bad_pressure_filter"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Data points with values of {1} for any of the {0} values have been removed.".format(filter_name,exclude)

    return ds_new

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
        This is correct for the oco2 Level 2 Product based on the CCI GHG Documentation for level-based extraction: 
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


def oco2_process_file(filename,site,species="co2",lat_bounds=[],lon_bounds=[],domain=None,
                       coord_bin=None,quality_filt=True,bad_pressure_filt=True,
                       name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],                       
                       mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31.,pressure_day_template=True,
                       write_nc=False,output_directory=obs_directory,
                       write_name=False,name_directory=name_csv_directory,
                       file_per_day=False,max_name_level=17,max_name_points=None,overwrite=False,
                       verbose=True):
    '''
    The oco2_process_file function processes input oco2 data for one file, applying designated filters, 
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
        - netCDF files (.nc) based on suitable format for acrg repository. See oco2_output() and oco2_output_filename() functions.
        - text files  (.csv) suitable for input into NAME. See output_name() function.
    
    Args:
        filename (str) : 
            Filename of oco2 "co2" Level 2 Data Product file (str)
        site (str) : 
            Specified sub-set defined for oco2 e.g. OCO2-INDIA. Should be defined within site_info.json
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
            Whether to remove data points based on the comparison between the oco2 surface pressure value
            and the NAME surface pressure value. 
            Exact conditions are specified by name_filters with associated values specified by cutoff
            and layer_range parameters if applicable.
            Default = False.
        mode (str, optional) : 
            Which retrieval mode (retr_flag) to filter by (if any). Options: 'land' or 'glint'.
            Set to None if no filtering by mode is required.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the input oco2 surface pressure
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
            processed oco2 Dataset
        
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
    oco2 = xray.open_dataset(filename)

    # Set 'time' as the index for 'sounding_id'
    oco2 = oco2.set_index({'sounding_id': 'time'})
    # Swap 'sounding_id' with 'time' to make 'time' the dimension
    oc2 = oco2.swap_dims({'sounding_id': 'time'})
    # Rename 'sounding_id' to 'time' (optional, if needed)
    oco2 = oco2.rename({'sounding_id': 'time'})
    
    # Rename levels to lev 
    oco2 = oco2.rename({'levels': 'lev'})
    oco2 = oco2.sortby(axis)
    #oco2 = ds_check_internal_unique(oco2,axis) # Check time values are unique and slightly modify if necessary

    if quality_filt:
        if verbose:
            print('Applying filter based on quality flag')
        oco2 = oco2_quality_filter(oco2)
    if bad_pressure_filt:
        if verbose:
            print('Applying filter based on bad pressure flag')
        oco2 = oco2_pressure_filter(oco2) # Removes points with any pressure values of -9999.99

    if name_sp_filt:
        if len(oco2[axis].values) > 0:
            if not name_filters:
                raise Exception("If name_sp_filt=True, name_filters must be specified.")
            if verbose:
                print('Applying filters based on NAME surface pressure.')
            
            oco2 = name_pressure_filter(oco2,filters=name_filters,pressure_domain=pressure_domain,
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
        lat_bounds = [min(oco2.latitude.values),max(oco2.latitude.values)]
    if not lon_bounds:
        lon_bounds = [min(oco2.longitude.values),max(oco2.longitude.values)]
    
    if coord_bin:
        if verbose:
            print('Binning data based on {0} degree bins'.format(coord_bin))
            print('Looking at area within latitude and longitude bounds: {0},{1}'.format(lat_bounds,lon_bounds))
        oco2 = binned_mean(oco2,lat_bounds,lon_bounds,domain,columns=columns,coord_bin=coord_bin)
        oco2 = oco2.sortby(axis)
        oco2 = ds_check_internal_unique(oco2,axis=axis) # After binning check time values are unique again
    else:
        if verbose:
            print('Looking at area within latitude and longitude bounds: {0},{1}'.format(lat_bounds,lon_bounds))
        oco2 = latlon_filter(oco2,lat_bounds,lon_bounds,columns=columns)

    oco2.attrs["input_file"] = os.path.split(filename)[1]

    if len(oco2[axis].values) > 0:
        if write_nc:
            output(oco2,site=site,species=species, network="oco2",output_directory=output_directory,
                   file_per_day=file_per_day,overwrite=overwrite)
        
        if write_name:
            output_name(oco2,network='oco2',site=site,max_level=max_name_level,use_name_pressure=use_name_pressure,
                        pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                        pressure_max_days=pressure_max_days,pressure_day_template=pressure_day_template,
                        name_directory=name_directory,file_per_day=file_per_day,
                        max_points=max_name_points,overwrite=overwrite)
    else:
        print('No points extracted from {} for specified parameters.'.format(filename))

    return oco2


def oco2_process(site,species="co2",input_directory=input_directory,start=None,end=None,
                  lat_bounds=[],lon_bounds=[],domain=None,coord_bin=None,
                  quality_filt=True,bad_pressure_filt=True,
                  name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],
                  mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                  pressure_domain=None,pressure_max_days=31,pressure_day_template=True,
                  write_nc=False,output_directory=obs_directory,
                  write_name=False,name_directory=name_csv_directory,file_per_day=False,
                  max_name_level=17,max_name_points=None,overwrite=False):
    '''
    The oco2_process function processes oco2 input files from a directory.
    As part of processing the oco2 input data there are multiple options including filtering based on
    various criteria and binning. See oco2_process_file() function for more details.
    
    Args:
        site (str) : 
            Specified sub-set defined for oco2 e.g. oco2-INDIA. Should be defined within site_info.json
        species (str, optional) : 
            Species of interest. Should be defined within species_info.json
            Default = "ch4"
        input_directory (str, optional) : 
            Top level directory containing oco2 CH4 Level 2 Data Product files only.
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
            Whether to remove data points based on the comparison between the oco2 surface pressure value
            and the NAME surface pressure value. 
            Exact conditions are specified by name_filters with associated values specified by cutoff
            and layer_range parameters if applicable.
            Default = False.
        mode (str, optional) : 
            Which retrieval mode (retr_flag) to filter by (if any). Options: 'land' or 'glint'.
            Set to None if no filtering by mode is required. 
            Default = None.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the input oco2 surface pressure
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
        xarray.Dataset : merged dataset for all processed oco2 files
        
        If write_nc is set to True:
            Set of .nc files will be written within path started by output_directory
        If write_name is set to True:
            Set of .csv files will be written within path started by name_directory
    '''

    
    if input_directory.find("$DATA_PATH"):
        input_directory = input_directory.replace("$DATA_PATH",str(data_path))
    
    input_directory_format = "monthly" # "year_split" (subdirectories are split into years) or None
    search_str = "*.nc"
    files = extract_files_oco2(input_directory,search_str,start=start,end=end,date_separator='')
    
    if domain and not (lat_bounds and lon_bounds):
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if (use_name_pressure or name_sp_filt) and domain and not pressure_domain:
        pressure_domain = domain
    
    if not write_nc and not write_name:
        print("\n**** WARNING: oco2_PROCESS IS SET TO NOT OUTPUT ANYTHING TO FILE. ****")
        print("Please cancel this run and restart with write_nc or write_name ")
        print("parameters set to True if you want to write the processed output to disk.")
        print("*************************************************************************")
    
    if len(files) > 0:
        for i,filename in enumerate(files):
            if i == 0:
                verbose = True
            else:
                verbose = False

            ds = oco2_process_file(filename,
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
                oco2 = ds
            else:
                oco2 = oco2.merge(ds)
                
    else:
        print("No oco2 files found within directory: {0} using search_str: {1} for dates {2} - {3}".format(input_directory,search_str,start,end))
        oco2 = None
    
    return oco2       
    

                 
        
