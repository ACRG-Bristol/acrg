#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 08:30:03 2019

@author: al18242

Functions for processing and using TROPOMI data
"""

import os
import socket
import numpy as np
import xarray as xr
import time as timing_module
from pathlib import Path
from collections import OrderedDict
import acrg_grid.regrid_xesmf as regrid
import acrg_satellite.gosat as gosat_fn

#from acrg_satellite.gosat import latlon_filter,extract_files,find_network

try:
    from acrg_config.paths import paths
except ImportError:
    data_path = os.getenv("DATA_PATH")
else:
    data_path = paths.data

if not data_path:
    server_name = socket.gethostname()
    if "bp1" in server_name:
        data_path = "/work/chxmr/shared/"
    elif "bc4" in server_name:
        data_path = "/mnt/storage/private/acrg/ACRG_Repository/data/"
    elif "snowy" in server_name:
        data_path = "/data/shared/"
    else:
        print("Unable to infer data_path")

home = Path(os.getenv("HOME"))

##TODO: Replace this with a central location when this has been sorted
input_directory = Path("/home/rt17603/shared/obs_raw/TROPOMI/SOUTHAMERICA/")

name_csv_directory = os.path.join(home,"NAME_files") # Where to write output NAME csv files
obs_directory = os.path.join(data_path,'obs/') # Where to write output nc files

def cornerGrid(geoData):
    '''Create lat and lon corner grids (for pcolormesh+xESMF etc) from geoData
    input
        geoData - dataframe from group "PRODUCT/SUPPORT_DATA/GEOLOCATIONS" of TROPOMI netCDF file
    
    return
        lat - 2d array of latitude grid corners
        lon - 2d array of longitude grid corners
    '''
    lat = np.zeros([geoData.scanline.size +1, geoData.ground_pixel.size +1])
    lon = np.zeros([geoData.scanline.size +1, geoData.ground_pixel.size +1])
    
    #lower left block is just lower left corner of all points
    lat[:-1, :-1] = geoData.latitude_bounds[0,:,:,0]
    lon[:-1, :-1] = geoData.longitude_bounds[0,:,:,0]
    
    #top row
    lat[-1, :-1] = geoData.latitude_bounds[0,-1,:,3]
    lon[-1, :-1] = geoData.longitude_bounds[0,-1,:,3]
    #right col
    lat[:-1, -1] = geoData.latitude_bounds[0,:,-1,1]
    lon[:-1, -1] = geoData.longitude_bounds[0,:,-1,1]
    #corner
    lat[-1,-1] = geoData.latitude_bounds[0, -1, -1, 2]
    lon[-1,-1] = geoData.longitude_bounds[0, -1, -1, 2]
       
    return lat, lon

def preProcessFile(filename,edge_coords="corners"):
    '''Combine important information from different netCDF groups into a single data frame
    
    input
        filename - string of full filename for TROPOMI file
        edge_coords - "corners" or "bounds"
    
    return
        tropomi_data - dataframe containing mixing ratio, quality flags, grid data, and column data
                        the "qa_pass" flag variable is True for good measurements and False for bad measurements
    '''
    tropomi_data = xr.open_dataset(filename, group="PRODUCT")
    tropomi_data_aux = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS")
    tropomi_data_geo = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS")
    tropomi_data_input = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA")
    
    #calculate boundaries of pressure bands from surface pressure and constant intervals
    pressure_data = (np.expand_dims(tropomi_data_input.surface_pressure,axis=3) - \
                                np.expand_dims(tropomi_data_input.pressure_interval,axis=3) * \
                                np.reshape(np.arange(0,13),newshape=(1,1,1,-1)))
    tropomi_data['pressure_levels'] = (["time", "scanline", "ground_pixel", "layer_bound"], pressure_data)
    
    tropomi_data['column_averaging_kernel'] = tropomi_data_aux.column_averaging_kernel
    tropomi_data['methane_profile_apriori']= tropomi_data_input.methane_profile_apriori
    
    if edge_coords == "corners":
        #calculate the meshgrid of corner points in the scan grid for further processing
        lat, lon = cornerGrid(tropomi_data_geo)
        #tropomi_data['lat_corners'] = (["time", "scanline_c", "ground_pixel_c"], np.expand_dims(lat,0))
        #tropomi_data['lon_corners'] = (["time", "scanline_c", "ground_pixel_c"],np.expand_dims(lon,0))
        tropomi_data = tropomi_data.assign_coords({'lat_corners':(["time", "scanline_c", "ground_pixel_c"], np.expand_dims(lat,0))})
        tropomi_data = tropomi_data.assign_coords({'lon_corners':(["time", "scanline_c", "ground_pixel_c"],np.expand_dims(lon,0))})
        
    elif edge_coords == "bounds":
        #Explicitly add bounds (may want to change this)
        tropomi_data = tropomi_data.assign_coords({'latitude_bounds':tropomi_data_geo['latitude_bounds']})
        tropomi_data = tropomi_data.assign_coords({'longitude_bounds':tropomi_data_geo['longitude_bounds']})
    
    #flag quality to remove outputs not flagged as 'success' and use recommended acceptable qa_values
    # Note: it may be the case that any value that passes the second condition necessarily passes the first
    quality_flag = (tropomi_data_aux.processing_quality_flags.values.astype(int) & 0xFF == 0) & \
                    (tropomi_data.qa_value.values > 0.5)
    
    tropomi_data["qa_pass"] = (["time", "scanline", "ground_pixel"], quality_flag)
    
    return tropomi_data


def check_in_grid(ds,output_lat,output_lon,latlon=["latitude","longitude"]):
    '''
    Check if any measurements are within the latitude-longitude range imposed.
    If the input grid is not within the output grid, an error is thrown when accessing 
    the regrid_xesmf functions.
    
    Args:
        ds (xarray.Dataset) :
            Dataset with latitudes and longitudes as data variables
        output_lat (numpy.array), output_lon (numpy.array) :
            Latitudes and longitudes for the values to use for filtering (each 1D)
        latlon (tuple, optional) :
            Names for latitude and longitude within ds.
            Default = ("latitude", "longitude")
           
    Returns:
        bool :
            Whether any tropomi orbit points are within the output_lat and output_lon range
    '''

    lat_bounds = [output_lat[0],output_lat[-1]]
    lon_bounds = [output_lon[0],output_lon[-1]]

    ds_check = gosat_fn.latlon_filter(ds,lat_bounds,lon_bounds,columns=latlon)

    if len(ds_check[latlon[0]]) == 0 and len(ds_check[latlon[1]]) == 0:
        inside_grid = False
    else:
        inside_grid = True
    
    return inside_grid

def all_nan(da):
    '''
    Check if all values within a given DataArray are NaNa.
    Returns: True (if all are NaN), False (if at least 1 value is not a NaN)
    '''
    return np.all(np.isnan(da.values))

def create_regrid_mask(da, pos_coords=("scanline", "ground_pixel")):
    '''
    Create mask for use with conservative_normed method. This should 
    match to the shape of the input array for regridding with boolean 
    values indicating which positions contain NaN values.
    
    Returns: numpy.array of boolean values
    '''
    dims = da.dims
    da_copy = da.copy(deep=True)
    
    for dim in dims:
        if dim not in pos_coords:
            da_copy = da_copy.isel({dim:0})
    
    input_values = da_copy.rename({pos_coords[0]:'x',pos_coords[1]:'y'})
    input_values = np.isfinite(input_values)
    
    return input_values


def regrid_da(da_tropomi,output_lat,output_lon,ds_tropomi_geo,
              method="conservative",latlon=["latitude","longitude"],
              pos_coords=("ground_pixel","scanline"),
              set_nan=True,create_corners=False,reuse_weights=False):
    '''
    The regrid_da function regrids one variable from a tropomi orbit onto a new 
    output grid.

    Inputs:
        da_tropomi (xarray.DataArray) :
            DataArray for one data variable from a Tropomi dataset. 
            From an opened tropomi file (netcdf group "PRODUCT") or created using 
            tropomi.preProcessFile() function (adds data from multiple groups).
        output_lat (numpy.array), output_lon (numpy.array) :
            Latitudes and longitudes for the output grid (each 1D)
        ds_tropomi_geo (xarray.Dataset) :
            Tropomi data from netcdf group "PRODUCT/SUPPORT_DATA/GEOLOCATIONS". 
            This dataset needs to contain "lat_bounds" and "lon_bounds" data variables.
        method (str, optional)
            string describing method to use when regridding. Should be one of:
                "conservative"
                "conservative_normed".
            Note that to use the conservative_normed method you need to have installed the 
            "masking" development branch of the xesmf package (contact Daniel for more details).
        latlon (tuple, optional) :
            Names for latitude and longitude within ds_tropomi dataset.
            Default = ("latitude", "longitude")
        pos_coords (tuple, optional) :
            Names for tropomi position co-ordinates within the orbit file.
            Default = ("scanline", "ground_pixel")
        set_nan (bool, optional) :
            After regridding, explicitly set values of 0 to np.nan.
            Default = True
        create_corners (bool, optional) :
            Re-create corner co-ordinates based on current "lat_bounds" and "lon_bounds"
            in the input file. This may be useful if the data has been filtered after
            the dataset was first created.
            Default = True
        reuse_weights (bool, optional) :
            Reuse previous weights created over the same grid area. This speeds
            up computation but should make sure this exactly matches the previous
            grid shape and any nan values.
            Default = False
    '''

    input_lat = da_tropomi[latlon[0]][0,:,:].values
    input_lon = da_tropomi[latlon[1]][0,:,:].values
    
    if ("lat_corners" not in ds_tropomi_geo and "lon_corners" not in ds_tropomi_geo) or (create_corners):
        input_lat_b, input_lon_b = cornerGrid(ds_tropomi_geo)
    else:
        input_lat_b = ds_tropomi_geo["lat_corners"].isel(time=0).values
        input_lon_b = ds_tropomi_geo["lon_corners"].isel(time=0).values

    input_grid = xr.Dataset({'lat': (['x', 'y'], input_lat),
                         'lon': (['x', 'y'], input_lon),
                         'lat_b': (['x_b', 'y_b'], input_lat_b),
                         'lon_b': (['x_b', 'y_b'], input_lon_b)})

    if method == "conservative_normed":
        ## mask is only needed for "conservative_normed" method
        #dims = da_tropomi.dims
        #da_copy = da_tropomi.copy(deep=True)
        ## Basing mask on first index for all extra dimensions (time, layer etc.)
        ## *For time this is fine but for other coords this assumption should be checked to be valid!*
        #for dim in dims:
        #    if dim not in pos_coords:
        #        da_copy = da_copy.isel({dim:0})
        #
        
        #
        #input_values = da_copy.rename({pos_coords[0]:'x',pos_coords[1]:'y'})
        
        ## mask is only needed for "conservative_normed" method
        nan_mask = create_regrid_mask(da_tropomi)
        input_grid = input_grid.assign(**{'mask': nan_mask})
        
    output_grid = regrid.create_xesmf_grid_uniform_cc(output_lat, output_lon)
    regridded = regrid.regrid_betweenGrids(da_tropomi, input_grid, output_grid,
                                           method=method,
                                           reuse_weights=False)#reuse_weights)
    if set_nan:
        regridded.values[regridded.values == 0.] = np.nan
    
    regridded = regridded.drop(labels=("lat","lon"))
    regridded = regridded.assign_coords(**{"x":output_lat,"y":output_lon})
    regridded = regridded.rename({"x":"lat","y":"lon"})

    return regridded


def regrid_orbit(ds_tropomi,output_lat,output_lon,names=None,ds_tropomi_geo=None,
                 method="conservative",latlon=("latitude", "longitude"),
                 pos_coords=("scanline", "ground_pixel"),set_nan=True,
                 filter_latlon=True,create_corners=True,verbose=False):
    '''
    The regrid_orbit function takes a tropomi orbit and regrids this onto a new 
    output grid.
    
    Inputs:
        ds_tropomi (xarray.Dataset) :
            Tropomi dataset. Opened tropomi file (netcdf group "PRODUCT") or 
            created using tropomi.preProcessFile() function (adds data from multiple groups).
        output_lat (numpy.array), output_lon (numpy.array) :
            Latitudes and longitudes for the output grid (each 1D)
        names (str/list, optional) :
            Name(s) of specific data variable(s) within dataset to regrid
            If this is set to None, all variables (except excluded variables)
        ds_tropomi_geo (xarray.Dataset, optional) :
            Tropomi data from netcdf group "PRODUCT/SUPPORT_DATA/GEOLOCATIONS". This dataset
            needs to contain "lat_bounds" and "lon_bounds" data variables.
            If tropomi.preProcessFile() has been used to create ds_tropomi this input is not
            necessary but is needed otherwise.
            Default = None
        method (str, optional)
            string describing method to use when regridding. Should be one of:
                "conservative"
                "conservative_normed".
            Note that to use the conservative_normed method you need to have installed the 
            "masking" development branch of the xesmf package (contact Daniel for more details).
        latlon (tuple, optional) :
            Names for latitude and longitude within ds_tropomi dataset.
            Default = ("latitude", "longitude")
        pos_coords (tuple, optional) :
            Names for tropomi position co-ordinates within the orbit file.
            Default = ("scanline", "ground_pixel")
        set_nan (bool, optional) :
            After regridding, explicitly set values of 0 to np.nan.
            Default = True
        filter_latlon (bool, optional) :
            Filter the input data to only include the range covered by the output grid.
            May make computation quicker.
            Default = True
        create_corners (bool, optional) :
            Re-create corner co-ordinates based on current "lat_bounds" and "lon_bounds"
            in the input file. This may be useful if the data has been filtered after
            the dataset was first created.
            Default = True
        verbose (bool, optional) :
            Print detailed output at each stage (verbose).
            Default = True
            
        
    Returns:
        xarray.Dataset :
            Regridded data for tropomi orbit
    '''

    # Check input grid overlap with the output grid either by explicitly filtering
    # the input data or by checking.
    if filter_latlon:
        lat_diff = np.abs(np.mean(output_lat[1:] - output_lat[:-1]))
        lon_diff = np.abs(np.mean(output_lon[1:] - output_lon[:-1]))
        lat_bounds = [output_lat[0]-lat_diff/2,output_lat[-1]+lat_diff/2]
        lon_bounds = [output_lon[0]-lon_diff/2,output_lon[-1]+lon_diff/2]
        
        ds_tropomi = gosat_fn.latlon_filter(ds_tropomi,lat_bounds,lon_bounds,columns=latlon)
        
        if len(ds_tropomi[latlon[0]]) == 0 and len(ds_tropomi[latlon[1]]) == 0:
            if verbose:
                print("No points in tropomi orbit within output grid")
            return None
    else:
        inside_grid = check_in_grid(ds_tropomi,output_lat,output_lon)
        if not inside_grid:
            if verbose:
                print("No points in tropomi orbit within output grid")
            return None

    # Create ds_tropomi_geo variable if it doesn't exist
    if not ds_tropomi_geo:
        ds_tropomi_geo = ds_tropomi
    
    # Expected dimensions in tropomi dataset both in order and name
    expected_dims = ("time",) + pos_coords
    
    if names:
        if isinstance(names,str):
            names = [names]
    else:
        # Create names if not defined
        names = ds_tropomi.data_vars
    
    excluded_dv = ["time_utc"]
    names = [name for name in names if name not in excluded_dv]
    
    # For each data variable apply the regridding algorithm
    for i,name in enumerate(names):
        ## TODO: Add way to check whether there are same nan values between 
        # one variable to the next. If so, you can reuse the weights (faster
        #computationally)
        
        dv = ds_tropomi[name]
        dims = dv.dims
        
        # If data variables has an additional dimension this needs to be transposed
        # to make sure the position coordinates are the last dimensions.
        # This is an expectation of the xesmf functionality
        if dims != expected_dims:
            if dims[0:len(expected_dims)] == expected_dims:
                extra_dimension = [dim for dim in dims if dim not in expected_dims][0]
            else:
                if verbose:
                    print(f"Cannot regrid variable {name} with dimensions {dims}")
                return None
                        
            dv = dv.transpose(*(extra_dimension,)+expected_dims)
            transpose = True
        else:
            transpose = False
        
        if all_nan(dv):
            if verbose:
                print(f"Data variable {name} has no values. Cannot regrid this orbit")
            return None

        # The same weights can be used for each data variable so can be calculated once
        # and then used again - this makes overall computation a lot faster.
        #if i == 0:
        #    reuse_weights = False
        #else:
        #    reuse_weights = True

        # The same weights can be used for each data variable so can be calculated once
        # and then used again - this makes overall computation a lot faster.
        if i == 0:
            mask_prev = create_regrid_mask(dv)
            reuse_weights = False
        else:
            mask = create_regrid_mask(dv)
            if np.array_equal(mask,mask_prev):
                print("***Reusing the same weights!***")
                reuse_weights = True
            else:
                reuse_weights = False
            mask_prev = mask

        
        regridded = regrid_da(dv,output_lat,output_lon,ds_tropomi_geo,
                              method=method,latlon=latlon,pos_coords=pos_coords,
                              set_nan=set_nan,create_corners=create_corners,
                              reuse_weights=reuse_weights)
        
        # Reorder dimensions to reflect original order
        if transpose:
            new_coord_dims = regridded.dims[-2:]
            new_dims = [new_coord_dims[0] if dim == pos_coords[0] else dim for dim in dims]
            new_dims = [new_coord_dims[1] if dim == pos_coords[1] else dim for dim in new_dims]         
            regridded = regridded.transpose(*new_dims)
        
        # Create / add to output dataset
        if i == 0:
            ds_tropomi_regridded = regridded.to_dataset(name=name)
        else:
            ds_tropomi_regridded = ds_tropomi_regridded.assign({name:regridded})
    
    ds_tropomi_regridded.attrs["regridded"] = f"Regridded to output grid of lat = {output_lat[0]} - {output_lat[-1]}, lon = {output_lon[0]} - {output_lon[-1]}"
    
    return ds_tropomi_regridded


def regrid_tropomi(ds_tropomi,lat_bounds,lon_bounds,coord_bin,
                   method="conservative",time_increment="10s"):
    '''
    The regrid_tropomi function regrids data for a given input tropomi dataset.
    
    Args:
        ds_tropomi (xarray.Dataset) :
            Tropomi dataset. Created using tropomi.preProcessFile() function 
            (adds data from multiple groups).
        lat_bounds (list) : 
            Upper and lower bounds for latitude. (two-item list e.g. [6.0,36.5])
        lon_bounds (list) : 
            Upper and lower bounds for longitude. (two-item list e.g. [70.0,90.5])
        coord_bin (float/list) : 
            The size of each bin in degrees.
            To specify the same bin for latitude and longitude a single values (float) can be used.
            To specify different bins include a two item list (e.g. [0.234,0.356]).
        method (str, optional)
            String describing regridding method to use.
            Should be one of "conservative" or "conservative_normed".
            Note that to use the conservative_normed method you need to have installed the "masking" development
            branch of the xesmf package (contact Daniel for more details).
            
            With the 'masking' branch of xESMF you can include a mask in the input_grid to ignore nan values
        time_increment (str/None, optional) :
            Time bins to split data into when regridding. This can be specified
            using pandas offset aliases:
                https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases
            For example:
                "10s" (10 seconds)
                "30min" (30 minutes)
                "1h" (1 hour)
            Should always be a day or less.
            Can set to None to not split the data by time.
            Default = "10s"
            
    Returns:
        xarray.Dataset :        
            Dataset variables re-gridded onto specified regular output grid.
    
    '''

    ## Create output grid from parameters
    if isinstance(coord_bin,float) or isinstance(coord_bin,int):
        coord_bin = [coord_bin,coord_bin]
    elif len(coord_bin) == 1:
        coord_bin = [coord_bin[0],coord_bin[0]]
    
    if len(coord_bin) != 2:
        print("WARNING: Did not recognise input for coord_bin into create_zbin function: {0}. Should be a two item object".format(coord_bin))
        return None

    dlat = coord_bin[0]
    dlon = coord_bin[1]

    output_lat = np.arange(lat_bounds[0],lat_bounds[-1]+dlat,dlat)
    output_lon = np.arange(lon_bounds[0],lon_bounds[-1]+dlon,dlon)
    
    ## Split by time
    time = "time"
    time_full = "delta_time" # Time as np.datetime objects
    
    ##TODO: Fix problem with this not working when time_increment is set
    # to None
    if time_increment is None:
        
        start_time = np.min(ds_tropomi[time_full].values)
        end_time = np.max(ds_tropomi[time_full].values)
        
        dt = start_time + (end_time - start_time)/2.
        data_split = [(dt,ds_tropomi)]
    
    else:
        dt0 = ds_tropomi[time_full].isel({time:0},drop=True)

        data_dt = ds_tropomi.assign_coords(coords={time_full:dt0})
        data_dt = data_dt.swap_dims({"scanline":time_full})
        data_split = data_dt.resample({time_full:time_increment})
    
    data_regridded = None
    for dt,data_t in data_split:
        
        if time_increment is not None:
            data_t = data_t.swap_dims({time_full:"scanline"})
        
        regridded = regrid_orbit(data_t,output_lat,output_lon,
                                 method=method)

        if regridded is not None:
            regridded[time] = [dt]
            if data_regridded is None:
                data_regridded = regridded #.expand_dims(time,axis=0)
            else:
                data_regridded = xr.concat([data_regridded,regridded],dim=time)
#        else:
#            print(f"No tropomi data within output grid in file {os.path.split(filename)[-1]} for datetime: {dt}")

    data_regridded.attrs["dlat"] = dlat
    data_regridded.attrs["dlon"] = dlon

    ##TODO: Here or somewhere else: write something to remove all the weight
    # file which get created after regridding is finished.
    # conservative_*.nc
    # conservative_normed_*.nc

    return data_regridded

def make_time_unique(ds,name="time"):
    '''
    For all times repeated within a dataset add small (nanosecond) increments
    to make these values unique.
    
    Args:
        ds (xarray.Dataset) :
            Dataset containing data variable for time (specified by name).
        name (str, optional) :
            Name of time data variable within dataset.
            Default = "time"
    
    Returns:
        xarray.Dataset :
            Original dataset with repeated time values modified.
    '''
    
    time = ds[name]
    time_group = time.groupby(name)
    
    time_unique = None
    for t,tg in time_group:
        repeats = len(tg)
        time_add = np.arange(0,repeats,step=1,dtype=int)
        # TODO: We are assuming this is in units of ns so may want to make that explicit
        time_update = tg + time_add
        
        dim = tg.dims[0]
        
        if time_unique is None:
            time_unique = time_update
        else:
            time_unique = xr.concat([time_unique,time_update],dim=dim)
    
    ds = ds.assign_coords({name:time_unique})
    
    ds.attrs["time_updated"] = "Repeated time values for each grid have had small arbitrary time increments added."

    return ds
                
   

def unravel_grid(ds_grid):
    '''
    Unravel regridded data file to create a timeseries with 1D time,lat,lon 
    values. Expects input data to be on a (time, lat, lon) 3D grid.
    
    Note this will also make time unique by adding small increments to
    any repeated time values.
    
    Args:
        ds_grid (xarray.Dataset) :
            Dataset containing data variables on a (time,lat,lon) 3D grid.
    
    Returns:
        xarray.Dataset :
            Original dataset flattened onto a 1D grid with time, lat, lon 
            values for each piece of data.
    '''
    
    temp_dim = "z"
    
    time = "time"
    latlon = ("lat","lon")
    
    # Stack time, latitude and longitude dimensions
    stack_dims = (time,) + latlon
    ds_flatten = ds_grid.stack({temp_dim:stack_dims})
    
    # Reset the temporary dimension and make lat, lon into data variables
    ds_flatten = ds_flatten.reset_index(temp_dim) # Remove MultiIndex
    ds_flatten = ds_flatten.reset_coords(names=latlon)
    
    # Remove nan values on the flattened dimension
    ds_flatten = ds_flatten.dropna(temp_dim)
    
    # Add small increments to repeated time values
    ds_flatten = make_time_unique(ds_flatten,name=time)
     
    # Make time the main dimension and co-ordinate and make sure this 
    # is the first dimension
    ds_flatten = ds_flatten.swap_dims({temp_dim:time})
    ds_flatten = ds_flatten.transpose(time,...)
    if temp_dim in ds_flatten.dims:
        ds_flatten = ds_flatten.drop(temp_dim)
    
    return ds_flatten

                # filename = site+"_"+date.replace('-','')
                # if max_points:
                #     if len(wh_date) > max_points:
                #         letter_split = chr(ord("A")+int(old_div(ID,max_points)))
                #         filename = "{}-{}".format(filename,letter_split)
                #     ID_1 = ID%max_points
                # else:
                #     ID_1 = ID
                # filename += '.csv'
                # filename = os.path.join(name_directory,filename)

    
def create_labels(number_of_points,max_level):
    '''
    Make ID_Level labels for inclusion in output NAME file.
    
    These are of of the form:
        N...LL
        e.g.
            100101 - points 1001, level 01
    
    N - will be equal to the number of digits for number_of_points
    L - always 2 digits for the level e.g. 01-17
    
    Args:
        number_of_points (int) :
            Total number of tropomi measurements
        max_level (int) :
            Number of pressure levels
    
    Returns:
        list :
            ID_Level labels for NAME output
    '''
    
    n_digits = len(str(number_of_points))
    
    labels = []
    for n in range(1,number_of_points+1):        
        point = f"{n:0{n_digits}d}"
        full_label = [point+f"{lev:02d}" for lev in range(1,max_level+1)]
        labels.extend(full_label)
    
    return labels
        

def write_tropomi_NAME(ds,site,max_level=None,max_points=50,
                       name_directory=name_csv_directory):
    '''
    Write NAME input files related to 1 day.
    For ds input expect a dataset with 1D time,lat,lon values.
    
    Args:
        ds (xarray.Dataset) :
            Dataset containing positions to write to NAME input files.
        site (str) :
            Name for the tropomi data selection e.g. TROPOMI-BRAZIL
        max_level (int/None, optional) :
            Maximum level to run to within NAME. Will only write 
            pressure levels up to and including this level if 
            specified.
            Default = None # All levels written to file
        max_points (int/None, optional):
            Maximum number of points to include in one NAME input
            file. If specified and number of points excees this number
            then multiple files will be created labelled as e.g. -A, =B
            etc. (see define_name_filenames() for more details).
        name_directory (str, optional) :
            Top level directory to write files for NAME. Full path will 
            be based on the "network" related to "site" 
            Default defined at the top of this module.            
         
    
    ##TODO: May want to generalise this to allow files to be written
    out for a longer time period rather than one day at a time.    
    '''
    # Note column mapping here is out:in values rather than in:out 
    # because some columns in output do not map to input data variables.
    
    time = "time"
    latlon = ("lat","lon")
    pressure = "pressure_levels"
    
    col_mapping = OrderedDict([('ID_Level',None),
                                ('Time','time'),
                                ('x','lon'),
                                ('dx',None),
                                ('y','lat'),
                                ('dy',None),
                                ('z',None),
                                ('dz',None)])    
    
    col_order = list(col_mapping.keys())[1:]
    
    # Create data frame for all data and then group?
    # How will NAME react to small time increment differences? Or should
    # we just cut off time for printing out to NAME file? And how will this
    # align with the footprints? With the current code this matches on time...
    
    #dx, dy will be based on whatever output grid was chosen so should be constant
    #dz can be based on the pressure_levels as these have layer_bounds rather than
    #just levels like with gosat.
    # - i and i+1 in pressure_levels
    # - Do we need to check these against the NAME surface again?

    if site is None:
        site = 'global'
    
    date = ds["time"].values[0]
    
    network = gosat_fn.find_network(site)[0] # Using first site as default

    output_directory = os.path.join(name_directory,network)
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    
    number_of_points = ds.dims[time]

    dpressure = np.abs(ds[pressure].diff(dim="layer_bound",label="lower"))
    pressure = ds[pressure].isel(layer_bound=slice(0,-1)) + dpressure/2.
    
    if max_level:
        dpressure = dpressure.isel({"layer_bound":slice(0,max_level+1)})
        pressure = pressure.isel({"layer_bound":slice(0,max_level+1)})
    else:
        max_level = dpressure["layer_bound"].size
    
    layer_bound = pressure["layer_bound"]
    
    ds_name = xr.Dataset()
    for new_name,name in col_mapping.items():
        if name is not None:
            dv = ds[name]
            dv = dv.expand_dims(dim={"layer_bound":layer_bound},axis=-1)
            ds_name = ds_name.assign({new_name:dv})
    
    dlat = xr.DataArray(np.array(ds.attrs["dlat"]))
    dlon = xr.DataArray(np.array(ds.attrs["dlon"])) 
    
    ## Can't do this because of missing points
    ##dlat = np.mean(ds[latlon[0]].diff(dim=time))#.values)
    ##dlon = np.mean(ds[latlon[1]].diff(dim=time))#.values)

    dlat = dlat.expand_dims(dim={"layer_bound":ds_name["layer_bound"]},axis=-1)
    dlon = dlon.expand_dims(dim={"layer_bound":ds_name["layer_bound"]},axis=-1)

    ds_name = ds_name.assign({"dx":dlon})
    ds_name = ds_name.assign({"dy":dlat})

    ds_name = ds_name.assign({"z":pressure})
    ds_name = ds_name.assign({"dz":dpressure})
    
    temp_dim = "stack"
    ds_name = ds_name.stack({temp_dim:("time","layer_bound")})
    ds_name = ds_name.reset_index(temp_dim)
    
    id_label = create_labels(number_of_points,max_level)
    ds_name = ds_name.assign_coords({"ID_Level":(temp_dim,id_label)})
    ds_name = ds_name.swap_dims({temp_dim:"ID_Level"})
    ds_name = ds_name.drop(labels=["time","layer_bound"])
    
    df = ds_name.to_dataframe()
    df = df[col_order]
    
    name_filenames = \
        gosat_fn.define_name_filenames(name_directory, site, date,
                                       number_of_points=number_of_points,
                                       max_points=max_points)
    
    if max_points:
        lines_per_file = max_points*max_level
    else:
        lines_per_file = number_of_points*max_level
    
    ##TODO: Check start:end is including all the correct lines, since this 
    # selection has been derived
    for i, name_fname in enumerate(name_filenames):
        
        start = i*lines_per_file
        end = (i+1)*lines_per_file

        print(f"Writing to {name_fname}")
        df.iloc[start:end].to_csv(name_fname)
  

def write_tropomi_output(ds,site,date,
                         #species="ch4",
                         output_directory=obs_directory):
    '''
    Write tropomi observation file for one day.
    For obs data would expect dataset with 1D time,lat,lon values.
    
    Args:
        ds (xarray.Dataset) :
            TROPOMI dataset to write to file. Should contain up to one 
            day's worth of data at present.
        site (str) :
            Name for the tropomi data selection e.g. TROPOMI-BRAZIL
        date (str) :
            Date related to this tropomi data.
        output_directory (str, optional) :
            Top-level directory to write output. Full path will be based 
            on the "network" related to "site".
    
    Returns:
        None
        
        Writes ds out as netcdf file.
    
    ##TODO: Also want to filter to max_level input? Or just write as an 
    attribute? Need to pass as a variable perhaps.
    
    ##TODO: May want to generalise this to allow files to be written
    out for a longer time period rather than one day at a time.
    '''
    
    instrument = "tropomi"
    satellite = "sentinel5p"
    
    ## TODO: Generalise to allow for different species than methane
    species = "ch4"
    
    network = gosat_fn.find_network(site)[0] # Using first site as default
    
    output_filename = \
        gosat_fn.define_obs_filename(output_directory,network,instrument,
                                     satellite,date,species,inlet=None,num=None)
    
    ##TODO: Add general attributes - time of creation, who created the
    # file etc.    
    
    print(f"Writing to file: {output_filename}")
    ds.to_netcdf(output_filename)
        


    

def tropomi_process(start_date,end_date,lat_bounds,lon_bounds,coord_bin,
                    max_level,site=None,max_points=50,
                    #species="ch4",
                    time_increment="10s",
                    #apply_qa_filter=True,
                    regrid_method="conservative",
                    input_directory=input_directory,
                    write_name=False,name_directory=name_csv_directory,
                    output_directory=obs_directory,verbose=False):
    '''
    Process tropomi data within date range including re-gridding onto 
    a regular latitute-longitude grid.
    
    Args:
        start_date (str) :
            YYYY-MM-DD
        end_date (str) :
            Up to but not including this date. YYYY-MM-DD
        lat_bounds (list) :
            Latitude lower and upper bounds to include in output (2-item list).
            Upper bound included.
        lon_bounds (list) :
            Longitude lower and upper bounds to include in output (2-item list).
            Upper bound included.
        coord_bin (float/list) :
            Bins to use for latitude and longitude dimensions.
            Specify one value to use for both or a 2-item list to
            use different values.
        time_increment (str, optional) :
            Time window to group tropomi points.
            Default = "10s" # 10 seconds
        regrid_method (str, optional) :
            Regridding method to use. Options include:
                "conservative"
                "conservative_normed"
            "conservative_normed" will retain more points but
            requires the installation of the "masking" branch for xesmf.
            Ask Rachel/Daniel for more details.
        input_directory (str, optional) :
            Top level directory for reading TROPOMI data.
            Full path will be based on "species".
        write_name (bool, optional) : 
            Write output .csv file for input into NAME
        name_directory (str, optional) : 
            Top level directory to write files for NAME. Full path will 
            be based on the "network" related to "site" 
            Default defined at the top of this module.
        output_directory (str, optional) :
            Top level directory to write output. Full path will be based 
            on the "network" related to "site"
        verbose (bool, optional) :
            Print more descriptive output as function executes.
            Default = False
    
    Returns:
        None
        
        Writes to netcdf file
        
        If write_name=True:
            Writes to csv NAME input file.
    
    
    ##TODO: Add (back in) filter for quality flag - may need to make
    # delta_time coordinate so we don't get NaT values
    ##TODO: Add species input and filter to choose correct search_species string

    '''
    timing_1 = timing_module.time()
    
    print(f"Processing tropomi files for date range: {start_date} - {end_date}")
    
    species = "ch4"
    ##TODO: Update input paths, currently hardwired to RT local BP1 area
    #base_directory = Path("/home/rt17603/shared/obs_raw/TROPOMI/SOUTHAMERICA/")
    #base_output_directory = Path("/home/rt17603/work/TROPOMI_processed/SOUTHAMERICA/")
    if species == "ch4":
        directory = input_directory / "CH4"
        #output_directory = base_output_directory / "CH4"
        search_species = "CH4____"
    analysis_mode = "OFFL"
    search_str = f"S5P_{analysis_mode}_L2__{search_species}*.nc"

    time = "time"
    #latlon=("latitude", "longitude")
    #pos_coords=("scanline", "ground_pixel")

    dates = np.arange(start_date,end_date,dtype=np.datetime64).astype(str)
    for date in dates:
        start = date
        end = (np.datetime64(date) + np.timedelta64(1,'D')).astype(str)
        
        filenames = gosat_fn.extract_files(directory,search_str,start,end)
        
        if len(filenames) == 0:
            current_analysis_mode = analysis_mode
            analysis_mode = "RPRO"
            if verbose:
                print(f"Couldn't find {current_analysis_mode} files, looking for {analysis_mode} (reprocessed)")
            search_str_alt = search_str.replace(current_analysis_mode,analysis_mode)
            filenames = gosat_fn.extract_files(directory,search_str_alt,start,end)
        
        if len(filenames) == 0:
            print(f"No files found for date: {start}")
            continue
        
        data_regridded = None
        for i,filename in enumerate(filenames):
            
            print(f"Processing input file: {filename}")
            ## TODO: Currently seems to be relying on edge_coords being
            # set to bounds - need to to update this and follow it through
            # the chain.
            data = preProcessFile(filename,edge_coords="bounds")
            
            ## TODO: Add in quality filter so this doesn't introduce
            # NaT objects in delta_time
            #if apply_qa_filter:
            #    data = data.where(data["qa_pass"])
            
            # Regrid data for each file
            regridded = regrid_tropomi(data,lat_bounds,lon_bounds,coord_bin,
                                       method=regrid_method,
                                       time_increment=time_increment)
            
            if regridded is not None:
                if data_regridded is None:
                    data_regridded = regridded
                else:
                    data_regridded = xr.concat([data_regridded,regridded],dim=time)
        
        ##TODO: Filter any points where some specific values are not defined
        # e.g. methane, pressure levels
        
        data_timeseries = unravel_grid(data_regridded)
        #TODO: Temporary line for now (update/remove)
        data_output = data_timeseries
        
        short_filenames = [os.path.split(fname)[1] for fname in filenames]
        short_filenames = ','.join(short_filenames)
        
        ## TODO: Add back in once quality filter is being used again
        #if apply_qa_filter:
        #    data_output.attrs["qa_filter"] = "Quality filter of > 0.5 applied"
        data_output.attrs["input_filename"] = short_filenames
        data_output.attrs["analysis_mode"] = f"Input files are from {analysis_mode} mode of analysis."
        
        write_tropomi_output(site,date,species=species,
                             output_directory=output_directory)
        
        # output_filename = output_directory / f"tropomi_sentinel5p_{start}_{species}-column.nc"
        # print(f"Writing to file: {output_filename}")
        # data_output.to_netcdf(output_filename)

        if write_name:
            ##TODO: Work out additional variables that need passing to this function
            write_tropomi_NAME(data_output,site=site,max_level=max_level,
                                max_points=max_points)

        print(f"\nTime to execute: {timing_module.time() - timing_1}\n")
        
    
    