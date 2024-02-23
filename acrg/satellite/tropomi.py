#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 08:30:03 2019

@author: al18242
@author: rt17603
@author: ky20893

This module provides a set of functions for processing and using TROPOMI data.

Within this module the main summary functions used for processing is:
    * tropomi_process(directory, ...) - process a directory of files

This applies the following processes to the input files:
 - Extract relevant details from the grouped netcdf TROPOMI data file
 - Filter by the recommended quality conditions (combined filter > 0.5). (optional)
 - Regrid the TROPOMI data to a regular 3D (time, lat, lon) grid (defined by lat_bounds,lon_bounds,coord_bin)
    - Note: this can include each time point along the scanline direction but by default
      this will be further grouped along the time axis to 1 minute increments (`time_increment`)
 - Unravel the created 3D (time, lat, lon) grid to create a 1D (time) array of time points each with an associated lat-lon value
   - Adds small time increments to make this unique for repeated time points (ground_pixel dimension)
 - Create netcdf output file (*format needs finalising*)
 - Creates file for using with NAME model (optional) (defined by `write_name` flag)
"""

import os
import glob
import socket
import numpy as np
import xarray as xr
import time as timing_module
from pathlib import Path
from collections import OrderedDict
import acrg.grid.regrid_xesmf as regrid
import acrg.satellite.gosat as gosat_fn

#from acrg_satellite.gosat import latlon_filter,extract_files,find_network

try:
    from acrg.config.paths import paths
except ImportError:
    data_path = os.getenv("DATA_PATH")
    acrg_path = os.getenv("ACRG_PATH")
else:
    data_path = paths.data
    acrg_path = paths.acrg

home = Path(os.getenv("HOME"))
if not acrg_path:
    acrg_path = os.path.join(home,"acrg")

if not data_path:
    server_name = socket.gethostname()
    if "bp1" in server_name:
        data_path = "/work/chxmr/shared/"
    elif "bc4" in server_name:
        data_path = "/mnt/storage/private/acrg/ACRG_Repository/data/"
    elif "snowy" in server_name:
        data_path = "/data/shared/"
    else:
        data_path = home
        print(f"Unable to infer data_path - setting to home directory: {data_path}")

##TODO: Replace this more general directory (not SOUTHAMERICA) once this is sorted
input_directory = os.path.join(data_path, "obs_raw/TROPOMI/ARCTIC/")

name_csv_directory = os.path.join(home,"NAME_files") # Where to write output NAME csv files
obs_directory = os.path.join(data_path,'obs/') # Where to write output nc files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure/")

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

def preProcessFile(filename,add_corners=False):
    '''Combine important information from different netCDF groups into a single data frame
    
    input
        filename - string of full filename for TROPOMI file
        add_corners (bool, optional) :
            Explicitly add the cooridinates for the corners to 
            the dataset. This will add the coordinates:
                "lat_corners", "lon_corners"
            Both with dimensions:
                ("time", "scanline_c", "ground_pixel_c")
            where the dimension of ground_pixel_c and scanline_c 
            are ground_pixel + 1 and scanline + 1 respectively.
            This is useful for plotting.
        
    return
        tropomi_data - dataframe containing mixing ratio, quality flags, grid data, and column data
                        the "qa_pass" flag variable is True for good measurements and False for bad measurements
    
    ## TODO: Add species input and allow this to choose which parameters
    # to extract based on the species.
    '''
    tropomi_data = xr.open_dataset(filename, group="PRODUCT")
    tropomi_data_aux = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS")
    tropomi_data_geo = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS")
    tropomi_data_input = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA")
    
    # If delta_time is a timedelta64 and not a datetime64 object
    # make sure the reference time is added to make this 
    # into a datetime64
    if tropomi_data["delta_time"].values.dtype == 'timedelta64[ns]':
        tropomi_data["delta_time"] = tropomi_data["delta_time"] + tropomi_data["time"]
    
    #calculate boundaries of pressure boundaries from surface pressure and constant intervals
    nlevel = tropomi_data["level"].shape[0]
    pressure_bounds = (np.expand_dims(tropomi_data_input.surface_pressure,axis=3) - \
                     np.expand_dims(tropomi_data_input.pressure_interval,axis=3) * \
                     np.reshape(np.arange(0,nlevel),newshape=(1,1,1,-1)))
    
    # Tropomi pressure is in Pa so /100 to convert to hPa, to be the same as in GOSAT obs files.
    tropomi_data['pressure_bounds'] = (["time", "scanline", "ground_pixel", "layer_bound"], pressure_bounds/100)
    tropomi_data['pressure_bounds'].attrs["units"] = "hPa"

    #calculate mid point of pressure layers from surface pressure and constant intervals
    nlayer = tropomi_data["layer"].shape[0]
    surface_layer_mid_pressure = tropomi_data_input.surface_pressure - tropomi_data_input.pressure_interval/2.
    pressure_data = (np.expand_dims(surface_layer_mid_pressure,axis=3) - \
                     np.expand_dims(tropomi_data_input.pressure_interval,axis=3) * \
                     np.reshape(np.arange(0,nlayer),newshape=(1,1,1,-1))) 

    # Tropomi pressure is in Pa so /100 to convert to hPa, to be the same as in GOSAT obs files.
    tropomi_data['pressure_levels'] = (["time", "scanline", "ground_pixel", "layer"], pressure_data/100)    
    tropomi_data['pressure_levels'].attrs["units"] = "hPa"
    tropomi_data['pressure_levels'].attrs['short_description'] = "Vertical altitude coordinate in pressure units as used for averaging kernels"
    
    #tropomi_data['column_averaging_kernel'] = tropomi_data_aux.column_averaging_kernel
    #tropomi_data['methane_profile_apriori']= tropomi_data_input.methane_profile_apriori
 
    ## Add additional parameters using same naming scheme as GOSAT
    # This is to be consistent with other processes within the repository

    # Check the averaging kernel is increasing with altitude across the layer dimension - this is reversed in the raw S5P_RPRO files 
    ### A much better check should be introduced at some point ###
    # If not, reverse the averaging kernel
    if np.nanmean(tropomi_data_aux.column_averaging_kernel, axis = (0,1,2))[0] < np.nanmean(tropomi_data_aux.column_averaging_kernel, axis = (0,1,2))[-1]:
        # print("Reversing the averaging kernel as this is reversed in the raw TROPOMI files.")
        tropomi_data['xch4_averaging_kernel'] = tropomi_data_aux.column_averaging_kernel[:,:,:,::-1]
    else:
        print("WARNING - Check raw input files thoroughly. Averaging kernel is in the correct order for use in the retrieval, which is not expected.")
        tropomi_data['xch4_averaging_kernel'] = tropomi_data_aux.column_averaging_kernel
    
    # Calculating the pressure weights from the dry air subcolumns
    tropomi_data['dry_air_subcolumns'] = tropomi_data_input.dry_air_subcolumns
    tropomi_data['pressure_weights'] = tropomi_data["dry_air_subcolumns"]/tropomi_data["dry_air_subcolumns"].sum(dim="layer")
    
    # Updating the a priori profile to match to GOSAT equations so we can
    # apply the pressure weights in the same way. Multiplying by 1e9 to be in units of ppb.

    # Check the ch4_profile_apriori is increasing with altitude across the layer dimension - this is reversed in the raw S5P_RPRO files 
    ### A much better check should be introduced at some point ###
    # If not, reverse the ch4_profile_apriori
    if np.nanmean(tropomi_data_aux.column_averaging_kernel, axis = (0,1,2))[0] < np.nanmean(tropomi_data_aux.column_averaging_kernel, axis = (0,1,2))[-1]:
        # print("Reversing the prior profile as this is reversed in the raw TROPOMI files.")
        tropomi_data['ch4_profile_apriori']= 1e9*tropomi_data_input.methane_profile_apriori[:,:,:,::-1]/tropomi_data["dry_air_subcolumns"]
    else:
        print("WARNING - Check raw input files thoroughly. Prior profile is in the correct order for use in the retrieval, which is not expected.")
        tropomi_data['ch4_profile_apriori']= 1e9*tropomi_data_input.methane_profile_apriori/tropomi_data["dry_air_subcolumns"]

    tropomi_data['ch4_profile_apriori'].attrs["units"] = "ppb"
    tropomi_data['ch4_profile_apriori'].attrs['short_description'] = "A priori mole fraction profile of atmospheric CH4."
    
    tropomi_data = tropomi_data.assign_coords({'latitude_bounds':tropomi_data_geo['latitude_bounds']})
    tropomi_data = tropomi_data.assign_coords({'longitude_bounds':tropomi_data_geo['longitude_bounds']})
    
    if add_corners:
        #calculate the meshgrid of corner points in the scan grid for further processing
        #lat, lon = cornerGrid(tropomi_data_geo)
        lat, lon = cornerGrid(tropomi_data)
        #tropomi_data['lat_corners'] = (["time", "scanline_c", "ground_pixel_c"], np.expand_dims(lat,0))
        #tropomi_data['lon_corners'] = (["time", "scanline_c", "ground_pixel_c"],np.expand_dims(lon,0))
        tropomi_data = tropomi_data.assign_coords({'lat_corners':(["time", "scanline_c", "ground_pixel_c"], np.expand_dims(lat,0))})
        tropomi_data = tropomi_data.assign_coords({'lon_corners':(["time", "scanline_c", "ground_pixel_c"],np.expand_dims(lon,0))})
        
    #flag quality to remove outputs not flagged as 'success' and use recommended acceptable qa_values
    # Note: it may be the case that any value that passes the second condition necessarily passes the first
    quality_flag = (tropomi_data_aux.processing_quality_flags.values.astype(int) & 0xFF == 0) & \
                    (tropomi_data.qa_value.values > 0.5)
    
    tropomi_data["qa_pass"] = (["time", "scanline", "ground_pixel"], quality_flag)

    # Rename methane mixing ratio variables to match GOSAT files
    tropomi_data = tropomi_data.rename({
        "methane_mixing_ratio":"xch4",
        "methane_mixing_ratio_precision":"xch4_uncertainty",
        "methane_mixing_ratio_bias_corrected":"xch4_bias_corrected"})
    
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
    
    Note: if the DataArray contains additional dimensions beyond
    those specified by pos_coords, the first value of that dimension
    is selected to create this mask.
    
    Args:
        da (xarray.DataArray) :
            Data values for one data variable (e.g. ds["xch4"]).
        pos_coords (tuple, optional) :
            The coordinate names describing the position within da.
            Default = ("scanline", "ground_pixel")
    
    Returns: numpy.array of boolean values
    '''
    dims = da.dims
    # Make a copy of the DataArray to make sure the underlying values
    # are not changed.
    da_copy = da.copy(deep=True)
    
    for dim in dims:
        if dim not in pos_coords:
            # For additional dimensions take the first value
            da_copy = da_copy.isel({dim:0})
    
    input_values = da_copy.rename({pos_coords[0]:'x',pos_coords[1]:'y'})
    input_values = np.isfinite(input_values)
    
    return input_values

def remove_weight_files(method, path=None):
    '''
    Remove weight files created by regrid method. Be careful not to
    do this if multiple regridding processes are being run in 
    parallel.
    
    Expect weight files have been created within the same folder as 
    this module with file naming convention:
        METHOD_NUMxNUM_NUMxNUM.nc
        e.g.
        conservative_normed_577x202_706x494.nc
    
    Args:
        method (str):
            Method name for regrid. Usually one of "conservative" or
            "conservative_normed".
        path (str/None, optional):
            Path to search for files. If not specified:
                 - looks within the same folder as this module file.
                 - if that can't be identified then it looks within 
                   ACRG_PATH/acrg_satellite folder.    
            
    Returns:
        None
        
        Files matching the path and search string are deleted.
    '''
    if path is None:
        try:
            filename = os.path.realpath(__file__)
            path = os.path.split(filename)[0]
        except NameError:
            path = os.path.join(acrg_path,"acrg_satellite")

    # Find all files within this folder than are netcdf files that
    # loosely match the method name and expected format (first pass)
    weight_files = glob.glob(os.path.join(path,f"{method}*.nc"))

    # Extra check to make sure the correct files are being selected
    # for deletion (using regular expressions).
    # e.g. don't delete "conservative_normed" files if running "conservative" method
    import re
    re_search = re.compile(os.path.join(path,f"{method}_\d+x\d+_\d+x\d+[.]nc"))    
    print(f"Removing regrid weighting files for {method} method from: {path}")
    for file in weight_files:
        if re_search.match(file):
            os.remove(file)


def regrid_da(da_tropomi,output_lat,output_lon,ds_tropomi_geo,
              method="conservative_normed",latlon=["latitude","longitude"],
              pos_coords=("ground_pixel","scanline"),
              set_nan=True,reuse_weights=False):
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
            This dataset needs to contain "latitude_bounds" and "longitude_bounds" data variables.
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
        reuse_weights (bool, optional) :
            Reuse previous weights created over the same grid area. This speeds
            up computation but should make sure this exactly matches the previous
            grid shape and any nan values.
            Default = False
    '''

    input_lat = da_tropomi[latlon[0]][0,:,:].values
    input_lon = da_tropomi[latlon[1]][0,:,:].values
    
    input_lat_b, input_lon_b = cornerGrid(ds_tropomi_geo)

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
    
    regridded = regridded.assign_coords(**{"x":output_lat,"y":output_lon})
    regridded = regridded.rename({"x":"lat","y":"lon"})
    
    # preserve attributes from input da
    regridded.attrs = da_tropomi.attrs

    return regridded


def regrid_subset(ds_tropomi,output_lat,output_lon,names=None,ds_tropomi_geo=None,
                 method="conservative_normed",latlon=("latitude", "longitude"),
                 pos_coords=("scanline", "ground_pixel"),
                 exclude=["time_utc","qa_pass","qa_value"],set_nan=True,
                 filter_latlon=True,clean_up_weights=False,
                 verbose=False):
    '''
    The regrid_subset function takes a subset from one tropomi orbit and 
    regrids this onto a new output grid.
    
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
        exclude (list, optional) :
            Data variables within the dataset to not include in the
            regridded output.
            Default = ["time_utc","qa_pass","qa_value"]
        set_nan (bool, optional) :
            After regridding, explicitly set values of 0 to np.nan.
            Default = True
        filter_latlon (bool, optional) :
            Filter the input data to only include the range covered by the output grid.
            May make computation quicker.
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
    
    #excluded_dv = ["time_utc"]
    names = [name for name in names if name not in exclude]
    
    # For each data variable apply the regridding algorithm
    for i,name in enumerate(names):
        
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
                              set_nan=set_nan,
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
    
    ds_tropomi_regridded.attrs["regridded"] = f"Regridded to output grid of lat = {output_lat[0]:.3f} - {output_lat[-1]:.3f}, lon = {output_lon[0]:.3f} - {output_lon[-1]:.3f}"
    
    if clean_up_weights:
        remove_weight_files(method)
    
    return ds_tropomi_regridded


def regrid_orbit(ds_tropomi,lat_bounds,lon_bounds,coord_bin,
                 method="conservative_normed",time_increment="1min",
                 exclude=["time_utc","qa_pass","qa_value"],
                 clean_up_weights=True):
    '''
    The regrid_orbit function regrids data for one tropomi orbit. This can 
    either be for the whole orbit at once or split into time windows based
    on the time_increment input.    
    
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
                - "10s" (10 seconds)
                - "1min" (1 minute)
                - "1h" (1 hour)
            Should always be a day or less.
            Can set to None to not split the data by time.
            
            Default = "1min"
        
        exclude (list, optional) :
            Data variables within the dataset to not include in the
            regridded output.
            
            Default = ["time_utc","qa_pass",""qa_value"]
            
        
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
        if time_full not in ds_tropomi.coords:
            dt0 = ds_tropomi[time_full].isel({time:0},drop=True)
            ds_tropomi = ds_tropomi.assign_coords(coords={time_full:dt0})
        data_dt = ds_tropomi.swap_dims({"scanline":time_full})
        data_split = data_dt.resample({time_full:time_increment})
    
    data_regridded = None
    for dt,data_t in data_split:
        
        if time_increment is not None:
            data_t = data_t.swap_dims({time_full:"scanline"})
        
        regridded = regrid_subset(data_t,output_lat,output_lon,
                                 method=method,exclude=exclude)

        if regridded is not None:
            regridded[time] = [dt]
            if data_regridded is None:
                data_regridded = regridded #.expand_dims(time,axis=0)
            else:
                data_regridded = xr.concat([data_regridded,regridded],dim=time)
        # else:
        #     print(f"No tropomi data within output grid in file for datetime: {dt}")
    
    if data_regridded is not None:
        data_regridded.attrs["dlat"] = dlat
        data_regridded.attrs["dlon"] = dlon
        
        if time_increment is None:
            comment = "Time has been resampled (to allow for regridding) to the midpoint for each orbit"
        else:
            comment = f"Time has been resampled (to allow for regridding) into bins of {time_increment}"
        
        data_regridded.attrs["time_resample"] = comment

    if clean_up_weights:
        remove_weight_files(method)

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
    time_group = np.unique(time)

    time_unique = None
    for t,tg in enumerate(time_group):
        # Adds a different number of nanoseconds to each repeat time to
        # make the times unique
        repeats = len(np.where(time == tg)[0])
        time_add = np.arange(0,repeats,step=1,dtype=int)

        if "[ns]" in str(time[t].dtype):
            # At the moment, needs datetime in ns
            time_update = tg + time_add
        else:
            raise ValueError("Do not recognise time units as nanoseconds")
            
        dim = time.dims[0]
        time_update_xarray = xr.DataArray(time_update, dims = dim)

        if time_unique is None:
            time_unique = time_update_xarray
        else:
            time_unique = xr.concat([time_unique,time_update_xarray],dim=dim)

    # print(f"Time unique unique values? : {len(np.unique(time_unique))}")
    # print(f"Original time unique values? : {len(np.unique(time))}")
    
    ds = ds.assign_coords({name:time_unique})
    
    ds.attrs["time_updated"] = "Repeated time values for each grid have had small arbitrary time increments added."

    return ds
                
   

def unravel_grid(ds_grid):
    '''
    Unravel regridded data file to create a timeseries with 1D time,lat,lon 
    values. Expects input data to be on a (time, lat, lon) 3D grid.
    
    Note this will also make time unique by adding small increments (ns) to
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
        point = f"{n:0{n_digits}d}" # Pad with zeros
        full_label = [point+f"{lev:02d}" for lev in range(1,max_level+1)]
        labels.extend(full_label)
    
    return labels

def use_NAME_surface_pressure(ds, pressure_domain, 
                              pressure_base_dir=name_pressure_directory,
                              pressure_NAME=None, p_column="pressure_levels",
                              p_coord="layer_bound",
                              columns=["latitude","longitude","time"],
                              day_template=True,max_days=31):
    '''
    Replace the surface pressure from within the TROPOMI file with the
    surface pressure from NAME for the NAME input (release point) files.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing pressure values to match the NAME surface pressure values to. Within the 
            dataset the pressure values should be related to latitude, longitude and time values.
        pressure_domain (str) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_NAME (np.array, optional) : 
            If pressure from NAME run has already been extracted, this can be specified explicitly to save 
            computing time.
            If not specified, columns from ds and pressure_dir will be used to extract matching pressure values.
        p_column (str, optional) : 
            Name for the pressure_levels data (str). Default = "pressure_levels"
        p_coord (str, optional) : 
            Name for the layer co-ordinate within pressure_levels data (str). 
            Default = "layer_bound"
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude, longitude and time values 
            (3 item list). Default = ["latitude","longitude","time"]
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        max_days (int, optional) : 
            Maximum number of days from time within ds to search for the relevant pressure data. 
            Default = 31 (days).

    Returns:
        xarray.Dataset:
            Original dataset with surface pressure levels replaced.
    
    '''
    
    if pressure_NAME is None:
        if pressure_domain is not None:
            pressure_NAME = \
                gosat_fn.name_pressure_match(ds,
                                             pressure_domain=pressure_domain,
                                             pressure_base_dir=pressure_base_dir,
                                             columns=columns,
                                             day_template=day_template,
                                             max_days=max_days)
        else:
            raise Exception("Pressure_domain must be specified if pressure values need to be \
                            extracted to use define_pressure_levels function.")
    # pressure from name_pressure_match() function is in hPa and needs to be
    # converted to Pa
    name_pressure_convert = 100.
    #print(pressure_NAME)
    pressure_NAME *= name_pressure_convert
    #print(pressure_NAME)
    
    # Reassign surface pressure to lower bound of layers
    surface_index = 0
    dims = ds[p_column].dims
    axis = dims.index(p_coord)
    if axis == len(dims)-1:
        ds[p_column][...,surface_index] = pressure_NAME
        ds[p_column].attrs["name_surface_pressure"] = f"Surface pressure ({p_coord}={surface_index}) matched to NAME surface pressure"
        #ds["surface_pressure"] = f"Surface pressure ({p_coord}={surface_index}) matched to NAME surface pressure"
        ds["surface_pressure"].attrs["updated"] = f"Surface pressure ({p_coord}={surface_index}) matched to NAME surface pressure"
    
    return ds    
      

def write_tropomi_NAME(ds,site,max_level=None,max_points=50,
                       network=None,
                       overwrite=False,
                       use_name_pressure=False,
                       pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31,
                       pressure_day_template=True,
                       name_directory=name_csv_directory):
    '''
    Write NAME input files related to 1 day.
    For ds input expect a dataset with 1D time,lat,lon values.
    
    See: define_name_filenames() function for details of NAME csv filename.

    Example file format (GOSAT example):
    ID_Level,Time,x,dx,y,dy,z,dz
    00101,2018-02-01 15:11:48,-33.506783,0.094621,-3.650032,0.094429,99193.199539,3614.003754
    00102,2018-02-01 15:11:48,-33.506783,0.094621,-3.650032,0.094429,93822.554398,7127.278900
    [...]
    00117,2018-02-01 15:11:48,-33.506783,0.094621,-3.650032,0.094429,4750.000000,3500.000000
    00201,2018-02-01 15:11:53,-33.561993,0.094649,-3.906265,0.094429,99196.054268,3611.492538
    00202,2018-02-01 15:11:53,-33.561993,0.094649,-3.906265,0.094429,93828.875732,7122.863770
    00203,2018-02-01 15:11:53,-33.561993,0.094649,-3.906265,0.094429,86756.071472,7022.74475
    [...]
    
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
        network (str/None, optional) : 
            Which network is being considered e.g. "GOSAT/GOSAT-INDIA"
            If present, this will be extend the output path e.g. "/shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/"
        overwrite (bool, optional) :
           Whether to allow overwiting of existing files, if present.
           Default = False.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the surface pressure value for each data 
            point.
        
        See: use_NAME_surface_pressure() function for arguments used when
        use_name_pressure is True
            
        name_directory (str, optional) :
            Top level directory to write files for NAME. Full path will 
            be based on the "network" related to "site" 
            Default defined at the top of this module.            
    
    Returns:
        None
        
        Dataset is reformatted and written output to a csv file.
    
    ##TODO: May want to generalise this to allow files to be written
    out for a longer time period rather than one day at a time.    
    '''
    ## Define parameters which can be directly extracted from input Dataset
    # Note column mapping here is out:in values rather than in:out 
    # because some columns in output do not map to input data variables.

    time = "time"
    lat = "lat"
    lon = "lon"
    pressure_column = "pressure_bounds"
    
    col_mapping = OrderedDict([('ID_Level',None),
                                ('Time','time'),
                                ('x','lon'),
                                ('dx',None),
                                ('y','lat'),
                                ('dy',None),
                                ('z',None),
                                ('dz',None)])    
    
    # Define order of columns in output
    col_order = list(col_mapping.keys())[1:]
    
    if site is None:
        site = 'global'
        
    ##TODO: Fix problem with extracting date here - coming out as a 
    #timedelta object?
    date = ds[time].values[0]
    
    if network is None:
        network = gosat_fn.find_network(site)[0] # Using first site as default

    ## Create output directiry if not present (based on network)
    output_directory = os.path.join(name_directory,network)
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    
    number_of_points = ds.dims[time]

    ## Replace TROPOMI surface pressure with NAME surface pressure (if specified)
    if use_name_pressure:
        columns = [lat, lon, time]
        ds = use_NAME_surface_pressure(ds, columns=columns,
                                       p_column=pressure_column,
                                       pressure_base_dir=pressure_base_dir,
                                       pressure_domain=pressure_domain,
                                       max_days=pressure_max_days,
                                       day_template=pressure_day_template)

    ## Extract pressure values (for z and dz outputs), these are in units of hPa, so convert to Pa by *100
    dpressure = 100*np.abs(ds[pressure_column].diff(dim="layer_bound",label="lower"))
    pressure = 100*ds[pressure_column].isel(layer_bound=slice(0,-1)) - dpressure/2.
    
    # Limit to max_level (if specified)
    if max_level:
        dpressure = dpressure.isel({"layer_bound":slice(0,max_level+1)})
        pressure = pressure.isel({"layer_bound":slice(0,max_level+1)})
    else:
        max_level = dpressure["layer_bound"].size
    
    ## Start to construct output - ds_name
    # Need repeated time, x, y for each layer (z)
    # so expand these to repeat along the layer_bound dimension as well
    # dims = (time x layer_bound)
    layer_bound = pressure["layer_bound"]
    
    ds_name = xr.Dataset()
    for new_name,name in col_mapping.items():
        if name is not None:
            dv = ds[name]
            dv = dv.expand_dims(dim={"layer_bound":layer_bound},axis=-1)
            ds_name = ds_name.assign({new_name:dv})
    
    ## Create dlon and dlat inputs (dx, dy) based on stored values
    # - "dlat" and "dlon" attributes are created when tropomi_regridded is applied
    # These will be the same values for all x, y inputs.
    dlat = xr.DataArray(np.array(ds.attrs["dlat"]))
    dlon = xr.DataArray(np.array(ds.attrs["dlon"])) 
    
    ### Can't calculate this using lat, lon values because of missing points
    ##dlat = np.mean(ds[latlon[0]].diff(dim=time))#.values)
    ##dlon = np.mean(ds[latlon[1]].diff(dim=time))#.values)

    dlat = dlat.expand_dims(dim={"layer_bound":ds_name["layer_bound"]},axis=-1)
    dlon = dlon.expand_dims(dim={"layer_bound":ds_name["layer_bound"]},axis=-1)

    ## Add remaining values to output Dataset
    ds_name = ds_name.assign({"dx":dlon})
    ds_name = ds_name.assign({"dy":dlat})

    ds_name = ds_name.assign({"z":pressure})
    ds_name = ds_name.assign({"dz":dpressure})
    
    ## Unravel (flatten) time x layer_bound onto one dimension
    temp_dim = "stack"
    ds_name = ds_name.stack({temp_dim:("time","layer_bound")})
    ds_name = ds_name.reset_index(temp_dim)
    
    ## Create ID values of the appropriate form and add to output
    id_label = create_labels(number_of_points,max_level)
    ds_name = ds_name.assign_coords({"ID_Level":(temp_dim,id_label)})
    ds_name = ds_name.swap_dims({temp_dim:"ID_Level"})
    ds_name = ds_name.drop(labels=["time","layer_bound"])
    
    ## Cast output to a simple indexed DataFrame and make sure columns 
    # are in the expected order
    df = ds_name.to_dataframe()
    df = df[col_order]
    
    ## Write data to one or more output files
    # If max_points is specified, the number of files is 
    # number_of_points//max_points + 1
    # i.e. files containing a maximum number of time points each
    name_filenames = \
        gosat_fn.define_name_filenames(output_directory, site, date,
                                       number_of_points=number_of_points,
                                       max_points=max_points)
    
    if max_points:
        lines_per_file = max_points*max_level
    else:
        lines_per_file = number_of_points*max_level
    
    ##TODO: Write in check for overwriting an existing file

    # Write correct number of points from DataFrame to each output file
    for i, name_fname in enumerate(name_filenames):
        
        start = i*lines_per_file
        end = (i+1)*lines_per_file

        # print("start,end",start,end)

        if not overwrite:
            if os.path.isfile(name_fname):
                print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(name_fname))
            else:
                print('Writing to filename:',name_fname)    
                df.iloc[start:end].to_csv(name_fname)
        else:
            print('Writing to filename:',name_fname)    
            df.iloc[start:end].to_csv(name_fname)



def write_tropomi_output(ds,site,date,species="ch4",
                         network=None,
                         overwrite=False,
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
        species (str, optional) :
            Species name for using within the filename.
            Default = "ch4"
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
    
    ##TODO: Add checks for overwriting files
    '''
    
    instrument = "tropomi"
    satellite = "sentinel5p"
    
    ## TODO: Generalise to allow for different species than methane
    species = "ch4"
    
    if network is None:
        network = gosat_fn.find_network(site)[0] # Using first site as default
    
    full_output_directory = os.path.join(output_directory,network)
    if not os.path.isdir(full_output_directory):
        os.makedirs(full_output_directory)

    output_filename = \
        gosat_fn.define_obs_filename(full_output_directory,instrument,
                                     satellite,date,species,inlet="column",num=None)
    
    ##TODO: Add general attributes - time of creation, who created the
    # file etc.    
    
    # Write file 
    gosat_fn.write_netcdf(ds, output_filename, overwrite)


def define_tropomi_search_str(species="ch4",analysis_mode="OFFL"):
    '''
    Define expected search string for tropomi files based on
    species and analysis mode.
    
    S5P_{analysis_mode}_L2__{search_species}*.nc
    
    e.g. search_species for "ch4" is "CH4____"

    Note: doesn't include path to filenames
    
    Args:
        species (str, optional):
            Currently only "ch4" for methane
        analysis_mode (str, optional):
            One of "OFFL" (offline mode) or "RPRO" (repreocessed).
    
    Returns:
        str:
            Filename search string including a wildcard (*)
    '''
    if species == "ch4":
        search_species = "CH4____"
    
    search_str = f"S5P_{analysis_mode}_L2__{search_species}*.nc"
    
    return search_str
    

def tropomi_regrid(start_date,end_date,lat_bounds,lon_bounds,coord_bin,
                   time_increment="1min",
                   quality_filt=True,
                   regrid_method="conservative_normed",
                   input_directory=input_directory,
                   allow_parallel=False,
                   verbose=False,
                   write = False, output_directory = None):
    '''
    The tropomi_regrid function combines tropomi data for a given date
    range onto a regular grid. This also involves grouping together
    measurement in time and the minimum time granularity is 
    determined by the the time_increment input.
    
    Args:
        start_date (str) :
            Start of the date range in the form "YYYY-MM-DD"
        end_date (str) :
            Up to but not including this date. "YYYY-MM-DD"
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
            Default = "1min" # 1 minute
        quality_filt (bool, optional) :
            Filter by the recommended quality conditions (combined
            filter > 0.5).
            Default = True
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
        verbose (bool, optional) :
            Print more descriptive output as function executes.
            Default = False
    
    Returns:
        xarray.Dataset :
            Grouped tropomi data on a regular grid with overall 
            dimensions of:
                'time', 'lat', 'lon, ['layer', 'layer_bound']
    '''
    
    species = "ch4"

    if allow_parallel:
        clean_up_weights = False
    else:
        clean_up_weights = True

    time = "time"
    time_full = "delta_time" # Time as np.datetime objects

    if species == "ch4":
        directory = os.path.join(input_directory, "CH4")
        #output_directory = base_output_directory / "CH4"
    analysis_mode = "OFFL"
    search_str = define_tropomi_search_str(species,analysis_mode)
        
    filenames = gosat_fn.extract_files(directory,search_str,start_date,end_date)
    
    if len(filenames) == 0:
        alt_analysis_mode = "RPRO"
        if verbose:
            print(f"Couldn't find {analysis_mode} files, looking for {alt_analysis_mode} (reprocessed)")
        search_str_alt = define_tropomi_search_str(species,alt_analysis_mode)
        #search_str_alt = search_str.replace(current_analysis_mode,analysis_mode)
        filenames = gosat_fn.extract_files(directory,search_str_alt,start_date,end_date)
    
    if len(filenames) == 0:
        print(f"No files found for date: {start_date} - {end_date}")
    
    data_regridded = None
    input_filenames = []
    for i,filename in enumerate(filenames):
        
        print(f"Processing input file: {filename}")
        ## TODO: Currently seems to be relying on edge_coords being
        # set to bounds - need to to update this and follow it through
        # the chain.
        data = preProcessFile(filename, add_corners=True)

        # Make delta_time into a coordinate. This will be used
        # if splitting my time increments and also means that
        # these values do not get turned into NaT values if using
        # filtering by the quality flag.
        dt0 = data[time_full].isel({time:0},drop=True)
        data = data.assign_coords(coords={time_full:dt0})

        if quality_filt:
            data = data.where(data["qa_pass"])
        
        # Regrid data for each file
        regridded = regrid_orbit(data,lat_bounds,lon_bounds,coord_bin,
                                   method=regrid_method,
                                   time_increment=time_increment,
                                   clean_up_weights=clean_up_weights)
        
        if regridded is not None:
            if data_regridded is None:
                data_regridded = regridded
            else:
                data_regridded = xr.concat([data_regridded,regridded],dim=time)
            # Record filenames which contain data
            input_filenames.append(filename)

    short_filenames = [os.path.split(fname)[1] for fname in input_filenames]
    short_filenames = ','.join(short_filenames)
    
    data_regridded.attrs["input_filename"] = short_filenames
    data_regridded.attrs["analysis_mode"] = f"Input files are from {analysis_mode} mode of analysis."

    if quality_filt:
        data_regridded.attrs["qa_filter"] = "Quality filter of > 0.5 applied"
    
    
    return data_regridded
    

def tropomi_process(site,start_date,end_date,lat_bounds,lon_bounds,
                    coord_bin,
                    network=None,
                    species="ch4",
                    time_increment="1min",
                    quality_filt=True,
                    regrid_method="conservative_normed",
                    input_directory=input_directory,
                    write_name=False,name_directory=name_csv_directory,
                    max_name_level=None,max_name_points=50,
                    file_per_day=True,
                    use_name_pressure=False,
                    pressure_base_dir=name_pressure_directory,
                    pressure_domain=None,pressure_max_days=31,
                    pressure_day_template=True,
                    output_directory=obs_directory,
                    overwrite=False,
                    allow_parallel=True,verbose=False):
    '''
    Process tropomi data within date range including re-gridding onto 
    a regular latitute-longitude grid.
    
    Args:
        site (str) :
            Site name for TROPOMI output e.g. "TROPOMI-BRAZIL"
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
        network (str/None, optional) :
            Specify the network for this site. By default this will
            use acrg_site_info file to extract this information
            based on the site name and will use the network to create
            the output file path for the obs and NAME data.
            Default = None    # Extract from definition file
        species (str, optional) :
            Species to extract. One of:
                - "ch4" (methane)
            TODO: This needs to be updated to accept other species
            Default = "ch4"
        time_increment (str, optional) :
            Time window to group tropomi points.
            Default = "10s" # 10 seconds
        quality_filt (bool, optional) :
            Filter by the recommended quality conditions (combined
            filter > 0.5).
            Default = True
        regrid_method (str, optional) :
            Regridding method to use. Options include:
                - "conservative"
                - "conservative_normed"
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
        max_name_level (int/None, optional) :
            Maximum level to extract from TROPOMI for running through 
            NAME.
            Default = None     # Uses all levels
        max_name_points (int/None, optional) :
            Maximum number of points within in each NAME input file.
            Default = 50
        file_per_day (bool, optional) :
            Group NAME points across a day into one file (subject to
            max_name_points) rather than outputting this as one point per
            file.
            Default = True
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
    
    
    ##TODO: Allow more species to be specified and filter to choose 
    # correct search_species string
    ##TODO: Make this return something as well as writing to file?
    '''
    timing_1 = timing_module.time()
    
    if allow_parallel:
        clean_up_weights = False
    else:
        clean_up_weights = True
    
    print(f"Processing tropomi files for date range: {start_date} - {end_date}")

    time = "time"
    time_full = "delta_time" # Time as np.datetime objects
    #latlon=("latitude", "longitude")
    #pos_coords=("scanline", "ground_pixel")
    
    if species == "ch4":
        directory = input_directory + "CH4"
        #output_directory = base_output_directory / "CH4"
        search_species = "CH4____"

    # analysis_mode = "OFFL"
    #### TEMPORARY CHANGE TO RPRO BECAUSE OF UNIDENTIFIED PROBLEM WITH IF STATEMENT BELOW ####
    analysis_mode = "RPRO"
    search_str = f"S5P_{analysis_mode}_L2__{search_species}*.nc"

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
        input_filenames = []
        for i,filename in enumerate(filenames):
            
            ## ONLY WORKS FOR METHANE FILES ("ch4") AT THE MOMENT!
            print(f"Processing input file: {filename}")
            data = preProcessFile(filename)


            # Making delta_time into a coordinate. This will be used
            # if splitting by time increments and also means that
            # these values do not get turned into NaT values if using
            # filtering by the quality flag.
            dt0 = data[time_full].isel({time:0},drop=True)
            data = data.assign_coords(coords={time_full:dt0})

            if quality_filt:
                data = data.where(data["qa_pass"])
            
            # Regrid data for each file
            regridded = regrid_orbit(data,lat_bounds,lon_bounds,coord_bin,
                                       method=regrid_method,
                                       time_increment=time_increment,
                                       clean_up_weights=clean_up_weights)
            
            if regridded is not None:
                if data_regridded is None:
                    data_regridded = regridded
                else:
                    data_regridded = xr.concat([data_regridded,regridded],dim=time)
                # Record filenames which contain data
                input_filenames.append(filename)

            # Check if there is no tropomi data in data_regridded, if so, skip this day

        if len(input_filenames) == 0:
            print(f"No tropomi data within output grid for {date}. Skipping processing of this date.")
        else:
        
            ##TODO: Filter any points where some specific values are not defined
            # e.g. methane, pressure levels       
            if verbose:
                print("Unravelling grid...") 
            
            data_timeseries = unravel_grid(data_regridded)
            
            short_filenames = [os.path.split(fname)[1] for fname in input_filenames]
            short_filenames = ','.join(short_filenames)
            
            if quality_filt:
                data_timeseries.attrs["qa_filter"] = "Quality filter of > 0.5 applied"
            data_timeseries.attrs["input_filename"] = short_filenames
            data_timeseries.attrs["analysis_mode"] = f"Input files are from {analysis_mode} mode of analysis."
            
            if max_name_level is None:
                print(f'Max name level: {data_timeseries["layer"].values[-1]}')
                data_timeseries.attrs["max_level"] = data_timeseries["layer"].values[-1]
            else:
                data_timeseries.attrs["max_level"] = max_name_level
            
            ## Rename "layer" output to match "lev" used within GOSAT files.
            # TROPOMI - level describes the bounds of each layer, layer = level-1
            # GOSAT - level describes the midpoint of each layer
            data_timeseries = data_timeseries.rename({"layer":"lev"})
            #if "layer_bound" in data_timeseries.dims:
            #    data_timeseries = data_timeseries.rename({"layer_bound":"lev_bound"})        
            
            write_tropomi_output(data_timeseries,site,date,species=species,
                                network=network,
                                output_directory=output_directory,
                                overwrite = overwrite)
                    
            if write_name:
                write_tropomi_NAME(data_timeseries,site=site,max_level=max_name_level,
                                    max_points=max_name_points,network=network,
                                    use_name_pressure=use_name_pressure,
                                    pressure_base_dir=pressure_base_dir,
                                    pressure_domain=pressure_domain,
                                    pressure_max_days=pressure_max_days,
                                    pressure_day_template=pressure_day_template,
                                    name_directory=name_directory,
                                    overwrite = overwrite)

        print(f"\nTime to execute: {timing_module.time() - timing_1}\n")

        