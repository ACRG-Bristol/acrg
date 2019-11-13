#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 08:30:03 2019

@author: al18242

Functions for processing and using TROPOMI data
"""
import numpy as np
import xarray as xr

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

def preProcessFile(filename):
    '''Combine important information from different netCDF groups into a single data frame
    
    input
        filename - string of full filename for TROPOMI file
    
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
    
    #calculate the meshgrid of corner points in the scan grid for further processing
    lat, lon = cornerGrid(tropomi_data_geo)
    tropomi_data['lat_corners'] = (["time", "scanline_c", "ground_pixel_c"], np.expand_dims(lat,0))
    tropomi_data['lon_corners'] = (["time", "scanline_c", "ground_pixel_c"],np.expand_dims(lon,0))
    
    #flag quality to remove outputs not flagged as 'success' and use recommended acceptable qa_values
    # Note: it may be the case that any value that passes the second condition necessarily passes the first
    quality_flag = (tropomi_data_aux.processing_quality_flags.values.astype(int) & 0xFF == 0) & \
                    (tropomi_data.qa_value.values > 0.5)
    
    tropomi_data["qa_pass"] = (["time", "scanline", "ground_pixel"], quality_flag)
    
    return tropomi_data