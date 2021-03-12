#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:08:56 2020

@author: rt17603
"""

import pytest
import os
import numpy as np
import xarray as xr
from pathlib import Path
from acrg_config.paths import paths
from acrg_satellite import tropomi

acrg_path = paths.acrg
test_data_dir = acrg_path / "tests/files/data/obs_raw/TROPOMI"

@pytest.fixture(scope="module")
def tropomi_filename():
    # Covers section of South America
    filename = "S5P_OFFL_L2__CH4____20190726T155905_20190726T174035_09240_01_010302_20190801T172237.nc"
    filename = test_data_dir / filename
    return filename

@pytest.fixture(scope="module")
def tropomi_raw_dataset(tropomi_filename):
    ds = tropomi.preProcessFile(tropomi_filename)
    return ds

def extent_Brazil():

    lat_min = -35.8 # -35.7525
    lat_max = +7.3  # +7.25138
    lon_min = -76.0 # -75.98444
    lon_max = -32.7 # -32.79305
    
    dlat = 0.5
    dlon = 0.5
    
    #output_lat = np.arange(lat_min,lat_max+dlat,dlat)
    #output_lon = np.arange(lon_min,lon_max+dlon,dlon)

    lat_bounds = [lat_min, lat_max]
    lon_bounds = [lon_min, lon_max]
    coord_bin = [dlat, dlon]

    return lat_bounds, lon_bounds, coord_bin

# Can use remove_weight_files function written in tropomi module
#def remove_weight_files(directory=".",method="conservative"):
#
#    weight_files = Path(directory).glob(f"{method}_*.nc")
#    for file in weight_files:
#        os.remove(file)

#%% Mocking dummy dataset

def dummy_coords():
    '''
    Creating dummy coordinates for a mocking.
    
    Dimensions:
        - time = 1
        - scanline = 4
        - ground_pixel = 2
        - level = 3
        - layer = 4 (level + 1)
    
    Creating offset latitude and longitude 2D values to allow us 
    to test regridding onto a regular grid.
    
    - latitude ("time", "scanline", "ground_pixel")
     - ground_pixel 1 - 9.8 - 10.4 in 0.2 degree increments
     - ground_pixel 2 - 10.0 - 10.6 in 0.2 degree increments
    
    - longitude ("time", "scanline", "ground_pixel")
     - ground_pixel 1 - 19.8 - 20.4 in 0.2 degree increments
     - ground_pixel 2 - 20.0 - 20.6 in 0.2 degree increments
    
    Returns:
        dict :
            In correct format for creation of an xarray Dataset.
    '''
    time = np.array(["2019-01-01"], dtype="datetime64[ns]")

    n_scanline = 4
    n_ground_pixel = 2
    n_level = 3

    scanline = np.arange(0.0, n_scanline, 1.0)
    ground_pixel = np.arange(0.0, n_ground_pixel, 1.0)
    corner = np.arange(0.0, 4.0, 1.0)
    layer = np.arange(0.0, n_level-1, 1.0)
    level = np.arange(0.0, n_level, 1.0)
    
    latitude = np.array([[[9.8, 10.0],
                          [10.0, 10.2],
                          [10.2, 10.4],
                          [10.4, 10.6]]])
    
    longitude = np.array([[[19.8, 20.0],
                           [20.0, 20.2],
                           [20.2, 20.4],
                           [20.4, 20.6]]])
    
    coords = {}
    coords["scanline"] = scanline
    coords["ground_pixel"] = ground_pixel
    coords["time"] = time
    coords["corner"] = corner
    coords["layer"] = layer
    coords["level"] = level
    
    coords["latitude"] = (("time", "scanline", "ground_pixel"), latitude)
    coords["longitude"] = (("time", "scanline", "ground_pixel"), longitude)
    
    return coords    


def dummy_main_product():
    '''
    Creating dummy Dataset for group="PRODUCT" within a tropomi file.
    
    Data variables:
     - delta_time ("time", "scanline")
        - Four times unevenly spaced between "2019-01-01T16:00:00"
     and "2019-01-01T17:00:00"
     - time_utc ("time", "scanline")
        - Same as delta_time but as a string
     - qa_value ("time", "scanline", "ground_pixel")
        - Values between 0.0 and 1.0
     - methane_mixing_ratio ("time", "scanline", "ground_pixel")
        - Values of 1, 0 or np.nan
     - methane_mixing_ratio_precision
        - methane_mixing_ratio * 0.1
     - methane_mixing_ratio_bias_corrected 
        - methane_mixing_ratio + 1.0
    
    Coords:
        Using values from dummy_coords() (though some assumptions
        made about the dimensionality)
    '''
    ### EXAMPLE
    # <xarray.Dataset>
    # Dimensions:                              (corner: 4, ground_pixel: 215, layer: 12, level: 13, scanline: 3246, time: 1)
    # Coordinates:
    #   * scanline                             (scanline) float64 0.0 ... 3.245e+03
    #   * ground_pixel                         (ground_pixel) float64 0.0 ... 214.0
    #   * time                                 (time) datetime64[ns] 2019-07-26
    #   * corner                               (corner) float64 0.0 1.0 2.0 3.0
    #   * layer                                (layer) float64 0.0 1.0 ... 10.0 11.0
    #   * level                                (level) float64 0.0 1.0 ... 11.0 12.0
    #     latitude                             (time, scanline, ground_pixel) float32 ...
    #     longitude                            (time, scanline, ground_pixel) float32 ...
    # Data variables:
    #     delta_time                           (time, scanline) datetime64[ns] ...
    #     time_utc                             (time, scanline) object ...
    #     qa_value                             (time, scanline, ground_pixel) float32 ...
    #     methane_mixing_ratio                 (time, scanline, ground_pixel) float32 ...
    #     methane_mixing_ratio_precision       (time, scanline, ground_pixel) float32 ...
    #     methane_mixing_ratio_bias_corrected  (time, scanline, ground_pixel) float32 ...
    ###
    
    data = {}
    
    times = [["2019-01-01T16:00:00.0000Z", 
             "2019-01-01T16:00:30.0000Z", 
             "2019-01-01T16:29:00.0000Z",
             "2019-01-01T17:00:00.0000Z"]]
    
    delta_time = np.array(times, dtype="datetime64[ns]")
    time_utc = np.array(times, dtype="object")
    qa_value = np.array([[[1.0, 1.0],
                          [1.0, 1.0],
                          [0.6, 0.6],
                          [0.0, 0.0]]])
    
    methane_mixing_ratio = np.array([[[1, 1],
                                      [1, 0],
                                      [1, 1],
                                      [1, np.nan]]])
                                
    methane_mixing_ratio_precision = methane_mixing_ratio*0.1
    methane_mixing_ratio_bias_corrected = methane_mixing_ratio + 1.0
    
    
    data = {}
    data["delta_time"] = (("time", "scanline"), delta_time)
    data["time_utc"] = (("time", "scanline"), time_utc)
    data["qa_value"] = (("time", "scanline", "ground_pixel"), qa_value)
    data["methane_mixing_ratio"] = (("time", "scanline", "ground_pixel"), methane_mixing_ratio)
    data["methane_mixing_ratio_precision"] = (("time", "scanline", "ground_pixel"), methane_mixing_ratio_precision)
    data["methane_mixing_ratio_bias_corrected"] = (("time", "scanline", "ground_pixel"), methane_mixing_ratio_bias_corrected)
    
    coords = dummy_coords()

    ## MAY NEED TO ADD ATTRS
    
    ds = xr.Dataset(data, coords)
    
    return ds

def dummy_geo_product():
    '''
    Creating dummy Dataset for group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS"
    within a tropomi file.

    latitude and longitude bounds derived to not lie on a grid but
    lie in a diamond pattern around the centre point.

    Data variables:
        latitude_bounds ("time", "scanline", "ground_pixel", "corner")
         - Based on dummy_coords() latitude values (with some assumptions
           about pixel size of 0.2)
        longitude_bounds ("time", "scanline", "ground_pixel", "corner")
         - Based on dummy_coords() latitude values (with some assumptions
           about pixel size of 0.2)
    
    Coords:
        None - not defined in tropomi input so not defined here.
    '''
    ### EXAMPLE
    #     <xarray.Dataset>
    # Dimensions:                (corner: 4, ground_pixel: 215, scanline: 3246, time: 1)
    # Dimensions without coordinates: corner, ground_pixel, scanline, time
    # Data variables:
    #     satellite_latitude     (time, scanline) float32 ...
    #     satellite_longitude    (time, scanline) float32 ...
    #     satellite_altitude     (time, scanline) float32 ...
    #     satellite_orbit_phase  (time, scanline) float32 ...
    #     solar_zenith_angle     (time, scanline, ground_pixel) float32 ...
    #     solar_azimuth_angle    (time, scanline, ground_pixel) float32 ...
    #     viewing_zenith_angle   (time, scanline, ground_pixel) float32 ...
    #     viewing_azimuth_angle  (time, scanline, ground_pixel) float32 ...
    #     latitude_bounds        (time, scanline, ground_pixel, corner) float32 ...
    #     longitude_bounds       (time, scanline, ground_pixel, corner) float32 ...
    #     geolocation_flags      (time, scanline, ground_pixel) float32 ...
    ###

    coords = dummy_coords()

    ## Create latitude bounds (edges) based on latitude values
    latitude_dims, latitude = coords["latitude"]
    lat_diff = 0.2 # Make more generic??
    #lat_size = (latitude[1:] - latitude[:-1]).mean()
    
    lat_corner_1 = latitude - lat_diff#/2
    lat_corner_2 = latitude
    lat_corner_3 = latitude + lat_diff#/2
    lat_corner_4 = latitude
    
    latitude_bounds = np.stack((lat_corner_1, lat_corner_2, 
                                lat_corner_3, lat_corner_4), axis=-1)

    ## Create longitude bounds (edges) based on longitude values
    longitude_dims, longitude = coords["latitude"]
    lon_diff = 0.2 # Make more generic??
    
    lon_corner_1 = longitude
    lon_corner_2 = longitude + lon_diff#/2
    lon_corner_3 = longitude
    lon_corner_4 = longitude - lon_diff#/2
    
    longitude_bounds = np.stack((lon_corner_1, lon_corner_2, 
                                lon_corner_3, lon_corner_4), axis=-1)    
    

    dims_corner = latitude_dims + ("corner",)

    ## Set up dictionary to contain correct format to convert to an
    # xarray Dataset
    # Note: Dimensions without coords for this group matches to tropomi 
    # data files.    
    data = {}
    data["latitude_bounds"] = (dims_corner, latitude_bounds)
    data["longitude_bounds"] = (dims_corner, longitude_bounds)

    ## MAY NEED TO ADD ATTRS
    
    ds = xr.Dataset(data)
    
    return ds

def dummy_detailed_results():
    '''
    Creating dummy Dataset for group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS"
    within a tropomi file.

    Data variables:
        column_averaging_kernel ("time", "scanline", "ground_pixel", "layer")
         - All values set to 0
        processing_quality_flags ("time", "scanline", "ground_pixel")
        - All values set to 1
    
    Coords:
        None - not defined in tropomi input so not defined here.    
    '''
    ### EXAMPLE
    #     <xarray.Dataset>
    # Dimensions:                                     (ground_pixel: 215, layer: 12, scanline: 3246, time: 1)
    # Dimensions without coordinates: ground_pixel, layer, scanline, time
    # Data variables:
    #     processing_quality_flags                    (time, scanline, ground_pixel) float64 ...
    #     number_of_spectral_points_in_retrieval      (time, scanline, ground_pixel) float32 ...
    #     number_of_spectral_points_in_retrieval_NIR  (time, scanline, ground_pixel) float32 ...
    #     column_averaging_kernel                     (time, scanline, ground_pixel, layer) float32 ...
    #     carbonmonoxide_total_column                 (time, scanline, ground_pixel) float32 ...
    #     carbonmonoxide_total_column_precision       (time, scanline, ground_pixel) float32 ...
    #     water_total_column                          (time, scanline, ground_pixel) float32 ...
    #     water_total_column_precision                (time, scanline, ground_pixel) float32 ...
    #     aerosol_size                                (time, scanline, ground_pixel) float32 ...
    #     aerosol_size_precision                      (time, scanline, ground_pixel) float32 ...
    #     aerosol_number_column                       (time, scanline, ground_pixel) float32 ...
    #     aerosol_number_column_precision             (time, scanline, ground_pixel) float32 ...
    #     aerosol_mid_altitude                        (time, scanline, ground_pixel) float32 ...
    #     aerosol_mid_altitude_precision              (time, scanline, ground_pixel) float32 ...
    #     surface_albedo_SWIR                         (time, scanline, ground_pixel) float32 ...
    #     surface_albedo_SWIR_precision               (time, scanline, ground_pixel) float32 ...
    #     surface_albedo_NIR                          (time, scanline, ground_pixel) float32 ...
    #     surface_albedo_NIR_precision                (time, scanline, ground_pixel) float32 ...
    #     aerosol_optical_thickness_SWIR              (time, scanline, ground_pixel) float32 ...
    #     aerosol_optical_thickness_NIR               (time, scanline, ground_pixel) float32 ...
    #     wavelength_calibration_offset_SWIR          (time, scanline, ground_pixel) float32 ...
    #     wavelength_calibration_offset_NIR           (time, scanline, ground_pixel) float32 ...
    #     chi_square                                  (time, scanline, ground_pixel) float32 ...
    #     chi_square_SWIR                             (time, scanline, ground_pixel) float32 ...
    #     chi_square_NIR                              (time, scanline, ground_pixel) float32 ...
    #     degrees_of_freedom                          (time, scanline, ground_pixel) float32 ...
    #     degrees_of_freedom_methane                  (time, scanline, ground_pixel) float32 ...
    #     degrees_of_freedom_aerosol                  (time, scanline, ground_pixel) float32 ...
    #     number_of_iterations                        (time, scanline, ground_pixel) float64 ...
    #     fluorescence                                (time, scanline, ground_pixel) float32 ...
    ###

    ## Grab coordinate definitions
    coords = dummy_coords()
    
    n_time = len(coords["time"])
    n_scanline = len(coords["scanline"])
    n_ground_pixel = len(coords["ground_pixel"])
    n_layer = len(coords["layer"])
    
    dims = ("time", "scanline", "ground_pixel")
    ndims = (n_time, n_scanline, n_ground_pixel)
    
    dims_layer = dims + ("layer",)
    ndims_layer = ndims + (n_layer,)    
    
    # Use one value for all entries of :
    # - column_averaging_kernel
    # - processing_quality_flags    
    
    column_averaging_kernel = np.ones(ndims_layer)
    processing_quality_flags = np.zeros(ndims)

    ## Set up dictionary to contain correct format to convert to an
    # xarray Dataset
    # Note: Dimensions without coords for this group matches to tropomi 
    # data files.
    data = {}
    data["column_averaging_kernel"] = (dims_layer, column_averaging_kernel)
    data["processing_quality_flags"] = (dims, processing_quality_flags)

    ds = xr.Dataset(data)

    ## MAY NEED TO ADD ATTRS
    
    return ds


def dummy_support_input():
    '''
    Creating dummy Dataset for group="PRODUCT/SUPPORT_DATA/INPUT_DATA"
    within a tropomi file.

    Data variables:
        surface_pressure ("time", "scanline", "ground_pixel")
            All values are 100000 (in Pa)
        pressure_interval ("time", "scanline", "ground_pixel")
            All values are 8000 (in Pa)
        methane_profile_apriori ("time", "scanline", "ground_pixel", "layer")
            All values of 1.0
        dry_air_subcolumns ("time", "scanline", "ground_pixel", "layer")
            30000 with intervals of 50 for each layer (in mol/m2??)
    
    Coords:
        None - not defined in tropomi input so not defined here.       
    '''
    ### EXAMPLE
    # <xarray.Dataset>
    # Dimensions:                                       (ground_pixel: 215, layer: 12, level: 13, scanline: 3246, time: 1)
    # Dimensions without coordinates: ground_pixel, layer, level, scanline, time
    # Data variables:
    #     surface_altitude                              (time, scanline, ground_pixel) float32 ...
    #     surface_altitude_precision                    (time, scanline, ground_pixel) float32 ...
    #     surface_classification                        (time, scanline, ground_pixel) float32 ...
    #     instrument_configuration_identifier           (time, scanline) float64 ...
    #     instrument_configuration_version              (time, scanline) float32 ...
    #     scaled_small_pixel_variance                   (time, scanline, ground_pixel) float32 ...
    #     eastward_wind                                 (time, scanline, ground_pixel) float32 ...
    #     northward_wind                                (time, scanline, ground_pixel) float32 ...
    #     methane_profile_apriori                       (time, scanline, ground_pixel, layer) float32 ...
    #     altitude_levels                               (time, scanline, ground_pixel, level) float32 ...
    #     dry_air_subcolumns                            (time, scanline, ground_pixel, layer) float32 ...
    #     surface_pressure                              (time, scanline, ground_pixel) float32 ...
    #     pressure_interval                             (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_SWIR_IFOV                (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_SWIR_OFOVa               (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_SWIR_OFOVb               (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_SWIR_OFOVc               (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_NIR_IFOV                 (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_NIR_OFOVa                (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_NIR_OFOVb                (time, scanline, ground_pixel) float32 ...
    #     cloud_fraction_VIIRS_NIR_OFOVc                (time, scanline, ground_pixel) float32 ...
    #     reflectance_cirrus_VIIRS_SWIR                 (time, scanline, ground_pixel) float32 ...
    #     reflectance_cirrus_VIIRS_NIR                  (time, scanline, ground_pixel) float32 ...
    #     apparent_scene_pressure                       (time, scanline, ground_pixel) float32 ...
    #     apparent_scene_pressure_standard_deviation    (time, scanline, ground_pixel) float32 ...
    #     methane_weak_twoband_total_column             (time, scanline, ground_pixel) float32 ...
    #     methane_strong_twoband_total_column           (time, scanline, ground_pixel) float32 ...
    #     methane_ratio_weak_strong_standard_deviation  (time, scanline, ground_pixel) float32 ...
    #     water_weak_twoband_total_column               (time, scanline, ground_pixel) float32 ...
    #     water_strong_twoband_total_column             (time, scanline, ground_pixel) float32 ...
    #     water_ratio_weak_strong_standard_deviation    (time, scanline, ground_pixel) float32 ...
    #     fluorescence_apriori                          (time, scanline, ground_pixel) float32 ...
    ###
    
    ## Grab coordinate definitions
    coords = dummy_coords()
    
    n_time = len(coords["time"])
    n_scanline = len(coords["scanline"])
    n_ground_pixel = len(coords["ground_pixel"])
    n_layer = len(coords["layer"])

    dims = ("time", "scanline", "ground_pixel")
    ndims = (n_time, n_scanline, n_ground_pixel)
    
    dims_layer = dims + ("layer",)
    ndims_layer = ndims + (n_layer,)

    ## Use one value for all entries of :
    # - surface_pressure
    # - pressure_interval
    # - methane_profile_apriori
    
    sp_value = 100000 # surface_pressure - Pa
    pi_value = 8000   # pressure_interval - Pa

    surface_pressure = np.zeros(ndims)
    surface_pressure += sp_value
    
    pressure_interval = np.zeros(ndims)
    pressure_interval += pi_value
    
    methane_profile_apriori = np.ones(ndims_layer)
    
    ## Create different dry_air_subcolumn value per layer
    # - From dry_air_bottom and decreasing in increments of dry_air_diff
    #dry_air_bottom = 28500 # mol/m2
    dry_air_bottom = 30000 # mol/m2
    dry_air_diff = 50 # mol/m2
    dry_air_top = dry_air_bottom - dry_air_diff*n_layer
    
    # Create range and then invert
    dry_air_subcolumns_l = np.arange(dry_air_top, dry_air_bottom, dry_air_diff)
    dry_air_subcolumns_l = np.flip(dry_air_subcolumns_l)
    
    # Tile function creates number of repeats across each specified axis.
    # - for last dimension (layer) don't want repeats so specify (1,)
    dry_air_subcolumns = np.tile(dry_air_subcolumns_l, ndims+(1,))
    
    ## Set up dictionary to contain correct format to convert to an
    # xarray Dataset
    # Note: Dimensions without coords for this group matches to tropomi 
    # data files.
    data = {}
    
    data["surface_pressure"] = (dims, surface_pressure)
    data["pressure_interval"] = (dims, pressure_interval)
    data["methane_profile_apriori"] = (dims_layer, methane_profile_apriori)
    data["dry_air_subcolumns"] = (dims_layer, dry_air_subcolumns)
    
    ds = xr.Dataset(data)
    
    ## MAY NEED TO ADD ATTRS
    
    return ds

def dummy_main_product_dt():
    '''
    Creating dummy Dataset for group="PRODUCT" within a tropomi file
    based on dummy_main_product() but replacing the delta_time value
    with a timedelta64 object rather than datetime64.
    
    This is to allow us to check that the preprocessFile() function
    can handle both cases and will produce the same output regardless.
    '''
    
    ds = dummy_main_product()
    
    delta_time = ds["delta_time"]
    time = ds["time"]
    
    delta_time_dt = delta_time - time
    ds = ds.assign({"delta_time":delta_time_dt})
    
    return ds
    

def mock_open_dataset(filename, group):
    '''
    Function to use for mocking xr.open_dataset call within
    tropomi.preProcessFile() function.
    
    Uses "group" to determine which dummy Dataset to return.
    
    Uses some "filename" inputs to determine special cases, otherwise
    uses the default:
        "test_delta_time.nc"
         - will use dummy_main_product_dt() function as a special
         case to check the delta_time output. All other inputs 
         are the same.
    '''
    
    ## Calls within code to mock:
    #tropomi_data = xr.open_dataset(filename, group="PRODUCT")
    #tropomi_data_aux = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS")
    #tropomi_data_geo = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS")
    #tropomi_data_input = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA")

    if group == "PRODUCT":
        if filename == "test_delta_time.nc":
            ds = dummy_main_product_dt()
        else:
            ds = dummy_main_product()
    elif group == "PRODUCT/SUPPORT_DATA/GEOLOCATIONS":
        ds = dummy_geo_product()
    elif group == "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS":
        ds = dummy_detailed_results()
    elif group == "PRODUCT/SUPPORT_DATA/INPUT_DATA":
        ds = dummy_support_input()
   
    return ds
    
#%% Test reading and pre-processing tropomi data

@pytest.mark.preprocess
def test_pre_process_xch4(mocker):
    '''
    Test xch4 output is as expected (matches defined input exactly)
    '''     
    mocker.patch.object(
            #'tropomi.preProcessFile.xr.open_dataset',
            tropomi.xr,
            "open_dataset",
            mock_open_dataset
    )

    ds = tropomi.preProcessFile("fake_name.nc")

    xch4 = ds["methane_mixing_ratio"]
    
    ds_expected_main = dummy_main_product()
    xch4_expected = ds_expected_main["methane_mixing_ratio"]
    
    xr.testing.assert_equal(xch4, xch4_expected)

@pytest.mark.preprocess
def test_pre_process_delta_time_1(mocker):
    '''
    Test delta time (datetime64) is included as a datetime object
    in the output.
    '''     
    mocker.patch.object(
            tropomi.xr,
            "open_dataset",
            mock_open_dataset
    )

    ds = tropomi.preProcessFile("fake_name.nc")

    delta_time = ds["delta_time"]
    
    ds_expected_main = dummy_main_product()
    delta_time_expected = ds_expected_main["delta_time"]
    
    xr.testing.assert_equal(delta_time, delta_time_expected)

@pytest.mark.preprocess
def test_pre_process_delta_time_2(mocker):
    '''
    Test delta time (timedelta64) is updated to be the correct datetime 
    object in the output.
    
    This uses a special input filename to preProcessFile of
    "test_delta_time.nc" to extract a special case where the delta_time
    in the input is defined as a timedelta64 object rather than
    datetime64.
    '''     
    mocker.patch.object(
            tropomi.xr,
            "open_dataset",
            mock_open_dataset
    )

    # Use specific filename to make sure correct input is used.
    ds = tropomi.preProcessFile("test_delta_time.nc")

    delta_time = ds["delta_time"]
    
    ds_expected_main = dummy_main_product()
    delta_time_expected = ds_expected_main["delta_time"]
    
    xr.testing.assert_equal(delta_time, delta_time_expected)

def calc_pressure_levels(surface_pressure, interval, nlevel=13):
    '''
    Calculate the pressure levels based on the surface pressure and
    the pressure interval.
    '''
    ndim = surface_pressure.ndim
    shape = (1,)*ndim + (-1,)    
    
    intervals = np.expand_dims(interval,-1)*np.reshape(np.arange(0,nlevel),newshape=shape)
    
    pressure_levels = np.expand_dims(surface_pressure,-1)
    pressure_levels = pressure_levels - intervals
    
    return pressure_levels    
    
@pytest.mark.preprocess
def test_pre_process_pressure(mocker):
    '''
    Test "pressure_levels" values (not coords) have been correctly 
    calculated based on inputs.
    '''     
    mocker.patch.object(
            tropomi.xr,
            "open_dataset",
            mock_open_dataset
    )

    ds = tropomi.preProcessFile("fake_name.nc")

    pressure_levels = ds["pressure_levels"]
    
    # Use known inputs to create expected pressure levels output
    ds_expected_support = dummy_support_input()
    surface_pressure = ds_expected_support["surface_pressure"]
    pressure_interval = ds_expected_support["pressure_interval"]
    nlevel = len(dummy_coords()["level"])
    
    pressure_levels_expect = calc_pressure_levels(surface_pressure, 
                                                  pressure_interval,
                                                  nlevel=nlevel)
    
    # Test numpy array values and shape (coords not important)
    np.testing.assert_allclose(pressure_levels, pressure_levels_expect)

##TODO: ADD TEST FOR PRESSURE_WEIGHTS WHEN WE'RE HAPPY THIS IS RIGHT

@pytest.mark.preprocess
def test_pre_process_quality_flag(mocker):
    '''
    Test "qa_pass" values (not coords) have been correctly included
    based on inputs.
    '''     
    mocker.patch.object(
            tropomi.xr,
            "open_dataset",
            mock_open_dataset
    )

    ds = tropomi.preProcessFile("fake_name.nc")

    qa_pass = ds["qa_pass"]
    
    
    # Use known inputs to create expected pressure levels output
    ds_expected_input = dummy_main_product()
    qa_value_expect = ds_expected_input["qa_value"]
    
    ## Overall qa_pass checks:
    # quality_flag = (tropomi_data_aux.processing_quality_flags.values.astype(int) & 0xFF == 0) & \
    #                 (tropomi_data.qa_value.values > 0.5)
    # but "processing_quality_flags" input have all been set to pass 
    # in our dummy data so should match to this simpler check
    qa_pass_expect = qa_value_expect > 0.5
    
    # Test numpy array values and shape (coords not important)
    np.testing.assert_allclose(qa_pass, qa_pass_expect)
    

#%% Test regridding

##TODO: Add tests using mock dataset - maybe replace Brazil inputs
# to do so

# @pytest.fixture(scope="module")
# def param_1_Brazil():
    
#     time_increment = "1d"
#     lat_bounds, lon_bounds, coord_bin = extent_Brazil()
    
#     param = {}
#     param["time_increment"] = time_increment
#     #param["output_lat"] = output_lat
#     #param["output_lon"] = output_lon
#     param["lat_bounds"] = lat_bounds
#     param["lon_bounds"] = lon_bounds
#     param["coord_bin"] = coord_bin
    
#     return param

# def test_regrid_time(tropomi_raw_dataset,param_1_Brazil):
    
#     out = tropomi.regrid_orbit(tropomi_raw_dataset,**param_1_Brazil)
    
#     #remove_weight_files()
#     #method = param_1_Brazil["method"]
#     #tropomi.remove_weight_files(method=method,path=".")
    
#     assert out

#%% Test reformatting to time, lat, lon

# ##TODO: Add tests using mock dataset - maybe replace Brazil inputs
# # to do so

# @pytest.fixture(scope="module")
# def tropomi_regridded_dataset(tropomi_raw_dataset,param_1_Brazil):
    
#     ds = tropomi.regrid_orbit(tropomi_raw_dataset,**param_1_Brazil)
    
#     #remove_weight_files()
#     #method = param_1_Brazil["method"]
#     #tropomi.remove_weight_files(method=method,path=".")
    
#     return ds


# def test_unravel_grid(tropomi_regridded_dataset):
    
#     out = tropomi.unravel_grid(tropomi_regridded_dataset)
    
#     assert out

#%%

# ##TODO: Add tests using mock dataset - maybe replace Brazil inputs
# # to do so

# @pytest.fixture(scope="module")
# def tropomi_timeseries_dataset(tropomi_regridded_dataset):
    
#     ds = tropomi.unravel_grid(tropomi_regridded_dataset)
    
#     return ds





    
    