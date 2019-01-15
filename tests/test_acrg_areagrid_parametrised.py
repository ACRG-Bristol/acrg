#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue 15 Jan 2019 

Simple tests for acrg_grid.py

This script should return a 2D array of areas based according to an array of 
inputs of lats by an array of input of lons.

Using a Parametrised test

@author mi19881
"""

import numpy as np
import netCDF4 as nc
import pytest
from acrg_grid.areagrid import areagrid

def test_acrg_areagrid_shape_size():
    '''
    Test if the shape and size of the output 2D array is correct
    '''
    lat_array = np.array([-31.0, -32.0, 33.0])
    lon_array = np.array([18.0, 19.0])
    temp_areagrid = areagrid(lat_array, lon_array)
    temp_areagrid_shape = temp_areagrid.shape
    temp_areagrid_size = temp_areagrid.size
    assert len(lat_array) == 3
    assert len(lon_array) == 2
    assert temp_areagrid_shape == (len(lat_array), len(lon_array))
    assert temp_areagrid_size == len(lat_array)*len(lon_array)


def test_acrg_areagrid_netcdf_shape_size():
    '''
    Test if the shape and size of the output 2D array is correct for 
    a netcdf file
    '''
    f = nc.Dataset( '/home/mi19881/programs/acrg/tests/files/NAME/emissions/EUROPE/ch4_EUROPE_2014.nc')
    lon_array = f.variables['lon'][:]
    lat_array = f.variables['lat'][:]
    temp_areagrid = areagrid(lat_array, lon_array)
    temp_areagrid_shape = temp_areagrid.shape
    temp_areagrid_size = temp_areagrid.size
    assert temp_areagrid_shape == (len(lat_array), len(lon_array))
    assert temp_areagrid_size == len(lat_array)*len(lon_array)
 

@pytest.mark.parametrize("lat_array, lon_array, expected", [
        (np.array([50.5, 51.5, 52.5]), np.array([-1.0, 0.0, 1.0]), [[7855.9, 7855.9, 7855.9], [7688.4, 7688.4, 7688.4], [7518.6, 7518.6, 7518.6]]),
        (np.array([46.5, 47.0, 47.5]), np.array([6.5, 7.0, 7.5]), [[2125.4, 2125.4, 2125.4], [2105.8, 2105.8, 2105.8], [2086. , 2086. , 2086. ]]),
        (np.array([-22.0, -20.0, -18.0]), np.array([-53.0, -51.0, -49.0]), [[45803.3, 45803.3, 45803.3], [46421.2, 46421.2, 46421.2], [46982.6, 46982.6, 46982.6]]),
        (np.array([-33.5, -32.5, -31.5]), np.array([18.5, 19.5, 20.5]), [[10299. , 10299. , 10299. ], [10416.4, 10416.4, 10416.4], [10530.6, 10530.6, 10530.6]]),
        (np.array([-37.5, -37.0, -36.5]), np.array([144.0, 144.5, 145.0]), [[2449.6, 2449.6, 2449.6], [2465.9, 2465.9, 2465.9], [2482. , 2482. , 2482. ]]),
        (np.array([64.0, 65.0, 66.0]), np.array([-109.0, -108.0, -107.0]), [[5414.1, 5414.1, 5414.1], [5219.6, 5219.6, 5219.6], [5023.4, 5023.4, 5023.4]]),
        ])
def test_acrg_areagrid_areas(lat_array, lon_array, expected):
    '''
    Test if the areas are calculated as for the benchmark code
    '''
    assert np.all(np.round(areagrid(lat_array, lon_array)/1000000.,1) == expected)

