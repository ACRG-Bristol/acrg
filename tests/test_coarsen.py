#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:28:12 2019

Tests for checking the basis function generation:
    coarsen

Testing the shape of the coarsened array.
Testing the values in a coarsened array.

@author: mi19881
"""

import os
import netCDF4 as nc
from builtins import range
import numpy as np
from acrg.grid.coarsen import coarsen
from acrg.config.paths import Paths

acrg_path = Paths.acrg


def test_acrg_coarsen_shape():
    '''
    Test if the shape of the regriddded 2D array is correct
    '''
    f = os.path.join(acrg_path,"tests/files/LPDM/emissions/EUROPE/ch4_EUROPE_reduced_2013.nc")
    f_emission = nc.Dataset(f)

    lon_array = f_emission.variables['lon'][:]
    lat_array = f_emission.variables['lat'][:]
    data_array = f_emission.variables['flux'][:]
    
    data_array_sub = data_array[:,:,0]
    
    coarsened_data, coarsend_lat, coarsened_lon = coarsen(arrayFine = data_array_sub, latFine = lat_array, lonFine = lon_array, factor = 2, mean = True)
    coarsened_data_shape = coarsened_data.shape
    
    assert coarsened_data_shape == (146, 195)

def test_acrg_coarsen():
    '''
    Test if the coarsen function is creating a new array with the correct values
    '''
    data_array = np.array(([1.,2.,3.,4.,5.,6.],[4.,5.,6.,7.,8.,9.],[7.,8.,9.,10.,11.,12.],[10.,11.,12.,13.,14.,15.],[13.,14.,15.,16.,17.,18.],[16.,17.,18.,19.,20.,21.]))
    lat_array = np.array([28.,29.,30.,31.,32.,33.])
    lon_array = np.array([25.,26.,27.,28.,29.,30.])
    
    coarsened_data, coarsend_lat, coarsened_lon = coarsen(arrayFine = data_array, latFine = lat_array, lonFine = lon_array, factor = 2, mean = True)
    coarsened_data_shape = coarsened_data.shape
    
    assert coarsened_data_shape == (3, 3)
    assert np.all(coarsened_data == np.array(([3.,5.,7.],[9.,11.,13.],[15.,17.,19.])))
