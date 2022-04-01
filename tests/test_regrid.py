#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue 15 Jan 2019 

Simple tests for regrid.py

This script should return a 2D or 3D array regridded from the original 

Using a Parametrised test

@author mi19881
"""

import numpy as np
import os
import netCDF4 as nc
import pytest

from acrg.grid.regrid import regrid2d
from acrg.grid.regrid import regrid3d
from acrg.grid.areagrid import areagrid
from acrg.config.paths import Paths

acrg_path = Paths.acrg

@pytest.fixture(scope="module")
def gen_netcdf_file_new():
    ''' Define directory and file containing footprint netcdf example. '''
    filename_path = os.path.join(acrg_path,"tests/files/LPDM/emissions/EUROPE/ch4_EUROPE_reduced_2013.nc")
    return filename_path

@pytest.fixture(scope="module")
def gen_netcdf_file_old():
    ''' Define directory and file containing emissions netcdf example. '''
    filename_path = os.path.join(acrg_path,"tests/files/data/ECMWF_CAMS/BC_CAMS_ch4_80.219.10.83_3x3_2017-12-01.nc")
    return filename_path


def test_acrg_regrid2d_netcdf_shape(gen_netcdf_file_old, gen_netcdf_file_new):
    '''
    Test if the shape of the regriddded 2D array is correct
    '''
    f_new = nc.Dataset(gen_netcdf_file_new)
    f_old = nc.Dataset(gen_netcdf_file_old)
    lon_array = f_old.variables['longitude'][:]
    lat_array = f_old.variables['latitude'][:]
    data_array = f_old.variables['ch4_c'][:]
    lon_out = f_new.variables['lon'][:]
    lat_out = f_new.variables['lat'][:]
    
    data_array_sub = data_array[0,0,:,:]
    
    regridded_data = regrid2d(data_array_sub, lat_array, lon_array, lat_out, lon_out)
    regridded_data_shape = regridded_data[0].shape
    
    assert regridded_data_shape == (len(lat_out), len(lon_out))


def test_acrg_regrid3d_netcdf_shape(gen_netcdf_file_old, gen_netcdf_file_new):
    '''
    Test if the shape of the regriddded 3D array is correct
    '''
    f_new = nc.Dataset(gen_netcdf_file_new)
    f_old = nc.Dataset(gen_netcdf_file_old)
    lon_array = f_old.variables['longitude'][:]
    lat_array = f_old.variables['latitude'][:]
    time_array = f_old.variables['time'][:]
    data_array = f_old.variables['ch4_c'][:]
    lon_out = f_new.variables['lon'][:]
    lat_out = f_new.variables['lat'][:]
    
    data_array_sub = data_array[:,0,:,:]
    data_array_sub = np.transpose(data_array_sub, (1, 2, 0))
    
    regridded_data = regrid3d(data_array_sub, lat_array, lon_array, lat_out, lon_out, time_array)
    regridded_data_shape = regridded_data.shape
    
    assert regridded_data_shape == (len(lat_out), len(lon_out), len(time_array))
    

def test_acrg_regrid2d_netcdf_conservation():
    '''
    Test if the regridding function conserves the emission over the domain when
    using regrid2d
    '''
    data_array = np.array(([1.,2.,3.],[4.,5.,6.],[7.,8.,9.]))
    lat_array = np.array([5.,15.,25.])
    lon_array = np.array([15.,25.,35.])
    lat_out = np.array([2.5,7.5,12.5,17.5,22.5,27.5])
    lon_out = np.array([12.5,17.5,22.5,27.5,32.5,37.5])
    
    regridded_data = regrid2d(data_array, lat_array, lon_array, lat_out, lon_out)
    regridded_data_shape = regridded_data[0].shape
    
    old_mean = np.round(np.mean(data_array),0)
    new_mean = np.round(np.mean(regridded_data[0]),0)
    
    areas_old = areagrid(lat_array, lon_array)
    areas_new = areagrid(lat_out, lon_out)
    
    total_grid_old = data_array*areas_old
    total_domain_emission_old = np.round(np.sum(total_grid_old),0)
    total_grid_new = regridded_data[0]*areas_new
    total_domain_emission_new = np.round(np.sum(total_grid_new),0)
    
    assert regridded_data_shape == (len(lat_out), len(lon_out))
    assert old_mean == new_mean
    assert total_domain_emission_old == total_domain_emission_new


def test_acrg_regrid3d_netcdf_conservation():
    '''
    Test if the regridding function conserves the emission over the domain when
    using regrid3d
    '''
    data_array = np.array(([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]], [[0.1,0.2,0.3],[0.4,0.5,0.6],[0.7,0.8,0.9]]))
    lat_array = np.array([5.,15.,25.])
    lon_array = np.array([15.,25.,35.])
    time_array = np.array([1.,2.])
    lat_out = np.array([2.5,7.5,12.5,17.5,22.5,27.5])
    lon_out = np.array([12.5,17.5,22.5,27.5,32.5,37.5])
    
    data_array_sub = np.transpose(data_array, (1,2,0))

    regridded_data = regrid3d(data_array_sub, lat_array, lon_array, lat_out, lon_out, time_array)
    regridded_data_shape = regridded_data.shape
    
    old_mean = np.round(np.mean(data_array_sub),0)
    new_mean = np.round(np.mean(regridded_data),0)
    
    areas_old = areagrid(lat_array, lon_array)
    areas_new = areagrid(lat_out, lon_out)
    
    regridded_data_sub = np.transpose(regridded_data, (2,0,1))
    
    total_grid_old = data_array*areas_old
    total_domain_emission_old = np.round(np.sum(total_grid_old),0)
    total_grid_new = regridded_data_sub*areas_new
    total_domain_emission_new = np.round(np.sum(total_grid_new),0)
    
    assert regridded_data_shape == (len(lat_out), len(lon_out), len(time_array))
    assert old_mean == new_mean
    assert total_domain_emission_old == total_domain_emission_new
