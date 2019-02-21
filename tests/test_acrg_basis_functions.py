 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs 17 Jan 2019 

Tests for checking the basis function generation:
    basis_blocks
    basis_transd
    basis_bc_blocks
    basis_bc_uniform

Testing the shape of the basis grids, the integers containted within the grid, 
and against a benchmarked output.

@author mi19881
"""
from __future__ import print_function

import numpy as np
import xarray as xray
import netCDF4 as nc
import glob
import pandas as pd

import os
from os.path import join

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH") 

if acrg_path is None:
    acrg_path = os.getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")

if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
          export DATA_PATH=/path/to/data/directory/ and restart python terminal")
    
import acrg_name.basis_functions
acrg_name.basis_functions.fields_file_path = join(acrg_path, 'tests/files/NAME/fp/')
acrg_name.basis_functions.basis_dir = join(acrg_path, 'tests/files/output/')
acrg_name.basis_functions.bc_basis_dir = join(acrg_path,'tests/files/output/')
if not os.path.exists(acrg_name.basis_functions.basis_dir):
    os.makedirs(acrg_name.basis_functions.basis_dir)
output_folder_to_create = join(acrg_name.basis_functions.basis_dir, 'EUROPE/')
if not os.path.exists(output_folder_to_create):
    os.makedirs(output_folder_to_create)
output_folder_to_create = join(acrg_name.basis_functions.bc_basis_dir, 'EUROPE/')
if not os.path.exists(output_folder_to_create):
    os.makedirs(output_folder_to_create)
        
import acrg_name.name
acrg_name.name.data_path = join(acrg_path, 'tests/files/')

def test_acrg_basis_blocks_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    acrg_name.basis_functions.basis_blocks(domain = "EUROPE", time = "2014-02-01", blocksize = 5, basis_case="Test_block")
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/NAME/fp/EUROPE/MHD-10magl_EUROPE_201402.nc")
    field_dataset = nc.Dataset(filename_field)
    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_block_EUROPE_2014.nc")
    test_basis_block = nc.Dataset(filename_basis)
    basis_array = test_basis_block.variables['basis'][:]
    basis_array_shape = basis_array.shape
    a = basis_array_shape[0]/5
    b = basis_array_shape[1]/5
    
    assert basis_array_shape == (lat_field_dim, lon_field_dim, 1)
    assert np.sum(np.remainder(basis_array,1.)) == 0.0
    assert np.min(basis_array) == 1.0
    assert np.max(basis_array) == a*b

def test_acrg_basis_blocks_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_blocks(domain = "EUROPE", time = "2014-02-01", blocksize = 5, basis_case ="Test_block")
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_block_EUROPE_2014.nc")
    test_basis_block = nc.Dataset(filename_basis)
    basis_array = test_basis_block.variables['basis'][:]
    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/basis_functions/EUROPE/Benchmark_block_EUROPE_2014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_basis_array = benchmark_dataset.variables['basis'][:]
    
    assert np.all(basis_array == benchmark_basis_array)

def test_acrg_basis_transd_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    acrg_name.basis_functions.basis_transd(domain = "EUROPE", time = "2014-02-01", basis_case="sub-transd", sub_lon_min = -80., sub_lon_max = 20., sub_lat_min = 30., sub_lat_max = 70.)
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/NAME/fp/EUROPE/MHD-10magl_EUROPE_201402.nc")
    field_dataset = nc.Dataset(filename_field)
    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/sub-transd_EUROPE_2014.nc")
    test_basis_transd = nc.Dataset(filename_basis)
    basis_array = test_basis_transd.variables['basis'][:]
    basis_array_shape = basis_array.shape
   
    assert basis_array_shape == (lat_field_dim, lon_field_dim, 1)
    assert np.sum(np.remainder(basis_array,1.)) == 0.0
    assert np.min(basis_array) == 0.0
    assert np.max(basis_array) == 4.0

def test_acrg_basis_transd_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_transd(domain = "EUROPE", time = "2014-02-01", basis_case="sub-transd", sub_lon_min = -80., sub_lon_max = 20., sub_lat_min = 30., sub_lat_max = 70.)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/sub-transd_EUROPE_2014.nc")
    test_basis_transd = nc.Dataset(filename_basis)
    basis_array = test_basis_transd.variables['basis'][:]
    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/basis_functions/EUROPE/Benchmark_sub-transd_EUROPE_2014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_basis_array = benchmark_dataset.variables['basis'][:]
    
    assert np.all(basis_array == benchmark_basis_array)


def test_acrg_basis_bc_blocks_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    acrg_name.basis_functions.basis_bc_blocks(domain = "EUROPE", basis_case='Test_NESW', time = "2014-02-01", vertical = 4)
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/NAME/fp/EUROPE/MHD-10magl_EUROPE_201402.nc")
    field_dataset = nc.Dataset(filename_field)
    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_NESW_EUROPE_2014.nc")
    test_basis_bc_block = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_block.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_block.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_block.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_block.variables['bc_basis_w'][:]
    bc_basis_n_array_shape =  bc_basis_n_array.shape
    bc_basis_e_array_shape =  bc_basis_e_array.shape
    bc_basis_s_array_shape =  bc_basis_s_array.shape
    bc_basis_w_array_shape =  bc_basis_w_array.shape
    
    assert bc_basis_n_array_shape == (20, lon_field_dim, 4*4, 1)
    assert bc_basis_e_array_shape == (20, lat_field_dim, 4*4, 1)
    assert bc_basis_s_array_shape == (20, lon_field_dim, 4*4, 1)
    assert bc_basis_w_array_shape == (20, lat_field_dim, 4*4, 1)
    assert np.sum(bc_basis_n_array) == lon_field_dim*20
    assert np.sum(bc_basis_e_array) == lat_field_dim*20
    assert np.sum(bc_basis_s_array) == lon_field_dim*20
    assert np.sum(bc_basis_w_array) == lat_field_dim*20


def test_acrg_basis_bc_blocks_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_bc_blocks(domain = "EUROPE", basis_case='Test_NESW', time = "2014-02-01", vertical = 4)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_NESW_EUROPE_2014.nc")
    test_basis_bc_block = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_block.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_block.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_block.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_block.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/bc_basis_functions/EUROPE/Benchmark_NESW_EUROPE_2014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)
    

def test_acrg_basis_bc_uniform_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    acrg_name.basis_functions.basis_bc_uniform(domain = "EUROPE", time = "2014-02-01", basis_case='Test_uniform')
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/NAME/fp/EUROPE/MHD-10magl_EUROPE_201402.nc")
    field_dataset = nc.Dataset(filename_field)
    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_uniform_EUROPE_022014.nc")
    test_basis_bc_block = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_block.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_block.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_block.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_block.variables['bc_basis_w'][:]
    bc_basis_n_array_shape =  bc_basis_n_array.shape
    bc_basis_e_array_shape =  bc_basis_e_array.shape
    bc_basis_s_array_shape =  bc_basis_s_array.shape
    bc_basis_w_array_shape =  bc_basis_w_array.shape
    
    assert bc_basis_n_array_shape == (20, lon_field_dim, 1, 1)
    assert bc_basis_e_array_shape == (20, lat_field_dim, 1, 1)
    assert bc_basis_s_array_shape == (20, lon_field_dim, 1, 1)
    assert bc_basis_w_array_shape == (20, lat_field_dim, 1, 1)
    assert np.sum(bc_basis_n_array) == lon_field_dim*20
    assert np.sum(bc_basis_e_array) == lat_field_dim*20
    assert np.sum(bc_basis_s_array) == lon_field_dim*20
    assert np.sum(bc_basis_w_array) == lat_field_dim*20


def test_acrg_basis_bc_uniform_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_bc_uniform(domain = "EUROPE", time = "2014-02-01", basis_case='Test_uniform')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_uniform_EUROPE_022014.nc")
    test_basis_bc_uniform = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/bc_basis_functions/EUROPE/Benchmark_uniform_EUROPE_022014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)
    

def test_acrg_basis_bc_all_gradients():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_bc_all_gradients(domain = "EUROPE", species = 'ch4', units = 'ppb', time = "2014-02-01", basis_case='Test_horiz-strat-grad')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_horiz-strat-grad_EUROPE_022014.nc")
    test_basis_bc_uniform = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/bc_basis_functions/EUROPE/Benchmark_horiz-strat-grad_EUROPE_022014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    assert bc_basis_n_array.shape == benchmark_bc_basis_n_array.shape
    assert bc_basis_e_array.shape == benchmark_bc_basis_e_array.shape
    assert bc_basis_s_array.shape == benchmark_bc_basis_s_array.shape
    assert bc_basis_w_array.shape == benchmark_bc_basis_w_array.shape
    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)


def test_acrg_basis_bc_horiz_gradient():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_bc_horiz_gradients(domain = "EUROPE", time = "2014-02-01", basis_case='Test_horiz-grad')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_horiz-grad_EUROPE_2014.nc")
    test_basis_bc_uniform = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/bc_basis_functions/EUROPE/Benchmark_horiz-grad_EUROPE_2014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    assert bc_basis_n_array.shape == benchmark_bc_basis_n_array.shape
    assert bc_basis_e_array.shape == benchmark_bc_basis_e_array.shape
    assert bc_basis_s_array.shape == benchmark_bc_basis_s_array.shape
    assert bc_basis_w_array.shape == benchmark_bc_basis_w_array.shape
    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)


def test_acrg_basis_bc_pca():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    acrg_name.basis_functions.basis_bc_pca(domain = "EUROPE", time = "2014-02-01", species ='ch4', basis_case='Test_pca', numregions = 4)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/ch4_Test_pca_EUROPE_022014.nc")
    test_basis_bc_uniform = nc.Dataset(filename_basis)
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/NAME/bc_basis_functions/EUROPE/Benchmark_ch4_pca_EUROPE_022014.nc")
    benchmark_dataset = nc.Dataset(filename_benchmark)
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    assert bc_basis_n_array.shape == benchmark_bc_basis_n_array.shape
    assert bc_basis_e_array.shape == benchmark_bc_basis_e_array.shape
    assert bc_basis_s_array.shape == benchmark_bc_basis_s_array.shape
    assert bc_basis_w_array.shape == benchmark_bc_basis_w_array.shape
    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)

