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
import numpy as np
import xarray as xray
import netCDF4 as nc
import glob
import pandas as pd
import os
from os.path import join
import sys
from acrg.config.paths import Paths
from acrg.name import basis_functions
from acrg.name import name


acrg_path = Paths.acrg

basis_functions.fields_file_path = join(acrg_path, 'tests/files/LPDM/fp_NAME/')
basis_functions.basis_dir = join(acrg_path, 'tests/files/output/')
basis_functions.bc_basis_dir = join(acrg_path,'tests/files/output/')

if not os.path.exists(basis_functions.basis_dir):
    os.makedirs(basis_functions.basis_dir)
output_folder_to_create = join(basis_functions.basis_dir, 'EUROPE/')
if not os.path.exists(output_folder_to_create):
    os.makedirs(output_folder_to_create)
output_folder_to_create = join(basis_functions.bc_basis_dir, 'EUROPE/')
if not os.path.exists(output_folder_to_create):
    os.makedirs(output_folder_to_create)

name.data_path = join(acrg_path, 'tests/files/')

def test_acrg_basis_transd_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    basis_functions.basis_transd(domain = "EUROPE", time = "2014-02-01", basis_case="sub-transd", sub_lon_min = -80., sub_lon_max = 20., sub_lat_min = 30., sub_lat_max = 70.)
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/LPDM/fp_NAME/EUROPE/MHD-10magl_EUROPE_201402.nc")
    with xray.open_dataset(filename_field) as temp:
        field_dataset = temp.load()

    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/sub-transd_EUROPE_2014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_transd = temp.load()
    basis_array = test_basis_transd.variables['basis'][:]
    basis_array_shape = basis_array.shape
   
    assert basis_array_shape == (lat_field_dim, lon_field_dim, 1)
    assert np.sum(np.remainder(basis_array,1.)) == 0.0
    assert np.min(basis_array) == 0.0
    assert np.max(basis_array) == 4.0
    
    os.remove(filename_basis)

def test_acrg_basis_transd_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    basis_functions.basis_transd(domain = "EUROPE", time = "2014-02-01", basis_case="sub-transd", sub_lon_min = -80., sub_lon_max = 20., sub_lat_min = 30., sub_lat_max = 70.)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/sub-transd_EUROPE_2014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_transd = temp.load()

    basis_array = test_basis_transd.variables['basis'][:]
    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/basis_functions/EUROPE/Benchmark_sub-transd_EUROPE_2014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
    benchmark_basis_array = benchmark_dataset.variables['basis'][:]
    
    assert np.all(basis_array == benchmark_basis_array)
    
    os.remove(filename_basis)

def test_acrg_basis_bc_blocks_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    basis_functions.basis_bc_blocks(domain = "EUROPE", basis_case='Test_NESW', time = "2014-02-01", vertical = 4)
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/LPDM/fp_NAME/EUROPE/MHD-10magl_EUROPE_201402.nc")
    with xray.open_dataset(filename_field) as temp:
        field_dataset = temp.load()

    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_NESW_EUROPE_2014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_block = temp.load()
    
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

    os.remove(filename_basis)
    
def test_acrg_basis_bc_blocks_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    basis_functions.basis_bc_blocks(domain = "EUROPE", basis_case='Test_NESW', time = "2014-02-01", vertical = 4)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_NESW_EUROPE_2014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_block = temp.load()
 
    bc_basis_n_array = test_basis_bc_block.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_block.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_block.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_block.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/EUROPE/Benchmark_NESW_EUROPE_2014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
 
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)

    os.remove(filename_basis)

def test_acrg_basis_bc_uniform_shape():
    '''
    Test if netcdf files called and saved correctly and that the elements of the basis function
    are integers, and the min and max are the correct values.
    '''
    basis_functions.basis_bc_uniform(domain = "EUROPE", time = "2014-02-01", basis_case='Test_uniform')
    ##Read in one of the field datasets
    filename_field = os.path.join(acrg_path,"tests/files/LPDM/fp_NAME/EUROPE/MHD-10magl_EUROPE_201402.nc")
    with xray.open_dataset(filename_field) as temp:
        field_dataset = temp.load()
    lon_field_dim = len(field_dataset['lon'][:])
    lat_field_dim = len(field_dataset['lat'][:])
    ## Read in the newly created netcdf file
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_uniform_EUROPE_022014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_block = temp.load()

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

    os.remove(filename_basis)

def test_acrg_basis_bc_uniform_output():
    '''
    Test if the basis function array corresponds to a benchmark output
    '''
    basis_functions.basis_bc_uniform(domain = "EUROPE", time = "2014-02-01", basis_case='Test_uniform')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_uniform_EUROPE_022014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_uniform = temp.load()

    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/EUROPE/Benchmark_uniform_EUROPE_022014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    
    assert np.all(bc_basis_n_array == benchmark_bc_basis_n_array)
    assert np.all(bc_basis_e_array == benchmark_bc_basis_e_array)
    assert np.all(bc_basis_s_array == benchmark_bc_basis_s_array)
    assert np.all(bc_basis_w_array == benchmark_bc_basis_w_array)
    
    os.remove(filename_basis)
    
def test_acrg_basis_bc_all_gradients():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    basis_functions.basis_bc_all_gradients(domain = "EUROPE", species = 'ch4', units = 'ppb', time = "2014-02-01", basis_case='Test_horiz-strat-grad')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_horiz-strat-grad_EUROPE_022014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_uniform = temp.load()
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/EUROPE/Benchmark_horiz-strat-grad_EUROPE_022014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
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

    os.remove(filename_basis)
    
def test_acrg_basis_bc_horiz_gradient():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    basis_functions.basis_bc_horiz_gradients(domain = "EUROPE", time = "2014-02-01", basis_case='Test_horiz-grad')
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/Test_horiz-grad_EUROPE_2014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_uniform = temp.load()
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/EUROPE/Benchmark_horiz-grad_EUROPE_2014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
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

    os.remove(filename_basis)
    
def test_acrg_basis_bc_pca():
    '''
    Test if the bc gradient array corresponds to a benchmark output
    '''
    basis_functions.basis_bc_pca(domain = "EUROPE", time = "2014-02-01", species ='ch4', basis_case='Test_pca', numregions = 4)
    filename_basis = os.path.join(acrg_path,"tests/files/output/EUROPE/ch4_Test_pca_EUROPE_022014.nc")
    with xray.open_dataset(filename_basis) as temp:
        test_basis_bc_uniform = temp.load()
    bc_basis_n_array = test_basis_bc_uniform.variables['bc_basis_n'][:]
    bc_basis_e_array = test_basis_bc_uniform.variables['bc_basis_e'][:]
    bc_basis_s_array = test_basis_bc_uniform.variables['bc_basis_s'][:]
    bc_basis_w_array = test_basis_bc_uniform.variables['bc_basis_w'][:]

    filename_benchmark = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/EUROPE/Benchmark_ch4_pca_EUROPE_022014.nc")
    with xray.open_dataset(filename_benchmark) as temp:
        benchmark_dataset = temp.load()
    benchmark_bc_basis_n_array = benchmark_dataset.variables['bc_basis_n'][:]
    benchmark_bc_basis_e_array = benchmark_dataset.variables['bc_basis_e'][:]
    benchmark_bc_basis_s_array = benchmark_dataset.variables['bc_basis_s'][:]
    benchmark_bc_basis_w_array = benchmark_dataset.variables['bc_basis_w'][:]

    assert bc_basis_n_array.shape == benchmark_bc_basis_n_array.shape
    assert bc_basis_e_array.shape == benchmark_bc_basis_e_array.shape
    assert bc_basis_s_array.shape == benchmark_bc_basis_s_array.shape
    assert bc_basis_w_array.shape == benchmark_bc_basis_w_array.shape
    
    os.remove(filename_basis)