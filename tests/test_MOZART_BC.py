#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:35 2018

Test suite for MOZART_BC.py script.
Checks the outputs of MOZART_BC for known curtains


To run this test suite only from within the tests/ directory use the syntax
>> pytest test_MOZART_BC.py


@author: ag12733
"""

import pytest
import os
import sys
import xarray
import numpy as np
import acrg_MOZART_BC as bc


if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH")
    data_path = os.getenv("DATA_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg
    data_path = paths.data


@pytest.fixture(scope="module")
def mozart_bc_benchmark_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path,'tests/files/LPDM/bc/EUROPE/Benchmark_ch4_EUROPE_201403.nc')
    return filename

@pytest.fixture(scope="module")
def mozart_bc_output_directory():
    ''' Define benchmark bc file '''
    output_dir = os.path.join(acrg_path,'tests/files/')
    return output_dir

@pytest.fixture(scope="module")
def mozart_bc_output_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path,'tests/files/LPDM/bc/EUROPE/ch4_EUROPE_201403.nc')
    return filename


@pytest.mark.long
def test_mozart_bc_outputs(mozart_bc_benchmark_file, mozart_bc_output_directory, mozart_bc_output_file):
    ''' Checks the netcdf output of MOZART_BC.py matches a benchmarked ch4_EUROPE_201403_benchmark.nc file '''
        
    with xarray.open_dataset(mozart_bc_benchmark_file) as temp:
        benchmark = temp.load()
    
    bc.MOZART_BC_nc(start = '2014-03-01', end = '2014-05-01', species = 'ch4', domain = 'EUROPE', freq = 'M', runname = 'NewEDGAR', output_dir = mozart_bc_output_directory)
    
    with xarray.open_dataset(mozart_bc_output_file) as temp:
        output = temp.load()
    
    assert np.array_equal(output["vmr_e"], benchmark["vmr_e"])
    assert np.array_equal(output["vmr_w"], benchmark["vmr_w"])
    assert np.array_equal(output["vmr_n"], benchmark["vmr_n"])   
    assert np.array_equal(output["vmr_s"], benchmark["vmr_s"])
