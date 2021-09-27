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
import xarray
import numpy as np
import glob

import acrg.BC.MOZART_BC as bc
from acrg.config.paths import Paths

acrg_path = Paths.acrg


@pytest.fixture(scope="module")
def mozart_bc_benchmark_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path,'tests/files/LPDM/bc/EUROPE/Benchmark_ch4_EUROPE_201403.nc')
    return filename

@pytest.fixture(scope="module")
def mozart_bc_output_directory():
    ''' Define benchmark bc file '''
    output_dir = os.path.join(acrg_path,'tests/files/output/')
    return output_dir

@pytest.fixture(scope="module")
def mozart_bc_output_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path,'tests/files/output/LPDM/bc/EUROPE/ch4_EUROPE_201403.nc')
    return filename


@pytest.mark.long
@pytest.mark.skipif(not glob.glob(os.path.join(Paths.data,"LPDM/bc")), reason="No access to files in data_path")
def test_mozart_bc_outputs(mozart_bc_benchmark_file, mozart_bc_output_directory, mozart_bc_output_file):
    ''' Checks the netcdf output of MOZART_BC.py matches a benchmarked ch4_EUROPE_201403_benchmark.nc file '''
        
    with xarray.open_dataset(mozart_bc_benchmark_file) as temp:
        benchmark = temp.load()
    
    bc.MOZART_BC_nc(start = '2014-03-01', end = '2014-05-01', species = 'ch4', domain = 'EUROPE', freq = 'M', runname = 'NewEDGAR', output_dir = mozart_bc_output_directory)
    
    with xarray.open_dataset(mozart_bc_output_file) as temp:
        output = temp.load()
    
    assert np.allclose(output["vmr_e"], benchmark["vmr_e"])
    assert np.allclose(output["vmr_w"], benchmark["vmr_w"])
    assert np.allclose(output["vmr_n"], benchmark["vmr_n"])   
    assert np.allclose(output["vmr_s"], benchmark["vmr_s"])
