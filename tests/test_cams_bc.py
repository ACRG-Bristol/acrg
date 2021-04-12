#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:35 2018

Test suite for MOZART_BC.py script.
Checks the outputs of MOZART_BC for known curtains


To run this test suite only from within the tests/ directory use the syntax
>> pytest test_MOZART_BC.py


@author: vf20487
"""

import pytest
import os
import glob
import copy
import numpy as np
import xarray as xr
from acrg_name.name import open_ds
from acrg_BC import BC_CAMS as bc
from acrg_config.paths import paths

acrg_path = paths.acrg

@pytest.fixture(scope="module")
def bc_benchmark():
    '''
    Define benchmark bc file
    
    vmr_n, vmr_s, vmr_e, and vmr_w are 2D arrays of ones (height and lat or lon)
    '''
    filename = os.path.join(acrg_path,'tests','files','LPDM','bc','EUROPE','Benchmark_ones_EUROPE_201603.nc')
    with xr.open_dataset(filename) as temp:
        benchmark = temp.load()
    return benchmark
    
@pytest.fixture(scope="module")
def bc_output_directory():
    '''
    Define output directory
    '''
    output_dir = os.path.join(acrg_path,'tests','files', 'LPDM', 'bc', 'EUROPE')
    return output_dir

@pytest.fixture(scope="module")
def bc_output_file():
    '''
    Define output filename
    '''
    filename = os.path.join(acrg_path,'tests','files','LPDM','bc','EUROPE','ch4_EUROPE_201603.nc')
    return filename

@pytest.fixture(scope="module")
def bc_output_climatology_file():
    '''
    Define output filename
    '''
    filename = os.path.join(acrg_path,'tests','files','LPDM','bc','EUROPE','ch4_EUROPE_201703_climatology.nc')
    return filename
    
@pytest.mark.long
def test_bc_outputs(bc_benchmark, bc_output_directory, bc_output_file):
    '''
    Use a dummy CAMS inversion file, which has vmr=1 and coordinates time, height, lat, & lon
    Create a boundary conditions file which should have contain 2D arrays of ones for 
    vmr_n, vmr_s, vmr_e, and vmr_w
    
    Compare to the benchmark file defined above
    '''
    cams_directory = bc_output_directory
    bc.makeCAMSBC(start          = '2016-03-01',
                  end            = '2016-04-01',
                  species        = 'ch4',
                  domain         = 'EUROPE', 
                  outdir         = bc_output_directory,
                  cams_directory = bc_output_directory,
                  cams_version     = 'latest',
                  overwrite      = True,
                  test           = True)

    with xr.open_dataset(bc_output_file) as temp:
        output = temp.load()
    
    os.remove(bc_output_file)

    assert np.array_equal(output["vmr_e"], bc_benchmark["vmr_e"])
    assert np.array_equal(output["vmr_w"], bc_benchmark["vmr_w"])
    assert np.array_equal(output["vmr_n"], bc_benchmark["vmr_n"])   
    assert np.array_equal(output["vmr_s"], bc_benchmark["vmr_s"])
    
@pytest.mark.long
def test_bc_climatology_outputs(bc_benchmark, bc_output_directory, bc_output_climatology_file):
    '''
    Use a dummy CAMS inversion file, which has vmr=1 and coordinates time, height, lat, & lon
    Create a boundary conditions file which should have contain 2D arrays of ones for 
    vmr_n, vmr_s, vmr_e, and vmr_w
    
    Compare to the benchmark file defined above
    '''
    cams_directory = bc_output_directory
    bc.makeCAMSBC(start            = '2017-03-01',
                  end              = '2017-04-01',
                  species          = 'ch4',
                  domain           = 'EUROPE', 
                  outdir           = bc_output_directory,
                  cams_directory   = bc_output_directory,
                  cams_version     = 'latest',
                  clim_start       = '2016-03-01',
                  clim_end         = '2016-04-01',
                  make_climatology = True,
                  overwrite        = True,
                  test             = True)
    
    with xr.open_dataset(bc_output_climatology_file) as temp:
        output = temp.load()
    
    os.remove(bc_output_climatology_file)

    assert np.array_equal(output["vmr_e"], bc_benchmark["vmr_e"])
    assert np.array_equal(output["vmr_w"], bc_benchmark["vmr_w"])
    assert np.array_equal(output["vmr_n"], bc_benchmark["vmr_n"])   
    assert np.array_equal(output["vmr_s"], bc_benchmark["vmr_s"])