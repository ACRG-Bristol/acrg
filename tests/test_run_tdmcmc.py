#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:35 2018

Test suite for run_tdmcmc.py script.
Checks the outputs of run_tdmcmc for known inputs.
Only checks the non-random components of the post_mcmc netcdf file

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_run_tdmcmc.py

Some tests have an additional decorator:
    @pytest.mark.long
which marks these tests as being part of a 'long' group of tests which can be excluded.

To run all tests except those labelled 'long' use the syntax
>> pytest test_tdmcmc_inputs.py -m 'not long'

@author: ag12733
"""

import pytest
import os
import sys
import subprocess
import xarray
import numpy as np

from acrg_config.paths import paths
acrg_path = paths.acrg

tdmcmc_path = os.path.join(acrg_path,"acrg_tdmcmc")
test_param_path = os.path.join(acrg_path,"tests/files/tdmcmc/")

@pytest.fixture(scope="module")
def tdmcmc_input_file():
    ''' Define tdmcmc inputs python script '''
    filename = os.path.join(tdmcmc_path,'tdmcmc_inputs.py')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_param_file():
    ''' Define tdmcmc param file input to run tdmcmc code '''
    filename = os.path.join(test_param_path,'param.ini')
    return filename

@pytest.fixture(scope="module")
def run_tdmcmc_benchmark_file():
    ''' Define benchmark post_mcmc.nc file '''
    filename = os.path.join(acrg_path,'tests/files/tdmcmc/output/output_AGAGE_ch4_2014-02-01_benchmark.nc')
    return filename

@pytest.fixture(scope="module")
def run_tdmcmc_output_file():
    ''' Define generated post_mcmc.nc file '''
    filename = os.path.join(acrg_path,'tests/files/output/output_AGAGE_ch4_2014-02-01.nc')
    return filename

@pytest.mark.skip(reason="Comparison no longer valid. Suspect that a new benchmark file with updated obs is needed.")
def test_run_tdmcmc_outputs(tdmcmc_input_file, tdmcmc_param_file, run_tdmcmc_benchmark_file, run_tdmcmc_output_file):
    ''' Checks the netcdf output of run_tdmcmc.py matches a benchmarked output_AGAGE_ch4_2014-02-01_benchmark.nc file '''
        
    with xarray.open_dataset(run_tdmcmc_benchmark_file) as temp:
        post_mcmc_benchmark = temp.load()
    
    result = subprocess.call(["python",tdmcmc_input_file,"-c{}".format(tdmcmc_param_file)])
    
    #first check run didn't fail:
    assert result == 0
    
    with xarray.open_dataset(run_tdmcmc_output_file) as temp:
        post_mcmc = temp.load()
    
    #assert np.array_equal(post_mcmc["h_v_all"], post_mcmc_benchmark["h_v_all"])
    assert np.isclose(post_mcmc["h_v_all"], post_mcmc_benchmark["h_v_all"]) # account for floating point differences
    assert np.array_equal(post_mcmc["y"], post_mcmc_benchmark["y"])
    assert np.array_equal(post_mcmc["sigma_measure"], post_mcmc_benchmark["sigma_measure"])    
    assert np.array_equal(post_mcmc["R_indices"], post_mcmc_benchmark["R_indices"])
    
    os.remove(run_tdmcmc_output_file)
