#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:35 2018

Test suite for tdmcmc_inputs.py script.
Checks whether script is able to run using subprocess module.

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_tdmcmc_inputs.py

Some tests have an additional decorator:
    @pytest.mark.long
which marks these tests as being part of a 'long' group of tests which can be excluded.

To run all tests except those labelled 'long' use the syntax
>> pytest test_tdmcmc_inputs.py -m 'not long'

@author: rt17603
"""

import pytest
import os
import subprocess
import glob

from acrg.config.paths import Paths


acrg_path = Paths.acrg
tdmcmc_path = os.path.join(acrg_path,"acrg/tdmcmc")
test_config_path = os.path.join(acrg_path,"tests/files/config")

@pytest.fixture(scope="module")
def tdmcmc_input_file():
    ''' Define tdmcmc inputs python script '''
    filename = os.path.join(tdmcmc_path,'tdmcmc_inputs.py')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_config_file():
    ''' Define tdmcmc config file input to run tdmcmc code '''
    filename = os.path.join(test_config_path,'tdmcmc_input_run.ini')
    return filename

@pytest.mark.skipif(not glob.glob(tdmcmc_path+"/*.so"), reason="Compiled Fortran programme unavailable.")
@pytest.mark.long
def test_tdmcmc_inputs(tdmcmc_input_file,tdmcmc_config_file):
    ''' Check that tdmcmc_inputs.py can be run with a standard tdmcmc config file '''
    result = subprocess.call(["python",tdmcmc_input_file,"-c{}".format(tdmcmc_config_file)])
    assert result == 0

@pytest.mark.skipif(not glob.glob(tdmcmc_path+"/*.so"), reason="Compiled Fortran programme unavailable.")
@pytest.mark.long
def test_tdmcmc_inputs_command_line(tdmcmc_input_file,tdmcmc_config_file):
    ''' Check that tdmcmc_inputs.py can be run with a standard tdmcmc config file '''
    result = subprocess.call(["python",tdmcmc_input_file, "2013-03-01", "2013-04-01", "-c{}".format(tdmcmc_config_file)])
    assert result == 0
    
