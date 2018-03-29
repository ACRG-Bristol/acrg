#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:35 2018

@author: rt17603
"""

import pytest
import os
import subprocess

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

tdmcmc_path = os.path.join(acrg_path,"acrg_tdmcmc")
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

def test_tdmcmc_run(tdmcmc_input_file,tdmcmc_config_file):
    ''' Check that tdmcmc_inputs.py can be run with a standard tdmcmc config file '''
    result = subprocess.call(["python",tdmcmc_input_file,"-c {}".format(tdmcmc_config_file)])
    assert result == 0
    