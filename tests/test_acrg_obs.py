#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:23:51 2018

Test suite for acrg_obs module.

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_acrg_obs.py

Some tests have an additional decorator:
    @pytest.mark.basic
which marks these tests as being part of a "basic" subset to be run if a full test run is not required.

To run only the basic tests use the syntax
>> pytest test_acrg_obs.py -m basic

@author: rt17603
"""

import pytest

#import acrg_agage as agage
import acrg_obs.read as read
from builtins import str

@pytest.fixture(scope="module")
def measurement_param():
    ''' Define set of measurement parameters to be used throughout the test suite '''
    param = {}
    param["sites"] = ["MHD","TAC"]
    param["start_date"] = str("2014-02-01")
    param["end_date"] = str("2014-03-01")
    param["inlets"] = ["10magl","100magl"]
    param["species"] = "ch4"
    
    return param

@pytest.fixture(scope="module")
def sat_measurement_param():
    ''' Define set of measurement parameters for satellite data '''
    param = {}
    param["sites"] = ["GOSAT-INDIA"]
    param["start_date"] = str("2013-07-01")
    param["end_date"] = str("2013-08-01")
    param["inlets"] = [None]
    param["species"] = "ch4"
    param["max_level"] = 17
    
    return param
 

#%%
#----------------------------

# TODO: Decide on other tests (if necessary) for functions within acrg_obs suite:
    # is_number(s)
    # synonyms(search_string, info)
    # listsearch(varnames, species, species_info, label="alt")
    # file_search_and_split(search_string)
    # quadratic_sum
    # ukmo_flags(site, site_info):
    # get_file_list(site, species, start, end, inlet,network = None, instrument = None, data_directory=None)
    # get(site_in, species_in, start = "1900-01-01", end = "2020-01-01",
    #        inlet=None, baseline=False, average=None, keep_missing=False,
    #        network = None, instrument = None,
    #        status_flag_unflagged = [0], data_directory=None)
    # get_gosat(site, species, max_level, start = "1900-01-01", end = "2020-01-01")
    # get_obs(sites, species, start = "1900-01-01", end = "2020-01-01",
    #            inlet = None, baseline = False, average = None, keep_missing=False,
    #            network = None, instrument = None, status_flag_unflagged = None,
    #            max_level = None):


#%%

@pytest.fixture()
def get_obs_sites_param(measurement_param):
    ''' Define set of input parameters for ground based sites for the get_obs() function. 
    Based on measurement_param '''
    # CAN'T SPECIFY FOLDER TO ACCESS OBSERVATIONS
    # CAN'T SPECIFY MULTIPLE HEIGHTS FOR SITES, SO HAVE TO SET TO NONE
    
    input_param = {}
    input_param["sites"] = measurement_param["sites"]
    input_param["species"] = measurement_param["species"]
    input_param["start_date"] = measurement_param["start_date"]
    input_param["end_date"] = measurement_param["end_date"]
    input_param["inlet"] = None

    return input_param

@pytest.fixture()
def get_obs_satellite_param(sat_measurement_param):
    ''' Define set of input parameters for satellite data for the get_obs() function. 
    Based on sat_measurement_param '''
    # CAN'T SPECIFY FOLDER TO ACCESS OBSERVATIONS
    
    input_param = {}
    input_param["sites"] = sat_measurement_param["sites"]
    input_param["species"] = sat_measurement_param["species"]
    input_param["start_date"] = sat_measurement_param["start_date"]
    input_param["end_date"] = sat_measurement_param["end_date"]
    input_param["inlet"] = sat_measurement_param["inlets"][0]
    input_param["max_level"] = sat_measurement_param["max_level"]

    return input_param

@pytest.mark.basic
def test_get_obs(get_obs_sites_param):
    '''
    Test get_obs() function returns a output in expected format when ground-based sites are specified.
    '''
    sites = get_obs_sites_param["sites"]
    out = read.get_obs(**get_obs_sites_param)
    assert isinstance(out,dict)
    
    for site in sites:
        assert site in out

@pytest.mark.basic  
def test_get_obs_sat(get_obs_satellite_param):
    '''
    Test get_obs() function returns output in expected format when satellite (GOSAT) sites are specified.
    '''
    sites = get_obs_satellite_param["sites"]
    out = read.get_obs(**get_obs_satellite_param)
    assert isinstance(out,dict)
    
    for site in sites:
        assert site in out
    
def test_get_obs_sat_nolevel(get_obs_satellite_param) :
    '''
    Test get_obs() raises ValueError when max_level is not specified for a satellite (GOSAT) site.
    '''
    with pytest.raises(ValueError) as e_info:
        get_obs_satellite_param["max_level"] = None
        out = read.get_obs(**get_obs_satellite_param)
        #assert out is None


