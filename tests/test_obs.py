#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:31:20 2019

@author: al18242
"""
import acrg_obs
import numpy as np
import pandas as pd
    
    
def checkUnits(data):
    '''
    Check the data units are numbers as needed in footprints_data_merge
    
    Inputs:
        data - output from get_obs
    '''
    assert isinstance(data[".units"], int) or isinstance(data[".units"], float)
    
def test_get_obs_site():
    '''
    call get_obs on an example data file and test the properties are as expected
    '''
    #test_data = pickle.load("path/to/data")
    start_date = "20160101"
    end_date = "20160201"
    recreated_data = acrg_obs.get_obs(["MHD"], "CH4", start_date, end_date,
                                      data_directory="files/obs/",
                                      keep_missing=True, average=["1H"])
    
    checkUnits(recreated_data)
    #assert recreated_data[".species"] == test_data[".species"]
    
    #test the date range is as expected
    assert np.amax(recreated_data["MHD"].index) == pd.to_datetime(end_date)-pd.Timedelta(hours=1)
    assert np.amin(recreated_data["MHD"].index) == pd.to_datetime(start_date)
    assert "mf" in recreated_data["MHD"].columns.values
    assert ("dmf" in recreated_data["MHD"].columns.values) or ("vmf" in recreated_data["MHD"].columns.values)
    
def test_get_obs_gosat():
    '''
    '''
    start_date = "20160602"
    end_date = "20160604"
    recreated_data = acrg_obs.get_obs(["GOSAT-UK"], "CH4", start_date, end_date,
                                      data_directory="files/obs/",
                                      max_level = 17)
    
    checkUnits(recreated_data)
    #assert recreated_data[".species"] == test_data[".species"]
    
    #test the date range is as expected
    assert np.amax(recreated_data["GOSAT-UK"].index) <= pd.to_datetime(end_date)-pd.Timedelta(hours=1)
    assert np.amin(recreated_data["GOSAT-UK"].index) >= pd.to_datetime(start_date)
    assert "mf" in recreated_data["GOSAT-UK"].columns.values
    assert ("dmf" in recreated_data["GOSAT-UK"].columns.values) or ("vmf" in recreated_data["GOSAT-UK"].columns.values)

#def test_process_raw_obs():