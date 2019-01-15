#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:31:20 2019

@author: al18242
"""
import acrg_obs
import numpy as np
import pandas as pd
import xarray as xr
import os

test_dir = os.path.dirname(os.path.abspath(__file__))

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

def test_process_utils_attributes():
    '''
    Test the acrg_obs.utils.attributes function
    
    Just makes sure that the function is returning a dataset with a few select 
    things changed. Could make this more comprehensive.
    '''
    
    attributes_test_file = os.path.join(test_dir,
                                        "files/obs/process_attributes_input.nc")

    with xr.open_dataset(attributes_test_file) as ds:
        ds.load()
    
    out = acrg_obs.utils.attributes(ds, "CFC-113", "MHD",
                                   global_attributes = {"test": "testing"},
                                   units = "ppt",
                                   scale = "TEST",
                                   sampling_period = 60,
                                   date_range = ["2000-01-01", "2000-01-10"])

    assert "cfc113" in out.keys()
    assert "time" in out.keys()
    assert out.time.attrs["sampling_period_seconds"] == 60
    assert "seconds since" in out.time.encoding["units"]
    assert out.attrs["Calibration_scale"] == "TEST"
    assert out.attrs['station_long_name'] == u'Mace Head, Ireland'
    assert out.attrs['test'] == u'testing'
    assert out.cfc113.units == 1e-12

