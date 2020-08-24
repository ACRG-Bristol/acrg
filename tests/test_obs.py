"""
Created on Mon Jan 14 11:31:20 2019

@author: al18242
"""
import acrg_obs
import pytest
import numpy as np
import pandas as pd
import xarray as xr
import os
import shutil
import glob
from pathlib import Path

test_dir = Path("__file__").parent


def test_get_obs_structure():
    '''
    Test that the get_obs function returns the correct output structure
    '''
    start_date = "20160101"
    end_date = "20170101"
    sites = ["MHD", "GOSAT-UK"]
    
    recreated_data = acrg_obs.get_obs(sites, "CH4", start_date, end_date,
                                      data_directory="files/obs/",
                                      keep_missing=True, average=["1H", None],
                                      max_level=17)
    
    assert isinstance(recreated_data, dict)
    
    for site in sites:
        assert site in recreated_data
    
    
def test_get_obs_site():
    '''
    call get_obs on an example data file and test the properties are as expected
    
    Called with keep_missing=True and average=["1H"] - these are then used to check that
    the start and end indicies within the output are as expected
    '''
    start_date = "20160101"
    end_date = "20160201"
    recreated_data = acrg_obs.get_obs(["MHD"], "CH4", start_date, end_date,
                                      data_directory="files/obs/",
                                      keep_missing=True, average=["1H"])
    
    #test the date range is as expected
    assert np.amax(recreated_data["MHD"][0].time) == pd.to_datetime(end_date)-pd.Timedelta(hours=1)
    assert np.amin(recreated_data["MHD"][0].time) == pd.to_datetime(start_date)
    assert "mf" in recreated_data["MHD"].variables
    assert ("dmf" in recreated_data["MHD"].variables) or ("vmf" in recreated_data["MHD"].variables)

    
def test_get_obs_gosat():
    '''
    call get_obs on an example data file and test the properties are as expected
    
    Date indicies are checked to be within the correct bounds as satellite measurements are not continous
    '''
    start_date = "20160602"
    end_date = "20160604"
    recreated_data = acrg_obs.get_obs(["GOSAT-UK"], "CH4", start_date, end_date,
                                      data_directory="files/obs/",
                                      max_level = 17)
    
    #test the date range is as expected
    assert np.amax(recreated_data["GOSAT-UK"][0].time) < pd.to_datetime(end_date)
    assert np.amin(recreated_data["GOSAT-UK"][0].time) >= pd.to_datetime(start_date)
    assert "mf" in recreated_data["GOSAT-UK"][0].variables
    assert ("dmf" in recreated_data["GOSAT-UK"][0].variables) or ("vmf" in recreated_data["GOSAT-UK"].variables)
    

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
    
    assert "cfc113" in out.variables
    assert "time" in out.variables
    assert out.time.attrs["sampling_period_seconds"] == 60
    assert "seconds since" in out.time.encoding["units"]
    assert out.attrs["Calibration_scale"] == "TEST"
    assert out.attrs['station_long_name'] == u'Mace Head, Ireland'
    assert out.attrs['test'] == u'testing'
    assert out.cfc113.units == "1e-12"


def test_obs_process_gc():
    '''
    Test GC processing
    '''
    
    gc_files_directory = test_dir / "files/obs/GC"
    
    acrg_obs.process_gcwerks.gc("CGO", "medusa", "AGAGE",
                                input_directory = gc_files_directory,
                                output_directory = gc_files_directory,
                                version = "TEST")

    # Test if CGO directory has been created
    assert (gc_files_directory / "CGO").exists()
    
    # Check that enough files have been created
    assert len(list((gc_files_directory / "CGO").glob("*.nc"))) == 56
    
    # As an example, get CF4 data
    cf4_file = gc_files_directory / "CGO/AGAGE-GCMSMedusa_CGO_20180101_cf4-70m-TEST.nc"
    
    # Check if file exists
    assert cf4_file.exists()
    
    # Open dataset
    with xr.open_dataset(cf4_file) as f:
        ds = f.load()
    
    # Check a particular value (note that time stamp is 10 minutes before analysis time,
    # because in GCWerks files, times are at the beginning of the sampling period)
    assert np.allclose(ds.sel(time = slice("2018-01-01 04:33", "2018-01-01 04:35")).cf4.values,
                       np.array(83.546))

    assert np.allclose(ds.sel(time = slice("2018-01-20", "2018-01-20"))["cf4_repeatability"].values[0:1],
                       np.array(0.03679))

    # clean up
    shutil.rmtree(gc_files_directory / "CGO")
    
    
def test_obs_process_crds():
    '''
    Test CRDS file processing from GCWerks script
    
    '''

    crds_files_directory = test_dir / "files/obs/CRDS"
    
    acrg_obs.process_gcwerks.crds("BSD", "DECC",
                                  input_directory = crds_files_directory,
                                  output_directory = crds_files_directory,
                                  version = "TEST")

    # Test if BSD directory has been created
    assert (crds_files_directory / "BSD").exists()
    
    # Check that enough files have been created
    assert len(list((crds_files_directory / "BSD").glob("*.nc"))) == 3
    
    # As an example, get CH4 data
    ch4_file = crds_files_directory / "BSD/DECC-picarro_BSD_20140130_ch4-TEST.nc"
    
    # Check if file exists
    assert ch4_file.exists()
    
    # Open dataset
    with xr.open_dataset(ch4_file) as f:
        ds = f.load()
    
    # Check a particular value (note that time stamp is 10 minutes before analysis time,
    # because in GCWerks files, times are at the beginning of the sampling period)
    assert np.allclose(ds.sel(time = slice("2014-01-30 14:00:00", "2014-01-30 14:01:00")).ch4.values,
                       np.array(1953.88))
    assert np.allclose(ds.sel(time = slice("2014-01-30 14:00:00", "2014-01-30 14:01:00"))["ch4_variability"].values,
                       np.array(0.398))

    # clean up
    shutil.rmtree(crds_files_directory / "BSD")
