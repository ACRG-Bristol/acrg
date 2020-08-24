#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 15:25:41 2018

Test suite for acrg_name.name module.

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_name.py

Some tests have an additional decorator:
    @pytest.mark.basic
OR
    @pytest.mark.long
which marks these tests as being part of a 'basic' subset to be run if a full test run is not required
or as part of a 'long' group of tests which can be excluded.

To run only the basic tests use the syntax
>> pytest test_name.py -m 'basic'

To run all tests except those labelled 'long' use the syntax
>> pytest test_name.py -m 'not long'

@author: rt17603
"""
from __future__ import division

from builtins import range
from past.utils import old_div
import pytest
import os
import sys
import glob
import numpy as np
import xarray as xray
import pandas as pd

import acrg_name.name as name
#import acrg_agage as agage
import acrg_obs.read as read

#%%

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg


@pytest.fixture(scope="module")
def fp_directory():
    ''' Define base directory containing footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/fp_NAME_minimal/")
    return directory

@pytest.fixture(scope="module")
def flux_directory():
    ''' Define base directory containing flux files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/emissions/")
    return directory

@pytest.fixture(scope="module")
def bc_directory():
    ''' Define base directory containing boundary condition files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/bc/")
    return directory

@pytest.fixture(scope="module")
def basis_directory():
    ''' Define base directory containing basis function files'''
    directory = os.path.join(acrg_path,"tests/files/LPDM/basis_functions/")
    return directory

@pytest.fixture(scope="module")
def bc_basis_directory():
    ''' Define base directory containing boundary condition basis function files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/bc_basis_functions/")
    return directory


@pytest.fixture(scope="module")
def measurement_param():
    ''' Define set of measurement parameters to be used throughout the test suite '''
    param = {}
    param["sites"] = ["MHD","TAC"]
    param["domain"] = "EUROPE"
    param["start"] = "2014-02-01"
    param["end"] = "2014-04-01"
    param["heights"] = ["10magl","100magl"]
    param["species"] = "ch4"
    
    return param

@pytest.fixture(scope="module")
def measurement_param_small():
    ''' Define set of measurement parameters to be used throughout the test suite '''
    param = {}
    param["sites"] = ["MHD"]
    param["domain"] = "EUROPE"
    param["start"] = "2014-02-01"
    param["end"] = "2014-03-01"
    param["height"] = "10magl"
    param["species"] = "ch4"
    
    return param

@pytest.fixture(scope="module")
def basis_function_param():
    ''' Define set of parameters for basis functions to be used throughout the test suite '''
    param = {}
    param["basis_case"] = "sub-transd"
    param["bc_basis_case"] = "NESW"
    
    return param

@pytest.fixture(scope="module")
def measurement_param_sat():
    ''' Define set of satellite measurement parameters to be used throughout the test suite '''
    param = {}
    param["site"] = "GOSAT-UK"
    param["domain"] = "EUROPE"
    param["start"] = "2014-02-20"
    param["end"] = "2014-02-21"
    param["height"] = "column"
    param["species"] = "ch4"
    param["max_level"] = 17
    
    return param

#%%
#------------------------------

@pytest.fixture()
def filenames_param(measurement_param,fp_directory):
    ''' Define set of input parameters for the filenames() function. Based on measurement_param '''
    
    input_param = {}
    input_param["site"] = measurement_param["sites"][0]
    input_param["domain"] = measurement_param["domain"]
    input_param["start"] = measurement_param["start"]
    input_param["end"] = measurement_param["end"]
    input_param["height"] = measurement_param["heights"][0]
    input_param["fp_directory"] = fp_directory

    return input_param

@pytest.mark.basic
def test_filenames(filenames_param):
    '''
    Test filenames function can find files of appropriate naming structure.
    Currently expect fp_directory/domain/site + "*" + "-" + height + "*" + domain + "*" + yearmonth + "*.nc"
    '''
    out = name.filenames(**filenames_param)
    assert out

def test_filenames_noheight(filenames_param):
    '''
    Test filenames() function can find height if it is not specified.
    Currently expect fp_directory/domain/site + "*" + "-" + height + "*" + domain + "*" + yearmonth + "*.nc"
    '''
    filenames_param["height"] = None
    out = name.filenames(**filenames_param)
    assert out

#%%
#----------------------------

@pytest.fixture()
def footprint_param(fp_directory,measurement_param):
    ''' Define set of input parameters for footprints() function. Based on measurement_param '''
    
    input_param = {}
    input_param["sitecode_or_filename"] = measurement_param["sites"][0]
    input_param["domain"] = measurement_param["domain"]
    input_param["start"] = measurement_param["start"]
    input_param["end"] = measurement_param["end"]
    input_param["height"] = measurement_param["heights"][0]
    input_param["species"] = measurement_param["species"]
    input_param["fp_directory"] = fp_directory
    
    return input_param

def test_footprints_from_file(fp_directory,measurement_param):
    '''
    Test dataset can be created from footprint filename with footprints() function
    '''
    domain = measurement_param["domain"]
    directory = os.path.join(fp_directory,domain)
    search_str = os.path.join(directory,"*.nc")
    fp_filename = glob.glob(search_str)[0]
    
    out = name.footprints(fp_filename)
    
    assert out

@pytest.mark.long
def test_footprints_from_site(footprint_param,flux_directory,bc_directory):
    '''
    Test dataset can be created from set of parameters with footprints() function.
    '''

    out = name.footprints(**footprint_param)
    
    assert out

#%%
#----------------------------

#read_netcdfs() - tested as part of following functions
#interp_times() - may need to add tests for this

@pytest.fixture()
def flux_param(flux_directory,measurement_param):
    ''' Define set of input parameters for flux() function. Based on measurement_param '''
    
    input_param = {}
    input_param["domain"] = measurement_param["domain"]
    input_param["species"] = measurement_param["species"]
    input_param["flux_directory"] = flux_directory
    
    return input_param

def test_flux(flux_param):
    '''
    Test dataset can be created by flux() function
    '''
    out = name.flux(**flux_param)
    assert out
 
#%%
#----------------------------

@pytest.fixture()
def bc_param(bc_directory,measurement_param):
    ''' Define set of input parameters for boundary_conditions() function. Based on measurement_param '''
    
    input_param = {}
    input_param["domain"] = measurement_param["domain"]
    input_param["species"] = measurement_param["species"]
    input_param["bc_directory"] = bc_directory
    
    return input_param

def test_boundary_conditions(bc_param):
    '''
    Test dataset can be created by boundary_conditions() function
    '''
    out = name.boundary_conditions(**bc_param)
    assert out

#%%
#----------------------------

@pytest.fixture()
def basis_param(basis_directory,measurement_param,basis_function_param):
    ''' Define set of input parameters for basis() function. Based on measurement_param '''
    
    input_param = {}
    input_param["domain"] = measurement_param["domain"]
    input_param["basis_case"] = basis_function_param["basis_case"]
    input_param["basis_directory"] = basis_directory
    
    return input_param

def test_basis(basis_param):
    '''
    Test dataset can be created by basis() function
    '''
    out = name.basis(**basis_param)
    assert out

#%%
#----------------------------

@pytest.fixture()
def bc_basis_param(bc_basis_directory,measurement_param,basis_function_param):
    ''' Define set of input parameters for basis_boundary_conditions() function. Based on measurement_param '''
    
    input_param = {}
    input_param["domain"] = measurement_param["domain"]
    input_param["basis_case"] = basis_function_param["bc_basis_case"]
    input_param["bc_basis_directory"] = bc_basis_directory
    
    return input_param

def test_bc_basis(bc_basis_param):
    '''
    Test dataset can be created by basis_boundary_conditions() function
    '''
    out = name.basis_boundary_conditions(**bc_basis_param)
    assert out

#%%
#----------------------------

@pytest.fixture()
def dsa():
    times = pd.date_range("2018-01-01", "2018-02-01", freq='4H')
    vals = np.ones_like(times, dtype='float')
    dsa = xray.Dataset({'vals_a': (['time'], vals)},
                          coords={'time':times})
    return dsa

@pytest.fixture()
def dsb():
    times = pd.date_range("2018-01-07", "2018-02-07", freq='2H')
    vals = np.ones_like(times, dtype='float')
    dsb = xray.Dataset({'vals_b': (['time'], vals)},
                          coords={'time':times})
    return dsb

@pytest.fixture()
def dsc():
    times = pd.date_range("2018-01-01", "2018-02-01", freq='4H')
    vals = np.ones_like(times, dtype='float')
    dsc = xray.Dataset({'vals_c': (['time'], vals)},
                          coords={'time':times})
    dsc.time.values[1] = dsc.time.values[1] + pd.Timedelta("1 second")
    return dsc

def test_combine_datasets(dsa, dsb, method = "nearest", tolerance = None):
    """
    Merge two datasets.   
    """
    ds_out = name.combine_datasets(dsa, dsb)
    
    errors = []
    if not np.sum(ds_out.time != dsa.time) == 0:
        errors.append("Output index does not match dsa index")
    if not ('vals_a' in ds_out.keys()) and ('vals_b' in ds_out.keys()):
        errors.append("Variables missing from combined dataset")
    
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
    
def test_align_datasets(dsa, dsb):
    dsa_out, dsb_out = name.align_datasets(dsa, dsb)
    
    dsa_period = dsa.time[1] - dsa.time[0]
    dsa_out_period = dsa_out.time[1] - dsa_out.time[0]
    dsb_period = dsb.time[1] - dsb.time[0]
    dsb_out_period = dsb_out.time[1] - dsb_out.time[0]
    errors = []
    if not (dsa_out.time[0] == dsb_out.time[0]):
        errors.append("Start time don't match")
    if not (dsa_out.time[-1] == dsb_out.time[-1]):
        errors.append("End time don't match")
    if not dsa_out_period == dsb_out_period:
        errors.append("Time periods don't match in outputs")
    if not dsa_out_period == max(dsa_period, dsb_period):
        errors.append("Output time period is not the courser period")
        
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
    
def test_indexesMatch(dsa, dsc):
    dsa2 = dsa.copy() * 2.
    
    errors = []
    if not name.indexesMatch(dsa, dsa2) == True:
        errors.append("Same indexes not seen as the same")
    if not name.indexesMatch(dsa, dsc) == False:
        errors.append("Different indexes not seen as different")
    assert not errors, "errors occured:\n{}".format("\n".join(errors))

#%%
# ----------------------------


#TODO:
#    Set up some tests for timeseries_HiTRes if needed

#def timeseries_HiTRes(fp_HiTRes_ds, domain, HiTRes_flux_name, Resid_flux_name,
#                      output_TS = True, output_fpXflux = True, flux_directory=flux_directory):
#    """
#    Compute flux * HiTRes footprints
#    
#    HiTRes footprints record the footprint at each 2 hour period back in time for the first 24 hours.
#    
#    Need a high time resolution (HiTRes) flux to multiply the first 24 hours back of footprints.
#    
#    Need a residual (Resid) flux to multiply the residual integrated footprint for the remainder of the 20 day period.
#    
#    Can output the timeseries (output_TS = True) and/or the sensitivity map (output_fpXflux = True)
#    """

#TODO:
#    Set up some tests for timeseries_boundary_conditions if needed

#def timeseries_boundary_conditions(ds):
#    """
#    Compute particle location * global model edges time series.
#    All that is required is that you input an xray
#    dataset with both the particle locations and vmr at domain edge fields present    
#    """ 

#%%
#----------------------------

@pytest.fixture(scope="module")
def data(measurement_param_small):
    ''' Define set of input parameters to use with footprint_data_merge() function. 
    Note: cannot specify data directory directly. '''
    input_param = {}
    input_param["sites"] = measurement_param_small["sites"]
    input_param["species"] = measurement_param_small["species"]
    input_param["start"] = measurement_param_small["start"]
    input_param["end"] = measurement_param_small["end"]
    # Can't specify data directory
    
    #measurement_data = read.get_obs(**input_param)
    time = pd.date_range(input_param["start"], input_param["end"], freq='2H')
    nt = len(time)
    obsdf = pd.DataFrame({"mf":np.random.rand(nt)*1000.,"dmf":np.random.rand(nt), "status_flag":np.zeros(nt)}, index=time)
    obsdf.index.name = 'time'
    obsds = xray.Dataset.from_dataframe(obsdf)
    obsds.attrs["inlet"] = "10m"
    obsds.mf.attrs["units"] = 1e-9
    obsds.attrs["scale"] = "testscale"
    obsds.attrs["species"] = "ch4"
    
    #measurement_data = {'.species' : 'ch4', '.units' : 1e-9, 'MHD' : obsdf}  
    measurement_data = {"MHD": [obsds]}
    

    return measurement_data

@pytest.fixture(scope="module")
def data_sat(measurement_param_sat):
    ''' Define set of input parameters to use with footprint_data_merge() function. Edited for use with satellite data.
    Note: cannot specify data directory directly. '''
    input_param = {}
    input_param["sites"] = measurement_param_sat["site"]
    input_param["species"] = measurement_param_sat["species"]
    input_param["start"] = measurement_param_sat["start"]
    input_param["end"] = measurement_param_sat["end"]
    input_param["max_level"] = measurement_param_sat["max_level"]
    # Can't specify data directory
    
    #measurement_data_sat = read.get_obs(**input_param)
    time = pd.date_range(input_param["start"], input_param["end"], freq='0.5H')
    nt = len(time)
    obsdf = pd.DataFrame({"mf":np.random.rand(nt)*1000.,"dmf":10*np.random.rand(nt), "mf_prior_factor":10*np.random.rand(nt), "mf_prior_upper_level_factor":15*np.random.rand(nt)}, index=time)
    obsdf.index.name = 'time'
    #obsdf.max_level = input_param["max_level"]
    obsds = xray.Dataset.from_dataframe(obsdf)
    obsds.mf.attrs["units"] = 1e-9
    obsds.attrs["max_level"] = input_param["max_level"]
    obsds.attrs["species"] = input_param["species"]
    measurement_data_sat = {input_param["sites"]: [obsds]}
    #measurement_data_sat = {'.species' : input_param["species"], '.units' : 1e-9, input_param["sites"] : obsdf}  
    
    return measurement_data_sat


@pytest.mark.basic
def test_fp_data_merge(data,measurement_param_small,fp_directory,flux_directory,bc_directory):
    '''
    Test footprints_data_merge() function (with one site).
    Check parameters within dictionary.
    Check data variables within dataset for site.
    '''
    site = measurement_param_small["sites"][0]
    expected_keys = [".species",".units",".flux",".bc",site]
    expected_data_var = ["mf","dmf","fp","particle_locations_n","particle_locations_e","particle_locations_s",
                         "particle_locations_w","bc"]
    
    out = name.footprints_data_merge(data,domain=measurement_param_small["domain"],fp_directory=fp_directory,
                                     flux_directory=flux_directory,bc_directory=bc_directory)
    ds = out[site]
    
    for key in expected_keys:
        assert key in out

    for data_var in expected_data_var:
        assert data_var in ds.data_vars

    return out

def test_fp_data_merge_long(data,measurement_param,fp_directory,flux_directory,bc_directory):
    '''
    Test footprints_data_merge() function (with one site).
    Check parameters within dictionary.
    Check data variables within dataset for site.
    '''
    site = measurement_param["sites"][0]
    expected_keys = [".species",".units",".bc",".flux",site]
    expected_data_var = ["mf","dmf","fp","particle_locations_n","particle_locations_e",
                         "particle_locations_s","particle_locations_w","bc"]#,"mf_mod"]
    
    out = name.footprints_data_merge(data,domain=measurement_param["domain"],fp_directory=fp_directory,
                                     flux_directory=flux_directory,bc_directory=bc_directory,
                                     calc_bc=True,calc_timeseries=True)
    ds = out[site]
    
    for key in expected_keys:
        assert key in out

    for data_var in expected_data_var:
        assert data_var in ds.data_vars

    return out

def test_fp_data_merge_sat(data_sat,measurement_param_sat,fp_directory,flux_directory,bc_directory):
    '''
    Test footprints_data_merge() function with GOSAT data.
    Check parameters within dictionary.
    Check data variables within dataset for site.
    '''
    site = measurement_param_sat["site"]
    expected_keys = [".species",".units",".bc",".flux",site]
    expected_data_var = ["mf","dmf","fp","particle_locations_n","particle_locations_e",
                         "particle_locations_s","particle_locations_w","bc"]#,"mf_mod"]
    
    out = name.footprints_data_merge(data=data_sat,domain=measurement_param_sat["domain"],fp_directory=fp_directory,
                                     flux_directory=flux_directory,bc_directory=bc_directory)
    
    ds = out[site]
    
    for key in expected_keys:
        assert key in out
        
    for data_var in expected_data_var:
        assert data_var in ds.data_vars
        
    return out
    
# TODO:
#    Add more complete tests for footprints_data_merge function as it is important to check the output

#%%
#----------------------------

@pytest.fixture(scope="module")
def fp_data_merge(data,measurement_param_small,fp_directory,flux_directory,bc_directory):
    ''' Create fp_and_data dictionary from footprints_data_merge() function '''
    out = name.footprints_data_merge(data,domain=measurement_param_small["domain"],fp_directory=fp_directory,
                                     flux_directory=flux_directory,bc_directory=bc_directory)
    return out

@pytest.fixture()
def fp_sensitivity_param(fp_data_merge,measurement_param_small,basis_function_param,basis_directory,flux_directory):
    ''' Define parameters for fp_sensitivity() function '''
    
    input_param = {}
    input_param["fp_and_data"] = fp_data_merge
    input_param["domain"] = measurement_param_small["domain"]
    input_param["basis_case"] = basis_function_param["basis_case"]
    input_param["basis_directory"] = basis_directory

    return input_param

def test_fp_sensitivity(fp_sensitivity_param):
    '''
    Test fp_sensitivity() function can create suitable output object from standard parameters.
    '''
    out = name.fp_sensitivity(**fp_sensitivity_param)
    assert out

#%%
#----------------------------

@pytest.fixture()
def bc_sensitivity_param(fp_data_merge,measurement_param_small,basis_function_param,bc_basis_directory):
    ''' Define parameters for fp_sensitivity() function '''
    
    input_param = {}
    input_param["fp_and_data"] = fp_data_merge
    input_param["domain"] = measurement_param_small["domain"]
    input_param["basis_case"] = basis_function_param["bc_basis_case"]
    input_param["bc_basis_directory"] = bc_basis_directory

    return input_param

def test_bc_sensitivity(bc_sensitivity_param):
    '''
    Test bc_sensitivity() function can create suitable output object from standard parameters.
    '''
    out = name.bc_sensitivity(**bc_sensitivity_param)
    assert out


# TODO:
#    Decide if tests are needed for merge_sensitivity function

#def merge_sensitivity(fp_data_H,
#                      out_filename = None,
#                      remove_nan = True):
#    """
#    Outputs y, y_site, y_time in a single array for all sites
#    (as opposed to a dictionary) and H and H_bc if present in dataset.
#    Assumes my default that NaN values in y will be removed.
#    """
#

#%%
#----------------------------

# LIST OF OPTIONS TAKEN FROM TDMCMC_INPUTS.PY SCRIPT
#; Options are:
#;  "daily_median": What it says
#;  "daytime": Only between 11:00 - 15:00 inclusive 
#;  "nighttime": Only b/w 23:00 - 03:00 inclusive
#;  "noon": Only 12:00 fp and obs used
#;  "pblh_gt_500":
#;  "pblh_gt_250": 
#;  "local_influence": Only keep times when localness is low
#;  "six_hr_mean":
## THESE OPTIONS ARE NOT INCLUDED IN THE FP_SENSITIVITY FUNCTION AT THE MOMENT
#;  "ferry_loc": GAUGE-FERRY specific - Used to filter out dodgy ferry locations
#;  "ferry_mf": GAUGE-FERRY specific - Used to filter out dodg ferry 
#;  "ferry_fp_zero": GAUGE-FERRY specific

@pytest.fixture()
def fp_data_H_merge(fp_data_merge,fp_sensitivity_param,bc_sensitivity_param):
    ''' '''
    fp_sensitivity_param["fp_and_data"] = fp_data_merge.copy()
    fp_data_H = name.fp_sensitivity(**fp_sensitivity_param)
    bc_sensitivity_param["fp_and_data"] = fp_data_H
    fp_data_H = name.bc_sensitivity(**bc_sensitivity_param)
    return fp_data_H

@pytest.fixture()
def fp_data_H_pblh_merge(fp_data_H_merge):
    ''' '''
    sites = [key for key in list(fp_data_H_merge.keys()) if key[0] != '.']
    fp_data_H_pblh = fp_data_H_merge.copy()
    for site in sites:
        fp_data_H_pblh = fp_data_H_pblh[site].assign(**{"pblh_threshold":500})
    return fp_data_H_pblh
    
def add_local_ratio(fp_data_H):
    #TODO: Why is this here?
    sites = [key for key in list(fp_data_H.keys()) if key[0] != '.']
    
    release_lons=np.zeros((len(sites)))
    release_lats=np.zeros((len(sites)))
    for si, site in enumerate(sites):
        release_lons[si]=fp_data_H[site].release_lon[0].values
        release_lats[si]=fp_data_H[site].release_lat[0].values
        dlon=fp_data_H[site].sub_lon[1].values-fp_data_H[site].sub_lon[0].values
        dlat=fp_data_H[site].sub_lat[1].values-fp_data_H[site].sub_lat[0].values
        local_sum=np.zeros((len(fp_data_H[site].mf)))
       
        for ti in range(len(fp_data_H[site].mf)):
            release_lon=fp_data_H[site].release_lon[ti].values
            release_lat=fp_data_H[site].release_lat[ti].values
            wh_rlon = np.where(abs(fp_data_H[site].sub_lon.values-release_lon) < dlon/2.)
            wh_rlat = np.where(abs(fp_data_H[site].sub_lat.values-release_lat) < dlat/2.)
            local_sum[ti] = old_div(np.sum(fp_data_H[site].sub_fp[
            wh_rlat[0][0]-2:wh_rlat[0][0]+3,wh_rlon[0][0]-2:wh_rlon[0][0]+3,ti].values),np.sum(
            fp_data_H[site].fp[:,:,ti].values))  
            
        local_ds = xray.Dataset({'local_ratio': (['time'], local_sum)},
                                        coords = {'time' : (fp_data_H[site].coords['time'])})
    
        fp_data_H[site] = fp_data_H[site].merge(local_ds)
    
    return fp_data_H

@pytest.fixture()
def fp_data_H_lr_merge(fp_data_H_merge):
    ''' '''
    fp_data_H_lr = add_local_ratio(fp_data_H_merge.copy())
    return fp_data_H_lr

def test_filtering_median(fp_data_H_merge):
    '''
    Test filtering() function can produce an output using "daily_median" filter
    '''
    filters = ["daily_median"]
    out = name.filtering(fp_data_H_merge,filters)
    assert out

def test_filtering_daytime(fp_data_H_merge):
    '''
    Test filtering() function can produce an output using "daytime" filter
    '''
    filters = ["daytime"]
    out = name.filtering(fp_data_H_merge,filters)
    assert out

def test_filtering_nighttime(fp_data_H_merge):
    '''
    Test filtering() function can produce an output using "nighttime" filter
    '''
    filters = ["nighttime"]
    out = name.filtering(fp_data_H_merge,filters)
    assert out

def test_filtering_noon(fp_data_H_merge):
    '''
    Test filtering() function can produce an output using "pblh_gt_threshold" filter
    '''
    filters = ["noon"]
    out = name.filtering(fp_data_H_merge,filters)
    assert out

def test_filtering_6hrmean(fp_data_H_merge):
    '''
    Test filtering() function can produce an output using "pblh_gt_threshold" filter
    '''
    filters = ["six_hr_mean"]
    out = name.filtering(fp_data_H_merge,filters)
    assert out

def test_filtering_local(fp_data_H_lr_merge):
    '''
    Test filtering() function can produce an output using "local_influence" filter
    '''
    filters = ["local_influence"]
    out = name.filtering(fp_data_H_lr_merge,filters)
    assert out

# TODO: 
#    Not working yet
#def test_filtering_pblh(fp_data_H_pblh_merge):
#    '''
#    Test filtering() function can produce an output using "noon" filter
#    '''
#    filters = ["pblh_gt_threshold"]
#    print "fp_data_H_pblh_merge",fp_data_H_pblh_merge
#    out = name.filtering(fp_data_H_pblh_merge,filters)
#    assert out

# TODO:
#    Not checked yet
#def test_filtering_local_lapse(fp_data_H_merge):
#    '''
#    Test filtering() function can produce an output using "local_lapse" filter
#    '''
#    filters = ["local_lapse"]
#    out = name.filtering(fp_data_H_merge,filters)
#    assert out

@pytest.mark.long
def test_filtering_mult(fp_data_H_lr_merge):
    '''
    Test filtering() function can produce an output using multiple (non-conflicting) filters ("nighttime and 
    "local_influence")
    '''
    filters = ["local_influence","nighttime"]
    out = name.filtering(fp_data_H_lr_merge,filters)
    assert out

#----------------------------

## SHOULD THE FOLLOWING PLOTTING FUNCTIONS/CLASSES BE SPLIT OUT INTO A DIFFERENT MODULE?

#TODO: Decide on any tests for 
#    plot_map_setup class
#    plot_map_zoom() function
#    plot() function
#    time_unique() function
#    plot3d() function
#    animate() function
#    get_country class
