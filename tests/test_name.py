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

import pytest
import os
import sys
import glob
import numpy as np
import xarray as xray
import pandas as pd
import pickle

import acrg.name.name as name
import acrg.obs.read as read
from acrg.config.paths import Paths


acrg_path = Paths.acrg
benchmarkdir = acrg_path / "tests/files/benchmark/"


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
def output_directory():
    ''' Define base directory containing output files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/benchmark/")
    return directory

@pytest.fixture()
def fs_mock(fs, fp_directory, flux_directory, bc_directory, basis_directory, bc_basis_directory):
    #add the real jsons to the fake file system:
    fs.add_real_file(os.path.join(acrg_path, "data/species_info.json"))
    fs.add_real_file(os.path.join(acrg_path, "data/site_info.json"))
    #create footprint files
    fs.create_file(os.path.join(fp_directory, "EUROPE", "MHD-10magl_EUROPE_201402.nc"))
    fs.create_file(os.path.join(fp_directory, "EUROPE", "MHD-10magl_EUROPE_201405.nc"))
    fs.create_file(os.path.join(fp_directory, "EUROPE", "TAC-100magl_EUROPE_201402.nc"))

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
def small_domain_param():
    ''' Define set of parameters to be used with tests on the SMALL-DOMAIN '''
    
    param = {}
    param["sites"] = ["MHD"]
    param["domain"] = "SMALL-DOMAIN"
    param["start"] = "2014-02-01"
    param["end"] = "2014-03-01"
    param["heights"] = ["10magl"]
    param["species"] = "ch4"
    param["basis_case"] = "simple-grid"
    
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

@pytest.fixture(scope="module")
def hitres_timeseries_benchmark_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path, 'tests', 'files', 'benchmark', 'HiTRes_timeseries_201801.nc')
    return filename

@pytest.fixture(scope="module")
def hitres_sens_benchmark_file():
    ''' Define benchmark bc file '''
    filename = os.path.join(acrg_path, 'tests', 'files', 'benchmark', 'HiTRes_sens_201801.nc')
    return filename

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
def test_filenames(fs_mock, filenames_param):
    '''
    Test filenames function can find files of appropriate naming structure.
    '''
    out = name.filenames(**filenames_param)
    assert out == [os.path.join(filenames_param["fp_directory"],"EUROPE/MHD-10magl_EUROPE_201402.nc")]
def test_filenames_noheight(fs_mock, filenames_param):
    '''
    Test filenames() function can find height if it is not specified.
    '''
    filenames_param["height"] = None
    out = name.filenames(**filenames_param)
    assert out == [os.path.join(filenames_param["fp_directory"],"EUROPE/MHD-10magl_EUROPE_201402.nc")]
    
def test_filenames_shortlived_notavailable(fs_mock, filenames_param):
    '''
    Test filenames() function refuses 30 day integrated footprints if species is shortlived
    '''
    filenames_param["species"] = "CHBr3"
    out = name.filenames(**filenames_param)
    assert len(out) == 0

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


@pytest.fixture()
def timeseries_periodic_1hr():
    #emulate a footprint timeseries
    times = pd.date_range("2018-01-01 00:00", "2018-02-01 00:00", freq='1H')
    values = np.ones_like(times, dtype='float')
    ds = xray.Dataset({'values': (['time'], values)},
                          coords={'time':times})
    return ds


@pytest.fixture()
def timeseries_periodic_30min():
    #emulate a faster periodic timeseries than footprints, perhaps for continuous data
    times = pd.date_range("2018-01-01", "2018-02-01", freq='30min')
    values = np.ones_like(times, dtype='float')
    ds = xray.Dataset({'values': (['time'], values)},
                          coords={'time':times})
    return ds


@pytest.fixture()
def timeseries_random():
    #a random timeseries, as may be produced by a satellite
    times = pd.to_datetime(["2018-01-01 12:35", "2018-01-01 12:55", "2018-01-05 09:33", "2018-01-17 16:01"])
    values = np.ones_like(times, dtype='float')
    ds = xray.Dataset({'values': (['time'], values)},
                          coords={'time':times})
    return ds


@pytest.fixture()
def timeseries_quasiperiodic_1day():
    #a timeseries that is approximately periodic, like a manually triggered flask sample
    times = pd.to_datetime(["2018-01-01 08:55", "2018-01-02 09:16", "2018-01-03 09:33", "2018-01-04 10:01"])
    values = np.ones_like(times, dtype='float')
    ds = xray.Dataset({'values': (['time'], values)},
                          coords={'time':times})
    return ds


def test_align_skip(timeseries_quasiperiodic_1day, timeseries_periodic_1hr):
    #test that aligning is skipped for certain platforms
    for platform in ("satellite", "flask"):
        ds, fp = name.align_datasets(timeseries_quasiperiodic_1day, timeseries_periodic_1hr, platform=platform)
        assert ds == timeseries_quasiperiodic_1day
        assert fp == timeseries_periodic_1hr
        

def test_align_periodic_datasets(timeseries_periodic_30min, timeseries_periodic_1hr):
    #test that two periodic timeseries will align correctly
    dsa = timeseries_periodic_30min
    dsb = timeseries_periodic_1hr
    
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
    
    
def test_combine_datasets_periodic(timeseries_periodic_30min, timeseries_periodic_1hr):
    dsa = timeseries_periodic_1hr
    dsa['a'] = xray.DataArray(np.arange(len(dsa.time)), dims=["time"])
    dsb = timeseries_periodic_1hr
    dsb['b'] = xray.DataArray(np.arange(len(dsb.time)), dims=["time"])
    ds_out = name.combine_datasets(dsa, dsb, method = "ffill")

    errors = []
    if not np.sum(ds_out.time != dsa.time) == 0:
        errors.append("Output index does not match dsa index")
    if not ('a' in ds_out.keys()) and ('b' in ds_out.keys()):
        errors.append("Variables missing from combined dataset")
    
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
    
    
def test_combine_datasets_random(timeseries_random, timeseries_periodic_1hr):
    dsa = timeseries_random
    dsa['a'] = xray.DataArray(np.arange(len(dsa.time)), dims=["time"])
    dsb = timeseries_periodic_1hr
    dsb['b'] = xray.DataArray(np.arange(len(dsb.time)), dims=["time"])
    ds_out = name.combine_datasets(dsa, dsb, method = "ffill")

    errors = []
    if not np.sum(ds_out.time != dsa.time) == 0:
        errors.append("Output index does not match dsa index")
    if not ('a' in ds_out.keys()) and ('b' in ds_out.keys()):
        errors.append("Variables missing from combined dataset")
    
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
    
    
def test_align_and_combine_quasiperiodic(timeseries_quasiperiodic_1day, timeseries_periodic_1hr):
    #Test for flask style data to be combined with regular hourly footprints
    dsa = timeseries_quasiperiodic_1day
    dsa['a'] = xray.DataArray(np.arange(len(dsa.time)), dims=["time"])
    dsb = timeseries_periodic_1hr
    dsb['b'] = xray.DataArray(np.arange(len(dsb.time)), dims=["time"])
    
    dsa_out, dsb_out = name.align_datasets(dsa, dsb, platform="flask")
    ds_out = name.combine_datasets(dsa_out, dsb_out, method = "ffill")
    
    expected_output_a =  dsa['a'].values
    expected_output_b = dsb['b'].sel({"time":dsa['a'].time.dt.floor("1H")}).values
    
    assert np.all(ds_out.a.values == expected_output_a)
    assert np.all(ds_out.b.values == expected_output_b)
    
    
def test_indexesMatch(timeseries_periodic_1hr, timeseries_quasiperiodic_1day):
    #test the index matching function for same and different indexes
    dsa2 = timeseries_periodic_1hr.copy() * 2.

    errors = []
    if not name.indexesMatch(timeseries_periodic_1hr, dsa2) == True:
        errors.append("Same indexes not seen as the same")
    if not name.indexesMatch(timeseries_periodic_1hr, timeseries_quasiperiodic_1day) == False:
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
    ''' 
    Define set of input parameters to use with footprint_data_merge() function. 
    Note: cannot specify data directory directly.
    '''
    
    input_param = {}
    input_param["sites"] = measurement_param_small["sites"]
    input_param["species"] = measurement_param_small["species"]
    input_param["start"] = measurement_param_small["start"]
    input_param["end"] = measurement_param_small["end"]
    # Can't specify data directory
    
    #measurement_data = read.get_obs(**input_param)
    time = pd.date_range(input_param["start"], input_param["end"], freq='2H')
    nt = len(time)
    obsdf = pd.DataFrame({"mf":np.arange(nt)*1000.,"mf_repeatability":np.arange(nt), "status_flag":np.zeros(nt)}, index=time)
    obsdf.index.name = 'time'
    obsds = xray.Dataset.from_dataframe(obsdf)
    obsds.attrs["inlet"] = "10m"
    obsds.mf.attrs["units"] = 1e-9
    obsds.attrs["scale"] = "testscale"
    obsds.attrs["species"] = "ch4"
    measurement_data =  {"MHD": [obsds]}

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
    time = pd.date_range(input_param["start"], input_param["end"], freq='0.5H')
    nt = len(time)
    obsdf = pd.DataFrame({"mf":np.random.rand(nt)*1000.,"mf_repeatability":10*np.random.rand(nt), "mf_prior_factor":10*np.random.rand(nt), "mf_prior_upper_level_factor":15*np.random.rand(nt)}, index=time)
    obsdf.index.name = 'time'
    #obsdf.max_level = input_param["max_level"]
    obsds = xray.Dataset.from_dataframe(obsdf)
    obsds.mf.attrs["units"] = 1e-9
    obsds.attrs["max_level"] = input_param["max_level"]
    obsds.attrs["species"] = input_param["species"]
    measurement_data_sat = {input_param["sites"]: [obsds]}
    
    return measurement_data_sat


@pytest.mark.basic
def test_fp_data_merge(data,measurement_param,fp_directory,flux_directory,bc_directory):
    '''
    Test footprints_data_merge() function (with one site).
    Compares the output of footprints_data_merge function against the benchmarked output.
    Output from footprints_data_merge is currently saved as a pickle, then reloaded here for
    comparison. 
    Check data variables within dataset for site.
    '''
    
    out = name.footprints_data_merge(data,domain=measurement_param["domain"],fp_directory=fp_directory,
                                     flux_directory=flux_directory,bc_directory=bc_directory,
                                     calc_bc=True,calc_timeseries=True)

    # How pickle file was created initially
    #   Before re-running this step, make sure all preceding tests pass
    # with open(benchmarkdir / "fp_data_merge_benchmark.pkl", "wb") as fpdm_file:
    #     pickle.dump(out, fpdm_file)

    with open(benchmarkdir / "fp_data_merge_benchmark.pkl", "rb") as fpdm_file:
        benchmark_out = pickle.load(fpdm_file)
    
    for key, value in out.items():
        if isinstance(value, xray.Dataset):
            xray.testing.assert_allclose(benchmark_out[key], value)
        else:
            assert benchmark_out[key] == value

            
def test_fp_data_merge_sat(data_sat,measurement_param_sat,fp_directory,flux_directory,bc_directory):
    '''
    Test footprints_data_merge() function with GOSAT data.
    Check parameters within dictionary.
    Check data variables within dataset for site.
    '''
    site = measurement_param_sat["site"]
    expected_keys = [".species",".units",".bc",".flux",site]
    expected_data_var = ["mf","mf_repeatability","fp","particle_locations_n","particle_locations_e",
                         "particle_locations_s","particle_locations_w","bc"]
    
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

@pytest.fixture(scope="module")
def fp_data_merge_small_domain(data,small_domain_param,fp_directory,flux_directory,bc_directory):
    ''' Create fp_and_data dictionary from footprints_data_merge() function '''
    out = name.footprints_data_merge(data,domain=small_domain_param["domain"],fp_directory=fp_directory,
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

@pytest.fixture()
def fp_sensitivity_small_domain_param(fp_data_merge_small_domain,small_domain_param,basis_directory,flux_directory):
    ''' Define parameters for fp_sensitivity() function '''
    
    input_param = {}
    input_param["fp_and_data"] = fp_data_merge_small_domain
    input_param["domain"] = small_domain_param["domain"]
    input_param["basis_case"] = small_domain_param["basis_case"]
    input_param["basis_directory"] = basis_directory

    return input_param

def test_fp_sensitivity(fp_sensitivity_param):
    '''
    Test fp_sensitivity() function can create suitable output object from standard parameters.
    '''
    out = name.fp_sensitivity(**fp_sensitivity_param)
    
    assert out
    
def test_fp_sensitivity_H(small_domain_param,fp_sensitivity_small_domain_param):
    '''
    Test fp_sensitivity() function can create suitable output object from standard parameters.
    '''
    fp_data_H = name.fp_sensitivity(**fp_sensitivity_small_domain_param)
    lat_shape = fp_data_H['.flux']['all'].flux.values.shape[0]
    lon_shape = fp_data_H['.flux']['all'].flux.values.shape[1]
    time_shape = fp_data_H[small_domain_param["sites"][0]].fp.values.shape[2]
    
    H_all = fp_data_H[small_domain_param["sites"][0]].fp.values * fp_data_H['.flux']['all'].flux.values
    H_all_flat = H_all.reshape(lat_shape*lon_shape,time_shape)
    
    with xray.open_dataset(glob.glob(os.path.join(fp_sensitivity_small_domain_param["basis_directory"],
                                                  fp_sensitivity_small_domain_param["domain"],
                                                  fp_sensitivity_small_domain_param["basis_case"]+'*.nc'))[0]) as basis:
        
        basis_flat = np.ravel(np.squeeze(basis.basis.values))
    
    H = np.zeros(((int(np.max(basis_flat)),time_shape)))
    
    for val in range(int(np.max(basis_flat))):
        H[val,:] = np.sum(H_all_flat[np.where(basis_flat == val+1)[0],:],axis=0)
    
    assert np.all(H == fp_data_H[small_domain_param["sites"][0]].H.values) 

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
    

@pytest.fixture()
def dummy_timeseries_dict_gen():
    ''' '''
#     time = pd.date_range("2010-01-01", "2010-01-02", freq="1H")
#     data = np.arange(len(time))
#     release_lon = np.repeat(10., len(time))
#     testds = xray.Dataset({"mf" : (["time"], data), "release_lon":(["time"],release_lon)}, coords={"time":(["time"],time)})
#     dummy_timeseries_dict = {"TEST":testds}
    
    time = pd.date_range("2010-01-01", "2010-01-02", freq="1H")
    data = np.arange(len(time))
    release_lon = np.repeat(10., len(time))
    release_lat = np.repeat(10., len(time))
    lat = np.arange(10)+5
    lon = np.arange(10)+5
    fp = np.ones((10,10,len(time)))
    fp[3:8,3:8,:] = 0.0
    fp[5,5,10:] = 100.
    testds = xray.Dataset({"mf" : (["time"], data), 
                           "release_lon":(["time"],release_lon),
                           "release_lat":(["time"],release_lat),
                           "fp":(["lat","lon","time"],fp)}, 
                          coords={"time":(["time"],time),
                                  "lat":(["lat"],lat),
                                  "lon":(["lon"],lon)})
    dummy_timeseries_dict = {"TEST":testds}
    
    return dummy_timeseries_dict

def test_filtering_median(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "daily_median" filter
    '''
    filters = ["daily_median"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([11.5,24.])).all()

def test_filtering_daytime(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "daytime" filter
    '''
    filters = ["daytime"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([11, 12, 13, 14, 15])).all()

def test_filtering_nighttime(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "nighttime" filter
    '''
    filters = ["nighttime"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([0, 1, 2, 3, 23, 24])).all()

def test_filtering_noon(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "pblh_gt_threshold" filter
    '''
    filters = ["noon"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([12])).all()

def test_filtering_6hrmean(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "pblh_gt_threshold" filter
    '''
    filters = ["six_hr_mean"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([2.5, 8.5, 14.5, 20.5, 24.0])).all()

def test_filtering_local(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using "local_influence" filter
    '''
    filters = ["local_influence"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])).all()

@pytest.mark.long
def test_hitres_timeseries(hitres_timeseries_benchmark_file, hitres_sens_benchmark_file):
    '''
    Checks the output of name.timeseries_HiTRes matches a benchmarked version saved to
    HiTRes_timeseries_201801.nc
    '''
    
    with xray.open_dataset(hitres_timeseries_benchmark_file) as benchmark:
        benchmark.load()
    with xray.open_dataset(hitres_sens_benchmark_file) as benchmark_sens:
        benchmark_sens.load()
    
    fp_file_name = os.path.join(acrg_path, 'tests', 'files', 'LPDM', 'fp_NAME_minimal', 'EUROPE', 'WAO-20magl_UKV_co2_EUROPE_201801.nc')
    with xray.open_dataset(fp_file_name) as footprint:
        footprint.load()
    
    emiss_file_name = os.path.join(acrg_path, 'tests', 'files', 'LPDM', 'emissions', 'SMALL-DOMAIN', 'co2-ukghg-total-1hr_EUROPE_2018.nc')
    with xray.open_dataset(emiss_file_name) as emiss:
        emiss.load()
    flux_dict = {'high_freq': emiss}
    
    timeseries, sens = name.timeseries_HiTRes(flux_dict = flux_dict,
                                              fp_HiTRes_ds = footprint,
                                              output_TS = True,
                                              output_fpXflux = True,
                                              output_type = 'Dataset',
                                              output_file = None,
                                              verbose = False,
                                              time_resolution = '1H')
    
    assert np.allclose(timeseries['total'], benchmark['total'])
    assert np.allclose(sens['total'], benchmark_sens['total'])
    
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
def test_filtering_mult(dummy_timeseries_dict_gen):
    '''
    Test filtering() function can produce an output using multiple (non-conflicting) filters ("nighttime and 
    "local_influence")
    '''
    filters = ["local_influence","nighttime"]
    out = name.filtering(dummy_timeseries_dict_gen,filters)
    assert np.isclose(out["TEST"].mf.values, np.array([0, 1, 2, 3])).all()

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
