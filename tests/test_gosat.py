#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:46:32 2018

Test cases for acrg_gosat module

@author: rt17603
"""
from __future__ import print_function

import pytest
import acrg_satellite.gosat as gosat
import os
import xarray as xray
import numpy as np

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

#%%
##### Create dummy data for small test case with known result #######

@pytest.fixture()
def ds_dimensions():
    dimension_1 = "time"
    dimension_2 = "level"
    
    return dimension_1,dimension_2

@pytest.fixture()
def ds_data_variables():
    
    dv = {"xch4":"xch4",
          "xch4_uncertainty":"xch4_uncertainty",
          "lat":"latitude",
          "lon":"longitude",
          "ak":"xch4_averaging_kernel",
          "qf_flag":"xch4_quality_flag"}
    
    return dv

@pytest.fixture()
def dummy_dataset(ds_dimensions,ds_data_variables):
    
    dimension_1,dimension_2 = ds_dimensions
    
    dim_1 = 5 # First order dimension - time
    dim_2 = 3  # Second order dimension - level
    
    data_dict = {}
    data_dict[ds_data_variables["xch4"]] = (dimension_1,np.arange(0.,dim_1,1))
    data_dict[ds_data_variables["xch4_uncertainty"]] = (dimension_1,np.zeros(dim_1)+0.1)
    data_dict[ds_data_variables["lat"]] = (dimension_1,np.linspace(30.,33,dim_1))
    data_dict[ds_data_variables["lon"]] = (dimension_1,np.linspace(70.,94.,dim_1))
    
    averaging_kernel = [range(1,dim_2+1,1)]*dim_1
    data_dict[ds_data_variables["ak"]] = ([dimension_1,dimension_2],np.array(averaging_kernel,dtype=np.float64))
    
    qf_name = ds_data_variables["qf_flag"]
    h = int(round(dim_1/2.0))
    #qf_0 = [0]*int(round(dim_1/2.0))
    #qf_1 = [1]*int(round(dim_1//2))
    qf_0 = [0]*h
    qf_1 = [1]*(dim_1-h)
    quality_flag = qf_0 + qf_1
    data_dict[qf_name] = (dimension_1,np.array(quality_flag))
    
    dt = [np.datetime64('2012-01-01T01:00:00','ns') + np.timedelta64(i,'h') for i in range(0,dim_1,1)] # Creating values on time axis
    datetimes = np.array(dt) 
    levels = np.arange(0,dim_2,1)
    
    coords = {dimension_1:datetimes,dimension_2:levels} # Setting up time and level coordinates for dataset
    ds_small = xray.Dataset(data_dict,coords=coords)
    
    return ds_small

######## Create new dataset with extra dimension ###############

@pytest.fixture()
def extra_dimension_ds(dummy_dataset):
    
    ds_small = dummy_dataset

    new_dv = "new"
    new_dim = "new_dim"
    new_coords = np.arange(0,10)
    dim_3 = 10
    da = xray.DataArray(np.arange(1,11),coords={new_dim:new_coords},dims={new_dim:dim_3})
    ds_extra = ds_small.assign(**{new_dv:da})
    
    return ds_extra,new_dim,new_dv

######## Create new dataset with inconsistent dimension ###############

@pytest.fixture()
def inconsistent_dimension_ds(ds_dimensions,dummy_dataset):

    dimension_1,dimension_2 = ds_dimensions
    ds_small = dummy_dataset

    datetimes = ds_small[dimension_1].values
    levels = ds_small[dimension_2].values
    dim_1 = ds_small.dims[dimension_1]
    dim_2 = ds_small.dims[dimension_2]
    
    dv = np.array([range(1,dim_1+1,1)]*dim_2)

    #da = xray.DataArray(dv,coords={dimension_1:datetimes,dimension_2:levels},dims={dimension_1:dim_1,dimension_2:dim_2})
    da = xray.DataArray(dv,coords={dimension_1:datetimes,dimension_2:levels},dims=(dimension_2,dimension_1))
    ds_inconsistent = ds_small.assign(extra=da)
    
    return ds_inconsistent

####### Create new dataset with extra error data variable ############

@pytest.fixture()
def extra_error_ds(ds_dimensions,dummy_dataset):
    
    dimension_1,dimension_2 = ds_dimensions
    ds_small = dummy_dataset
    dim_1 = ds_small.dims[dimension_1]
    new_dv="xch4_error"
    
    da_extra_err = xray.DataArray(np.zeros(dim_1)+0.5,coords={dimension_1:ds_small[dimension_1]},dims={dimension_1:dim_1})
    ds_extra_err = ds_small.assign(**{new_dv:da_extra_err})
    
    return ds_extra_err,new_dv

####### Create new dataset with exposure_id data variable ############

@pytest.fixture()
def str_dv_ds(ds_dimensions,dummy_dataset):
    
    dimension_1 = ds_dimensions[0]
    ds_small = dummy_dataset
    dim_1 = ds_small.dims[dimension_1]
    str_dv = "exposure_id"
    
    da_str_dv = xray.DataArray(np.array(["Name{}".format(i) for i in range(dim_1)]),coords={dimension_1:ds_small[dimension_1]},
                              dims={dimension_1:dim_1})
    ds_str_dv = ds_small.assign(**{str_dv:da_str_dv})
    
    return ds_str_dv,str_dv

####### Create new dataset with only one value the first axis ############

@pytest.fixture()
def single_first_axis_ds(ds_dimensions,ds_data_variables,dummy_dataset):
    
    dimension_1 = ds_dimensions[0]
    data_var = ds_data_variables["xch4"]
    ds_single = dummy_dataset.copy(deep=True)
    
    ds_single[data_var] = ds_single[data_var][np.array([0])]
    ds_single = ds_single.dropna(dimension_1)
    
    return ds_single

####### Create new dataset with a range of time values on axis ############

@pytest.fixture()
def time_ds():
    
    dimension = "time"
    date = '2012-01-01'
    datetimes = np.arange('{}T00:33:45'.format(date),'{}T12:03:45'.format(date),10,dtype="datetime64[m]")
    values = np.arange(len(datetimes))
    
    ds_time = xray.Dataset({"values":(dimension,values)},coords={dimension:datetimes})
    
    return ds_time,dimension,date
    
##### Create data from test files for following test cases #####

@pytest.fixture(scope="module")
def gosat_data_folder():
    ''' Define directory containing example gosat data. '''
    directory = os.path.join(acrg_path,"tests/files/data/obs_raw/GOSAT/CH4/")
    return directory

@pytest.fixture(scope="module")
def gosat_dataset_raw(gosat_data_folder):
    ''' Create dataset from gosat data with no processing. '''
    directory = gosat_data_folder
    filename = "ESACCI-GHG-L2-CH4-GOSAT-OCPR-20120101-fv7.2.nc"
    filename = os.path.join(directory,filename)
    ds_raw = xray.open_dataset(filename)
    return ds_raw
    
@pytest.fixture(scope="module")
def gosat_dataset(gosat_dataset_raw):
    ''' Create gosat dataset with coordinate values added (necessary for most additional processing steps). '''
    ds_raw = gosat_dataset_raw
    ds = gosat.gosat_add_coords(ds_raw)
    return ds

###### Tests for apply_filter ################
def test_str_filter(dummy_dataset):
    '''
    Test Exception is raised when a str is specified as a filter_array
    '''
    filter_array = 'SPAM'
    with pytest.raises(KeyError) as e_info:
        gosat.apply_filter(dummy_dataset,filter_array)
    #self.assertRaises("KeyError",gosat.apply_filter(ds,filter_array))
    print(e_info)
    
def test_out_of_range_filter(ds_dimensions,dummy_dataset):
    '''
    Test Exception is raised when indices within filter_array are out of range
    '''
    dim_apply = ds_dimensions[0]
    filter_array = np.arange(len(dummy_dataset[dim_apply]),len(dummy_dataset[dim_apply])+10,1)
    with pytest.raises(IndexError) as e_info:
        gosat.apply_filter(dummy_dataset,filter_array)
    #self.assertRaises("IndexError",gosat.apply_filter(ds,filter_array))
    print(e_info)

def test_first_dim_filter(ds_dimensions,extra_dimension_ds):
    '''
    Test dimensions of dataset are correct after applying filter to first dimension
    '''
    dimension_1,dimension_2 = ds_dimensions
    ds_extra,dimension_3 = extra_dimension_ds[0:2]
    
    filter_array = np.array([0,1,2,3])
    dim_2 = ds_extra.dims[dimension_2]
    dim_3 = ds_extra.dims[dimension_3]
    
    out = gosat.apply_filter(ds_extra,filter_array,dim_apply=dimension_1)    
    
    assert out[dimension_1].size == len(filter_array)
    assert out[dimension_2].size == dim_2
    assert out[dimension_3].size == dim_3

def test_other_filter(ds_dimensions,extra_dimension_ds):
    '''
    Test dimensions of dataset are correct when applying filter to a different, first order dimension
    '''
    dimension_1,dimension_2 = ds_dimensions
    ds_extra,dimension_3 = extra_dimension_ds[0:2]
    
    dim_1 = ds_extra.dims[dimension_1]
    dim_2 = ds_extra.dims[dimension_2]
    filter_array = np.array([0])
    
    out = gosat.apply_filter(ds_extra,filter_array,dim_apply=dimension_3)
    
    assert out[dimension_1].size == dim_1
    assert out[dimension_2].size == dim_2
    assert out[dimension_3].size == len(filter_array)

def test_empty_filter(ds_dimensions,extra_dimension_ds):
    '''
    Test dimensions of dataset are correct when applying a zero length filter_array (should return empty on time dimension - i.e. not entries to keep)
    '''
    dimension_1,dimension_2 = ds_dimensions
    ds_extra,dimension_3 = extra_dimension_ds[0:2]
    
    filter_array = []
    dim_2 = ds_extra.dims[dimension_2]
    dim_3 = ds_extra.dims[dimension_3]
    
    out = gosat.apply_filter(ds_extra,filter_array=filter_array,dim_apply=dimension_1)
    
    assert out[dimension_1].size == 0
    assert out[dimension_2].size == dim_2
    assert out[dimension_3].size == dim_3

def test_full_gosat_filter(ds_dimensions,gosat_dataset):
    '''
    Test dimensions of dataset are correct for full gosat dataset
    '''
    dim_apply = ds_dimensions[0]
    filter_array = np.array([0,1,2,3])
    
    out = gosat.apply_filter(gosat_dataset,filter_array,dim_apply)
    
    assert out.time.size == len(filter_array)


###### Tests for gosat_quality_filter #######

def test_no_quality_filter(gosat_dataset):
    '''
    Test Exception is raised when "xch4_quality_flag" column is not included within dataset
    '''
    qf_name = "xch4_quality_flag"
    ds_no_qf = gosat_dataset.copy(deep=True)
    ds_no_qf = ds_no_qf.drop(qf_name)
        
    with pytest.raises(KeyError) as e_info:
        gosat.gosat_quality_filter(ds_no_qf)
    print(e_info)

def test_full_gosat_qf(ds_dimensions,gosat_dataset):
    '''
    Test dimensions of dataset are correct after quality filter has been applied to full gosat dataset.
    '''
    qf_name = "xch4_quality_flag"
    dimension_1 = ds_dimensions[0]
    out = gosat.gosat_quality_filter(gosat_dataset)
    wh = np.where(gosat_dataset[qf_name] == 0)[0]
    
    assert out[dimension_1].size == len(wh)

###### Tests for gosat_mode_filter ##########

###### Tests for gosat_pressure_filter ######


###### Tests for latlon_filter ##############

@pytest.fixture(scope="module")
def coord_boundaries():
    lat_bounds = [6.0,36.5]
    lon_bounds = [68.0,93.0]
    return lat_bounds,lon_bounds

def test_no_lat(coord_boundaries,gosat_dataset):
    '''
    Test Exception is raised if latitude column is not within dataset
    '''
    lat_bounds,lon_bounds = coord_boundaries
    columns = ["wrong_col","longitude"]
    with pytest.raises(KeyError) as e_info:
        gosat.latlon_filter(gosat_dataset,lat_bounds,lon_bounds,columns=columns)
    print(e_info)

def test_no_lon(coord_boundaries,gosat_dataset):
    '''
    Test Exception is raised if longitude column is not within dataset
    '''
    lat_bounds,lon_bounds = coord_boundaries
    columns = ["latitude","wrong_col"]
    with pytest.raises(KeyError) as e_info:
        gosat.latlon_filter(gosat_dataset,lat_bounds,lon_bounds,columns=columns)
    print(e_info)

def test_latlon(coord_boundaries,ds_dimensions,gosat_dataset):
    '''
    Test dimensions of dataset are correct after lat-lon filter has been applied
    '''
    lat_bounds,lon_bounds = coord_boundaries
    dimension_1 = ds_dimensions[0]
    columns=["latitude","longitude"]
    out = gosat.latlon_filter(gosat_dataset,lat_bounds,lon_bounds,columns=columns)
    
    wh = np.where((out[columns[0]] >= lat_bounds[0]) & (out[columns[0]] < lat_bounds[1]) & (out[columns[1]] >= lon_bounds[0]) & (out[columns[1]] < lon_bounds[1]))[0]
    
    assert out[dimension_1].size == len(wh)

###### Tests for coord_order ################

def test_wrong_data_var(dummy_dataset):
    '''
    Test KeyError is raised when specified data variable is not within dataset
    '''
    with pytest.raises(KeyError) as e_info:
        gosat.coord_order(dummy_dataset,data_vars=["wrong_name"])

def test_check_inconsistent(ds_dimensions,inconsistent_dimension_ds):
    '''
    Test function check for consistency of the order of coordinates returns False when coord order is not consistent for all data variables.
    '''
    dimension_1,dimension_2 = ds_dimensions
    out,check = gosat.coord_order(inconsistent_dimension_ds,check_consistent=True)
    
    assert check == False
    assert out["1"] == [dimension_1,dimension_2]
    assert out["2"] == [dimension_2,dimension_1]

def test_check_consistent(ds_dimensions,dummy_dataset):
    '''
    Test function check for consistency of the order of coordinates returns True when coord order is consistent for all data variables.
    '''
    dimension_1,dimension_2 = ds_dimensions
    out,check = gosat.coord_order(dummy_dataset,check_consistent=True)
    
    assert check == True
    assert out["1"] == [dimension_1]
    assert out["2"] == [dimension_2]

def test_no_check(dummy_dataset):
    '''
    Test consistency check is not performed and no value returned when check_consistent option is not selected
    '''
    out = gosat.coord_order(dummy_dataset,check_consistent=False)
    
    assert not isinstance(out,tuple)

def test_data_var(extra_dimension_ds):
    '''
    Test only dimensions of specified data variable are returned
    '''
    ds_extra,dimension_3,new_dv = extra_dimension_ds
    out = gosat.coord_order(ds_extra,data_vars=[new_dv])

    assert out["1"] == [dimension_3]

def test_full_gosat_coord(gosat_dataset):
    '''
    Test coord_order function can be run successfully over full gosat dataset
    '''
    out = gosat.coord_order(gosat_dataset)

###### Tests for midpoint ###################

def test_midpoint():
    '''
    Test output of midpoint is consistent with expected value
    '''
    datetimes = [np.datetime64('2010-01-01T00:00:00'),np.datetime64('2010-01-02T00:00:00')]
    dt = gosat.midpoint(datetimes)
    
    dt_expect = np.datetime64('2010-01-01T00:00:00') + (np.datetime64('2010-01-02T00:00:00') - np.datetime64('2010-01-01T00:00:00'))/2.
    
    #assert dt == dt_expect
    assert np.abs((dt - dt_expect)/np.timedelta64(1,'s')) < 1

###### Tests for calc_mid ###################

def test_mid_wrong_name(dummy_dataset):
    '''
    Test Exception is raised if data variable name is not within dataset
    '''
    with pytest.raises(KeyError) as e_info:
        name = "SPAM"
        out = gosat.calc_mid(dummy_dataset,name)
    
def test_mid_datetime(ds_dimensions,dummy_dataset):
    '''
    Test correct format (np.datetime64) is returned for "time"
    '''
    
    name = ds_dimensions[0] # Should be "time" based dimension
    out = gosat.calc_mid(dummy_dataset,name)
    
    assert isinstance(out,np.datetime64)

def test_mid_float(ds_data_variables,dummy_dataset):
    '''
    Test correct format is returned for a float based variable ("xch4")
    '''
    name = ds_data_variables["xch4"]
    out = gosat.calc_mid(dummy_dataset,name)
    
    if isinstance(out,xray.core.dataarray.DataArray):
        assert out.dtype == float
        assert out.size == 1
    else:
        assert isinstance(out,float)

def test_mid_higherdim(ds_data_variables,dummy_dataset):
    '''
    Test higher dimensions are retained when calculating the mid point on first order dimension
    '''
    name = ds_data_variables["ak"]
    dim_2 = dummy_dataset[name].shape[1]
    out = gosat.calc_mid(dummy_dataset,name)
    
    # xarray drops top level dimension when this value is returned
    assert out.shape[0] == dim_2

def test_full_gosat_mid(ds_data_variables,gosat_dataset):
    '''
    Test calc_mid function can be successfully run over full gosat dataset
    '''
    name = ds_data_variables["xch4"]
    out = gosat.calc_mid(gosat_dataset,name)
    
    if isinstance(out,xray.core.dataarray.DataArray):
        assert out.dtype == np.float32 or out.dtype == np.float64
        assert out.size == 1
    else:
        assert isinstance(out,float)

###### Tests for associated_error_name ######

def test_wrong_error_ident(ds_data_variables,dummy_dataset):
    '''
    Test associated_error_name function will return None if variable cannot be found from stripping error_ident from error_name.
    In this case the new variable name will be the same as error_name and so will be rejected.
    '''
    error_name = ds_data_variables["xch4_uncertainty"]
    error_ident = "false"
    out = gosat.associated_error_name(dummy_dataset,error_name,error_ident)
    
    assert out == None

def test_no_data_variable(ds_data_variables,dummy_dataset):
    '''
    Test associated_error_name function will return None if variable cannot be found from stripping error_ident from error_name.
    In this case new variable name will not be the same as error_name but will be rejected as it is not a data variable which exists within the dataset.
    '''
    error_name = ds_data_variables["xch4_uncertainty"]
    error_ident = "certainty"
    out = gosat.associated_error_name(dummy_dataset,error_name,error_ident)
    
    assert out == None
    
def test_associated_error_name_found(ds_data_variables,extra_error_ds):
    '''
    Test name can be found using a different error_ident and error_name than default.
    '''
    name = ds_data_variables["xch4"]
    ds_extra_err,error_name = extra_error_ds
    error_ident = error_name.split('_')[-1]
    out = gosat.associated_error_name(ds_extra_err,error_name,error_ident)
   
    assert out == name 

###### Tests for calc_mid_err ###############

def test_full_gosat_error_mid(ds_data_variables,gosat_dataset):
    '''
    Test calc_mid function can be successfully run over full gosat dataset
    '''
    name = ds_data_variables["xch4"]
    error_name = ds_data_variables["xch4_uncertainty"]
    out = gosat.calc_mid_err(gosat_dataset,name,error_name)  

###### Tests for concat_str ##################

def test_concat_str(str_dv_ds):
    '''
    Test concat_str function returns the expected string concatenation from a data variable 
    within a dataset
    '''
    ds,name = str_dv_ds
    separator = ','
    values = ds[name].values
    
    out = gosat.concat_str(ds,name,separator)
    out_split = out.split(separator)
    
    assert out.find(separator) != -1
    assert len(out_split) == len(values)
    
def test_keyerror_concat_str(dummy_dataset):
    '''
    Test concat_str function returns a KeyError when the data variable name cannot be found
    '''
    ds_small = dummy_dataset
    name = "unknown"
    with pytest.raises(KeyError) as e_info:
        gosat.concat_str(ds_small,name)

###### Tests for distance_lat,distance_lon #

def test_compare_dist_lat_lon():
    '''
    Test distance_lat and distance_lon functions produce sensible values when outputs are compared to each other.
    i.e. > 0 in all cases, equal at the equator and not equal away from the equator.
    '''
    distances = np.arange(10,110,10)
    for distance in distances:
        out1 = gosat.distance_lat(distance)
        out2 = gosat.distance_lon(distance,0) # set to latitude = 0
        assert out1 > 0
        assert out2 > 0
        assert out1 == out2
    
    for distance in distances:
        out1 = gosat.distance_lat(distance)
        out2 = gosat.distance_lon(distance,10) # set to latitude = 10
        assert out1 > 0
        assert out2 > 0
        assert out1 != out2
  
def test_check_symmetry_dist_lon():
    '''
    Test distance_lon function is symmetrical around the equator for give latitude values.
    '''
    latitudes_below0 = np.arange(-80,0,10)
    latitudes_above0 = np.arange(10,90,10)[::-1]
    
    for above,below in zip(latitudes_above0,latitudes_below0):
        out1 = gosat.distance_lon(10,above)
        out2 = gosat.distance_lon(10,below)
        assert out1 > 0
        assert out2 > 0
        assert out1 == out2

###### Tests for calc_dlat,calc_dlon #######

def test_calc_dlat(ds_data_variables,dummy_dataset):
    '''
    Test calc_dlat function can return the correct value when a range of latitude values is contained within the dataset.
    i.e. should return the maximum distance (max(lat) - min(lat)) + sounding in degrees of latitude
    '''
    ds_small = dummy_dataset
    name = ds_data_variables["lat"]
    sounding = 10.5
    
    min_lat = np.min(dummy_dataset[name])
    max_lat = np.max(dummy_dataset[name])
    
    diff = max_lat - min_lat
    
    out = gosat.calc_dlat(ds_small,sounding=sounding,lat_column=name)
    dlat = gosat.distance_lat(sounding) + diff
    
    assert out > 0
    assert out == dlat

def test_calc_dlat_single(ds_data_variables,single_first_axis_ds):
    '''
    Test calc_dlat function can return the correct value when only one value exists in the time dimension.
    i.e. should just return the dlat equivalent to 10.5 km sounding
    '''
    name = ds_data_variables["lat"]
    sounding = 10.5
    
    out = gosat.calc_dlat(single_first_axis_ds,sounding=sounding,lat_column=name)
    dlat = gosat.distance_lat(sounding)
    
    assert out > 0
    assert out == dlat

def test_calc_dlat_sounding(ds_data_variables,single_first_axis_ds):
    '''
    Test calc_dlat function can return the correct value when a different sounding value is provided.
    '''
    name = ds_data_variables["lat"]
    sounding = 100
    
    out = gosat.calc_dlat(single_first_axis_ds,sounding=sounding,lat_column=name)
    dlat = gosat.distance_lat(sounding)
    
    assert out > 0
    assert out == dlat 

def test_calc_dlon(ds_data_variables,dummy_dataset):
    '''
    Test calc_dlon function can return the correct value when a range of longitude values is contained within the dataset.
    i.e. should return the maximum distance (max(lon) - min(lon)) + sounding in degrees of longitude
    '''
    ds_small = dummy_dataset
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    sounding = 10.5
    
    min_lon = np.min(ds_small[columns[1]])
    max_lon = np.max(ds_small[columns[1]])
    
    lat = ds_small[columns[0]]
    
    diff = max_lon - min_lon
    
    out = gosat.calc_dlon(ds_small,sounding=sounding,columns=columns)
    dlon = gosat.distance_lon(sounding,np.mean(lat)) + diff
    
    assert out > 0
    assert np.array_equal(out,dlon)

def test_calc_dlon_single(ds_data_variables,single_first_axis_ds):
    '''
    Test calc_dlon function can return the correct value when only one value exists in the time dimension.
    i.e. should just return the dlon equivalent to 10.5 km sounding at 0 degrees latitude
    '''
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    sounding = 10.5
    
    lat = single_first_axis_ds[columns[0]]
    
    out = gosat.calc_dlon(single_first_axis_ds,sounding=sounding,columns=columns)
    dlon = gosat.distance_lon(sounding,np.mean(lat))
    
    assert out > 0
    assert out == dlon
    

def test_calc_dlon_sounding(ds_data_variables,single_first_axis_ds):
    '''
    Test calc_dlon function can return the correct value when a different sounding value is provided.
    '''
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    sounding = 100
    
    lat = single_first_axis_ds[columns[0]]
    
    out = gosat.calc_dlon(single_first_axis_ds,sounding=sounding,columns=columns)
    dlon = gosat.distance_lon(sounding,np.mean(lat))
    
    assert out > 0
    assert out == dlon

###### Tests for mean_ds ####################

# Multiple uncertainty values in error_ident list
# Consistent check
# Full dataset check
# Returns None when given multiple first order coords

def test_mult_error_ident_mean(ds_data_variables,extra_error_ds):
    '''
    Test multiple error identifiers can be specified and identified when calculating means within the function.
    '''
    ds_extra_err,error_name = extra_error_ds
    name = ds_data_variables["xch4"]
    error_names = [ds_data_variables["xch4_uncertainty"],error_name]
    error_ident = [ename.split('_')[-1] for ename in error_names]
    
    mid_errors = []
    for ename in error_names:
        mid_error = np.sqrt(np.mean(ds_extra_err[ename])**2. + np.std(ds_extra_err[name])**2.)
        mid_errors.append(mid_error)

    out = gosat.mean_ds(ds_extra_err,error_ident)
    
    for ename,mid_error in zip(error_names,mid_errors):
        assert out[ename].values == mid_error.values

def test_inconsistent_mean_ds(inconsistent_dimension_ds):
    '''
    Test mean_ds will return None if coordinate order is inconsistent
    '''
    out = gosat.mean_ds(inconsistent_dimension_ds)
    
    assert out == None
    
def test_mult_first_order_mean_ds(extra_dimension_ds):
    '''
    Test mean_ds will return None if there are multiple first order coordinates
    '''
    ds_extra = extra_dimension_ds[0]
    out = gosat.mean_ds(ds_extra)
    
    assert out == None

###### Tests for gosat_add_coords ###########

def test_incorrect_gosat_add_coords(dummy_dataset):
    '''
    Test gosat_add_coords will return None if the dimensions are not recognised as m,n
    '''
    out = gosat.gosat_add_coords(dummy_dataset)
    
    assert out == None

def test_data_vars_gosat_add(ds_data_variables,gosat_dataset_raw):
    '''
    Test data variables can be specified within the gosat_add_coords
    '''
    data_vars = [ds_data_variables["xch4"]]
    
    out = gosat.gosat_add_coords(gosat_dataset_raw,data_vars=data_vars)
    out_data_vars = out.data_vars
    
    for in_var,out_var in zip(data_vars,out_data_vars):
        assert out_var == in_var
 
def test_add_coords(gosat_dataset_raw):
    '''
    Test add_coords function can be run successfully over a full gosat dataset
    '''
    out = gosat.gosat_add_coords(gosat_dataset_raw)
    
    assert out != None

###### Tests for create_zbin ################

def test_zbin_1bin(coord_boundaries,ds_data_variables,dummy_dataset):
    '''
    Test output is produced when coordinate bin is specified as a float (should apply to both lat and lon)
    '''
    lat_bounds,lon_bounds=coord_boundaries
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    coord_bin = 0.5
    out = gosat.create_zbin(dummy_dataset,lat_bounds,lon_bounds,columns=columns,coord_bin=coord_bin)
    assert out

def test_zbin_1binlist(coord_boundaries,ds_data_variables,dummy_dataset):
    '''
    Test output is produced when coordinate bin is specified as a one item list (should apply to both lat and lon)
    '''
    lat_bounds,lon_bounds=coord_boundaries
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    coord_bin = [0.5]
    out = gosat.create_zbin(dummy_dataset,lat_bounds,lon_bounds,columns=columns,coord_bin=coord_bin)
    assert out

def test_zbin_2bin(coord_boundaries,ds_data_variables,dummy_dataset):
    '''
    Test output is produced when coordinate bins are specified as a two item list with the first value applied to lat and the second to lon.
    '''
    lat_bounds,lon_bounds=coord_boundaries
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    coord_bin = [0.1,0.5]
    out = gosat.create_zbin(dummy_dataset,lat_bounds,lon_bounds,columns=columns,coord_bin=coord_bin)
    assert out

def test_zbin_3binlist(coord_boundaries,ds_data_variables,dummy_dataset):
    '''
    Test None is returned if list for coordinate bin is too long
    '''
    #with pytest.raises(KeyError) as e_info:
    lat_bounds,lon_bounds=coord_boundaries
    columns = [ds_data_variables["lat"],ds_data_variables["lon"]]
    coord_bin = [0.1,0.5,0.5]
    out = gosat.create_zbin(dummy_dataset,lat_bounds,lon_bounds,columns=columns,coord_bin=coord_bin)
    assert out == None

###### Tests for zbin_filter ################

###### Tests for binned_mean ################

###### Tests for extract_dates ##############

def test_extract_dates(ds_dimensions,dummy_dataset):
    '''
    Test dates can be extracted from dummy dataset using extract_dates() function.
    '''
    # Assuming dimension 1 is related to time
    time_dim = ds_dimensions[0]
    dates = gosat.extract_dates(dummy_dataset,time_dim)

def test_extract_set_dates(time_ds):
    '''
    Test correct date can be extracted for all time points when times are all on the same day.
    '''
    ds_time,dimension,date = time_ds
    
    out = gosat.extract_dates(ds_time,dimension)

    date_out = np.unique(out)

    assert np.array_equal(date_out,np.array([date]))
    
    

#def extract_dates(ds,dim="time"):
#    '''
#    The extract_dates function converts datetime objects into date strings for all values along a given dimension.
#    Note: dimension should contain an array of np.datetime64 objects.
#    Args:
#        ds  : xarray.Dataset object. Should contain dimension specified by dim
#        dim : string of dimension containing datetime objects. Default = "time"
#    Returns:
#        np.array (str) : dates as an array of strings
#    '''
#    dtype = 'M8[D]'
#    return ds[dim].values.astype(dtype).astype(str) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)

# extract_files
# ds_check_unique

###### Tests for name_pressure_file #########

def test_open_name_pressure_file():
    '''
    Test pressure file can be opened using name_pressure_file() function.
    '''
    filename = os.path.join(acrg_path,"tests/files/LPDM/surface_pressure/SOUTHAMERICA/Pressure_C1_20120101_1d.txt")
    out = gosat.name_pressure_file(filename)
    assert out != None

###### Tests for name_pressure ##############

def test_name_pressure():
    '''
    Test pressure file can be opened and extracted with name_pressure() function.
    '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/surface_pressure/SOUTHAMERICA")
    start_date = "2012-01-01"
    end_date = "2012-01-02"
    max_days = 1
    
    out = gosat.name_pressure(directory,start_date,end_date,max_days=max_days)
    
###### Tests for name_pressure_match #######

def test_name_match_diff_columns(gosat_dataset):
    '''
    Test different column names can be used when calling name_pressure_match() function suitable output
    is produced.
    '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/surface_pressure")
    pressure_domain="SOUTHAMERICA"
    
    dimension_1 = "time"
    
    ds_diff_col_name = gosat_dataset.copy(deep=True)
    ds_diff_col_name = ds_diff_col_name.rename({"latitude":"lat","longitude":"lon"})
    
    out = gosat.name_pressure_match(ds_diff_col_name,columns=["lat","lon","time"],
                                    pressure_base_dir=directory,pressure_domain=pressure_domain)
    
    assert len(out) == ds_diff_col_name.dims[dimension_1]
    
@pytest.fixture(scope="module")
def gosat_brazil_day_offset_dataset(gosat_data_folder):
    ''' Create a dataset with a date offset from the first of the month to check pressure values are
    being selected correctly from the first day. '''
    directory = gosat_data_folder
    filename = "ESACCI-GHG-L2-CH4-GOSAT-OCPR-20120105-fv7.2.nc"
    filename = os.path.join(directory,filename)
    ds_raw = xray.open_dataset(filename)
    ds = gosat.gosat_add_coords(ds_raw)
    
    lat_bounds = [-35.753,7.251] # UPDATE
    lon_bounds = [-75.984,-32.793] # UPDATE

    ds = gosat.latlon_filter(ds,lat_bounds,lon_bounds)
    
    return ds

def test_pressure_match_offset(gosat_brazil_day_offset_dataset,ds_dimensions,ds_data_variables):
    '''
    Test offset dates match the expected pressure values for two example times (first and last times 
    within data) when day_template option is used.
    Note: day_template option is used if max_days > 1, using the first date as a template for the change
    over the course of a day.
    '''
    
    directory = os.path.join(acrg_path,"tests/files/LPDM/surface_pressure/")
    pressure_domain="SOUTHAMERICA"
    full_directory = os.path.join(directory,pressure_domain)

    lat_column = ds_data_variables["lat"]
    lon_column = ds_data_variables["lon"]
    time_column = ds_dimensions[0]

    index_1 = 0
    index_2 = len(gosat_brazil_day_offset_dataset[time_column])-1

    filter_array = np.array([index_1,index_2])
    gosat_brazil_day_offset_dataset_cut = gosat.apply_filter(gosat_brazil_day_offset_dataset.copy(deep=True),filter_array=filter_array,dim_apply=time_column)
   
    day_template=True
    max_days=31
    pressure_convert=1.

    time = gosat_brazil_day_offset_dataset[time_column].values
    
    start_date = min(time).astype('M8[D]').astype(str)
    end_date = max(time).astype('M8[D]').astype(str)
    
    name_pressure = gosat.name_pressure(full_directory,start_date=start_date,end_date=end_date,
                                        max_days=max_days)
    
    time_add = [np.timedelta64(t - t.astype('M8[D]'),'s') for t in time[filter_array]]
    #name_start_date = np.datetime64("2012-01-01")
    name_start_date = max(name_pressure["time"]).astype('M8[D]')-np.timedelta64(1,'D')
    offset_time = [name_start_date+ta for ta in time_add]
    
    lon = gosat_brazil_day_offset_dataset_cut[lon_column].values
    lat = gosat_brazil_day_offset_dataset_cut[lat_column].values

    expected_values = [name_pressure.sel(method="nearest",**{lat_column:la,lon_column:lo,time_column:t})["surface_pressure"].values for lo,la,t in zip(lon,lat,offset_time)]
    expected_values = np.array(expected_values)

    out = gosat.name_pressure_match(gosat_brazil_day_offset_dataset_cut,pressure_base_dir=directory,
                                    pressure_domain=pressure_domain,
                                    day_template=day_template,max_days=max_days,
                                    pressure_convert=pressure_convert)
   
    assert np.array_equal(out,expected_values)

def test_pressure_match_offset_day(gosat_brazil_day_offset_dataset,ds_dimensions):
    '''
    Test one day can be successfully used as a template for all days of a month when looking
    at the pressure profile.
    '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/surface_pressure/")
    pressure_domain="SOUTHAMERICA"

    gosat_brazil_day_offset_dataset_dummy = gosat_brazil_day_offset_dataset.copy(deep=True)
    gosat_brazil_day_offset_dataset_dummy["time"] += np.timedelta64(1,'D')
    
    day_template=True
    max_days=31
    
    dimension_1 = ds_dimensions[0]
    time = gosat_brazil_day_offset_dataset[dimension_1].values

    out = gosat.name_pressure_match(gosat_brazil_day_offset_dataset,pressure_base_dir=directory,
                                    pressure_domain=pressure_domain,
                                    day_template=day_template,max_days=max_days)
    out_ref = gosat.name_pressure_match(gosat_brazil_day_offset_dataset_dummy,
                                        pressure_base_dir=directory,pressure_domain=pressure_domain,
                                        day_template=day_template,max_days=max_days)
    
#    for p,t,pr in zip(out,time,out_ref):
#        print t,p,pr
    
    assert np.array_equal(out,out_ref)


# name_topog_match
# name_pressure_filter
# midpoint_bounds
# define_pressure_levels

###### Tests for ds_check_internal_unique ###

@pytest.fixture(scope="module")
def gosat_dataset_repeat():
    ''' Create gosat dataset object from data with repeated entries (2014-12-05). '''
    directory = os.path.join(acrg_path,"tests/files/data/obs_raw/GOSAT/CH4/")
    filename = "ESACCI-GHG-L2-CH4-GOSAT-OCPR-20141205-fv7.2.nc"
    filename = os.path.join(directory,filename)
    ds_raw = xray.open_dataset(filename)
    ds = gosat.gosat_add_coords(ds_raw)
    return ds

def test_internal_unique_repeat(gosat_dataset_repeat,gosat_dataset):
    '''
    Test processed data from file containing repeated entries can be merged with 
    another gosat dataset after ds_check_internal_unique() function has been called.
    Note: Error would not raised when dataset is created but when xarray attempts
    to merge that dataset with another dataset.
    '''
    gosat_dataset_repeat = gosat_dataset_repeat.sortby("time")
    gosat_dataset_2 = gosat.ds_check_internal_unique(gosat_dataset_repeat)
    out = gosat_dataset.merge(gosat_dataset_2)
    


###### Tests for gosat_output_filename ######

# gosat_construct_output - probably tested in other function more

###### Tests for gosat_output ###############

###### Tests for gosat_process_file #########

###### Tests for gosat_process ##############


   