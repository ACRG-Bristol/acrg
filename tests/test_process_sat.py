#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 10:15:58 2020

@author: rt17603
"""

import pytest
import os
import glob
import gzip
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr

import acrg.name.process as process
from acrg.name.name import open_ds
from acrg.grid.areagrid import areagrid
from . import process_helper as ph

from acrg.config.paths import Paths


acrg_path = Paths.acrg

#%% Setting up directory fixtures

@pytest.fixture()
def sat_byday_dummy_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay/")
    return directory

@pytest.fixture()
def sat_byday_mf_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay_MF/")
    return directory

@pytest.fixture()
def sat_byday_conc_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay_Conc/")
    return directory

@pytest.fixture()
def sat_bypoint_directory():
    ''' Define base directory containing satellite footprint files with separated data points '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByPoint/")
    return directory

@pytest.fixture()
def benchmark_output_directory():
    ''' Define base directory containing benchmark files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/benchmark/")
    return directory

@pytest.fixture()
def sat_param():
    '''
    Define input parameters for satellite data to test
    Keys include:
        "domain" (str),"site" (str),"height" (str),"year" (int),"month" (int),"satellite" (bool), "upper_level" (bool), "max_level" (bool)
    '''
    param = {}
    param["domain"] = "SOUTHAMERICA"
    param["site"] = "GOSAT-BRAZIL"
    param["height"] = "column"
    #param["year"] = 2012
    param["year"] = 2013
    param["month"] = 1
    param["satellite"] = True
    param["upper_level"] = 17
    param["max_level"] = 17
    
    return param

@pytest.fixture()
def subfolder_sat(sat_param):
    ''' Define subfolder name for satellite data'''
    domain = sat_param["domain"]
    site = sat_param["site"]
    height = sat_param["height"]
    
    subfolder = ph.define_subfolder(domain,site,height)
    return subfolder

@pytest.fixture()
def sat_param_dummy(sat_param):
    '''
    Define input parameters for dummy satellite data (based on sat_param)
    
    Want simplified
     - upper_level
     - max_level
    values
    '''
    param = sat_param
    param["upper_level"] = 2
    param["max_level"] = 2
    return param

@pytest.fixture()
def subfolder_sat_dummy(sat_param_dummy):
    ''' Define subfolder name for satellite data'''
    domain = sat_param_dummy["domain"]
    site = sat_param_dummy["site"]
    height = sat_param_dummy["height"]
    
    subfolder = ph.define_subfolder(domain,site,height)
    return subfolder

@pytest.fixture()
def subfolder_sat_byday_dummy(sat_byday_dummy_directory,subfolder_sat_dummy):
    ''' Define full subfolder path for satellite data with points grouped by day'''
    subfolder = os.path.join(sat_byday_dummy_directory, subfolder_sat_dummy)  
    return subfolder

@pytest.fixture()
def subfolder_sat_byday_mf(sat_byday_mf_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data (explicit mf data) with points grouped by day'''
    subfolder = os.path.join(sat_byday_mf_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def subfolder_sat_byday_conc(sat_byday_conc_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data (explicit conc data) with points grouped by day'''
    subfolder = os.path.join(sat_byday_conc_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def subfolder_sat_bypoint(sat_bypoint_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data with separate points'''
    subfolder = os.path.join(sat_bypoint_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def folder_names():
    '''
    Define subfolders within NAME output directory structure
    Keys includes:
        "fields_folder","particles_folder","met_folder","processed_folder","obs_folder"
    '''
    param = {}
    param["conc_folder"] = "Fields_files"
    param["mixr_folder"] = "MixR_files"
    param["particles_folder"] = "Particle_files"
    param["met_folder"] = "Met"
    param["processed_folder"] = "Processed_Fields_files"
    param["obs_folder"] = "Observations"

    return param  

#%% Testing read_fields() function

@pytest.fixture()
def get_fields_files_sat_dummy_byday(subfolder_sat_byday_dummy,folder_names,sat_param_dummy):
    '''
    Get filenames of fields files for satellite run with points grouped by day.
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = ph.create_datestr(sat_param_dummy)
    fields_folder = folder_names["mixr_folder"]
    
    fields_files = ph.get_fields_files(subfolder_sat_byday_dummy,fields_folder,datestr)
    return fields_files


@pytest.fixture()
def get_fields_files_sat_mf_byday(subfolder_sat_byday_mf,folder_names,sat_param):
    '''
    Get filenames of fields files for satellite run with points grouped by day.
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = ph.create_datestr(sat_param)
    fields_folder = folder_names["mixr_folder"]
    
    fields_files = ph.get_fields_files(subfolder_sat_byday_mf,fields_folder,datestr)
    
    return fields_files

@pytest.fixture()
def get_fields_files_sat_bypoint(subfolder_sat_bypoint,folder_names,sat_param):
    '''
    Get filenames of fields files for satellite run with separate points.
    Finds field files based on datestr (see create_datestr), multiple files per day
    (one for each point, filename syntax *YYYYMMDD-NNN*).
    '''
    datestr = ph.create_datestr(sat_param)
    fields_folder = folder_names["conc_folder"] # Previous fields folder name
    
    fields_files = ph.get_fields_files(subfolder_sat_bypoint,fields_folder,datestr)
    
    return fields_files

#@pytest.fixture()
#def get_fields_files_sat_byday_mf(subfolder_sat_byday_mf,folder_names,sat_param):
#    '''
#    Get filenames of fields files for satellite run with points grouped by day.
#    Finds field files based on datestr (see create_datestr), one file per day.
#    '''
#    datestr = ph.create_datestr(sat_param)
#    fields_folder = folder_names["fields_folder"]
#    
#    fields_files = ph.get_fields_files(subfolder_sat_byday_mf,fields_folder,datestr)
#    
#    return fields_files

@pytest.fixture()
def read_fields_file_sat_dummy_byday(get_fields_files_sat_dummy_byday):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_dummy_byday[0]
    header, column_headings, data_arrays, namever = process.read_file(fields_file)
    
    return header, column_headings, data_arrays, namever

@pytest.fixture()
def read_fields_file_sat_mf_byday(get_fields_files_sat_mf_byday):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_mf_byday[0]
    header, column_headings, data_arrays, namever  = process.read_file(fields_file)
    
    return header, column_headings, data_arrays, namever 

@pytest.fixture()
def read_fields_file_sat_bypoint(get_fields_files_sat_bypoint):
    ''' Read first satellite by point field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_bypoint[0]
    header, column_headings, data_arrays, namever = process.read_file(fields_file)
    
    return header, column_headings, data_arrays, namever

@pytest.mark.long
def test_read_field_file_sat_mf_byday(get_fields_files_sat_mf_byday):
    '''  Test read_file function can run for satellite data separated by day. '''
    fields_file = get_fields_files_sat_mf_byday[0]
    header, column_headings, data_arrays, namever = process.read_file(fields_file)

    assert header is not None
    assert column_headings is not None
    assert data_arrays is not None
    assert namever == 3

@pytest.mark.long
def test_read_field_file_sat_bypoint(get_fields_files_sat_bypoint):
    ''' Test read_file function can run for satellite data separated by point. '''
    point = 10
    fields_file = get_fields_files_sat_bypoint[point-1]
    header, column_headings, data_arrays, namever = process.read_file(fields_file)

    assert header is not None
    assert column_headings is not None
    assert data_arrays is not None
    assert namever == 3

#    point = 10
#    da_range = [(point-1)*17,(point-1)*17+17]
#    data_arrays_point = data_arrays[da_range[0]:da_range[1]]
#    for i,array in enumerate(data_arrays_point):
#        print i,np.min(array),np.max(array),np.mean(array),np.sum(array)
#    print "Over 17 elements of point {0}".format(point),np.mean(data_arrays_point)

@pytest.fixture()
def sat_byday_dummy_mixr():
    '''
    Defining dummy values for fields files. Should be consistent with
    the dummy fields file data which we can check against.

    Matching to values within:
      files/LPDM/raw_output/Satellite_ByDay/SOUTHAMERICA_GOSAT-BRAZIL_column/MixR_files/Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_synthetic_20130101.txt.gz

    Data include 4 fields entries for 2 lat, lon points.
    Overall grid is 100 x 100
    '''
    # Matching to values within:
    #  files/LPDM/raw_output/Satellite_ByDay/SOUTHAMERICA_GOSAT-BRAZIL_column/MixR_files/Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_synthetic_20130101.txt.gz

    # Expected header values
    header = {
    "Run name":                   "SOUTHAMERICA_GOSAT-BRAZIL_column_20130101_Synthetic",
    "Run time":                   "26/04/2019 09:30:38.040 UTC+01:00",
    "Met data":                   "NWP Flow.Global_Mk6_PT1 Flow; NWP Flow.Global_Mk6_PT2 Flow; NWP Flow.Global_Mk6_PT3 Flow; NWP Flow.Global_Mk6_PT4 Flow; NWP Flow.Global_Mk6_PT5 Flow; NWP Flow.Global_Mk6_PT6 Flow; NWP Flow.Global_Mk6_PT7 Flow; NWP Flow.Global_Mk6_PT8 Flow; NWP Flow.Global_Mk6_PT9 Flow; NWP Flow.Global_Mk6_PT10 Flow; NWP Flow.Global_Mk6_PT11 Flow; NWP Flow.Global_Mk6_PT12 Flow; NWP Flow.Global_Mk6_PT13 Flow; NWP Flow.Global_Mk6_PT14 Flow",
    "Start of release":           "01/01/2013 16:00 UTC",
    "End of release":             "01/01/2013 14:00 UTC",
    "Source strength":            "Multiple Sources",
    "Release location":           "Multiple Sources",
    "Release height":             "Multiple Sources",
    "Run duration":               "30day 1hr 46min",
    "X grid origin":                -90.00000,
    "Y grid origin":                -50.00000,  
    "X grid size":                         100,
    "Y grid size":                         100,
    "X grid resolution":            0.5000000,  
    "Y grid resolution":            1.0000000, 
    "Number of preliminary cols":            4,
    "Number of field cols":                  4,
    }

    # Expected column headings
    column_headings = {'category': ["INERT", "INERT", "INERT", "INERT"],
                       'name':  ["Req_H1Z1_1", "Req_H1Z1_2", "Req_H1Z1_3", "Req_H1Z1_4"],
                       'quantity':  ["Mixing Ratio", "Mixing Ratio", "Mixing Ratio", "Mixing Ratio"],
                       'species': ["INERT-TRACER", "INERT-TRACER", "INERT-TRACER", "INERT-TRACER"],
                       'unit': ["ppm s","ppm s","ppm s","ppm s"],
                       'source': ["SourceID_1", "SourceID_2", "SourceID_3", "SourceID_4"],     
                       'ensemble': ["No ensemble averaging", "No ensemble averaging", "No ensemble averaging", "No ensemble averaging"],
                       'time_averaging': ["30day 1hr 0min integral", "30day 1hr 0min integral", "30day 1hr 0min integral", "30day 1hr 0min integral"],
                       'horiz_averaging': ["No horizontal averaging", "No horizontal averaging", "No horizontal averaging", "No horizontal averaging"],
                       'vert_averaging': ["No vertical averaging", "No vertical averaging", "No vertical averaging", "No vertical averaging"],
                       'prob_percentile': ["",     "",     "",     ""],
                       'prob_percentile_ensemble': ["",     "",     "",     ""],
                       'prob_percentile_time': ["",     "",     "",     ""],
                       'time': ["02/12/2012 14:58 UTC",          "02/12/2012 14:58 UTC",          "02/12/2012 14:58 UTC",          "02/12/2012 14:58 UTC"],
                       'z_level': ["Z = 20.00000 m agl",            "Z = 20.00000 m agl",            "Z = 20.00000 m agl",            "Z = 20.00000 m agl"],
                       'D': ["",     "",     "",     ""],
                    }

    # Update time values to be in Python datetime format
    column_headings["time"] = [dt.datetime.strptime(t, '%d/%m/%Y  %H:%M %Z') for t in column_headings["time"]]

    # Add empty columns to the start of each entry
    starter = ['']*4 # 4 empty columns at the start for each row
    column_headings = {key:starter + value for key, value in column_headings.items()}

    # Create expected data_arrays object
    # Take expected values and expand them into a 2D nlat x nlon grid
    xgrid = header["X grid size"]
    ygrid = header["Y grid size"]
    grid = np.zeros((ygrid, xgrid), dtype=np.float32)
    data_arrays = [grid.copy(), grid.copy(), grid.copy(), grid.copy()]
    
    # Should match exactly to relevant values in "Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_synthetic_20130101.txt.gz"
    xindex = np.array([1, 2]) # X index
    yindex = np.array([2, 1]) # Y index
    fields = [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0]] # Field columns x 4

    xindex = xindex - 1 # namever3 index correction
    yindex = yindex - 1 # namever3 index correction

    for data_array, field in zip(data_arrays, fields):
        for y, x, val in zip(yindex, xindex, field):
            data_array[y, x] = val
    
    namever = 3

    return header, column_headings, data_arrays, namever


def test_read_file_sat_dummy_byday(get_fields_files_sat_dummy_byday,
                                   sat_byday_dummy_mixr):
    '''
    Test read_file function against expected outputs. Based on synthetic data.
    '''
    fields_file = get_fields_files_sat_dummy_byday[0]
    header, column_headings, data_arrays, namever = process.read_file(fields_file)

    exp_header, exp_columns_headings, exp_data_arrays, exp_namever = \
                                                        sat_byday_dummy_mixr

    assert header.keys() == exp_header.keys()

    for key, value in header.items():
        exp_value = exp_header[key]
        if isinstance(value, float):
            assert np.isclose(value, exp_value)
        else:
            assert value == exp_value
    
    assert column_headings.keys() == exp_columns_headings.keys()

    for key, value_list in column_headings.items():
        exp_value_list = exp_columns_headings[key]
        if isinstance(value_list[0], float):
            assert np.allclose(value_list, exp_value_list)
        else:
            assert value_list == exp_value_list

    for data_array, exp_data_array in zip(data_arrays, exp_data_arrays):
        assert np.allclose(data_array, exp_data_array)

    assert namever == exp_namever

#%% Test define grid function

@pytest.mark.long
def test_define_grid_sat_mf_byday(read_fields_file_sat_mf_byday,sat_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with points 
    grouped by day.
    '''
    header,column_headings,data_arrays,namever = read_fields_file_sat_mf_byday
    out = process.define_grid(namever,header,column_headings,satellite=True,upper_level=sat_param["upper_level"])
    
    assert out

@pytest.mark.long
def test_define_grid_sat_bypoint(read_fields_file_sat_bypoint,sat_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with separate 
    points.
    '''
    header,column_headings,data_arrays,namever = read_fields_file_sat_bypoint
    out = process.define_grid(namever,header,column_headings,satellite=True,upper_level=sat_param["upper_level"])

    assert out

#%% Test particle_locations()

@pytest.fixture()
def define_grid_sat_dummy_byday(read_fields_file_sat_dummy_byday,sat_param_dummy):
    ''' Create output from process.define_grid() function for one satellite by day field file '''
    header,column_headings,data_arrays,namever = read_fields_file_sat_dummy_byday
    upper_level=sat_param_dummy["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(namever,header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def define_grid_sat_mf_byday(read_fields_file_sat_mf_byday,sat_param):
    ''' Create output from process.define_grid() function for one satellite by day field file '''
    header,column_headings,data_arrays,namever = read_fields_file_sat_mf_byday
    upper_level=sat_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(namever,header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def define_grid_sat_bypoint(read_fields_file_sat_bypoint,sat_param):
    ''' Create output from process.define_grid() function for one satellite by point field file '''
    header,column_headings,data_arrays,namever = read_fields_file_sat_bypoint
    upper_level=sat_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(namever,header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def get_particle_files_sat_dummy_byday(subfolder_sat_byday_dummy,folder_names,
                                       get_fields_files_sat_dummy_byday):
    ''' 
    Get filenames of particle files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = ph.get_particle_files(subfolder_sat_byday_dummy,particles_folder,
                                        get_fields_files_sat_dummy_byday)
    return particle_files


@pytest.fixture()
def get_particle_files_sat_mf_byday(subfolder_sat_byday_mf,folder_names,
                                       get_fields_files_sat_mf_byday):
    ''' 
    Get filenames of particle files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = ph.get_particle_files(subfolder_sat_byday_mf,particles_folder,
                                        get_fields_files_sat_mf_byday)
    return particle_files

@pytest.fixture()
def get_particle_files_sat_bypoint(subfolder_sat_bypoint,folder_names,
                                         get_fields_files_sat_bypoint):
    ''' 
    Get filenames of particle files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = ph.get_particle_files(subfolder_sat_bypoint,particles_folder,
                                        get_fields_files_sat_bypoint)
    return particle_files

@pytest.fixture()
def define_heights():
    ''' Define height range for input into particle_locations function. '''
    dheights = 1000
    heights = np.arange(0, 19001, dheights) + dheights/2.
    return heights

@pytest.fixture()
def sat_byday_dummy_part():
    '''
    Defining dummy values for particle histograms. Should be consistent with
    the dummy particle location data which we can check against.
    
    Data include 4 id values, 2 points with 2 levels for each, 16 particles (2x2x4).
    Point 001 (id 1, 2): 2 particles at each of the south and west boundary (none at north and east)
    Point 002 (id 3, 4): 1 particle at each boundary (n, e, s, w)
    '''

    slice_dict = [{"time":[0],"lev":[0]},
                  {"time":[0],"lev":[1]},
                  {"time":[1],"lev":[0]},
                  {"time":[1],"lev":[1]}]

    part_hist = [{"pl_n":0,"pl_e":0,"pl_s":0.5,"pl_w":0.5},
                 {"pl_n":0,"pl_e":0,"pl_s":0.5,"pl_w":0.5},
                 {"pl_n":0.25,"pl_e":0.25,"pl_s":0.25,"pl_w":0.25},
                 {"pl_n":0.25,"pl_e":0.25,"pl_s":0.25,"pl_w":0.25}]

    return part_hist,slice_dict

def test_particle_locations_sat_dummy_byday(define_grid_sat_dummy_byday,
                                         get_particle_files_sat_dummy_byday,
                                         define_heights,sat_param_dummy,
                                         sat_byday_dummy_part):
    '''
    Test particle_locations() function can produce the correct histogram output for
    a known input when satellite points have been grouped by day.
    '''
    lons, lats, levs, time, timeStep = define_grid_sat_dummy_byday
    particle_file = get_particle_files_sat_dummy_byday[0]
    heights =  define_heights   
    upper_level = sat_param_dummy["upper_level"]
    
    out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
    expected_part_hist, slice_dict = sat_byday_dummy_part
    
    for id_slice,exp_hist_per_id in zip(slice_dict,expected_part_hist):
        for pl_key,value in exp_hist_per_id.items():
            pl_per_id = out[pl_key].isel(id_slice)
            assert pl_per_id.sum() == value

@pytest.mark.long
def test_particle_locations_sat_mf_byday(define_grid_sat_mf_byday,
                                         get_particle_files_sat_mf_byday,
                                         define_heights,sat_param):
    '''
    Test particle_locations() function can produce the correct output for satellite data when points have been
    grouped by day.
    '''
    lons, lats, levs, time, timeStep = define_grid_sat_mf_byday
    particle_file = get_particle_files_sat_mf_byday[0]
    heights =  define_heights   
    upper_level = sat_param["upper_level"]
    
    out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
    assert len(out["time"].values) > 1
    assert out
    ### TODO: ADD MORE STRINGENT TEST

@pytest.mark.long
def test_particle_locations_sat_bypoint(define_grid_sat_bypoint,get_particle_files_sat_bypoint,
                                            define_heights,sat_param):
    '''
    Test particle_locations() function can produce the correct output for satellite data when points are separate.
    '''
    lons, lats, levs, time, timeStep = define_grid_sat_bypoint
    particle_files = get_particle_files_sat_bypoint
    heights =  define_heights   
    upper_level = sat_param["upper_level"]
    
    for particle_file in particle_files:
        out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
        assert len(out["time"].values) == 1
        assert out
        ### TODO: ADD MORE STRINGENT TEST

#%% Test read_met()

@pytest.fixture()   
def get_met_files_sat_dummy_byday(subfolder_sat_byday_dummy,folder_names,get_fields_files_sat_dummy_byday):
    '''
    Get filenames of met files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    met_folder = folder_names["met_folder"]
    met_files = ph.get_met_files(subfolder_sat_byday_dummy,met_folder,get_fields_files_sat_dummy_byday,satellite=True)
    
    return met_files

@pytest.fixture()   
def get_met_files_sat_mf_byday(subfolder_sat_byday_mf,folder_names,get_fields_files_sat_mf_byday):
    '''
    Get filenames of met files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    met_folder = folder_names["met_folder"]
    met_files = ph.get_met_files(subfolder_sat_byday_mf,met_folder,get_fields_files_sat_mf_byday,satellite=True)
    
    return met_files

@pytest.fixture()
def get_met_files_sat_bypoint(subfolder_sat_bypoint,folder_names,get_fields_files_sat_bypoint):
    '''
    Get filenames of met files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    
    met_folder = folder_names["met_folder"]
    met_files = ph.get_met_files(subfolder_sat_bypoint,met_folder,get_fields_files_sat_bypoint,satellite=True)

    return met_files

@pytest.fixture()
def sat_byday_dummy_met():
    '''
    Defining dummy values for the satellite met data. Should be consistent with
    the dummy satellite met data which we can check against
    
    TODO: May want to read this in in some way rather than just define values
    '''
    values = {}
    
    values["release_lon"] = [-50.0,-50.0,-51.0,-51.0]
    values["release_lat"] = [-10.0,-10.0,-11.0,-11.0]
    values["temp"] = [25.0,20.0,25.0,20.0]
    values["PBLH"] = [600.0,590.0,730.0,730.0]
    values["press"] = [99000.0,95000.0,99863.93,94433.27]
    values["wind"] = [7.0,7.5,7.361012,7.119352]
    values["wind_direction"] = [100.0,110.0,43.85613,41.78423]
    values["label"] = ["00101","00102","00201","00202"]

    df = pd.DataFrame(values)

    return df

    #                    ,                        ,                        ,               Temperature (C),          Boundary layer depth,                 Pressure (Pa),                    Wind speed,      Wind direction (degrees),    
    #01/01/2013 14:58 UTC,               -50.00000,               -10.00000,                      25.00000,                      600.0000,                      99000.00,                      7.00000,                      100.0000,    
    #01/01/2013 14:58 UTC,               -50.00000,               -10.00000,                      20.00000,                      590.0000,                      95000.00,                      7.500000,                      110.0000,    
    #01/01/2013 15:01 UTC,               -51.00000,               -11.00000,                      25.00000,                      730.0000,                      99863.93,                      7.361012,                      43.85613,
    #01/01/2013 15:01 UTC,               -51.00000,               -11.00000,                      20.00000,                      730.0000,                      94433.27,                      7.119352,                      41.78423,

def test_read_met_sat_dummy_byday(get_met_files_sat_dummy_byday,sat_byday_dummy_met):
    '''
    Test output from read_met function against expected values.
    '''
    met_files_1 = get_met_files_sat_dummy_byday[0] # Met files for 1 day
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    expected = sat_byday_dummy_met
    compare_columns = expected.columns
    
    for col in compare_columns:
        assert np.array_equal(out[col],expected[col])
    
def test_read_met_sat_mf_byday(get_met_files_sat_mf_byday):
    '''
    Test read_met function can read in one met file when satellite data is grouped by day.
    '''
    met_files_1 = get_met_files_sat_mf_byday[0] # Met files for 1 day
    
    total_file_num = len(met_files_1)
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    assert len(out.index) == total_file_num

def test_read_met_sat_bypoint(get_met_files_sat_bypoint,sat_param):
    '''
    Test read_met function can read in one met file when satellite data is separated into points.
    '''
    met_files_1 = get_met_files_sat_bypoint[0]
    total_file_num = len(met_files_1)
   
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)

    assert len(out.index) == total_file_num

#%% Test footprint_array()

@pytest.fixture()
def read_met_sat_dummy_byday(get_met_files_sat_dummy_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Returns: list of dataframes.
    '''
    return ph.read_met_files(get_met_files_sat_dummy_byday,satellite=True)

@pytest.fixture()
def read_met_sat_mf_byday(get_met_files_sat_mf_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Returns: list of dataframes.
    '''
    return ph.read_met_files(get_met_files_sat_mf_byday,satellite=True)

@pytest.fixture()
def read_met_sat_bypoint(get_met_files_sat_bypoint):
    '''
    Read met data for satellite data is separated into points.
    Returns: list of dataframes.
    '''
    return ph.read_met_files(get_met_files_sat_bypoint,satellite=True)

@pytest.fixture()
def read_met_separate_sat_byday(get_met_files_sat_mf_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Separates points into individual dataframes.??
    Returns: nested list of dataframes.
    '''
    return ph.read_met_files_separate(get_met_files_sat_mf_byday,satellite=True)

@pytest.fixture()
def read_met_separate_sat_bypoint(get_met_files_sat_bypoint):
    '''
    Read met data for satellite run is separated into points.
    Separates points into individual dataframes.
    This is the input format expected by the original process script.
    Returns: nested list of dataframes.
    '''
    return ph.read_met_files_separate(get_met_files_sat_bypoint,satellite=True)

def test_met_satellite_split_dummy_byday(read_met_sat_dummy_byday,sat_byday_dummy_met):
    '''
    Test that output from met_satellite_split() function is as expected for 
    satellite data grouped by day.
    Using dummy dataset to check values and split is correct.
    '''
    # Expect output of this function to be a 2-item list containing dataframes
    # 1 dataframe per point - should have been grouped according to the level 
    met = process.met_satellite_split(read_met_sat_dummy_byday)
    
    # Need to break expected values out into each point for comparison as well
    compare = sat_byday_dummy_met
    compare_columns = compare.columns # Grab columns before adding extra point column.
    
    num_col = compare["label"].apply(lambda x: x[:-2])
    compare["point"] = num_col

    point_num = np.unique(num_col.values) # Get unique point values (should be 001,002)
    for num,data in zip(point_num,met):
        filt = compare["point"]==num
        expected = compare[filt] # Filter df of expected values to only contain one point
        for col in compare_columns:
            assert np.array_equal(data[col],expected[col])

def get_obs_prefix(subfolder,obs_folder):
    '''
    Create prefix for obs_folder.
    This is a string of the path based on the subfolder and obs_folder name.
    Example: /data/rt17603/NAME_output/SOUTHAMERICA_GOSAT-BRAZIL_column/Observations/
    Note: This is a new folder which has to be created by the user after the NAME run is completed.
    Args:
        subfolder (str):
            Path up to subfolder of NAME output. The final folder is based on the domain_site_height/ syntax
            e.g. SOUTHAMERICA_GOSAT-BRAZIL_column
        met_folder (str) :
            Folder containing the met files or folders. e.g. "Observations"
    Returns:
        str:
            Full path to observations files.
    '''
    return os.path.join(subfolder,obs_folder)

def get_obs_files(subfolder,obs_folder,field_files):
    '''
    Get filenames of all obs files matching date strings extracted from fields files.
    Only relevant to satellite data where we have to extract the times from the observations files rather 
    than work it out from the NAME input.
    Returns a list of filenames.
    '''
    
    obs_prefix = get_obs_prefix(subfolder,obs_folder)
    
    file_datestrs = ph.get_field_file_datestr(field_files) # Can be YYYYMMDD or YYYYMMDD-NNN
    
    obs_files = []
    for file_datestr in file_datestrs:
        obs_search_string = os.path.join(obs_prefix,"*{datestr}_*.nc".format(datestr=file_datestr))
        obs_file_search = glob.glob(obs_search_string)
        if obs_file_search:
            obs_files.extend(obs_file_search)

    return obs_files

@pytest.fixture()
def get_obs_files_sat_dummy_byday(subfolder_sat_byday_dummy,folder_names,
                                  get_fields_files_sat_dummy_byday):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_sat_byday_dummy,
                              obs_folder,get_fields_files_sat_dummy_byday)
    
    return obs_files

@pytest.fixture()
def get_obs_files_sat_mf_byday(subfolder_sat_byday_mf,folder_names,
                                  get_fields_files_sat_mf_byday):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_sat_byday_mf,
                              obs_folder,get_fields_files_sat_mf_byday)
    
    return obs_files

@pytest.fixture()
def get_obs_files_sat_bypoint(subfolder_sat_bypoint,folder_names,
                                  get_fields_files_sat_bypoint):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_sat_bypoint,
                              obs_folder,get_fields_files_sat_bypoint)
    
    return obs_files

@pytest.fixture()
def get_time_step_sat_dummy_byday(subfolder_sat_byday_dummy):
    ''' Extract time step value from within satellite_byday subfolder '''
    return ph.get_time_step(subfolder_sat_byday_dummy)

@pytest.fixture()
def get_time_step_sat_mf_byday(subfolder_sat_byday_mf):
    ''' Extract time step value from within satellite_byday subfolder '''
    return ph.get_time_step(subfolder_sat_byday_mf)

@pytest.fixture()
def get_time_step_sat_bypoint(subfolder_sat_bypoint):
    ''' Extract time step value from within satellite_byday subfolder '''
    return ph.get_time_step(subfolder_sat_bypoint)

def calc_area(df,lonlat_col=["X (Lat-Long)","Y (Lat-Long)"]):
    ''' Calculate the area using the dataframe extracted from a fields file '''
    lon = df[lonlat_col[0]].values
    lat = df[lonlat_col[1]].values
    
    areas = areagrid(lat,lon)
    areas_values = np.diagonal(areas)
    
    return areas_values

def convert_units_fp(value,area,timeStep):
    ''' Convert units from ppm s to (mol/mol)/(mol/m2/s)'''
    converted_value = value*area*1e-6*1./(3600.*timeStep*1.)
    return converted_value

@pytest.fixture()
def sat_byday_dummy_fields(get_fields_files_sat_dummy_byday,
                           get_time_step_sat_dummy_byday,
                           sat_param_dummy):
    '''
    Extract fields files values from dummy file and calculate unit conversion.
    Outputs a dataframe which has been stacked to include all values in one column
    and has the associated converted values.
    
    Also adds the correct indices for the associated fp dataset (time,lev,lat,lon)
    as:
        coord"_index" i.e. "lat_index","lon_index","time_index","lev_index"
    
    Returns:
        pandas.DataFrame,str:
            Dataframe containing the values, name of the new column containing 
            the converted values.
    '''
    fields_file_1 = get_fields_files_sat_dummy_byday[0]
    time_step = get_time_step_sat_dummy_byday
    
    header_rows = 36
    
    ## Extract values from fields files
    df = pd.read_csv(fields_file_1,skiprows=header_rows)
    
    prelim_header = 16
    fields_header = 17
    title_line2 = 21
    title_line1 = 36
    with gzip.open(fields_file_1) as temp:
        
        lines = temp.readlines()

        prelim_col = lines[prelim_header].split(b':')[-1].strip()
        num_prelim_col = int(prelim_col)
        
        fields_col = lines[fields_header].split(b':')[-1].strip()
        num_fields_col = int(fields_col)
        total_col = num_prelim_col+num_fields_col
        
        titles2 = lines[title_line2].split(b',')[num_prelim_col:total_col]
        titles2 = [l.strip().decode("utf-8") for l in titles2]
        
        titles1 = lines[title_line1].split(b',')[:num_prelim_col]
        titles1 = [l.strip().decode("utf-8") for l in titles1]

    # Create title line
    titles = titles1 + titles2
    
    df = df.dropna(how="all",axis=1)
    df.columns = titles

    # Calculate areas based on latitude and longitdue values in fields file
    areas = calc_area(df,lonlat_col=titles1[2:4])
    
    # Rearrange dataframe to put all values in one column
    var_name = "Run"
    value_name = "MixR"
    df_match = df.melt(id_vars=titles1,value_vars=titles2,
                         var_name=var_name,value_name=value_name)
    df_match = df_match.reset_index(drop=True)

    # Make sure there is a correct area value for each row in new dataframe
    areas_reshape = list(areas)*num_fields_col

    upper_level = sat_param_dummy["upper_level"]
    num_cells = len(df.index)

    ## Find correct indexes to relate to final fp dataset (time,lev,lat,lon)
    time_index = []
    lev_index = []
    conv_values = []
    for i in range(num_fields_col):
        
        if i == 0:
            time_id=0
            lev_id=0
        elif (i)%upper_level: # Loop through levels while i is not exactly divisable by upper_level
            lev_id+=1
        else: # Set time value to next value and reset level to 0
            time_id+=1
            lev_id=0
        
        for j in range(num_cells):
            
            time_index.append(time_id)
            lev_index.append(lev_id)
            
            # Calculate converted value (from ppm s to (mol/mol)/(mol/m2/s))
            k = i*num_cells+j
            a = areas_reshape[k]
            v = df_match[value_name].iloc[k]
            conv_value = convert_units_fp(v,a,time_step)
            conv_values.append(conv_value)

    # Add converted value to dataframe as a new column
    conv_name = "MixR_conv"
    df_match[conv_name] = np.array(conv_values)

    # Add index values so we can map to the fp dataset
    df_match["lon_index"] = df_match[titles1[0]] - 1
    df_match["lat_index"] = df_match[titles1[1]] - 1
    df_match["time_index"] = np.array(time_index)
    df_match["lev_index"] = np.array(lev_index)

    return df_match,conv_name

def test_footprint_array_sat_dummy_byday(get_fields_files_sat_dummy_byday,
                                         get_particle_files_sat_dummy_byday,
                                         read_met_sat_dummy_byday,
                                         get_obs_files_sat_dummy_byday,
                                         get_time_step_sat_dummy_byday,
                                         sat_param_dummy,
                                         sat_byday_dummy_fields):
    '''
    Test footprint_array() function against dummy satellite data.
    Test that conversion has been completed correctly and that values are in the
    correct positions within the derived lat/lon grid.
    '''
    fields_file_1 = get_fields_files_sat_dummy_byday[0]
    particle_file_1 = get_particle_files_sat_dummy_byday[0]
    met = read_met_sat_dummy_byday[0]
    #obs_file_1 = get_obs_files_sat_byday[0]
    time_step = get_time_step_sat_dummy_byday
    upper_level = sat_param_dummy["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,
                                  time_step=time_step,
                                  upper_level=upper_level)
    
    fp = out["fp"]
    
    df_expected,name = sat_byday_dummy_fields
    
    for index,row in df_expected.iterrows():
        slice_dict = {"time":[row["time_index"]],
                      "lat":row["lat_index"],
                      "lon":row["lon_index"],          
                      "lev":row["lev_index"]}
        value = fp[slice_dict].values[0]
        expected_value = row[name]
        assert value == expected_value   

@pytest.mark.long
def test_footprint_array_sat_mf_byday(get_fields_files_sat_mf_byday,
                                      get_particle_files_sat_mf_byday,
                                      read_met_sat_mf_byday,
                                      #get_obs_files_sat_mf_byday,
                                      get_time_step_sat_mf_byday,
                                      sat_param):
    '''
    Test footprint_array function can return an output for satellite files when grouped by day.
    '''
    fields_file_1 = get_fields_files_sat_mf_byday[0]
    particle_file_1 = get_particle_files_sat_mf_byday[0]
    met = read_met_sat_mf_byday[0]
    #obs_file_1 = get_obs_files_sat_byday[0]
    time_step = get_time_step_sat_mf_byday
    upper_level = sat_param["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,upper_level=upper_level)

@pytest.mark.long
def test_footprint_array_sat_bypoint(get_fields_files_sat_bypoint,
                                     get_particle_files_sat_bypoint,
                                     read_met_sat_bypoint,
                                     get_obs_files_sat_bypoint,
                                     get_time_step_sat_bypoint,
                                     sat_param):
    '''
    Test footprint_array function can return an output for satellite files when separated by point.
    '''
    point = 11
    
    fields_file_1 = get_fields_files_sat_bypoint[point-1]
    particle_file_1 = get_particle_files_sat_bypoint[point-1]
    met = read_met_sat_bypoint[point-1]
    #obs_file_1 = get_obs_files_sat_bypoint[0]
    time_step = get_time_step_sat_bypoint
    upper_level = sat_param["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,upper_level=upper_level,
                                  time_step=time_step)
    
    print(out["wind_direction"])

#%% Test footprint_concatenate()

@pytest.fixture()
def fc_param_sat_mf_byday(subfolder_sat_byday_mf,read_met_sat_mf_byday,
                             folder_names,sat_param,
                             get_time_step_sat_mf_byday):
                             #get_obs_files_sat_byday):
    '''
    Define parameters for input into footprint_concatenate function for satellite files when grouped by day.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = ph.footprint_concatenate_param(subfolder_sat_byday_mf,read_met_sat_mf_byday,
                                        folder_names,sat_param,satellite=True,
                                        fields="mixr",
                                        get_time_step=get_time_step_sat_mf_byday)
                                        #get_obs_files=get_obs_files_sat_byday)
    return param

@pytest.fixture()
def fc_param_sat_bypoint(subfolder_sat_bypoint,read_met_sat_bypoint,
                             folder_names,sat_param,
                             get_time_step_sat_bypoint):
                             #get_obs_files_sat_bypoint):
    '''
    Define parameters for input into footprint_concatenate function for satellite files when separated by point.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = ph.footprint_concatenate_param(subfolder_sat_bypoint,read_met_sat_bypoint,
                                           folder_names,sat_param,satellite=True,
                                           fields="conc",
                                           get_time_step=get_time_step_sat_bypoint)
                                           #get_obs_files=get_obs_files_sat_bypoint)
    
    return param

@pytest.mark.long
def test_footprint_concatenate_sat_mf_byday(fc_param_sat_mf_byday):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are grouped by day.
    '''
    param = fc_param_sat_mf_byday
    param["datestr"] = param["datestr"][0]
    
    out = process.footprint_concatenate(**param)
  
@pytest.mark.long
def test_footprint_concatenate_sat_bypoint(fc_param_sat_bypoint):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are separated by point.
    '''
    param = fc_param_sat_bypoint
    param["datestr"] = param["datestr"][0]
    param["met"] = param["met"][0]

    out = process.footprint_concatenate(**param)

#%% Test satellite_vertical_profile()
    
@pytest.fixture()
def footprint_concatenate_sat_bypoint(fc_param_sat_bypoint):
    '''
    Create output from footprint_concatenate function for satellite data seperated by point.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per point.
    '''
    param = fc_param_sat_bypoint.copy()
    met_points = fc_param_sat_bypoint["met"][:]
    date_strings = fc_param_sat_bypoint["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process.footprint_concatenate(**param))
    return fp_all

@pytest.fixture()
def footprint_concatenate_sat_mf_byday(fc_param_sat_mf_byday):
    '''
    Create output from footprint_concatenate function for satellite data grouped by day.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per day.
    '''
    param = fc_param_sat_mf_byday.copy()
    met_points = fc_param_sat_mf_byday["met"][:]
    date_strings = fc_param_sat_mf_byday["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process.footprint_concatenate(**param))
    
    return fp_all

@pytest.mark.long
def test_satellite_vertical_profile_bypoint(footprint_concatenate_sat_bypoint,
                                          get_obs_files_sat_bypoint,sat_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is separated by point.
    '''
    sat_obs_file = get_obs_files_sat_bypoint[0]
    fp = footprint_concatenate_sat_bypoint[0]
    out = process.satellite_vertical_profile(fp,sat_obs_file,max_level=sat_param["max_level"])

    assert out is not None

@pytest.mark.long
def test_satellite_vertical_profile_mf_byday(footprint_concatenate_sat_mf_byday,
                                          get_obs_files_sat_mf_byday,sat_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is grouped by day.
    '''
    sat_obs_file = get_obs_files_sat_mf_byday[0]
    fp = footprint_concatenate_sat_mf_byday[0]
    out = process.satellite_vertical_profile(fp,sat_obs_file,max_level=sat_param["max_level"])

    assert out is not None

#%% Test process()

@pytest.fixture()
def process_sat_byday_mf_param(sat_param,folder_names,sat_byday_mf_directory):
    '''
    Define input parameters for process.process function for satellite data grouped by day.
    Additional "base_dir" parameter added to point to Satellite_ByDay folder.
    '''
    param = ph.process_param(sat_param,folder_names,satellite=True, fields="mixr")
    param["base_dir"] = sat_byday_mf_directory
    param["process_dir"] = os.path.join(sat_byday_mf_directory,"Processed_Fields_files")

    return param

@pytest.fixture()
def process_sat_bypoint_param(sat_param,folder_names,sat_bypoint_directory):
    '''
    Define input parameters for process.process function for satellite data separated by point.
    Additional "base_dir" parameter added to point to Satellite_ByPoint folder.
    '''    
    param = ph.process_param(sat_param,folder_names,satellite=True, fields="conc")
    param["base_dir"] = sat_bypoint_directory
    param["process_dir"] = os.path.join(sat_bypoint_directory,"Processed_Fields_files")
    
    return param

@pytest.mark.long
def test_process_sat_bypoint(process_sat_bypoint_param):
    '''
    Test process function produces an output for satellite data separated by point.
    '''
    
    ph.remove_processed_file(process_sat_bypoint_param)    
    out = process.process(**process_sat_bypoint_param)

    assert out is not None

@pytest.mark.long
def test_process_sat_mf_byday(process_sat_byday_mf_param):
    '''
    Test process function produces an output for satellite data grouped by day in mole fraction output.
    '''
    
    ph.remove_processed_file(process_sat_byday_mf_param)
    out = process.process(**process_sat_byday_mf_param)

    assert out is not None

#@pytest.mark.skip(reason="Possible comparison error due to new short-lived footprints. Requires updating.")
def test_process_sat_byday_mf_benchmark(process_sat_byday_mf_param,
                                        benchmark_output_directory): 
    '''
    Test mol fraction output against benchmarked file (created 2019-05-20). This marks
    a code update to check the units within the fields file and then use the correct unit 
    conversion. NAME output should now primarily be in ppm s units rather than gs/m3 but 
    this should ensure backwards compatibility with previous runs.
    '''
    
    ph.remove_processed_file(process_sat_byday_mf_param)    
    process.process(**process_sat_byday_mf_param)

    process.process(**process_sat_byday_mf_param)
    output_file = ph.find_processed_file(process_sat_byday_mf_param)
    print(output_file)
    out = open_ds(output_file)

    # Benchmark was run on 15/09/2021
    # GOSAT-BRAZIL, SOUTHAMERICA domain, January, 2013
    # *If this does not match to your current inputs - benchmark must be updated.*

    benchmark_file = os.path.join(benchmark_output_directory, "GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc")
    benchmark = xr.open_dataset(benchmark_file)

    xr.testing.assert_equal(out, benchmark)
