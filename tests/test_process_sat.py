#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 10:15:58 2020

@author: rt17603
"""

import pytest
import os
import sys
import glob
import gzip
import numpy as np
import pandas as pd
import xarray as xray
import acrg_name.process as process
import deprecated.process as process_org
from acrg_name.name import open_ds
import process_general as pg
import pdb

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg

#%% Setting up directory fixtures

@pytest.fixture()
def sat_byday_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay/")
    return directory

@pytest.fixture()
def sat_bypoint_directory():
    ''' Define base directory containing satellite footprint files with separated data points '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByPoint/")
    return directory

@pytest.fixture()
def sat_byday_conc_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay_MF/")
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
    #param["max_level"] = 17
    
    return param

@pytest.fixture()
def sat_param_dummy(sat_param):
    '''
    Define input parameters for dummy satellite data (based on sat_param)
    '''
    param = sat_param
    sat_param["max_level"] = 2
    return param

@pytest.fixture()
def subfolder_sat(sat_param):
    ''' Define subfolder name for satellite data'''
    domain = sat_param["domain"]
    site = sat_param["site"]
    height = sat_param["height"]
    
    subfolder = pg.define_subfolder(domain,site,height)
    return subfolder

@pytest.fixture()
def subfolder_sat_byday(sat_byday_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data with points grouped by day'''
    subfolder = os.path.join(sat_byday_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def subfolder_sat_bypoint(sat_bypoint_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data with separate points'''
    subfolder = os.path.join(sat_bypoint_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def subfolder_sat_byday_mf(sat_byday_mf_directory,subfolder_sat):
    ''' Define full subfolder path for satellite data with points grouped by day'''
    subfolder = os.path.join(sat_byday_mf_directory,subfolder_sat)  
    return subfolder

@pytest.fixture()
def folder_names():
    '''
    Define subfolders within NAME output directory structure
    Keys includes:
        "fields_folder","particles_folder","met_folder","processed_folder","obs_folder"
    '''
    param = {}
    param["fields_folder"] = "Fields_files"
    param["mixr_folder"] = "MixR_files"
    param["particles_folder"] = "Particle_files"
    param["met_folder"] = "Met"
    param["processed_folder"] = "Processed_Fields_files"
    param["obs_folder"] = "Observations"

    return param  

#%% Testing read_fields() function

@pytest.fixture()
def get_fields_files_sat_byday(subfolder_sat_byday,folder_names,sat_param_dummy):
    '''
    Get filenames of fields files for satellite run with points grouped by day.
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = pg.create_datestr(sat_param_dummy)
    #fields_folder = folder_names["fields_folder"]
    mixr_folder = folder_names["mixr_folder"]
    
    fields_files = pg.get_fields_files(subfolder_sat_byday,mixr_folder,datestr)
    
    return fields_files

@pytest.fixture()
def get_fields_files_sat_bypoint(subfolder_sat_bypoint,folder_names,sat_param):
    '''
    Get filenames of fields files for satellite run with separate points.
    Finds field files based on datestr (see create_datestr), multiple files per day
    (one for each point, filename syntax *YYYYMMDD-NNN*).
    '''
    datestr = pg.create_datestr(sat_param)
    fields_folder = folder_names["fields_folder"]
    
    fields_files = pg.get_fields_files(subfolder_sat_bypoint,fields_folder,datestr)
    
    return fields_files

#@pytest.fixture()
#def get_fields_files_sat_byday_mf(subfolder_sat_byday_mf,folder_names,sat_param):
#    '''
#    Get filenames of fields files for satellite run with points grouped by day.
#    Finds field files based on datestr (see create_datestr), one file per day.
#    '''
#    datestr = pg.create_datestr(sat_param)
#    fields_folder = folder_names["fields_folder"]
#    
#    fields_files = pg.get_fields_files(subfolder_sat_byday_mf,fields_folder,datestr)
#    
#    return fields_files

@pytest.fixture()
def read_fields_file_sat_byday(get_fields_files_sat_byday):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_byday[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    #print "header,column_headings",header,column_headings
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_fields_file_sat_bypoint(get_fields_files_sat_bypoint):
    ''' Read first satellite by point field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_bypoint[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_fields_file_sat_byday_mf(get_fields_files_sat_byday_mf):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_sat_byday_mf[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

def test_read_field_file_sat_bypoint(get_fields_files_sat_bypoint):
    ''' Test read_file function can run for satellite data separated by point. '''
    point = 10
    fields_file = get_fields_files_sat_bypoint[point-1]
    header, column_headings, data_arrays = process.read_file(fields_file)

def test_read_field_file_sat_byday(get_fields_files_sat_byday):
    '''  Test read_file function can run for satellite data separated by day. '''
    fields_file = get_fields_files_sat_byday[0]
    header, column_headings, data_arrays = process.read_file(fields_file)

#    point = 10
#    da_range = [(point-1)*17,(point-1)*17+17]
#    data_arrays_point = data_arrays[da_range[0]:da_range[1]]
#    for i,array in enumerate(data_arrays_point):
#        print i,np.min(array),np.max(array),np.mean(array),np.sum(array)
#    print "Over 17 elements of point {0}".format(point),np.mean(data_arrays_point)


    ### TODO: ADD MORE STRINGENT TEST    

def test_define_grid_sat_byday(read_fields_file_sat_byday,sat_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with points 
    grouped by day.
    '''
    header,column_headings,data_arrays = read_fields_file_sat_byday
    out = process.define_grid(header,column_headings,satellite=True,upper_level=sat_param["upper_level"])
    
    assert out

    ### TODO: ADD MORE STRINGENT TEST

def test_define_grid_sat_bypoint(read_fields_file_sat_bypoint,sat_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with separate 
    points.
    '''
    header,column_headings,data_arrays = read_fields_file_sat_bypoint
    out = process.define_grid(header,column_headings,satellite=True,upper_level=sat_param["upper_level"])

    assert out
    ### TODO: ADD MORE STRINGENT TEST

#%% Test particle_locations()

@pytest.fixture()
def define_grid_sat_byday(read_fields_file_sat_byday,sat_param):
    ''' Create output from process.define_grid() function for one satellite by day field file '''
    header,column_headings,data_arrays = read_fields_file_sat_byday
    upper_level=sat_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def define_grid_sat_bypoint(read_fields_file_sat_bypoint,sat_param):
    ''' Create output from process.define_grid() function for one satellite by point field file '''
    header,column_headings,data_arrays = read_fields_file_sat_bypoint
    upper_level=sat_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def get_particle_files_sat_byday(subfolder_sat_byday,folder_names,
                                       get_fields_files_sat_byday):
    ''' 
    Get filenames of particle files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = pg.get_particle_files(subfolder_sat_byday,particles_folder,
                                        get_fields_files_sat_byday)
    return particle_files

@pytest.fixture()
def get_particle_files_sat_bypoint(subfolder_sat_bypoint,folder_names,
                                         get_fields_files_sat_bypoint):
    ''' 
    Get filenames of particle files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = pg.get_particle_files(subfolder_sat_bypoint,particles_folder,
                                        get_fields_files_sat_bypoint)
    return particle_files

def test_particle_locations_sat_byday(define_grid_sat_byday,get_particle_files_sat_byday,
                                            define_heights,sat_param):
    '''
    Test particle_locations() function can produce the correct output for satellite data when points have been
    grouped by day.
    '''
    lons, lats, levs, time, timeStep = define_grid_sat_byday
    particle_file = get_particle_files_sat_byday[0]
    heights =  define_heights   
    upper_level = sat_param["upper_level"]
    
    out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
    assert len(out["time"].values) > 1
    assert out
    ### TODO: ADD MORE STRINGENT TEST
   
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

@pytest.fixture()
def define_grid_org_sat_bypoint(read_fields_file_sat_bypoint):
    ''' 
    Create output from original process script: process_org.define_grid() function for one satellite by point 
    field file.
    For comparison with new process script.
    '''
    header,column_headings,data_arrays = read_fields_file_sat_bypoint
    lons, lats, levs, time, timeStep = process_org.define_grid(header,column_headings,satellite=True)
    
    return lons, lats, levs, time, timeStep

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_particle_locations_sat_bypoint_against_org(define_grid_sat_bypoint,
                                                         define_grid_org_sat_bypoint,
                                                         get_particle_files_sat_bypoint,
                                                         define_heights,sat_param):
    '''
    Test particle locations output is identical between particle files processed using the original
    process script and the new process script. 
    This is for satellite output separated by point since original script cannot process by day output.
    '''
    lons1, lats1, levs1, time1, timeStep1 = define_grid_sat_bypoint
    lons2, lats2, levs2, time2, timeStep2 = define_grid_org_sat_bypoint
    particle_files = get_particle_files_sat_bypoint
    heights =  define_heights   
    upper_level = sat_param["upper_level"]
    
    pl_data_vars = ["pl_n","pl_e","pl_s","pl_w"]
    
    for i,particle_file in enumerate(particle_files):
        out = process.particle_locations(particle_file,time1,lats1,lons1,levs1,heights,id_is_lev=False,satellite=True,
                                         upper_level=upper_level)
        out_org = process_org.particle_locations(particle_file,time2,lats2,lons2,levs2,heights,id_is_lev=True)
    
        for pl in pl_data_vars:
            assert np.array_equal(out[pl].values,out_org[pl].values)

#%% Test read_met()

@pytest.fixture()   
def get_met_files_sat_byday(subfolder_sat_byday,folder_names,get_fields_files_sat_byday):
    '''
    Get filenames of met files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    met_folder = folder_names["met_folder"]
    met_files = pg.get_met_files(subfolder_sat_byday,met_folder,get_fields_files_sat_byday,satellite=True)
    
    return met_files

@pytest.fixture()
def get_met_files_sat_bypoint(subfolder_sat_bypoint,folder_names,get_fields_files_sat_bypoint):
    '''
    Get filenames of met files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    
    met_folder = folder_names["met_folder"]
    met_files = pg.get_met_files(subfolder_sat_bypoint,met_folder,get_fields_files_sat_bypoint,satellite=True)

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
    
def test_read_met_sat_byday(get_met_files_sat_byday):
    '''
    Test read_met function can read in one met file when satellite data is grouped by day.
    '''
    met_files_1 = get_met_files_sat_byday[0] # Met files for 1 day
    
    total_file_num = len(met_files_1)
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    assert len(out.index) == total_file_num

def test_read_met_sat_byday_values(get_met_files_sat_byday,sat_byday_dummy_met):
    '''
    Test output from read_met function against expected values.
    '''
    met_files_1 = get_met_files_sat_byday[0] # Met files for 1 day
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    expected = sat_byday_dummy_met
    compare_columns = expected.columns
    
    for col in compare_columns:
        assert np.array_equal(out[col],expected[col])

def test_read_met_sat_bypoint(get_met_files_sat_bypoint,sat_param):
    '''
    Test read_met function can read in one met file when satellite data is separated into points.
    '''
    met_files_1 = get_met_files_sat_bypoint[0]
    total_file_num = len(met_files_1)
   
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)

    assert len(out.index) == total_file_num

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_read_met_sat_bypoint_against_org(get_met_files_sat_bypoint,sat_param):
    '''
    Test values within output from original process script: process_org.read_met() function is identical to
    met data when processed using new process script.
    '''
    met_files_1 = get_met_files_sat_bypoint[0]
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    out_org_list = []
    for met_file in met_files_1:
        out_org_list.append(process_org.read_met(met_file))
    
    if "label" in out.columns:
        out.drop("label",inplace=True,axis=1)
    
    for t,out_org in enumerate(out_org_list):
        assert np.all(out.iloc[t] == out_org)

#%% Test footprint_array()

@pytest.fixture()
def read_met_sat_byday(get_met_files_sat_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Returns: list of dataframes.
    '''
    return pg.read_met_files(get_met_files_sat_byday,satellite=True)

@pytest.fixture()
def read_met_sat_bypoint(get_met_files_sat_bypoint):
    '''
    Read met data for satellite data is separated into points.
    Returns: list of dataframes.
    '''
    return pg.read_met_files(get_met_files_sat_bypoint,satellite=True)

@pytest.fixture()
def read_met_separate_sat_byday(get_met_files_sat_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Separates points into individual dataframes.??
    Returns: nested list of dataframes.
    '''
    return pg.read_met_files_separate(get_met_files_sat_byday,satellite=True)

@pytest.fixture()
def read_met_separate_sat_bypoint(get_met_files_sat_bypoint):
    '''
    Read met data for satellite run is separated into points.
    Separates points into individual dataframes.
    This is the input format expected by the original process script.
    Returns: nested list of dataframes.
    '''
    return pg.read_met_files_separate(get_met_files_sat_bypoint,satellite=True)

def test_met_satellite_split_byday(read_met_sat_byday,sat_byday_dummy_met):
    '''
    Test that output from met_satellite_split() function is as expected for 
    satellite data grouped by day.
    Using dummy dataset to check values and split is correct.
    '''
    # Expect output of this function to be a 2-item list containing dataframes
    # 1 dataframe per point - should have been grouped according to the level 
    met = process.met_satellite_split(read_met_sat_byday)
    
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
    
    file_datestrs = pg.get_field_file_datestr(field_files) # Can be YYYYMMDD or YYYYMMDD-NNN
    
    obs_files = []
    for file_datestr in file_datestrs:
        obs_search_string = os.path.join(obs_prefix,"*{datestr}_*.nc".format(datestr=file_datestr))
        obs_file_search = glob.glob(obs_search_string)
        if obs_file_search:
            obs_files.extend(obs_file_search)

    return obs_files

@pytest.fixture()
def get_obs_files_sat_byday(subfolder_sat_byday,folder_names,
                                  get_fields_files_sat_byday):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_sat_byday,
                              obs_folder,get_fields_files_sat_byday)
    
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
def get_time_step_sat_byday(subfolder_sat_byday):
    ''' Extract time step value from within satellite_byday subfolder '''
    return pg.get_time_step(subfolder_sat_byday)

@pytest.fixture()
def get_time_step_sat_bypoint(subfolder_sat_bypoint):
    ''' Extract time step value from within satellite_byday subfolder '''
    return pg.get_time_step(subfolder_sat_bypoint)

@pytest.fixture()
def sat_byday_dummy_fields(get_fields_files_sat_byday):
    '''
    '''

    fields_file_1 = get_fields_files_sat_byday[0]
    
    header_rows = 37
    
    df = pd.read_csv(fields_file_1,skiprows=header_rows)
    
    prelim_header = 17
    title_line2 = 21
    title_line1 = 37
    with gzip.open(fields_file_1) as temp:
        
        lines = temp.readlines()

        prelim_col = lines[prelim_header].split(b':')[-1].strip()
        num_prelim_col = int(prelim_col)
        
        titles2 = lines[title_line2].split(b',')[num_prelim_col:]
        titles2 = [l.strip().decode("utf-8") for l in titles2]
        
        titles1 = lines[title_line1].split(b',')[:num_prelim_col]
        titles1 = [l.strip().decode("utf-8") for l in titles1]

    titles = titles1 + titles2
    
    df.columns = titles
    df = df.dropna(how="all",axis=1)

    return df,titles1,titles2

    # values = {}
    
    # values["release_lon"] = [-50.0,-50.0,-51.0,-51.0]
    # values["release_lat"] = [-10.0,-10.0,-11.0,-11.0]
    # values["temp"] = [25.0,20.0,25.0,20.0]
    # values["PBLH"] = [600.0,590.0,730.0,730.0]
    # values["press"] = [99000.0,95000.0,99863.93,94433.27]
    # values["wind"] = [7.0,7.5,7.361012,7.119352]
    # values["wind_direction"] = [100.0,110.0,43.85613,41.78423]
    # values["label"] = ["00101","00102","00201","00202"]

    # df = pd.DataFrame(values)

    # return df

    # #                    ,                        ,                        ,               Temperature (C),          Boundary layer depth,                 Pressure (Pa),                    Wind speed,      Wind direction (degrees),    
    # #01/01/2013 14:58 UTC,               -50.00000,               -10.00000,                      25.00000,                      600.0000,                      99000.00,                      7.00000,                      100.0000,    
    # #01/01/2013 14:58 UTC,               -50.00000,               -10.00000,                      20.00000,                      590.0000,                      95000.00,                      7.500000,                      110.0000,    
    # #01/01/2013 15:01 UTC,               -51.00000,               -11.00000,                      25.00000,                      730.0000,                      99863.93,                      7.361012,                      43.85613,
    # #01/01/2013 15:01 UTC,               -51.00000,               -11.00000,   

def test_footprint_array_sat_byday(get_fields_files_sat_byday,
                                         get_particle_files_sat_byday,
                                         read_met_sat_byday,
                                         get_obs_files_sat_byday,
                                         get_time_step_sat_byday,
                                         sat_param_dummy):
    '''
    Test footprint_array function can return an output for satellite files when grouped by day.
    '''
    fields_file_1 = get_fields_files_sat_byday[0]
    particle_file_1 = get_particle_files_sat_byday[0]
    met = read_met_sat_byday[0]
    #obs_file_1 = get_obs_files_sat_byday[0]
    time_step = get_time_step_sat_byday
    upper_level = sat_param["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,upper_level=upper_level)
    

#TODO: ADD THIS TEST OF VALUES HERE!
def test_footprint_array_sat_byday_values(get_fields_files_sat_byday,
                                        get_particle_files_sat_byday,
                                        read_met_sat_byday,
                                        get_obs_files_sat_byday,
                                        get_time_step_sat_byday,
                                        sat_param_dummy,
                                        sat_byday_dummy_fields):
    '''
    
    '''
    fields_file_1 = get_fields_files_sat_byday[0]
    particle_file_1 = get_particle_files_sat_byday[0]
    met = read_met_sat_byday[0]
    #obs_file_1 = get_obs_files_sat_byday[0]
    time_step = get_time_step_sat_byday
    upper_level = sat_param_dummy["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,
                                  upper_level=upper_level)
    
    expected,axis_cols,field_cols = sat_byday_dummy_fields
    
    pdb.set_trace()

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

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_footprint_array_sat_bypoint_against_org(get_fields_files_sat_bypoint,
                                           get_particle_files_sat_bypoint,
                                           read_met_sat_bypoint,
                                           read_met_separate_sat_bypoint,
                                           get_obs_files_sat_bypoint,
                                           get_time_step_sat_bypoint,
                                           sat_param):
    '''
    Test output from original process script: process_org.footprint_array() produces an identical output
    to the new process script over all files extracted based on field files datestr values.
    '''
    fields_files_1 = get_fields_files_sat_bypoint
    particle_files_1 = get_particle_files_sat_bypoint
    met_1 = read_met_sat_bypoint
    met_2 = read_met_separate_sat_bypoint
    time_step = get_time_step_sat_bypoint
    upper_level = sat_param["upper_level"]
    
    for fields_file,particle_file,met_data_1,met_data_2 in zip(fields_files_1,particle_files_1,met_1,met_2):
        
        out = process.footprint_array(fields_file,particle_file,met_data_1,satellite=True,time_step=time_step,upper_level=upper_level)
        out_org = process_org.footprint_array(fields_file,particle_file,met_data_2,satellite=True,time_step=time_step)
        
        data_vars = out.data_vars
        
        for dv in data_vars:
            assert np.array_equal(out[dv].values,out_org[dv].values)

#%% Test footprint_concatenate()

@pytest.fixture()
def fc_param_sat_byday(subfolder_sat_byday,read_met_sat_byday,
                             folder_names,sat_param,
                             get_time_step_sat_byday):
                             #get_obs_files_sat_byday):
    '''
    Define parameters for input into footprint_concatenate function for satellite files when grouped by day.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = pg.footprint_concatenate_param(subfolder_sat_byday,read_met_sat_byday,
                                        folder_names,sat_param,satellite=True,
                                        get_time_step=get_time_step_sat_byday)
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
    param = pg.footprint_concatenate_param(subfolder_sat_bypoint,read_met_sat_bypoint,
                                        folder_names,sat_param,satellite=True,
                                        get_time_step=get_time_step_sat_bypoint)
                                        #get_obs_files=get_obs_files_sat_bypoint)
    
    return param

@pytest.fixture()
def fc_param_org_sat_bypoint(fc_param_sat_bypoint,read_met_separate_sat_bypoint):
    '''
    Define parameters for input into original process script: process_org.footprint_concatenate() for
    satellite files when separated by point.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = fc_param_sat_bypoint.copy()
    if "upper_level" in list(param.keys()):
        param.pop("upper_level")
    param["met"] = read_met_separate_sat_bypoint # Re-define based on orginal script expected input.
    
    return param
  
def test_footprint_concatenate_sat_byday(fc_param_sat_byday):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are grouped by day.
    '''
    param = fc_param_sat_byday
    param["datestr"] = param["datestr"][0]
    
    out = process.footprint_concatenate(**param)
  
def test_footprint_concatenate_sat_bypoint(fc_param_sat_bypoint):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are separated by point.
    '''
    param = fc_param_sat_bypoint
    param["datestr"] = param["datestr"][0]
    param["met"] = param["met"][0]

    out = process.footprint_concatenate(**param)

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_footprint_concatenate_sat_bypoint_against_org(fc_param_sat_bypoint,
                                                             fc_param_org_sat_bypoint):
    '''
    Test output from original process script: process_org.footprint_concatenate against output from new
    process script for a satellite run where output is separated by point.
    '''
    
    param = fc_param_sat_bypoint.copy() 
    met_points = fc_param_sat_bypoint["met"][:]
    date_strings = fc_param_sat_bypoint["datestr"][:]
    
    param_org = fc_param_org_sat_bypoint.copy()
    met_points_org = fc_param_org_sat_bypoint["met"][:]
    date_strings_org = fc_param_sat_bypoint["datestr"][:]

    out_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        out_all.append(process.footprint_concatenate(**param))
    
    out_org_all = []
    for met,datestr in zip(met_points_org,date_strings_org):
        print("inputs",param_org)
        param_org["met"] = met
        param_org["datestr"] = datestr
        out_org_all.append(process_org.footprint_concatenate(**param_org))
    
    data_vars = out_all[0].data_vars
    
    for out, out_org in zip(out_all,out_org_all):
    
        for dv in data_vars:
            assert np.array_equal(out[dv].values,out_org[dv].values)

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
def footprint_concatenate_sat_byday(fc_param_sat_byday):
    '''
    Create output from footprint_concatenate function for satellite data grouped by day.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per day.
    '''
    param = fc_param_sat_byday.copy()
    met_points = fc_param_sat_byday["met"][:]
    date_strings = fc_param_sat_byday["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process.footprint_concatenate(**param))
    return fp_all

def test_satellite_vertical_profile_bypoint(footprint_concatenate_sat_bypoint,
                                          get_obs_files_sat_bypoint,sat_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is separated by point.
    '''
    sat_obs_file = get_obs_files_sat_bypoint[0]
    fp = footprint_concatenate_sat_bypoint[0]
    out = process.satellite_vertical_profile(fp,sat_obs_file,max_level=sat_param["max_level"])

def test_satellite_vertical_profile_byday(footprint_concatenate_sat_byday,
                                          get_obs_files_sat_byday,sat_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is grouped by day.
    '''
    sat_obs_file = get_obs_files_sat_byday[0]
    fp = footprint_concatenate_sat_byday[0]
    out = process.satellite_vertical_profile(fp,sat_obs_file,max_level=sat_param["max_level"])

@pytest.fixture()
def footprint_concatenate_org_sat_bypoint(fc_param_org_sat_bypoint):
    '''
    Create output from original process script: process_org.footprint_concatenate function for satellite data 
    seperated by point.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per point.
    '''
    param = fc_param_org_sat_bypoint.copy()
    met_points = fc_param_org_sat_bypoint["met"][:]
    date_strings = fc_param_org_sat_bypoint["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process_org.footprint_concatenate(**param))
    
    return fp_all

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_satellite_vertical_profile_bypoint_against_org(footprint_concatenate_sat_bypoint,
                                                        footprint_concatenate_org_sat_bypoint,
                                                        get_obs_files_sat_bypoint,sat_param):
    '''
    Test output from original process script: process_org.satellite_vertical_profile matches the output 
    from new process script for satellite data separated into points.
    '''
   
    sat_obs_files = get_obs_files_sat_bypoint
    fp_all = footprint_concatenate_sat_bypoint
    fp_all_org = footprint_concatenate_org_sat_bypoint
    
    for fp,fp_org,sat_obs_file in zip(fp_all,fp_all_org,sat_obs_files):
        out = process.satellite_vertical_profile(fp,sat_obs_file,max_level=sat_param["max_level"])
        out_org = process_org.satellite_vertical_profile(fp_org,sat_obs_file,
                                                                  max_level=sat_param["max_level"])
        data_vars = out.data_vars
        for dv in data_vars:
            assert np.array_equal(out[dv].values,out_org[dv].values)

#%% Test process()

@pytest.fixture()
def process_sat_byday_param(sat_param,folder_names,sat_byday_directory):
    '''
    Define input parameters for process.process function for satellite data grouped by day.
    Additional "base_dir" parameter added to point to Satellite_ByDay folder.
    '''
    param = pg.process_param(sat_param,folder_names,satellite=True)
    param["base_dir"] = sat_byday_directory
    param["process_dir"] = os.path.join(sat_byday_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_sat_bypoint_param(sat_param,folder_names,sat_bypoint_directory):
    '''
    Define input parameters for process.process function for satellite data separated by point.
    Additional "base_dir" parameter added to point to Satellite_ByPoint folder.
    '''    
    param = pg.process_param(sat_param,folder_names,satellite=True)
    param["base_dir"] = sat_bypoint_directory
    param["process_dir"] = os.path.join(sat_bypoint_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_sat_org_bypoint_param(sat_param,folder_names,sat_bypoint_directory):
    '''
    Define input parameters for process.process function for satellite data separated by point.
    Additional "base_dir" parameter added to point to Satellite_ByPoint folder.
    '''    
    param = pg.process_param(sat_param,folder_names,satellite=True)
    param["base_dir"] = sat_bypoint_directory
    if "upper_level" in list(param.keys()):
        param.pop("upper_level")
    param["process_dir"] = os.path.join(sat_bypoint_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_sat_byday_mf_param(sat_param,folder_names,sat_byday_mf_directory):
    '''
    Define input parameters for process.process function for satellite data grouped by day.
    Additional "base_dir" parameter added to point to Satellite_ByDay folder.
    '''
    param = pg.process_param(sat_param,folder_names,satellite=True)
    param["base_dir"] = sat_byday_mf_directory
    param["process_dir"] = os.path.join(sat_byday_mf_directory,"Processed_Fields_files")
    
    return param

       
def test_process_sat_byday(process_sat_byday_param,
                                 subfolder_sat_byday,folder_names,sat_param):
    '''
    Test process function produces an output for satellite data grouped by day.
    '''
    
    processed_folder = folder_names["processed_folder"]
    pg.remove_processed_file(subfolder_sat_byday,processed_folder,sat_param)
    
    out = process.process(**process_sat_byday_param)

def test_process_sat_bypoint(process_sat_bypoint_param,
                                   subfolder_sat_bypoint,folder_names,sat_param):
    '''
    Test process function produces an output for satellite data separated by point.
    '''
    
    processed_folder = folder_names["processed_folder"]
    pg.remove_processed_file(subfolder_sat_bypoint,processed_folder,sat_param)
    
    out = process.process(**process_sat_bypoint_param)

def test_process_sat_mf_byday(process_sat_byday_mf_param,
                                 subfolder_sat_byday_mf,folder_names,sat_param):
    '''
    Test process function produces an output for satellite data grouped by day in mole fraction output.
    '''
    
    processed_folder = folder_names["processed_folder"]
    pg.remove_processed_file(subfolder_sat_byday_mf,processed_folder,sat_param)
    
    out = process.process(**process_sat_byday_mf_param)

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_process_sat_bypoint_against_org(process_sat_bypoint_param,
                                               process_sat_org_bypoint_param,
                                               subfolder_sat_bypoint,folder_names,sat_param): 
    '''
    Test output of original process function: process_org.process against new process script for satellite
    data separated by point.
    '''
    processed_folder = folder_names["processed_folder"]
    pg.remove_processed_file(subfolder_sat_bypoint,processed_folder,sat_param)

    out = process.process(**process_sat_bypoint_param)

    processed_folder_org = "Processed_Fields_files_Org"    
    org_folder = os.path.join(subfolder_sat_bypoint,processed_folder_org)
    if not os.path.exists(org_folder):
        os.makedirs(org_folder)
    
    pg.remove_processed_file(subfolder_sat_bypoint,processed_folder_org,sat_param)
    
    process_sat_org_bypoint_param["processed_folder"] = processed_folder_org
    out_org = process_org.process(**process_sat_org_bypoint_param)
    
    data_vars = out.data_vars
    
    for dv in data_vars:
        assert np.array_equal(out[dv].values,out_org[dv].values)   

@pytest.mark.skip(reason="Possible comparison error due to new short-lived footprints. Requires updating.")
def test_process_sat_byday_mf_against_bench(process_sat_byday_mf_param,
                                                  subfolder_sat_byday_mf,
                                                  folder_names,sat_param): 
    '''
    Test mol fraction output against benchmarked file (created 2019-05-20). This marks
    a code update to check the units within the fields file and then use the correct unit 
    conversion. NAME output should now primarily be in ppm s units rather than gs/m3 but 
    this should ensure backwards compatibility with previous runs.
    '''
    processed_folder = folder_names["processed_folder"]
    pg.remove_processed_file(subfolder_sat_byday_mf,processed_folder,sat_param)

    process.process(**process_sat_byday_mf_param)
    output_file = pg.find_processed_file(subfolder_sat_byday_mf,
                                      processed_folder,
                                      sat_param)
    out = open_ds(output_file)

    processed_folder_bench = "Processed_Fields_files_Benchmark"    
    bench_file = pg.find_processed_file(subfolder_sat_byday_mf,
                                     processed_folder_bench,
                                     sat_param)
    out_bench = open_ds(bench_file) # Need to open comparison file and put in the same format

    data_vars = out.data_vars
    
    for dv in data_vars:
        assert np.array_equal(out[dv].values,out_bench[dv].values) 
