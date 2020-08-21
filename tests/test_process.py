#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 15:52:02 2018

Simple tests (for now) for the process script after being updated to process satellite
data grouped by day (as well as by point).

NOTE:
Tests which compare the current output to the output from the previous process script which is contained 
with the deprecated/ folder have now been skipped as they do not adhere to the new unit checking 
and conversion (g/m2/s or ppm s output).
The previous process script was originally used as a benchmark as this output has been compared
with the output produced by Alistair.

The comparison tests compare:
    - outputs from several functions for satellite by point processing
    - output from process function for site data

To run all tests:
    $ pytest test_process.py
To only run tests which check against a benchmark file:
    $ pytest test_process.py -m "bench"

@author: rt17603
"""
from __future__ import print_function

from builtins import zip
from builtins import str
import pytest
import os
import sys
import glob
import numpy as np
import xarray as xray
import acrg_name.process as process
import deprecated.process as process_org
from acrg_name.name import open_ds
import datetime
import pandas as pd
import pdb

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg

#%%  Setting up directory fixtures

@pytest.fixture()
def satellite_byday_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay/")
    return directory

@pytest.fixture()
def satellite_bypoint_directory():
    ''' Define base directory containing satellite footprint files with separated data points '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByPoint/")
    return directory

@pytest.fixture()
def site_directory():
    ''' Define base directory containing satellite footprint files with separated data points '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Site/")
    return directory

@pytest.fixture()
def satellite_byday_mf_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Satellite_ByDay_MF/")
    return directory

@pytest.fixture()
def benchmark_output_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/output/")
    return directory

@pytest.fixture()
def satellite_param():
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
def site_param():
    '''
    Define input parameters for site data to test.
    Keys include:
        "domain" (str),"site" (str),"height" (str),"year" (int),"month" (int),"satellite" (bool)
    '''
    param = {}
    param["domain"] = "CARIBBEAN"
    param["site"] = "RPB"
    param["height"] = "10magl"
    param["year"] = 2010
    param["month"] = 1
    param["satellite"] = False
   
    return param

def define_subfolder(domain,site,height):
    ''' Define subfolder name based on domain, site and height '''
    subfolder = "{domain}_{site}_{height}/".format(domain=domain,site=site,height=height)
    return subfolder

@pytest.fixture()
def subfolder_satellite(satellite_param):
    ''' Define subfolder name for satellite data'''
    domain = satellite_param["domain"]
    site = satellite_param["site"]
    height = satellite_param["height"]
    
    subfolder = define_subfolder(domain,site,height)
    return subfolder

@pytest.fixture()
def subfolder_satellite_byday(satellite_byday_directory,subfolder_satellite):
    ''' Define full subfolder path for satellite data with points grouped by day'''
    subfolder = os.path.join(satellite_byday_directory,subfolder_satellite)  
    return subfolder

@pytest.fixture()
def subfolder_satellite_bypoint(satellite_bypoint_directory,subfolder_satellite):
    ''' Define full subfolder path for satellite data with separate points'''
    subfolder = os.path.join(satellite_bypoint_directory,subfolder_satellite)  
    return subfolder

@pytest.fixture()
def subfolder_satellite_byday_mf(satellite_byday_mf_directory,subfolder_satellite):
    ''' Define full subfolder path for satellite data with points grouped by day'''
    subfolder = os.path.join(satellite_byday_mf_directory,subfolder_satellite)  
    return subfolder

@pytest.fixture()
def subfolder_site(site_param,site_directory):
    ''' Define full subfolder path for site data. '''
    domain = site_param["domain"]
    site = site_param["site"]
    height = site_param["height"]
    
    subfolder_site = define_subfolder(domain,site,height)
    subfolder = os.path.join(site_directory,subfolder_site)  
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
    param["mixr_fields_folder"] = "MixR_files"
    param["mixr_hourly_folder"] = "MixR_hourly"
    param["particles_folder"] = "Particle_files"
    param["met_folder"] = "Met_daily"
    param["processed_folder"] = "Processed_Fields_files"
    param["obs_folder"] = "Observations"

    return param    

#%% Testing read_fields() function

def get_fields_prefix(subfolder,fields_folder):
    '''
    Create prefix for fields_folder. This is a string of the path based on the subfolder and 
    fields_folder name.
    Example: /data/rt17603/NAME_output/SOUTHAMERICA_GOSAT-BRAZIL_column/Fields_files/
    Args:
        subfolder (str):
            Path to subfolder of NAME output. The final folder is based on the domain_site_height/ syntax
            e.g. SOUTHAMERICA_GOSAT-BRAZIL_column
        fields_folder (str) :
            Folder containing the fields files. e.g. "Fields_folder"
    Returns:
        str:
            Full path to fields files
    '''
    prefix = os.path.join(subfolder,fields_folder)
    return os.path.join(prefix,'')

def get_fields_files(subfolder,fields_folder,datestr):
    ''' Get filenames of all fields files matching datestr (year + month) '''
    fields_prefix = get_fields_prefix(subfolder,fields_folder)
    search_str = "*{date}*.txt*".format(date=datestr)
    
    fields_files = sorted(glob.glob(os.path.join(fields_prefix,search_str)))
    
    return fields_files

def create_datestr(param):
    ''' Create basic date string based on year and month of the format "YYYYMM" '''
    return str(param["year"]) + str(param["month"]).zfill(2)

@pytest.fixture()
def get_fields_files_satellite_byday(subfolder_satellite_byday,folder_names,satellite_param):
    '''
    Get filenames of fields files for satellite run with points grouped by day.
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(satellite_param)
    fields_folder = folder_names["fields_folder"]
    
    fields_files = get_fields_files(subfolder_satellite_byday,fields_folder,datestr)
    
    return fields_files

@pytest.fixture()
def get_fields_files_satellite_bypoint(subfolder_satellite_bypoint,folder_names,satellite_param):
    '''
    Get filenames of fields files for satellite run with separate points.
    Finds field files based on datestr (see create_datestr), multiple files per day
    (one for each point, filename syntax *YYYYMMDD-NNN*).
    '''
    datestr = create_datestr(satellite_param)
    fields_folder = folder_names["fields_folder"]
    
    fields_files = get_fields_files(subfolder_satellite_bypoint,fields_folder,datestr)
    
    return fields_files

@pytest.fixture()
def get_fields_files_site(subfolder_site,folder_names,site_param):
    ''' 
    Get filenames of fields files for satellite run with separate points. 
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(site_param)
    fields_folder = folder_names["fields_folder"]
    
    fields_files = get_fields_files(subfolder_site,fields_folder,datestr)
    return fields_files

@pytest.fixture()
def get_mixr_files_site(subfolder_site,folder_names,site_param):
    ''' 
    Get filenames of fields files for satellite run with separate points. 
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(site_param)
    fields_folder = folder_names["mixr_fields_folder"]
    
    fields_files = get_fields_files(subfolder_site,fields_folder,datestr)
    return fields_files

@pytest.fixture()
def get_mixr_hourly_site(subfolder_site,folder_names,site_param):
    ''' 
    Get filenames of fields files for satellite run with separate points. 
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(site_param)
    fields_folder = folder_names["mixr_hourly_folder"]
    
    fields_files = get_fields_files(subfolder_site,fields_folder,datestr)
    return fields_files

@pytest.fixture()
def get_fields_files_satellite_byday_mf(subfolder_satellite_byday_mf,folder_names,satellite_param):
    '''
    Get filenames of fields files for satellite run with points grouped by day.
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(satellite_param)
    fields_folder = folder_names["fields_folder"]
    
    fields_files = get_fields_files(subfolder_satellite_byday_mf,fields_folder,datestr)
    
    return fields_files

@pytest.fixture()
def mixr_file_information():
    ''' Define values in the test fields file '''
    
    header = {'Title': 'Back_General',
             'Run time': '1128UTC 17/05/2020',
             'Met data': 'NWP Flow.Regional_flow',
             'Start of release': '0000UTC 02/01/2010',
             'End of release': '0000UTC 01/01/2010',
             'Release rate': 'Multiple Sources',
             'Release location': 'Multiple Sources',
             'Release height': 'Multiple Sources',
             'Forecast duration': '744 hours',
             'X grid origin': -100.188,
             'Y grid origin': -26.126,
             'X grid size': 391,
             'Y grid size': 340,
             'X grid resolution': 0.352,
             'Y grid resolution': 0.234,
             'Number of fields': 2}
    
    column_headings = {'species_category': ['', '', '', '', 'INERT', 'INERT'],
                         'species': ['', '', '', '', 'INERT_C_00', 'INERT_C_01'],
                         'cell_measure': ['','','','','744 hr time integrated','744 hr time integrated'],
                         'quantity': ['', '', '', '', 'Mixing ratio', 'Mixing ratio'],
                         'unit': ['', '', '', '', 'ppms', 'ppms'],
                         'z_level': ['','','','','From     0 -    40m agl','From     0 -    40m agl'],
                         'time': ['X grid','Y grid','Longitude','Latitude',
                          datetime.datetime(2009, 12, 2, 0, 0),
                          datetime.datetime(2009, 12, 2, 0, 0)]}
    
    data_arrays = []
    data_arrays.append(np.zeros((340,391)))
    data_arrays[0][131,0] = 1
    data_arrays[0][135,1] = 2
    data_arrays.append(np.zeros((340,391)))
    data_arrays[1][131,0] = 3
    data_arrays[1][135,1] = 4
    
    namever = 2
    
    lons = np.arange(-100.188, -100.188+0.352*391-0.01, 0.352) + 0.352/2
    lats = np.arange(-26.126 , -26.126 +0.234*340-0.01, 0.234) + 0.234/2
    
    levs = ['From     0 -    40m agl']
    
    time = [datetime.datetime(2010, 1, 1, 0, 0), datetime.datetime(2010, 1, 1, 12, 0)]
    
    timeStep = 12
    
    return header, column_headings, data_arrays, namever, lons, lats, levs, time, timeStep

@pytest.fixture()
def particle_file_benchmark():
    '''Benchmark infromation from the dummy particle file'''
    pl_n = np.zeros((2,1,391,20))
    pl_w = np.zeros((2,1,340,20))
    mean_age_n = np.zeros((2,1,391,20))
    mean_age_w = np.zeros((2,1,340,20))
    
    pl_n[0,0,66,1] = 0.25
    pl_n[0,0,89,0] = 0.25
    pl_n[1,0,87,0] = 0.25
    pl_n[1,0,95,1] = 0.25
    
    pl_w[0,0,197,1] = 0.25
    pl_w[0,0,207,0] = 0.25
    pl_w[1,0,197,0] = 0.25
    pl_w[1,0,213,1] = 0.25
    
    mean_age_n[0,0,66,1] = 100
    mean_age_n[0,0,89,0] = 120
    mean_age_n[1,0,87,0] = 100
    mean_age_n[1,0,95,1] = 120
    
    mean_age_w[0,0,197,1] = 100
    mean_age_w[0,0,207,0] = 120
    mean_age_w[1,0,197,0] = 100
    mean_age_w[1,0,213,1] = 120
    
    return pl_n, pl_w, mean_age_n, mean_age_w

@pytest.fixture()
def met_file_benchmark():
    '''
    Defining dummy values for the satellite met data. Should be consistent with
    the dummy satellite met data which we can check against
    
    TODO: May want to read this in in some way rather than just define values
    '''
    values = {}
    
    values["release_lon"] = [-59.43330,-59.43330]
    values["release_lat"] = [13.16670,13.16670]
    values["temp"] = [27.07520,26.98370]
    values["PBLH"] = [778.9576,782.4219]
    values["press"] = [101161.4,101163.7]
    values["wind"] = [5.697351,5.695338]
    values["wind_direction"] = [72.90572,76.33878]

    df = pd.DataFrame(values)

    return df
    
def test_read_mixr_file(mixr_file_information, get_mixr_files_site):
    ''' Test read_file function '''

    # true values defined based on input file
    header_bench, column_headings_bench, data_arrays_bench, namever_bench, lons_bench, lats_bench, levs_bench, time_bench, timeStep_bench = mixr_file_information

    fields_file_1 = get_mixr_files_site[0]

    header, column_headings, data_arrays, namever = process.read_file(fields_file_1)
    
    assert np.array_equal(data_arrays, data_arrays_bench)
    assert header == header_bench
    
    for ii in range(len(data_arrays_bench)):
        assert np.array_equal(data_arrays[ii],data_arrays_bench[ii])
        
    assert namever == namever_bench

@pytest.fixture()
def read_fields_file_satellite_byday(get_fields_files_satellite_byday):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_satellite_byday[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    #print "header,column_headings",header,column_headings
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_fields_file_satellite_bypoint(get_fields_files_satellite_bypoint):
    ''' Read first satellite by point field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_satellite_bypoint[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_fields_file_site(get_fields_files_site):
    ''' Read first site field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_site[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_mixr_file_site(get_fields_files_site):
    ''' Read first site field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_site[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

@pytest.fixture()
def read_fields_file_satellite_byday_mf(get_fields_files_satellite_byday_mf):
    ''' Read first satellite by day field file using process.read_file() function and create output. '''
    fields_file = get_fields_files_satellite_byday_mf[0]
    header, column_headings, data_arrays = process.read_file(fields_file)
    
    return header, column_headings, data_arrays

def test_read_field_file_satellite_bypoint(get_fields_files_satellite_bypoint):
    ''' Test read_file function can run for satellite data separated by point. '''
    point = 10
    fields_file = get_fields_files_satellite_bypoint[point-1]
    header, column_headings, data_arrays = process.read_file(fields_file)

def test_read_field_file_satellite_byday(get_fields_files_satellite_byday):
    '''  Test read_file function can run for satellite data separated by day. '''
    fields_file = get_fields_files_satellite_byday[0]
    header, column_headings, data_arrays = process.read_file(fields_file)

#    point = 10
#    da_range = [(point-1)*17,(point-1)*17+17]
#    data_arrays_point = data_arrays[da_range[0]:da_range[1]]
#    for i,array in enumerate(data_arrays_point):
#        print i,np.min(array),np.max(array),np.mean(array),np.sum(array)
#    print "Over 17 elements of point {0}".format(point),np.mean(data_arrays_point)


    ### TODO: ADD MORE STRINGENT TEST    

def test_define_grid_satellite_byday(read_fields_file_satellite_byday,satellite_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with points 
    grouped by day.
    '''
    header,column_headings,data_arrays = read_fields_file_satellite_byday
    out = process.define_grid(header,column_headings,satellite=True,upper_level=satellite_param["upper_level"])
    
    assert out

    ### TODO: ADD MORE STRINGENT TEST

def test_define_grid_satellite_bypoint(read_fields_file_satellite_bypoint,satellite_param):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over satellite data with separate 
    points.
    '''
    header,column_headings,data_arrays = read_fields_file_satellite_bypoint
    out = process.define_grid(header,column_headings,satellite=True,upper_level=satellite_param["upper_level"])

    assert out
    ### TODO: ADD MORE STRINGENT TEST

def test_define_grid_site(mixr_file_information):
    '''
    Test that grid can be defined correctly when extracted from a NAME run over site data.
    '''
    
    header_bench, column_headings_bench, data_arrays_bench, namever_bench, lons_bench, lats_bench, levs_bench, time_bench, timeStep_bench = mixr_file_information
    
    lons, lats, levs, time, timeStep = process.define_grid(namever_bench, header_bench,column_headings_bench,satellite=False)

    assert np.array_equal(lons,lons_bench)
    assert np.array_equal(lats,lats_bench)
    assert levs == levs_bench
    assert time == time_bench
    assert timeStep == timeStep_bench 
    

def test_satellite_byday_units(read_fields_file_satellite_byday):
    '''
    Test units column heading for satellite by day input is being read in correctly (expect g s/m^3)
    '''    
    header,column_headings,data_arrays = read_fields_file_satellite_byday
    units_column = column_headings["unit"]
    units = units_column[4:][0] # First 4 headings are empty
    print("units",units)
    print(units.replace(' ',''))
    assert units.replace(' ','') == "gs/m^3" or units.replace(' ','') == "gs/m3"

def test_satellite_byday_mf_units(read_fields_file_satellite_byday_mf):
    '''
    Test units column heading for mol fraction input is being read in correctly (expect ppm s)
    '''    
    header,column_headings,data_arrays = read_fields_file_satellite_byday_mf
    units_column = column_headings["unit"]
    units = units_column[4:][0] # First 4 headings are empty
    assert units.replace(' ','') == "ppms"
 
#%%  Test particle_locations()
  
@pytest.fixture()
def define_grid_satellite_byday(read_fields_file_satellite_byday,satellite_param):
    ''' Create output from process.define_grid() function for one satellite by day field file '''
    header,column_headings,data_arrays = read_fields_file_satellite_byday
    upper_level=satellite_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def define_grid_satellite_bypoint(read_fields_file_satellite_bypoint,satellite_param):
    ''' Create output from process.define_grid() function for one satellite by point field file '''
    header,column_headings,data_arrays = read_fields_file_satellite_bypoint
    upper_level=satellite_param["upper_level"]
    lons, lats, levs, time, timeStep = process.define_grid(header,column_headings,satellite=True,upper_level=upper_level)
    
    return lons, lats, levs, time, timeStep

@pytest.fixture()
def define_grid_site(read_fields_file_site):
    ''' Create output from process.define_grid() function for one site field file '''
    header,column_headings,data_arrays = read_fields_file_site
    lons, lats, levs, time, timeStep = process.define_grid(header,column_headings,satellite=False)
    
    return lons, lats, levs, time, timeStep

def get_field_file_datestr(field_files):
    ''' 
    Get date strings for file searching from fields files. This takes the last section of the filename when
    split by "_" and extension and path are removed.
    e.g. for Site files "Fields_MACE_HEAD_20120331.txt.gz" --> "20120331"
         for ByDay files "Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_20120101.txt.gz" --> "20120101"
         for ByPoint files "Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_ByPoint_20120101-029.txt.gz" --> "20120101-029"
    Args:
        fields_files (list):
            List of field files. Will extract date strings from these files.
    Returns:
        list:
            All date strings for the input field_files
    '''
    file_datestrs = [os.path.split(f)[1].split(".txt")[0].split("_")[-1] for f in field_files]
    return file_datestrs

def get_particle_prefix(subfolder,particles_folder):
    '''
    Create prefix for particles_folder.
    This is a string of the path based on the subfolder and particles_folder name.
    Example: /data/rt17603/NAME_output/SOUTHAMERICA_GOSAT-BRAZIL_column/Particle_files/
    Args:
        subfolder (str):
            Path up to subfolder of NAME output. The final folder is based on the domain_site_height/ syntax
            e.g. SOUTHAMERICA_GOSAT-BRAZIL_column
        particles_folder (str) :
            Folder containing the particle location files. e.g. "Particle_files"
    Returns:
        str:
            Full path to particle location files
    '''
    prefix = os.path.join(subfolder,particles_folder)
    return os.path.join(prefix,'')

def get_particle_files(subfolder,particles_folder,field_files):
    '''
    Get filenames of all particle files matching date strings extracted from fields files.
    Returns a list of filenames.
    '''
    particle_prefix = get_particle_prefix(subfolder,particles_folder)
    
    file_datestrs = get_field_file_datestr(field_files) # Can be YYYYMMDD or YYYYMMDD-NNN
    
    particle_files = []
    for file_datestr in file_datestrs:
        particle_file_search_string = os.path.join(particle_prefix,"*{datestr}*.txt*".format(datestr=file_datestr))
        particle_file_search = glob.glob(particle_file_search_string)
        if particle_file_search:
            particle_files.extend(particle_file_search)

    return particle_files

@pytest.fixture()
def get_particle_files_satellite_byday(subfolder_satellite_byday,folder_names,
                                       get_fields_files_satellite_byday):
    ''' 
    Get filenames of particle files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = get_particle_files(subfolder_satellite_byday,particles_folder,
                                        get_fields_files_satellite_byday)
    return particle_files

@pytest.fixture()
def get_particle_files_satellite_bypoint(subfolder_satellite_bypoint,folder_names,
                                         get_fields_files_satellite_bypoint):
    ''' 
    Get filenames of particle files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = get_particle_files(subfolder_satellite_bypoint,particles_folder,
                                        get_fields_files_satellite_bypoint)
    return particle_files

@pytest.fixture()
def get_particle_files_site(subfolder_site,folder_names,
                                         get_mixr_files_site):
    ''' 
    Get filenames of particle files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = get_particle_files(subfolder_site,particles_folder,
                                        get_mixr_files_site)
    return particle_files

@pytest.fixture()
def define_heights():
    ''' Define height range for input into particle_locations function. '''
    dheights = 1000
    heights = np.arange(0, 19001, dheights) + dheights/2.
    return heights

def test_particle_locations_satellite_byday(define_grid_satellite_byday,get_particle_files_satellite_byday,
                                            define_heights,satellite_param):
    '''
    Test particle_locations() function can produce the correct output for satellite data when points have been
    grouped by day.
    '''
    lons, lats, levs, time, timeStep = define_grid_satellite_byday
    particle_file = get_particle_files_satellite_byday[0]
    heights =  define_heights   
    upper_level = satellite_param["upper_level"]
    
    out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
    assert len(out["time"].values) > 1
    assert out
    ### TODO: ADD MORE STRINGENT TEST
   
def test_particle_locations_satellite_bypoint(define_grid_satellite_bypoint,get_particle_files_satellite_bypoint,
                                            define_heights,satellite_param):
    '''
    Test particle_locations() function can produce the correct output for satellite data when points are separate.
    '''
    lons, lats, levs, time, timeStep = define_grid_satellite_bypoint
    particle_files = get_particle_files_satellite_bypoint
    heights =  define_heights   
    upper_level = satellite_param["upper_level"]
    
    for particle_file in particle_files:
        out = process.particle_locations(particle_file,time,lats,lons,levs,heights,id_is_lev=False,
                                     satellite=True,upper_level=upper_level)
    
        assert len(out["time"].values) == 1
        assert out
        ### TODO: ADD MORE STRINGENT TEST

def test_particle_locations_site(particle_file_benchmark,mixr_file_information,get_particle_files_site,define_heights):
    '''
    Test particle_locations() function can produce the correct output for a. dummy particle location file.
    '''

    header_bench, column_headings_bench, data_arrays_bench, namever_bench, lons_bench, lats_bench, levs_bench, time_bench, timeStep_bench = mixr_file_information

    particle_file_1 = get_particle_files_site[0]
    heights =  define_heights
    
    pl_n_bench, pl_w_bench, mean_age_n_bench, mean_age_w_bench = particle_file_benchmark
    

    hist = process.particle_locations(particle_file_1,time_bench,lats_bench,lons_bench,levs_bench, heights,id_is_lev=False,
                                 satellite=False)

    assert np.array_equal(hist["pl_n"].values,pl_n_bench)
    assert np.array_equal(hist["pl_w"].values,pl_w_bench)
    assert np.array_equal(hist["mean_age_n"].values,mean_age_n_bench)
    assert np.array_equal(hist["mean_age_w"].values,mean_age_w_bench)        

@pytest.fixture()
def define_grid_org_satellite_bypoint(read_fields_file_satellite_bypoint):
    ''' 
    Create output from original process script: process_org.define_grid() function for one satellite by point 
    field file.
    For comparison with new process script.
    '''
    header,column_headings,data_arrays = read_fields_file_satellite_bypoint
    lons, lats, levs, time, timeStep = process_org.define_grid(header,column_headings,satellite=True)
    
    return lons, lats, levs, time, timeStep

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_particle_locations_satellite_bypoint_against_org(define_grid_satellite_bypoint,
                                                         define_grid_org_satellite_bypoint,
                                                         get_particle_files_satellite_bypoint,
                                                         define_heights,satellite_param):
    '''
    Test particle locations output is identical between particle files processed using the original
    process script and the new process script. 
    This is for satellite output separated by point since original script cannot process by day output.
    '''
    lons1, lats1, levs1, time1, timeStep1 = define_grid_satellite_bypoint
    lons2, lats2, levs2, time2, timeStep2 = define_grid_org_satellite_bypoint
    particle_files = get_particle_files_satellite_bypoint
    heights =  define_heights   
    upper_level = satellite_param["upper_level"]
    
    pl_data_vars = ["pl_n","pl_e","pl_s","pl_w"]
    
    for i,particle_file in enumerate(particle_files):
        out = process.particle_locations(particle_file,time1,lats1,lons1,levs1,heights,id_is_lev=False,satellite=True,
                                         upper_level=upper_level)
        out_org = process_org.particle_locations(particle_file,time2,lats2,lons2,levs2,heights,id_is_lev=True)
    
        for pl in pl_data_vars:
            assert np.array_equal(out[pl].values,out_org[pl].values)

#%% Test read_met()

def get_met_prefix(subfolder,met_folder):
    '''
    Create prefix for met_folder.
    This is a string of the path based on the subfolder and met_folder name.
    Example: /data/rt17603/NAME_output/SOUTHAMERICA_GOSAT-BRAZIL_column/Met/
    Args:
        subfolder (str):
            Path up to subfolder of NAME output. The final folder is based on the domain_site_height/ syntax
            e.g. SOUTHAMERICA_GOSAT-BRAZIL_column
        met_folder (str) :
            Folder containing the met files or folders. e.g. "Met"
    Returns:
        str:
            Full path to met file or folders containing separated met files.
    '''
    prefix = os.path.join(subfolder,met_folder)
    return os.path.join(prefix,'')

def get_met_files(subfolder,met_folder,field_files=None,satellite=False):
    '''
    Get filenames of met files within full path to the met folder.
    For site data all files from this folder are extracted.
    For satellite data file names are exracted from within folder with matching date strings from fields files.
    Args:
        subfolder (str) :
            See get_met_prefix
        met_folder (str) :
            See get_met_prefix
        fields_files (list / None) :
            Only needed if satellite = True.
            List of filenames (str) for fields files. Used to extract date string values (datestr).
        satellite (bool) :
            If satellite is True, met data is extracted based on search str "*datestr/*.txt*" 
            If satellite is False, met data is extracted based on search str *.txt*
    Returns:
        If satellite is True:
            nested list:
                list of filenames, one list for each date string.
        If satellite is False:
            list:
                list of filenames
    '''
    met_prefix = get_met_prefix(subfolder,met_folder)
    
    if satellite:
        file_datestrs = get_field_file_datestr(field_files) # Can be YYYYMMDD or YYYYMMDD-NNN

        all_met_files = []
        for file_datestr in file_datestrs:
            met_search_str = os.path.join(met_prefix,"*{datestr}/*.txt*".format(datestr=file_datestr))
            met_files = glob.glob(met_search_str)
            all_met_files.append(met_files)
    else:
        met_search_str = subfolder + met_folder + "/*.txt*"
        met_files = glob.glob(met_search_str)
        all_met_files = [met_files]
    
    return all_met_files

@pytest.fixture()   
def get_met_files_satellite_byday(subfolder_satellite_byday,folder_names,get_fields_files_satellite_byday):
    '''
    Get filenames of met files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    met_folder = folder_names["met_folder"]
    met_files = get_met_files(subfolder_satellite_byday,met_folder,get_fields_files_satellite_byday,satellite=True)
    
    return met_files

@pytest.fixture()
def get_met_files_satellite_bypoint(subfolder_satellite_bypoint,folder_names,get_fields_files_satellite_bypoint):
    '''
    Get filenames of met files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    
    met_folder = folder_names["met_folder"]
    met_files = get_met_files(subfolder_satellite_bypoint,met_folder,get_fields_files_satellite_bypoint,satellite=True)
    
    return met_files

@pytest.fixture()
def get_met_files_site(subfolder_site,folder_names):
    ''' Get filenames of met files for site run. '''
    
    met_folder = folder_names["met_folder"]
    met_files = get_met_files(subfolder_site,met_folder,satellite=False)
    
    return met_files


def test_read_met_satellite_byday(get_met_files_satellite_byday,satellite_param):
    '''
    Test read_met function can read in one met file when satellite data is grouped by day
    '''
    met_files_1 = get_met_files_satellite_byday[0]
    total_file_num = len(met_files_1)
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    assert len(out.index) == total_file_num

def test_read_met_satellite_bypoint(get_met_files_satellite_bypoint,satellite_param):
    '''
    Test read_met function can read in one met file when satellite data is separated into points.
    '''
    met_files_1 = get_met_files_satellite_bypoint[0]
    total_file_num = len(met_files_1)
   
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)

    assert len(out.index) == total_file_num

def test_read_met_site(get_met_files_site, met_file_benchmark, read_met_site):
    '''
    Test read_met function can read in site met data successfully.
    '''
    met_files_1 = get_met_files_site[0]
    
    out = process.read_met(met_files_1,vertical_profile=False, satellite=False)
    
    expected = met_file_benchmark
    compare_columns = expected.columns
    
    for col in compare_columns:
        assert np.array_equal(out[col],expected[col])

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_read_met_satellitebypoint_against_org(get_met_files_satellite_bypoint,satellite_param):
    '''
    Test values within output from original process script: process_org.read_met() function is identical to
    met data when processed using new process script.
    '''
    met_files_1 = get_met_files_satellite_bypoint[0]
    
    out = process.read_met(met_files_1,satellite=True,vertical_profile=False)
    
    out_org_list = []
    for met_file in met_files_1:
        out_org_list.append(process_org.read_met(met_file))
    
    if "label" in out.columns:
        out.drop("label",inplace=True,axis=1)
    
    for t,out_org in enumerate(out_org_list):
        assert np.all(out.iloc[t] == out_org)
        

#%% Test footprint_array()

def read_met_files(met_files,satellite=False):
    '''
    Reads in met files based on extracted met filenames. Expects nested list of met files and processes each
    list of met_files together so they are output as one dataframe per set of met files.
    Args:
        met_files (nested list):
            List of each set of met files. One set of met files extracted for each input date string 
            (date string is based on fields filenames).
        satellite (bool, optional):
            Whether to process as satellite data or site data.
    Returns:
        list(pandas.Dataframe):
            List of dataframes, one for each set of met_files
    '''
    df_all = [process.read_met(fnames,satellite=satellite) for fnames in met_files]
    return df_all

def read_met_files_separate(met_files,satellite=False):
    '''
    Reads in met files based on extracted met filenames. Expects nested list of met files but processes each
    met_file separately so they are output as one dataframe per met file.
    Args:
        See read_met_files function
    Returns:
        nested list(pandas.Dataframe):
            List of dataframes, one list for set of met files, one dataframe for every met file.
    '''
    df_sep = [[process.read_met(fname,satellite=satellite) for fname in fnames] for fnames in met_files]
    return df_sep

@pytest.fixture()
def read_met_satellite_byday(get_met_files_satellite_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Returns: list of dataframes.
    '''
    return read_met_files(get_met_files_satellite_byday,satellite=True)

@pytest.fixture()
def read_met_satellite_bypoint(get_met_files_satellite_bypoint):
    '''
    Read met data for satellite data is separated into points.
    Returns: list of dataframes.
    '''
    return read_met_files(get_met_files_satellite_bypoint,satellite=True)

@pytest.fixture()
def read_met_separate_satellite_byday(get_met_files_satellite_byday):
    '''
    Read met data for satellite run with points grouped by day.
    Separates points into individual dataframes.??
    Returns: nested list of dataframes.
    '''
    return read_met_files_separate(get_met_files_satellite_byday,satellite=True)

@pytest.fixture()
def read_met_separate_satellite_bypoint(get_met_files_satellite_bypoint):
    '''
    Read met data for satellite run is separated into points.
    Separates points into individual dataframes.
    This is the input format expected by the original process script.
    Returns: nested list of dataframes.
    '''
    return read_met_files_separate(get_met_files_satellite_bypoint,satellite=True)

@pytest.fixture()
def read_met_site(get_met_files_site):
    ''' Read met data for site data. '''
    return read_met_files(get_met_files_site, satellite=False)

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
    Get filenames of all particle files matching date strings extracted from fields files.
    Only relevant to satellite data where we have to extract the times from the observations files rather 
    than work it out from the NAME input.
    Returns a list of filenames.
    '''
    
    obs_prefix = get_obs_prefix(subfolder,obs_folder)
    
    file_datestrs = get_field_file_datestr(field_files) # Can be YYYYMMDD or YYYYMMDD-NNN
    
    obs_files = []
    for file_datestr in file_datestrs:
        obs_search_string = os.path.join(obs_prefix,"*{datestr}_*.nc".format(datestr=file_datestr))
        obs_file_search = glob.glob(obs_search_string)
        if obs_file_search:
            obs_files.extend(obs_file_search)

    return obs_files

@pytest.fixture()
def get_obs_files_satellite_byday(subfolder_satellite_byday,folder_names,
                                  get_fields_files_satellite_byday):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_satellite_byday,
                              obs_folder,get_fields_files_satellite_byday)
    
    return obs_files

@pytest.fixture()
def get_obs_files_satellite_bypoint(subfolder_satellite_bypoint,folder_names,
                                  get_fields_files_satellite_bypoint):
    '''
    Get filenames of observation files for satellite run with points grouped by day.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    obs_folder = folder_names["obs_folder"]
    
    obs_files = get_obs_files(subfolder_satellite_bypoint,
                              obs_folder,get_fields_files_satellite_bypoint)
    
    return obs_files
 
def get_time_step(subfolder):
    '''
    Extract value from time step file, if it exists. Should be a "time_step.txt" file within the NAME subfolder.
    File should only contain a number which is the number of seconds of the retrieval_time.
    Only relevant to satellite data where we may have a release time < 1 minute. NAME maximum time resolution
    of one minute means the time step cannot be extracted from the output files.
    '''
    time_step_file = os.path.join(subfolder,"time_step.txt")

    if os.path.exists(time_step_file):
        with open(time_step_file) as f:
            timeStep = float(f.read())/3600.
    else:
        timeStep = None
    
    return timeStep

@pytest.fixture()
def get_time_step_satellite_byday(subfolder_satellite_byday):
    ''' Extract time step value from within satellite_byday subfolder '''
    return get_time_step(subfolder_satellite_byday)

@pytest.fixture()
def get_time_step_satellite_bypoint(subfolder_satellite_bypoint):
    ''' Extract time step value from within satellite_byday subfolder '''
    return get_time_step(subfolder_satellite_bypoint)

def test_footprint_array_satellite_byday(get_fields_files_satellite_byday,
                                         get_particle_files_satellite_byday,
                                         read_met_satellite_byday,
                                         get_obs_files_satellite_byday,
                                         get_time_step_satellite_byday,
                                         satellite_param):
    '''
    Test footprint_array function can return an output for satellite files when grouped by day.
    '''
    fields_file_1 = get_fields_files_satellite_byday[0]
    particle_file_1 = get_particle_files_satellite_byday[0]
    met = read_met_satellite_byday[0]
    #obs_file_1 = get_obs_files_satellite_byday[0]
    time_step = get_time_step_satellite_byday
    upper_level = satellite_param["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,upper_level=upper_level,
                                  time_step=time_step)

def test_footprint_array_satellite_bypoint(get_fields_files_satellite_bypoint,
                                           get_particle_files_satellite_bypoint,
                                           read_met_satellite_bypoint,
                                           get_obs_files_satellite_bypoint,
                                           get_time_step_satellite_bypoint,
                                           satellite_param):
    '''
    Test footprint_array function can return an output for satellite files when separated by point.
    '''
    point = 11
    
    fields_file_1 = get_fields_files_satellite_bypoint[point-1]
    particle_file_1 = get_particle_files_satellite_bypoint[point-1]
    met = read_met_satellite_bypoint[point-1]
    #obs_file_1 = get_obs_files_satellite_bypoint[0]
    time_step = get_time_step_satellite_bypoint
    upper_level = satellite_param["upper_level"]
    
    out = process.footprint_array(fields_file_1,particle_file_1,met,satellite=True,upper_level=upper_level,
                                  time_step=time_step)
    
    print(out["wind_direction"])

def test_footprint_array_site(benchmark_output_directory,get_mixr_files_site,
                              get_particle_files_site,
                              get_met_files_site):
    '''
    Test footprint_array function can return an output for site files.
    '''
    
    fp_bench_file = os.path.join(benchmark_output_directory,'Test_process_footprint_array.nc')
    
    with xray.open_dataset(fp_bench_file) as f: 
        fp_bench = f.load()
    
    fields_file_1 = get_mixr_files_site[0]
    particle_file_1 = get_particle_files_site[0]
    
    met_files_1 = get_met_files_site[0]
    
    met = process.read_met(met_files_1,vertical_profile=False, satellite=False)
    
    out = process.footprint_array(fields_file_1, particle_file_1, met)
    
    for key in fp_bench.keys():
        assert np.array_equal(fp_bench[key], out[key])

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_footprint_array_satellite_bypoint_against_org(get_fields_files_satellite_bypoint,
                                           get_particle_files_satellite_bypoint,
                                           read_met_satellite_bypoint,
                                           read_met_separate_satellite_bypoint,
                                           get_obs_files_satellite_bypoint,
                                           get_time_step_satellite_bypoint,
                                           satellite_param):
    '''
    Test output from original process script: process_org.footprint_array() produces an identical output
    to the new process script over all files extracted based on field files datestr values.
    '''
    fields_files_1 = get_fields_files_satellite_bypoint
    particle_files_1 = get_particle_files_satellite_bypoint
    met_1 = read_met_satellite_bypoint
    met_2 = read_met_separate_satellite_bypoint
    time_step = get_time_step_satellite_bypoint
    upper_level = satellite_param["upper_level"]
    
    for fields_file,particle_file,met_data_1,met_data_2 in zip(fields_files_1,particle_files_1,met_1,met_2):
        
        out = process.footprint_array(fields_file,particle_file,met_data_1,satellite=True,time_step=time_step,upper_level=upper_level)
        out_org = process_org.footprint_array(fields_file,particle_file,met_data_2,satellite=True,time_step=time_step)
        
        data_vars = out.data_vars
        
        for dv in data_vars:
            assert np.array_equal(out[dv].values,out_org[dv].values)
        
#%% Test footprint_concatenate()

def footprint_concatenate_param(subfolder,read_met,folder_names,parameters,satellite,
                                get_time_step=None,get_obs_files=None):
    '''
    Define parameters for input into footprint_concatenate() function. 
    Args:
        subfolder (str) :
            Path to subfolder of NAME output. The final folder is based on the domain_site_height/ syntax
            e.g. SOUTHAMERICA_GOSAT-BRAZIL_column.
        read_met (list) :
            Output from any read_met_files() or read_met_files_separate function. Should contain list of
            dataframes (or nested list of dataframes)
        folder_names (dict) :
            Dictionary of folder_names - see folder_names() fixture
        parameters (dict) :
            Dictionary of parameters relevant to run - see satellite_param() or site_param() fixtures
        satellite (bool) :
            If the input data is satellite data
        get_time_step (float / None, optional) :
            Optional if satellite=True; not needed if satellite=False.
            Time step value - see get_time_step function
        #get_obs_files (list, optional)
        #    Only needed if satellite=True.
        #    List of observations files.
    
    Returns:
        dict :
            Key-value pairs matching input into footprint_concatenate function.
            **NOTE: At the moment "datestr" key is linked to a list of datestr rather than one datestr value.
            Need to re-assign "datestr" key to one datestr value to param dictionary before passing to function **
    '''
    param = {}
    param["fields_prefix"] = get_fields_prefix(subfolder,folder_names["fields_folder"])
    param["particle_prefix"] = get_particle_prefix(subfolder,folder_names["particles_folder"])
    param["met"] = read_met
    param["satellite"] = satellite
    
    datestr = create_datestr(parameters)
    field_files = get_fields_files(subfolder,folder_names["fields_folder"],datestr)
    param["datestr"] = get_field_file_datestr(field_files)
    
    if satellite:
        #param["obs_file"] = get_obs_files[0]
        param["time_step"] = get_time_step
        param["upper_level"] = parameters["upper_level"]

    return param

@pytest.fixture()
def fc_param_satellite_byday(subfolder_satellite_byday,read_met_satellite_byday,
                             folder_names,satellite_param,
                             get_time_step_satellite_byday,
                             get_obs_files_satellite_byday):
    '''
    Define parameters for input into footprint_concatenate function for satellite files when grouped by day.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = footprint_concatenate_param(subfolder_satellite_byday,read_met_satellite_byday,
                                        folder_names,satellite_param,satellite=True,
                                        get_time_step=get_time_step_satellite_byday,
                                        get_obs_files=get_obs_files_satellite_byday)
    
    return param

@pytest.fixture()
def fc_param_satellite_bypoint(subfolder_satellite_bypoint,read_met_satellite_bypoint,
                             folder_names,satellite_param,
                             get_time_step_satellite_bypoint,
                             get_obs_files_satellite_bypoint):
    '''
    Define parameters for input into footprint_concatenate function for satellite files when separated by point.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = footprint_concatenate_param(subfolder_satellite_bypoint,read_met_satellite_bypoint,
                                        folder_names,satellite_param,satellite=True,
                                        get_time_step=get_time_step_satellite_bypoint,
                                        get_obs_files=get_obs_files_satellite_bypoint)
    
    return param

@pytest.fixture()
def fc_param_site(subfolder_site,read_met_site,folder_names,site_param):
    '''
    Define parameters for input into footprint_concatenate function for site files.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = footprint_concatenate_param(subfolder_site,read_met_site,
                                        folder_names,site_param,satellite=False)
    
    return param

@pytest.fixture()
def fc_param_org_satellite_bypoint(fc_param_satellite_bypoint,read_met_separate_satellite_bypoint):
    '''
    Define parameters for input into original process script: process_org.footprint_concatenate() for
    satellite files when separated by point.
    **NOTE: At the moment "datestr" and "met" keys are linked to a lists of datestr and met data rather 
    than one value. Need to re-assign "datestr" and "met" keys to singular values in param dictionary before 
    passing to function **
    '''
    param = fc_param_satellite_bypoint.copy()
    if "upper_level" in list(param.keys()):
        param.pop("upper_level")
    param["met"] = read_met_separate_satellite_bypoint # Re-define based on orginal script expected input.
    
    return param
  
def test_footprint_concatenate_satellite_byday(fc_param_satellite_byday):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are grouped by day.
    '''
    param = fc_param_satellite_byday
    param["datestr"] = param["datestr"][0]
    
    out = process.footprint_concatenate(**param)
  
def test_footprint_concatenate_satellite_bypoint(fc_param_satellite_bypoint):
    '''
    Test footprint_concatenate function produces an output when files for satellite data are separated by point.
    '''
    param = fc_param_satellite_bypoint
    param["datestr"] = param["datestr"][0]
    param["met"] = param["met"][0]

    out = process.footprint_concatenate(**param)
    
def test_footprint_concatenate_site(fc_param_site):
    '''
    Test footprint_concatenate function produces an output when files are for site data.
    '''
    param = fc_param_site
    param["datestr"] = param["datestr"][0]

    out = process.footprint_concatenate(**param)
    
@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_footprint_concatenate_satellite_bypoint_against_org(fc_param_satellite_bypoint,
                                                             fc_param_org_satellite_bypoint):
    '''
    Test output from original process script: process_org.footprint_concatenate against output from new
    process script for a satellite run where output is separated by point.
    '''
    
    param = fc_param_satellite_bypoint.copy() 
    met_points = fc_param_satellite_bypoint["met"][:]
    date_strings = fc_param_satellite_bypoint["datestr"][:]
    
    param_org = fc_param_org_satellite_bypoint.copy()
    met_points_org = fc_param_org_satellite_bypoint["met"][:]
    date_strings_org = fc_param_satellite_bypoint["datestr"][:]

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
def footprint_concatenate_satellite_bypoint(fc_param_satellite_bypoint):
    '''
    Create output from footprint_concatenate function for satellite data seperated by point.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per point.
    '''
    param = fc_param_satellite_bypoint.copy()
    met_points = fc_param_satellite_bypoint["met"][:]
    date_strings = fc_param_satellite_bypoint["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process.footprint_concatenate(**param))
    return fp_all

@pytest.fixture()
def footprint_concatenate_satellite_byday(fc_param_satellite_byday):
    '''
    Create output from footprint_concatenate function for satellite data grouped by day.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per day.
    '''
    param = fc_param_satellite_byday.copy()
    met_points = fc_param_satellite_byday["met"][:]
    date_strings = fc_param_satellite_byday["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process.footprint_concatenate(**param))
    return fp_all

def test_satellite_vertical_profile_bypoint(footprint_concatenate_satellite_bypoint,
                                          get_obs_files_satellite_bypoint,satellite_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is separated by point.
    '''
    satellite_obs_file = get_obs_files_satellite_bypoint[0]
    fp = footprint_concatenate_satellite_bypoint[0]
    out = process.satellite_vertical_profile(fp,satellite_obs_file,max_level=satellite_param["max_level"])

def test_satellite_vertical_profile_byday(footprint_concatenate_satellite_byday,
                                          get_obs_files_satellite_byday,satellite_param):
    '''
    Test satellite_vertical_profile funcion can produce an output when data from a satellite run 
    is grouped by day.
    '''
    satellite_obs_file = get_obs_files_satellite_byday[0]
    fp = footprint_concatenate_satellite_byday[0]
    out = process.satellite_vertical_profile(fp,satellite_obs_file,max_level=satellite_param["max_level"])

@pytest.fixture()
def footprint_concatenate_org_satellite_bypoint(fc_param_org_satellite_bypoint):
    '''
    Create output from original process script: process_org.footprint_concatenate function for satellite data 
    seperated by point.
    Returns list of xarray.Dataset objects. One for each footprint, one footprint per point.
    '''
    param = fc_param_org_satellite_bypoint.copy()
    met_points = fc_param_org_satellite_bypoint["met"][:]
    date_strings = fc_param_org_satellite_bypoint["datestr"][:]
    
    fp_all = []
    for met,datestr in zip(met_points,date_strings):
        param["met"] = met
        param["datestr"] = datestr
        fp_all.append(process_org.footprint_concatenate(**param))
    
    return fp_all

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_satellite_vertical_profile_bypoint_against_org(footprint_concatenate_satellite_bypoint,
                                                        footprint_concatenate_org_satellite_bypoint,
                                                        get_obs_files_satellite_bypoint,satellite_param):
    '''
    Test output from original process script: process_org.satellite_vertical_profile matches the output 
    from new process script for satellite data separated into points.
    '''
   
    satellite_obs_files = get_obs_files_satellite_bypoint
    fp_all = footprint_concatenate_satellite_bypoint
    fp_all_org = footprint_concatenate_org_satellite_bypoint
    
    for fp,fp_org,satellite_obs_file in zip(fp_all,fp_all_org,satellite_obs_files):
        out = process.satellite_vertical_profile(fp,satellite_obs_file,max_level=satellite_param["max_level"])
        out_org = process_org.satellite_vertical_profile(fp_org,satellite_obs_file,
                                                                  max_level=satellite_param["max_level"])
        data_vars = out.data_vars
        for dv in data_vars:
            assert np.array_equal(out[dv].values,out_org[dv].values)

#%% Test process()

def process_param(input_param,folder_names,satellite):
    '''
    Define input parameters for process.process function.
    Args:
        input_param (dict) :
            Dictionary of input parameters - see satellite_param and site_param fixtures
        folder_names (dict) :
            Dictionary of folder names within subfolder - see folder_names fixture
        satellite (bool) :
            If input parameters are related to satellite data.
            If so, additional parameters for "max_level" and "upper_level" will be added.
    Returns:
        dict :
            Key-value pairs to pass to process function
    '''
    param = {}
    param["domain"] = input_param["domain"]
    param["site"] = input_param["site"]
    param["height"] = input_param["height"]
    param["year"] = input_param["year"]
    param["month"] = input_param["month"]
    if satellite:
        param["max_level"] = input_param["max_level"]
        param["upper_level"] = input_param["upper_level"]
    
    param["fields_folder"] = folder_names["fields_folder"]
    param["particles_folder"] = folder_names["particles_folder"]
    param["met_folder"] = folder_names["met_folder"]
    param["processed_folder"] = folder_names["processed_folder"]

    param["satellite"] = satellite
    
    return param
    
@pytest.fixture()
def process_satellite_byday_param(satellite_param,folder_names,satellite_byday_directory):
    '''
    Define input parameters for process.process function for satellite data grouped by day.
    Additional "base_dir" parameter added to point to Satellite_ByDay folder.
    '''
    param = process_param(satellite_param,folder_names,satellite=True)
    param["base_dir"] = satellite_byday_directory
    param["process_dir"] = os.path.join(satellite_byday_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_satellite_bypoint_param(satellite_param,folder_names,satellite_bypoint_directory):
    '''
    Define input parameters for process.process function for satellite data separated by point.
    Additional "base_dir" parameter added to point to Satellite_ByPoint folder.
    '''    
    param = process_param(satellite_param,folder_names,satellite=True)
    param["base_dir"] = satellite_bypoint_directory
    param["process_dir"] = os.path.join(satellite_bypoint_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_satellite_org_bypoint_param(satellite_param,folder_names,satellite_bypoint_directory):
    '''
    Define input parameters for process.process function for satellite data separated by point.
    Additional "base_dir" parameter added to point to Satellite_ByPoint folder.
    '''    
    param = process_param(satellite_param,folder_names,satellite=True)
    param["base_dir"] = satellite_bypoint_directory
    if "upper_level" in list(param.keys()):
        param.pop("upper_level")
    param["process_dir"] = os.path.join(satellite_bypoint_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_site_param(site_param,folder_names,site_directory):
    '''
    Define input parameters for process.process function for site data.
    Additional "base_dir" parameter added to point to Site folder.
    '''    
    param = process_param(site_param,folder_names,satellite=False)
    param["base_dir"] = site_directory
    param["process_dir"] = os.path.join(site_directory,"Processed_Fields_files")
    
    return param

@pytest.fixture()
def process_satellite_byday_mf_param(satellite_param,folder_names,satellite_byday_mf_directory):
    '''
    Define input parameters for process.process function for satellite data grouped by day.
    Additional "base_dir" parameter added to point to Satellite_ByDay folder.
    '''
    param = process_param(satellite_param,folder_names,satellite=True)
    param["base_dir"] = satellite_byday_mf_directory
    param["process_dir"] = os.path.join(satellite_byday_mf_directory,"Processed_Fields_files")
    
    return param

def find_processed_file(subfolder,processed_folder,param):
    '''
    Find NAME processed file based on search string:
        "{site}-{height}_{domain}_{year}{month:02}*.nc"
    Args:
        subfolder (str) :
            Main sub-folder containing NAME output files
        processed_folder (str) :
            Directory containing the processed file.
        param (dict) :
            Parameter dictionary containing
            'site','height','domain','year','month'
    Returns
        str / list:
            Filename
        None :
            If no files or more than one file are found
    '''
    directory = os.path.join(subfolder,processed_folder)
    search_str = "{site}-{height}_{domain}_{year}{month:02}*.nc".format(site=param["site"],height=param["height"],domain=param["domain"],year=param["year"],month=param["month"])
    full_search_str = os.path.join(directory,search_str)
    print("search string",full_search_str)
    processed_file = glob.glob(full_search_str)
    if len(processed_file) == 1:
        processed_file = processed_file[0]
    elif len(processed_file) == 0:
        print("No processed file found")
        processed_file = None
    elif len(processed_file) > 1:
        print("Too many processed files found. Not unique")
        processed_file = None
    
    return processed_file

def remove_processed_file(subfolder,processed_folder,param):
    '''
    Check and remove any files matching datestr extracted from param input from processed files folder.
    This is to allow test of process function to be run multiple times.
    Processed files are found within processed files folder based on datestr from 
    param["year"] + param["month"].zfill(2)
    '''
    
    files_in_folder = find_processed_file(subfolder,processed_folder,param)
    
#    folder = os.path.join(subfolder,processed_folder)
#    
#    search_str = "*{}*.nc".format(create_datestr(param))
#    search_str = os.path.join(folder,search_str)
#    files_in_folder = glob.glob(search_str)
    
    if files_in_folder is not None:
        os.remove(files_in_folder)
    
       
def test_process_satellite_byday(process_satellite_byday_param,
                                 subfolder_satellite_byday,folder_names,satellite_param):
    '''
    Test process function produces an output for satellite data grouped by day.
    '''
    
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_satellite_byday,processed_folder,satellite_param)
    
    out = process.process(**process_satellite_byday_param)

def test_process_satellite_bypoint(process_satellite_bypoint_param,
                                   subfolder_satellite_bypoint,folder_names,satellite_param):
    '''
    Test process function produces an output for satellite data separated by point.
    '''
    
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_satellite_bypoint,processed_folder,satellite_param)
    
    out = process.process(**process_satellite_bypoint_param)

def test_process_site(process_site_param,subfolder_site,folder_names,site_param):
    '''
    Test process function produces an output for site data.
    '''
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_site,processed_folder,site_param)
    
    out = process.process(**process_site_param)

def test_process_satellite_mf_byday(process_satellite_byday_mf_param,
                                 subfolder_satellite_byday_mf,folder_names,satellite_param):
    '''
    Test process function produces an output for satellite data grouped by day in mole fraction output.
    '''
    
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_satellite_byday_mf,processed_folder,satellite_param)
    
    out = process.process(**process_satellite_byday_mf_param)

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_process_site_against_org(process_site_param,subfolder_site,folder_names,site_param):
    '''
    Test output of original process function: process_org.process against new process script for 
    site data.
    '''
    
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_site,processed_folder,site_param)
    process_site_param_org = process_site_param.copy()
    process_site_param["use_surface_conditions"] = False
    
    out = process.process(**process_site_param)
    
    processed_folder_org = "Processed_Fields_files_Org"
    org_folder = os.path.join(subfolder_site,processed_folder_org)
    if not os.path.exists(org_folder):
        os.makedirs(org_folder)
    
    remove_processed_file(subfolder_site,processed_folder_org,site_param)
    
    process_site_param_org["processed_folder"] = processed_folder_org
    out_org = process_org.process(**process_site_param_org)
    
    print("out_org",out_org)
    
    data_vars = out.data_vars
    
    for dv in data_vars:
        assert np.array_equal(out[dv].values,out_org[dv].values)

@pytest.mark.skip(reason="No longer able to compare after update to use mole fraction output")
def test_process_satellite_bypoint_against_org(process_satellite_bypoint_param,
                                               process_satellite_org_bypoint_param,
                                               subfolder_satellite_bypoint,folder_names,satellite_param): 
    '''
    Test output of original process function: process_org.process against new process script for satellite
    data separated by point.
    '''
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_satellite_bypoint,processed_folder,satellite_param)

    out = process.process(**process_satellite_bypoint_param)

    processed_folder_org = "Processed_Fields_files_Org"    
    org_folder = os.path.join(subfolder_satellite_bypoint,processed_folder_org)
    if not os.path.exists(org_folder):
        os.makedirs(org_folder)
    
    remove_processed_file(subfolder_satellite_bypoint,processed_folder_org,satellite_param)
    
    process_satellite_org_bypoint_param["processed_folder"] = processed_folder_org
    out_org = process_org.process(**process_satellite_org_bypoint_param)
    
    data_vars = out.data_vars
    
    for dv in data_vars:
        assert np.array_equal(out[dv].values,out_org[dv].values)   

@pytest.mark.skip(reason="Possible comparison error due to new short-lived footprints. Requires updating.")
def test_process_satellite_byday_mf_against_bench(process_satellite_byday_mf_param,
                                                  subfolder_satellite_byday_mf,
                                                  folder_names,satellite_param): 
    '''
    Test mol fraction output against benchmarked file (created 2019-05-20). This marks
    a code update to check the units within the fields file and then use the correct unit 
    conversion. NAME output should now primarily be in ppm s units rather than gs/m3 but 
    this should ensure backwards compatibility with previous runs.
    '''
    processed_folder = folder_names["processed_folder"]
    remove_processed_file(subfolder_satellite_byday_mf,processed_folder,satellite_param)

    process.process(**process_satellite_byday_mf_param)
    output_file = find_processed_file(subfolder_satellite_byday_mf,
                                      processed_folder,
                                      satellite_param)
    out = open_ds(output_file)

    processed_folder_bench = "Processed_Fields_files_Benchmark"    
    bench_file = find_processed_file(subfolder_satellite_byday_mf,
                                     processed_folder_bench,
                                     satellite_param)
    out_bench = open_ds(bench_file) # Need to open comparison file and put in the same format

    data_vars = out.data_vars
    
    for dv in data_vars:
        assert np.array_equal(out[dv].values,out_bench[dv].values)   
