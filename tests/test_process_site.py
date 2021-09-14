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
import pytest
import os
import glob
import numpy as np
import datetime
import pandas as pd
import xarray as xr

import acrg.name.process as process
from acrg.config.paths import Paths

acrg_path = Paths.acrg

#%%  Setting up directory fixtures

@pytest.fixture()
def site_directory():
    ''' Define base directory containing satellite footprint files with separated data points '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/Site/")
    return directory

@pytest.fixture()
def benchmark_output_directory():
    ''' Define base directory containing grouped satellite footprint files '''
    directory = os.path.join(acrg_path,"tests/files/LPDM/raw_output/benchmarks/")
    return directory


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

@pytest.fixture()
def site_param_xunbounded():
    '''
    Define input parameters for site data to test.
    Keys include:
        "domain" (str),"site" (str),"height" (str),"year" (int),"month" (int),"satellite" (bool)
    '''
    param = {}
    param["domain"] = "ARCTIC"
    param["site"] = "CBR"
    param["height"] = "10magl"
    param["year"] = 2010
    param["month"] = 1
    param["satellite"] = False
   
    return param

@pytest.fixture()
def site_param_species():
    '''
    Define input parameters for site data to test.
    Keys include:
        "domain" (str),"site" (str),"height" (str),"year" (int),"month" (int),"satellite" (bool), "species" (str)
    '''
    param = {}
    param["domain"] = "CARIBBEAN"
    param["site"] = "RPB"
    param["height"] = "10magl"
    param["year"] = 2010
    param["month"] = 1
    param["satellite"] = False
    param["species"] = "HFO-1234zee"
   
    return param

@pytest.fixture()
def site_param_co2():
    '''
    Define input parameters for site data to test.
    Keys include:
        "domain" (str),"site" (str),"height" (str),"year" (int),"month" (int),"satellite" (bool), "species" (str)
    '''
    param = {}
    param["domain"] = "CARIBBEAN"
    param["site"] = "RPB"
    param["height"] = "10magl"
    param["year"] = 2010
    param["month"] = 1
    param["satellite"] = False
    param["species"] = "co2"
   
    return param

#%%

def define_subfolder(domain,site,height):
    ''' Define subfolder name based on domain, site and height '''
    subfolder = "{domain}_{site}_{height}/".format(domain=domain,site=site,height=height)
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
def subfolder_xunbounded_site(site_param_xunbounded,site_directory):
    ''' Define full subfolder path for site data. '''
    domain = site_param_xunbounded["domain"]
    site = site_param_xunbounded["site"]
    height = site_param_xunbounded["height"]
    
    subfolder_site = define_subfolder(domain,site,height)
    subfolder = os.path.join(site_directory,subfolder_site)  
    return subfolder

@pytest.fixture()
def folder_names():
    '''
    Define subfolders within NAME output directory structure
    Keys includes:
        "fields_folder", "fields_folder_hourly", "particles_folder","met_folder","obs_folder", 
        (also "mixr_fields_folder","mixr_hourly_folder" but largely deprecated)
    '''
    param = {}
    
    param["fields_folder"] = "MixR_files"
    param["fields_folder_hourly"] = "MixR_hourly"
    param["particles_folder"] = "Particle_files"
    param["met_folder"] = "Met_daily"
    param["obs_folder"] = "Observations"
    
    param["mixr_fields_folder"] = param["fields_folder"] # supporting legacy code (for now)
    param["mixr_hourly_folder"] = param["fields_folder_hourly"]  # supporting legacy code (for now)

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
def get_mixr_files_xunbounded_site(subfolder_xunbounded_site,folder_names,site_param_xunbounded):
    ''' 
    Get filenames of fields files for satellite run with separate points. 
    Finds field files based on datestr (see create_datestr), one file per day.
    '''
    datestr = create_datestr(site_param_xunbounded)
    fields_folder = folder_names["mixr_fields_folder"]
    
    fields_files = get_fields_files(subfolder_xunbounded_site,fields_folder,datestr)
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
def mixr_file_xunbounded_information():
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
             'X grid origin': -98.07578,
             'Y grid origin': 35.97656,
             'X grid size': 1024,
             'Y grid size': 230,
             'X grid resolution': 0.3515625,
             'Y grid resolution': 0.2343750,
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
    data_arrays.append(np.zeros((230,1024)))
    data_arrays[0][131,0] = 1
    data_arrays[0][135,1] = 2
    data_arrays.append(np.zeros((230,1024)))
    data_arrays[1][131,0] = 3
    data_arrays[1][135,1] = 4
    
    namever = 2
    
    lons = np.arange(-98.07578, -98.07578+0.3515625*1024-0.01, 0.3515625) + 0.3515625/2
    lats = np.arange(35.97656 , 35.97656 +0.2343750*230-0.01, 0.2343750) + 0.2343750/2
    
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
def particle_file_xunbounded_benchmark():
    '''Benchmark infromation from the dummy particle file, replacing values on E and W directions with 0 to test this feature'''
    pl_n = np.zeros((2,1,1024,20))
    pl_w = np.zeros((2,1,230,20))
    pl_e = np.zeros((2,1,230,20))
    mean_age_n = np.zeros((2,1,1024,20))
    mean_age_w = np.zeros((2,1,230,20))
    mean_age_e = np.zeros((2,1,230,20))
    
    pl_n[0,0,2,0] = 0.5   
    pl_n[0,0,3,0] = 0.5   
    pl_w[0,0,9,0] = 0
    pl_e[0,0,14,0] = 0
    
    mean_age_n[0,0,2,0] = 120
    mean_age_n[0,0,3,0] = 100
    mean_age_w[0,0,9,0] = 0
    mean_age_e[0,0,14,0] = 0

    
    return pl_n, pl_w, pl_e, mean_age_n, mean_age_w, mean_age_e

@pytest.fixture()
def met_file_benchmark():
    '''
    Defining dummy values for the satellite met data. Should be consistent with
    the dummy satellite met data which we can check against
    
    TODO: May want to read this in in some way rather than just define values
    '''
    values = {}
    
    #values["release_lon"] = [-59.43330,-59.43330]
    #values["release_lat"] = [13.16670,13.16670]
    #values["temp"] = [27.07520,26.98370]
    #values["PBLH"] = [778.9576,782.4219]
    #values["press"] = [101161.4,101163.7]
    #values["wind"] = [5.697351,5.695338]
    #values["wind_direction"] = [72.90572,76.33878]

    # Define values for each parameter
    # Repeated for hours 0-11 and 12-23
    release_lon_H_0_11 = [-59.43330]*12 # Repeat the same value 12 times for each hour
    release_lat_H_0_11 = [13.16670]*12
    temp_H_0_11 = [27.07520]*12
    PBLH_H_0_11 = [778.9576]*12
    press_H_0_11 = [101161.4]*12
    wind_H_0_11 = [5.697351]*12
    wind_direction_H_0_11 = [72.90572]*12

    release_lon_H_12_23 = [-59.43330]*12 # Repeat the same value 12 times for each hour
    release_lat_H_12_23 = [13.16670]*12
    temp_H_12_23 = [26.98370]*12
    PBLH_H_12_23 = [782.4219]*12
    press_H_12_23 = [101163.7]*12
    wind_H_12_23 = [5.695338]*12
    wind_direction_H_12_23 = [76.33878]*12

    # Concatenate two sets of values together to make a table
    values["release_lon"] = release_lon_H_0_11 + release_lon_H_12_23
    values["release_lat"] = release_lat_H_0_11 + release_lat_H_12_23
    values["temp"] = temp_H_0_11 + temp_H_12_23
    values["PBLH"] = PBLH_H_0_11 + PBLH_H_12_23
    values["press"] = press_H_0_11 + press_H_12_23
    values["wind"] = wind_H_0_11 + wind_H_12_23
    values["wind_direction"] = wind_direction_H_0_11 + wind_direction_H_12_23

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
    
 
#%%  Test particle_locations()
  

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
def get_particle_files_xunbounded_site(subfolder_xunbounded_site,folder_names,
                                         get_mixr_files_xunbounded_site):
    ''' 
    Get filenames of particle files for satellite run separated by point.
    Note: filenames found are based on datestr extracted from the input field files.
    '''
    particles_folder = folder_names["particles_folder"]
    particle_files = get_particle_files(subfolder_xunbounded_site,particles_folder,
                                        get_mixr_files_xunbounded_site)
    return particle_files

@pytest.fixture()
def define_heights():
    ''' Define height range for input into particle_locations function. '''
    dheights = 1000
    heights = np.arange(0, 19001, dheights) + dheights/2.
    return heights

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

def test_particle_locations_xunbounded_site(particle_file_xunbounded_benchmark,mixr_file_xunbounded_information,get_particle_files_xunbounded_site,define_heights):
    '''
    Test particle_locations() function can produce the correct output for a. dummy particle location file.
    '''

    header_bench, column_headings_bench, data_arrays_bench, namever_bench, lons_bench, lats_bench, levs_bench, time_bench, timeStep_bench = mixr_file_xunbounded_information

    particle_file_1 = get_particle_files_xunbounded_site[0]
    heights =  define_heights
    
    pl_n_bench, pl_w_bench, pl_e_bench, mean_age_n_bench, mean_age_w_bench, mean_age_e_bench = particle_file_xunbounded_benchmark

    hist = process.particle_locations(particle_file_1,time_bench,lats_bench,lons_bench,levs_bench, heights,id_is_lev=False,
                                 satellite=False)

    assert np.array_equal(hist["pl_n"].values,pl_n_bench)
    assert np.array_equal(hist["pl_w"].values,pl_w_bench)
    assert np.array_equal(hist["pl_e"].values,pl_e_bench)
    assert np.array_equal(hist["mean_age_n"].values,mean_age_n_bench)
    assert np.array_equal(hist["mean_age_w"].values,mean_age_w_bench)        
    assert np.array_equal(hist["mean_age_e"].values,mean_age_e_bench) 


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
def get_met_files_site(subfolder_site,folder_names):
    ''' Get filenames of met files for site run. '''
    
    met_folder = folder_names["met_folder"]
    met_files = get_met_files(subfolder_site,met_folder,satellite=False)
    
    return met_files

def test_read_met_site(get_met_files_site, met_file_benchmark):
    '''
    Test read_met function can read in site met data successfully.
    '''
    met_files_1 = get_met_files_site[0]
    
    out = process.read_met(met_files_1,vertical_profile=False, satellite=False)
    
    expected = met_file_benchmark
    compare_columns = expected.columns
    
    for col in compare_columns:
        assert np.array_equal(out[col],expected[col])

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

def test_footprint_array_site(benchmark_output_directory,get_mixr_files_site,
                              get_particle_files_site,
                              read_met_site):
    '''
    Test footprint_array function can return an output for site files.
    '''
    
    fp_bench_file = os.path.join(benchmark_output_directory,'Test_process_footprint_array.nc')
    
    with xr.open_dataset(fp_bench_file) as f: 
        fp_bench = f.load()
    
    fields_file_1 = get_mixr_files_site[0]
    particle_file_1 = get_particle_files_site[0]
    
#     met_files_1 = get_met_files_site[0]
    
#     met = process.read_met(met_files_1,vertical_profile=False, satellite=False)

    met = read_met_site
    
    out = process.footprint_array(fields_file_1, particle_file_1, met)
    
    for key in fp_bench.keys():
        assert np.array_equal(fp_bench[key], out[key])

       
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

    # footprint_concatenate expects str rather than list input for "datestr" and "met"
    if isinstance(param["datestr"], list):
        print("datestr", param["datestr"])
        param["datestr"] = param["datestr"][0]
    
    if isinstance(param["met"], list):
        print("met", param["met"])
        param["met"] = param["met"][0]

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

#@pytest.mark.skip(reason="Needs to be updated.")
def test_footprint_concatenate_site(fc_param_site):
    '''
    Test footprint_concatenate function produces an output when files are for site data.
    '''
    param = fc_param_site

    out = process.footprint_concatenate(**param)
    assert out is not None
    
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
    
    if "species" in input_param:
        param["species"] = input_param["species"]
        param["fields_folder"] = folder_names["fields_folder_hourly"]
    else:
        param["fields_folder"] = folder_names["fields_folder"]
    
    param["particles_folder"] = folder_names["particles_folder"]
    param["met_folder"] = folder_names["met_folder"]
    #param["processed_folder"] = folder_names["processed_folder"]

    param["satellite"] = satellite
    
    return param
    
@pytest.fixture()
def process_site_param(site_param,folder_names,site_directory):
    '''
    Define input parameters for process.process function for site data.
    Additional "base_dir" parameter added to point to Site folder.
    '''    
    param = process_param(site_param,folder_names,satellite=False)
    param["base_dir"] = site_directory
    
    process_dir = os.path.join(site_directory, "Processed_Fields_files")
    param["process_dir"] = process_dir
    
    return param

def find_processed_file(param):
    '''
    Find NAME processed file based on search string:
        "{site}-{height}[_{species}]_{domain}_{year}{month:02}*.nc"
    Args:
        #subfolder (str) :
        #    Main sub-folder containing NAME output files
        #processed_folder (str) :
        #    Directory containing the processed file.
        param (dict) :
            Parameter dictionary containing
            'site','height','domain','year','month', 'species', 'process_dir'
    Returns
        str / list:
            Filename
        None :
            If no files or more than one file is found
    '''
    #directory = os.path.join(subfolder,processed_folder)
    directory = os.path.join(param["process_dir"], param["domain"])

    if "species" in param:
        optional_species_str = "_"+param["species"].lower()
    else:
        optional_species_str = ""

    search_str = f"{param['site']}-{param['height']}{optional_species_str}_{param['domain']}_{param['year']}{param['month']:02}*.nc"

    full_search_str = os.path.join(directory,search_str)
    
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

#def remove_processed_file(subfolder,processed_folder,param):
def remove_processed_file(param):
    '''
    Check and remove any files matching datestr extracted from param input from processed files folder.
    This is to allow test of process function to be run multiple times.
    Processed files are found within processed files folder based on datestr from 
    param["year"] + param["month"].zfill(2)
    '''
    
    files_in_folder = find_processed_file(param)
    
#    folder = os.path.join(subfolder,processed_folder)
#    
#    search_str = "*{}*.nc".format(create_datestr(param))
#    search_str = os.path.join(folder,search_str)
#    files_in_folder = glob.glob(search_str)
    
    if files_in_folder is not None:
        os.remove(files_in_folder)
    
# Rudimentary check that this runs (for now)
def test_process_site(process_site_param):
    '''
    Test process function produces an output for site data.
    '''
    remove_processed_file(process_site_param)
    out = process.process(**process_site_param)

    assert out is not None

def test_process_site_benchmark(process_site_param):
    '''
    Test process function output against a benchmark for the same data.
    '''
    remove_processed_file(process_site_param)
    process.process(**process_site_param)

    out_file = find_processed_file(process_site_param)
    out = xr.open_dataset(out_file)

    # Benchmark was run on 14/09/2021
    # RPB, 10magl, CARIBBEAN domain, January, 2010 data
    # *If this does not match to your current inputs - benchmark must be updated.*

    benchmark_file = "files/LPDM/raw_output/benchmarks/RPB-10magl_CARIBBEAN_201001.nc"
    benchmark = xr.open_dataset(benchmark_file)

    print(out["particle_locations_n"], benchmark["particle_locations_n"])

    xr.testing.assert_equal(out, benchmark)


#%%

@pytest.fixture()
def process_site_species_param(site_param_species,folder_names,site_directory):
    '''
    Define input parameters for process.process function for site data.
    Additional "base_dir" parameter added to point to Site folder.
    '''    
    param = process_param(site_param_species, folder_names, satellite=False)
    
    param["base_dir"] = site_directory
    
    process_dir = os.path.join(site_directory, "Processed_Fields_files")
    param["process_dir"] = process_dir

    #param["unzip"] = True
    
    return param

# Rudimentary check that this runs (for now)
def test_process_site_species(process_site_species_param):
    '''
    Test process function can work with hourly input netcdf data.

    Additional input parameters:
     - "species" = "HFO-1234zee"
     - "fields_folder" = "MixR_hourly"

    This is run against a paired down netcdf hourly file which contains a 
    4-hour back period (OutputTime) for 1 release (ReleaseTime) at 04:00 
    '''
    remove_processed_file(process_site_species_param)
    
    print(process_site_species_param)
    out = process.process(**process_site_species_param)

    assert out is not None
    assert out is not False


def test_process_site_species_benchmark(process_site_species_param):
    '''
    Test process function output against a benchmark for the same data.
    '''
    remove_processed_file(process_site_species_param)
    
    process.process(**process_site_species_param)

    out_file = find_processed_file(process_site_species_param)
    out = xr.open_dataset(out_file)

    # Benchmark was run on 14/09/2021
    # RPB, 10magl, CARIBBEAN domain, January, 2010 data for species hfo-1234zee
    # *If this does not match to your current inputs - benchmark must be updated.*

    benchmark_file = "files/LPDM/raw_output/benchmarks/RPB-10magl_hfo-1234zee_CARIBBEAN_201001.nc"
    benchmark = xr.open_dataset(benchmark_file)

    xr.testing.assert_equal(out, benchmark)


@pytest.fixture()
def process_site_co2_param(site_param_co2,folder_names,site_directory):
    '''
    Define input parameters for process.process function for site data.
    Additional "base_dir" parameter added to point to Site folder.
    '''    
    param = process_param(site_param_co2,folder_names,satellite=False)
    
    param["base_dir"] = site_directory

    process_dir = os.path.join(site_directory, "Processed_Fields_files")
    param["process_dir"] = process_dir

    return param


# Rudimentary check that this runs (for now)
def test_process_site_co2(process_site_co2_param):
    '''
    Test process function work with hourly input netcdf data when species is co2.

    Additional input parameters:
     - "species" = "co2"
     - "fields_folder" = "MixR_hourly"

    This is run against a paired down netcdf hourly file which contains a 
    4-hour back period (OutputTime) for 1 release (ReleaseTime) at 04:00 
    '''
    remove_processed_file(process_site_co2_param)
    
    out = process.process(**process_site_co2_param)
    
    assert out is not None
    assert out is not False

def test_process_site_co2_benchmark(process_site_co2_param):
    '''
    Test process function output against a benchmark for the same data.
    '''
    remove_processed_file(process_site_co2_param)
    process.process(**process_site_co2_param)

    out_file = find_processed_file(process_site_co2_param)
    out = xr.open_dataset(out_file)

    # Benchmark was run on 14/09/2021
    # RPB, 10magl, CARIBBEAN domain, January, 2010 for co2 data
    # *If this does not match to your current inputs - benchmark must be updated.*

    benchmark_file = "files/LPDM/raw_output/benchmarks/RPB-10magl_co2_CARIBBEAN_201001.nc"
    benchmark = xr.open_dataset(benchmark_file)

    xr.testing.assert_equal(out, benchmark)

