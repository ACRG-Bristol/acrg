#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 12:39:12 2020

@author: rt17603
"""

import os
import glob

import acrg_name.process as process

def define_subfolder(domain,site,height):
    ''' Define subfolder name based on domain, site and height '''
    subfolder = "{domain}_{site}_{height}/".format(domain=domain,site=site,height=height)
    return subfolder

def create_datestr(param):
    ''' Create basic date string based on year and month of the format "YYYYMM" '''
    return str(param["year"]) + str(param["month"]).zfill(2)

#%% Fields files / MixR

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

#%% Particle files

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

#%% Met files

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

#%% Footprint concatenate

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

def footprint_concatenate_param(subfolder,read_met,folder_names,parameters,satellite,
                                get_time_step=None):
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
        param["time_step"] = get_time_step
        param["upper_level"] = parameters["upper_level"]

    return param

#%% Process

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

