# -*- coding: utf-8 -*-
"""
Get AGAGE data, average, truncate and filter for baseline

Examples:

Get Mace Head CH4 observations and average into 3 hourly periods:

    import acrg_agage as agage
    time, ch4, sigma = agage.get("MHD", "CH4", average="3H")
    
Get Mace Head CH4 observations, truncate for 2010-2012 inclusive
and average into 6 hourly periods:

    import acrg_agage as agage
    time, ch4, sigma = agage.get("MHD", "CH4", startY=2010, endY=2013,average="6H")

Calculate Cape Grim monthly means, with baseline filtering:
    
    import acrg_agage as agage
    time, hfc134a, sigma = 
        agage.get("CGO", "HFC-134a", baseline=True, average="1M")

Created on Sat Dec 27 17:17:01 2014
@author: chxmr
"""

import numpy as np
import pandas as pd
import glob
from os.path import split, join
from os import getenv
import re
from netCDF4 import Dataset
from acrg_time import convert
import json
import datetime as dt
import xarray as xr
from collections import OrderedDict
import os

acrg_path = getenv("ACRG_PATH")
data_path = getenv("DATA_PATH")

if acrg_path is None:
    acrg_path = getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")
if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

root_directory= join(data_path, "obs/")

#Get site info and species info from JSON files
with open(join(acrg_path, "acrg_species_info.json")) as f:
    species_info=json.load(f)

with open(join(acrg_path, "acrg_site_info_2018.json")) as f:
    site_info=json.load(f, object_pairs_hook=OrderedDict)

def open_ds(path):
    """
    Function efficiently opens xray datasets.
    """
    # use a context manager, to ensure the file gets closed after use
    with xray.open_dataset(path) as ds:
        ds.load()
    return ds    

#def is_number(s):
#    try:
#        float(s)
#        return True
#    except ValueError:
#        return False
##
def synonyms(search_string, info, alternative_label = "alt"):
    '''
    Check to see if there are other names that we should be using for
    a particular input. E.g. If CFC-11 or CFC11 was input, go on to use cfc-11,
    as this is used in species_info.json
    
    Args:
        search_string (str): Input string that you're trying to match
        info (dict): Dictionary whose keys are the "default" values, and an
             variable that contains other possible names
    Returns:
        corrected string
    '''
    
    keys=info.keys()
    
    #First test whether site matches keys (case insensitive)
    out_strings = \
        [k for k in keys if k.upper() == search_string.upper()]

    #If not found, search synonyms
    if len(out_strings) is 0:
        for k in keys:
            matched_strings = \
                [s for s in info[k][alternative_label] \
                    if s.upper() == search_string.upper()]
            if len(matched_strings) is not 0:
                out_strings = [k]
                break

    if len(out_strings) is 1:
        out_string = out_strings[0]
    else:
        out_string = None

    return out_string

def listsearch(possible_strings, correct_string,
               further_info, further_info_alternative_label = "alt"):
    '''
    Given a list of strings, check whether any of them match some correct string.
    They may differ by case, or there may be synonyms. 
    
    Args:
        possible_strings (list of strings): list of strings to search through
        correct_string (str): "correct" string that we are trying to match
        further_info (dict): dictionary containing "correct" keys,
            along with potential synonyms under the "further_info_alternative_label" key

    Returns:
        string that can be used as a search term to identify matches in
        "possible_strings"
    '''
    
    for v in possible_strings:
        if correct_string.upper() == v.upper():
            out_string = v
            break
        else:
            matched_strings=[s for s in further_info[correct_string][further_info_alternative_label] \
                if s.upper() == v.upper()]
            if len(matched_strings) is not 0:
                out_string = v
                break
            else:
                out_string = None

    return out_string


def file_search_and_split(search_string):
    '''
    Search for files using a particular search string, 
    then split wherever an underscore is found
    
    Args:
        search_string (str): A string of characters that will be searched
                             using glob
    
    Returns:
        fnames (list): list of matching filenames
        file info (list): for each filename a list of strings of text that was
                        separated by underscores
    '''

    files=glob.glob(search_string)
    fnames = sorted([split(f)[1] for f in files])
    return fnames, [f.split("_") for f in fnames]


def quadratic_sum(x):
    return np.sqrt(np.sum(x**2))/float(len(x))
    

def get_file_list(site, species,
                  inlet = None,
                  network = None,
                  instrument = None,
                  data_directory=None):
    '''
    Find all files that fit the site, species, start and end dates, inlet
    network and instrument.
    '''
    
    # Look for data in data directionry + site name folder
    if data_directory is None:
        data_directory = join(root_directory, site)
    else:
        data_directory = join(data_directory, site)

    # Test if file path exists
    if not os.path.exists(data_directory):
        print("Directory %s doesn't exist." %data_directory)
        return data_directory, None

    #Get file list and break down file name
    fnames, file_info = file_search_and_split(join(data_directory, "*.nc"))

    # Test if there are any files in the directory
    if len(fnames) == 0:
        print("Can't find any netCDF files in %s." % data_directory)
        return data_directory, None
    
    # Double check that sites match
    file_site = [f[1] for f in file_info]
    fnames = [f for (s, f) in zip(file_site, fnames) if s.upper() == site.upper()]
    file_info = [f for (s, f) in zip(file_site, file_info) if s.upper() == site.upper()]
    if len(fnames) == 0:
        print("Can't find any %s files in %s." %(site, data_directory) )
        return data_directory, None
    
    # Test if species is in folder
    file_species = [re.split("-|\.", f[-1])[0] for f in file_info]
    file_species_string = listsearch(file_species, species, species_info)
    if file_species_string is None:
        print("Can't find files for %s in %s" %(species, data_directory) )
        return data_directory, None
    # Subset of matching files
    fnames = [f for (sp, f) in zip(file_species, fnames) if sp == file_species_string]
    file_info = [f for (sp, f) in zip(file_species, file_info) if sp == file_species_string]

    # If instrument is specified, only return that instrument
    if instrument is not None:
        # Subset of matching files
        file_instrument = [re.split("-", f[0])[1] for f in file_info]
        fnames = [f for (inst, f) in zip(file_instrument, fnames) if inst == instrument]
        file_info = [f for (inst, f) in zip(file_instrument, file_info) if inst == instrument]
        if len(fnames) == 0:
            print("Cant find file matching instrument %s in %s" %(instrument, data_directory))
            return data_directory, None
 
    # See if any sites are separated by height
    file_inlet = [re.split("-|\.", f[-1])[1] for f in file_info]
    file_inlet = [None if fh == "nc" else fh for fh in file_inlet]
    # Are there any files in the folder for which no inlet height is defined?
    if any([f is None for f in file_inlet]):
        # If so, do none of the matching files have a height? If so, we need more info
        if not all([f is None for f in file_inlet]):
            print("It looks like there are some files in %s" %data_directory)
            print(" for which an inlet height has been defined, and others where none has been defined.")
            print("Make sure an inlet height is defined for all files, or change this function!")
            return data_directory, None
        # Else, we don't need to do anything            
    else:
        # If no inlet specified, pick first element in height list
        if inlet is None:
            file_inlet_string = site_info[site][network]["height"][0]
        else:
            file_inlet_string = inlet
        # Subset of matching files
        fnames = [f for (ht, f) in zip(file_inlet, fnames) if ht == file_inlet_string]
        file_info = [f for (ht, f) in zip(file_inlet, file_info) if ht == file_inlet_string]
        if len(fnames) == 0:
            print("Can't find any matching inlets: %s" %inlet)
            return data_directory, None

    # Return output
    if len(fnames) == 0:
        print("For some reason I got to this stage and there aren't any files here...")
        return data_directory, None
    else:
        return data_directory, sorted(fnames)



def get(site, species_in, network = None, start = None, end = None,
        inlet = None, average = None, keep_missing = False,
        instrument = None,
        status_flag_unflagged = [0], data_directory = None):
    '''
    Get measurements from one site
    
    '''
    
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)
    
    if site not in site_info.keys():
        print("No site called %s." % site)
        print("Either try a different name, or add name to acrg_site_info.json.")
        return

    species = synonyms(species_in, species_info)
    if species is None:
        print("No species called %s." % species_in)
        print("Either try a different name, or add name to species_info.json.")
        return
    
    # If no network specified, pick the first network in site_info dictionary
    if network is None:
        network = site_info[site].keys()[0]
        print("Assuming network is %s" %network)
    
    
    data_directory, files = get_file_list(site, species,
                                          network = network,
                                          inlet = inlet,
                                          instrument = instrument,
                                          data_directory = data_directory)

    # Iterate through files
    if files is not None:
        
        data_frames = []
        cal = []
        
        for f in files:
            
            print("Reading: " + f)
            
            # Read files using xarray
            with xr.open_dataset(join(data_directory, f)) as fxr:
                ds = fxr.load()
                
            # Record calibration scales
            if "Calibration_scale" in ds.attrs.keys():
                cal.append(ds.attrs["Calibration_scale"])            
                
                ncvarname=listsearch(ds.keys(), species, species_info)
                
                if ncvarname is None:
                    print("Can't find mole fraction variable name '" + species + 
                          "' or alternatives in file. Either change " + 
                          "variable name, or add alternative to " + 
                          "acrg_species_info.json")
                    return None

                # Create single site data frame
                df = pd.DataFrame({"mf": ds[ncvarname][:]},
                                  index = ds.time)
                
                if ds[ncvarname].units.isdigit():
                    units = float(ds[ncvarname].units)
                else:
                    units = str(ds[ncvarname].units)
                
                #Get repeatability
                if ncvarname + "_repeatability" in ds.keys():
                    file_dmf=ds[ncvarname + "_repeatability"].values
                    if len(file_dmf) > 0:
                        df["dmf"] = file_dmf[:]
        
                #Get variability
                if ncvarname + "_variability" in ds.keys():
                    file_vmf=ds[ncvarname + "_variability"]
                    if len(file_vmf) > 0:
                        df["vmf"] = file_vmf[:]
                
                # If ship read lat and lon data
                # Check if site info has a keyword called platform
                if 'platform' in site_info[site].keys():
                    if site_info[site][network]["platform"] == 'ship':
                        file_lat=ds["latitude"].values
                        if len(file_lat) > 0:                                
                            df["meas_lat"] = file_lat
                            
                        file_lon=ds["longitude"].values
                        if len(file_lon) > 0:
                            df["meas_lon"] = file_lon
                           
                    
                    #If platform is aircraft, get altitude data
                    if site_info[site]["platform"] == 'aircraft':
                        if "alt" in ds.keys():
                            file_alt=ds["alt"].values
                            if len(file_alt) > 0:
                                df["altitude"] = file_alt    
                        
                #Get status flag
                if ncvarname + "_status_flag" in ds.keys():
                    file_flag=ds[ncvarname + "_status_flag"].values
                    if len(file_flag) > 0:
                        df["status_flag"] = file_flag                
                        # Flag out multiple flags
                        flag = [False for _ in range(len(df.index))]
                        for f in status_flag_unflagged:
                            flag = flag | (df.status_flag == f)
                        df = df[flag]
               
                        
                if units != "permil":
                    df = df[df.mf > 0.]

                if len(df) > 0:
                    data_frames.append(df)
    
        # If multiple files, merge
        if len(data_frames) > 0:

            # Check if all calibration scales are the same
            if not all([c == cal[0] for c in cal]):
                print("ERROR: Can't merge the following files, as calibration scales are different:")
                for f in files:
                    print(f)
                return None
            
            data_frame = pd.concat(data_frames).sort_index()
            data_frame.index.name = 'time'
            data_frame = data_frame[start_time : end_time]
            
            if len(data_frame) == 0:
                return None

        else:

            return None


        # If averaging is set, resample
        if average is not None:
            how = {}
            for key in data_frame.columns:
                if key == "dmf":
                    how[key] = quadratic_sum
                elif key == "vmf":
                    # Calculate std of 1 min mf obs in av period as new vmf 
                    how[key] = "std"
                    data_frame["vmf"] = data_frame["mf"]
                else:
                    how[key] = "mean"
            
            if keep_missing == True:
                if min(data_frame.index) > start_time:
                    dum_frame = pd.DataFrame({"status_flag": float('nan')},
                                         index = np.array([start_time]))   
                    dum_frame["mf"] =  float('nan')
                    dum_frame["dmf"] =  float('nan')
                    dum_frame.index.name = 'time'                                                                               
                    data_frame = data_frame.append(dum_frame)
            
                if max(data_frame.index) < end_time:
                    dum_frame2 = pd.DataFrame({"status_flag": float('nan')},
                                         index = np.array([end_time]))   
                    dum_frame2["mf"] =  float('nan')
                    dum_frame2["dmf"] =  float('nan')
                    dum_frame2.index.name = 'time'
                    data_frame = data_frame.append(dum_frame2)
            
            data_frame=data_frame.resample(average, how=how)
            if keep_missing == True:
                data_frame=data_frame.drop(data_frame.index[-1])
              
        # Drop NaNs
        if keep_missing == False:  
            data_frame.dropna(inplace = True)
        
        data_frame.mf.units = units
        data_frame.files = files
    
        return data_frame

    else:
        
        return None


#def get_gosat(site, species, max_level, start = "1900-01-01", end = "2020-01-01", data_directory = None):
#    
#    if max_level is None:
#        print "ERROR: MAX LEVEL REQUIRED FOR SATELLITE OBS DATA"
#        return None
#        
#    start_time = convert.reftime(start)
#    end_time = convert.reftime(end)
#
#    data_directory, files = get_file_list(site, species, start_time, end_time,
#                                          None, data_directory = data_directory)
#    files_date = [convert.reftime(f.split("_")[2][0:8]) for f in files]
#    files_keep = [True if d >= start_time and d < end_time else False for d in files_date]
#    files = [f for (f, k) in zip(files, files_keep) if k == True]
#
#    data = []
#    for i,f in enumerate(files):
#        #data.append(xray.open_dataset(join(data_directory,f)))
#        # rt17603: 06/04/2018 Changed this to use function containing with statement. Was possible causing problems with too many files being open at once..
#        data.append(open_ds(join(data_directory,f)))
#        
#    data = xray.concat(data, dim = "time")
#
#    lower_levels =  range(0,max_level)
#
#    prior_factor = (data.pressure_weights[dict(lev=list(lower_levels))]* \
#                    (1.-data.xch4_averaging_kernel[dict(lev=list(lower_levels))])* \
#                    data.ch4_profile_apriori[dict(lev=list(lower_levels))]).sum(dim = "lev")
#                    
#    upper_levels = range(max_level, len(data.lev.values))            
#    prior_upper_level_factor = (data.pressure_weights[dict(lev=list(upper_levels))]* \
#                    data.ch4_profile_apriori[dict(lev=list(upper_levels))]).sum(dim = "lev")
#                
#    data["mf_prior_factor"] = prior_factor
#    data["mf_prior_upper_level_factor"] = prior_upper_level_factor
#    data["mf"] = data.xch4 - data.mf_prior_factor - data.mf_prior_upper_level_factor
#    data["dmf"] = data.xch4_uncertainty
#
#    # rt17603: 06/04/2018 Added drop variables to ensure lev and id dimensions are also dropped, Causing problems in footprints_data_merge() function
#    drop_data_vars = ["xch4","xch4_uncertainty","lon","lat","ch4_profile_apriori","xch4_averaging_kernel",
#                      "pressure_levels","pressure_weights","exposure_id"]
#    drop_coords = ["lev","id"]
#    
#    for dv in drop_data_vars:
#        if dv in data.data_vars:
#            data = data.drop(dv)
#    for coord in drop_coords:
#        if coord in data.coords:
#            data = data.drop(coord)
#
#    #data = data.drop("lev")
#    #data = data.drop(["xch4", "xch4_uncertainty", "lon", "lat"])
#    data = data.to_dataframe()
#    # rt17603: 06/04/2018 Added sort because some data was not being read in time order. Causing problems in footprints_data_merge() function
#    data = data.sort_index()
#    
#    data.max_level = max_level
#    if species.upper() == "CH4":
#        data.mf.units = 1e-9
#    if species.upper() == "CO2":
#        data.mf.units = 1e-6
#
#    
#    return data
#

def get_obs(sites, species, start = "1900-01-01", end = "2020-01-01",
            height = None, baseline = False, average = None, keep_missing=False,
            network = None, instrument = None, status_flag_unflagged = None,
            max_level = None, data_directory = None):

    # retrieves obervations for a set of sites and species between start and end dates
    # max_level only pertains to satellite data

    def check_list_and_length(var, sites, error_message_string):
        if var is not None:
            if type(var) is not list:
                var = [var]
            if len(var) != len(var):
                print("If you're going to specify " + error_message_string +
                      ", make sure the length of the height list is " +
                      "the same length as sites list. Returning.")
                return None
        else:
            var = [None for i in sites]
        return var

    if status_flag_unflagged:
        if len(status_flag_unflagged) != len(sites):
            print("Status_flag_unflagged must be a LIST OF LISTS with the same dimension as sites")
            return None
    else:
        status_flag_unflagged = [[0] for _ in sites]
        print("Assuming status flag = 0 for all sites")

    # Prepare inputs
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)

    if end_time < start_time:
        print("End time is before start time")
        return None

    if type(sites) is not list:
        sites = [sites]

    height = check_list_and_length(height, sites, "height")
    average = check_list_and_length(average, sites, "average")
    network = check_list_and_length(network, sites, "network")
    instrument = check_list_and_length(instrument, sites, "instrument")
    
    if height == None or average == None or network == None or instrument == None:
        return None

    # Get data
    obs = {}
    for si, site in enumerate(sites):

        if "GOSAT" in site.upper():
            data = get_gosat(site, species,
                       start = start_time, end = end_time, max_level = max_level,
                       data_directory = data_directory)
            if data is None:
                return
                
        else:
            data = get(site, species, height = height[si],
                       start = start_time, end = end_time,
                       average = average[si],
                       network = network[si],
                       instrument = instrument[si],
                       keep_missing = keep_missing,
                       status_flag_unflagged = status_flag_unflagged[si],
                       data_directory = data_directory)                     
                       
        if data is not None:
            obs[site] = data.copy()            
            if "GOSAT" in site.upper():
                obs[site].max_level = data.max_level
                
        # Add some attributes
        if data is not None:
            obs[".species"] = species
            obs[".units"] = data.mf.units

    
    if len(obs) == 0:
        return None
    else:
        return obs
    