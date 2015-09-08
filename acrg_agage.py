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
import xray

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

with open(join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def synonyms(search_string, info):
    
    keys=info.keys()
    
    #First test whether site matches keys (case insensitive)
    out_strings = \
        [k for k in keys if k.upper() == search_string.upper()]

    #If not found, search synonyms
    if len(out_strings) is 0:
        for k in keys:
            matched_strings = \
                [s for s in info[k]["alt"] \
                    if s.upper() == search_string.upper()]
            if len(matched_strings) is not 0:
                out_strings = [k]
                break

    if len(out_strings) is 1:
        out_string = out_strings[0]
    else:
        out_string = None

    return out_string

def listsearch(varnames, species, species_info, label="alt"):

    for v in varnames:
        if species.upper() == v.upper():
            out_string = v
            break
        else:
            matched_strings=[s for s in species_info[species][label] \
                if s.upper() == v.upper()]
            if len(matched_strings) is not 0:
                out_string = v
                break
            else:
                out_string = None

    return out_string


def file_search_and_split(search_string):
    flag_files=glob.glob(search_string)
    fnames = sorted([split(f)[1] for f in flag_files])
    return fnames, [f.split("_") for f in fnames]


def quadratic_sum(x):
    return np.sqrt(np.sum(x**2))/float(len(x))


#Get Met Office baseline flags
def ukmo_flags(site, site_info):
    
    flag_directory=join(root_directory, "flags/")
    fnames, file_info=file_search_and_split(
        flag_directory + "*.txt")
    file_site = [f[0] for f in file_info]
    file_site_string=listsearch(file_site, site, site_info)

    if file_site_string is None:

        return None
        
    else:

        files=sorted([f for f in fnames if file_site_string in f])
    
        flag=[]
        flag_time=[]
        
        for f in files:
            flag_data=pd.io.parsers.read_csv(join(flag_directory, f), 
                                             delim_whitespace=True, skiprows=6)
            flag_time = flag_time + [dt.datetime(y, m, d, h, mi)
                for y, m, d, h, mi in 
                zip(flag_data["YY"], flag_data["MM"], flag_data["DD"], 
                    flag_data["HH"], flag_data["MI"])]
            flag.append(flag_data["B"])
        
        flag=np.hstack(flag)
        
        return pd.DataFrame(flag, index=flag_time, columns=("flags",))


def get_file_list(site, species, start, end, height,
                  network = None, instrument = None):
        
    if network is None:
        file_network_string = site_info[site]["network"]
    else:
        file_network_string = network

    data_directory=join(root_directory, file_network_string)
    
    if height is None:
        file_height_string = site_info[site]["height"][0]
    else:
        if type(height) is not list:
            height = [height]

        file_height_string = listsearch(height, site, site_info, 
                                        label="height")
        if file_height_string is None:
            print("Height " + height + " doesn't exist in site_info.json. "
                + "Available heights are " + str(site_info[site]["height"])
                + ". Leave blank for default height.")
            return data_directory, None
    
    #Get file info
    fnames, file_info = file_search_and_split(join(data_directory, "*.nc"))

    if len(fnames) == 0:
        print("Can't find any data files: " + join(data_directory, "*.nc"))
        return data_directory, None
        
    file_site = [f[1] for f in file_info]
    file_species = [re.split("-|\.", f[-1])[0] for f in file_info]
    file_height = [re.split("-|\.", f[-1])[1] for f in file_info]
    
    #Get file list
    file_species_string = listsearch(file_species, species, species_info)
    if file_species_string is None:
        print("Can't find files for species " + species + " in " + data_directory)
        return data_directory, None
        
    file_site_string = listsearch(file_site, site, site_info)
    if file_site_string is None:
        print("Can't find files for site " + site + " in " + data_directory)
        return data_directory, None

    files = [f for f, si, sp, hi in \
             zip(fnames, file_site, file_species, file_height) if \
                 si.upper() == file_site_string.upper() and
                 sp.upper() == file_species_string.upper() and
                 hi.upper() == file_height_string.upper()]

    if instrument is not None:
        files = [f for f in files if instrument.upper() in f.upper()]

    if len(files) == 0:
        return data_directory, None
    else:
        files.sort()
        return data_directory, files


def get(site_in, species_in, start = "1900-01-01", end = "2020-01-01",
        height=None, baseline=False, average=None,
        network = None, instrument = None, status_flag_unflagged = 0):
    
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)
    
    site = synonyms(site_in, site_info)
    if site is None:
        print("No site called " + site_in +
            ". Either try a different name, or add name to site_info.json.")
        return

    species = synonyms(species_in, species_info)
    if species is None:
        print("No species called " + species_in +
            ". Either try a different name, or add name to species_info.json.")
        return
        
    data_directory, files = get_file_list(site, species, start_time, end_time,
                                          height, network = network,
                                          instrument = instrument)

    #Get files
    #####################################
    
    if files is not None:
    
        data_frames = []
        
        for f in files:
            
            print("Reading: " + f)
    
            skip = False
            
            ncf=Dataset(join(data_directory, f), 'r')
    
            if "time" not in ncf.variables:
                print("Skipping: " + f + ". No time variable")
                skip = True
    
            else:
                if ("seconds" in ncf.variables["time"].units) is True:
                    time = convert.sec2time(ncf.variables["time"][:], 
                                            ncf.variables["time"].units[14:])
                elif ("minutes" in ncf.variables["time"].units) is True:
                    time = convert.min2time(ncf.variables["time"][:], 
                                            ncf.variables["time"].units[14:]) 
                elif ("days" in ncf.variables["time"].units) is True:
                    time = convert.day2time(ncf.variables["time"][:], 
                                            ncf.variables["time"].units[14:])
                else: 
                    print("Time unit is not a recognized unit (seconds, minuties or days since")
                            
                if max(time) < start_time:
                    skip = True
                if min(time) > end_time:
                    skip = True
    
            if not skip:
    
                ncvarname=listsearch(ncf.variables.keys(), species, species_info)
                
                if ncvarname is None:
                    print("Can't find mole fraction variable name '" + species + 
                          "' or alternatives in file. Either change " + 
                          "variable name, or add alternative to " + 
                          "acrg_species_info.json")
                    ncf.close()
                    return None

                df = pd.DataFrame({"mf": ncf.variables[ncvarname][:]},
                                  index = time)
                
                if is_number(ncf.variables[ncvarname].units):
                    units = float(ncf.variables[ncvarname].units)
                else:
                    units = str(ncf.variables[ncvarname].units)
                
                #Get repeatability
                if ncvarname + "_repeatability" in ncf.variables.keys():
                    file_dmf=ncf.variables[ncvarname + "_repeatability"]
                    if len(file_dmf) > 0:
                        df["dmf"] = file_dmf[:]
        
                #Get variability
                if ncvarname + "_variability" in ncf.variables.keys():
                    file_vmf=ncf.variables[ncvarname + "_variability"]
                    if len(file_vmf) > 0:
                        df["vmf"] = file_vmf[:]
        
                #Get status flag
                if ncvarname + "_status_flag" in ncf.variables.keys():
                    file_flag=ncf.variables[ncvarname + "_status_flag"]
                    if len(file_flag) > 0:
                        df["status_flag"] = file_flag[:]
                        df = df[df.status_flag == status_flag_unflagged]

                if units != "permil":
                    df = df[df.mf > 0.]

                data_frames.append(df)
    
            ncf.close()
    
        if len(data_frames) > 0:
            data_frame = pd.concat(data_frames).sort()
            data_frame.index.name = 'time'
            data_frame = data_frame[start_time : end_time]
        else:
            return None

        #Do baseline filtering
        if baseline:
            #Get flags
            flagdf=ukmo_flags(site, site_info)
            
            if flagdf is None:
                print("No baseline flags, sorry!")
            else:
                #Truncate time series if no flags available
                data_frame=data_frame[min(flagdf.index) : max(flagdf.index)]
                #Re-index flag time series to mf Data frame index
                data_frame["flag"]=flagdf.reindex(index=data_frame.index,
                                                  method='pad')
                #Only retain background
                data_frame=data_frame[data_frame['flag']==10]
        
        if average is not None:
            how = {}
            for key in data_frame.columns:
                if key == "dmf":
                    how[key] = quadratic_sum
                else:
                    how[key] = "median"
            data_frame=data_frame.resample(average, how=how)
    
        # Drop NaNs
        data_frame.dropna(inplace = True)
    
        data_frame.mf.units = units
    
        return data_frame

    else:
        
        return None


def get_gosat(site, species, start = "1900-01-01", end = "2020-01-01"):
    
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)

    data_directory, files = get_file_list(site, species, start_time, end_time,
                                          None)
    files_date = [convert.reftime(f.split("_")[2][0:8]) for f in files]
    files_keep = [True if d >= start_time and d < end_time else False for d in files_date]
    files = [f for (f, k) in zip(files, files_keep) if k == True]

    data = []
    for f in files:
        data.append(xray.open_dataset(join(data_directory,f)))
        
    data = xray.concat(data, dim = "time")

    prior_factor = (data.pressure_weights* \
                    (1.-data.xch4_averaging_kernel)* \
                    data.ch4_profile_apriori).sum(dim = "lev")

    data["mf_prior_factor"] = prior_factor
    data["mf"] = data.xch4 - data.mf_prior_factor
    data["dmf"] = data.xch4_uncertainty

    data = data.drop("lev")
    data = data.drop(["xch4", "xch4_uncertainty", "lon", "lat"])
    data = data.to_dataframe()
    
    if species.upper() == "CH4":
        data.mf.units = 1e-9
    if species.upper() == "CO2":
        data.mf.units = 1e-6

    return data


def get_obs(sites, species, start = "1900-01-01", end = "2020-01-01",
            height = None, baseline = False, average = None,
            network = None, instrument = None, status_flag_unflagged = 0):


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
        print(site)
        if "GOSAT" in site.upper():
            data = get_gosat(site, species,
                       start = start_time, end = end_time)
        else:
            data = get(site, species, height = height[si],
                       start = start_time, end = end_time,
                       average = average[si],
                       network = network[si],
                       instrument = instrument[si],
                       status_flag_unflagged = status_flag_unflagged)
                       
        if data is not None:
            obs[site] = data.copy()
    
    # Add some attributes
    if data is not None:
        obs[".species"] = species
        obs[".units"] = data.mf.units
    
    if len(obs) == 0:
        return None
    else:
        return obs
    