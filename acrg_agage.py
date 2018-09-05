# -*- coding: utf-8 -*-
"""
Get AGAGE data, average, truncate and filter for baseline

Expected format for filenames is:
	"network-instrument_site_date_species-height.nc"
	e.g. AGAGE-GC-FID_MHD_19940101_ch4-10m.nc

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
import os
from os import getenv
import re
from netCDF4 import Dataset
from acrg_time import convert
import json
import datetime as dt
import xarray as xray
import pdb

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

root_directory = os.path.join(data_path, "obs/")

#Get site info and species info from JSON files
with open(os.path.join(acrg_path, "acrg_species_info.json")) as f:
    species_info=json.load(f)

with open(os.path.join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f)

def open_ds(path):
    
    """
    Function efficiently opens xray datasets.
    
    Example:
        ds = open_ds("path/file.nc")
    """
    # use a context manager, to ensure the file gets closed after use
    with xray.open_dataset(path) as ds:
        ds.load()
    return ds    

def is_number(s):
    """
    Is it a number?
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def synonyms(search_string, info):
    """
    Todo:
        Make doc string
    """
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
    """Search over species info for variable.
    
    Args:
        varnames (list):
            variable names to search for
        species (str):
            Species of interest
        species_info (dict):
            species info from JSON value
        label (str):
            Entry label.
            Default='alt'
            
    Returns:
            out_string (str):
                Matches to search
                
    Note:
        I'm not actually sure if this is a correct description of that it's 
        doing. So if somebody actually has to use it then it's probably good 
        to check.
    """
    
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
    fnames = sorted([os.path.split(f)[1] for f in flag_files])
    return fnames, [f.split("_") for f in fnames]


def quadratic_sum(x):
    """
    Computes np.sqrt(np.sum(x**2))/float(len(x))
    
    Args: 
        x (array):
            
    """
    return np.sqrt(np.sum(x**2))/float(len(x))
    

#Get Met Office baseline flags
def ukmo_flags(site, site_info):
    """Function to get Met Office baseline flags
    
    Args:
       site (str):
           Name of site of interest as acrg_site_info.json.
       site_info (dict): 
          site info from JSON file.
          
    Returns:
        (pandas.Dataframe) Dataframe containing Met Office baseline flags for 
        site.
    """
    
    flag_directory=os.path.join(root_directory, "flags/")
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
            flag_data=pd.io.parsers.read_csv(os.path.join(flag_directory, f), 
                                             delim_whitespace=True, skiprows=6)
            flag_time = flag_time + [dt.datetime(y, m, d, h, mi)
                for y, m, d, h, mi in 
                zip(flag_data["YY"], flag_data["MM"], flag_data["DD"], 
                    flag_data["HH"], flag_data["MI"])]
            flag.append(flag_data["B"])
        
        flag=np.hstack(flag)
        
        return pd.DataFrame(flag, index=flag_time, columns=("flags",))


def get_file_list(site, species, start, end, height,
                  network = None, instrument = None, data_directory=None):
    
    """Gets list of relevent files for site, species, time and height.
    
    Args:
        site (str) :
            Obs site. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        start (str) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
        end (str) : 
            End date in same format as start for range of files to find.
        height (str/list) : 
            Height of inlet for input data (must match number of sites).
        network (str/list, optional) : 
            Network for the site/instrument (must match number of sites).
            Default = None.
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites).
            Default = None.
        data_directory (str, optional) :
            flux_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by network.
            Default=None.
            
    Returns:
        data_directory (str):
            If data_directory is specified on input then this is returned.
            Else returns relevent directory.
        files (list):
            List of relevent files.
            
    """
             
    if network is None:
        file_network_string = site_info[site]["network"]
    else:
        file_network_string = network
    
    if data_directory is None:
        data_directory=os.path.join(root_directory,file_network_string)
    else:
        data_directory=os.path.join(data_directory,file_network_string)
        
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
    fnames, file_info = file_search_and_split(os.path.join(data_directory,"*.nc"))

    if len(fnames) == 0:
        print("Can't find any data files: " + os.path.join(data_directory,"*.nc"))
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
        height=None, baseline=False, average=None, keep_missing=False,
        network = None, instrument = None,
        status_flag_unflagged = [0], data_directory=None):
    """
    Args:    
        site_in (str) :
            Site of interest. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species_in (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        start (str, optional) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
            Default = "1900-01-01".
        end (str, optional) : 
            End date in same format as start for range of files to find.
            Default="2020-01-01".
        height (str/list, optional) : 
            Height of inlet for input data (must match number of sites).
            Default=None
        baseline (bool, optional) : 
            *** Not actually used in this function, at present? ***
            Default=False.
        average (str/list, optional) :
            Averaging period for each dataset (for each site) ((must match number of sites)).
            Each value should be a string of the form e.g. "2H", "30min" (should match pandas offset 
            aliases format).
            Default=None.
        keep_missing (bool, optional) :
            Whether to keep missing data points or drop them.
            default=False.
        network (str/list, optional) : 
            Network for the site/instrument (must match number of sites).
            Default=None.
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites). 
            Default=None.
                status_flag_unflagged (list, optional) : 
            The value to use when filtering by status_flag. 
            Default = [0]
        data_directory (str, optional) :
            flux_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by network.
            Default=None.
            
    Returns:
        (xarray dataframe):
            Get timeseries data frame for observations at site of species.
        
        
        
    """
    
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)

    site = synonyms(site_in, site_info)
    
    if site is None:
        print("No site called {}. Either try a different name, or add name to site_info.json.".format(site_in))
        return

    species = synonyms(species_in, species_info)
    if species is None:
        print("No species called {}. Either try a different name, or add name to species_info.json.".format(species_in))
        return
        
    data_directory, files = get_file_list(site, species, start_time, end_time,
                                          height, network = network,
                                          instrument = instrument, data_directory=data_directory)                                 
    
    #Get files
    #####################################
    if files is not None:
    
        data_frames = []
        
        for f in files:
            
            print("Reading: " + f)
    
            skip = False
            
            ncf=Dataset(os.path.join(data_directory, f), 'r')
            
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
                elif ("hours" in ncf.variables["time"].units) is True:
                    time = convert.hours2time(ncf.variables["time"][:], 
                                            ncf.variables["time"].units[14:]) 
                elif ("days" in ncf.variables["time"].units) is True:
                    time = convert.day2time(ncf.variables["time"][:], 
                                            ncf.variables["time"].units[11:])
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
                
#                If ship read lat and lon data
#               Check if site info has a keyword called platform
                if 'platform' in site_info[site].keys():
                    if site_info[site]["platform"] == 'ship':
                        file_lat=ncf.variables["latitude"]
                        if len(file_lat) > 0:                                
                            df["meas_lat"] = file_lat[:]
                            
                        file_lon=ncf.variables["longitude"]
                        if len(file_lon) > 0:
                            df["meas_lon"] = file_lon[:]
                           
                    
                    #If platform is aircraft, get altitude data
                    if site_info[site]["platform"] == 'aircraft':
                        if "alt" in ncf.variables.keys():
                            file_alt=ncf.variables["alt"]
                            if len(file_alt) > 0:
                                df["altitude"] = file_alt[:]        
                        
                #Get status flag
                if ncvarname + "_status_flag" in ncf.variables.keys():
                    file_flag=ncf.variables[ncvarname + "_status_flag"]
                    if len(file_flag) > 0:
                        df["status_flag"] = file_flag[:]                  
                        # Flag out multiple flags
                        flag = [False for _ in range(len(df.index))]
                        for f in status_flag_unflagged:
                            flag = flag | (df.status_flag == f)
                        df = df[flag]
               
                        
                if units != "permil":
                    df = df[df.mf > 0.]

                if len(df) > 0:
                    data_frames.append(df)
    
            ncf.close()
    
        if len(data_frames) > 0:
            data_frame = pd.concat(data_frames).sort_index()
            data_frame.index.name = 'time'
            data_frame = data_frame[start_time : end_time]            
            if len(data_frame) == 0:
                return None
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


def get_gosat(site, species, max_level, start = "1900-01-01", end = "2020-01-01", data_directory = None):
    """retrieves obervations for a set of sites and species between start and 
    end dates for GOSAT 
    
    Args:    
        site (str) :
            Site of interest. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        max_level (int) : 
            Required for satellite data only. Maximum level to extract up to from within satellite data.
        start (str, optional) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
            Default = "1900-01-01".
        end (str, optional) : 
            End date in same format as start for range of files to find.
            Default="2020-01-01".
        data_directory (str, optional) :
            flux_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by network.
            Default=None.
            
    Returns:
        (xarray dataframe):
            xarray data frame for GOSAT observations.
            
    """
    if max_level is None:
        print "ERROR: MAX LEVEL REQUIRED FOR SATELLITE OBS DATA"
        return None
        
    start_time = convert.reftime(start)
    end_time = convert.reftime(end)

    data_directory, files = get_file_list(site, species, start_time, end_time,
                                          None, data_directory = data_directory)
    files_date = [convert.reftime(f.split("_")[2][0:8]) for f in files]
    files_keep = [True if d >= start_time and d < end_time else False for d in files_date]
    files = [f for (f, k) in zip(files, files_keep) if k == True]

    data = []
    for i,f in enumerate(files):
        #data.append(xray.open_dataset(os.path.join(data_directory,f)))
        # rt17603: 06/04/2018 Changed this to use function containing with statement. Was possible causing problems with too many files being open at once..
        data.append(open_ds(os.path.join(data_directory,f)))
        
    data = xray.concat(data, dim = "time")

    lower_levels =  range(0,max_level)

    prior_factor = (data.pressure_weights[dict(lev=list(lower_levels))]* \
                    (1.-data.xch4_averaging_kernel[dict(lev=list(lower_levels))])* \
                    data.ch4_profile_apriori[dict(lev=list(lower_levels))]).sum(dim = "lev")
                    
    upper_levels = range(max_level, len(data.lev.values))            
    prior_upper_level_factor = (data.pressure_weights[dict(lev=list(upper_levels))]* \
                    data.ch4_profile_apriori[dict(lev=list(upper_levels))]).sum(dim = "lev")
                
    data["mf_prior_factor"] = prior_factor
    data["mf_prior_upper_level_factor"] = prior_upper_level_factor
    data["mf"] = data.xch4 - data.mf_prior_factor - data.mf_prior_upper_level_factor
    data["dmf"] = data.xch4_uncertainty

    # rt17603: 06/04/2018 Added drop variables to ensure lev and id dimensions are also dropped, Causing problems in footprints_data_merge() function
    drop_data_vars = ["xch4","xch4_uncertainty","lon","lat","ch4_profile_apriori","xch4_averaging_kernel",
                      "pressure_levels","pressure_weights","exposure_id"]
    drop_coords = ["lev","id"]
    
    for dv in drop_data_vars:
        if dv in data.data_vars:
            data = data.drop(dv)
    for coord in drop_coords:
        if coord in data.coords:
            data = data.drop(coord)

    #data = data.drop("lev")
    #data = data.drop(["xch4", "xch4_uncertainty", "lon", "lat"])
    data = data.to_dataframe()
    # rt17603: 06/04/2018 Added sort because some data was not being read in time order. 
    # Causing problems in footprints_data_merge() function
    data = data.sort_index()
    
    data.max_level = max_level
    if species.upper() == "CH4":
        data.mf.units = 1e-9
    if species.upper() == "CO2":
        data.mf.units = 1e-6

    
    return data


def get_obs(sites, species, start, end,
            height = None, baseline = False, average = None, keep_missing=False,
            network = None, instrument = None, status_flag_unflagged = None,
            max_level = None, data_directory = None):
    """
    The get_obs function retrieves obervations for a set of sites and species between start and end dates
    Note: max_level only pertains to satellite data
    
    TODO: 
        At the moment satellite data is identified by finding "GOSAT" in the site name. This will be 
        changed to include a check by "platform" label from acrg_site_info.json to identify the type 
        of observation (e.g. satellite, ferry, aircraft)
    
    If height, network, instrument or average inputs are specified they should have 
    the same number of entries as there are sites. Format:
        - For one site, can include as a str or a one item list.
        - For multiple sites, must be a list matching the number of sites.
    The status_flag_unflagged must also match the number of sites but must always be a list.
 
    If not explicitly specified, height and network values can be extracted from acrg_site_info.json. If
    there are multiple values are present, the first one will be used by default 
    (*** May change with Matt's new acrg_site_info.json format ***)
    
    For example if we wanted to read in daily averaged methane data for Mace Head and Tacolneston for 2012 
    we could include:
    
        get_obs(sites=["MHD","TAC"],species="ch4",start="2012-01-01",end="2013-01-01",height=["10m","100m"],
                average=["24H","24H"],network=["AGAGE","DECC"])
    
    Args:
        sites (list) :
            Site list. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        start (str) : 
            Start date in format "YYYY-MM-DD" for range of files to find.
        end (str) : 
            End date in same format as start for range of files to find.
        height (str/list, optional) : 
            Height of inlet for input data (must match number of sites).
        baseline (bool, optional) : 
            *** Not actually used in this function, at present? ***
        average (str/list, optional) :
            Averaging period for each dataset (for each site) ((must match number of sites)).
            Each value should be a string of the form e.g. "2H", "30min" (should match pandas offset 
            aliases format). 
        keep_missing (bool, optional) :
            Whether to keep missing data points or drop them.
        network (str/list, optional) : 
            Network for the site/instrument (must match number of sites). (optional)
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites). (optional)
        status_flag_unflagged (list, optional) : 
            The value to use when filtering by status_flag. Default = [0]
        max_level (int) : 
            Required for satellite data only. Maximum level to extract up to from within satellite data.
        data_directory (str, optional) :
            flux_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by network.
    
    Returns:
        dict(pandas.DataFrame) : 
            pandas Dataframe for every site, keywords of ".species" and ".units" are also included in
            the dictionary.    
    """

    def check_list_and_length(var, sites, error_message_string):
        if var is not None:
            if type(var) is not list:
                var = [var]
            if len(var) != len(sites):
                print("If you're going to specify {}, make sure the length of ".format(error_message_string) +
                      "the height list is the same length as sites list. Returning.")
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
        #return None
        raise Exception("All of height, average, network and instrument parameters must match "+
                        "the dimensionality of sites parameter. "+
                        "height : {0}; average: {1}; network: {2}; instrument: {3}\n".format(height,average,network,instrument)+
                        "sites: {}, number of sites: {}".format(sites,len(sites)))

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
    
