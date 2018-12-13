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
import json
import xarray as xr
from collections import OrderedDict

acrg_path = os.path.dirname(os.path.realpath(__file__))
data_path = getenv("DATA_PATH")

if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

# Set default obs folder
obs_directory = os.path.join(data_path, "obs_2018/")

#Get site info and species info from JSON files
with open(os.path.join(acrg_path, "acrg_species_info.json")) as f:
    species_info=json.load(f)

with open(os.path.join(acrg_path, "acrg_site_info_2018.json")) as f:
    site_info=json.load(f, object_pairs_hook=OrderedDict)

def open_ds(path):
    """
    Function efficiently opens xray datasets.
    
    Example:
        ds = open_ds("path/file.nc")
    """
    # use a context manager, to ensure the file gets closed after use
    with xr.open_dataset(path) as ds:
        ds.load()
    return ds    


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
    fnames = sorted([os.path.split(f)[1] for f in files])

    return fnames, [f.split("_") for f in fnames]


def quadratic_sum(x):
    """
    Computes np.sqrt(np.sum(x**2))/float(len(x))
    
    Args: 
        x (array):
            
    """
    return np.sqrt(np.sum(x**2))/len(x)
    

def file_list(site, species,
              inlet = None,
              network = None,
              instrument = None,
              version = None,
              data_directory=None):
    '''
    Find all files that fit the site, species, start and end dates, inlet
    network and instrument.
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
    '''
    
    # Look for data in data directionry + site name folder
    if data_directory is None:
        data_directory = os.path.join(obs_directory, site)
    else:
        data_directory = os.path.join(data_directory, site)

    # Test if file path exists
    if not os.path.exists(data_directory):
        print("Directory %s doesn't exist." %data_directory)
        return data_directory, None

    #Get file list and break down file name
    fnames, file_info = file_search_and_split(os.path.join(data_directory,
                                                           "*.nc"))

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
    
    # If network is specified, only return that network
    if network is not None:
        # Subset of matching files
        file_network = [re.split("-", f[0])[0] for f in file_info]
        fnames = [f for (nwork, f) in zip(file_network, fnames) if nwork == network]
        file_info = [f for (nwork, f) in zip(file_network, file_info) if nwork == network]
        if len(fnames) == 0:
            print("Cant find file matching network %s in %s" %(network, data_directory))
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
    # Test to see whether there are 3 hyphen-delimited entries in file suffix
    # If only 2, there isn't an inlet
    # If there is an inlet, it'll be in the second position
    # This will need to change if we ever change that suffix
    file_suffix_length = [len(re.split("-", f[-1])) for f in file_info]
    
    # Check if there are a mix of inlet and non-inlet files (can't handle this)
    if len(set(file_suffix_length)) > 1:
        print("ERROR: Can't have a mix of files with and without inlets defined in %s." %data_directory)
        return data_directory, None

    # Are there any files in the folder for which no inlet height is defined?
    # If so, don't need to do anything. Otherwise, choose inlet
    if file_suffix_length[0] > 2:
        
        # If all inlets are selected, do nothing
        if inlet != "all":
        
            file_inlet = [re.split("-|\.", f[-1])[1] for f in file_info]
            
            # If no inlet specified, pick first element in parameters height list
            if inlet is None:
                file_inlet_string = site_info[site][network]["height"][0]
            else:
                file_inlet_string = inlet
            # Subset of matching files
            fnames = [f for (ht, f) in zip(file_inlet, fnames) if ht == file_inlet_string]
            file_info = [f for (ht, f) in zip(file_inlet, file_info) if ht == file_inlet_string]
            
            if len(fnames) == 0:
                print("Can't find any matching inlets: %s, or use 'all'" %inlet)
                return data_directory, None


    # specify version number
    file_version = [re.split("-|\.", f[-1])[-2] for f in file_info]
    
    # If specific version is requested, use that, otherwise, take the max version
    if version is None:
        file_version_string = max(file_version)
    else:
        file_version_string = version
        
    # Subset of matching files
    fnames = [f for (v, f) in zip(file_version, fnames) if v == file_version_string]
    file_info = [f for (v, f) in zip(file_version, file_info) if v == file_version_string]
    
    if len(fnames) == 0:
        print("Can't find any matching versions: %s" %version)
        return data_directory, None
    else:
        # Return file list
        return data_directory, sorted(fnames)



def get_single_site(site, species_in,
                    network = None,
                    start_date = None, end_date = None,
                    inlet = None, average = None,
                    keep_missing = False,
                    instrument = None,
                    status_flag_unflagged = [0],
                    version = None,
                    data_directory = None):
    '''
    Get measurements from one site
    
    Args:    
        site_in (str) :
            Site of interest. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species_in (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        start_date (str, optional) : 
            Output start date in a format that Pandas can interpret
            Default = None.
        end_date (str, optional) : 
            Output end date in a format that Pandas can interpret
            Default=None.
        inlet (str/list, optional) : 
            Inlet label. If you want to merge all inlets, use "all"
            Default=None
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
            directpry can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by site.
            Default=None.
            
    Returns:
        (Pandas dataframe):
            Timeseries data frame for observations at site of species.

    '''
    
    # Check that site is in acrg_site_info.json
    if site not in site_info.keys():
        print("No site called %s." % site)
        print("Either try a different name, or add name to acrg_site_info.json.")
        return

    species = synonyms(species_in, species_info)
    if species is None:
        print("No species called %s." % species_in)
        print("Either try a different name, or add name to species_info.json.")
        return
    if species != species_in:
        print("... changing species from %s to %s" % (species_in, species))
    
    # If no network specified, pick the first network in site_info dictionary
    if network is None:
        network = site_info[site].keys()[0]
        print("... assuming network is %s" %network)
    elif network not in site_info[site].keys():
        print("Error: Available networks for site %s are %s." % (site, site_info[site].keys()))
        print("       You'll need to add network %s to acrg_site_info.json" %network)
        return
        
    data_directory, files = file_list(site, species,
                                      network = network,
                                      inlet = inlet,
                                      instrument = instrument,
                                      version = version,
                                      data_directory = data_directory)

    # Iterate through files
    if files is not None:
        
        data_frames = []
        cal = []
        
        for f in files:
            
            print("... reading: " + f)
            
            # Read files using xarray
            with xr.open_dataset(os.path.join(data_directory, f)) as fxr:
                ds = fxr.load()
                
            # Record calibration scales
            if "Calibration_scale" in ds.attrs.keys():
                cal.append(ds.attrs["Calibration_scale"])            
                
                ncvarname = listsearch(ds.keys(), species, species_info)
                
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
                if ncvarname + " repeatability" in ds.keys():
                    file_dmf=ds[ncvarname + " repeatability"].values
                    if len(file_dmf) > 0:
                        df["dmf"] = file_dmf[:]
        
                #Get variability
                if ncvarname + " variability" in ds.keys():
                    file_vmf=ds[ncvarname + " variability"]
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
                if ncvarname + " status_flag" in ds.keys():
                    file_flag=ds[ncvarname + " status_flag"].values
                    if len(file_flag) > 0:
                        df["status_flag"] = file_flag                
                        # Flag out multiple flags
                        flag = [False for _ in range(len(df.index))]
                        for f in status_flag_unflagged:
                            flag = flag | (df.status_flag == f)
                        df = df[flag]
                               
#                if units != "permil" and units != "per meg":
#                    df = df[df.mf > 0.]

                if len(df) > 0:
                    data_frames.append(df)
    
        if len(data_frames) == 0:
            warningMessage = '''For some reason, there is no valid data 
                                within files %s 
                                (e.g. all flagged)''' % ",".join(files)
            print(warningMessage)
            return None
            
        # Check if all calibration scales are the same
        if not all([c == cal[0] for c in cal]):
            errorMessage = '''Can't merge the following files,
                              as calibration scales 
                              are different: ''' % ",".join(files)
            raise ValueError(errorMessage)

        # Merge files        
        data_frame = pd.concat(data_frames).sort_index()
        data_frame.index.name = 'time'

        # Cut out requested time period
        data_frame = data_frame[start_date : end_date]
        if len(data_frame) == 0:
            print("Warning, no data within range")
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
                # Pad with an empty entry at the start date
                if min(data_frame.index) > pd.to_datetime(start_date):
                    data_frame.loc[pd.to_datetime(start_date)] = \
                        [np.nan for col in data_frame.columns]           
                # Pad with an empty entry at the end date
                if max(data_frame.index) < pd.to_datetime(end_date):
                    data_frame.loc[pd.to_datetime(end_date)] = \
                        [np.nan for col in data_frame.columns]
                # Now sort to get everything in the right order
                data_frame.sort_index(inplace = True)
            
            # Resample
            data_frame=data_frame.resample(average, how = how)
            if keep_missing == True:
                data_frame=data_frame.drop(data_frame.index[-1])
             
        # Drop NaNs
        if keep_missing == False:  
            data_frame.dropna(inplace = True)
        
        data_frame.mf.units = units
        data_frame.mf.scale = cal[0]
        data_frame.files = files

        return data_frame

    else:

        return None


def get_gosat(site, species, max_level,
              start_date = None, end_date = None,
              data_directory = None):
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
        start_date (str, optional) : 
            Start date in Pandas-readable format. E.g. "2000-01-01 00:00"
            Default = None.
        end_date (str, optional) : 
            End date in Pandas-readable format.
            Default = None.
        data_directory (str, optional) :
            directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by site.
            Default=None.
            
    Returns:
        (xarray dataframe):
            xarray data frame for GOSAT observations.
            
    """
    if max_level is None:
        raise ValueError("MAX LEVEL REQUIRED FOR SATELLITE OBS DATA")
            
    data_directory = os.path.join(obs_directory, "GOSAT", site)
    files = glob.glob(os.path.join(data_directory, '*.nc'))
    files = [os.path.split(f)[-1] for f in files]

    files_date = [pd.to_datetime(f.split("_")[2][0:8]) for f in files]

    data = []
    for (f, d) in zip(files, files_date):
        if d >= pd.to_datetime(start_date) and d < pd.to_datetime(end_date):
            with xr.open_dataset(os.path.join(data_directory,f)) as fxr:
                data.append(fxr.load())
    
    if len(data) == 0:
        return None
    
    data = xr.concat(data, dim = "time")

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
        data.mf.units = '1e-9'
    if species.upper() == "CO2":
        data.mf.units = '1e-6'

    data.mf.scale = "GOSAT"

    return data
    
def get_obs(sites, species,
            start_date = None, end_date = None,
            inlet = None,
            average = None,
            keep_missing=False,
            network = None,
            instrument = None,
            version = None,
            status_flag_unflagged = None,
            max_level = None,
            data_directory = None):
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
    
        get_obs(sites=["MHD","TAC"],species="ch4",start_date="2012-01-01",end_date="2013-01-01",height=["10m","100m"],
                average=["24H","24H"],network=["AGAGE","DECC"])
    
    Args:
        sites (list) :
            Site list. All sites should be defined within acrg_site_info.json. 
            E.g. ["MHD"] for Mace Head site.
        species (str) :
            Species identifier. All species names should be defined within acrg_species_info.json. 
            E.g. "ch4" for methane.
        start_date (str) : 
            Output start date in a format that Pandas can interpret
        end_date (str) : 
            Output end date in a format that Pandas can interpret
        inlet (str/list, optional) : 
            Height of inlet for input data (must match number of sites).
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
        if var != None:
            if type(var) is not list:
                raise ValueError("Variable '%s' must be a list" % error_message_string)
            if len(var) != len(sites):
                raise ValueError("Length of variable '%s' is different to 'sites'" % error_message_string)
        else:
            var = [None for i in sites]
        return var

    if status_flag_unflagged:
        if len(status_flag_unflagged) != len(sites):
            errorMessage = '''Status_flag_unflagged must be a 
                              LIST OF LISTS with the same 
                              first dimension as sites'''
            raise ValueError(errorMessage)
    else:
        status_flag_unflagged = [[0] for _ in sites]
        print("Assuming status flag = 0 for all sites")

    if type(sites) is not list:
        raise ValueError("Sites variable must be a list")

    inlet = check_list_and_length(inlet, sites, "inlet")
    average = check_list_and_length(average, sites, "average")
    network = check_list_and_length(network, sites, "network")
    instrument = check_list_and_length(instrument, sites, "instrument")
    version = check_list_and_length(version, sites, "version")
    
    # Get data
    obs = {}
    units = []
    scales = []
    
    for si, site in enumerate(sites):
        print("Getting %s data for %s..." %(species, site))
        if "GOSAT" in site.upper():
            data = get_gosat(site, species,
                       start_date = start_date, end_date = end_date,
                       max_level = max_level,
                       data_directory = data_directory)
        else:
            data = get_single_site(site, species, inlet = inlet[si],
                                   start_date = start_date, end_date = end_date,
                                   average = average[si],
                                   network = network[si],
                                   instrument = instrument[si],
                                   version = version[si],
                                   keep_missing = keep_missing,
                                   status_flag_unflagged = status_flag_unflagged[si],
                                   data_directory = data_directory)
        
        if data is None:
            obs[site] = None
            units.append(None)
            scales.append(None)
        else:    
            obs[site] = data.copy()
            if "GOSAT" in site.upper():
                if data is not None:
                    obs[site].max_level = data.max_level
            units.append(data.mf.units)
            scales.append(data.mf.scale)
    
    # Raise error if units don't match
    if len(set([u for u in units if u != None])) > 1:
        print(set(units))
        siteUnits = [':'.join([site, u]) for (site, u) in zip(sites, str(units)) if u is not None]
        errorMessage = '''Units don't match for these sites: %s''' % ', '.join(siteUnits)        
        raise ValueError(errorMessage)

    # Warning if scales don't match
    if len(set([s for s in scales if s != None])) > 1:
        siteScales = [':'.join([site, scale]) for (site, scale) in zip(sites, scales) if scale is not None]
        warningMessage = '''WARNING: scales don't match for these sites:
                            %s''' % ', '.join(siteScales)
        print(warningMessage)

    # Add some attributes
    obs[".species"] = species
    obs[".units"] = units[0]    # since all units are the same
    obs[".scales"] = {si:s for (si, s) in zip(sites, scales)}
    
    return obs
    


def plot(data_dict):
    '''
    Plot the data dictionary created by get_obs
    
    Args:
        data_dict (dict) :
            Dictionary of Pandas DataFrames output by get_obs. 
            Keys are site names.
        
    '''
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    plots = []
    
    for site, df in data_dict.iteritems():
        if site[0] != '.':
            if df is None:
                plots.append(mpatches.Patch(label="%s (no data)" %site, color = "white"))
            else:
                if "vmf" in df.columns:
                    plots.append(plt.errorbar(df.index, df.mf, df.vmf,
                                 linewidth = 0,
                                 marker = '.', markersize = 3.,
                                 label = site))
                if "dmf" in df.columns:
                    plots.append(plt.errorbar(df.index, df.mf, df.dmf,
                                 linewidth = 0, 
                                 marker = '.', markersize = 3.,
                                 label = site))

    plt.legend(handles = plots)
    plt.ylabel("%s (%s)" %(data_dict[".species"], data_dict[".units"]))
    
    plt.show()



