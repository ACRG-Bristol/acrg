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
from __future__ import print_function
from __future__ import division

from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import pandas as pd
import glob
import os
import re
import json
import xarray as xr
from collections import OrderedDict
import sys
import sqlite3

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH")
    data_path = os.getenv("DATA_PATH")
    obs_directory = os.path.join(data_path,"obs")    
else:
    from acrg_config.paths import paths

    acrg_path = paths.acrg
    obs_directory = paths.obs

#Get site info and species info from JSON files
#with open(acrg_path / "acrg_species_info.json") as f:
with open(os.path.join(acrg_path,"acrg_species_info.json")) as f:
    species_info=json.load(f)

with open(os.path.join(acrg_path,"acrg_site_info.json")) as f:
    site_info=json.load(f, object_pairs_hook=OrderedDict)

def is_number(s):
    """
    Is it a number?
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

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
    
    keys=list(info.keys())
    
    #First test whether site matches keys (case insensitive)
    out_strings = \
        [k for k in keys if k.upper() == search_string.upper()]

    #If not found, search synonyms
    if len(out_strings) == 0:
        for k in keys:
            matched_strings = \
                [s for s in info[k][alternative_label] \
                    if s.upper() == search_string.upper()]
            if len(matched_strings) != 0:
                out_strings = [k]
                break

    if len(out_strings) == 1:
        out_string = out_strings[0]
    else:
        out_string = None

    return out_string


def scale_convert(df, species, to_scale):
    '''
    Convert to a new calibration scale, based on conversions in acrg_obs_scale_convert.csv
    
    Args:
        df (Pandas Dataframe): Must contain an mf variable (mole fraction), with a .scale attribute
        species (str): species name
        to_scale (str): Calibration scale to convert to
    
    '''
    
    # If scale is already correct, return
    df_scale = df.scale
    if df_scale == to_scale:
        return(df)
    else:
        print(f"... converting scale to {to_scale}")
    
    scale_converter = pd.read_csv("acrg_obs_scale_convert.csv")
    scale_converter_scales = scale_converter[scale_converter.isin([species, df_scale, to_scale])][["species", "scale1", "scale2"]].dropna(axis=0, how = "any")
    
    if len(scale_converter_scales) == 0:
        errorMessage = f'''Scales {df_scale} and {to_scale} are not both in any one row in acrg_obs_scale_convert.csv for species {species}'''
        raise ValueError(errorMessage)
    elif len(scale_converter_scales) > 1:
        errorMessage = f'''Duplicate rows in acrg_obs_scale_convert.csv?'''
        raise ValueError(errorMessage)
    else:
        row = scale_converter_scales.index[0]
    
    converter = scale_converter.loc[row]

    if to_scale == converter["scale1"]:
        direction = "2to1"
    else:
        direction = "1to2"

    df.mf = eval(converter[direction].replace("X", "df.mf"))
    df.scale = to_scale
    
    return(df)


def get_single_site(site, species_in,
                    network = None,
                    start_date = "1900-01-01", end_date = "2100-01-01",
                    inlet = None, average = None,
                    keep_missing = False,
                    instrument = None,
                    status_flag_unflagged = [0],
                    version = None,
                    data_directory = None,
                    calibration_scale = None,
                    verbose = False):
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
    if site not in list(site_info.keys()):
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
    
    
    
    # Open defaults file
    df_defaults = pd.read_csv(paths.acrg / "acrg_obs/acrg_obs_defaults.csv",
                             parse_dates = ["startDate", "endDate"])
    df_defaults.dropna(axis = 0, how = "all", inplace = True)

    # Set early start date and late end date, if empty
    df_defaults["startDate"] = df_defaults["startDate"].fillna(pd.Timestamp("1900-01-01"))
    df_defaults["endDate"] = df_defaults["endDate"].fillna(pd.Timestamp("2100-01-01"))

    df_defaults.replace(np.nan, "%", inplace = True)

    
    
    # Create database in memory
    conn = sqlite3.connect(":memory:")

    dtypes = {"site": "text",
              "species": "text",
              "startDate": "timestamp",
              "endDate": "timestamp",
              "instrument": "text",
              "inlet": "text"}

    # Convert to database
    df_defaults.to_sql("defaults", conn, if_exists="replace", dtype = dtypes)

    c = conn.cursor()

    # Attach the obs database
    c.execute("ATTACH DATABASE ? as obs", (str(paths.obs / "obs.db"),))

    # Change some variables for query
    start_date_query = pd.Timestamp(start_date).to_pydatetime()
    end_date_query = pd.Timestamp(end_date).to_pydatetime()
    species_query = species.replace("-", "").lower()
    
    
    # Run a couple of initial queries to see whether there are any defaults defined for a particular site
    df_defaults_for_site = pd.read_sql_query("SELECT * FROM defaults WHERE site=?",
                                             conn, params = (site,))
    df_defaults_for_site_species = pd.read_sql_query("SELECT * FROM defaults WHERE site=? AND species=?", 
                                                     conn, params = (site, species))

    # Check if defaults need to be over-written
    if (inlet != None or instrument != None or network != None) and len(df_defaults_for_site) >= 1:
        print("... You've set either an inlet, instrument or network, overriding any defaults.")
        print("... Best to set all three of these, if you want to avoid ambiguity.")
        override_defaults = True
    else:
        override_defaults = False

        
        
    # Read filenames from database
    if len(df_defaults_for_site) == 0 or override_defaults:
        # Query only the 'files' table table to determine which files to read

        print("... no defaults set")
        query = '''
                SELECT files.filename
                FROM obs.files
                WHERE files.species=? COLLATE NOCASE AND
                      files.site=? COLLATE NOCASE AND
                      (files.endDate > date(?) AND files.startDate < date(?))
                '''
        params = [species_query, site, start_date_query, end_date_query]

        # If an inlet, network or instrument are specified, append to the query
        if inlet != None:
            query += " AND files.inlet=?"
            params.append(inlet)

        if network != None:
            query += " AND files.network=?"
            params.append(network)

        if instrument != None:
            query += " AND files.instrument=?"
            params.append(instrument)

    else:
        # Create an inner join of the defaults and files table and filter based on defaults
        # Note that the LIKE statement below allows us to have a wildcard (%) in some columns
        # The wildcard is set wherever there is a missing value in the defaults table, and
        #  where there's only one inlet in the files database

        query = '''
                SELECT filename, defaultStartDate, defaultEndDate FROM
                (
                 SELECT files.filename,
                        files.startDate,
                        files.endDate,
                        files.site,
                        files.inlet,
                        files.species,
                        files.instrument,
                        files.network,
                        defaults.inlet,
                        defaults.site,
                        defaults.instrument,
                        defaults.network,
                        defaults.startDate AS defaultStartDate,
                        defaults.endDate AS defaultEndDate
                FROM obs.files
                INNER JOIN defaults
                ON files.site = defaults.site
                WHERE files.species = ? COLLATE NOCASE AND
                        files.species LIKE defaults.species COLLATE NOCASE AND
                        files.site = ? COLLATE NOCASE AND
                        files.inlet LIKE defaults.inlet AND
                        files.instrument LIKE defaults.instrument AND
                        files.network LIKE defaults.network COLLATE NOCASE AND
                        (files.endDate > date(?) AND files.startDate < date(?)) AND
                        (defaults.endDate > date(?) AND defaults.startDate < date(?))
                        )
                '''


        # If species explicitly appears in the defaults file, enforce matching to that value
        # This is needed in the case that a species is measured on two instruments, 
        #  but one instrument measures a whole load of stuff and therefore a wildcard is set 
        #  for that instrument (e.g. SF6 on the ECD and Medusa at TAC)
        if len(df_defaults_for_site_species) > 0:
            query = query.replace("files.species LIKE defaults.species",
                                  "files.species = defaults.species")

        params = [species_query, site,
                  start_date_query, end_date_query, start_date_query, end_date_query]

    if verbose:
        print(query)

    # Run query and get list of files
    files_to_get = c.execute(query, params)
    
    obs_files = []

    # Retrieve files
    for f in files_to_get:

        print(f"... reading {f[0]}")
        with xr.open_dataset(f[0]) as f_ds:
            ds = f_ds.load()

        # If 3 elements, it means the query returned a start and end date
        # Otherwise, return whole dataset
        if len(f) == 3:
            if pd.Timestamp(start_date) > pd.Timestamp(f[1]):
                slice_start = pd.Timestamp(start_date)
            else:
                slice_start = pd.Timestamp(f[1])
            if pd.Timestamp(end_date) < pd.Timestamp(f[2]):
                slice_end = pd.Timestamp(end_date) - pd.Timedelta("1 ns")
            else:
                slice_end = pd.Timestamp(f[2])

        else:

            slice_start = pd.Timestamp(start_date)
            slice_end = pd.Timestamp(end_date)

        if slice_start.round(freq = "T") != pd.Timestamp("1900-01-01") or  slice_end.round(freq = "T") != pd.Timestamp("2100-01-01"):
            print(f"... slicing from {slice_start} to {slice_end}")
            ds = ds.loc[dict(time = slice(slice_start, slice_end))]

        # If averaging is set, resample
        if average != None:

            # First, just do a mean resample on all variables
            print(f"... resampling to {average}")
            ds_resampled = ds.resample(time = average, keep_attrs = True
                                       ).mean()
            # keep_attrs doesn't seem to work for some reason, so manually copy
            ds_resampled.attrs = ds.attrs.copy()

            
            # For some variables, need a different type of resampling
            # For all variables, copy attributes
            for var in ds.variables:
                if "repeatability" in var:
                    ds_resampled[var] = np.sqrt((ds[var]**2).resample(time = average).sum()) / \
                                                 ds[var].resample(time = average).count()

                elif "variability" in var:
                    # Calculate std of 1 min mf obs in av period as new vmf 
                    ds_resampled[var] = ds[var].resample(time = average,
                                                         keep_attrs = True).std()

                # Copy over some attributes
                if "long_name" in ds[var].attrs:
                    ds_resampled[var].attrs["long_name"] = ds[var].attrs["long_name"]
                if "units" in ds[var].attrs:
                    ds_resampled[var].attrs["units"] = ds[var].attrs["units"]

                ds = ds_resampled.copy()

        obs_files.append(ds)

    conn.close()
    
    return(obs_files)


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
        raise ValueError("'max_level' ARGUMENT REQUIRED FOR SATELLITE OBS DATA")
            
    data_directory = os.path.join(data_directory, "GOSAT", site)
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

    lower_levels =  list(range(0,max_level))

    prior_factor = (data.pressure_weights[dict(lev=list(lower_levels))]* \
                    (1.-data.xch4_averaging_kernel[dict(lev=list(lower_levels))])* \
                    data.ch4_profile_apriori[dict(lev=list(lower_levels))]).sum(dim = "lev")
                    
    upper_levels = list(range(max_level, len(data.lev.values)))            
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

    data = data.to_dataframe()
    data = MetaDataFrame(data)
    # rt17603: 06/04/2018 Added sort because some data was not being read in time order. 
    # Causing problems in footprints_data_merge() function
    data = data.sort_index()
    
    data.max_level = max_level
    if species.upper() == "CH4":
        data.units = 1e-9
    if species.upper() == "CO2":
        data.units = 1e-6

    data.scale = "GOSAT"

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
            data_directory = None,
            calibration_scale = None):
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
            If you want to merge all inlets, use "all"
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

    if type(average) is not list:
        average = [average]*len(sites)

    if (data_directory is None):
        data_directory = obs_directory

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
                                   data_directory = data_directory,
                                   calibration_scale = calibration_scale)
        
        if data is None:
            obs[site] = None
            units.append(None)
            scales.append(None)
        else:
            # reset mutli index into the expected standard index
            obs[site] = data.reset_index().set_index("time")
            if "GOSAT" in site.upper():
                if data is not None:
                    obs[site].max_level = data.max_level
            units.append(data.units)
            scales.append(data.scale)
    
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
    
def label_species(species):
    '''
    Write species label in correct format for plotting with subscripts for number of atoms.
    e.g. "CH4" becomes "CH$_{4}$" or "CHCl3" becomes "CHCl$_{3}$"
    '''
    import re
    num = re.findall(r"\d+",species)
    species_str = species
    if num:
        for n in num:
            species_str = species_str.replace(n,"$_{{{num}}}$".format(num=n)) # {{{}}} needed to both use string formatting and print literal curly brackets.
    
    return species_str

def plot(data_dict, output_file = None):
    '''
    Plot the data dictionary created by get_obs
    
    Args:
        data_dict (dict) :
            Dictionary of Pandas DataFrames output by get_obs. 
            Keys are site names.
    kwargs:
        output_file (str): specify a file path, if you want to save the figure
    '''
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    plots = []
    
    for site, df in data_dict.items():
        if site[0] != '.':
            if df is None:
                plots.append(mpatches.Patch(label="%s (no data)" %site, color = "white"))
            else:
                if "vmf" in df.columns:
                    error_col = "vmf"
                    errors = df[error_col]
                elif "dmf" in df.columns:
                    error_col = "dmf"
                    errors = df[error_col]
                else:
                    #TODO: Make this faster by duing a line plot, if there are no error bars
                    errors = df.mf*0.
                
                plots.append(plt.errorbar(df.index, df.mf,
                                          yerr = errors,
                                          linewidth = 0, 
                                          marker = '.', markersize = 6.,
                                          label = site))

    plt.legend(handles = plots)
    #plt.ylabel("%s (%s)" %(data_dict[".species"], data_dict[".units"]))
    plt.ylabel("%s (%s)" %(label_species(data_dict[".species"]), data_dict[".units"]))

    if output_file:
        plt.savefig(output_file)
    
    plt.show()
