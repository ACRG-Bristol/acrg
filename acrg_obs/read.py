# -*- coding: utf-8 -*-
"""
@author: chxmr
"""
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
from acrg_config.paths import paths
from acrg_utils import is_number
from acrg_obs.utils import obs_database
import numexpr as ne
from pathlib import Path

acrg_path = paths.acrg
obs_directory = paths.obs

#Get site info and species info from JSON files
with open(acrg_path / "acrg_species_info.json") as f:
    species_info=json.load(f)

with open(acrg_path / "acrg_site_info.json") as f:
    site_info=json.load(f, object_pairs_hook=OrderedDict)

    
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


def scale_convert(ds, species, to_scale):
    '''
    Convert to a new calibration scale, based on conversions in acrg_obs_scale_convert.csv
    
    Args:
        ds (xarray Dataset): Must contain an mf variable (mole fraction), and scale must be in global attributes
        species (str): species name
        to_scale (str): Calibration scale to convert to
    
    '''
    
    # If scale is already correct, return
    ds_scale = ds.attrs["scale"]
    if ds_scale == to_scale:
        return(ds)
    else:
        print(f"... converting scale to {to_scale}")
    
    scale_converter = pd.read_csv(acrg_path / "acrg_obs/acrg_obs_scale_convert.csv")
    scale_converter_scales = scale_converter[scale_converter.isin([species.upper(), ds_scale, to_scale])][["species", "scale1", "scale2"]].dropna(axis=0, how = "any")
    
    if len(scale_converter_scales) == 0:
        errorMessage = f'''Scales {ds_scale} and {to_scale} are not both in any one row in acrg_obs_scale_convert.csv for species {species}'''
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

    # scale_convert file has variable X in equations, so let's create it
    X = 1.
    scale_factor = ne.evaluate(converter[direction])
    ds["mf"].values *= scale_factor

    ds.attrs["scale"] = to_scale
    
    return(ds)


def get_single_site(site, species_in,
                    network = None,
                    start_date = "1900-01-01", end_date = "2100-01-01",
                    inlet = None, average = None,
                    instrument = None,
                    status_flag_unflagged = [0],
                    keep_missing = None,
                    file_path = None,
                    data_directory = None,
                    calibration_scale = None,
                    verbose = False):
    '''
    Get measurements from one site as a list of xarray datasets.
    If there are multiple instruments and inlets at a particular site, 
    note that the acrg_obs_defaults.csv file may be referenced to determine which instrument and inlet to use for each time period.
    If an inlet or instrument changes at some point during time period, multiple datasets will be returned,
    one for each inlet/instrument.
    
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
        file_path (str, optional) :
            Path to file. If this is used, network, inlet and instrument are ignored. site and species are still required.
            Default=None.
        data_directory (str or pathlib.Path, optional):
            User-defined obs directory
        calibration_scale (str, optional) :
            Convert to this calibration scale (original scale and new scale must both be in acrg_obs_scale_convert.csv)
    Returns:
        (list of xarray datasets):
            Mole fraction time series data as an xarray dataset, returned in a list. 
            Each list element is for a unique instrument and inlet.
            If either of these changes at some point during the timeseries, they are added as separate list elements.
            
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

        
    if file_path is not None:
        
        files_to_get = [(file_path, "na", "na"),]
        species_query = species

    else:
        
        # Open defaults file
        df_defaults = pd.read_csv(paths.acrg / "acrg_obs/acrg_obs_defaults.csv",
                                 parse_dates = ["startDate", "endDate"])
        df_defaults.dropna(axis = 0, how = "all", inplace = True)

        # Set early start date and late end date, if empty
        df_defaults["startDate"] = df_defaults["startDate"].fillna(pd.Timestamp("1900-01-01"))
        df_defaults["endDate"] = df_defaults["endDate"].fillna(pd.Timestamp("2100-01-01"))

        df_defaults.replace(np.nan, "%", inplace = True)

        # Read defaults database into memory
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
        if data_directory is None:
            # Attach the obs database
            c.execute("ATTACH DATABASE ? as obs", (str(paths.obs / "obs.db"),))
        else:
            # Attach the obs database
            c.execute("ATTACH DATABASE ? as obs", (str(Path(data_directory) / "obs.db"),))
            override_defaults = True

        # Change some variables for query
        start_date_query = pd.Timestamp(start_date).to_pydatetime()
        end_date_query = pd.Timestamp(end_date).to_pydatetime()
        species_query = species.replace("-", "").lower()

        if data_directory is None:
            # Run a couple of initial queries to see whether there are any defaults defined for a particular site
            df_defaults_for_site = pd.read_sql_query("SELECT * FROM defaults WHERE site = ? COLLATE NOCASE",
                                                     conn, params = (site,))
            df_defaults_for_site_species = pd.read_sql_query('''
                                                             SELECT * FROM defaults 
                                                             WHERE site = ? COLLATE NOCASE AND 
                                                             species = ? COLLATE NOCASE
                                                             ''', 
                                                             conn, params = (site, species_query))

            # Check if defaults need to be over-written
            if (inlet != None or instrument != None or network != None) and len(df_defaults_for_site) >= 1:
                print("... You've set either an inlet, instrument or network, overriding any defaults.")
                print("... Best to set all three of these, if you want to avoid ambiguity.")
                override_defaults = True
            else:
                override_defaults = False
        else:
            df_defaults_for_site = []

            
        # For debugging, get some info on the site and species in the obs.db
        df_site_species = pd.read_sql_query('''
                            SELECT files.inlet, files.instrument, files.startDate, files.endDate
                            FROM obs.files
                            WHERE files.species=? COLLATE NOCASE AND
                                  files.site=? COLLATE NOCASE''', 
                            conn, params = (species_query, site))
        
        # Read filenames from database
        if len(df_defaults_for_site) == 0 or override_defaults:
            # Query only the 'files' table table to determine which files to read

            print("... no defaults set")
            query = '''
                    SELECT files.filename, files.inlet, files.instrument
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
                    SELECT filename, inlet, instrument, defaultStartDate, defaultEndDate FROM
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

        # Store query with parameters inserted, for error messages
        query_params = query
        for p in params:
            query_params = query_params.replace("?", str(p), 1)

        if verbose:
            print(query_params)

        # Run query and get list of files
        files_to_get = list(c.execute(query, params))
        
        # close database
        conn.close()


    # Create empty output list (will contain xarray datasets)
    obs_files = []
    
    # If no files were found, print warning message explaining why
    if len(files_to_get) == 0 and len(df_defaults_for_site) == 0:
        
        print(f'''No files were found for this combination of inputs (site, species, date range, etc.)
                Have a look in {[data_directory, str(paths.obs)][data_directory == None]} and check that the expected files are there.

                ---------------------------------------------------
                obs.db entries for this site and species:
                '''.replace("  ", ""))
        print(df_site_species.to_string())

        return obs_files
    
    elif len(files_to_get) == 0 and len(df_defaults_for_site) > 1:
        
        print(f'''No files were found for this combination of inputs.
                
                Note that there are some defaults set for this site.
                
                1. Make sure that you've run acrg_obs.utils.obs_database() and that it executes correctly.
                This will make the obs.db file up-to-date with the filesystem.

                2. Check that there is some data for your selected inputs (e.g. within date range)
                in {[data_directory, str(paths.obs)][data_directory == None]}
                
                3. Take a look at acrg_obs_defaults.csv. A common issue is that defaults
                are added for some species at a particular site, but there's no
                instruction for the remaining species. If that's the case, 
                add a row to the file, which is empty, except for the site name. 
                Be careful that this doesn't lead to any ambiguity in the files retrieved.
                
                4. Another common issue is that the species name in the defaults file doesn't match
                the name in the files database. Make sure the species is the same as in the filename
                e.g. "cfc11", rather than "CFC-11"
                
                ---------------------------------------------------
                obs.db entries for this site:
                '''.replace("  ", ""))
                
        print(df_site_species.to_string())
        print('''---------------------------------------------------
                Defaults for this site:
                '''.replace("  ", ""))
        print(df_defaults_for_site.to_string())

        # Conditional statements from query
        # print((query_params.split("WHERE ")[1]).replace("  ", "")[:-2])
        
        return obs_files
        
    
    # Retrieve files
    for f in files_to_get:

        print(f"... reading {f[0]}")
        with xr.open_dataset(f[0]) as f_ds:
            ds = f_ds.load()

        # Remove any spaces from variable names:
        var_rename = {}
        for var in ds.variables:
            if " " in var:
                var_rename[var] = var.replace(" ", "_")
        ds = ds.rename_vars(var_rename)

        # If 5 elements, it means the query returned a start and end date
        # Otherwise, return whole dataset
        if len(f) == 5:
            if pd.Timestamp(start_date) > pd.Timestamp(f[3]):
                slice_start = pd.Timestamp(start_date)
            else:
                slice_start = pd.Timestamp(f[3])
            if pd.Timestamp(end_date) < pd.Timestamp(f[4]):
                slice_end = pd.Timestamp(end_date) - pd.Timedelta("1 ns")
            else:
                slice_end = pd.Timestamp(f[4]) - pd.Timedelta("1 ns")

        else:
            slice_start = pd.Timestamp(start_date)
            slice_end = pd.Timestamp(end_date) - pd.Timedelta("1 ns")

        if slice_start.round(freq = "T") != pd.Timestamp("1900-01-01") or  slice_end.round(freq = "T") != pd.Timestamp("2100-01-01"):
            print(f"... slicing from {slice_start} to {slice_end}")
            ds = ds.loc[dict(time = slice(slice_start, slice_end))]

        # If averaging is set, resample
        if average != None:
            
            if keep_missing == True:
                
                # Create a dataset with one element and NaNs to prepend or append
                ds_single_element = ds[dict(time = 0)]
                for v in ds_single_element.variables:
                    if v != "time":
                        ds_single_element[v].values = np.nan
                
                ds_concat = []
                
                # Pad with an empty entry at the start date
                if min(ds.time) > pd.Timestamp(start_date):
                    ds_single_element_start = ds_single_element.copy()
                    ds_single_element_start.time.values = pd.Timestamp(start_date)
                    ds_concat.append(ds_single_element_start)

                ds_concat.append(ds)
                
                # Pad with an empty entry at the end date
                if max(ds.time) < pd.Timestamp(end_date):
                    ds_single_element_end = ds_single_element.copy()
                    ds_single_element_end.time.values = pd.Timestamp(end_date) - pd.Timedelta("1ns")
                    ds_concat.append(ds_single_element_end)
                
                ds = xr.concat(ds_concat, dim="time")
                    
                # Now sort to get everything in the right order
                ds = ds.sortby("time")      
            
            # First, just do a mean resample on all variables
            print(f"... resampling to {average}")
            ds_resampled = ds.resample(time = average, keep_attrs = True
                                       ).mean(skipna=True)
            # keep_attrs doesn't seem to work for some reason, so manually copy
            ds_resampled.attrs = ds.attrs.copy()
            
            # For some variables, need a different type of resampling
            for var in ds.variables:
                if "repeatability" in var:
                    ds_resampled[var] = np.sqrt((ds[var]**2).resample(time = average).sum()) / \
                                                 ds[var].resample(time = average).count()

                # Copy over some attributes
                if "long_name" in ds[var].attrs:
                    ds_resampled[var].attrs["long_name"] = ds[var].attrs["long_name"]
                if "units" in ds[var].attrs:
                    ds_resampled[var].attrs["units"] = ds[var].attrs["units"]

            # Create a new variability variable, containing the standard deviation within the resampling period
            ds_resampled[f"{species_query}_variability"] = ds[species_query].resample(time = average,
                                                                                      keep_attrs = True).std(skipna=False)
            # If there are any periods where only one measurement was resampled, just use the median variability
            ds_resampled[f"{species_query}_variability"][ds_resampled[f"{species_query}_variability"] == 0.] = \
                                                                ds_resampled[f"{species_query}_variability"].median()
            # Create attributes for variability variable
            ds_resampled[f"{species_query}_variability"].attrs["long_name"] = f"{ds[species_query].attrs['long_name']}_variability"
            ds_resampled[f"{species_query}_variability"].attrs["units"] = ds[species_query].attrs["units"]

            # Resampling may introduce NaNs, so remove, if not keep_missing
            if keep_missing == False:
                ds_resampled = ds_resampled.dropna(dim = "time")
                    
            ds = ds_resampled.copy()

        # Rename variables
        rename = {}

        for var in ds.variables:
            if var.lower() == species_query.lower():
                rename[var] = "mf"
            if "repeatability" in var:
                rename[var] = "mf_repeatability"
            if "variability" in var:
                rename[var] = "mf_variability"
            if "number_of_observations" in var:
                rename[var] = "mf_number_of_observations"
            if "status_flag" in var:
                rename[var] = "status_flag"
            if "integration_flag" in var:
                rename[var] = "integration_flag"

        ds = ds.rename_vars(rename)

        # Append inlet and filename to attributes
        ds.attrs["filename"] = f[0]
        if f[1] == "%":
            first_network = next(iter(site_info[site]))
            ds.attrs["inlet"] = site_info[site][first_network]["height"][0]
        else:
            ds.attrs["inlet"] = f[1]
        ds.attrs["instrument"] = f[2]
        ds.attrs["species"] = species

        # Find calibration scale in file
        scale_count = 0
        for attr in ds.attrs:
            if "calibration" in attr.lower() or "scale" in attr.lower():
                scale = ds.attrs[attr]
                scale_count += 1
        
        if scale_count > 1:
            raise Exception("Ambiguous calibration scale: Check if the file has multiple global attributs containing the words calibration and scale")
        ds.attrs["scale"] = scale

        # Convert calibration scale, if needed
        if calibration_scale != None:
            ds = scale_convert(ds, species, calibration_scale)
        
        # Store dataset
        obs_files.append(ds)

    # Check if units match
    units = [f.mf.attrs["units"] for f in obs_files]
    if len(set(units)) > 1:
        errorMessage = f'''Units don't match for these files: {[(f.mf.attrs["units"],f.attrs["filename"]) for f in obs_files]}'''
        raise ValueError(errorMessage)

    scales = [f.attrs["scale"] for f in obs_files]
    if len(set(scales)) > 1: 
        print(f"WARNING: scales don't match for these files: {[(f.attrs['scale'],f.attrs['filename']) for f in obs_files]}")
        print("... suggest setting calibration_scale to convert")

    return obs_files


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
    
    if data_directory is None:
        gosat_directory = obs_directory / "GOSAT" / site
    else:
        gosat_directory = data_directory / "GOSAT" / site
    files = [f.name for f in gosat_directory.glob("*.nc")]

    files_date = [pd.to_datetime(f.split("_")[2][0:8]) for f in files]

    data = []
    for (f, d) in zip(files, files_date):
        if d >= pd.to_datetime(start_date) and d < pd.to_datetime(end_date):
            with xr.open_dataset(gosat_directory / f) as fxr:
                data.append(fxr.load())
    
    if len(data) == 0:
        return []
    
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
    data["mf_repeatability"] = data.xch4_uncertainty

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

    data = data.sortby("time")
    
    data.attrs["max_level"] = max_level
    if species.upper() == "CH4":
        data.mf.attrs["units"] = '1e-9'
        data.attrs["species"] = "CH4"
    if species.upper() == "CO2":
        data.mf.attrs["units"] = '1e-6'
        data.attrs["species"] = "CO2"

    data.attrs["scale"] = "GOSAT"

    # return single element list
    return [data,]


def get_obs(sites, species,
            start_date = "1900-01-01", end_date = "2100-01-01",
            inlet = None,
            average = None,
            keep_missing=False,
            network = None,
            instrument = None,
            status_flag_unflagged = None,
            max_level = None,
            data_directory = None,
            file_paths = None,
            calibration_scale = None):
    """
    The get_obs function retrieves obervations for a set of sites and species between start and end dates.
    This is essentially a wrapper function for get_single_site, to read in multiple sites.    
    
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
 
    For some sites where species are measured on multiple inlets and/or instruments,
    the acrg_obs_defaults.csv file may be read to determine which inlet, instrument to be returned for each time period.
    
    
    Examples:
    If we wanted to read in daily averaged methane data for Mace Head and Tacolneston for 2012:
    
        get_obs(sites=["MHD","TAC"],species="ch4",start_date="2012-01-01",end_date="2013-01-01",height=["10m","100m"],
                average=["24H","24H"],network=["AGAGE","DECC"])
    
    To just go with the defaults at these sites, do:
        get_obs(sites=["MHD","TAC"],species="ch4",start_date="2012-01-01",end_date="2013-01-01")
    
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
            Whether to keep missing data points or drop them when averaging.
        network (str/list, optional) : 
            Network for the site/instrument (must match number of sites). (optional)
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites). (optional)
        status_flag_unflagged (list, optional) : 
            The value to use when filtering by status_flag. Default = [0]
        max_level (int) : 
            Required for satellite data only. Maximum level to extract up to from within satellite data.
        data_directory (str, optional) :
            User-defined directory where data is stored. Note that if this is specified, an obs database will be created in this directory.
        file_paths (list of str, optional):
            Paths for specific files to read
        calibration_scale (str, optional) :
            Convert to this calibration scale (original scale and new scale must both be in acrg_obs_scale_convert.csv)            
    
    Returns:
        dict(list(xarray.Dataset)) : 
            Returns a dictionary with site codes as keys. For each site, a list of xarray datasets is returned (see get_single_site)
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

    inlet = check_list_and_length(inlet, sites, "inlet")
    average = check_list_and_length(average, sites, "average")
    network = check_list_and_length(network, sites, "network")
    instrument = check_list_and_length(instrument, sites, "instrument")

    if data_directory is not None:
        if not (Path(obs_directory) / "obs.db").exists():
            # Need to create a database in data_directory
            obs_database(data_directory = data_directory)
    
    # Get data
    obs = {}
    
    for si, site in enumerate(sites):
        print("Getting %s data for %s..." %(species, site))
        if "GOSAT" in site.upper():
            obs[site] = get_gosat(site, species,
                                   start_date = start_date, end_date = end_date,
                                   max_level = max_level,
                                   data_directory = data_directory)
        else:
            if file_paths is not None:
                file_path = file_paths[si]
            else:
                file_path = None
                        
            obs[site] = get_single_site(site, species, inlet = inlet[si],
                                   start_date = start_date, end_date = end_date,
                                   average = average[si],
                                   network = network[si],
                                   instrument = instrument[si],
                                   keep_missing = keep_missing,
                                   status_flag_unflagged = status_flag_unflagged[si],
                                   file_path = file_path,
                                   calibration_scale = calibration_scale,
                                   data_directory = data_directory)

    # Raise error if units don't match
    units = []
    for s, obs_site in obs.items():
        if len(obs_site) > 0:
            units.append(obs_site[0].mf.attrs["units"])
        else:
            units.append(None)

    if len(set(filter(None, units))) > 1:
        siteUnits = [': '.join([site, str(u)]) for (site, u) in zip(sites, units) if u is not None]
        errorMessage = f'''Units don't match for these sites: {siteUnits}'''
        raise ValueError(errorMessage)

    # Warning if scales don't match
    scales = []
    for s, obs_site in obs.items():
        if len(obs_site) > 0:
            scales.append(obs_site[0].attrs["scale"])
        else:
            scales.append(None)
    
    if len(set(filter(None, scales))) > 1:
        siteScales = [': '.join([site, scale]) for (site, scale) in zip(sites, scales) if scale is not None]
        warningMessage = '''WARNING: scales don't match for these sites:
                            %s''' % ', '.join(siteScales)
        print(warningMessage)

    return obs
    
