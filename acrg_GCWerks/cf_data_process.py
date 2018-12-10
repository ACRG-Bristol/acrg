# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 14:08:07 2015

@author: chxmr
"""

import numpy as np
import pandas as pd
import os
from os.path import join, split
from datetime import datetime as dt
from datetime import timedelta as td
import glob
import xarray as xray
import json
from os import getenv, stat
import time
import getpass
import fnmatch

# Read site info file
acrg_path = getenv("ACRG_PATH")
info_file = join(acrg_path,
                 "acrg_GCWerks/cf_data_process_parameters.json")
with open(info_file) as sf:
    params = json.load(sf)

site_info_file = join(acrg_path, "acrg_site_info.json")
with open(site_info_file) as sf:
    site_params = json.load(sf)

# Output unit strings (upper case for matching)
unit_species = {"CO2": "1e-6",
                "CH4": "1e-9",
                "C2H6": "1e-9",
                "N2O": "1e-9",
                "CO": "1e-9",
                "DCH4C13": "1",
                "DCH4D": "1",
                "DCO2C13": "1",
                "DCO2C14": "1",
                "DO2N2" : "1",
                "APO" : "1"}

# If units are non-standard, this attribute can be used
unit_species_long = {"DCH4C13": "permil",
                     "DCH4D": "permil",
                     "DCO2C13": "permil",
                     "DCO2C14": "permil",
                     "DO2N2" : "per meg",
                     "APO" : "per meg"}

unit_interpret = {"ppm": "1e-6",
                  "ppb": "1e-9",
                  "ppt": "1e-12",
                  "ppq": "1e-15",
                  "else": "unknown"}

# Default calibration scales
# TODO: Remove this? seems dangerous
scales = {"CO2": "NOAA-2007",
          "CH4": "NOAA-2004A",
          "N2O": "SIO-98",
          "CO": "Unknown"}


# For species which need more than just a hyphen removing or changing to lower case
# First element of list is the output variable name,
# second is the long name for variable standard_name and long_name
# Keys are upper case
species_translator = {"CO2": ["co2", "carbon_dioxide"],
                      "CH4": ["ch4", "methane"],
                      "ETHANE": ["c2h6", "ethane"],
                      "PROPANE": ["c3h8", "propane"],
                      "C-PROPANE": ["cc3h8", "cpropane"],
                      "BENZENE": ["c6h6", "benzene"],
                      "TOLUENE": ["c6h5ch3", "methylbenzene"],
                      "ETHENE": ["c2f4", "ethene"],
                      "ETHYNE": ["c2f2", "ethyne"],
                      "N2O": ["n2o", "nitrous_oxide"],
                      "CO": ["co", "carbon_monoxide"],
                      "H-1211": ["halon1211", "halon1211"],
                      "H-1301": ["halon1301", "halon1301"],
                      "H-2402": ["halon2402", "halon2402"],
                      "PCE": ["c2cl4", "tetrachloroethylene"],
                      "TCE": ["c2hcl3", "trichloroethylene"],
                      "PFC-116": ["c2f6", "hexafluoroethane"],
                      "PFC-218": ["c3f8", "octafluoropropane"],
                      "PFC-318": ["c4f8", "cyclooctafluorobutane"],
                      "F-113": ["cfc113", "cfc113"],
                      "H2_PDD": ["h2", "hydrogen"],
                      "NE_PDD": ["Ne", "neon"],                      
                      "DO2N2": ["do2n2", "ratio_of_oxygen_to_nitrogen"],
                      "DCH4C13": ["dch4c13", "delta_ch4_c13"],
                      "DCH4D": ["dch4d", "delta_ch4_d"],
                      "DCO2C13": ["dco2c13", "delta_co2_c13"],
                      "DCO2C14": ["dco2c14", "delta_co2_c14"],
                      "APO": ["apo", "atmospheric_potential_oxygen"]
                      }

# Translate header strings
crds_header_string_interpret = {"C": "",
                                "stdev": " variability",
                                "N": " number_of_observations"}

def parser_YYMMDD(yymmdd):
    return dt.strptime(yymmdd, '%y%m%d')


def get_directories(default_input_directory,
                    default_output_directory,
                    user_specified_input_directory = None,
                    user_specified_output_directory = None):

    # If an output directory directory is set, use that, otherwise from json file
    if user_specified_output_directory:
        output_folder = user_specified_output_directory
    else:
        output_folder = default_output_directory

    # If an input directory directory is set, use that, otherwise from json file
    if user_specified_input_directory:
        if user_specified_output_directory is None:
            print("EXITING: You must also specify an output folder if you specify an input")
            return None
        data_folder = user_specified_input_directory
    else:
        data_folder = default_input_directory

    return data_folder, output_folder


def site_info_attributes(site):

    attributes = {}
    attributes_list = {"longitude": "station_longitude",
                       "latitude": "station_latitude",
                       "long_name": "station_long_name",
                       "height_station_masl": "station_height_masl"}

    if site in site_params.keys():
        for at in attributes_list.keys():
            if at in site_params[site].keys():
                attributes[attributes_list[at]] = site_params[site][at]
        return attributes
    else:
        return None

def attributes(ds, species, site,
               global_attributes = None,
               units = None,
               scale = None,
               sampling_period = None):
    """
    Format attributes for netCDF file
    Attributes of xarray DataSet are modified, and variable names are changed
    
    If the species is a standard mole fraction then either:
        - species name will used in lower case in the file and variable names
            but with any hyphens taken out
        - name will be changed according to the species_translator dictionary
    
    If the species is isotopic data or a non-standard variable (e.g. APO):
        - Isotopes species names should begin with a "D"
            (Annoyingly, the code currently picks up "Desflurane" too. I've
             fixed this for now, but if we get a lot of other "D" species, we 
             should make this better)
        - I suggest naming for isotopologues should be d<species><isotope>, e.g.
            dCH4C13, or dCO2C14
        - Any non-standard variables should be listed in the species_translator
            dictionary
    
    Args:
        ds (xarray dataset): Should contain variables such as "ch4", "ch4 repeatability".
            Must have a "time" dimension.
        species (string): Species name. e.g. "CH4", "HFC-134a", "dCH4C13"
        site (string): Three-letter site code
        
        global_attribuates (dict, optional): Dictionary containing any info you want to 
            add to the file header (e.g. {"Contact": "Matt Rigby"})
        units (string, optional): This routine will try to guess the units
            unless this is specified. Options are in units_interpret
        scale (string, optional): Calibration scale for file header.
        sampling_period (int, optional): Number of seconds for which air 
            sample is taken. Only for time variable attribute
    """

    # Rename species
    for key in ds.keys():
        if species in key:
            if species.upper() in species_translator.keys():
                # Rename based on species_translator, if available
                species_out = species_translator[species.upper()][0]
            else:
                # Rename species to be lower case and without hyphens
                species_out = species.lower().replace("-", "")
            ds.rename({key: key.replace(species, species_out)}, inplace = True)

    # Global attributes
    #############################################
    if global_attributes is None:
        global_attributes = {}

    # Add some defaults
    for key, value in params["global_attributes"].iteritems():
        global_attributes[key] = value
    global_attributes["File created"] = str(dt.now())
    global_attributes["Conventions"] = "CF-1.6"

    # Add user
    global_attributes["Processed by"] = "%s@bristol.ac.uk" % getpass.getuser()    


    for key, values in global_attributes.iteritems():
        ds.attrs[key] = values

    # Add some site attributes
    global_attributes_site = site_info_attributes(site.upper())
    if global_attributes_site is not None:
        for key, values in global_attributes_site.iteritems():
            ds.attrs[key] = values

    # Add calibration scale
    if scale:
        ds.attrs["Calibration_scale"] = scale
    else:
        if species.upper() in scales.keys():
            ds.attrs["Calibration_scale"] = scales[species.upper()]
        else:
            ds.attrs["Calibration_scale"] = "unknown"

    # Add species name
    ds.attrs["species"] = species_out

    # Species-specific attributes
    #############################################

    # Long name
    if (species.upper()[0] == "D" and species.upper() != "DESFLURANE") or species.upper() == "APO":
        sp_long = species_translator[species.upper()][1]
    elif species.upper() in species_translator.keys():
        sp_long = "mole_fraction_of_" + species_translator[species.upper()][1] + "_in_air"
    else:
        sp_long = "mole_fraction_of_" + species_out + "_in_air"

    ancillary_variables = ""

    for key in ds.keys():

        if species_out in key:

            # Standard name attribute
            #ds[key].attrs["standard_name"]=key.replace(species_out, sp_long)
            ds[key].attrs["long_name"]=key.replace(species_out, sp_long)

            # If units are required for variable, add attribute
            if (key == species_out) or \
                ("variability" in key) or \
                ("repeatability" in key):
                if units is None:
                    ds[key].attrs["units"] = unit_species[species.upper()]
                else:
                    if units in unit_interpret.keys():
                        ds[key].attrs["units"] = unit_interpret[units]
                    else:
                        ds[key].attrs["units"] = unit_interpret["else"]

                # if units are non-standard, add explanation
                if species.upper() in unit_species_long.keys():
                    ds[key].attrs["units_description"] = unit_species_long[species.upper()]

            # Add to list of ancilliary variables                    
            if key != species_out:
                ancillary_variables += " " + key

    # Write ancilliary variable list
    ds[species_out].attrs["ancilliary_variables"] = ancillary_variables.strip()

    # Add quality flag attributes
    ##################################

    flag_key = [key for key in ds.keys() if " status_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = unflagged, 1 = flagged",
                              "long_name":
                              ds[species_out].attrs["long_name"] + " status_flag"}

    # Add integration flag attributes
    ##################################

    flag_key = [key for key in ds.keys() if " integration_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = area, 1 = height",
                              "standard_name":
                              ds[species_out].attrs["long_name"] + " integration_flag",
                              "comment":
                              "GC peak integration method (by height or by area). " +
                              "Does not indicate data quality"}

    # Set time encoding
    #########################################

    # Check if there are duplicate time stamps
    if len(set(ds.time.values)) < len(ds.time.values):
        print("WARNING. Dupliate time stamps")

    first_year = str(ds.time.to_pandas().index.to_pydatetime()[0].year)

    ds.time.encoding = {"units": "seconds since " + \
                        first_year + "-01-01 00:00:00"}
    ds.time.attrs["label"] = "left"
    ds.time.attrs["comment"] = "Time stamp corresponds to beginning of sampling period. " + \
                               "Time since midnight UTC of reference date. " + \
                               "Note that sampling periods are approximate."
    if sampling_period:
        ds.time.attrs["sampling_period_seconds"] = sampling_period

    return ds


def output_filename(output_directory,
                    network,
                    instrument,
                    site,
                    time_start,
                    species,
                    inlet = None,
                    version = None):
    """
    Create an output filename in the format 
    output_directory/site/network-instrument_site_YYYYMMDD_species[-inlet]-version.nc
    
    Also creates a site directory in output_directory, if one doesn't exist
    
    Args:
        
        inlet (string, optional): Label for inlet. If none supplied, assumes that there
            is only one inlet, and no inlet info is added to file name
        version (string, optional): Version number for file. If none supplied, will
            use YYYYMMDD for day processed
    
    """
    
    # Check if directory exists. If not, create it
    if not os.path.exists(join(output_directory, site)):
        os.makedirs(join(output_directory, site))
    
    # Create suffix. If inlet is specified, append to species name
    suffix = species
    if inlet is not None:
        suffix += "-" + inlet
    if version is not None:
        suffix += "-" + version
    else:
        suffix += "-" + time.strftime("%Y%m%d")
    
    # Return output filename
    return join(output_directory,
                "%s/%s-%s_%s_%04i%02i%02i_%s.nc" %(site,
                                                   network,
                                                   instrument,
                                                   site,
                                                   time_start.year,
                                                   time_start.month,
                                                   time_start.day,
                                                   suffix))


# ICOS
########################################################

def icos_data_read(data_file, species):

    print("Reading " + data_file)

    # Find out how many header lines there are
    nheader = 0
    with open(data_file, "rb") as f:
        for l in f:
            if l[0] != "#":
                break
            nheader += 1

    # Read CSV file
    df =  pd.read_csv(data_file,
                      skiprows = nheader-1,
                      parse_dates = {"time": ["Year", "Month", "Day", "Hour", "Minute"]},
#                      date_parser = lambda s: pd.to_datetime(s, format = "%Y %m %d %H %M"),
                      index_col = "time",
                      sep = ";",
                      usecols = ["Day", "Month", "Year", "Hour", "Minute",
                                 str(species.lower()), "SamplingHeight",
                                 "Stdev", "NbPoints"],
                      dtype = {"Day": np.int,
                               "Month": np.int,
                               "Year": np.int,
                               "Hour": np.int,
                               "Minute": np.int,
                               species.lower(): np.float,
                               "Stdev": np.float,
                               "SamplingHeight": np.float,
                               "NbPoints": np.int},
                      na_values = "-999.99")

    # Format time
    df.index = pd.to_datetime(df.index, format = "%Y %m %d %H %M")

    df = df[df[species.lower()] >= 0.]

    # Remove duplicate indices
    df.reset_index(inplace = True)
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')

    # Rename columns
    df.rename(columns = {species.lower(): species.upper(),
                         "Stdev": species.upper() + " variability",
                         "NbPoints": species.upper() + " number_of_observations"},
               inplace = True)

    df.index.name = "time"

    # Convert to Dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())

    return ds


def icos(site, network = "ICOS",
         input_directory = None,
         output_directory = None):

    # Get directories and site strings
    params_icos = params["ICOS"]
    site_string = params_icos[site]["gcwerks_site_name"]

    data_folder, output_folder = \
            get_directories(params_icos["directory"].replace("%site", site_string),
                            params_icos["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    # Search for species and inlets from file names
    data_file_search = join(data_folder, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    data_file_names = [split(f)[1] for f in data_files]
    species_and_inlet = [(f.split(".")[1], f.split(".")[-2]) \
                         for f in data_file_names]

    for i, (species, inlet) in enumerate(species_and_inlet):

        if stat(data_files[i]).st_size > 0:

            # Create Pandas dataframe
            ds = icos_data_read(data_files[i], species.upper())

            # Sort out attributes
            global_attributes = params_icos[site.upper()]["global_attributes"]
            global_attributes["inlet_height_magl"] = float(params_icos[site]["inlet_rename"][inlet][:-1])

            ds = attributes(ds,
                            species.upper(),
                            site.upper(),
                            global_attributes = global_attributes,
                            sampling_period = 60)

            # Write file
            nc_filename = output_filename(output_folder,
                                          network,
                                          "CRDS",
                                          site.upper(),
                                          ds.time.to_pandas().index.to_pydatetime()[0],
                                          ds.species,
                                          params_icos[site]["inlet_rename"][inlet])

            ds.to_netcdf(nc_filename)

            print("Written " + nc_filename)

        else:
            print("Skipping empty file: %s" % data_files[i])

# GC FUNCTIONS
###############################################################

def gc_data_read(dotC_file, scale = {}, units = {}):

    species = []

    # Read header
    header = pd.read_csv(dotC_file,
                         skiprows=2,
                         nrows=2,
                         header = None,
                         sep=r"\s+")

    # Read data
    df = pd.read_csv(dotC_file,
                     skiprows=4,
                     sep=r"\s+")

    # Time index
    time = [dt(df.yyyy[i], df.mm[i], df.dd[i], df.hh[i], df.mi[i]) \
            for i in range(len(df))]
    df.index = time

    # Drop duplicates
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')

    # Rename flag column with species name
    for i, key in enumerate(df.keys()):
        if key[0:4] == "Flag":
            quality_flag = []
            area_height_flag = []
            for flags in df[key].values:

                # Quality flag
                if flags[0] == "-":
                    quality_flag.append(0)
                else:
                    quality_flag.append(1)

                # Area/height
                if flags[1] == "-":
                    area_height_flag.append(0)  # Area
                else:
                    area_height_flag.append(1)  # Height

            df = df.rename(columns = {key: df.keys()[i-1] + "_flag"})
            df[df.keys()[i-1] + " status_flag"] = quality_flag
            df[df.keys()[i-1] + " integration_flag"] = area_height_flag
            scale[df.keys()[i-1]] = header[i-1][0]
            units[df.keys()[i-1]] = header[i-1][1]
            species.append(df.keys()[i-1])

    return df, species, units, scale


def gc_precisions_read(precisions_file):

    # Read precision species
    precision_species = list(pd.read_csv(precisions_file,
                                         skiprows=3,
                                         nrows = 1,
                                         header = None,
                                         sep=r"\s+").values[0][1:])

    # Read precisions
    precision = pd.read_csv(precisions_file,
                            skiprows=5,
                            header = None,
                            sep=r"\s+", dtype = str,
                            index_col = 0,
                            parse_dates = True,
                            date_parser = parser_YYMMDD)

    # Rename index column
    precision.index.names = ["index"]

    # Drop duplicates
    precision = precision.reset_index().drop_duplicates(subset='index').set_index('index')

    return precision, precision_species


def gc(site, instrument, network,
       input_directory = None,
       output_directory = None):
    """
    Process GC data per site and instrument
    Instruments can be:
        "GCMD": GC multi-detector (output will be labeled GC-FID or GC-ECD)
        "GCMS": GC ADS (output GC-ADS)
        "medusa": GC medusa (output GC-MEDUSA)

    Network is the network name for output file.
    """

    site_gcwerks = params["GC"][site]["gcwerks_site_name"]
    instrument_gcwerks = params["GC"]["instruments"][instrument]

    data_folder, output_folder = \
            get_directories(params["GC"]["directory"][instrument],
                            params["GC"]["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    search_strings = []
    for suffix in params["GC"]["instruments_suffix"][instrument]:
        # Search string
        search_string = join(data_folder,
                             site_gcwerks + \
                             instrument_gcwerks + \
                             suffix + ".??.C")
        search_strings.append(search_string)

        data_files = sorted(glob.glob(search_string))
        if len(data_files) > 0:
            break

    # Error if can't find files
    if len(data_files) == 0.:
        print("ERROR: can't find any files: " + \
              ",\r".join(search_strings))
        return None

    precision_files = [data_file[0:-2] + ".precisions.C" \
                        for data_file in data_files]

    dfs = []
    scale = {}
    units = {}

    for fi, data_file in enumerate(data_files):

        print("Reading " + data_file)

        # Get observations
        df, species, units, scale = gc_data_read(data_file,
                                                 scale = scale,
                                                 units = units)

        # Get precision
        precision, precision_species = gc_precisions_read(precision_files[fi])

        # Merge precisions into dataframe
        for sp in species:
            precision_index = precision_species.index(sp)*2+1
            df[sp + " repeatability"] = precision[precision_index].\
                                            astype(float).\
                                            reindex_like(df, "pad")

        dfs.append(df)

    # Concatenate
    dfs = pd.concat(dfs).sort_index()

    # Apply timestamp offset so that timestamp reflects start of sampling
    time = dfs.index.values
    time_offset = np.timedelta64(td(seconds = params["GC"]["timestamp_correct_seconds"][instrument]))
    time = [t + time_offset for t in time]
    dfs.index = time

    # Label time index
    dfs.index.name = "time"

    # Convert to xray dataset
    ds = xray.Dataset.from_dataframe(dfs)

    # Get species from scale dictionary
    species = scale.keys()

    inlets = params["GC"][site]["inlets"]

    # Process each species in file
    for sp in species:

        global_attributes = params["GC"][site.upper()]["global_attributes"]
        global_attributes["comment"] = params["GC"]["comment"][instrument]


        # Now go through each inlet (if required)
        for inleti, inlet in enumerate(inlets):

            # There is only one inlet, just use all data, and don't lable inlet in filename
            if (inlet == "any") or (inlet == "air"):
                
                print("Processing %s, assuming single inlet..." %sp)
                
                ds_sp = ds[[sp,
                            sp + " repeatability",
                            sp + " status_flag",
                            sp + " integration_flag"]]
                
                # No inlet label in file name
                inlet_label = None
                
                # Add list of all inlet lables, for record
                global_attributes["inlet_height_magl"] = \
                                    ", ".join(set(ds["Inlet"].values))

            else:
                # Get specific inlet
                
                print("Processing " + sp + ", " + inlet + "...")
                
                # if inlet is in the format "date_YYYYMMDD_YYYYMMDD", split by date
                if inlet[0:4] == "date":
                    dates = inlet.split("_")[1:]
                    slice_dict = dict(time = slice(dates[0], dates[1]))
                    ds_sliced = ds.loc[slice_dict]
                    ds_sp = ds_sliced[[sp,
                                       sp + " repeatability",
                                       sp + " status_flag",
                                       sp + " integration_flag"]]
    
                    global_attributes["inlet_height_magl"] = \
                        ", ".join(set(ds_sliced.Inlet.values))
                
                else:
                    
                    # Use UNIX pattern matching to find matching inlets
                    # select_inlet is a list of True or False
                    select_inlet = [fnmatch.fnmatch(i, inlet) for i in ds.Inlet.values]
                    # now create a DataArray of True or False
                    select_ds = xray.DataArray(select_inlet, coords = [ds.time],
                                               dims = ["time"])
                    
                    # sub-set ds
                    ds_sp = ds.where(select_ds, drop = True)[[sp,
                                                              sp + " repeatability",
                                                              sp + " status_flag",
                                                              sp + " integration_flag"]]

                    # Keep a record of which inlets were included
                    global_attributes["inlet_height_magl"] = \
                        ", ".join(set(ds.where(select_ds, drop = True).Inlet.values))
    
                # re-label inlet if required
                if "inlet_label" in params["GC"][site].keys():
                    inlet_label = params["GC"][site]["inlet_label"][inleti]
                else:
                   inlet_label = inlet
                  
    

            # Drop NaNs
            ds_sp = ds_sp.dropna("time")

            if len(ds_sp.time) == 0:

                print("... no data in file, skipping " + sp)

            else:

                # Sort out attributes
                ds_sp = attributes(ds_sp, sp, site.upper(),
                                   global_attributes = global_attributes,
                                   units = units[sp],
                                   scale = scale[sp],
                                   sampling_period = params["GC"]["sampling_period"][instrument])

                # Get instrument name for output
                if sp.upper() in params["GC"]["instruments_out"][instrument]:
                    instrument_out = params["GC"]["instruments_out"][instrument][sp]
                else:
                    instrument_out = params["GC"]["instruments_out"][instrument]["else"]

                # Write file
                print(output_folder)
                nc_filename = output_filename(output_folder,
                                              network,
                                              instrument_out,
                                              site.upper(),
                                              ds_sp.time.to_pandas().index.to_pydatetime()[0],
                                              ds_sp.species,
                                              inlet = inlet_label)
                print("Writing... " + nc_filename)
                ds_sp.to_netcdf(nc_filename)
                print("... written.")



def crds_data_read(data_file):

    print("Reading " + data_file)

    # Read file header
    df_header = pd.read_csv(data_file,
                         skiprows=1,
                         nrows = 2,
                         header = None,
                         sep=r"\s+")

    header = []
    species = []

    # Create header list
    for i in df_header.columns:
        if df_header[i][0] != '-':
            header.append(df_header[i][0].upper() + \
                          crds_header_string_interpret[df_header[i][1]])
            if df_header[i][1] == "C":
                species.append(df_header[i][0].upper())
        else:
            header.append(df_header[i][1].upper())

    # Read data
    df = pd.read_csv(data_file,
                     skiprows=4,
                     header = None,
                     sep=r"\s+",
                     names = header,
                     dtype = {"DATE": str, "TIME": str})

    # Interpret time
    time = [dt(2000 + int(date[0:2]),
                      int(date[2:4]),
                      int(date[4:]),
                      int(time[0:2]),
                      int(time[2:4]),
                      int(time[4:])) \
            for date, time in zip(df["DATE"].values, df["TIME"].values)]
    df.index = time

    # Remove duplicate indices
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')

    # Convert to Dataset
    df.index.name = "time"
    ds = xray.Dataset.from_dataframe(df.sort_index())

    return ds, species


def crds(site, network,
         input_directory = None,
         output_directory = None,
         version = None):
    """
    Process CRDS data

    site : Three letter site code
    network : Network string only for output
    """
    params_crds = params["CRDS"]

    site_string = params_crds[site]["gcwerks_site_name"]

    data_folder, output_folder = \
            get_directories(params_crds["directory"].replace("%site", site_string),
                            params["CRDS"]["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    # Search for species and inlets from file names
    data_file_search = join(data_folder, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    inlets = [f.split(".")[-2] for f in data_files]

    for i, inlet in enumerate(inlets):

        # Create Pandas dataframe
        ds, species = crds_data_read(data_files[i])

        # Write netCDF file for each species
        for sp in species:

            # Species-specific dataset
            ds_sp = ds[[sp,
                        sp + " variability",
                        sp + " number_of_observations"]]
            ds_sp = ds_sp.dropna("time")

            global_attributes = params_crds[site]["global_attributes"]
            global_attributes["inlet_height_magl"] = float(inlet[0:-1])
            global_attributes["comment"] = params_crds["comment"]

            ds_sp = attributes(ds_sp, sp, site.upper(),
                               global_attributes = global_attributes,
                               scale = scales[sp],
                               sampling_period=60)

            # Write file
            nc_filename = output_filename(output_folder,
                                          network,
                                          "CRDS",
                                          site.upper(),
                                          ds_sp.time.to_pandas().index.to_pydatetime()[0],
                                          ds_sp.species,
                                          inlet = inlet,
                                          version = version)
            print("Writing " + nc_filename)
            ds_sp.to_netcdf(nc_filename)
            print("... written.")


def ale_gage(site, network):
    """
    Process Georgia Tech ALE or GAGE observations
    
    Args:
        site (str): ADR, RPB, ORG, SMO, CGO or MHD
        network (str): ALE or GAGE

    """
    
    import fortranformat as ff

    ale_directory = "/dagage2/agage/summary/git/ale_new/complete/"
    gage_directory = "/dagage2/agage/summary/git/gage_new/complete/"

    output_directory = "/data/shared/obs_2018/"

    site_translate = {"ADR": "adrigole",
                      "RPB": "barbados",
                      "ORG": "oregon",
                      "SMO": "samoa",
                      "CGO": "tasmania",
                      "MHD": "macehead"}

    if network == "ALE":
        data_directory = ale_directory
    elif network == "GAGE":
        data_directory = gage_directory
    else:
        print("Network needs to be ALE or GAGE")
        return None


    fnames = sorted(glob.glob(join(data_directory,
                                   site_translate[site] + "/" + site + "*.dap")))


    formatter = ff.FortranRecordReader('(F10.5, 2I4,I6, 2I4,I6,1X,10(F10.3,a1))')


    dfs = []
    for fname in fnames:

        print("Reading... " + fname)

        header = []

        with open(fname) as f:
            for i in range(6):
                header.append(f.readline())

            lines = f.readlines()

        scales = header[-3].split()
        units = header[-2].split()
        species = header[-1].split()

        dayi = species.index("DD")
        monthi = species.index("MM")
        yeari = species.index("YYYY")
        houri = species.index("hh")
        mini = species.index("min")

        data = []
        time = []

        for line in lines:
            data_line = formatter.read(line)

            if data_line[mini] < 60 and data_line[houri] < 24:
                data.append([d for d in data_line if d != " " and \
                                                     d != None and \
                                                     d != "P"])
                time.append(dt(data_line[yeari],
                               data_line[monthi],
                               data_line[dayi],
                               data_line[houri],
                               data_line[mini]))

        data = np.vstack(data)
        data = data[:, 7:]

        df = pd.DataFrame(data = data, columns = species[7:], index = time)
        df.replace(to_replace = 0., value=np.NaN, inplace = True)
        dfs.append(df)

    scales = scales[7:]
    units = units[7:]
    species = species[7:]

    df = pd.concat(dfs)

    # Write netCDF file for each species
    for si, sp in enumerate(species):

        # Remove duplicate indices
        df_sp = df[sp].reset_index().drop_duplicates(subset='index').set_index('index')

        # Convert to Dataset
        df_sp.index.name = "time"
        ds = xray.Dataset.from_dataframe(df_sp.sort_index())
        ds = ds.dropna("time")

        ds = attributes(ds, sp, site.upper(),
                       scale = scales[si],
                       sampling_period=60,
                       units = units[si])

        # Write file
        nc_filename = output_filename(output_directory,
                                      network,
                                      "GCECD",
                                      site.upper(),
                                      ds.time.to_pandas().index.to_pydatetime()[0],
                                      ds.species,
                                      None)
        print("Writing " + nc_filename)
        ds.to_netcdf(nc_filename)
        print("... written.")


def mhd_o3():

    channels = ["channel1", "channel0", "channel2"]
    base_directory = "/dagage2/agage/macehead-ozone/results/export/"

    df = []

    for channel in channels:

        files_channel = sorted(glob.glob(join(base_directory, channel, "*.csv")))

        for f in files_channel:

            df.append(pd.read_csv(f, sep=",",
                                  names = ["datetime",
                                           "ozone",
                                           "ozone_variability",
                                           "ozone_number_samples"],
                                  na_values = "NA",
                                  index_col = "datetime",
                                  parse_dates = ["datetime"]))
            df[-1].dropna(inplace = True)

    df = pd.concat(df)
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')
    df.sort_index(inplace = True)

    # Convert to Dataset
    df.index.name = "time"
    ds = xray.Dataset.from_dataframe(df)

    ds = attributes(ds,
                    "ozone",
                    "MHD",
                    scale = "SCALE",
                    sampling_period=60*60,
                    units = "ppb")

    # Write file
    nc_filename = output_filename("/dagage2/agage/metoffice/processed_observations_2018",
                                  "AURN",
                                  "thermo",
                                  "MHD",
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  site_params["MHD"]["height"][0])
    print("Writing " + nc_filename)
    ds.to_netcdf(nc_filename)
    print("... written.")


def decc_data_freeze():

    input_directory = "/dagage2/agage/summary/gccompare-net/snapshot/current-frozendata/data-net/"
    output_directory = "/dagage2/agage/summary/gccompare-net/snapshot/current-frozendata/data-net/processed/"

    # ICOS
    icos("MHD", network = "ICOS", input_directory = input_directory, output_directory = output_directory)

    # GAUGE CRDS data
    crds("HFD", "GAUGE", input_directory = input_directory, output_directory = output_directory)
    crds("BSD", "GAUGE", input_directory = input_directory, output_directory = output_directory)

    # GAUGE GC data
    gc("BSD", "GCMD", "GAUGE", input_directory = input_directory, output_directory = output_directory)
    gc("HFD", "GCMD", "GAUGE", input_directory = input_directory, output_directory = output_directory)

    # DECC CRDS data
    crds("TTA", "DECC", input_directory = input_directory, output_directory = output_directory)
    crds("RGL", "DECC", input_directory = input_directory, output_directory = output_directory)
    crds("TAC", "DECC", input_directory = input_directory, output_directory = output_directory)

    # DECC GC data
    gc("TAC", "GCMD", "DECC", input_directory = input_directory, output_directory = output_directory)
    gc("RGL", "GCMD", "DECC", input_directory = input_directory, output_directory = output_directory)

    # DECC Medusa
    gc("TAC", "medusa", "DECC", input_directory = input_directory, output_directory = output_directory)

    # AGAGE GC data
    gc("MHD", "GCMD", "AGAGE", input_directory = input_directory, output_directory = output_directory)

    # AGAGE GCMS data
    gc("MHD", "GCMS", "AGAGE", input_directory = input_directory, output_directory = output_directory)

    # AGAGE Medusa
    gc("MHD", "medusa", "AGAGE", input_directory = input_directory, output_directory = output_directory)





if __name__ == "__main__":

#    # AGAGE Medusa
#    gc("MHD", "medusa", "AGAGE")
#    gc("CGO", "medusa", "AGAGE")
#    gc("GSN", "medusa", "AGAGE")
#    gc("SDZ", "medusa", "AGAGE")
#    gc("THD", "medusa", "AGAGE")
#    gc("RPB", "medusa", "AGAGE")
#    gc("SMO", "medusa", "AGAGE")
#    gc("SIO", "medusa", "AGAGE")
#    gc("JFJ", "medusa", "AGAGE")
#    gc("CMN", "medusa", "AGAGE")
#    gc("ZEP", "medusa", "AGAGE")
#
#    # AGAGE GC data
#    gc("RPB", "GCMD", "AGAGE")
#    gc("CGO", "GCMD", "AGAGE")
#    gc("MHD", "GCMD", "AGAGE")
#    gc("SMO", "GCMD", "AGAGE")
#    gc("THD", "GCMD", "AGAGE")
#
#    # AGAGE GCMS data
#    gc("CGO", "GCMS", "AGAGE")
#    gc("MHD", "GCMS", "AGAGE")
#    gc("RPB", "GCMS", "AGAGE")
#    gc("SMO", "GCMS", "AGAGE")
#    gc("THD", "GCMS", "AGAGE")
#    gc("JFJ", "GCMS", "AGAGE")
#    gc("CMN", "GCMS", "AGAGE")
#    gc("ZEP", "GCMS", "AGAGE")
#
    # AGAGE CRDS data
    crds("RPB", "AGAGE")

    # ICOS
#    icos("TTA", network = "DECC")
    icos("MHD", network = "ICOS")

    # GAUGE CRDS data
    crds("HFD", "DECC")
    crds("BSD", "DECC")

    # GAUGE GC data
    gc("BSD", "GCMD", "DECC")
    gc("HFD", "GCMD", "DECC")

    # DECC CRDS data
    crds("TTA", "DECC")
    crds("RGL", "DECC")
    crds("TAC", "DECC")

    # DECC GC data
    gc("TAC", "GCMD", "DECC")
    gc("RGL", "GCMD", "DECC")

    # DECC Medusa
    gc("TAC", "medusa", "DECC")


#    # Copy files
#    networks = ["AGAGE", "GAUGE", "DECC", "ICOS"]
#    src_dir = "/dagage2/agage/metoffice/processed_observations_2018"
#    dest_dir = "/data/shared/obs_2018"
#
#    for network in networks:
#        files = glob.glob(join(src_dir, network, "*.nc"))
#        for f in files:
#            print("Copying %s..." % (split(f)[-1]))
#            shutil.copy(f, join(dest_dir, network))
