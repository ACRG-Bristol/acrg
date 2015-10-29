# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 14:08:07 2015

@author: chxmr
"""

import numpy as np
import pandas as pd
from os.path import join
from datetime import datetime as dt
import glob
import xray
import json
from os import getenv


# Read site info file
acrg_path = getenv("ACRG_PATH")
info_file = join(acrg_path,
                 "acrg_GCWerks/cf_data_process_parameters.json")
with open(info_file) as sf:
    params = json.load(sf)

# Output unit strings
unit_species = {"CO2": "1e-6",
                "CH4": "1e-9",
                "N2O": "1e-9",
                "CO": "1e-6"}

unit_interpret = {"ppm": "1e-6",
                 "ppb": "1e-9",
                 "ppt": "1e-12",
                 "else": "unknown"}

# Calibration scales
scales = {"CO2": "NOAA-2007",
          "CH4": "NOAA-2004",
          "N2O": "SIO-98",
          "CO": "Unknown"}

# Species long-names for output
species_long = {"CO2": "carbon_dioxide",
                "CH4": "methane",
                "N2O": "nitrous_oxide",
                "CO": "carbon_monixide"}

# Translate header strings
crds_header_string_interpret = {"C": "",
                                "stdev": "_variability",
                                "N": "_number_of_observations"}


def attributes(ds, species, global_attributes = None,
               units = None, scale = None):
    """
    Format attributes for netCDF file
    """

    # Rename species to be lower case and without hyphens
    species_out = species.lower().replace("-", "")
    for key in ds.keys():
        if species in key:
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
    
    for key, values in global_attributes.iteritems():
        ds.attrs[key] = values

    if scale is not None:
        ds.attrs["Calibration_scale"] = scale
    else:
        if species.upper() in scales.keys():
            ds.attrs["Calibration_scale"] = scales[species.upper()]
        else:
            ds.attrs["Calibration_scale"] = "unknown"

    ds.attrs["species"] = species_out

    # Species-specific attributes
    #############################################

    # Long name
    if species.upper() in species_long.keys():
        sp_long = "mole_fraction_of_" + species_long[species.upper()] + "_in_air"
    else:
        sp_long = "mole_fraction_of_" + species_out + "_in_air"
    
    ancillary_variables = ""
    
    for key in ds.keys():
        
        if species_out in key:
            
            # Standard name attribute
            ds[key].attrs["standard_name"]=key.replace(species_out, sp_long)
            ds[key].attrs["long_name"]=key.replace(species_out, sp_long)
    
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
            
            ancillary_variables += " " + key

    ds[species_out].attrs["ancilliary_variables"] = ancillary_variables.strip()

    # Add quality flag attributes
    ##################################

    flag_key = [key for key in ds.keys() if "_status_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = unflagged, 1 = flagged",
                              "standard_name":
                              ds[species_out].attrs["standard_name"] + "_status_flag"}
    
    # Set time encoding
    #########################################

    first_year = str(ds.time.to_pandas().index.to_pydatetime()[0].year)
    ds.time.encoding = {"units": "seconds since " + \
                        first_year + "-01-01 00:00"}
        
    return ds

def output_filename(output_directory,
                    network,
                    instrument,
                    site,
                    year,
                    species,
                    inlet):
    
    return join(output_directory,
                network + "-" + \
                instrument + "_" + \
                site + "_" + \
                year + "0101_" + \
                species + "-" + \
                inlet + ".nc")

# ICOS
########################################################

def icos_data_read(data_file, species):

    print("Reading " + data_file)

    df =  pd.read_csv(data_file,
                      sep = r"\s+")
    
    # Sometimes header appears in middle of file
    if not "int" in str(df.Year.dtype):
        df = df[df.Year != "Year"]

    df[[species.lower(), "Stdev"]] = df[[species.lower(), "Stdev"]].astype(float)

    df = df[df[species.lower()] >= 0.]
    
    df.reset_index(inplace = True)

    df[["Year", "Month", "Day", "Hour", "Minute"]] = \
        df[["Year", "Month", "Day", "Hour", "Minute"]].astype(int)
    
    time = []
    for y, m, d, h, mi in zip(df.Year, df.Month, df.Day, df.Hour, df.Minute):
        time.append(dt(y, m, d, h, mi))

    df.index = time
    
    # Remove unwanted columns
    df.drop(["Year", "Month", "Day", "Hour", "Minute", "index",
             "Site", "Date", 'Flag', 'QualityId',
             u'LastModifiedYear', u'LastModifiedMonth', u'LastModifiedDay',
             u'LastModifiedHour', u'LastModifiedMinute', u'AutoDescriptiveFlag',
             u'ManualDescriptiveFlag'],
            1,
            inplace = True)

    # Remove duplicate indices
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')
    
    # Rename columns
    df.rename(columns = {"Stdev": species.upper() + "_variability",
                         "NbPoints": species.upper() + "_number_of_observations"},
               inplace = True)

    df.index.name = "time"

    # Convert to Dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())

    return ds


def icos(site, species, network):
    
    params_icos = params["ICOS"]
    site_string = params_icos[site]["gcwerks_site_name"]
    inlets = params_icos[site]["inlets"].keys()
    instruments = params_icos[site]["instruments"]
    
    for instrument in instruments:
        
        for inlet in inlets:

            data_directory = params_icos["directory"].replace("%site", site_string)
            out_directory = params_icos["directory_output"]
            
            # Find data file
            data_file_search = join(data_directory, site.lower() + "." + \
                                    species.lower() + ".1minute." + inlet + ".dat")
            data_file = glob.glob(data_file_search)
            
            if len(data_file) == 0:
                print("Can't find file " + data_file_search)
                return None
            else:
                data_file = data_file[0]
            
            # Create Pandas dataframe
            ds = icos_data_read(data_file, species)
            
            # Sort out attributes
            global_attributes = params_icos[site.upper()]["global_attributes"]
            ds = attributes(ds, species, global_attributes = global_attributes)
            
            # Write file
            nc_filename = output_filename(out_directory,
                                          network,
                                          instrument,
                                          site.upper(),
                                          str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                          ds.species,
                                          params_icos[site]["inlets"][inlet])
            
            ds.to_netcdf(nc_filename)
            
            print("Written " + nc_filename)


# GC FUNCTIONS
###############################################################

def gc_data_read(dotC_file, scale = {}, units = {}):
    
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

    species = []

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
                    area_height_flag.append(0)  # Height

            df = df.rename(columns = {key: df.keys()[i-1] + "_flag"})
            df[df.keys()[i-1] + "_status_flag"] = quality_flag
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
                            parse_dates = True)

    return precision, precision_species


def gc(site, instrument, network):
    """
    Process GC data per site and instrument
    Instruments can be:
        "GCMD": GC multi-detector (output will be labeled GC-FID or GC-ECD)
        "GCMS": GC ADS (output GC-ADS)
        "medusa": GC medusa (output GC-MEDUSA)
    
    Network can be:
        "AGAGE"
        "DECC"
        "GAUGE"
    """

    site_gcwerks = params["GC"][site]["gcwerks_site_name"]
    instrument_gcwerks = params["GC"]["instruments"][instrument]
    if instrument_gcwerks != "":
        instrument_gcwerks = "-" + instrument_gcwerks
    
    data_folder = params["GC"]["directory"][instrument]
    data_files = sorted(glob.glob(join(data_folder,
                                  site_gcwerks + \
                                  instrument_gcwerks + ".??.C")))[0:2]
    
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
            df[sp + "_repeatability"] = precision[precision_index].\
                                            astype(float).\
                                            reindex_like(df, "pad")
        
        dfs.append(df)
    
    # Concatenate
    dfs = pd.concat(dfs).sort_index()
    
    # Drop duplicate indices
    dfs = dfs.reset_index().drop_duplicates(subset='index').set_index('index')

    # Label time index
    dfs.index.name = "time"

    # Convert to xray dataset
    ds = xray.Dataset.from_dataframe(dfs)

    # Get species from scale dictionary
    species = scale.keys()
    
    # Get inlets
    inlets = set(ds["Inlet"].values)
    inlets = [inlet for inlet in inlets if inlet[-1] is "m"]

    for sp in species:
        
        for inlet in inlets:        
    
            ds_sp = ds.where(ds.Inlet == inlet)[[sp,
                                                 sp + "_repeatability",
                                                 sp + "_status_flag"]]

            # Drop NaNs
            ds_sp = ds_sp.dropna("time")

            # Sort out attributes
            global_attributes = params["GC"][site.upper()]["global_attributes"]
            ds_sp = attributes(ds_sp, sp,
                               global_attributes = global_attributes,
                               units = units[sp],
                               scale = scale[sp])

            # Get instrument name for output
            if sp in params["GC"]["instruments_out"][instrument]:
                instrument_out = params["GC"]["instruments_out"][instrument][sp]
            else:
                instrument_out = params["GC"]["instruments_out"][instrument]["else"]

            # Write file
            nc_filename = output_filename(params["GC"]["directory_output"],
                                          network,
                                          instrument_out,
                                          site.upper(),
                                          str(ds_sp.time.to_pandas().index.to_pydatetime()[0].year),
                                          ds_sp.species,
                                          inlet)

            ds_sp.to_netcdf(nc_filename)
            print("Written " + nc_filename)
            


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


def crds(site, network):
    """
    Process CRDS data
    
    site : Three letter site code
    network : Network string only for output 
    """
    params_crds = params["CRDS"]

    site_string = params_crds[site]["gcwerks_site_name"]
    inlets = params_crds[site]["inlets"]
    
    for inlet in inlets:
    
        data_file = glob.glob(join(params_crds["directory"].replace("%site",
                                                                    site_string),
                                   site.lower() + \
                                   ".picarro.1minute." + inlet + ".dat"))[0]
        
        # Create Pandas dataframe
        ds, species = crds_data_read(data_file)
        
        # Write netCDF file for each species
        for sp in species:
                        
            # Species-specific dataset
            ds_sp = ds[[sp, sp + "_variability", sp + "_number_of_observations"]]
            ds_sp.dropna("time")
            
            global_attributes = params_crds[site]["global_attributes"]

            ds_sp = attributes(ds_sp, sp,
                               global_attributes = global_attributes,
                               scale = scales[sp])
            
            # Write file
            nc_filename = output_filename(params["CRDS"]["directory_output"],
                                          network,
                                          "CRDS",
                                          site.upper(),
                                          str(ds_sp.time.to_pandas().index.to_pydatetime()[0].year),
                                          ds_sp.species,
                                          inlet)

            ds_sp.to_netcdf(nc_filename)
            print("Written " + nc_filename)
