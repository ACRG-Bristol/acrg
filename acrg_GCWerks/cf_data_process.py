# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 14:08:07 2015

@author: chxmr
"""

import numpy as np
import pandas as pd
from os.path import join, split
from datetime import datetime as dt
from datetime import timedelta as td
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

site_info_file = join(acrg_path, "acrg_site_info.json")
with open(site_info_file) as sf:
    site_params = json.load(sf)

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
                "CO": "carbon_monoxide"}

# Translate header strings
crds_header_string_interpret = {"C": "",
                                "stdev": "_variability",
                                "N": "_number_of_observations"}

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

def attributes(ds, species, site, global_attributes = None,
               units = None, scale = None,
               integration_period = None):
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
        
    # Add some site attributes
    global_attributes_site = site_info_attributes(site.upper())
    if global_attributes_site is not None:
        for key, values in global_attributes_site.iteritems():
            ds.attrs[key] = values
    
    # Add calibration scale
    if scale is not None:
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

    # Add integration flag attributes
    ##################################

    flag_key = [key for key in ds.keys() if "_integration_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = area, 1 = height",
                              "standard_name":
                              ds[species_out].attrs["standard_name"] + "_integration_flag",
                              "comment":
                              "GC peak integration method (by height or by area). " + 
                              "Does not indicate data quality"}
    
    # Set time encoding
    #########################################

    first_year = str(ds.time.to_pandas().index.to_pydatetime()[0].year)
    ds.time.encoding = {"units": "seconds since " + \
                        first_year + "-01-01 00:00:00"}
    ds.time.attrs["label"] = "left"
    ds.time.attrs["comment"] = "Time stamp corresponds to beginning of integration period. " + \
                               "Time since midnight UTC of reference date."
    if integration_period:
        ds.time.attrs["period"] = integration_period
    
    return ds

def output_filename(output_directory,
                    network,
                    instrument,
                    site,
                    year,
                    species,
                    inlet):
    
    return join(output_directory,
                network + "/" + \
                network + "-" + \
                instrument + "_" + \
                site + "_" + \
                year + "0101_" + \
                species + "-" + \
                inlet + ".nc")


# UCAM
########################################################

def ucam(site, species):
    '''
    Process University of Cambridge data files

    Inputs are site code and species

    Assumes that file names start with a date, routine will pick the latest one
    '''

    params_ucam = params["UCAM"]
    ucam_site = params_ucam[site]["ucam_name"]

    fnames = glob.glob(join(params_ucam["directory"],
                            "*_" + site.lower() + "_" + \
                            species.lower() + \
                            "*_ucam.csv"))

    #Pick most recent file
    fname = fnames[-1]
    
    print("Reading " + fname + "... this can take a while")
    
    if params_ucam[site]["instrument"] == "CRDS":
        date_col = ucam_site + "_pic_date"
    else:
        date_col = ucam_site + "_date"
        
    df = pd.read_csv(fname,
                     parse_dates = [date_col],
                     index_col = [date_col]).sort_index()
    
    rename_dict_all = {ucam_site + "_data_obs_scaled": species.upper(),
                       ucam_site + "_obs_repeatability": species.upper() + "_repeatability",
                       ucam_site + "_cal_uncertainty": species.upper() + "_calibration_uncertainty",
                       ucam_site + "_pic_" + species.upper(): species.upper(),
                       ucam_site + "_pic_SD": species.upper() + "_repeatability"}

    rename_dict = {}
    for key in rename_dict_all.keys():
        if key in df.keys():
            rename_dict[key] = rename_dict_all[key]
                   
    df.rename(columns = rename_dict,
              inplace = True)

    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"
    
    ds = xray.Dataset.from_dataframe(df.sort_index())
    
    global_attributes = params_ucam[site]["global_attributes"]
    
    ds = attributes(ds,
                    species.upper(),
                    site.upper(),
                    global_attributes=global_attributes)
    
    # Write file
    nc_filename = output_filename(params_ucam["directory_output"],
                                  "UCAM",
                                  params_ucam[site]["instrument"],
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  params_ucam[site]["inlet"])
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)
    

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


def icos(site, network = "ICOS"):
    
    # Get directories and site strings
    params_icos = params["ICOS"]
    site_string = params_icos[site]["gcwerks_site_name"]
    data_directory = params_icos["directory"].replace("%site", site_string)
    out_directory = params_icos["directory_output"]

    # Search for species and inlets from file names
    data_file_search = join(data_directory, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    data_file_names = [split(f)[1] for f in data_files]
    species_and_inlet = [(f.split(".")[1], f.split(".")[-2]) \
                         for f in data_file_names]
    
    for i, (species, inlet) in enumerate(species_and_inlet):
        
        # Create Pandas dataframe
        ds = icos_data_read(data_files[i], species.upper())
        
        # Sort out attributes
        global_attributes = params_icos[site.upper()]["global_attributes"]
        global_attributes["inlet_height_magl"] = float(inlet[:-1])

        ds = attributes(ds,
                        species.upper(),
                        site.upper(),
                        global_attributes = global_attributes,
                        integration_period = "1 minute")

        # Write file
        nc_filename = output_filename(out_directory,
                                      network,
                                      "CRDS",
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      params_icos[site]["inlet_rename"][inlet])
        
        ds.to_netcdf(nc_filename)
        
        print("Written " + nc_filename)


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
            df[df.keys()[i-1] + "_status_flag"] = quality_flag
            df[df.keys()[i-1] + "_integration_flag"] = area_height_flag
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
    
    # Rename index column
    precision.index.names = ["index"]

    # Drop duplicates
    precision = precision.reset_index().drop_duplicates(subset='index').set_index('index')

    return precision, precision_species


def gc(site, instrument, network):
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
    if instrument_gcwerks != "":
        instrument_gcwerks = "-" + instrument_gcwerks
    
    data_folder = params["GC"]["directory"][instrument]

    # Search string
    search_string = join(data_folder,
                         site_gcwerks + \
                         instrument_gcwerks + ".??.C")
    # Search string with -md attached
    search_string_md = join(data_folder,
                            site_gcwerks + \
                            instrument_gcwerks + "-md.??.C")

    # Search for data files
    # First try to find MD file called site.YY.C
    data_files = sorted(glob.glob(search_string))
    if instrument == "GCMD":
        if len(data_files) == 0:
            # Else try to look for site-md.YY.C
            data_files = sorted(glob.glob(search_string_md))
    
    # Error if can't find files
    if len(data_files) == 0.:
        print("ERROR: can't find any files: " + search_string + " or " + \
              search_string_md)
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
            df[sp + "_repeatability"] = precision[precision_index].\
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
        
    for sp in species:

        global_attributes = params["GC"][site.upper()]["global_attributes"]
        global_attributes["comment"] = params["GC"]["comment"][instrument]

        for inlet in inlets:        
            
            print("Processing " + sp + ", " + inlet + "...")
            
            if inlet == "any":
                ds_sp = ds[[sp,
                            sp + "_repeatability",
                            sp + "_status_flag",
                            sp + "_integration_flag"]]
                inlet_label = params["GC"][site.upper()]["inlet_label"][0]
                global_attributes["inlet_height_magl"] = ", ".join(set(ds["Inlet"].values))
                
            else:
                ds_sp = ds.where(ds.Inlet == inlet)[[sp,
                                                     sp + "_repeatability",
                                                     sp + "_status_flag",
                                                     sp + "_integration_flag"]]
                global_attributes["inlet_height_magl"] = float(inlet[:-1])
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
                                   integration_period = params["GC"]["integration_period"][instrument])
    
                # Get instrument name for output
                if sp.upper() in params["GC"]["instruments_out"][instrument]:
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
                                              inlet_label)
                                              
                ds_sp.to_netcdf(nc_filename)
                print("... written " + nc_filename)
            


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
    data_directory = params_crds["directory"].replace("%site", site_string)

    # Search for species and inlets from file names
    data_file_search = join(data_directory, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    inlets = [f.split(".")[-2] for f in data_files]
    
    for i, inlet in enumerate(inlets):
        
        # Create Pandas dataframe
        ds, species = crds_data_read(data_files[i])
        
        # Write netCDF file for each species
        for sp in species:
            
            # Species-specific dataset
            ds_sp = ds[[sp,
                        sp + "_variability",
                        sp + "_number_of_observations"]]
            ds_sp.dropna("time")
            
            global_attributes = params_crds[site]["global_attributes"]
            global_attributes["inlet_height_magl"] = float(inlet[0:-1])
            global_attributes["comment"] = params_crds["comment"]
    
            ds_sp = attributes(ds_sp, sp, site.upper(),
                               global_attributes = global_attributes,
                               scale = scales[sp],
                               integration_period="1 minute")
            
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


if __name__ == "__main__":

    # ICOS
    icos("TTA")
    icos("MHD", network = "LSCE")

    # GAUGE CRDS data
    crds("HFD", "GAUGE")
    crds("BSD", "GAUGE")

    # GAUGE GC data
    gc("BSD", "GCMD", "GAUGE")
    gc("HFD", "GCMD", "GAUGE")

    # DECC CRDS data
    crds("TTA", "DECC")
    crds("RGL", "DECC")
    crds("TAC", "DECC")

    # DECC GC data    
    gc("TAC", "GCMD", "DECC")
    gc("RGL", "GCMD", "DECC")
    gc("BSD", "GCMD", "DECC")

    # DECC Medusa
    gc("TAC", "medusa", "DECC")

    # AGAGE GC data    
    gc("CGO", "GCMD", "AGAGE")
    gc("MHD", "GCMD", "AGAGE")
    gc("RPB", "GCMD", "AGAGE")
    gc("SMO", "GCMD", "AGAGE")
    gc("THD", "GCMD", "AGAGE")

    # AGAGE Medusa
    gc("MHD", "medusa", "AGAGE")
    gc("CGO", "medusa", "AGAGE")
    gc("GSN", "medusa", "AGAGE")
    gc("SDZ", "medusa", "AGAGE")
    gc("THD", "medusa", "AGAGE")
    gc("RPB", "medusa", "AGAGE")
    gc("SMO", "medusa", "AGAGE")
    gc("SIO", "medusa", "AGAGE")
    
    # University of Cambridge GC
    ucam("WAO", "CH4")
    ucam("HAD", "CH4")
    ucam("TIL", "CH4")
    