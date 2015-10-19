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
                 "acrg_GCWerks/data_process_parameters.json")
with open(info_file) as sf:
    params = json.load(sf)

# Output unit strings
unit_strings = {"CO2": "1e-6",
                "CH4": "1e-9",
                "N2O": "1e-9",
                "CO": "1e6"}
                
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
header_string_translate = {"C": "",
                           "stdev": "_variability",
                           "N": "_number_of_observations"}


def data_read(data_file):

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
                          header_string_translate[df_header[i][1]])
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

    return df, species


def run(site):
    
    site_string = params[site]["gcwerks_site_name"] + "-picarro"
    inlets = params[site]["inlets"]
    
    for inlet in inlets:
    
        data_file = glob.glob("/dagage2/agage/" + site_string + \
                              "/results-gcwerks/" + site.lower() + \
                              ".picarro.1minute." + inlet + ".dat")[0]
        
        print("Reading " + data_file)
        
        # Create Pandas dataframe
        df, species = data_read(data_file)
        
        # Remove duplicate indices
        df = df.reset_index().drop_duplicates(subset='index').set_index('index')
        
        # Convert to Dataset
        df.index.name = "time"
        ds = xray.Dataset.from_dataframe(df.sort_index())
        
        
        # Write netCDF file for each species
        for sp in species:
            
            # Long name
            sp_long_name = "mole_fraction_of_" + species_long[sp] + "_in_air"
            
            # Species-specific dataset
            ds_sp = ds[[sp, sp + "_variability", sp + "_number_of_observations"]]
            ds_sp.dropna("time")
            
            # Set time units
            ds_sp.time.encoding = {"units": "seconds since 1970-01-01 00:00"}
            
            # Variable attributes
            ds_sp[sp].attrs = {"units": unit_strings[sp],
                               "standard_name": sp_long_name,
                               "ancilliary_variables": sp_long_name + "_variability " + \
                                                       sp_long_name + "_number_of_observations"}
            ds_sp[sp + "_variability"].attrs = {"units": unit_strings[sp],
                               "standard_name": sp_long_name + "_variability"}    
            ds_sp[sp + "_number_of_observations"].attrs = {"units": "",
                               "standard_name": sp_long_name + "_number_of_observations"}
        
            # Global attributes
            ds_sp.attrs = {"Contact": \
                           params[site]["species_defaults"]["contact"],
                           "Contact email": \
                           params[site]["species_defaults"]["contact_email"],
                           "Calibration scale": scales[sp]}
        
            # Earlist year
            first_year = str(ds_sp.time.to_pandas().index.to_pydatetime()[0].year)
        
            # Write file
            out_filename = join(params["output_directory"],
                                 params[site]["species_defaults"]["network"] + \
                                 "-" + "CRDS_" + \
                                 site + "_" + first_year + \
                                 "0101_" + sp.replace("-", "").lower() + "-" + \
                                 inlet + ".nc")
            ds_sp.to_netcdf(out_filename)

            print("... written " + out_filename)