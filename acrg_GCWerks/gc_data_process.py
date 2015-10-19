# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:42:24 2015

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
unit_strings = {"ppm": "1e-6",
                "ppb": "1e-9",
                "ppt": "1e-12",
                "--": "1."}

# Output instrument strings
instrument_strings = {"": "GCMD"}


def dotC_read(dotC_file, scale = {}, units = {}):
    
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
            df = df.rename(columns = {key: df.keys()[i-1] + "_flag"})
            scale[df.keys()[i-1]] = header[i-1][0]
            units[df.keys()[i-1]] = header[i-1][1]
            species.append(df.keys()[i-1])

    return df, species, units, scale


def precisions_read(precisions_file):

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


def run(site):

    for instrument in params[site]["instruments"]:
    
        data_files = sorted(glob.glob(join(params["agage_directory"],
                                      params[site]["gcwerks_site_name"] + \
                                      instrument + \
                                      ".??.C")))
        precision_files = [data_file[0:-2] + ".precisions.C" \
                            for data_file in data_files]
        
        dfs = []
        scale = {}
        units = {}

        for fi, data_file in enumerate(data_files):
        
            print("Reading " + data_file)
        
            # Get observations
            df, species, units, scale = dotC_read(data_file,
                                                  scale = scale,
                                                  units = units)
        
            # Get precision
            precision, precision_species = precisions_read(precision_files[fi])

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
        inlets = [inlet for inlet in inlets if "pump" not in inlet and \
                                               "tank" not in inlet and\
                                               "Sib" not in inlet]
        
        for sp in species:
            for inlet in inlets:        
        
                ds_sp = ds.where(ds.Inlet == inlet)[[sp, sp + "_repeatability"]]
                
                # Drop NaNs
                ds_sp = ds_sp.dropna("time")

                # Set units and attributes
                ds_sp.time.encoding = {"units": "seconds since 1970-01-01 00:00"}
                ds_sp[sp].attrs = {"units": unit_strings[units[sp]]}
                
                ds_sp.attrs = {"Contact": \
                               params[site]["species_defaults"]["contact"],
                               "Contact email": \
                               params[site]["species_defaults"]["contact_email"],
                               "Calibration scale": scale[sp]}
                
                ds_sp.to_netcdf(join(params["output_directory"],
                                params[site]["species_defaults"]["network"] + "-" + \
                                instrument_strings[instrument] + "_" + \
                                site + "_" + inlet + "agl-" + 
                                sp.replace("-", "") + ".nc"))
