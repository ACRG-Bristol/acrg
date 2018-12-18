#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 14:02:12 2018

@author: chxmr
"""

import json
import pandas as pd
import glob
from os import getenv
from os.path import join
import xarray as xr
from acrg_obs.utils import attributes, output_filename


# Site info file
acrg_path = getenv("ACRG_PATH")
data_path = getenv("DATA_PATH")
site_info_file = join(acrg_path, "acrg_site_info.json")
with open(site_info_file) as sf:
    site_params = json.load(sf)

# Set default obs folder
obs_directory = join(data_path, "obs_2018/")


def ufrank(site, species):
    '''
    Process Goethe Frankfurt University data files for Taunus

    Inputs are site code and species

    Assumes that file names start with a date, routine will pick the latest one
    '''

    params = {
        "directory" : "/data/shared/obs_raw/UFrank/",
        "directory_output" : "/dagage2/agage/temp/DECC/taunus_obs/",
        "scale": {
            "PFC-318": "SIO-14",
            "SO2F2": "SIO-07",
            "HFC-32": "SIO-07",
            "HFC-125": "UB-98",
            "HFC-134a": "SIO-05",
            "HFC-143a": "SIO-07",
            "HFC-152a": "SIO-05",
            "HFC-227ea": "Empa-2005",
            "HFC-236fa": "Empa-2009-p",
            "HFC-245fa": "Empa-2005",
            "HCFC-22": "SIO-05",
            "HCFC-141b": "SIO-05",
            "HCFC-142b": "SIO-05",
            "CFC11": "SIO-05",
            "CFC12": "SIO-05",
            "CFC113": "SIO-05",
            "CFC114": "SIO-05",
            "CFC115": "SIO-05",
            "H-1301": "SIO-05",
            "H-1211": "SIO-05",
            "CH3Cl": "SIO-05",
            "CH3Br": "SIO-05",
            "CH3I": "NOAA-Dec09",
            "CH2Cl2": "UB-98",
            "CHBr3": "NOAA-Dec09",
            "CH3CCl3": "SIO-05",
            "PCE": "NOAA-2003B",
            "COS": "NOAA-SIO-p1"},
        "TAU" : {
            "ufrank_name": "tau",
            "instrument": "Medusa",
            "inlet": "8m",
            "global_attributes": {
                "data_owner": "Tanja Schuck",
                "data_owner_email": "schuck@iau.uni-frankfurt.de"
            }
            }
    }


    ufrank_site = params[site]["ufrank_name"]

    fnames = sorted(glob.glob(join(params["directory"],
                            "*_" + site.lower() + "_" + \
                            species.lower() + \
                            "*_ufrank.csv")))

    #Pick most recent file
    fname = fnames[-1]
    
    print("Reading " + fname + "... this can take a while")
    
    if params[site]["instrument"] == "Medusa":
        date_col = "date"
    else:
        date_col = ufrank_site + "_date"
        
    df = pd.read_csv(fname,
                     parse_dates = [date_col],
                     index_col = [date_col],
                     dayfirst = True).sort_index()
    
    rename_dict_all = {ufrank_site + "_data_obs_scaled": species.upper(),
                       ufrank_site + "_" + species.upper(): species.upper(),
                       ufrank_site + "_SD": species.upper() + " repeatability"}

    rename_dict = {}
    for key in rename_dict_all.keys():
        if key in df.keys():
            rename_dict[key] = rename_dict_all[key]
                   
    df.rename(columns = rename_dict,
              inplace = True)

    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"
    
    ds = xr.Dataset.from_dataframe(df.sort_index())
    
    global_attributes = params[site]["global_attributes"]
    global_attributes["inlet_magl"] = params[site]["inlet"]
    
    ds = attributes(ds,
                    species.upper(),
                    site.upper(),
                    global_attributes = global_attributes,
                    scale = params["scale"][species.upper()],
                    sampling_period = 60,
                    units = "ppt")
    
    # Write file
    nc_filename = output_filename(params["directory_output"],
                                  "UFRANK",
                                  params[site]["instrument"],
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  None)
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)

