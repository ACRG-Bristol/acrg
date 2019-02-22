#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 14:02:12 2018

@author: chxmr
"""
from __future__ import print_function

import json
import pandas as pd
import glob
from os import getenv
from os.path import join, split
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


def ufrank(site):
    '''
    Process Goethe Frankfurt University data files for Taunus

    Inputs are site code and species

    Assumes that file names start with a date, routine will pick the latest one
    '''

    params = {
        "directory" : "/data/shared/obs_raw/UFrank/",
        "directory_output" : "/data/shared/obs_2018",
        "scale": {
            "C4F8": "SIO-14",
            "SO2F2": "SIO-07",
            "HFC32": "SIO-07",
            "HFC125": "UB-98",
            "HFC134A": "SIO-05",
            "HFC143A": "SIO-07",
            "HFC152A": "SIO-05",
            "HFC227EA": "Empa-2005",
            "HFC236FA": "Empa-2009-p",
            "HFC245FA": "Empa-2005",
            "HCFC22": "SIO-05",
            "HCFC141B": "SIO-05",
            "HCFC142B": "SIO-05",
            "CFC11": "SIO-05",
            "CFC12": "SIO-05",
            "CFC113": "SIO-05",
            "CFC114": "SIO-05",
            "CFC115": "SIO-05",
            "HALON1301": "SIO-05",
            "HALON1211": "SIO-05",
            "CH3CL": "SIO-05",
            "CH3BR": "SIO-05",
            "CH3I": "NOAA-Dec09",
            "CH2CL2": "UB-98",
            "CHBR3": "NOAA-Dec09",
            "CH3CCL3": "SIO-05",
            "C2CL4": "NOAA-2003B",
            "COS": "NOAA-SIO-p1"},
        "TAU" : {
            "ufrank_name": "tau",
            "instrument": "GCMS",
            "inlet": "8m",
            "global_attributes": {
                "data_owner": "Tanja Schuck",
                "data_owner_email": "schuck@iau.uni-frankfurt.de"
            }
            }
    }

    # Find name of species
    fnames = glob.glob(join(params["directory"], "*%s*_ufrank.csv" % site.lower()))
    species_ufrank = set([split(f)[-1].split("_")[-2] for f in fnames])
    
    for species in species_ufrank:

        #Pick most recent file
        fname = sorted(glob.glob(join(params["directory"], "*%s*%s*_ufrank.csv" % 
                                      (site.lower(), species))))[-1]

        print("Reading %s ..." %fname)
        
        date_col = "date"
        
        df = pd.read_csv(fname,
                         parse_dates = [date_col],
                         index_col = [date_col]).sort_index()

        df.rename(columns = {species + "_SD": species + " repeatability"},
                  inplace = True)
        
        # Drop duplicates
        df.index.name = "index"
        df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
        df.index.name = "time"
        
        ds = xr.Dataset.from_dataframe(df.sort_index())
        
        global_attributes = params[site]["global_attributes"]
        global_attributes["inlet_magl"] = params[site]["inlet"]
        
        ds = attributes(ds,
                        species,
                        site,
                        global_attributes = global_attributes,
                        scale = params["scale"][species.upper()],
                        sampling_period = 60,
                        units = "ppt")
        
        # Write file
        nc_filename = output_filename(params["directory_output"],
                                      "UFRANK",
                                      params[site]["instrument"],
                                      site.upper(),
                                      ds.time.to_pandas().index.to_pydatetime()[0],
                                      ds.species,
                                      None)
        
        print(" ... writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)

