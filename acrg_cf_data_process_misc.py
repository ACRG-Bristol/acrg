# -*- coding: utf-8 -*-
"""
Data processing routines for various networks/stations

Created on Thu Jan 14 10:47:45 2016

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
from acrg_GCWerks.cf_data_process import attributes, output_filename
import calendar
import re


# Site info file
acrg_path = getenv("ACRG_PATH")
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

# UCAM
########################################################

def ucam(site, species):
    '''
    Process University of Cambridge data files

    Inputs are site code and species

    Assumes that file names start with a date, routine will pick the latest one
    '''

    params = {
        "directory" : "/data/shared/obs_raw/UCAM/",
        "directory_output" : "/data/shared/obs/",
        "scale": {
            "CH4": "WMOx2004a",
            "CO2": "WMOx2007"},
        "EHL" : {
            "ucam_name": "Eh",
            "instrument": "CRDS",
            "inlet": "25m",
            "global_attributes": {
                "data_owner": "Neil Harris",
                "data_owner_email": "nrh1000@cam.ac.uk"
            }
        },
        "HAD" : {
            "ucam_name": "Had",
            "instrument": "GC-FID",
            "inlet": "25m",
            "global_attributes": {
                "data_owner": "Neil Harris",
                "data_owner_email": "nrh1000@cam.ac.uk"
            }
        },
        "TIL" : {
            "ucam_name": "Til",
            "instrument": "GC-FID",
            "inlet": "27m",
            "global_attributes": {
                "data_owner": "Neil Harris",
                "data_owner_email": "nrh1000@cam.ac.uk"
            }
        },
        "WAO" : {
            "ucam_name": "Wao",
            "instrument": "GC-FID",
            "inlet": "25m",
            "global_attributes": {
                "data_owner": "Neil Harris",
                "data_owner_email": "nrh1000@cam.ac.uk"
            }
        }
    }


    ucam_site = params[site]["ucam_name"]

    fnames = sorted(glob.glob(join(params["directory"],
                            "*_" + site.lower() + "_" + \
                            species.lower() + \
                            "*_ucam.csv")))

    #Pick most recent file
    fname = fnames[-1]
    
    print("Reading " + fname + "... this can take a while")
    
    if params[site]["instrument"] == "CRDS":
        date_col = ucam_site + "_pic_date"
    else:
        date_col = ucam_site + "_date"
        
    df = pd.read_csv(fname,
                     parse_dates = [date_col],
                     index_col = [date_col],
                     dayfirst = True).sort_index()
    
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
    
    global_attributes = params[site]["global_attributes"]
    
    ds = attributes(ds,
                    species.upper(),
                    site.upper(),
                    global_attributes = global_attributes,
                    scale = params["scale"][species.upper()])
    
    # Write file
    nc_filename = output_filename(params["directory_output"],
                                  "UCAM",
                                  params[site]["instrument"],
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  params[site]["inlet"])
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)

def ucam_run():

    # University of Cambridge GC
    ucam("WAO", "CH4")
    ucam("HAD", "CH4")
    ucam("TIL", "CH4")
    ucam("EHL", "CH4")
    ucam("EHL", "CO2")



def cbw():
    '''
    Process Cabauw

    '''

    def cbw_write(df, instrument, inlet):
        
        ds = xray.Dataset.from_dataframe(df.sort_index())

        ds = attributes(ds,
                        species.upper(),
                        site.upper(),
                        global_attributes=params["global_attributes"])
        
        # Write file
        nc_filename = output_filename(params["directory_output"],
                                      "ECN",
                                      instrument,
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      inlet)
    
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)


    params = {
        "directory" : "/data/shared/obs_raw/ECN/",
        "directory_output" : "/data/shared/obs/",
        "inlets": ["20m", "60m", "120m", "200m"],
        "global_attributes" : {
            "contact": "Danielle van Dinther"
        }
    }


    site = "CBW"
    species = "CH4"
    
    fnames = sorted(glob.glob(join(params["directory"],
                                   "*" + site.upper() + "*.csv")))
    print("Reading...")
    print(fnames)

    rename_dict = {"CH4_std_20m": "STD_CH4_20m",
                   "CH4_std_60m": "STD_CH4_60m",
                   "CH4_std_120m": "STD_CH4_120m",
                   "CH4_std_200m": "STD_CH4_200m"}

    dfs = []
    for fname in fnames:
        dff = pd.read_csv(fname, sep = ";",
                          parse_dates = ["Time"],
                          index_col = ["Time"],
                          dayfirst = True)
        dff.rename(columns = rename_dict, inplace = True)
        dfs.append(dff)

    # Combine
    if len(dfs) > 1:
        df = dfs[0].combine_first(dfs[1])
    else:
        df = dfs[0].copy()
        del dfs
    
    # remove duplicates
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"

    # output for each inlet
    for inlet in params["inlets"]:
        
        df_inlet = df[[species + "_" + inlet,
                       "STD_" + species + "_" + inlet]]

        df_inlet.rename(columns = {species + "_" + inlet : species,
                                   "STD_" + species + "_" + inlet : species + "_variability"},
                        inplace = True)
        
        df_inlet = df_inlet[np.isfinite(df_inlet[species])]

        # GC data
        df_inlet_gc = df_inlet[np.isnan(df_inlet[species + "_variability"])]
        df_inlet_gc.rename(columns = {species + "_variability":
                                      species + "_repeatability"},
                           inplace = True)
        
        df_inlet_gc["1992":"1997"][species + "_repeatability"] = 10.
        df_inlet_gc["2000":"2004-10"][species + "_repeatability"] = 3.
        df_inlet_gc["2004-11":"2009-1"][species + "_repeatability"] = 2.
        df_inlet_gc["2009-2":][species + "_repeatability"] = 1.
        
        cbw_write(df_inlet_gc, "GC-FID", inlet)

        # CRDS data
        df_inlet_crds = df_inlet[np.isfinite(df_inlet[species + "_variability"])]
        cbw_write(df_inlet_crds, "CRDS", inlet)


def ec(site):
    """
    Process Environment Canada data
    """
    params = {
        "directory" : "/data/shared/obs_raw/EC/",
        "directory_output" : "/data/shared/obs/",
        "global_attributes" : {
            "contact": "Doug Worthy",
            "averaging": "hourly averages of GC data"
        }
    }
    
    def ec_correct_time(df):
        time = []
        for t in df.index:
            year, month, day, hour, minute = map(int, re.split(" |:|-", t))
            if hour == 24:
                hour = 0
                if day == calendar.monthrange(year, month)[1]:
                    day = 1
                    if month == 12:
                        month = 1
                        year = year + 1
                    else:
                        month = month + 1
                else:
                    day = day + 1
            time.append(dt(year,
                           month,
                           day,
                           hour,
                           minute))
        df.index = time
        df.index.name = "time"
        return df
        
    fnames = sorted(glob.glob(join(params["directory"],
                                   "*" + site.upper() + "*/*/*/*.dat")))

    df = []
    
    for fname in fnames:
        
        header_count = 0
        with open(fname, "r") as f:
            for line in f:
                if line[0] == "C":
                    header_count+=1
                    last_line = line
                else:
                    break
        
        header = last_line[5:].split()
        header[2] += "x"
        header[3] += "x"
    
        dff = pd.read_csv(fname, sep = r"\s+",
                         skiprows = header_count,
                         header = None,
                         names = header,
                         parse_dates = {"time": ["DATE", "TIME"]},
                         index_col = "time")
        
        dff = dff[["CH4", "ND", "SD"]]
    
        dff.rename(columns = {"ND": "CH4_number_of_observations",
                              "SD": "CH4_variability"},
                   inplace = True)

        dff = ec_correct_time(dff)

        df.append(dff)

    df = pd.concat(df)
    
    # remove duplicates
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"
 
    # Convert to xray dataset
    ds = xray.Dataset.from_dataframe(df)

    # Add attributes
    ds = attributes(ds,
                    "CH4",
                    site.upper(),
                    global_attributes = params["global_attributes"],
                    scale = "NOAA-2004")
    
    # Write file
    nc_filename = output_filename(params["directory_output"],
                                  "EC",
                                  "GC-FID",
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  site_params[site]["height"][0])
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)

def ec_run():
    
    ec("EGB")
    ec("ESP")
    ec("ETL")
    ec("FSD")
    ec("LLB")
    ec("WSA")
    



# UCAM
########################################################

def uea(site, species):
    '''
    Process University of East Anglia data files

    Inputs are site code and species
    '''

    params = {
        "directory" : "/data/shared/obs_raw/UEA/",
        "directory_output" : "/data/shared/obs/",
        "scale": {
            "CH4": "NOAA2004",
            "CO2": "WMOx2007"},
        "WAO" : {
            "instrument": "GC-FID",
            "inlet": "21m",
            "global_attributes": {
                "data_owner": "Grant Forster",
                "data_owner_email": "g.forster@uea.ac.uk"
            }
        }
    }


    fnames = sorted(glob.glob(join(params["directory"],
                              site.lower() + "*" + \
                              species.lower() + \
                              "*.na")))
    
    data = []
    time = []
    
    # Read files
    ###############################
    
    for fname in fnames:

        f = open(fname)    
        line = f.readline()
        header_lines, dummy = [int(l) for l in line.split()]
        header = []
        for i in range(header_lines-2):
            header.append(f.readline())
        columns = f.readline().split()
        lines = f.readlines()
        f.close()

        time_file = []
        data_file = []
        
        for line in lines:
            
            data_line = line.split()
            
            time_file.append(dt(int(data_line[columns.index("year")]),
                                int(data_line[columns.index("month")]),
                                int(data_line[columns.index("day")]),
                                int(data_line[columns.index("hour")]),
                                int(data_line[columns.index("minute")]),
                                int(data_line[columns.index("second")])))

            data_file.append(data_line)

        data_file = np.array(data_file)

        time.append(time_file)
        data.append(data_file)
    
    time = np.concatenate(time)
    data = np.concatenate(data)


    # Store file in dataframes
    ch4i = columns.index("CH4_ppb")
    ch4_flagi = columns.index("CH4_ppb_flag")
    ch4_df = pd.DataFrame(data = data[:, [ch4i, ch4_flagi]].astype(float),
                 columns = ["CH4", "CH4_status_flag"],
                 index = time)
    ch4_df = ch4_df[ch4_df.CH4_status_flag <= 1]

    ch4i = columns.index("CH4_WT_ppb")
    ch4_flagi = columns.index("ch4_WT_flag")
    ch4_wt_df = pd.DataFrame(data = data[:, [ch4i, ch4_flagi]].astype(float),
                 columns = ["CH4_WT", "CH4_WT_flag"],
                 index = time)
    ch4_wt_df = ch4_wt_df[ch4_wt_df.CH4_WT_flag <= 1]

    ch4i = columns.index("CH4_TT_ppb")
    ch4_flagi = columns.index("ch4_TT_flag")
    ch4_tt_df = pd.DataFrame(data = data[:, [ch4i, ch4_flagi]].astype(float),
                 columns = ["CH4_TT", "CH4_TT_flag"],
                 index = time)
    ch4_tt_df = ch4_tt_df[ch4_tt_df.CH4_TT_flag <= 1]

    # Calculate sigma and append to dataframe
    ch4_df["CH4_repeatability"] = \
        ch4_wt_df.CH4_WT.resample("1D", how = np.std).reindex_like(ch4_df, method = "nearest")


    # Drop any duplicates
    ch4_df.index.name = "index"
    ch4_df = ch4_df.reset_index().drop_duplicates(subset='index').set_index('index')              
    ch4_df.index.name = "time"
    
    # Sort and convert to dataset
    ds = xray.Dataset.from_dataframe(ch4_df.sort_index())
    
    global_attributes = params[site]["global_attributes"]
    global_attributes["inlet_height_magl"] = int(params[site]["inlet"][:-1])
    
    # Add attributes
    ds = attributes(ds,
                    species.upper(),
                    site.upper(),
                    global_attributes = global_attributes,
                    scale = params["scale"][species.upper()],
                    sampling_period = 60)
    
    # Write file
    nc_filename = output_filename(params["directory_output"],
                                  "UEA",
                                  params[site]["instrument"],
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  params[site]["inlet"])
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)
