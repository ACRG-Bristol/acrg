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
import xarray as xray
import json
from os import getenv
from acrg_GCWerks.cf_data_process import attributes, output_filename
import calendar
import re
import dateutil
from acrg_time import convert
import pytz


# Site info file
acrg_path = getenv("ACRG_PATH")
data_path = getenv("DATA_PATH")
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
            "inlet": "21m",
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
    global_attributes["inlet_height_magl"] = float(params[site]["inlet"][:-1])
    
    ds = attributes(ds,
                    species.upper(),
                    site.upper(),
                    global_attributes = global_attributes,
                    scale = params["scale"][species.upper()],
                    sampling_period = 60)
    
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

        global_attributes_inlet = params["global_attributes"].copy()
        global_attributes_inlet["inlet_height_magl"] = inlet[:-1]

        ds = attributes(ds,
                        species.upper(),
                        site.upper(),
                        global_attributes=global_attributes_inlet,
                        sampling_period = 60)
        
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
            "Data_owner": "Danielle van Dinther"
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
                          parse_dates = ["Date"],
                          index_col = ["Date"],
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
        with odffpen(fname, "r") as f:
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


def noaa_montzka_usa(species, network = "NOAA"):
    '''
    Process Steve Montzka's USA files

    '''

    def date_parser(year, month, day, hour, minute, second):
        return dt(year, month, day, hour, minute, second)

    params = {
        "directory" : "/data/shared/obs_raw/NOAA/montzka/",
        "directory_output" : "/data/shared/obs/",
        "scale": {
            "HFC-134a": "NOAA1995"
            },
        "species_noaa": {
            "HFC-134a": "HFC_134a"
            },
        "global_attributes": {
            "data_owner": "Steve Montzka",
            "data_owner_email": "stephen.a.montzka@noaa.gov"
            }
        }

    fnames = glob.glob(join(params["directory"],
                            params["species_noaa"][species] + "*.csv"))
    fname = fnames[0]

    with open(fname) as f:
        header = f.readline().split(",")
        lines = f.readlines()
    
    section_description_start = []
    section_description_end = []
    section_description = []

    header_section = False
    for li, l in enumerate(lines):
        if l[0] == "*" and not header_section:
            header_section = True
            section_description_start.append(li)
        if l[0:2] == "* " and header_section:
            section_description.append(l[1:])
        if l[0] != "*" and header_section:
            section_description_end.append(li)
            header_section = False

    platforms = []
    for desc in section_description:
        if "Daily" in desc and "Ground" in desc:
            platforms.append("-tower")
        elif "Biweekly" in desc and "Aircraft" in desc:
            platforms.append("-aircraft")
        elif "Weekly" in desc and "Ground" in desc:
            platforms.append("")
    
    for pi, platform in enumerate(platforms):

        if pi == len(platforms)-1:
            lastrow = None
        else:
            lastrow = section_description_start[pi+1] - section_description_end[pi]
            
        df = pd.read_csv(fname, skiprows = section_description_end[pi] + 1,
                         names = header, nrows = lastrow,
                         parse_dates = {"time": ["sample_year", "sample_month", "sample_day",
                                        "sample_hour", "sample_minute", "sample_second"]},
                         dtype = {"sample_year": np.int,
                                  "sample_month": np.int,
                                  "sample_day": np.int,
                                  "sample_hour": np.int,
                                  "sample_minute": np.int,
                                  "sample_second": np.int},
                         index_col = "time",
                         date_parser = date_parser,
                         skipinitialspace = True)
    
        df = df[["site", "sample_latitude", "sample_longitude", "sample_altitude",
                 "sample_terrain_height", "mole_fraction"]]
    
        df.rename(columns = {"mole_fraction": species}, inplace = True)
        
        for site in set(df.site.values):
            
            dfs = df[df.site == site]
            dfs.drop(["site"], axis = 1, inplace = True)
    
            dfs.index.name = "index"
            dfs = dfs.reset_index().drop_duplicates(subset='index').set_index('index')              
            dfs.index.name = "time"
            
            global_attributes = params["global_attributes"]
            inlet_heights = set(np.round((dfs["sample_altitude"] - \
                                          dfs["sample_terrain_height"]).values/10.)*10)
    
            if len(inlet_heights) > 4:
    
                global_attributes_site = global_attributes.copy()
                global_attributes_site["inlet_height_magl"] = int(-999)
    
                # Sort and convert to dataset
                ds = xray.Dataset.from_dataframe(dfs.sort_index())
                
                # Add attributes
                ds = attributes(ds,
                                species,
                                site.upper(),
                                global_attributes = global_attributes_site,
                                scale = params["scale"][species],
                                sampling_period = 60,
                                units = "ppt")
                
                # Write file
                nc_filename = output_filename(params["directory_output"],
                                              "NOAA",
                                              "GCMS" + platform,
                                              site.upper(),
                                              str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                              ds.species,
                                              "various")
                
                print("Writing " + nc_filename)
                
                ds.to_netcdf(nc_filename)
    
            else:
                
                inlet_heights_rounded = np.round((dfs["sample_altitude"] - \
                                          dfs["sample_terrain_height"]).values/10.)*10
                
                inlet_str = [str(int(ih)) + "m" for ih in inlet_heights]
                
                for inlet_i, inlet in enumerate(inlet_heights):
    
                    dfs_inlet = dfs[inlet_heights_rounded == inlet]
    
                    # Sort and convert to dataset
                    ds = xray.Dataset.from_dataframe(dfs_inlet.sort_index())
                    
                    global_attributes_site = global_attributes.copy()
                    global_attributes_site["inlet_height_magl"] = int(inlet)
                    
                    # Add attributes
                    ds = attributes(ds,
                                    species,
                                    site.upper(),
                                    global_attributes = global_attributes_site,
                                    scale = params["scale"][species],
                                    sampling_period = 60,
                                    units = "ppt")
                    
                    # Write file
                    nc_filename = output_filename(params["directory_output"],
                                                  "NOAA",
                                                  "GCMS" + platform,
                                                  site.upper(),
                                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                                  ds.species,
                                                  inlet_str[inlet_i])
                    
                    print("Writing " + nc_filename)
                    
                    ds.to_netcdf(nc_filename)
                    

def noaa_ccgg(species):

    def date_parser(year, month, day, hour, minute, second):
        return dt(year, month, day, hour, minute, second)
    
    params = {
        "directory" : "/data/shared/obs_raw/NOAA/CCGG/surface",
        "directory_output" : "/data/shared/obs/",
        "units": {"CH4": "ppb",
                  "C2H6": "ppb",
                  "CO2": "ppm",
                  "CH4C13": "permil"},
        "instrument": {"CH4": "GC-FID",
                       "C2H6": "GC-FID",
                       "CO2": "NDIR",
                       "CH4C13": "IRMS"},
        "scale": {"CH4": "NOAA04",
                  "C2H6": "NOAA12",
                  "CO2": "WMO_X2007",
                  "CH4C13": "NOAA-INSTAAR"},
        "global_attributes": {
            "data_owner": "Ed Dlugokencky",
            "data_owner_email": "ed.dlugokencky@noaa.gov"
            }
        }

    
    fnames=glob.glob(join(params["directory"],
                          species.upper(),
                          "event",
                          '*event.txt'))

    for fname in fnames:
    
        header = []
        with open(fname, "r") as f:
            for line in f:
                if line[0] == "#":
                    header.append(line)
                else:
                    break
    
        columns = header[-1][14:].split()
        
        df = pd.read_csv(fname, skiprows = len(header), sep = r"\s+",
                         names = columns,
                         parse_dates = {"time": ["sample_year", "sample_month", "sample_day",
                                        "sample_hour", "sample_minute", "sample_seconds"]},
                         dtype = {"sample_year": np.int,
                                  "sample_month": np.int,
                                  "sample_day": np.int,
                                  "sample_hour": np.int,
                                  "sample_minute": np.int,
                                  "sample_seconds": np.int},
                         index_col = "time",
                         date_parser = date_parser,
                         skipinitialspace = True)
        
        site = df["sample_site_code"][0]
        
        flag = []
        selection_flag = []
        for flag_str in df.analysis_flag:
            flag.append(flag_str[0] == '.')
            selection_flag.append(int(flag_str[1] != '.'))
            
        df[species + "_status_flag"] = flag
        df[species + "_selection_flag"] = selection_flag
    
        df = df[df[species + "_status_flag"]]
        
        df = df[["sample_latitude", "sample_longitude", "sample_altitude",
                 "analysis_value", "analysis_uncertainty",
                 species + "_selection_flag"]]
        
        df.rename(columns = {"analysis_value": species,
                             "analysis_uncertainty": species + "_repeatability",
                             "sample_longitude": "longitude",
                             "sample_latitude": "latitude",
                             "sample_altitude": "altitude"}, inplace = True)
        
        df.index.name = "index"
        df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
        df.index.name = "time"
        
        global_attributes = params["global_attributes"]
    
        # Sort and convert to dataset
        ds = xray.Dataset.from_dataframe(df.sort_index())
        
        # Add attributes
        ds = attributes(ds,
                        species,
                        site.upper(),
                        global_attributes = global_attributes,
                        scale = params["scale"][species],
                        sampling_period = 60)
    
        # Write file
        nc_filename = output_filename(params["directory_output"],
                                      "NOAA-CCGG",
                                      params["instrument"][species],
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      "various")
        
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)
    

def gla(species):
    
    fnames = glob.glob("/data/shared/obs_raw/GAUGE/*GLA*" + species.lower() + "*10m.nc")
    
    ds = []
    for fname in fnames:
        with xray.open_dataset(fname) as f:
            dsf = f.load()
        if species.lower() in dsf.variables.keys():
            dsf.rename({species.lower(): species.upper()}, inplace = True)
        dsf.rename({"uncertainty": species.upper() + "_variability"}, inplace = True)
        ds.append(dsf)
    
    ds = xray.concat(ds, dim = "time")

    # Change time stamp (currently the middle)
    ds['time'] -= np.median(ds.time.values[1:] - ds.time.values[0:-1])/2.
    ds.attrs["comment"] = "Time stamp used to be middle of averaging period, " + \
                    "now subtracted half the " + \
                    "median to approximate the start"

    # Add attributes
    ds = attributes(ds,
                    species,
                    "GLA",
                    global_attributes = {"inlet_height_magl": 10.},
                    sampling_period = 3*60)

    # Write file
    nc_filename = output_filename("/data/shared/obs/",
                                  "GAUGE",
                                  "FTS",
                                  "GLA",
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  "10m")
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)
    

def wdcgg_read(fname, species,
               permil = False,
               repeatability_column = None):
        
    skip = 0
    with open(fname, 'r') as f:
        for li in range(200):
            line = f.readline()
            if line[0] == "C":
                skip += 1
                header_line = line
            else:
                break
    
    columns = header_line[4:].split()
    columns[2] = columns[2] + "_2"
    columns[3] = columns[3] + "_2"
    
    df = pd.read_csv(fname, skiprows = skip,
                     sep = r"\s+", names = columns,
                     parse_dates = {"time": ["DATE", "TIME"]},
                     index_col = "time")
    
    sp_file = species[:-1].upper() + species[-1].lower()
    
    df = df.rename(columns = {sp_file: species})
    
    output_columns = {species: species}
    if repeatability_column:
        output_columns[repeatability_column] = species.lower() + "_repeatability"
    
    df = df[output_columns.keys()]

    df.rename(columns = output_columns, inplace = True)

    if not permil:
        df = df[df[species] > 0.]
    
    # Drop duplicates
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"

    if type(df.index) != pd.tseries.index.DatetimeIndex:
        # Index is not of type time. This is sometimes because time is 99:99
        keep_row = []
        for i in df.index:
            d, t = i.split(" ")
            if t[0:2] == "99":
                keep_row.append(False)
            else:
                keep_row.append(True)

        df = df[keep_row]
        df.index = pd.to_datetime(df.index)

    return df


def nies_read(network, site,
              global_attributes = {},
              instrument = "",
              assume_repeatability = None):

    params = {
        "directory_output" : "/data/shared/obs/"
        }

    
    directories = glob.glob("/data/shared/obs_raw/" + network + \
                            "/" + site + "/*/*")

    species = []
    fnames = []
    for directory in directories:
        species.append(directory.split("/")[-1])
        fnames.append(glob.glob("/data/shared/obs_raw/" + network + \
                                "/" + site + "/*/" + species[-1] + "/event/*.dat")[0])
                                
    for sp, fname in zip(species, fnames):

        df = wdcgg_read(fname, sp)
        
        if assume_repeatability:
            df[sp + "_repeatability"] = df[sp]*assume_repeatability
            global_attributes["Assumed_repeatability_%"] = int(assume_repeatability*100.)
        
        # Sort and convert to dataset
        ds = xray.Dataset.from_dataframe(df.sort_index())
        
        # Add attributes
        ds = attributes(ds,
                        sp,
                        site.upper(),
                        global_attributes = global_attributes,
                        sampling_period = 60,
                        units = "ppt")

        if assume_repeatability:
            ds[sp.lower() + "_repeatability"].attrs["Comment"] = \
                "NOTE: This is an assumed value. Contact data owner."
    
        # Write file
        nc_filename = output_filename(params["directory_output"],
                                      network,
                                      instrument,
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      "various")
        
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)


def nies_wdcgg():

    global_attributes = {"data_owner": "NIES",
                         "data_owner_email": "lnmukaih@nies.go.jp"
                         }
    
    nies_read("NIES", "HAT",
              global_attributes = global_attributes,
              instrument = "GCMS",
              assume_repeatability = 0.03)
    nies_read("NIES", "COI",
              global_attributes = global_attributes,
              instrument = "GCMS",
              assume_repeatability = 0.03)



def nies(fname, species, site, units = "ppt"):
    '''
    Examples of files that this will process:
        fname = "/data/shared/obs_raw/NIES/HAT/HAT_20170804_PFC-318.xlsx"
        species = "PFC-318"
        site = "HAT"
        
        fname = "/data/shared/obs_raw/NIES/HAT/HAT_20170628_CHCl3.txt"
        species = "CHCl3"
        site = "HAT"
    '''

    global_attributes = {"data_owner": "Takuya Saito",
                         "data_owner_email": "saito.takuya@nies.go.jp"
                         }    
    
    
    if fname.split(".")[1] == "xlsx":
        df = pd.read_excel(fname, parse_dates = [0], index_col = [0])
    else:
        df = pd.read_csv(fname, sep = "\t", parse_dates = [0], index_col = [0])
    
    print("Assuming data is in JST. Check input file. CONVERTING TO UTC.")
    
   # df.index = df.index.tz_localize(pytz.timezone("Japan")).tz_convert(pytz.utc).tz_localize(None) # Previous solution
    df.index = df.index.tz_localize(pytz.timezone("Japan")).tz_convert(None) # Simpler solution

    # Sort
    df.sort_index(inplace = True)

    # Rename columns to species
    df.rename(columns = {df.columns[0]: species}, inplace = True)
    
    # Drop duplicates and rename index
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"

    
    # Add a repeatability column
    df[species + "_repeatability"] = df[species]*0.05

    # Convert to dataset
    ds = xray.Dataset.from_dataframe(df)

    ds[species + "_repeatability"].attrs["Comment"] = \
        "NOTE: This is an assumed value. Contact data owner."
    
    
    # Add attributes
    ds = attributes(ds,
                    species,
                    site,
                    scale = "Check raw data file or contact data owner",
                    global_attributes = global_attributes,
                    sampling_period = 60,
                    units = units)

    # Write file
    #output_path = "/data/shared/obs/"
    output_path = join(data_path,"obs/")
    nc_filename = output_filename(output_path,
                                  "NIES",
                                  "GCMS",
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  "various")
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)



def niwa(site, species):

    global_attributes = {"data_owner": "NIWA",
                         "data_owner_email": "sylvia.nichol@niwa.co.nz gordon.brailsford@niwa.co.nz",
                         "calibration_scale": "V-PDB derived from IAEA value for NBS19"
                         }

    fname = glob.glob("/data/shared/obs_raw/NIWA/" + \
                      site.lower() + "*" + species.lower() + "*.txt")[0]
    
    df = wdcgg_read(fname, species,
                    permil = True,
                    repeatability_column="SD")

    # Sort and convert to dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())

    # Only post-1992 data is good
    ds = ds.loc[dict(time = slice("1992-01-01", "2020-01-01"))]

    # Add attributes
    ds = attributes(ds,
                    species,
                    site.upper(),
                    global_attributes = global_attributes,
                    sampling_period = 60,
                    units = "permil")

    # Write file
    nc_filename = output_filename("/data/shared/obs/",
                                  "NIWA",
                                  "GCMS",
                                  site.upper(),
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  "various")
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)

def niwa_run():
    
    niwa("ARH", "13ch4")
    niwa("BHD", "13ch4")


def uci_13ch4():


    for site in ["NWR", "MDO"]:
        
        if site == "NWR":
            df = pd.read_csv("/data/shared/obs_raw/UCI/nwrch4_edited.dat",
                             skiprows = 29,
                             sep = r"\s+",
                             names = ["date", "time_local", "ch4", "13ch4", "ch3d", "unknown1", "w1", "w2", "w3"],
                             parse_dates = {"time": ["date", "time_local"]},
                             index_col = "time",
                             na_values = "NM")
        elif site == "MDO":
            df = pd.read_csv("/data/shared/obs_raw/UCI/mdoch4_edited.dat",
                             skiprows = 24,
                             sep = r"\s+",
                             names = ["date", "time_local", "am_pm", "n", "ch4", "13ch4", "ch3d", "unknown1", "w1", "w2", "w3", "w4"],
                             parse_dates = {"time": ["date", "time_local", "am_pm"]},
                             index_col = "time",
                             na_values = "NM")

        df.index = df.index.tz_localize("US/Pacific")
        df.index = [np.datetime64(t) for t in df.index.tz_convert("UTC")]
        df.index.name = "time"

        df = df[["13ch4"]]
        df["13ch4_repeatability"] = 0.1
        df = df[np.isfinite(df["13ch4"])]

        # Sort and convert to dataset
        ds = xray.Dataset.from_dataframe(df.sort_index())
        
        # Add attributes
        ds = attributes(ds,
                        "13ch4",
                        site.upper(),
                        global_attributes = {"data_owner": "S. Tyler (UC. Irvine)"},
                        sampling_period = 60,
                        units = "permil")
    
        # Write file
        nc_filename = output_filename("/data/shared/obs/",
                                      "UCI",
                                      "GCMS",
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      "various")
        
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)



def uw_13ch4():
    
    sites = {"CGO": "CG",
             "OPW": "CP",
             "MLO": "ML",
             "BHD": "NZ",
             "BRW": "PB",
             "SMO": "SM"}
    
    for site in sites.keys():
        
        fname = glob.glob("/data/shared/obs_raw/UW/" + sites[site] + "Extract.prn.txt")[0]

        with open(fname, 'r') as f:
            lines = f.readlines()
                
        data = []
        
        data_started = False
        for li, line in enumerate(lines):
            if "Date" in line and "SIL ID #" in line:
                data_start = li
                data_started = True
            if data_started == True and len(line) < 5:
                data_end = li
                break
    
        data = lines[data_start+1:data_end]
    
        date = []
        ch4c13 = []
        
        for d in data:
            dsp = d.split()
            if len(dsp) > 1:
                date.append(dsp[1])
                ch4c13.append(float(dsp[-1]))
    
        df = pd.DataFrame(data = {"13ch4": np.array(ch4c13)}, index= pd.to_datetime(date))
        df.index.name = "time"
    
        df = df.groupby(level = 0, axis = 0).mean()
    
        df["13ch4_repeatability"] = 0.1
    
        # Sort and convert to dataset
        ds = xray.Dataset.from_dataframe(df.sort_index())
    
        # Add attributes
        ds = attributes(ds,
                        "13ch4",
                        site.upper(),
                        global_attributes = {"data_owner": "Quay, P.D. (UW)"},
                        sampling_period = 60,
                        units = "permil")
    
        # Write file
        nc_filename = output_filename("/data/shared/obs/",
                                      "UW",
                                      "GCMS",
                                      site.upper(),
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      "various")
        
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)
    
    
def uhei_13ch4():

    
    sites = {"ALT": "/data/shared/obs_raw/UHei/d13Cdata_Alert.txt",
             "IZA": "/data/shared/obs_raw/UHei/d13Cdata_Izana.txt"}

    for site, fname in sites.iteritems():
        with open(fname, "r") as f:
            lines = f.readlines()

        data_started = False
        for li, line in enumerate(lines):
            if "start" in line and "stop" in line and "CH4" in line:
                data_start = li
                data_started = True
            if data_started == True and line[0:4] == "----":
                data_end = li
                break

        columns = lines[data_start].split()
        columns[3] = columns[3] + "CH4"
        data_lines = lines[data_start+1:data_end]
        
        data = []
        date = []
        for line in data_lines:
            dsp = line.split()
            date.append(convert.decimal2time(float(dsp[0])))
            data.append([float(d) if d != "n.m." else np.NaN for d in dsp[2:]])
    
        df = pd.DataFrame(data = data, index = date, columns=columns[2:], dtype = float)

        df.rename(columns = {"d13CH4": "13ch4", "1sig": "13ch4_repeatability"},
                  inplace = True)
        df = df[["13ch4", "13ch4_repeatability"]]
        df.index.name = "time"
    
        # Sort and convert to dataset
        ds = xray.Dataset.from_dataframe(df.sort_index())
    
        # Add attributes
        ds = attributes(ds,
                        "13ch4",
                        site,
                        global_attributes = {"data_owner": "Ingeborg Levin (U. Heidelberg)"},
                        sampling_period = 60,
                        units = "permil")
    
        # Write file
        nc_filename = output_filename("/data/shared/obs/",
                                      "UHei",
                                      "GCMS",
                                      site,
                                      str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                      ds.species,
                                      "various")
        
        print("Writing " + nc_filename)
        
        ds.to_netcdf(nc_filename)


    # Neumayer 
    #########################################################################

    with open("/data/shared/obs_raw/UHei/d13Cdata_Neumayer.txt", "r") as f:
        lines = f.readlines()

    data_started = False
    for li, line in enumerate(lines):
        if "year" in line and "month" in line and "day" in line:
            data_start = li
            data_started = True
        if data_started == True and line[0:4] == "----":
            data_end = li
            break

    columns = lines[data_start].split()
    columns[4] = columns[4] + "CH4"
    data_lines = lines[data_start+1:data_end]

    data = []
    date = []
    for line in data_lines:
        dsp = line.split()
        date.append(dt(int(dsp[0]), int(dsp[1]), int(dsp[2])))
        data.append([float(d) if d != "n.m." else np.NaN for d in dsp[3:]])

    df = pd.DataFrame(data = data, index = date, columns=columns[3:], dtype = float)

    df.rename(columns = {"d13CH4": "13ch4", "1sig": "13ch4_repeatability"},
              inplace = True)
    df = df[["13ch4", "13ch4_repeatability"]]
    df.index.name = "time"

    # Sort and convert to dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())

    # Add attributes
    ds = attributes(ds,
                    "13ch4",
                    "NMY",
                    global_attributes = {"data_owner": "Ingeborg Levin (U. Heidelberg)"},
                    sampling_period = 60,
                    units = "permil")

    # Write file
    nc_filename = output_filename("/data/shared/obs/",
                                  "UHei",
                                  "GCMS",
                                  "NMY",
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  "various")
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)



def saws():
    ''' 
    South African Weather Service observations
    
    '''
    
    
    fname = "/data/shared/obs_raw/SAWS/CPT_CH4_1983_2015_All data.txt"
    df = pd.read_csv(fname, skiprows=10,
                     delimiter="\t", names = ["Time", "CH4"],
                     index_col = "Time", parse_dates=["Time"],
                     dayfirst=True)
    
    # remove duplicates
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"

    # remove NaN
    df = df[np.isfinite(df["CH4"])]
    
    # Sort and convert to dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())
    
    
    # Add attributes
    ds = attributes(ds,
                    "CH4",
                    "CPT",
                    global_attributes = {"inlet_height_magl": 30.,
                                         "data owner": "Casper Labuschagne (casper.labuschagne@weathersa.co.za)"},
                    sampling_period = 30)

    # Write file
    nc_filename = output_filename("/data/shared/obs/",
                                  "SAWS",
                                  "GC-FID-CRDS",
                                  "CPT",
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  "30m")
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)
    
    return(ds, df)


def uex():
    params = {
            "site" : "CVO",
            "directory" : "/data/shared/obs_raw/UEX/",
            "directory_output" : "/data/shared/obs/",
            "global_attributes" : {
                    "contact": "Elena Kozlova",
                    "averaging": "minute averaged OA-ICOS"
                    }
            }
    
    fnames = sorted(glob.glob(join(params["directory"],"*.txt")))

    df = []
    
    for fname in fnames:
        
        header_count = 0
        with open(fname, "r") as f:
            for line in f:
                if line[0] == "D":
                    header_count+=1
                    header = line.split()
                else:
                    break
    
        dff = pd.read_csv(fname, sep = r"\s+",
                         skiprows = header_count,
                         header = None,
                         names = header,
                         parse_dates = {"time": ["DATE", "TIME"]},
                         index_col = "time")
        
        dff = dff[["DATA", "ND", "SD"]]
    
        dff.rename(columns = {"DATA" : "CH4",
                                "ND": "CH4_number_of_observations",
                                "SD": "CH4_variability"},
                   inplace = True)

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
                    params['site'].upper(),
                    global_attributes = params["global_attributes"],
                    scale = "NOAA-2004")
    
    # Write file
    nc_filename = output_filename(params["directory_output"],
                                  "UEX",
                                  "OA-ICOS",
                                  params["site"],
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  site_params[params["site"]]["height"][0])
    
    print("Writing " + nc_filename)
    
    ds.to_netcdf(nc_filename)


def globalview_co2_obspack(site, height):

    height = str(height)
    site = str(site)
    
    if height == 'surface':
        fname = glob.glob("/data/shared/obs_raw/EUROCOM/ObsPack/co2_" + site.lower() + "*" + ".nc" )
    else:
        fname = glob.glob("/data/shared/obs_raw/EUROCOM/ObsPack/co2_" + site.lower() + "*" + "-" + height +"magl.nc" )
    
    instrument_dict = {'MHD':'CRDS',
                       'RGL':'CRDS',
                       'TAC':'CRDS',
                       'TTA':'CRDS',
                       'CBW':'NDIR',
                       'HUN':'NDIR'}

    if len(fname) == 0:
        print "Can't find file for site %s and height %s" %(site, height)
    elif len(fname) > 1:
        print "Ambiguous filename for site %s and height %s" %(site, height)
    elif len(fname) == 1:
        ds = xray.open_dataset(fname[0])
        
        species = ds.dataset_parameter
        if site in ['ces', 'CES']:
            site = 'CBW'
        else:
            site = ds.site_code
            
        if ds.value.units == "mol mol-1":
            print ds.value.values[0], float(unit_species[species.upper()])
            values = ds.value.values/float(unit_species[species.upper()])
            unc_values = ds.value_unc.values/float(unit_species[species.upper()])
            print values[0], unc_values[0]
        else:
            print "You need to create a unit conversion for the input units"
        
        ds2 = xray.Dataset({species.upper(): (['time'],values),
                            species.upper() + "_repeatability": (['time'],unc_values)},
#                            species.upper() + "_status_flag": (['time'],ds.obs_flag.values)},
                            coords = {'time': ds.time.values})

    
        global_attributes = {'origin': ds.obspack_name,
                             'original_filename': ds.dataset_name,
                             'citation': ds.obspack_citation,
                             'inlet_height_%s' %ds.dataset_intake_ht_unit : ds.dataset_intake_ht,
                             'main_provider_name' : ds.provider_1_name,
                             'main_provider_affiliation': ds.provider_1_affiliation,
                             'main_provider_email': ds.provider_1_email,
                             'fair_usage': ds.obspack_fair_use}
    
        ds2 = attributes(ds2,
                        species.upper(),
                        site.upper(),
                        global_attributes = global_attributes,
                        units = 'ppm',
                        scale = ds.dataset_calibration_scale,
                        sampling_period = None)
    
        # Write file
        directory_output = "/data/shared/obs/"
        instrument = instrument_dict[site.upper()]
        if site == 'MHD':
            inlet = "10m"
        else:
            inlet = "%im" %(int(float(ds.dataset_intake_ht)))
        
        nc_filename = output_filename(directory_output,
                                  "EUROCOM",
                                  instrument,
                                  site.upper(),
                                  str(ds2.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds2.species,
                                  inlet)
    
        print("Writing " + nc_filename)
        
        ds2.to_netcdf(nc_filename)

