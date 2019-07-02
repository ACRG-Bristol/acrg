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
import pytz
import numpy as np
import datetime as dt
from collections import OrderedDict

# Site info file
acrg_path = getenv("ACRG_PATH")
data_path = getenv("DATA_PATH")
site_info_file = join(acrg_path, "acrg_site_info.json")
with open(site_info_file) as sf:
    site_params = json.load(sf)

# Set default obs folder
obs_directory = join(data_path, "obs_2018/")



def wdcgg_read(fname, species,
               permil = False,
               repeatability_column = None):
    '''
    Read data from World Data Centre for Greenhouse Gases
    '''
    
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
        output_columns[repeatability_column] = species.lower() + " repeatability"
    
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
              instrument = "",
              global_attributes = {},
              assume_repeatability = None):
    '''
    Read NIES data files
    '''
    
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
            df[sp + " repeatability"] = df[sp]*assume_repeatability
            global_attributes["Assumed_repeatability_%"] = int(assume_repeatability*100.)
        
        # Sort and convert to dataset
        ds = xr.Dataset.from_dataframe(df.sort_index())
        
        # Add attributes
        ds = attributes(ds,
                        sp,
                        site.upper(),
                        global_attributes = global_attributes,
                        sampling_period = 60,
                        units = "ppt")

        if assume_repeatability:
            ds[sp.lower() + " repeatability"].attrs["Comment"] = \
                "NOTE: This is an assumed value. Contact data owner."
    
        # Write file
        nc_filename = output_filename(obs_directory,
                                      network,
                                      instrument,
                                      site.upper(),
                                      ds.time.to_pandas().index.to_pydatetime()[0],
                                      ds.species)
        
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
    
    
    repeatability = {"CFC-11": 0.008}
    
    if fname.split(".")[1] == "xlsx":
        df = pd.read_excel(fname, parse_dates = [0], index_col = [0])
    elif fname.split(".")[1] == "csv":
        df = pd.read_csv(fname, sep = ",", parse_dates = [0], index_col = [0],
                         skipinitialspace=True)
    else:
        df = pd.read_csv(fname, sep = "\t", parse_dates = [0], index_col = [0])
    
    print("Assuming data is in JST. Check input file. CONVERTING TO UTC.")

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
    df[species + " repeatability"] = df[species]*repeatability[species]

    # Convert to dataset
    ds = xr.Dataset.from_dataframe(df)

    ds[species + " repeatability"].attrs["Comment"] = \
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
    nc_filename = output_filename(obs_directory,
                                  "NIES",
                                  "GCMS",
                                  site.upper(),
                                  ds.time.to_pandas().index.to_pydatetime()[0],
                                  ds.species)

    print("Writing " + nc_filename)
    ds.to_netcdf(nc_filename)


def nies_n2o_ch4(species):
    '''
    N2O files from NIES
    '''            
    params = {
        "site" : "COI",
        "scale": {
            "CH4": "NIES-94",
            "N2O": "NIES-96"},
        "instrument": {
                "CH4": "GCFID",
                "N2O": "GCECD"},
        "global_attributes" : {
                "contact": "Yasunori Tohjima (tohjima@nies.go.jp)" ,
                "averaging": "20 minutes"
                }
        }
    if species.lower() == 'ch4':
        fname = "/data/shared/obs_raw/NIES/COI/COICH4_Hourly_withSTD.TXT"
        df = pd.read_csv(fname, skiprows=1,
                 delimiter=",", names = ["Time", species.upper(), 'sd', 'N'],
                 index_col = "Time", parse_dates=["Time"],
                 dayfirst=True)
    elif species.lower() == 'n2o':
        fname = "/data/shared/obs_raw/NIES/COI/COIN2O_Hourly_withSTD.txt"
        df = pd.read_csv(fname, skiprows=1,
                 delimiter=",", names = ["Time", species.upper(), 'STD', 'n'],
                 index_col = "Time", parse_dates=["Time"],
                 dayfirst=True)
              
    
    print("Assuming data is in JST. Check input file. CONVERTING TO UTC.")
    
    df.index = df.index.tz_localize(pytz.timezone("Japan")).tz_convert(None) # Simpler solution

    # Sort
    df.sort_index(inplace = True)

    # Rename columns to species
    df.rename(columns = {df.columns[0]: species.upper()}, inplace = True)

    df.rename(columns = {df.columns[1]: species.upper() + " repeatability"}, inplace = True)
    df.rename(columns = {df.columns[2]: species.upper() + " number_of_observations"}, inplace = True)
    
    # Drop duplicates and rename index
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"
    
    # remove 9999
    df = df[df[species.upper()]<9999]

    # Convert to dataset
    ds = xr.Dataset.from_dataframe(df)
    

    ds = ds.where((ds[species.upper() + " repeatability"] < 9000), drop = True)
    
    # Add attributes

    ds = attributes(ds,
                    species.upper(),
                    params['site'].upper(),
                    global_attributes = params["global_attributes"],
                    scale = params["scale"][species.upper()])
   
    # Write file
    nc_filename = output_filename(obs_directory,
                                  "NIES",
                                  params["instrument"][species.upper()],
                                  params["site"],
                                  ds.time.to_pandas().index.to_pydatetime()[0],
                                  ds.species)
    
    print("Writing " + nc_filename)
    ds.to_netcdf(nc_filename)



def sio_carbon_cycle(species):
    ''' 
    Observations from the SIO carbon cycle group
    
    '''    
    params = {
        "site" : "THD",
        "scale": {
            "CO2": "WMO X2007",
            "DO2N2": "SIO",
            "APO": "SIO"},
        "instrument": {
                "CO2": "LICOR",
                "DO2N2": "DFCA",
                "APO": "DFCA"},
        "directory" : "/data/shared/obs_raw/SIOcarboncycle/",
        "global_attributes" : {
                "contact": "Timothy Lueker <tlueker@ucsd.edu>" ,
                "averaging": "8 minute"
                }
        }
            
    fname = "/data/shared/obs_raw/SIOcarboncycle/Trinidad_full.csv"
    
    if species.lower() == 'do2n2':
        species = 'do2n2'

    df = pd.read_csv(fname, skiprows=3,
             delimiter=",", names = ["YYYY","MM", "D", "HH", "mm", "VLV", "FG-C", "FG-O", 'CO2','DO2N2', 'APO'],
             parse_dates = [[0,1,2,3,4]], index_col = False,  engine='python')

    df["YYYY_MM_D_HH_mm"] = pd.to_datetime(df["YYYY_MM_D_HH_mm"], format = '%Y %m %d %H %M')     
    df = df.set_index("YYYY_MM_D_HH_mm", drop = True)

    # remove duplicates
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')              
    df.index.name = "time"

    # hack to filter spurious CO2 and APO data but need more permanent fix!!!! #
    df['CO2'] = np.array(df['CO2'], dtype=float)
    df['DO2N2'] = np.array(df['DO2N2'], dtype=float)
    df['APO'] = np.array(df['APO'], dtype=float)
    df = df[df['CO2']<450]
    df = df[df['CO2']>350]
    df = df[df['APO']<0]
    df = df[df['APO']>-400]

#    df.index = df.index.tz_localize(pytz.timezone("Africa/Johannesburg")).tz_convert(None) # Simpler solution
    if species.lower() == 'co2':
        df = df.drop(["VLV", "FG-C", "FG-O",'DO2N2', 'APO'], axis=1)
    elif species.lower() == 'do2n2':
        df = df.drop(["VLV", "FG-C", "FG-O", 'CO2', 'APO'], axis=1)   
    elif species.lower() == 'apo':
        df = df.drop(["VLV", "FG-C", "FG-O", 'CO2','DO2N2'], axis=1)  
       
    # remove NaN
    df[species.upper()] = np.array(df[species.upper()], dtype=float)
    df = df[np.isfinite(df[species.upper()])]    
    
    # Sort and convert to dataset
    ds = xr.Dataset.from_dataframe(df.sort_index())
    
    
    # Add attributes
    ds = attributes(ds,
                    species.upper(),
                    params['site'].upper(),
                    global_attributes = params["global_attributes"],
                    scale = params["scale"][species.upper()],
                    sampling_period = 8*60)


    # Write file
    nc_filename = output_filename(obs_directory,
                                  "SIOCC",
                                  params["instrument"][species.upper()],
                                  params["site"],
                                  ds.time.to_pandas().index.to_pydatetime()[0],
                                  ds.species)

    print("Writing " + nc_filename)
    ds.to_netcdf(nc_filename)




def ufrank(site = "TAU"):
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

def atto(species="CH4",inlet=1,version="v1"):
    '''
    Process Amazon Tall Tower Observatory data.

    Search string used is "*{version}*.dat" and currently reprocesses all matching files in
    folder.
    
    General notes:
       - For the period of 2013-06-01 to 2013-10-29 one of the two CRDS instruments (D9) was 
       substituted and then replaced. Multiple files are created if this time period is covered.
       - CO2 was measured by two instruments (one of which was substituted) as described above,
       so multiple files for CO2 covering the same time period are generated.
       - Due to the buffer system and the length of the inlets, there is some time delay 
       between the arrival of the air mass at the inlet and the time of the measurement. This
       has *not* been factored into the time points in the output.

    Example of input file (ATTO_Pic89_db30_v1_f201203_v1.dat) (line numbers added in []) 
    [1]   60 1001
    [...]
    [57]  20120301	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
    [58]	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999	9999
    [59]  [s]	[ppb]	[ppb]	[ppb]	[ppb]	[ppb]	[ppm]	[ppm]	[ppm]	[ppm]	[ppm]	[ppb]	[ppb]	[ppb]	[ppb]	[ppb]	[ppm]	[ppm]	[ppm]	[ppm]	[ppm]
    [60]  TimeUTC	D8CO1_1	D8CO1_3	D8CO1_4	D8CO1_5	D8CO1_6	D8CO2_1	D8CO2_3	D8CO2_4	D8CO2_5	D8CO2_6	D9CH4_1	D9CH4_3	D9CH4_4	D9CH4_5	D9CH4_6	D9CO2_1	D9CO2_3	D9CO2_4	D9CO2_5	D9CO2_6
    [61]  518400	117.39	113.91	118.48	120.00	123.70	399.57	397.99	402.17	430.73	431.23	1841.28	1841.43	1840.63	1839.96	1839.38	399.51	398.05	402.19	430.69	431.20
    [62]  520200	115.15	117.60	123.04	124.25	125.12	396.45	396.08	399.77	420.84	425.70	1839.65	1840.56	1841.10	1841.07	1841.22	396.42	396.14	399.87	420.83	425.62

    Args:
        species (str) :
            Options are:
                "CH4","CO2" and "CO".
        inlet (int) :
            5 different inlet heights are available corresponding to:
                1 - 79m
                3 - 53m
                4 - 38m
                5 - 24m
                6 - 4m
        version (str, optional) :
            Version of the data. At the moment only "v1" is available but other data may become
            available with new versions in future.
    
    TODO: There is a BUG that when splitting between the different instruments there is an 
    overlap between the end date of the first file (2013-06-01) and the start date of the secof
    file. Xarray selection here doesn't seem to be doing what we expect. May cause issues when
    using the data over this period.
    '''
    
    ## Time input notes:
    # TimeUTC column contains seconds since UTC 0
    #  Data have been averaged over 30 minute intervals.
    #  The time of measurement by the Picarro instruments has been taken.
    #  Due to the buffer system and the length of the inlets, there is some time delay 
    #  between the arrival of the air mass at the inlet and the time of the measurement.
    ## Instrument notes:
    # The data have been measured by CRDS instruments 'Picarro G1302 CKADS-18' (column names starting with 'D8') and 'Picarro G1301 CFADS-109' (column names starting with 'D9').
    # except for 2013-06-01 to 2013-10-29, where a 'Los Gatos OA-ICOS' analyser was used as replacement ('Los Gatos OA-ICOS' (column names starting with 'D9').
    # The data is calibrated on the scales NOAA-2004 (CH4), WMOX2007 (CO2), WMO CO X2004 (CO), which are implemented at MPI-BGC, via primary standards from NOAA, and are propagated to ATTO.
    
    params = {
        "site" : "ATO",
        "network":"ATTO",
        "scale": {
            "CO2": "WMO X2007",
            "CH4": "NOAA-2004",
            "CO": "WMO CO X2004"},
        "instrument": {
                "CO2": ["CRDS-CKADS-18","CRDS-CFADS-109"],
                "CH4": ["CRDS-CFADS-109"],
                "CO": ["CRDS-CKADS-18"]},
        "substitute" : {
                "start_date":"2013-06-01",
                "end_date":"2013-11-01",
                "org_instrument":'CRDS-CFADS-109',
                "sub_instrument":"OA-ICOS"},
        "inlet" :
            {1:"79m",
             3:"53m",
             4:"38m",
             5:"24m",
             6:"4m"},
        "units":
            {"CO2":"ppm",
             "CH4":"ppb",
             "CO":"ppb"},
        "measurement_period":30,
        "directory" : "/data/shared/obs_raw/ATTO/",
        "directory_output" : "/data/shared/obs_2018",
        "global_attributes" : OrderedDict([
                ("data_owner","Jost Lavric"),
                ("data_owner_email", "jost.lavric@bgc-jena.mpg.de"),
                ("data_creator","David Walter"),
                ("data_creator_email","david.walter@mpic.de")])
               }

    search_str = join(params["directory"],"*{}*.dat".format(version))
    fnames = glob.glob(search_str)
    fnames.sort()

    title_line = 0
    na_values = "9999"
    time_col = "TimeUTC"

    def define_meas_time(t,date,date_fmt="%Y%m%d"):
        if isinstance(date,str):
            date = dt.datetime.strptime(date,date_fmt)
            date = np.datetime64(date)
        datetime = date + np.timedelta64(t,'s')# + np.timedelta64(int(meas_period/2.),'s')
        return datetime
    
    for instrument in params["instrument"][species]:
        
        ds_list = []
        for fname in fnames:
    
            if instrument == "CRDS-CKADS-18":
                col_start = "D8"
            elif instrument == "CRDS-CFADS-109":
                col_start = "D9"
            
            if species.upper() == "CO":
                species_str = species.upper()+"1"
            else:
                species_str = species.upper()
                
            data_col = "{col}{species}_{inlet}".format(col=col_start,species=species_str,inlet=inlet)#,"D9CH4_3","D9CH4_4","D9CH4_5","D9CH4_6"]
            
            ## Number of header lines is included in first line of file and date is included
            # two before the end of the header lines. Extract these values from each file.
            with open(fname) as fl:
                header_lines = int(fl.readline().split(' ')[0])-1
                fl.seek(0) # Reset to top of file
                line_containing_date = header_lines-2
                date_str = fl.readlines()[line_containing_date-1].split('\t')[0]

            df = pd.read_csv(fname,header=title_line,delim_whitespace=True,
                             skiprows=header_lines,usecols=[time_col,data_col],na_values=na_values)
            
       
            ## Remove any NaN values
            df = df[np.isfinite(df[data_col])]    
        
            ## Rename data column to species name
            df = df.rename({data_col:species.upper()},axis="columns")
            
            ## Construct UTC time
            df["time"] = df[time_col].apply(define_meas_time,date=date_str)
            df = df.drop(time_col,axis=1)
            df = df.set_index("time")
    
            # TODO: No error included so could arbitrarily set to e.g. 5% of data?
            #df[species+"_repeatability"] = df[species]*0.05
            
            if not df.empty:
                ds_list.append(xr.Dataset.from_dataframe(df.sort_index()))
            else:
                print("No data for {} for {} instrument from file: {}".format(species,params["site"],split(fname)[1]))
    
        ds = xr.concat(ds_list,dim="time")
    
        ## One instrument was substituted for a period of time so need to create different
        # file for these dates.
        ds_mult = []
        if params["substitute"]["org_instrument"] == instrument:
            instrument_sub = params["substitute"]["sub_instrument"]
            
            sub_start = np.datetime64(params["substitute"]["start_date"])
            sub_end = np.datetime64(params["substitute"]["end_date"])
            
            start = ds.time.values[0]
            end = ds.time.values[-1]
    
            date_range = [start,sub_start,sub_end,end]        
            instruments = [instrument,instrument_sub,instrument]
    
            for i,ins in enumerate(instruments):
    
                if i != len(instruments)-1:
                    s = date_range[i]
                    e = date_range[i+1]-np.timedelta64(1,"D")
                else:
                    s = date_range[i]
                    e = date_range[i+1]
                
                ds_1 = ds.sel(**{"time":slice(s,e)})
                #print("Sub-set dataset start and end time",ds_1["time"].values[0],ds_1["time"].values[-1])
                num_t = ds_1.dims["time"]
    
                if num_t > 0:
                    ds_mult.append((ins,ds_1))
        else:
            ds_mult.append((instrument,ds))
                
        for inst,ds in ds_mult:        
            
            ## Add attributes including global attributes
            global_attributes = params["global_attributes"]
            global_attributes["averaging"] = "{} minutes".format(params["measurement_period"])
            global_attributes["inlet_magl"] = params["inlet"][inlet]
            
            ds = attributes(ds,
                            species,
                            params["site"],
                            global_attributes = global_attributes,
                            scale = params["scale"][species],
                            sampling_period = params["measurement_period"]*60,
                            units = params["units"][species])
         
            ## Define filename and write to netcdf
            nc_filename = output_filename(params["directory_output"],
                                              params["network"],
                                              inst,
                                              params["site"].upper(),
                                              ds.time.to_pandas().index.to_pydatetime()[0],
                                              ds.species,
                                              inlet=params["inlet"][inlet],
                                              version="v1")
        
            print(" ... writing " + nc_filename)
                
            ds.to_netcdf(nc_filename)
    