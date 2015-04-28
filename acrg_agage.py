# -*- coding: utf-8 -*-
"""
Get AGAGE data, average, truncate and filter for baseline

Examples:

Get Mace Head CH4 observations and average into 3 hourly periods:

    import acrg_agage as agage
    time, ch4, sigma = agage.get("MHD", "CH4", average="3H")
    
Get Mace Head CH4 observations, truncate for 2010-2012 inclusive
and average into 6 hourly periods:

    import acrg_agage as agage
    time, ch4, sigma = agage.get("MHD", "CH4", startY=2010, endY=2013,average="6H")

Calculate Cape Grim monthly means, with baseline filtering:
    
    import acrg_agage as agage
    time, hfc134a, sigma = 
        agage.get("CGO", "HFC-134a", baseline=True, average="1M")

Created on Sat Dec 27 17:17:01 2014
@author: chxmr
"""

import numpy as np
import pandas as pd
import glob
from os.path import split, realpath
import re
from netCDF4 import Dataset
from acrg_time import convert
import json
import datetime as dt

root_directory="/data/shared/AGAGE_UKMO"

def synonyms(search_string, info):
    
    keys=info.keys()
    
    #First test whether site matches keys (case insensitive)
    out_strings = \
        [k for k in keys if k.upper() == search_string.upper()]

    #If not found, search synonyms
    if len(out_strings) is 0:
        for k in keys:
            matched_strings = \
                [s for s in info[k]["alt"] \
                    if s.upper() == search_string.upper()]
            if len(matched_strings) is not 0:
                out_strings = [k]
                break

    if len(out_strings) is 1:
        out_string = out_strings[0]
    else:
        out_string = None

    return out_string

def listsearch(varnames, species, species_info, label="alt"):
    
    for v in varnames:
        if species.upper() == v.upper():
            out_string = v
            break
        else:
            matched_strings=[s for s in species_info[species][label] \
                if s.upper() == v.upper()]
            if len(matched_strings) is not 0:
                out_string = v
                break
            else:
                out_string = None

    return out_string


def file_search_and_split(search_string):
    flag_files=glob.glob(search_string)
    fnames = sorted([split(f)[1] for f in flag_files])
    return fnames, [f.split("_") for f in fnames]


def quadratic_sum(x):
    return np.sqrt(np.sum(x**2))/float(len(x))


#Get Met Office baseline flags
def ukmo_flags(site, site_info):
    
    flag_directory=root_directory + "/flags/"
    fnames, file_info=file_search_and_split(
        flag_directory + "*.txt")
    file_site = [f[0] for f in file_info]
    file_site_string=listsearch(file_site, site, site_info)
    files=sorted([f for f in fnames if file_site_string in f])
    
    flag=[]
    flag_time=[]
    
    for f in files:
        flag_data=pd.io.parsers.read_csv(flag_directory + "/" + f, 
                                         delim_whitespace=True, skiprows=6)
        flag_time = flag_time + [dt.datetime(y, m, d, h, mi)
            for y, m, d, h, mi in 
            zip(flag_data["YY"], flag_data["MM"], flag_data["DD"], 
                flag_data["HH"], flag_data["MI"])]
        flag.append(flag_data["B"])
    
    flag=np.hstack(flag)
    
    return pd.DataFrame(flag, index=flag_time, columns=("flags",))


def get(site_in, species_in, startY=None,endY=None,
        height=None, baseline=False, average=None, 
        output_variability=False):
        
        #Get site info and species info from JSON files
        acrg_path=split(realpath(__file__))
        
        with open(acrg_path[0] + "/acrg_species_info.json") as f:
            species_info=json.load(f)

        with open(acrg_path[0] + "/acrg_site_info.json") as f:
            site_info=json.load(f)
        
        site = synonyms(site_in, site_info)
        if site is None:
            print("No site called " + site_in +
                ". Either try a different name, or add name to site_info.json.")
            return

        species = synonyms(species_in, species_info)
        if species is None:
            print("No species called " + species_in +
                ". Either try a different name, or add name to species_info.json.")
            return
        
        data_directory=root_directory + "/data/"

        if height is None:
            file_height_string = site_info[site]["height"][0]
        else:
            file_height_string = listsearch(height, site, site_info, 
                                            label="height")
            if file_height_string is None:
                print("Height " + height + " doesn't exist in site_info.json. "
                    + "Available heights are " + str(site_info[site]["height"])
                    + ". Leave blank for default height.")
                return
        
        #Get file info
        fnames, file_info=file_search_and_split(data_directory + "*.nc")
        if len(fnames) == 0:
            print("Can't find data files")
        file_inst = [f[0] for f in file_info]
        file_site = [f[1] for f in file_info]
        file_species = [re.split("-|\.", f[-1])[0] for f in file_info]
        file_height = [re.split("-|\.", f[-1])[1] for f in file_info]
        
        
        #Get file list
        file_species_string = listsearch(file_species, species, species_info)
        file_site_string = listsearch(file_site, site, site_info)
        files = [f for f, si, sp, hi in \
            zip(fnames, file_site, file_species, file_height) if \
            si == file_site_string and
            sp == file_species_string and
            hi == file_height_string]
#        files = sorted([f for f in fnames if 
#        #    inst in f and 
#            file_site_string in f and 
#            file_species_string in f and 
#            file_height_string in f])

        #Get files
        time=[]
        mf=[]
        dmf=[]
        vmf=[]
        status_flag=[]
        
        for f in files:
            ncf=Dataset(data_directory + f, 'r')
            time = time + (
                convert.sec2time(ncf.variables["time"][:], 
                                 ncf.variables["time"].units[14:]))
            ncvarname=listsearch(ncf.variables, species, species_info)
            #Get mole fraction
            mf.append(ncf.variables[ncvarname][:])
            #Get repeatability
            file_dmf=ncf.variables[ncvarname + "_repeatability"]
            if len(file_dmf) > 0:
                dmf.append(file_dmf[:])
            #Get variability
            file_vmf=ncf.variables[ncvarname + "_variability"]
            if len(file_vmf) > 0:
                vmf.append(file_vmf[:])
            #Get status flag
            file_flag=ncf.variables[ncvarname + "_status_flag"]
            if len(file_flag) > 0:
                status_flag.append(file_flag[:])
                
            ncf.close()
        
        time=np.array(time)
        mf=np.hstack(mf)
        if len(dmf) > 0:
            dmf=np.hstack(dmf)
        else:
            dmf=np.zeros((len(mf)))
        
        if len(vmf) > 0:
            vmf=np.hstack(vmf)
        else:
            vmf=np.zeros((len(mf)))

        status_flag = np.hstack(status_flag)
        
        #Flag out bad data
        wh=np.where(status_flag < 3)
        time=time[wh]
        mf=mf[wh]
        dmf=dmf[wh]
        vmf=vmf[wh]

        #Put into data frame        
        mfdf=(pd.DataFrame(zip(mf, dmf, vmf), 
                           index=time, columns=("mf", "dmf", "vmf"))).sort()
        if startY is not None:
            if endY is None:
                print 'Need to specify an end year as well'
                return
            else:
                if endY == startY:
                    print 'End year needs to be different from start year'
                    return
                mfdf=mfdf.truncate(dt.datetime(startY, 1, 1), dt.datetime(endY, 1, 1))
        #Do baseline filtering
        if baseline:
            #Get flags
            flagdf=ukmo_flags(site, site_info)
            #Truncate time series if no flags available
            mfdf=mfdf[min(flagdf.index):max(flagdf.index)]
            #Re-index flag time series to mf Data frame index
            mfdf["flag"]=flagdf.reindex(index=mfdf.index, method='pad')
            #Only retain background
            mfdf=mfdf[mfdf['flag']==10]
            #Remove flags
            mfdf.drop('flag', 1)
        
        if average is not None:
            mfdf=mfdf.resample(average, 
                               how={"mf": "median", 
                                    "dmf": quadratic_sum,
                                    "vmf": "mean"})

        #Output with either repeatability or variability
        if output_variability:
            return mfdf.index.to_pydatetime(), np.array(mfdf['mf']), \
                np.array(mfdf['vmf'])
        else:
            return mfdf.index.to_pydatetime(), np.array(mfdf['mf']), \
                np.array(mfdf['dmf'])
