# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 16:43:34 2014

Extract NOAA observations stored in CF convention format

Examples:

Get all site names from stored filenames:
    
    import acrg_noaa as noaa
    noaa.sites("CH4")
        ['ABP',
         'ALT',
         'AMS',
         ...]

Read monthly NOAA methane data for 2005 at Mace Head:

    import acrg_noaa as noaa
    data=noaa.read("MHD", "CH4", startYear=2005, endYear=2006, monthly=True)
    print data.time
        [datetime.datetime(2005, 1, 31, 0, 0),
         datetime.datetime(2005, 2, 28, 0, 0),
         datetime.datetime(2005, 3, 31, 0, 0),
         ...]

    print data.mf
        array([ 1847.11364746,  1850.75964355,  1843.13830566,  1983.51013184,
                1855.42993164,  1864.32507324,  1831.22009277,  1815.51501465,
                1836.64611816,  1873.51403809,  1869.26098633,  1856.56835938], dtype=float32)

sitelocs returns the FIRST lat/lon listed in each NOAA flask file.

      
@author: Matt Rigby
"""
from __future__ import print_function

import numpy as np
import acrg_cf as cf
import glob
import pandas
import datetime as dt

def directory(species):
    
    baseDirectory='/data/shared/GAUGE/'
    speciesDirectory={"CH4":"CH4/NOAA_flask/processed/", \
                        "CO2":"CO2/NOAA_flask/processed/"}
    return baseDirectory + speciesDirectory[species.upper()]

def filename(site, species):

    dirname=directory(species)
    filenames=glob.glob(dirname + "*")
    filename=[f for f in filenames if ('_' + site.upper() )in f.upper()]
    
    if len(filename) > 1:
        print("More than one matching filename " + site.upper())
        return 0
    if len(filename) == 0:
        print("No matching filename " + site.upper())
        return 0

    return filename[0]

def varname(species):
    
    varnames={"CH4":"ch4", \
                "CO2": "co2"}
    varname=varnames[species.upper()]
    return varname

def sites(species):
    
    filenames=glob.glob(directory(species) + '*')
    site=[]
    for f in filenames:
        s=f[f.rfind('/')+1:]
        s=s[s.find('_')+1:]
        site.append(s[0:s.find('_')].upper())
    
    site.sort()

    return site

class sitelocs:
    def __init__(self, species):
    
        sitenames = sites(species)    
    
        filenames=glob.glob(directory(species) + '*')
    
        lat=np.empty(np.shape(sitenames))
        lon=np.empty(np.shape(sitenames))
        alt=np.empty(np.shape(sitenames))
        
        for f in np.arange(len(filenames)):
            data_f = read(sitenames[f], species)
            lat[f] = data_f.lat[0]
            lon[f] = data_f.lon[0]
            alt[f] = data_f.alt[0]

        self.lat = lat
        self.lon = lon
        self.alt = alt


class read:
    def __init__(self, site, species, \
        startYear=2000, endYear=2014, monthly=True):
        
        #Read (CF-compliant) netCDF file
        fileOut = cf.read(filename(site, species), varname(species))
        
        time=fileOut.time
        mf=fileOut.mf
        
        #If there isn't an element at the start or end of the desired
        # time period, add a null mole fraction. This helps create a 
        # regular array if resampling
        if monthly:
            if min(time) > dt.datetime(startYear, 1, 31):
                time.append(dt.datetime(startYear, 1, 1))
                
                mf=np.append(mf, np.array(float('nan')))
            if max(time) < dt.datetime(endYear-1, 12, 31):
                time.append(dt.datetime(endYear-1, 12, 31))
                mf=np.append(mf, np.array(float('nan')))

        #Create Pandas time series
        ts=pandas.Series(mf, index=time)

        #Extract monthly means, if requested
        if monthly:
            ts=ts.resample('M', how='mean')

        ts=ts.truncate(dt.datetime(startYear, 1, 1), \
            dt.datetime(endYear, 1, 1))
        
        self.time=(ts.index.to_pydatetime()).tolist()
        self.mf=np.array(ts[:])
        self.lon=fileOut.lon
        self.lat=fileOut.lat
        self.alt=fileOut.alt
        