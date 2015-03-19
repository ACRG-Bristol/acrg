# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 11:37:12 2014

Read CF-compliant netCDF file containing time, mole fractions, lon, lat, alt

This program just gets the variables and interprets the time co-ordinate

@author: chxmr
"""

from netCDF4 import Dataset
import datetime as dt

class read:
    def __init__(self, filename, varname):
    
        f=Dataset(filename, 'r')
        mf=f.variables[varname][:]
        secs=f.variables['time'][:]
        lon=f.variables['longitude'][:]
        lat=f.variables['latitude'][:]
        alt=f.variables['altitude'][:]
        
        time_ref=dt.datetime.strptime(f.variables['time'].units[14:], \
            "%Y-%m-%d %H:%M:%S")
        time=[time_ref + dt.timedelta(seconds=long(s)) for s in secs]
        
        self.time=time
        self.mf=mf
        self.lon=lon
        self.lat=lat
        self.alt=alt