# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:13:48 2015

@author: chxmr
"""

import numpy as np
import datetime as dt
import os
import netCDF4 as nc
import getpass
import acrg_time
import re
import numpy as np 
import xarray as xray
import collections as c
import sys
import pandas as pd

def write(lat, lon, time, flux, species, domain,
          source, title, prior_info_dict,
          regridder_used = 'acrg_grid.regrid.regrid_3D',
          copy_from_year = None, climatology = False, flux_comments = None,
          output_directory = '/data/shared/NAME/emissions/'):
    '''
    Write a flux file
    
    -time- should be in numpy.datetime64 format in either an array (format from xarray) or a datetime index (format from pandas).
    See 'Creating datetime64 data' : http://xarray.pydata.org/en/stable/time-series.html
    
    -flux- should be in mol/m2/s and -units- should be mol/m2/s.
    
    -source- is the sources in file. E.g. 'ff' (fossil fuel), 'agriculture'.
    source = None if the file contains all sources of a species. 
    If multiple sources: -source- is a chain of sources: 'waste-and-agriculture' (with hyphens between words).
    
    -title- gives more information about what is in the file, e.g. Fossil Fuel CO2.
    
    prior_info_dict = {'NAME_OF_PRIOR' : ['VERSION','RAW RESOLUTION', 'REFERENCE']}
    
    -copy_from_year- if the data is the same as another year but with a different timestamp give the original year here as a string
    
    -climatology- if the data is a climatology set this to True and give detail using flux_comments 
    
    TODO: Add some error checking (e.g. check that domain is correct)
    '''
    
    print "WARNING: Make sure time stamp is start of time period (i.e. 1st of month\
            for monthly data or 1st January for yearly data)."
    print "WARNING: Make sure coordinates are centre of the gridbox."
    print "WARNING: Make sure fluxes are in mol/m2/s."
        
    if source == None:
        file_source = species
        source_name = species + '-total'
    else:
        file_source = species + '-' + source
        source_name = file_source
    
    file_source = file_source.lower()
    species = species.lower()  
        
    # Check that the flux is in the correct shape
    if np.shape(flux) != tuple((np.shape(lat)[0], np.shape(lon)[0], np.shape(time)[0])):
        print "Flux doesn't have dimensions lat x lon x time"
        print "Reshape your flux array and try again"
        return
        
    if type(time[0]) == np.datetime64:
        time=time
    else:
        sys.exit('Time format not correct, needs to be type numpy.datetime64')

        
    #Open netCDF file
    year = pd.DatetimeIndex([time[0]]).year[0]
    if copy_from_year != None:
        ncname = output_directory + '%s/%s_%s_%s_copy-from-%s.nc' %(domain, file_source, domain, year, copy_from_year)
    if climatology == True:
        ncname = output_directory + '%s/%s_%s_%s_climatology.nc' %(domain, file_source, domain, year)
    else:
        ncname = output_directory + '%s/%s_%s_%s.nc' %(domain, file_source, domain, year)

    if os.path.isfile(ncname) == True:
        answer = raw_input("You are about to overwrite an existing file, do you want to continue? Y/N")
        if answer == 'N':
            sys.exit()
        elif answer == 'Y':
            pass
    
    flux_attrs = {"source" : source_name,
                  "units" : 'mol/m2/s',
                  "species" : species} 
    
    lat_attrs = {"long_name" : "latitude",
                 "units" : "degrees_north",
                 "notes" : "centre of cell"}
    
    lon_attrs = {"long_name" : "longitude",
                 "units" : "degrees_east",
                 "notes" : "centre of cell"}
    

    glob_attrs = {"title":title,
                  "author" : getpass.getuser(),
                  "date_created" : np.str(dt.datetime.today()),
                  "prior_used" : ["%s: Version %s, raw resolution %s, reference %s. " %(i, prior_info_dict[i][0], prior_info_dict[i][1], prior_info_dict[i][2]) for i in prior_info_dict.keys()],
                  "regridder" : "Created using %s" %regridder_used}

    if flux_comments != None:
        glob_attrs['comments'] = flux_comments
        if copy_from_year != None:
            glob_attrs['comments'] = "Fluxes copied from year %s." %copy_from_year + glob_attrs['comments']
    
    if copy_from_year != None:
        glob_attrs['comments'] = "Fluxes copied from year %s." %copy_from_year

    flux_ds = xray.Dataset({'flux':(['lat','lon','time'], flux, flux_attrs)},
                              coords = {'lat' : lat,
                                        'lon' : lon,
                                        'time' : time},
                              attrs = c.OrderedDict(glob_attrs))
    
    flux_ds.lat.attrs = lat_attrs
    flux_ds.lon.attrs = lon_attrs
    flux_ds.time.attrs['notes'] = "Start of time period"
    

    flux_ds.flux.encoding = {'zlib':True}                        
    flux_ds.to_netcdf(ncname, mode='w')    


class EDGARread:
    def __init__(self, filename_of_EDGAR_emissions):

        f = nc.Dataset(filename_of_EDGAR_emissions, 'r')
    
        #Get grid
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get flux of species
#        variables = f.variables.keys()
#        species = str(variables[2])
        variables = [str(i) for i in f.variables.keys()]
        for i in variables:
            while i not in ['lat','lon','time']:
                species = i
                if species is not None:
                    break
        flux = f.variables[species][:,:]
        units = f.variables[species].units

        f.close()
        
        #Get year and datetime date of Edgar file
        filename = os.path.split(filename_of_EDGAR_emissions)[1]
        match = re.compile(r'_\d{4}_')
        m = match.findall(filename)
        if len(m) == 1:
            year = m[0].strip('_')
            date = dt.datetime.strptime(year, '%Y')
        elif len(m) > 1:
            print "Can't find correct date."
            year = None
            date = None
        
        species = species.split('_')[1]

        self.lon = lon
        self.lat = lat
        self.flux = flux
        self.species = species
        self.units = units
        self.year = year
        self.date = date