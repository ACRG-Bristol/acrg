# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:13:48 2015

@author: chxmr
"""

import datetime as dt
import os
import netCDF4 as nc
import getpass
import acrg_time

def write(lat, lon, flux, species, domain, year,
          comments = ""):
    '''
    Write a flux file
    
    TODO: Add some error checking (e.g. check that domain is correct)
    '''
    
    if type(year) == dt.date:
        time = year
        year = str(year.year)
    elif type(year) == str:
        time = dt.datetime.strptime(year, '%Y')
        year = year
    
    species = species.lower()    
    
    #Open netCDF file
    ncname = '/data/shared/NAME/emissions/%s/%s_%s_%s.nc' %(domain, species, domain, year)

    if os.path.isfile(ncname) == True:
        answer = raw_input("You are about to overwrite an existing file, do you want to continue? Y/N")
        if answer == 'N':
            return
        elif answer == 'Y':
            f=nc.Dataset(ncname, 'w')
   
    elif os.path.isfile(ncname) == False:
        f=nc.Dataset(ncname, 'w')
        
    #Create dimensions
    f.createDimension('time', len([time]))
    f.createDimension('lat', len(lat))
    f.createDimension('lon', len(lon))
    
    #Header attributes
    f.title = "EDGAR emissions"
    f.author = getpass.getuser()
    f.date_created = str(dt.datetime.today())
    f.comments = comments

    #Time variable
    time_seconds, time_reference = acrg_time.convert.time2sec(time)
    nctime=f.createVariable('time', 'i', ('time',))
    nctime[:]= time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + str(time_reference)
    nctime.calendar='gregorian'

    #Longitude variable
    nclon=f.createVariable('lon', 'd', ('lon',))
    nclon[:]=lon
    nclon.long_name='longitude'
    nclon.standard_name='longitude'
    nclon.units='degrees_east'
    
    #Latitude variable
    nclat=f.createVariable('lat', 'd', ('lat',))
    nclat[:]=lat
    nclat.long_name='latitude'
    nclat.standard_name='latitude'
    nclat.units='degrees_north'
    
    #Mole fraction variable
    ncflux=f.createVariable('flux', \
        'd', ('lat', 'lon', 'time'))
    ncflux[:,:,:]=flux
    ncflux.units='mol/m2/s'

    f.close() 
