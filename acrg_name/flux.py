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

def write(lat, lon, time, flux, species, domain, year,
          comments = "", title="EDGAR emissions"):
    '''
    Write a flux file
    
    TODO: Add some error checking (e.g. check that domain is correct)
    '''
    if type(year) == dt.date:
        #time = year
        year = str(year.year)
    elif type(year) == str:
        #time = dt.datetime.strptime(year, '%Y')
        year = year
    elif isinstance(year[0], dt.date):
        time = year
        year = str(year[0].year)
    else:
        raise ValueError('Year needs to be in datetime.date or string format')
        
    if type(time) == dt.date:
        time = time
    elif type(time) == str:
        time = dt.datetime.strptime(time, '%Y')
    elif type(time) == np.ndarray:
        if type(time[0]) == dt.date:
            time = time
        else:
            raise ValueError('Time format not correct')
    else:
        time = time
       # raise ValueError('Time needs to be in datetime.date or string format')
    
        

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
    f.createDimension('time', len(time))
    f.createDimension('lat', len(lat))
    f.createDimension('lon', len(lon))
    
    #Header attributes
    f.title = title
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
    # Check that the flux is in the correct shape
    if np.shape(flux) == tuple((np.shape(lat)[0], np.shape(lon)[0], np.shape(time)[0])):
        ncflux=f.createVariable('flux', \
            'd', ('lat', 'lon', 'time'))
        ncflux[:,:,:]=flux
        ncflux.units='mol/m2/s'
    else:
        print 'The flux needs to be lat x lon x time which it is not'
        print 'Reshape your flux array and try again'
    f.close() 


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