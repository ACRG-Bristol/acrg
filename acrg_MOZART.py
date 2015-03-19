# -*- coding: utf-8 -*-
"""


File to read in standard MOZART output file
Output will contain the timestamp, lat, lon, surface emissions and concentrations
If you want any other fields then you'll need to add them in
    import acrg_MOZART as mozart    
    data=mozart.read('filename')
    print data.time
    
Can also extract data for sites given lats and lons
sitefile = a text file with the site acronym, lat and then lon separated by spaces or tabs
    import acrg_MOZART as mozart
    data=mozart.filter('filename', 'sitefile')

    

"""
import netCDF4
import numpy as np
import datetime as dt
from acrg_grid.hybridcoords import hybridcoords as hybrid_coords
import pdb

 # Class to read in the data
class read:
    def __init__(self, filename):
    
        if type(filename) == tuple:
            filename = filename[0]
        
        data=netCDF4.Dataset(filename)
        
        conc_varname = (str(data.__getattribute__('title'))).strip()+'_VMR_avrg'
        emiss_varname = (str(data.__getattribute__('title'))).strip()+'_SRF_EMIS_avrg'
        
        
        conc = data.variables[conc_varname][:]
        emis = data.variables[emiss_varname][:]
        date = data.variables['date'][:] # Date YYYYMMDD
        secs = data.variables['datesec'][:] # Seconds to be added to above date to make date-time variable
        lon = data.variables['lon'][:]
        lat = data.variables['lat'][:]
        PS=data.variables['PS'][:]
        P0=data.variables['P0'][:]
        hyai=data.variables['hyai'][:]        
        hybi=data.variables['hybi'][:]
        
        # Split up the date based on position        
        year = np.asarray([int(((date.astype('str'))[i])[0:4]) for i in np.arange(len(date))])
        month = np.asarray([int(((date.astype('str'))[i])[4:6]) for i in np.arange(len(date))])
        day = np.asarray([int(((date.astype('str'))[i])[6:8]) for i in np.arange(len(date))])

        
        dt_time = [dt.timedelta(seconds=(secs[i]).astype('int')) for i in np.arange(len(date))]
        dt_date = [dt.datetime(year[i],month[i],day[i]) for i in np.arange(len(date))]

        time_t = [dt_time[i] + dt_date[i] for i in np.arange(len(date))]

        P = np.empty((len(date), len(hyai)-1, len(lat), len(lon)))

        for i in np.arange(len(date)):         
            P_i = hybrid_coords(hyai, PS[i,:,:],  B=hybi, P0=P0, half=1)
            P[i,:,:,:] = P_i
         
        self.time = time_t
        self.emis = emis
        self.conc = conc
        self.lon = lon
        self.lat = lat
        self.P0 =P0
        self.PS =PS
        self.hyai=hyai
        self.hybi=hybi
        self.pressure = P
        self.date = date
        self.filename = filename         
        self.species = str(data.__getattribute__('title')).strip()
        self.case = str(data.__getattribute__('case')).strip()
        self.concunits = data.variables[conc_varname].getncattr('units')
        self.emissunits = data.variables[emiss_varname].getncattr('units')
        self.pressureunits = data.variables['P0'].getncattr('units')
                  
 # Class to read in the site file GONZI netcdf version       
class read_sitefile_nc:
    def __init__(self, sitefile):
        
        if type(sitefile) == tuple:
            data=netCDF4.Dataset(sitefile[0], 'r')
        elif  type(sitefile) == str:
            data=netCDF4.Dataset(sitefile, 'r')
        
        # Set flag for ferry which has no ALT
        if type(sitefile) == tuple:        
            if 'acf' in sitefile[0]:
                ferry_flag =0
            if 'ferry' in sitefile[0]:
                ferry_flag = 1
        
        if type(sitefile) == str:        
            if 'acf' in sitefile:
                ferry_flag = 0
            if 'ferry' in sitefile:
                ferry_flag = 1
      
        
        
        # Extract the data    
        lat = data.variables['LAT'][:]
        lon = data.variables['LONG'][:]
        if ferry_flag == 0:
            alt = data.variables['ALTI'][:] # Date YYYYMMDD
        if ferry_flag == 1:
             alt = np.empty(len(lat))
             
        Year = data.variables['YYYY'][:] # Seconds to be added to above date to make date-time variable
        Month = data.variables['MM'][:]
        Day = data.variables['DD'][:]
        SecSinceMidnight=data.variables['Time_Slots'][:]
        
        
        # Create the time variable
        time_t = [dt.datetime(Year[i],Month[i],Day[i]) + dt.timedelta(seconds=(SecSinceMidnight[i]).astype('int')) for i in np.arange(len(Day))]

        # Remove the superfluous points
        good_i = (np.where(lat[:,0] != -9999))[0]
        
        self.time = [time_t[good_i[i]] for i in np.arange(len(good_i))]
        self.lat = lat[good_i]
        self.lon = lon[good_i]
        self.alt = alt[good_i]
        self.filename = sitefile

 # Class to estimate the pressure levels from a given altitude
 # This uses the scale height and a simple atmospheric model
 # P = P0 * e^(-Z/H)
 # Where P = Pressure
 # P0 = surface pressure in Pa
 # -Z = altitude (km)
 # H = scale height (km) this is hard coded to be 7.64km which is the global mean
 # This assumes that the altitude is in m unless other units are given
units = 'm'
class calc_pressure:
    def __init__(self, altitude, P0, units=units):

        from math import exp
        
        if units == 'm':
            altitude = altitude/1000.0
        elif units =='km':
            altitude = altitude
        elif units =='miles':
            altitude = altitude*1.609344
        elif (units in {'m','km','miles'}) == False:
            print 'the units you have given are not a listed option'
            print 'please give the altitude in m, km or miles'
        
        pressure = P0*exp((-1*altitude)/7.64) 
        
        self.pressure = pressure


 # Class to read in the site file txt version
class read_sitefile_txt:
    def __init__(self, sitefile):
        
        if type(sitefile) == tuple:
            sitedata=np.genfromtxt(sitefile[0], dtype=str, skip_header=1)
        elif type(sitefile) == str:
            sitedata=np.genfromtxt(sitefile, dtype=str, skip_header=1)
        
        sitenames = sitedata[:,0]
        lat = sitedata[:,1]
        lon = sitedata[:,2]
        alt = sitedata[:,3]
        
        self.sitenames = sitenames
        self.lat = lat.astype('float')
        self.lon = lon.astype('float')
        self.alt = alt.astype('float')
        self.filename = sitefile
        
 # Class to extract the lat lon for the correct site from the text site files  
class extract_site_info:
    def __init__(self, sitefile, sitename):
        
        if sum([sitefile[0].find('acf'), sitefile[0].find('ferry')])  < -1:
            sitedata=read_sitefile_txt(sitefile)
        else:
            sitedata=read_sitefile_nc(sitefile)
            
        site_i = (np.where(sitedata.sitenames==sitename))[0]      
        
        self.sitename = sitedata.sitenames[site_i]        
        self.lat = sitedata.lat[site_i]
        self.lon = sitedata.lon[site_i]
        self.alt = sitedata.alt[site_i]
         
        
# Class to filter the data
# This uses an individual MOZART history file
# and a site file which contains the lat, lon and alt for fixed sites
Ave = 0
class data_filter_fixed:
    def __init__(self, mzfile, sitefile, Ave=Ave):
        
       
        import acrg_grid.haversine as haver
        
    
        # Read MOZART file name
        data = read(mzfile)
        
        
        # Read site info file
        siteinfo = read_sitefile_txt(sitefile)
            
        
        # create parameters to store output
        pressures = np.empty((len(siteinfo.sitenames), len(data.time)))      
        concs = np.empty((len(siteinfo.sitenames), len(data.time)))
        concs_sd = np.empty((len(siteinfo.sitenames), len(data.time)))
        
        emissions = np.empty((len(siteinfo.sitenames), len(data.time)))
        
        
        # loop through each lat/lon in the siteinfo
        for j in np.arange(len(siteinfo.lat)):
            
            # Use the haverseine code to determine the closest surface grid point to the given location
            haversine =  haver.multipledistances([siteinfo.lat[j], siteinfo.lon[j]] , data.lat, data.lon)     
            
            # As it's tall tower then we'll want the data for all time stamps
            # Extract column pressure at the correct lat/lon for all timestamps
            # This should be a time by lev array
            column_P = np.squeeze(data.pressure[:,[np.arange(len(data.hyai)-1)],haversine.mindist_index[0], haversine.mindist_index[1]])        
            
            if Ave != 0:
                # Just take the average of the bottom 7 levels that's all the data roughly < 1 km
    
                # Put data for each lat/lon (i.e. j value) into output arrays 
                # I want the lowest 7 levels which are actually the LAST 7 not the first
                index = np.empty(7)
                index[:] = len(data.hyai)-1
                index = index - np.arange(7) -1
                pressures[j,:] = np.mean(np.squeeze(column_P[:,[index]]), axis=1)
                
                # Extract the corresponding concentrations and emissions
                conc_j = np.squeeze(data.conc[:,:, haversine.mindist_index[0], haversine.mindist_index[1]])
                concs[j,:] = np.mean(np.squeeze(conc_j[:,[index]]), axis=1)
                concs_sd[j,:] = np.std(np.squeeze(conc_j[:,[index]]), axis=1)
                
    
            else:
                
                # Convert the altitude to a pressure using P0 from the model
                site_pressure = calc_pressure(siteinfo.alt[j], data.P0)             
                print 'Assuming that the altitudes are given in m'
                  
                # Loop through each time step and match site pressure to column level
                lev_i = np.empty(len(data.time))
                
                for i in np.arange(len(data.time)):
                    lev_i[i] = np.where(abs(column_P[i,:] - site_pressure.pressure) == min(abs(column_P[i,:] - site_pressure.pressure)))[0]

                # Extract the corresponding concentrations and emissions
                # Put data for each lat/lon (i.e. j value) into output arrays            
                pressures[j,:] = np.squeeze(column_P[:,lev_i.astype('int')])
                #pdb.set_trace()
                concs[j,:] = np.squeeze(data.conc[:,lev_i.astype('int'), haversine.mindist_index[0], haversine.mindist_index[1]])
            
            
            
            emissions[j,:] = np.squeeze(data.emis[[np.arange(len(data.time))], haversine.mindist_index[0], haversine.mindist_index[1]])
                        
                    
                    
        self.sitenames = siteinfo.sitenames
        self.sitetype = 'TT'
        self.lat = siteinfo.lat
        self.lon = siteinfo.lon
        self.alt = siteinfo.alt
        self.pressure = np.squeeze(pressures)
        self.time = data.time
        self.conc = np.squeeze(concs)
        self.conc_sd = np.squeeze(concs_sd)
        self.emis = np.squeeze(emissions)
        self.case = data.case.strip()
        self.species = data.species.strip()
        self.mzfile = data.filename.strip()
        self.sitefile = sitefile        
        self.concunits = data.concunits
        self.emissunits = data.emissunits
        self.pressureunits = data.pressureunits         
        
# Class to filter the data
# This uses an individual MOZART history file
# and a site file which contains the lat, lon, alt and time for moving sites
class data_filter_moving:
    def __init__(self, mzfile, sitefile):   
        
       
        import acrg_grid.haversine as haver
        
    
        # Read MOZART file name
        data = read(mzfile)
        
        
        # Read site info file
        siteinfo = read_sitefile_nc(sitefile)
        
        # Define site_names as are they're not a given parameter in the site info file
        if type(siteinfo.filename) == tuple:        
            if 'acf' in siteinfo.filename[0]:
                sitename = 'acf'
            if 'ferry' in siteinfo.filename[0]:
                sitename = 'ferry'
        
        if type(siteinfo.filename) == str:        
            if 'acf' in siteinfo.filename:
                sitename = 'acf'
            if 'ferry' in siteinfo.filename:
                sitename = 'ferry'
        
        # Only want to match times that are within the range of the file
        modeltime_range = [data.time[0], data.time[1] - data.time[0] + data.time[-1]]
        
        if sum(np.asarray(siteinfo.time) >= modeltime_range[0])/len(siteinfo.time) + sum(np.asarray(siteinfo.time) <= modeltime_range[1]) >= 2:
            starttime_data = np.where((np.asarray(siteinfo.time) >= modeltime_range[0]))[0][0]      
            endtime_data = np.where((np.asarray(siteinfo.time) <= modeltime_range[1]))[0][-1]      
        
        
            # extract times, lats, lons and alts that correspond to the time range of the model output file
            times = siteinfo.time[starttime_data:endtime_data+1]
            lats = siteinfo.lat[starttime_data:endtime_data+1, :]
            lons = siteinfo.lon[starttime_data:endtime_data+1, :]
            alts = siteinfo.alt[starttime_data:endtime_data+1, :]
            
            
            # create parameters to store output
            pressures = np.empty(len(times))      
            concs = np.empty(len(times))
            emissions = np.empty(len(times))
            
            
            # loop through each lat/lon in the same time range as the model output
            for j in np.arange(len(lats)):
                
                # Use the haverseine code to determine the closest surface grid point to the given location
                haversine =  haver.multipledistances([lats[j], lons[j]] , data.lat, data.lon)     
                
                
                # Convert the altitude to a pressure using P0 from the model
                site_pressure = calc_pressure(alts[j], data.P0)             
                print 'Assuming that the altitudes are given in m'
                
                
                # Match obs timestamp to the closest model timestamp
                timeindex_j = np.where(abs(np.asarray(data.time) - times[j]) == min(abs(np.asarray(data.time) - times[j])))[0]
    
                    
                # Extract column pressure at the correct lat/lon for all timestamps
                column_P = np.squeeze(data.pressure[timeindex_j,:,haversine.mindist_index[0], haversine.mindist_index[1]])
                    
                    
                # Match site pressure to column level
                lev_i = np.where(abs(column_P - site_pressure.pressure) == min(abs(column_P - site_pressure.pressure)))[0]
                         
                         
                # Extract the corresponding concentrations and emission
                # Put data for each lat/lon (i.e. j value) into output arrays            
                pressures[j] = np.squeeze(column_P[lev_i])
                concs[j] = np.squeeze(data.conc[timeindex_j,lev_i, haversine.mindist_index[0], haversine.mindist_index[1]])
                emissions[j] = np.squeeze(data.emis[timeindex_j,haversine.mindist_index[0], haversine.mindist_index[1]])
    
                
                
            self.sitenames = sitename
            self.sitetype = sitename
            self.lat = lats
            self.lon = lons
            self.alt = alts
            self.pressure = np.squeeze(pressures)
            self.time = times
            self.conc = np.squeeze(concs)
            self.emis = np.squeeze(emissions)
            self.case = data.case.strip()
            self.species = data.species.strip()
            self.mzfile = mzfile
            self.sitefile = sitefile
            self.concunits = data.concunits
            self.emissunits = data.emissunits
            self.pressureunits = data.pressureunits
            
        else:
            print 'There were no time matched ' + sitename + ' points for file: ' + data.filename
            self.time = 0
            
# Class to write out the data from a moving platform
filename = 0
outdir = 0
mzfile = 0
class write_ncdf_moving:
    def __init__(self, filtereddata, outdir = outdir, filename=filename, mzfile=mzfile):
        
        import os
        
        if mzfile == 0:
            mzfile = filtereddata.mzfile
        
        if type(outdir) == int:
            outdir = os.path.dirname(mzfile)
    
        if type(filename) == int:
            if type(mzfile) == tuple:  
                filetimestamp = mzfile[0][mzfile[0].find('h0')+ 2: mzfile[0].find('.nc')]
            if type(mzfile) == str:
                filetimestamp = mzfile[mzfile.find('h0')+ 2: mzfile.find('.nc')]
            filename = filtereddata.species + '_' + filtereddata.case + '.MZT.h0.'+  filetimestamp + '_' + filtereddata.sitetype+'.nc'
        
	print 'writing file: ' + outdir + '/'+ filename
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(outdir + '/'+ filename, 'w')
        
        # Create the dimensions
        ncF.createDimension('time', len(filtereddata.time))
        ncF.createDimension('sitenames', len(filtereddata.sitenames))
        
        # Make some global attributes
        ncF.mzfile = filtereddata.mzfile
        ncF.sitefile = filtereddata.sitefile
        ncF.case = filtereddata.case
        ncF.species = filtereddata.species
        ncF.sitetype = filtereddata.sitetype
        
        # Create the variables
        #ncsitenames = ncF.createVariable('sitename','s', ('sitenames',))
        ncsitenames = ncF.createVariable('sitename',type(filtereddata.sitenames[0]), ('sitenames',))        
        
        ncyear = ncF.createVariable('year', 'i', ('time',))
        ncmonth = ncF.createVariable('month', 'i', ('time',))
        ncday = ncF.createVariable('day', 'i', ('time',))
        nchour = ncF.createVariable('hour', 'i', ('time',))
        ncminute = ncF.createVariable('minute', 'i', ('time',))    
        
        nclon = ncF.createVariable('lon', 'f', ('time',))
        nclat = ncF.createVariable('lat', 'f', ('time',))
        ncalt = ncF.createVariable('alt', 'f', ('time',))
        ncpress = ncF.createVariable('pressure', 'f', ('time',))
        ncconc = ncF.createVariable('conc', 'f', ('time'))
        ncemiss = ncF.createVariable('emiss', 'f', ('time',))
        
        # Fill the variables
        ncsitenames[:] = filtereddata.sitenames        
        
        ncyear[:] = [filtereddata.time[t].year for t in np.arange(len(filtereddata.time))]
        ncmonth[:] = [filtereddata.time[t].month for t in np.arange(len(filtereddata.time))]        
        ncday[:] = [filtereddata.time[t].day for t in np.arange(len(filtereddata.time))] 
        nchour[:] = [filtereddata.time[t].hour for t in np.arange(len(filtereddata.time))]
        ncminute[:] = [filtereddata.time[t].minute for t in np.arange(len(filtereddata.time))]        
   
        nclon[:] = filtereddata.lon
        nclat[:] = filtereddata.lat
        ncalt[:] = filtereddata.alt
        ncpress[:] = filtereddata.pressure
        ncconc[:] = filtereddata.conc
        ncemiss[:] = filtereddata.emis
        
        
        # Give the variables some attributes        
        nclon.units = 'Degrees east'
        nclat.units = 'Degrees north'
        ncconc.units = filtereddata.concunits
        ncemiss.units = filtereddata.emissunits
        ncpress.units = filtereddata.pressureunits

        
        ncF.close()
        print "Written " + outdir + '/'+ filename


# Class to write out the data from a stationary platform
filename = 0
outdir = 0

class write_ncdf_TT:
    def __init__(self, filtereddata, outdir = outdir, filename=filename):
        
        import os
        
        if type(outdir) == int:
            outdir = os.path.dirname(filtereddata.mzfile)
    
        if type(filename) == int:
            if type(filtereddata.mzfile) == tuple:            
                filetimestamp = filtereddata.mzfile[0][filtereddata.mzfile[0].find('h0')+ 2: filtereddata.mzfile[0].find('.nc')]
            if type(filtereddata.mzfile) == str:
                filetimestamp = filtereddata.mzfile[filtereddata.mzfile.find('h0')+ 2: filtereddata.mzfile.find('.nc')]
            filename = filtereddata.species + '_' + filtereddata.case + '.mzt.h0'+  filetimestamp + '_' + filtereddata.sitetype+'.nc'

        print 'writing file: ' + outdir + '/'+ filename
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(outdir + '/'+ filename, 'w')
        
        # Create the dimensions
        ncF.createDimension('time', len(filtereddata.time))
        ncF.createDimension('sitenames', len(filtereddata.sitenames))
        
        # Make some global attributes
        ncF.mzfile = filtereddata.mzfile
        ncF.sitefile = filtereddata.sitefile
        ncF.case = filtereddata.case
        ncF.species = filtereddata.species
        ncF.sitetype = filtereddata.sitetype
        
        # Create the variables
        #ncsitenames = ncF.createVariable('sitename','s', ('sitenames',))
        ncsitenames = ncF.createVariable('sitenames',type(filtereddata.sitenames[0]), ('sitenames',))        
 
        ncyear = ncF.createVariable('year', 'i', ('time',))
        ncmonth = ncF.createVariable('month', 'i', ('time',))
        ncday = ncF.createVariable('day', 'i', ('time',))
        nchour = ncF.createVariable('hour', 'i', ('time',))
        ncminute = ncF.createVariable('minute', 'i', ('time',))    
        
        nclon = ncF.createVariable('lon', 'f', ('sitenames',))
        nclat = ncF.createVariable('lat', 'f', ('sitenames',))
        
        ncalt = ncF.createVariable('alt', 'f', ('sitenames',))
        ncpress = ncF.createVariable('pressure', 'f', ('sitenames', 'time',))
        ncconc = ncF.createVariable('conc', 'f', ('sitenames', 'time',))
        ncconc_SD = ncF.createVariable('conc_sd', 'f', ('sitenames', 'time',))
        ncemiss = ncF.createVariable('emiss', 'f', ('sitenames', 'time',))
        
        # Fill the variables
        ncsitenames[:] = filtereddata.sitenames        
        
        ncyear[:] = [filtereddata.time[t].year for t in np.arange(len(filtereddata.time))]
        ncmonth[:] = [filtereddata.time[t].month for t in np.arange(len(filtereddata.time))]        
        ncday[:] = [filtereddata.time[t].day for t in np.arange(len(filtereddata.time))] 
        nchour[:] = [filtereddata.time[t].hour for t in np.arange(len(filtereddata.time))]
        ncminute[:] = [filtereddata.time[t].minute for t in np.arange(len(filtereddata.time))]        
   
        nclon[:] = filtereddata.lon
        nclat[:] = filtereddata.lat
        ncalt[:] = filtereddata.alt
        ncpress[:,:] = filtereddata.pressure
        ncconc[:,:] = filtereddata.conc
        ncconc_SD[:,:] = filtereddata.conc_sd
        ncemiss[:,:] = filtereddata.emis
        
        
        # Give the variables some attributes        
        nclon.units = 'Degrees east'
        nclat.units = 'Degrees north'
        ncconc.units = filtereddata.concunits
        ncemiss.units = filtereddata.emissunits
        ncpress.units = filtereddata.pressureunits

        
        ncF.close()
        print "Written " + outdir + '/'+ filename

# Class to read in the filtered output
class read_ncdf_TT:
    def __init__(self, filenames):
        
        # need to sort the files
        filenames.sort()
        
        for j in np.arange(len(filenames)):
            
            data=netCDF4.Dataset(filenames[j])
            
            year_j = data.variables['year'][:]
            month_j = data.variables['month'][:] 
            day_j = data.variables['day'][:] 
            hour_j = data.variables['hour'][:]
            minute_j = data.variables['minute'][:] 
            
            lon = data.variables['lon'][:]
            lat = data.variables['lat'][:]
            alt = data.variables['alt'][:]            
            
            sitenames = data.variables['sitenames'][:] 
            
            conc_j = np.transpose(data.variables['conc'][:])
            conc_SD_j = np.transpose(data.variables['conc_sd'][:])
            emis_j = np.transpose(data.variables['emiss'][:])
            pressure_j = np.transpose(data.variables['pressure'][:])
            
            dt_date_j = [dt.datetime(year_j[i],month_j[i],day_j[i],hour_j[i], minute_j[i]) for i in np.arange(len(year_j))]
    
            if j == 0:
            
                year = year_j
                month = month_j
                day = day_j
                hour = hour_j
                minute = minute_j
                conc = conc_j
                conc_SD = conc_SD_j                
                emis = emis_j
                pressure =  pressure_j
                dt_date = dt_date_j              
                
            else:
                
                year = np.concatenate((year, year_j))
                month = np.concatenate((month, month_j))
                day = np.concatenate((day, day_j))
                hour = np.concatenate((hour, hour_j))                
                minute = np.concatenate((minute, minute_j))
                conc = np.concatenate((conc, conc_j))
                conc_SD = np.concatenate((conc_SD, conc_SD_j))
                emis = np.concatenate((emis, emis_j))
                pressure = np.concatenate((pressure, pressure_j))                
                dt_date = np.concatenate((dt_date, dt_date_j))
      
                
        self.time = dt_date
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.emis = emis
        self.conc = conc
        self.conc_SD = conc_SD
        self.pressure = pressure
        self.alt = alt
        self.lon = lon
        self.lat = lat
        self.sitenames = sitenames
        self.filenames = filenames         
        self.species = str(data.__getattribute__('species')).strip()
        self.case = str(data.__getattribute__('case')).strip()
        self.sitetype = str(data.__getattribute__('sitetype')).strip()
        self.sitefile = str(data.__getattribute__('sitefile')).strip()
        self.concunits = data.variables['conc'].getncattr('units')
        self.emissunits = data.variables['emiss'].getncattr('units')
        self.pressureunits = data.variables['pressure'].getncattr('units')
        

