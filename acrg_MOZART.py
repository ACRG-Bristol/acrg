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
import bisect
import os

# ___________________________________________________________________________________________________________________
# CODE TO READ THE DIFFERENT DATA TYPES
# ___________________________________________________________________________________________________________________

 # Class to read in the data
class read:
    def __init__(self, filename):
    
        print 'Reading file : ' + filename
        
        if type(filename) == tuple:
            filename = filename[0]
        
        data=netCDF4.Dataset(filename)
        
        if 'h0' in filename:
            conc_varname = (str(data.__getattribute__('title'))).strip()+'_VMR_avrg'
        
        if 'h1' in filename:        
            conc_varname = (str(data.__getattribute__('title'))).strip()+'_13:30_LT'
          
        conc = data.variables[conc_varname][:]
        date = data.variables['date'][:] # Date YYYYMMDD
        secs = data.variables['datesec'][:] # Seconds to be added to above date to make date-time variable
        lon = data.variables['lon'][:].astype('float')
        lat = data.variables['lat'][:].astype('float')
        lev = data.variables['lev'][:].astype('int')
        PS=data.variables['PS'][:].astype('float')
        P0=data.variables['P0'][:].astype('float')
        hyai=data.variables['hyai'][:].astype('float')
        hybi=data.variables['hybi'][:].astype('float')
        
        
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
            print 'date ' + str(date[i]) + ' processed'
         
        self.time = time_t
        self.conc = conc
        self.lon = lon
        self.lat = lat
        self.lev = lev
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
        self.pressureunits = data.variables['P0'].getncattr('units')
        
        if 'h0' in filename:
            emiss_varname = (str(data.__getattribute__('title'))).strip()+'_SRF_EMIS_avrg'
            emis = data.variables[emiss_varname][:]
            self.emis = emis             
            self.emissunits = data.variables[emiss_varname].getncattr('units')
             
             
 # Class to read in the site file GONZI netcdf version       
class read_sitefile_GONZI_nc:
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
      
        
        # Extract the variable names
        keys = np.array((data.variables.keys()))
       
        
        # Extract the data   
        lat = data.variables['LAT'][:]
        lon = data.variables['LONG'][:]
        if ferry_flag == 0:
            alt = data.variables['ALTI'][:] # Date YYYYMMDD
            alt_units = data.variables['ALTI'].units
        if ferry_flag == 1:
             alt = np.empty(len(lat))
             alt_units = ''
             
        Year = data.variables['YYYY'][:] # Seconds to be added to above date to make date-time variable
        Month = data.variables['MM'][:]
        Day = data.variables['DD'][:]
        SecSinceMidnight=data.variables['Time_Slots'][:]
        
        CH4 = data.variables['CH4'][:]        
        CO2 = data.variables['CO2'][:]     
        
        # only do this if N2O is given
        if len(np.where(keys == 'N2O')[0]) > 0:
            N2O = data.variables['N2O'][:]          
        else:
            N2O = np.empty((len(CO2),3))
            N2O[:] = 'nan'
        
        model_CH4 = data.variables['Model_CH4'][:]
        model_CO2 = data.variables['Model_CO2'][:]        
        model_N2O = data.variables['Model_N2O'][:]
        
        # Create the time variable
        time_t = [dt.datetime(Year[i],Month[i],Day[i]) + dt.timedelta(seconds=(SecSinceMidnight[i]).astype('int')) for i in np.arange(len(Day))]

        # Remove the superfluous points
        good_i = (np.where(lat[:,0] != -9999))[0]
        
        self.time = [time_t[good_i[i]] for i in np.arange(len(good_i))]
        self.lat = lat[good_i]
        self.lon = lon[good_i]
        self.alt = alt[good_i]
        self.filename = sitefile
        self.conc = np.squeeze(np.array([[CH4[good_i]],[CO2[good_i]],[N2O[good_i]]]))
        self.species_tags = np.array(['CH4','CO2','N2O'])
        self.good_index = good_i
        self.alt_units = alt_units
        self.model_CH4 = model_CH4[good_i]
        self.model_CO2 = model_CO2[good_i]
        self.model_N2O= model_N2O[good_i]


 # Class to read in the netcdf site file written by matt listing the fixed sites    
class read_fixed_sitefile_nc:
    def __init__(self, sitefile = 0, species = 'CH4', dir = '/data/shared/GAUGE/'):
        
        if type(sitefile) == int:
            sitefile = dir + species + '/mozart_obs_stationary.nc'
        
        print 'Using site file : ' + sitefile
        
        data=netCDF4.Dataset(sitefile, 'r')

        if species == 'CH4':
            species_lc = 'ch4'
        
        if species == 'CO2':
            species_lc = 'co2'
        
        if species == 'N2O':
            species_lc = 'n2o'
        
        
        # Extract the data    
        time = data.variables['time'][:] # "seconds since 2004-01-01 00:00:00" ;
        site = data.variables['site'][:]
        network = data.variables['network'][:]
        lon = data.variables['longitude'][:]
        lat = data.variables['latitude'][:]
        alt = data.variables['altitude'] # "m  above 1.9x2.5 MOZART surface level (m)" ;
        conc = data.variables[species_lc]
        repeatability = data.variables[species_lc+'_repeatability']
	      
        # Create the time variable
        dateunits = data.variables['time'].getncattr('units')
        sincedate = dateunits[dateunits.find('seconds since ') +14:-1] # seconds since 2014-01-01 00:00:00
                
        time_dt = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in time]
        #time_dt = [dt.datetime(2009,1,1,0,0,0) + dt.timedelta(seconds=(i).astype('int')) for i in time]
        

        self.time = time_dt
        self.lat = lat
        self.lon = lon
        self.alt = alt[:]
        self.alt_units = alt.units
        self.site = site
        self.species = species
        self.species_lc = species_lc
        self.conc = conc[:]
        self.repeatability = repeatability[:]
        self.units = conc.units
        self.network = network
        self.sitefile = sitefile


 # Class to read in the netcdf site file written by matt listing the column sites (satelite + TCON)  
class read_column_sitefile_nc:
    def __init__(self, sitefile = 0, species = 'CH4', dir = '/data/shared/GAUGE/mozart_obs_column/', month = 1, year = 2009):
        
        if type(sitefile) == int:
            sitefile = dir + species + '/mozart_obs_column' + str(month).zfill(2) +str(year)+ '.nc'
        
        data=netCDF4.Dataset(sitefile, 'r')

        if species == 'CH4':
            species_lc = 'ch4'
        
        if species == 'CO2':
            species_lc = 'co2'
        
        if species == 'N2O':
            species_lc = 'n2o'
        
        # Extract the data    
        time = data.variables['time'][:] # "seconds since 2004-01-01 00:00:00" ;
        network = data.variables['network'][:]
        lon = data.variables['longitude'][:]
        lat = data.variables['latitude'][:]
        gas_data = data.variables[species_lc]
        gas_repeatability = data.variables[species_lc+'_repeatability']
        averaging_kernel = data.variables['averaging_kernel'][:]
        pressure = data.variables['pressure']
       
        # Create the time variable
        dateunits = data.variables['time'].getncattr('units')
        sincedate = dateunits[dateunits.find('seconds since ') +14:-1] # seconds since 2009-01-01 00:00:00
                
        time_dt = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in time]
        #time_dt = [dt.datetime(2009,1,1,0,0,0) + dt.timedelta(seconds=(i).astype('int')) for i in time]

        self.time = time_dt
        self.time_secs = time
        self.lat = lat
        self.lon = lon
        self.no_lev = np.shape(pressure[:])[1]
        self.gasname = species_lc
        self.conc = gas_data[:]
        self.repeatability = gas_repeatability[:]
        self.units = gas_data.units
        self.averaging_kernel = averaging_kernel
        self.pressure = pressure[:]
        self.pressure_units = pressure.units
        self.network = network
        self.sitefile = sitefile
 

# Class to read in the netcdf site file written by matt listing the mobile GAUGE sites (ferry and aircraft)  
class read_mobile_sitefile_nc:
    def __init__(self, sitefile = 0, species = 'CH4', dir = '/data/shared/GAUGE/', month = 1, year = 2003):
        
        if type(sitefile) == type(0):
            sitefile = dir + species + '/mozart_obs_mobile/mozart_obs_mobile_' +str(year)+ str(month).zfill(2)+'.nc'
        
        # Check if the site file exists
        exists = os.path.isfile(sitefile)
                
        if exists :
            
            print 'Site file being read:'
            print sitefile
            
            data=netCDF4.Dataset(sitefile, 'r')
    
            if species == 'CH4':
                species_lc = 'ch4'
            
            if species == 'CO2':
                species_lc = 'co2'
            
            if species == 'N2O':
                species_lc = 'n2o'
            
            # Extract the data    
            time = data.variables['time'][:] # "seconds since 2004-01-01 00:00:00" ;
            site = data.variables['site'][:]
            network = data.variables['network'][:]
            scale = data.variables['scale'][:]
            lon = data.variables['longitude'][:]
            lat = data.variables['latitude'][:]
            alt = data.variables['altitude']
            gas_data = data.variables[species_lc]
            gas_repeatability = data.variables[species_lc+'_repeatability']
            pressure = data.variables['pressure']
            
            # Create the time variable
            dateunits = data.variables['time'].getncattr('units')
            sincedate = dateunits[dateunits.find('seconds since ') +14:-1] # seconds since 2014-01-01 00:00:00
                    
            time_dt = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in time]
            #time_dt = [dt.datetime(2009,1,1,0,0,0) + dt.timedelta(seconds=(i).astype('int')) for i in time]
    
            self.time = time_dt
            self.time_secs = time
            self.lat = lat
            self.lon = lon
            self.alt = alt[:]
            
            self.alt_units = alt.units
            self.scale = scale
            
            self.gasname = species_lc
            self.conc = gas_data[:]
            self.repeatability = gas_repeatability[:]
            self.units = gas_data.units
            self.pressure = pressure[:]
            self.pressure_units = pressure.units
            self.network = network
            self.site = site
            self.sitefile = sitefile
            self.fileexists = exists
    
        else:
            print 'Sitefile ' + sitefile + ' does not exist'
            self.fileexists = exists




# ___________________________________________________________________________________________________________________
# CODE TO DO USEFUL STUFF
# ___________________________________________________________________________________________________________________
        
 # Class to estimate the pressure levels from a given altitude
 # This uses the scale height and a simple atmospheric model
 # P = P0 * e^(-Z/H)
 # Where P = Pressure
 # P0 = surface pressure in Pa
 # -Z = altitude (km)
 # H = scale height (km) this is hard coded to be 7.64km which is the global mean
 # This assumes that the altitude is in m unless other units are given
class calc_pressure:
    def __init__(self, altitude, P0, units='m'):

        from math import exp
        
        if units == 'm':
            altitude = altitude/1000.0
        elif units =='metres':
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
            sitedata=read_sitefile_GONZI_nc(sitefile)
            
        site_i = (np.where(sitedata.sitenames==sitename))[0]      
        
        self.sitename = sitedata.sitenames[site_i]        
        self.lat = sitedata.lat[site_i]
        self.lon = sitedata.lon[site_i]
        self.alt = sitedata.alt[site_i]

# Class to match to the closeest lat/lon using bisect 
class match_latlon:
    def __init__(self, lat, lon, lat_array, lon_array):
        
        import bisect
        
        # NB: Assuming evenly spaced grid
        lat_spacing = lat_array[1] - lat_array[0]
        lon_spacing = lon_array[1] - lon_array[0]
        
        
        # If the lon array is from 0 to 360
        # then need to check if the lon point is < 0
        # if it's < 0 then convert it
        if min(lon_array) >= 0 and lon < 0 :
            lon = lon + 360
        
        # If the lon array is from -180 to 180
        # then need to check that the lon point is > 180
        # if it's > 180 then convert it
        if min(lon_array) >= 0 and lon < 0 :
            lon = 360 - lon
        
        
        lat_index = bisect.bisect(lat_array+(lat_spacing/2), lat) # nb: the adjustment shifts the grid so that it returns the closest point
        lon_index = bisect.bisect(lon_array+(lon_spacing/2), lon) # nb: the adjustment shifts the grid so that it returns the closest point


        # Check the edges of the lon grid
        if lon_index == len(lon_array):
            lon_index = 0
            
        
        self.inputlocation = np.array((lat,lon))
        self.lat_array = lat_array
        self.lon_array = lon_array
        self.closestpoint = np.array((lat_array[lat_index], lon_array[lon_index]))
        self.closestindex = np.array((lat_index, lon_index))
        








# ___________________________________________________________________________________________________________________
# CODE TO MATCH THE DIFFERENT DATA TYPES
# ___________________________________________________________________________________________________________________

        
# Class to filter the data
# This uses an individual MOZART history file
# and matt's netcdf site file which contains the lat, lon and alt for fixed sites
class data_filter_fixed:
    def __init__(self, mzfile, sitefile, Ave=0):
        
     
        # Read MOZART file name
        data = read(mzfile)
        
        
        # Read site info file
        # siteinfo = read_sitefile_txt(sitefile)
        species = data.species
        siteinfo = read_fixed_sitefile_nc(species=species)    
        
        # create parameters to store output
        # initialise to nans
        pressures = np.empty((len(siteinfo.site), len(data.time), 3, 3, 3))*np.nan  
        concs = np.empty((len(siteinfo.site), len(data.time), 3, 3, 3))*np.nan 
        concs_sd = np.empty((len(siteinfo.site), len(data.time), 3, 3, 3))*np.nan 
        
        levs = np.empty((len(siteinfo.site),len(data.time),3))*np.nan 
        
        emissions = np.empty((len(siteinfo.site), len(data.time), 3, 3))*np.nan 
        
        matched_lats = np.empty((len(siteinfo.site), 3))*np.nan 
        matched_lons = np.empty((len(siteinfo.site), 3))*np.nan 
        matched_levs = np.empty((len(siteinfo.site), len(data.time), 3))*np.nan 
        
        site_pressure = np.empty(len(siteinfo.site))*np.nan         
        
        
        # loop through each lat/lon in the siteinfo
        for j in np.arange(len(siteinfo.lat)):
            
            
            # use the match_latlot code which is based on bisect rather than calculating the exact distance
            latlon_index = match_latlon(siteinfo.lat[j],siteinfo.lon[j],data.lat, data.lon)      
            
            # As it's tall tower then we'll want the data for all time stamps
            # Extract column pressure at the correct lat/lon for all timestamps
            # This should be a time by lev array
            # column_P = np.squeeze(data.pressure[:,[np.arange(len(data.hyai)-1)],haversine.mindist_index[0], haversine.mindist_index[1]])        

            # Changed 16/1/2015 to extract
            # the closest match and then a cube surrounding that point, 
            # - north, south, east, west, above and below

            # find the lat and lon range
            print 'site: ' + siteinfo.site[j] 
            print 'site postion: ' + str(siteinfo.lat[j]) + ', ' + str(siteinfo.lon[j])
              
            
            # Adjustments for hitting the edge of the grid
            # lat index = 0 i.e. south pole
            if latlon_index.closestpoint[0] == 0:
                lat_range = np.squeeze([latlon_index.closestindex[0], latlon_index.closestindex[0]+1])
                matched_lats[j,0:2] = data.lat[lat_range]
                
                
            # lat index  = 95 i.e. north pole
            elif latlon_index.closestpoint[0] == 95:
                lat_range = np.squeeze([latlon_index.closestindex[0]-1, latlon_index.closestindex[0]])
                matched_lats[j,0:2] = data.lat[lat_range]
                
            else:            
                lat_range = np.squeeze([latlon_index.closestindex[0]-1, latlon_index.closestindex[0], latlon_index.closestindex[0]+1])
                matched_lats[j,:] = data.lat[lat_range]
                

            
            # Adjustments for hitting the edge of the grid
            # lon index = 0 i.e. weastern edge
            if latlon_index.closestindex[1] == 0:
                lon_range = [143,0,1]
            
            # lon index  = 143 i.e. eastern edge
            elif latlon_index.closestindex[1] == 143:
                lon_range = [142,143,0]
                
            else:            
                lon_range = np.squeeze([latlon_index.closestindex[1]-1, latlon_index.closestindex[1], latlon_index.closestindex[1]+1])
            
            matched_lons[j,:] = data.lon[lon_range]
            
            
            print 'model lats: ' + str(matched_lats[j,:])
            print 'model lons: ' + str(matched_lons[j,:])

            

            
            # Extract the column pressure at each time point for the matching (central) point of the cube
            # Pressure is time x lev x lat x lon
            column_P = np.squeeze(data.pressure[:,:,latlon_index.closestindex[0], latlon_index.closestindex[1]])        
            
            
            
            if Ave != 0:
                # Just take the average of the bottom 7 levels that's all the data roughly < 1 km
    
                # Put data for each lat/lon (i.e. j value) into output arrays 
                # I want the lowest 7 levels which are actually the LAST 7 not the first
                index = np.empty(7)
                index[:] = len(data.hyai)-1
                index = index - np.arange(7) -1
                pressures[j,:] = np.mean(np.squeeze(column_P[:,[index]]), axis=1)
                
                # Extract the corresponding concentrations and emissions
                conc_j = np.squeeze(data.conc[:,:, latlon_index.closestindex[0], latlon_index.closestindex[1]])
                concs[j,:] = np.mean(np.squeeze(conc_j[:,[index]]), axis=1)
                concs_sd[j,:] = np.std(np.squeeze(conc_j[:,[index]]), axis=1)
                
    
            else:
                
                
                
                # Convert the altitude to a pressure using P0 from the model
                site_pressure[j] = (calc_pressure(siteinfo.alt[j], data.P0, units=siteinfo.alt_units)).pressure             
                print 'site pressure: ' + str(site_pressure[j])               
                
                # Loop through each time step and match site pressure to column level
                # Extract the cube of data
                # Store the positions
                for i in np.arange(len(data.time)):
                    # Find the index of the closest level based on pressure
                    # column_P is time x lev
                    lev_index = np.where(abs(column_P[i,:] - site_pressure[j]) == min(abs(column_P[i,:] - site_pressure[j])))[0]
                    
                    # As I'm using "fancy' indexing rather than slicing I need to do each dimension separately 
                    # Extract data for time = i
                    pressures_i = np.squeeze(data.pressure[i,:,:,:])
                    concs_i = np.squeeze(data.conc[i,:,:,:])
                    emissions_i = np.squeeze(data.emis[i,:,:])
                    
                    # Extract data for lat range
                    pressures_lat = pressures_i[:,lat_range,:]
                    concs_lat = concs_i[:,lat_range,:]
                    emissions_lat = emissions_i[lat_range,:]
                    
                    # Extract data for lon range  
                    pressures_lon = pressures_lat[:,:,lon_range]
                    concs_lon = concs_lat[:,:,lon_range]
                    emissions_lon = emissions_lat[:,lon_range]
                    
                    #pdb.set_trace()
                    
                    # Define the range of the height of the cube
                    # NB: the levels are 1 based while the levels are 1 based
                    # if using the bottom level
                    if lev_index == 55:
                        lev_i = np.squeeze([54,55])
                                                
                        # store the level numbers used                    
                        levs[j,i,0:2] =  [55, 56]                  
                        
                        #print 'lev_i: ' + str(lev_i)                        
                        
                        # extract and store the pressure levels used and corresponding concs
                        if len(lat_range) == 2:
                            pressures[j,i,0:2,0:2,:] = pressures_lon[lev_i,:,:]
                            concs[j,i,0:2,0:2,:] = concs_lon[lev_i,:,:]
                        else:
                            pressures[j,i,0:2,:,:] = pressures_lon[lev_i,:,:]
                            concs[j,i,0:2,:,:] = concs_lon[lev_i,:,:]
                    
                        
                    else:
                        lev_i = np.squeeze([lev_index-1, lev_index, lev_index+1])
                    
                        
                        # store the level numbers used  
                        # NB: levels are 1 based while the index is 0 based
                        levs[j,i,:] =  [lev_index, lev_index+1, lev_index+2]                    
                        
                        #print 'lev_i: ' + str(lev_i)
                        
                        # extract and store the pressure levels used and corresponding concs
                        if len(lat_range) == 2:
                            pressures[j,i,:,0:2,:] = np.squeeze(pressures_lon[lev_i,:,:])
                            concs[j,i,:,0:2,:]  = np.squeeze(concs_lon[lev_i,:,:])
                        else:
                            pressures[j,i,:,:,:] = np.squeeze(pressures_lon[lev_i,:,:])
                            concs[j,i,:,:,:]  = np.squeeze(concs_lon[lev_i,:,:])
                    
                    
                if len(lat_range) == 2 :
                    emissions[j,:,0:2,:] = emissions_lon
                else:
                    emissions[j,:,:,:] = emissions_lon
                       
               
                    
        self.site = siteinfo.site
        self.sitetype = 'TT'
        self.sitenames = siteinfo.site
        self.site_lat = siteinfo.lat
        self.site_lon = siteinfo.lon
        self.site_alt = siteinfo.alt
        self.site_pressure = site_pressure
        self.time = data.time
        
        self.model_pressure = np.squeeze(pressures)
        self.model_lat = matched_lats
        self.model_lon = matched_lons
        self.model_levs = matched_levs
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
# and Matt's netcdf site file which contains the lat, lon and alt for column sites
class data_filter_column:
    def __init__(self, mzfile, sitefile, month=1, year=2009):
        
        # Read MOZART file name
        data = read(mzfile)
        
        
        # Read site info file
        species = data.species
        columndata = read_column_sitefile_nc(species=species, sitefile=sitefile, month=month, year=year)
        
        # Only want to look at time points that are within the range of the data file
        # Convert model time to seconds since 1/1/2009 00:00:00
        # compare that to Matt's timestamps which are already seconds since 1/1/2009 00:00:00
        model_secs = [(i - dt.datetime(2009,1,1,0,0,0)).total_seconds() for i in data.time]
        # model_secs = [(i - dt.datetime(2003,1,1,0,0,0)).total_seconds() for i in data.time]
        
        good_index = np.where( (columndata.time_secs >= model_secs[0]) & (columndata.time_secs <= model_secs[-1]) )[0]
        
        if len(good_index) != 0:
            
            # create parameters to store output
            # initialise to nans
            model_pressures = np.empty((len(columndata.time), len(data.lev)))*np.nan  
            
            print "model pressures shape"
            print np.shape(model_pressures)
            
            model_concs = np.empty((len(good_index), len(data.lev)))*np.nan 
            model_times = []
            matched_lats = np.empty(len(good_index))*np.nan 
            matched_lons = np.empty(len(good_index))*np.nan 
            
            
            # loop through each time and lat/lon in the siteinfo
            for j in np.arange(len(good_index)):
                
                # Use the haverseine code to determine the closest surface grid point to the given location
                # haversine =  haver.multipledistances([columndata.lat[good_index[j]], columndata.lon[good_index[j]]] , data.lat, data.lon)     
                
                # Use the match_latlon which is based on bisect
                latlon_index = match_latlon(columndata.lat[j],columndata.lon[j],data.lat, data.lon)      

                    
                # Match the column timestamp to the closest model timestamp
                time_index = np.where( abs(model_secs - columndata.time_secs[good_index[j]]) == np.nanmin(abs(model_secs - columndata.time_secs[good_index[j]] ) ) )[0]
                
                # Sometimes can be halfway between two timestamps default to the first one
                if len(time_index) > 1:
                    time_index = time_index[0]
                                
                print 'time index: ' + str(time_index)
    
                print 'j: ' + str(j)
                
                print 'column time: ' + str(columndata.time[good_index[j]])
                print 'model time: ' + str(data.time[time_index])
                print 'column postion: ' + str(columndata.lat[good_index[j]]) + ', ' + str(columndata.lon[good_index[j]])
                  
                
                matched_lats[j] = data.lat[latlon_index.closestindex[0]]
                matched_lons[j] = data.lon[latlon_index.closestindex[1]]
                
                
                print 'model lats: ' + str(matched_lats[j])
                print 'model lons: ' + str(matched_lons[j])
    
                
                
                # Extract the column pressure at the matching time point for the matching lat/lon
                # Pressure is time x lev x lat x lon
                model_times.append(data.time[time_index])         
                model_pressures[j,:] = np.squeeze(data.pressure[time_index,:, latlon_index.closestindex[0], latlon_index.closestindex[1]])        
                model_concs[j,:] = np.squeeze(data.conc[time_index,:, latlon_index.closestindex[0], latlon_index.closestindex[1]])
        
                            
           
                        
            self.sitetype = 'Column'
            self.site_lat = columndata.lat[good_index]
            self.site_lon = columndata.lon[good_index]
            self.site_pressure = columndata.pressure[good_index]
            self.site_time = [columndata.time[k] for k in good_index]
            
            self.model_time =np.array(model_times)
            self.model_pressure = np.squeeze(model_pressures)
            self.model_lat = matched_lats
            self.model_lon = matched_lons
            self.conc = np.squeeze(model_concs)
            
            self.case = data.case.strip()
            self.species = species
            self.mzfile = data.filename.strip()
            self.sitefile = sitefile        
            self.concunits = data.concunits
            self.pressureunits = data.pressureunits    
            self.good_index = good_index
        
        else:
            self.good_index = 0                             



# Class to filter the data
# This uses an individual MOZART history file
# and a site file which contains the lat, lon, alt and time for moving sites
class data_filter_mobile:
    def __init__(self, mzfile, sitefile = '/data/shared/GAUGE/CH4/mozart_obs_mobile/mozart_obs_mobile_200301.nc'):   
                
    
        # Read MOZART file name
        data = read(mzfile)
        
        species = data.species
        #print 'reading site file : ' 
        print sitefile
        
        # Read site info file        
        siteinfo = read_mobile_sitefile_nc(species = species, month = data.time[0].month, year = data.time[0].year)
        
        # only proceed if the site file exists
        if siteinfo.fileexists :
            
            matched_data = data_match_mobile(species = species, \
                                            model_conc = data.conc, \
                                            model_lat = data.lat, \
                                            model_lon = data.lon, \
                                            model_time = data.time, \
                                            model_P0 = data.P0, \
                                            model_pressure = data.pressure, \
                                            model_emission = 0, \
                                            obs_conc = siteinfo.conc, \
                                            obs_lat = siteinfo.lat, \
                                            obs_lon = siteinfo.lon, \
                                            obs_time = siteinfo.time, \
                                            obs_pressure = siteinfo.pressure, \
                                            obs_alt = 0, \
                                            obs_alt_units = 0, \
                                            quiet = 1)
         
                
            self.network = siteinfo.network
            self.site = siteinfo.site
            self.sitetype = 'mobile'
            self.case = data.case.strip()
            self.species = species
            self.mzfile = mzfile
            self.sitefile = sitefile
            self.concunits = data.concunits
            self.pressureunits = data.pressureunits
            
            self.site_lat = siteinfo.lat
            self.site_lon = siteinfo.lon
            self.site_alt = siteinfo.alt
            self.site_pressure = siteinfo.pressure
            self.site_conc = siteinfo.conc
            self.site_repeat = siteinfo.repeatability
            self.site_units = siteinfo.units
            self.site_time = siteinfo.time
            
            self.model_pressure = matched_data.matched_pressure
            self.model_lat = matched_data.matched_lat
            self.model_lon = matched_data.matched_lon
            self.model_time = matched_data.matched_time
            
            self.conc = matched_data.matched_conc
            
            self.fileexists = True

            
        else:
            print 'There is no matching mobile obs data for : ' + mzfile
            self.fileexists = False
 

    
# Generalised code to match model 4D array (time/lat/lon/P) to given obs vector (dim = time) with associated lat, lons and P (alt) 
# Assumes that the model and obs times are given as datetimes
# Set quiet = 0 to stop matched info being printed
class data_match_mobile:
    def __init__(self, species = 0, \
                model_conc = 0, \
                model_lat = 0, \
                model_lon = 0, \
                model_time = 0, \
                model_P0 = 0, \
                model_pressure = 0, \
                model_emission = 0, \
                obs_conc = 0, \
                obs_lat = 0, \
                obs_lon = 0, \
                obs_time = 0, \
                obs_pressure = 0, \
                obs_alt = 0, \
                obs_alt_units = 0, \
                quiet = 1):   
                
            # check if the matching will be using alt or P
            if type(obs_pressure) == type(0):
                if type(obs_alt) == type(0):
                    print 'Not obs pressure or altitude were given'
                    print 'Assuming ground level'
 
 
            # Put in checks of the inputs
            if type(model_time) == type(0):
                print 'You need to give model time'
            if type(model_conc) == type(0):
                print 'You need to give model concentrations'
            if type(model_lat) == type(0):
                print 'You need to give model lat'
            if type(model_lon) == type(0):
                print 'You need to give model lon'
            if type(model_pressure) == type(0):
                print 'You need to give model pressure'
            if type(model_P0) == type(0):
                print 'You need to give model P0'
                
            if type(obs_time) == type(0):
                print 'You need to give obs time'
            if type(obs_conc) == type(0):
                print 'You need to give obs conc'
            if type(obs_lat) == type(0):
                print 'You need to give obs lat'
            if type(obs_lon) == type(0):
                print 'You need to give obs lon'
            


             
            # Convert the model times and the obs times to seconds
            # Assumes that these are given as datetimes
            
            model_secs = np.asarray([(i - dt.datetime(2009,1,1,0,0,0)).total_seconds() for i in model_time])
            
            obs_secs = np.asarray([(i - dt.datetime(2009,1,1,0,0,0)).total_seconds() for i in obs_time])

                      
 
            # create parameters to store output
            matched_pressure = np.empty(len(obs_time))      
            matched_conc = np.empty(len(obs_time))
            matched_time = []
            matched_lat = np.empty(len(obs_time))
            matched_lon = np.empty(len(obs_time))
            matched_emissions = np.empty(len(obs_time))
            
            # loop through each lat/lon in the same time range as the model output
            for j in np.arange(len(obs_time)):
                
                # Use the haverseine code to determine the closest model grid point to the obs location
                latlon_index = (match_latlon(obs_lat[j],obs_lon[j],model_lat, model_lon)).closestindex      
                matched_lat[j] = model_lat[latlon_index[0]]                
                matched_lon[j] = model_lon[latlon_index[1]] 
                
                
                # Match obs timestamp to the closest model timestamp
                time_gap = model_secs[1] - model_secs[0]
                timeindex_j = bisect.bisect(model_secs + (time_gap/2), obs_secs[j])
                
                
                # if it's AFTER the last model point then bisect returns the number of elements in the data.time array
                if timeindex_j == len(model_time):
                    timeindex_j = len(model_time) - 1
                    print 'Switched to last element of time array'
                
                matched_time.append(model_time[timeindex_j])
                
                
                # Determine the closest level
                # if there's no alt or pressure given assuming it's the surface (which for mozart is the last one)
                print j                
                if type(obs_pressure) == type(0) and type(obs_alt) == type(0):
                    
                    lev_i = -1
                
                else:
                    # if obs pressure isn't given use alt
                    if type(obs_pressure) == type(0):

                        obs_pressure_j = calc_pressure(obs_alt[j], model_P0, units=obs_alt_units).pressure                    
                    
                    else:
                        
                        obs_pressure_j = obs_pressure[j]
                    
                    # Extract column pressure at the correct lat/lon for the  matching timestamp
                    
                    column_P = np.squeeze(model_pressure[timeindex_j,:,latlon_index[0], latlon_index[1]])
                    
                    # Match site pressure to column level
                    # as the pressure levels aren't evenly matched i'd need to find the gap and then compare to 
                    # the pressures  on either side which is unlikely to be much faster than this anyway
                    lev_i = np.where(abs(column_P - obs_pressure_j) == min(abs(column_P - obs_pressure_j)))[0]
             
                    matched_pressure[j] = np.squeeze(column_P[lev_i])
                     
                # Extract the corresponding concentrations and emission
                # Put data for each lat/lon (i.e. j value) into output arrays            
                
                matched_conc[j] = np.squeeze(model_conc[timeindex_j,lev_i, latlon_index[0],latlon_index[1]])
                
                if type(model_emission) != type(1):
                    matched_emissions[j] = np.squeeze(model_emission[timeindex_j,lev_i, latlon_index[0],latlon_index[1]])
    
    
                if quiet != 0 :
                    np.shape(obs_lat)
                    type(obs_lat)
                    print 'Obs lat: ' + str(obs_lat[j])
                    print 'Model lat: ' + str(matched_lat[j])
                    print 'Obs lon: ' + str(obs_lon[j])
                    print 'Model lon: ' + str(matched_lon[j])
                       
                    print 'Obs time: ' + str(obs_time[j])
                    print 'Model time: ' + str(matched_time[j])
                    
                    if lev_i != -1 :
                        print 'Obs pressure: ' + str(obs_pressure_j)
                    print 'Model pressure: ' + str(matched_pressure[j])
                    
            self.obs_time = obs_time    
            self.obs_lat = obs_lat
            self.obs_lon = obs_lon
            self.obs_alt = obs_alt
            self.obs_pressure = obs_pressure
            self.obs_conc = obs_conc

            self.model_time = model_time
            self.model_lat = model_lat
            self.model_lon = model_lon
            self.model_pressure = model_pressure
            self.model_conc = model_conc
            
            self.matched_time = matched_time
            self.matched_lat = matched_lat
            self.matched_lon = matched_lon
            self.matched_pressure = matched_pressure
            self.matched_conc = matched_conc
            self.matched_emissions = matched_emissions
            
            self.species = species
          
 
 
 
 
# ___________________________________________________________________________________________________________________
# CODE TO WRITE THE DIFFERENT DATA TYPES
# ___________________________________________________________________________________________________________________
            
            
# Class to write out the data from a moving platform
# Uses the h0 file 
class write_ncdf_mobile:
    def __init__(self, filtereddata, outdir = 0, filename=0, mzfile=0):
        
        import os
        
        if mzfile == 0:
            mzfile = filtereddata.mzfile
        
        if type(outdir) == int:
            outdir = os.path.dirname(mzfile)
    
        
        if type(filename) == int:
            if type(mzfile) == tuple:  
                filetimestamp = (filtereddata.mzfile[0])[(filtereddata.mzfile[0]).find('h0')+ 2: (filtereddata.mzfile[0]).find('.nc')]
            if type(mzfile) == str:
                filetimestamp = filtereddata.mzfile[filtereddata.mzfile.find('h0')+ 2: filtereddata.mzfile.find('.nc')]
            filename = filtereddata.species + '_' + filtereddata.case + '.mzt.h0'+  filetimestamp + '_' + filtereddata.sitetype+'.nc'
        
        print 'writing file: ' + outdir + '/'+ filename
        
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(outdir + '/'+ filename, 'w')
        
        # Create the dimensions
        ncF.createDimension('time', len(filtereddata.site_time))
        
        # Make some global attributes
        ncF.mzfile = filtereddata.mzfile
        ncF.sitefile = filtereddata.sitefile
        ncF.case = filtereddata.case
        ncF.species = filtereddata.species
        ncF.sitetype = filtereddata.sitetype
        
        # Create the variables
        ncsitenames = ncF.createVariable('sitename',type(filtereddata.site[0]), ('time',))        
        
        ncobs_time = ncF.createVariable('obs_time', 'i', ('time',))
        ncobs_lon = ncF.createVariable('obs_lon', 'f', ('time',))
        ncobs_lat = ncF.createVariable('obs_lat', 'f', ('time',))
        ncobs_alt = ncF.createVariable('obs_alt', 'f', ('time',))
        ncobs_press = ncF.createVariable('obs_pressure', 'f', ('time',))
        ncobs_conc = ncF.createVariable('obs_conc', 'f', ('time',))
        ncobs_repeat = ncF.createVariable('obs_repeatability', 'f', ('time',))
        
        ncmodel_time = ncF.createVariable('model_time', 'i', ('time',))
        ncmodel_lon = ncF.createVariable('model_lon', 'f', ('time',))
        ncmodel_lat = ncF.createVariable('model_lat', 'f', ('time',))
        ncmodel_press = ncF.createVariable('model_pressure', 'f', ('time',))

        ncmodel_conc = ncF.createVariable('model_conc', 'f', ('time'))
        
        
        # Fill the variables
        ncsitenames[:] = filtereddata.site        
                
        # times as seconds since 1/1/2009
        ncobs_time[:] = [(t - dt.datetime(2009,1,1,0,0,0)).total_seconds() for t in filtereddata.site_time]     
    
        ncobs_lon[:] = filtereddata.site_lon
        ncobs_lat[:] = filtereddata.site_lat
        ncobs_alt[:] = filtereddata.site_alt
        ncobs_press[:] = filtereddata.site_pressure
        ncobs_conc[:] = filtereddata.site_conc
        ncobs_repeat[:] = filtereddata.site_repeat
        
        # times as seconds since 1/1/2009
        ncmodel_time[:] = [(t - dt.datetime(2009,1,1,0,0,0)).total_seconds() for t in filtereddata.model_time]     
    
        ncmodel_lon[:] = filtereddata.model_lon
        ncmodel_lat[:] = filtereddata.model_lat
        ncmodel_press[:] = filtereddata.model_pressure

        ncmodel_conc[:] = filtereddata.conc
        
        
        # Give the variables some attributes        
        ncmodel_lon.units = 'Degrees east'
        ncmodel_lat.units = 'Degrees north'
        ncmodel_conc.units = filtereddata.concunits
        ncobs_conc.units = filtereddata.site_units
        ncmodel_press.units = filtereddata.pressureunits
        ncmodel_time.units = 'seconds since 2009-01-01 00:00:00'
        
        ncF.close()
        print "Written " + outdir + '/'+ filename



# Class to write out the data from a stationary platform
class write_ncdf_fixed:
    def __init__(self, filtereddata, outdir = 0, filename=0):
        
        import os
        
        mzfile = filtereddata.mzfile
        
        if type(outdir) == int:
            outdir = os.path.dirname(mzfile)
    
        if type(filename) == int:
            if type(mzfile) == tuple:  
                filetimestamp = mzfile[0][mzfile[0].find('h0')+ 2: mzfile[0].find('.nc')]
            if type(mzfile) == str:
                filetimestamp = mzfile[mzfile.find('h0')+ 2: mzfile.find('.nc')]
            filename = filtereddata.species + '_' + filtereddata.case + '.mzt.h0'+  filetimestamp + '_' + filtereddata.sitetype+'.nc'
        
        print 'writing file: ' + outdir + '/'+ filename
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(outdir + '/'+ filename, 'w')
        
        # Create the dimensions
        ncF.createDimension('time', len(filtereddata.time))
        ncF.createDimension('sitenames', len(filtereddata.site))
        ncF.createDimension('boxwidth', 3)
        
        # Make some global attributes
        ncF.mzfile = filtereddata.mzfile
        ncF.sitefile = filtereddata.sitefile
        ncF.case = filtereddata.case
        ncF.species = filtereddata.species
        ncF.sitetype = filtereddata.sitetype
        
        # Create the variables
        #ncsitenames = ncF.createVariable('sitename','s', ('sitenames',))
        ncsitenames = ncF.createVariable('sitename',type(filtereddata.sitenames[0]), ('sitenames',))        
        
        nctime = ncF.createVariable('time', 'i', ('time',))

        nclon = ncF.createVariable('lon', 'f', ('sitenames','boxwidth'))
        nclat = ncF.createVariable('lat', 'f', ('sitenames','boxwidth'))
        ncpress = ncF.createVariable('pressure', 'f', ('sitenames','time','boxwidth','boxwidth','boxwidth',))
        ncconc = ncF.createVariable('conc', 'f', ('sitenames','time','boxwidth','boxwidth','boxwidth',))
        ncemiss = ncF.createVariable('emiss', 'f', ('sitenames','time','boxwidth','boxwidth',))
        
        # Fill the variables
        ncsitenames[:] = filtereddata.sitenames        
                
        # times as seconds since 1/1/2009
        nctime[:] = [(t - dt.datetime(2009,1,1,0,0,0)).total_seconds() for t in filtereddata.time]     
    
   
        nclon[:] = filtereddata.model_lon
        nclat[:] = filtereddata.model_lat
        ncpress[:] = filtereddata.model_pressure
        ncconc[:] = filtereddata.conc
        ncemiss[:] = filtereddata.emis
        
        
        # Give the variables some attributes        
        nclon.units = 'Degrees east'
        nclat.units = 'Degrees north'
        ncconc.units = filtereddata.concunits
        ncemiss.units = filtereddata.emissunits
        ncpress.units = filtereddata.pressureunits
        nctime.units = 'seconds since 2009-01-01 00:00:00'
        
        ncF.close()
        print "Written " + outdir + '/'+ filename


# Class to write out the data from a stationary platform
class write_ncdf_column:
    def __init__(self, filtereddata, outdir = 0, filename=0):
        
        import os
        
        if type(outdir) == int:
            outdir = os.path.dirname(filtereddata.mzfile)
    
        if type(filename) == int:
            if type(filtereddata.mzfile) == tuple:            
                filetimestamp = filtereddata.mzfile[0][filtereddata.mzfile[0].find('h1')+ 2: filtereddata.mzfile[0].find('.nc')]
            if type(filtereddata.mzfile) == str:
                filetimestamp = filtereddata.mzfile[filtereddata.mzfile.find('h1')+ 2: filtereddata.mzfile.find('.nc')]
            filename = filtereddata.species + '_' + filtereddata.case + '.mzt.h1'+  filetimestamp + '_' + filtereddata.sitetype+'.nc'

        print 'writing file: ' + outdir + '/'+ filename
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(outdir + '/'+ filename, 'w')
        
        # Create the dimensions
        ncF.createDimension('time', len(filtereddata.model_time))
        ncF.createDimension('lev', np.shape(filtereddata.model_pressure)[1])
        
        
        # Make some global attributes
        ncF.mzfile = filtereddata.mzfile
        ncF.sitefile = filtereddata.sitefile
        ncF.case = filtereddata.case
        ncF.species = filtereddata.species
        ncF.sitetype = filtereddata.sitetype
        
        # Create the variables 
        nctime = ncF.createVariable('time', 'i', ('time',))   
        
        nclon = ncF.createVariable('lon', 'f', ('time',))
        nclat = ncF.createVariable('lat', 'f', ('time',))
        
        ncpress = ncF.createVariable('pressure', 'f', ('time', 'lev',))
        ncconc = ncF.createVariable('conc', 'f', ('time', 'lev',))
        
        # Fill the variables               
        # times as seconds since 1/1/2009
        nctime[:] = [(t - dt.datetime(2009,1,1,0,0,0)).total_seconds() for t in filtereddata.model_time]     
   
        
   
        nclon[:] = filtereddata.model_lon
        nclat[:] = filtereddata.model_lat
        ncpress[:,:] = filtereddata.model_pressure
        ncconc[:,:] = filtereddata.conc
        
        
        # Give the variables some attributes        
        nclon.units = 'Degrees east'
        nclat.units = 'Degrees north'
        ncconc.units = filtereddata.concunits
        nctime.units = 'seconds since 2009-01-01 00:00:00'
        ncpress.units = filtereddata.pressureunits

        
        ncF.close()
        print "Written " + outdir + '/'+ filename


# Class to write out the data for the mobile platforms in the GONZI format
# this format has a lot of -9999's in the obs 
# need to filter these out first
# Extracts from the h0 files and fills in the Model_species variable in the GONZI file

class write_GONZI:
    def __init__(self, species = 'N2O', \
        modeloutput_dir = '/data/as13988/MOZART/', \
        obs_dir = '/shared_data/snowy/shared/GAUGE/', \
        outdir = 0, filename=0):
        
        import os
        import shutil
        import fnmatch
        
        
        gonzifile_acf = '/data/as13988/GAUGE/model_protocol_acf_v3.nc'      
        gonzifile_ferry = '/data/as13988/GAUGE/model_protocol_ferry_v3.nc'
         
        # Copy the original gonzi files
        # check if a copy already exists
        if os.path.exists(gonzifile_acf[:-3] + '_MZT.nc') != True:
            shutil.copy2(gonzifile_acf,  gonzifile_acf[:-3] + '_MZT.nc')

        if os.path.exists(gonzifile_ferry[:-3] + '_MZT.nc') != True:
            shutil.copy2(gonzifile_ferry,  gonzifile_ferry[:-3] + '_MZT.nc')         
         
         
        # Read in the gonzi aircraft file 
        gonzi_acf = read_sitefile_GONZI_nc(gonzifile_acf[:-3] + '_MZT.nc')        
        gonzi_ferry = read_sitefile_GONZI_nc(gonzifile_ferry[:-3] + '_MZT.nc')  
        
        # Find all the model h0 output files
        filepattern = 'FWDModelComparison_NewEDGAR*h0*'
             
        matches = []
        for root, dirnames, filenames in os.walk(modeloutput_dir+species+'/output/FWDModelComparison_NewEDGAR/'):
            for filename in fnmatch.filter(filenames, filepattern):
                matches.append(os.path.join(root, filename))   
                
        filenames = matches
        filenames.sort()
        
        
        # determine the years and months of these files and calculate a corresponding date
        years = [int(i[-19:-15]) for i in filenames]
        months = [int(i[-14:-12]) for i in filenames]
                
        filedates = [dt.datetime(years[j], months[j], 1) for j in np.arange(len(years))]
        
        
        # Only want to read in the model output for months/years in the obs
        # use bisect to find where the min and max dates intersect with the filedates
        start = bisect.bisect_left(filedates, min(gonzi_acf.time)) -1
        finish = bisect.bisect_right(filedates, max(gonzi_acf.time))
        
        modelfiles = filenames[start:finish]
        
        
        # create lists to store the model data
        model_time = []
        model_conc = []
        model_lat = []
        model_lon = []      
        model_pressure = []
        
        
        # Read in the model output files
        for i in modelfiles:
            modeltemp = read(i)
         
            model_time.extend(modeltemp.time)
            model_conc.extend(modeltemp.conc)
            model_pressure.extend(modeltemp.pressure)
        
        model_lat = modeltemp.lat
        model_lon = modeltemp.lon
        model_P0 = modeltemp.P0
        
         
        conc_index = np.where(gonzi_acf.species_tags == species)[0]      
         
        # Use the data matching code to extract the matching data
        matched_data = data_match_mobile(species = species, \
                                        model_conc = np.array(model_conc), \
                                        model_lat = np.array(model_lat), \
                                        model_lon = np.array(model_lon), \
                                        model_time = model_time, \
                                        model_P0 = model_P0, \
                                        model_pressure = np.array(model_pressure), \
                                        obs_conc = np.squeeze(np.array(gonzi_acf.conc[conc_index,:,0])), \
                                        obs_lat = np.squeeze(np.array(gonzi_acf.lat[:,0])), \
                                        obs_lon = np.squeeze(np.array(gonzi_acf.lon[:,0])), \
                                        obs_time = gonzi_acf.time, \
                                        obs_pressure = 0, \
                                        obs_alt = np.squeeze(np.array(gonzi_acf.alt[:,0])), \
                                        obs_alt_units = gonzi_acf.alt_units, \
                                        quiet = 1)   
                

        
       
        
        print 'Writing to file: ' + gonzifile_acf[:-3] + '_MZT.nc'
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(gonzifile_acf[:-3] + '_MZT.nc', 'a')
        
        # Extract the species variable
        species_var = ncF.variables['Model_'+species][:]
        
        # check that the length of the variable matches the length of the extracted data
        if len(species_var[gonzi_acf.good_index]) != len(matched_data.matched_conc):
            print 'Something has gone horribly wrong with the matching'
            print 'Data not inserted into the gonzi file'
        
        else:
            
            # Put the matched data into the non NAN sections of the variables 
            species_var[[gonzi_acf.good_index],0] = matched_data.matched_conc
            ncF.variables['Model_' + species][:] = species_var
        
            ncF.close()
            print'"Written file: ' + gonzifile_acf[:-3] + '_MZT.nc'


        # Match to the ferry file

        # Only want to read in the model output for months/years in the obs
        # use bisect to find where the min and max dates intersect with the filedates
        start = bisect.bisect_left(filedates, min(gonzi_ferry.time)) -1
        finish = bisect.bisect_right(filedates, max(gonzi_ferry.time))
        
        modelfiles = filenames[start:finish]
        
        
        # create lists to store the model data
        model_time = []
        model_conc = []
        model_lat = []
        model_lon = []      
        model_pressure = []
        
        
        # Read in the model output files
        for i in modelfiles:
            modeltemp = read(i)
         
            model_time.extend(modeltemp.time)
            model_conc.extend(modeltemp.conc)
            model_pressure.extend(modeltemp.pressure)
        
        model_lat = modeltemp.lat
        model_lon = modeltemp.lon
        model_P0 = modeltemp.P0
        
         
        conc_index = np.where(gonzi_ferry.species_tags == species)[0]      
         
        # Use the data matching code to extract the matching data
        matched_data = data_match_mobile(species = species, \
                                        model_conc = np.array(model_conc), \
                                        model_lat = np.array(model_lat), \
                                        model_lon = np.array(model_lon), \
                                        model_time = model_time, \
                                        model_P0 = model_P0, \
                                        model_pressure = np.array(model_pressure), \
                                        obs_conc = np.squeeze(np.array(gonzi_ferry.conc[conc_index,:,0])), \
                                        obs_lat = np.squeeze(np.array(gonzi_ferry.lat[:,0])), \
                                        obs_lon = np.squeeze(np.array(gonzi_ferry.lon[:,0])), \
                                        obs_time = gonzi_ferry.time, \
                                        obs_pressure = 0, \
                                        obs_alt = 0, \
                                        obs_alt_units = gonzi_ferry.alt_units, \
                                        quiet = 1)   
                

        
       
        
        print 'Writing to file: ' + gonzifile_ferry[:-3] + '_MZT.nc'
        
        #Write NetCDF file
        ncF = netCDF4.Dataset(gonzifile_ferry[:-3] + '_MZT.nc', 'a')
        
        # Extract the species variable
        species_var = ncF.variables['Model_'+species][:]
        
        # check that the length of the variable matches the length of the extracted data
        if len(species_var[gonzi_ferry.good_index]) != len(matched_data.matched_conc):
            print 'Something has gone horribly wrong with the matching'
            print 'Data not inserted into the gonzi file'
        
        else:
            
            # Put the matched data into the non NAN sections of the variables 
            species_var[[gonzi_ferry.good_index],0] = matched_data.matched_conc
            ncF.variables['Model_' + species][:] = species_var
        
            ncF.close()
            print'"Written file: ' + gonzifile_ferry[:-3] + '_MZT.nc'
            
            
# ___________________________________________________________________________________________________________________
# CODE TO READ THE DIFFERENT DATA TYPES
# ___________________________________________________________________________________________________________________
# Class to read in the filtered output
class read_ncdf_fixed:
    def __init__(self, filepattern = '*_TT.nc', filenames = 0, species='CH4', directory = 1):
        
        if type(filenames) == type(0):
            import fnmatch
            import os
            
            if type(directory) == type(1):
                directory = '/data/as13988/MOZART/'+species+'/output/FWDModelComparison_NewEDGAR/'
            
            filepattern = '*'+species+filepattern
            
            matches = []
            for root, dirnames, filenames in os.walk(directory):
                for filename in fnmatch.filter(filenames, filepattern):
                    matches.append(os.path.join(root, filename))   
                    
        filenames = matches
        # need to sort the files
        filenames.sort()
        
        if len(filenames) == 0 :
            print 'There are no files matching the file pattern in the given directory'

        else:
            
            for j in np.arange(len(filenames)):
                
                print 'Reading file : ' + str(filenames[j])
                
                data=netCDF4.Dataset(filenames[j])
                
                lon = data.variables['lon'][:]
                lat = data.variables['lat'][:]
                
                #sitenames = data.variables['sitename'][:] 
                sitenames = data.variables['sitename'][:] 
                
                time_j = np.transpose(data.variables['time'][:])
                conc_j = np.transpose(data.variables['conc'][:])
                emis_j = np.transpose(data.variables['emiss'][:])
                pressure_j = np.transpose(data.variables['pressure'][:])
                
                # Create the time variable
                dateunits = data.variables['time'].getncattr('units')
                sincedate = dateunits[dateunits.find('seconds since ') +14:-1] # seconds since 2014-01-01 00:00:00
                    
                dt_date_j = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in time_j]
                
                
                if j == 0:
                
                    conc = conc_j              
                    emis = emis_j
                    pressure =  pressure_j
                    dt_date = dt_date_j              
                    
                else:
                    
                    
                    conc = np.concatenate((conc, conc_j), axis=3)
                    emis = np.concatenate((emis, emis_j), axis=2)
                    pressure = np.concatenate((pressure, pressure_j), axis=3)                
                    dt_date = np.concatenate((dt_date, dt_date_j))
          
                    
            self.time = dt_date
            self.emis = emis
            self.conc = conc
            self.pressure = pressure
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
 

# Class to read in the filtered output
class read_ncdf_mobile:
    def __init__(self, species = 'CH4', filepattern = '*_mobile.nc', filenames = 0, \
        directory = 1):
        
        if type(filenames) == type(0):
            import fnmatch
            import os
            
            if type(directory) == type(1):
                directory = '/data/as13988/MOZART/'+species+'/output/FWDModelComparison_NewEDGAR/'
            
            filepattern = '*'+species+filepattern
             
            
            matches = []
            for root, dirnames, filenames in os.walk(directory):
                for filename in fnmatch.filter(filenames, filepattern):
                    matches.append(os.path.join(root, filename))   
                    
            filenames = matches
            
        # need to sort the files
        filenames.sort()
        
        if len(filenames) == 0:        
            print "There are no files in the given directory matching the file pattern"
        
        else: 
            for j in np.arange(len(filenames)):
                
                print 'Reading file : ' + str(filenames[j])
                
                data=netCDF4.Dataset(filenames[j])
                            
                model_lon_j = data.variables['model_lon'][:]
                model_lat_j = data.variables['model_lat'][:]
                obs_lon_j = data.variables['obs_lon'][:]
                obs_lat_j = data.variables['obs_lat'][:]
               
                #sitenames = data.variables['sitename'][:] 
                sitenames_j = data.variables['sitename'][:] 
                
                model_time_j = np.transpose(data.variables['model_time'][:])
                obs_time_j = np.transpose(data.variables['obs_time'][:])

                model_conc_j = np.transpose(data.variables['model_conc'][:])
                obs_conc_j = np.transpose(data.variables['obs_conc'][:])

                model_pressure_j = np.transpose(data.variables['model_pressure'][:])
                obs_pressure_j = np.transpose(data.variables['obs_pressure'][:])
                
                # Create the time variable
                dateunits = data.variables['model_time'].getncattr('units')
                sincedate = dateunits[dateunits.find('seconds since ') +14:-1] # seconds since 2014-01-01 00:00:00
                    
                model_dt_date_j = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in model_time_j]
                obs_dt_date_j = [ dt.datetime.strptime(sincedate, "%Y-%m-%d %H:%M:%S") + dt.timedelta(seconds=(i).astype('int')) for i in obs_time_j]
               
                
                if j == 0:
    
                    sitenames = sitenames_j            
                
                    model_conc = model_conc_j      
                    obs_conc = obs_conc_j                    

                    model_pressure =  model_pressure_j
                    obs_pressure =  obs_pressure_j
                    
                    model_dt_date = model_dt_date_j         
                    obs_dt_date = obs_dt_date_j         
                    
                    model_lat = model_lat_j
                    model_lon = model_lon_j
                    obs_lat = obs_lat_j
                    obs_lon = obs_lon_j
                    
    
                else:
                    
                    sitenames = np.concatenate((sitenames, sitenames_j))               
                    
                    model_conc = np.concatenate((model_conc, model_conc_j))
                    obs_conc = np.concatenate((obs_conc, obs_conc_j))

                    model_pressure = np.concatenate((model_pressure, model_pressure_j))    
                    obs_pressure = np.concatenate((obs_pressure, obs_pressure_j)) 
                    
                    model_dt_date = np.concatenate((model_dt_date, model_dt_date_j))
                    obs_dt_date = np.concatenate((obs_dt_date, obs_dt_date_j))
                    
                    model_lat = np.concatenate((model_lat, model_lat_j))
                    model_lon = np.concatenate((model_lon, model_lon_j))
                    obs_lat = np.concatenate((obs_lat, obs_lat_j))
                    obs_lon = np.concatenate((obs_lon, obs_lon_j))
    
                    
                    
            self.model_time = model_dt_date
            self.obs_time = obs_dt_date
            self.model_conc = model_conc
            self.obs_conc = obs_conc
            self.model_pressure = model_pressure
            self.obs_pressure = obs_pressure
            self.model_lon = model_lon
            self.model_lat = model_lat
            self.obs_lon = obs_lon
            self.obs_lat = obs_lat
            self.sitenames = sitenames
            self.filenames = filenames         
            self.species = str(data.__getattribute__('species')).strip()
            self.case = str(data.__getattribute__('case')).strip()
            self.sitetype = str(data.__getattribute__('sitetype')).strip()
            self.sitefile = str(data.__getattribute__('sitefile')).strip()
            
            self.concunits = data.variables['model_conc'].getncattr('units')
            self.pressureunits = data.variables['model_pressure'].getncattr('units')
       
# ___________________________________________________________________________________________________________________
# CODE TO PLOT THE DIFFERENT DATA TYPES
# ___________________________________________________________________________________________________________________

# Class to plot the filtered output
# Plots the output of read_ncdf_fixed
class plot_ncdf_fixed:
    def __init__(self, data, sitename = 'mhd', scaling = 1e06, x_range = 1, save_plot = 0):
        
        import matplotlib.ticker as ticker
        import matplotlib.pyplot as plt
        
        # Extract the model data for the given site
        sites = data.sitenames
        matched = np.where(sites == sitename)[0]
        
          
        if len(matched) > 1:
            # This means that the site was listed twice in the site file as two measurement networks share it.
            # As we're pulling out model data at the same location it doesn't matter which version we use
            # Default to using the first one
            matched = [matched[0]]
            
        
        if len(matched) == 0:
            print 'There is no sitename matching : ' + sitename
            print 'Check the sitename or try using lowercase'
        else: 
            time = data.time
            
            # conc is lon x lat x lev x time x site
            conc = np.squeeze(data.conc[:,:,:,:,matched])
            #pressure = np.squeeze(data.pressure[:,:,:,:,matched])
                       
            
            reshapedconc = np.reshape(conc,(27,np.shape(conc)[-1]))
            
            # Indicies for lats, lons and levs
            lat_i = [1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3]
            lon_i = [1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3]
            lev_i = [1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]  
            
            n_colours = 27        
            
            colours = generatecolours(n_colours).RGB
            
            # Plot model ooutput for fixed site data
            fig = plt.figure()
            fig.subplots_adjust(right = 0.8)
            
            legend_spacing = np.arange(n_colours)   
            
                        
            
            for i in np.arange(n_colours):      
                
                conc_i = reshapedconc[i]        
    
                # Plot the data
                plt1 = plt.subplot()
                            
                plt1.plot(time, conc_i*scaling, "-", color = colours[i], markersize = 3)
                    
                y_formatter = ticker.ScalarFormatter(useOffset=False)
                plt1.yaxis.set_major_formatter(y_formatter)
                
                x_tickno_formatter = ticker.MaxNLocator(5)
                plt1.xaxis.set_major_locator(x_tickno_formatter)
                
                plt1.set_title('Model output at '+ sitename)
                plt1.set_ylabel(data.species + '(' + data.concunits + '*' + str(scaling) + ')')
                plt1.set_xlabel('Time')
    
                legend_i = str(data.lat[matched,lat_i[i]-1][0]) + ', ' +str(data.lon[matched,lon_i[i]-1][0]) + ', ' + str(lev_i[i]-1)
    
                plt.figtext(0.82, 0.85-(0.03*legend_spacing[i]), legend_i, verticalalignment='bottom', \
                horizontalalignment='left', color=colours[i], fontsize=8)
                    
            plt1.plot(time, reshapedconc[13]*scaling, "--", color = 'black', markersize = 3)
    
    
            plt.figtext(0.82, 0.88, 'Lat, Lon, Lev', verticalalignment='bottom', horizontalalignment='left', color='black', fontsize=8)
            plt.figtext(0.6, 0.2, 'Central point ---', verticalalignment='bottom', horizontalalignment='left', color='black', fontsize=8)
               
            
            if type(x_range) != type(1):
                plt1.set_xlim(x_range)            
            
            if save_plot != 0:
                outdir = os.path.dirname(os.path.dirname(data.filenames[0]))
                fig.savefig(outdir + '/plots/' + data.species+ '_'+ sitename+ '_Model.png', dpi=100)
                print'Figure saved as : ' + outdir + '/plots/' + data.species+ '_'+ sitename+ '_Model.png'
            
            
            plt.show()
        
            plt.close()
            
            
            #pdb.set_trace()
            
            # Plot model output and obs    
            # Read in the obs
            obs = read_fixed_sitefile_nc(species = data.species, dir = '/data/shared/GAUGE/')
        
            # Extract the data for the given site
            index = np.where(obs.site == sitename)[0] 
            
            # loop through each network at each site
            for j in index:
                
                print 'Plotting site ' + obs.site[j]
                print 'Plotting data from ' + obs.network[j]                
                
                # find where the obs exist i.e. != nan   
                obs_index = np.where(np.isfinite(np.squeeze(obs.conc[j,:])))[0]      
                
                print 'Number of finite elements in obs data = ' + str(len(obs_index))
                
                if len(obs_index) != 0:   
                
                    obs_conc = np.squeeze(obs.conc[j,:])[obs_index]
                    obs_time = [obs.time[i] for i in obs_index]
                    
                    # Only plot model output for the same time range as the obs
                    start = bisect.bisect_left(time, min(obs_time)) -1
                    finish = bisect.bisect_right(time, max(obs_time))
                
                    model_time = time[start:finish]
                    model_conc = reshapedconc[13][start:finish]
    
                else:
                    model_time = time
                    model_conc = reshapedconc[13]
                    obs_conc = []
                    obs_time = []
                    
                    
                #pdb.set_trace()  
                
                fig = plt.figure()
                 
                # Plot the data
                plt1 = plt.subplot()
                            
                plt1.plot(model_time, model_conc*scaling, "bo-", markersize = 3, markeredgecolor = 'blue')
                plt1.plot(obs_time, np.squeeze(obs_conc), "ro-", markersize = 3, markeredgecolor = 'red')
                 
                y_formatter = ticker.ScalarFormatter(useOffset=False)
                plt1.yaxis.set_major_formatter(y_formatter)
                
                x_tickno_formatter = ticker.MaxNLocator(5)
                plt1.xaxis.set_major_locator(x_tickno_formatter)
                
                plt1.set_title('Model output and obsservations at '+ sitename + '(' + obs.network[j] + ')' )
                plt1.set_ylabel(data.species + '(' + data.concunits + '*' + str(scaling) + ')')
                plt1.set_xlabel('Time')
    
                
                if save_plot != 0:
                    outdir = os.path.dirname(os.path.dirname(data.filenames[0]))
                    fig.savefig(outdir + '/plots/' + data.species+ '_'+ obs.network[j] + '_'+ sitename+ '.png', dpi=100)
                    print'Figure saved as : ' + outdir + '/plots/' + data.species+ '_'+ obs.network[j] + '_'+ sitename+ '.png'
                
                
                plt.show()
            
                plt.close()
         
            
    
           

        
        
        
# Class to plot the filtered output
# Plots the output of read_ncdf_mobile
# defaults to not plotting the obs concentrations as we only have these for CH4 at the moment
class plot_ncdf_mobile:
    def __init__(self, data, sitename = 'ferry', save_plot = 0, scaling = 1e06):
        
        import matplotlib.ticker as ticker
        import matplotlib.pyplot as plt
        
        # Extract the data for the given site
        sites = data.sitenames
        matched = np.where(sites == sitename)[0]
        
        if len(matched) == 0 :        
            print 'There are no sitenames matching ' + str(sitename)
            
        else:
            # conc is of dimension time 
            # sitenames are also of dimension time        
        
            time = np.squeeze(data.obs_time[matched])
            
            model_conc = np.squeeze(data.model_conc[matched])
            obs_conc = np.squeeze(data.obs_conc[matched])
                
            
            model_pressure = np.squeeze(data.model_pressure[matched])
            obs_pressure = np.squeeze(data.obs_pressure[matched])
            
            model_lat = np.squeeze(data.model_lat[matched])
            obs_lat = np.squeeze(data.obs_lat[matched])
            
            model_lon = np.squeeze(data.model_lon[matched])
            obs_lon = np.squeeze(data.obs_lon[matched])
            
            model_GT180 = np.where(model_lon > 180)[0]
            obs_GT180 = np.where(obs_lon > 180)[0]
            
                   
            
            if len(model_GT180) != 0:
                model_lon[model_GT180] = model_lon[model_GT180] - 360
            
            if len(obs_GT180) != 0:
                obs_lon[obs_GT180] = obs_lon[obs_GT180] - 360


            # Plot of 
            fig = plt.figure()
            
            # Plot the data
            plt1 = plt.subplot(6,1,1)
                        
            plt1.plot(time, model_conc*scaling, "bo-",  markersize = 3)
            
            if no_obs !=0:
                plt1.plot(time, obs_conc, "ro-",  markersize = 3)
                
            # plot the model and the obs time/space matched output
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
             
            plt1.set_title(data.species + 'model output (blue) and obs (red) for '+ sitename)
            plt1.set_ylabel(data.species + '(' + data.concunits + '*' + str(scaling) + ')')

            
            # plot the difference between the model and the obs
            plt5 = plt.subplot(6,1,2)
                        
            plt5.plot(time, model_conc*scaling - obs_conc, "go-", markersize = 3)
            
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt5.yaxis.set_major_formatter(y_formatter)
             
            plt5.set_ylabel('Conc diff')
            
            
            # plot the model and obs pressures
            plt2 = plt.subplot(6,1,3)
                        
            plt2.plot(time, model_pressure, "bo-", markersize = 3)
            plt2.plot(time, obs_pressure, "ro-", markersize = 3)
                
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt2.yaxis.set_major_formatter(y_formatter)
             
            plt2.set_ylabel('Pressure (Pa)')


            # plot the model and obs pressure difference
            plt6 = plt.subplot(6,1,4)
                        
            plt6.plot(time, model_pressure- obs_pressure, "go-", markersize = 3)
                
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt6.yaxis.set_major_formatter(y_formatter)
             
            plt6.set_ylabel('Pressure diff (Pa)')
            
            
            # plot the model and obs lats
            plt3 = plt.subplot(6,1,5)
                        
            plt3.plot(time, model_lat, "bo-", markersize = 3)
            plt3.plot(time, obs_lat, "ro-", markersize = 3)
                
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt3.yaxis.set_major_formatter(y_formatter)
             
            plt3.set_ylabel('lat')


            # plot the model and obs lons
            plt4 = plt.subplot(6,1,6)
                        
            plt4.plot(time, model_lon, "bo-", markersize = 3)
            plt4.plot(time, obs_lon, "ro-", markersize = 3)
                
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt4.yaxis.set_major_formatter(y_formatter)
             
            plt4.set_ylabel('lon')

            fig.set_size_inches(6,8)

            if save_plot != 0:
                outdir = os.path.dirname(os.path.dirname(data.filenames[0]))
                fig.savefig(outdir + '/plots/' + data.species+ '_'+ sitename+ '.png', dpi=100)
                print'Figure saved as : ' + outdir + '/plots/' + data.species+ '_'+ sitename+ '.png'
            
            
            plt.show()      



# Generates a given number of HSV or RGB tuples spread over colour space 
class generatecolours:
    def __init__(self, N):
        
        import colorsys        
        
        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        
        self.RGB = RGB_tuples
        self.HSV = HSV_tuples


