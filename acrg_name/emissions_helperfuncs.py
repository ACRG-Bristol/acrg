#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 09:50:09 2018

@author: lw13938

This script contains functions to read emissions data sets, convert them to 
mol/m2/s and regrid them onto a desired lat lon grid.
This is primarily useful for creating priors for inversions.

The functions for the GFED and EDGAR data are general. 
Currently the others are just CH4. 

These functions regrid on each time step in a loop which is very slow. 
This should probably be changed at some point to regrid over all timesteps in 
one go.

"""

import numpy as np
import glob 
import h5py
from acrg_grid.regrid import regrid2d
import xarray as xr
from acrg_grid import areagrid
import datetime
import pandas as pd
import os

def getGFED(year, lon_out, lat_out, timeframe='monthly', monthrange = [1,2,3,4,5,6,7,8,9,10,11,12], soi='CH4', incagr=False):
    """
    Gets GFED 4.1s data, puts into mol/m2/s and regrids to desired 
    lats and lons for a given year and months.
    Will output monthly, daily or 3 hourly emissions.
    
    Args:
        year (int): 
            Year of interest
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        timeframe (str):
            Time over which data is aggregated.
            Can be either 'monthly', 'daily' or '3hourly'.
            Only 'monthly' data available before 2003.
            Default = 'monthly'
        monthrange (list):
            The months that you want data for. NB these don't have to be 
            sequential.
            Default = [1,2,3,4,5,6,7,8,9,10,11,12] (all months)
        soi (str):
            Which species you want to look at as defined in 
            GFED4_Emission_Factors.csv.
            e.g. soi = 'CO2'
            Default = 'CH4'
        incagr (bool):
            To include agricultural waste burning (which is otherwise in EDGAR) 
            set incagr=True.
            Default = False
        
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    
    Example: 
        emissions  = getGFED(2012, lat, lon)
    
    Note:
        The regridding in this function is really slow as it does it in a loop
        over all time steps.
        There's definitely a faster way of doing this.
    """    
    months       = '01','02','03','04','05','06','07','08','09','10','11','12'
    if incagr == False:
       sources      = 'SAVA','BORF','TEMF','DEFO','PEAT'
       sourceindex = [7,5,3,1,11] #
    else:
       sources      = 'SAVA','BORF','TEMF','DEFO','PEAT','AGRI'
       sourceindex = [7,5,3,1,11,9] #
        
    
    directory    = '/data/shared/Gridded_fluxes/GFED4/fire_emissions_v4_R1_1293/data'
    soi = soi.upper()
    
    ###
    ###  Read in emission factors
    ###
    species = [] # names of the different gas and aerosol species
    EFs     = np.zeros((41, 23)) # Species and sources
    
    k = 0
    f = open(directory+'/GFED4_Emission_Factors.csv')
    while 1:
        line = f.readline()
        line = line.rstrip()
    
        if line.startswith('carbon') == False:
            continue
        else:
            contents = line.split(",")
            species.append(contents[0])
            xvals = [i for i,x in enumerate(contents) if x=='x']
            for j in xvals:
               contents[j] = np.nan
            EFs[k,:] = contents[1:]
            k += 1
            while 1:
                line = f.readline()
                line = line.rstrip()
                if line.startswith(','):
                    break
                contents = line.split(",")
                species.append(contents[0])
                xvals = [i for i,x in enumerate(contents) if (x=='x') or (x=='')]
                for j in xvals:
                    contents[j] = np.nan
                EFs[k,:] = np.asarray(contents[1:])
                k += 1
        break
          
    f.close()
    
    # Find species of interest in file
    kc = 0
    for k in species:
        if k.find(soi) > -1:
            soix = kc
        kc += 1
    EF = EFs[soix,:]
    
    #Check that there are daily/3 hourly data available
    if (year < 2003) & (timeframe != 'monthly'):
        print "No %s data before 2003. Switching to monthly emissions." % timeframe
        timeframe = 'monthly'
    
    #Check to see if desired year has a dataset
    possyears = np.empty(shape=[0,0],dtype=int)
    for f in glob.glob(directory+'/GFED4.1s_*'+'.hdf5'):
        fname = f.split('/')[-1]
        fyear = fname[9:13]      #Extract year from filename
        possyears = np.append(possyears, int(fyear))
    if year > max(possyears):
        print "%s is later than latest year in GFED v4.1 database" % str(year)
        print "Using %s as the closest year" % str(max((possyears)))
        year = max(possyears)
    if year < min(possyears):
        print "%s is earlier than earliest year in GFED v4.1 database" % str(year)
        print "Using %s as the closest year" % str(min((possyears)))
        year = min(possyears)
    
    #Get emissions
    string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    f = h5py.File(string, 'r')   
    lat = np.flipud(np.unique(np.array(f['lat']) - 0.125))
    lon = np.unique(np.array(f['lon']) + 0.125)
    leapyears = np.arange(1900, 2100, 4)    
    #dim = np.array([31,28,31,30,31,30,31,30,31,30,31,31])  
    dim = np.array([31,28,31,30,31,30,31,31,30,31,30,31]) #rt17603: Updated to align with days in a month 
    if (np.min(abs(year-leapyears)) == 0):
        dim[1] = 29
    
    monthrange = [m-1 for m in monthrange] # rt17603: Reset to zero-indexed list (jan=0,feb=1 etc)
    dim = dim[monthrange]    #Use only desired months
    
    if timeframe == 'monthly':
        emissions = np.zeros((len(dim), len(lat), len(lon)))
        for month in range(len(dim)):
           convert2secs = dim[month]*24*3600
            # read in DM emissions
           string = '/emissions/'+months[month]+'/DM'
           DM_emissions = f[string][:]
           for source in range(len(sources)):
               # read in the fractional contribution of each source
               string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
               contribution = f[string][:]
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), and the emission factor (g per kg DM burned)
               #emissions += DM_emissions * contribution * EF[sourceindex[source]]
               #Then convert from g/m2/month to mol/m2/s
               emissions[month, :,:] += (DM_emissions * contribution * EF[sourceindex[source]]) / (EF[0] * convert2secs)
    elif timeframe == 'daily':
        emissions = np.zeros((np.sum(dim), len(lat), len(lon)))
        convert2secs = 24*3600
        d = 0
        for month in range(len(dim)):
           convert2secs = dim[month]*24*3600
            # read in DM emissions
           string = '/emissions/'+months[month]+'/DM'
           DM_emissions = f[string][:]
           contribution = f[string][:]
           for day in range(dim[month]):
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), the emission factor (g per kg DM burned), and fractional
               #daily contribution.
               #Then convert from g/m2/day to mol/m2/s
               daystr = '/emissions/'+months[month]+'/daily_fraction/day_'+str(day+1)
               dayfrac = f[daystr][:]
               for source in range(len(sources)):
                   # read in the fractional contribution of each source
                   string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
                   emissions[d, :,:] += (DM_emissions * contribution * dayfrac * EF[sourceindex[source]]) / (EF[0] * convert2secs)
               d += 1
    elif timeframe == '3hourly':
        emissions = np.zeros((np.sum(dim)*8, len(lat), len(lon)))
        convert2secs = 3*3600
        h = 0
        #rt17603: Made into a list rather than a list containing one tuple.
        diurnalcyclenames = ["UTC_0-3h", "UTC_3-6h", "UTC_6-9h", "UTC_9-12h", "UTC_12-15h", "UTC_15-18h", "UTC_18-21h", "UTC_21-24h"]
        for month in range(len(dim)):
            # read in DM emissions
           string = '/emissions/'+months[month]+'/DM'
           DM_emissions = f[string][:]
           contribution = f[string][:]
           for day in range(dim[month]):
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), the emission factor (g per kg DM burned), fractional
               #daily contribution and 3 hourly fractional contribution.
               #Then convert from g/m2/3 hour to mol/m2/s
               daystr = '/emissions/'+months[month]+'/daily_fraction/day_'+str(day+1)
               dayfrac = f[daystr][:]
               for dc in range(len(diurnalcyclenames)):
                   dcstr = '/emissions/'+months[month]+'/diurnal_cycle/'+diurnalcyclenames[dc]
                   dcfrac = f[dcstr][:]
                   for source in range(len(sources)):
                       # read in the fractional contribution of each source
                       string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
                       emissions[h, :,:] += (DM_emissions * contribution * dayfrac * dcfrac * EF[sourceindex[source]]) / (EF[0] * convert2secs)
                   h += 1
        
    #Regrid the data to desired region
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = emissions.shape[0]
    
    narr = np.zeros((nlat, nlon, nt))
   
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(emissions[i,:,:], lat, lon,
                             lat_out, lon_out)
    return(narr)
    
def getedgarannualtotals(year, lon_out, lat_out, soi='CH4'):
    """
    Get annual emission totals for species of interest from EDGAR v4.3.2 data.
    At the time of making this the only species available are CH4 and N20.
    Regrids to the desired lats and lons.
    
    Args:
        year (int): 
            Year of interest
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        soi (str):
            Which species you want to look at. 
            e.g. soi = 'CH4'
            Default = 'CH4'
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon]
    
    """
    
    soi = soi.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = '/data/shared/Gridded_fluxes/'+soi+'/EDGAR_v4.3.2/v432_'+soi+'_TOTALS_nc/'
    
    #Check to see range of years. If desired year falls outside of this range 
    #then take closest year
    possyears = np.empty(shape=[0,0],dtype=int)
    for f in glob.glob(edpath+'v432_'+soi+'_*'):
        fname = f.split('/')[-1]
        fyear = fname[9:13]      #Extract year from filename
        possyears = np.append(possyears, int(fyear))
    if year > max(possyears):
        print "%s is later than max year in EDGAR database" % str(year)
        print "Using %s as the closest year" % str(max((possyears)))
        year = max(possyears)
    if year < min(possyears):
        print "%s is earlier than min year in EDGAR database" % str(year)
        print "Using %s as the closest year" % str(min((possyears)))
        year = min(possyears)
        
    
    #Read in EDGAR data of annual mean CH4 emissions
    #units are in kg/m2/s
    edgar = edpath+'v432_'+soi+'_'+str(year)+'.0.1x0.1.nc'
    
    ds = xr.open_dataset(edgar)
    soiname = 'emi_'+soi.lower()    
    tot = ds[soiname].values*1e3/16.04
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon))    
    narr, reg = regrid2d(tot, lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)    
    
    
def getothernaturalCH4(lon_out, lat_out):
    """
    Gets monthly climatology for 'other natural' CH4 contributions: 
    (volcanoes, termites, hydrates) and regrids these to desired lat and lons.
    From Fung et al 1987
    File orginally created by Anne (I think)
    These are in kg/m2/s
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
        
        
    """
    path = '/data/shared/GAUGE/CH4/'
    otherfn = 'CH4_flux_natural_global_climatology.nc'
    ods = xr.open_dataset(path+otherfn)
    lat = ods.latitude.values
    lon = ods.longitude.values
    
    # Total hydrates should be scaled down to ~5 Tg/year 
    # (https://doi.org/10.1002/2016RG000534)
    area = areagrid(lat, lon)
    ods.HYDR.values *= 5./np.sum(ods.HYDR.values*area*1e-9*3600*24*365)
    
    oemiss =  (ods.HYDR.values + ods.OCEAN.values + ods.VOLC.values + ods.TERMITE.values)*1e3/16.04
    
    #Regrid the data to desired region
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = len(ods.date)
    
    narr = np.zeros((nlat, nlon,nt))
    
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(oemiss[i,:,:], lat, lon,
                             lat_out, lon_out)
    return(narr)
    
def getsoilsinkCH4(lon_out, lat_out):
    """
    This gets the monthly climatology of CH4 soil sinks
    From Bousquet et al 2006
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
            N.B. These are negative as it's a sink.
    """
    path = '/data/shared/GAUGE/CH4/'
    ssinkfn = 'CH4_flux_soilssink_global_climatology.nc'
    
    ds = xr.open_dataset(path+ssinkfn)
    
    soilsink = ds.SOILS.values*1e3/16.04
    lat_in = ds.latitude.values
    lon_in = ds.longitude.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = len(ds.date)
    
    narr = np.zeros((nlat, nlon, nt))
    
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(soilsink[i,:,:], lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)
    
def getbloomwetlandsCH4(year, lon_out, lat_out, monthly=True):
    """
    Global wetlands CH4 emissions from Bloom et al, regridded to desired 
    lats and lons for year of interest.
    
    Args:
        year (int):
            Year of interest.
        lon_out (array): 
            Longitudes to output the data on.
        lat_out (array):
            Latitudes to output the data on.
        monthly (bool):
            Should the data be averaged by month? 
            True (default) = monthly averaged
            False = daily averaged
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    """
    path = '/data/shared/GAUGE/CH4/'
    bloomwetlands = 'CH4_flux_wetlands_rice_global_2003_2009.nc'
    
    possyears = np.arange(7) + 2003
    if year > max(possyears):
        print "%s is later than max year in Bloom et al wetlands data" % str(year)
        print "Using %s as the closest year" % str(max((possyears)))
        year = max(possyears)
    if year < min(possyears):
        print "%s is earlier than min year in Bloom et al data" % str(year)
        print "Using %s as the closest year" % str(min((possyears)))
        year = min(possyears)
    
    ds = xr.open_dataset(path+bloomwetlands)
    
    d = ds.date.values.astype(str)
    ddt = np.empty_like(d)
    for i in range(len(ddt)):
       ddt[i] = datetime.datetime.strptime(d[i], '%Y%m%d')
    ds.date.values = pd.to_datetime(ddt)
    ds = ds.sel(date=str(year))
    
    if monthly == True:
        ds = ds.resample('M', 'date')
    
    bloomch4 = ds.CH4_FLUX.values*1e3/16.04
    lat_in = ds.latitude.values
    lon_in = ds.longitude.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = len(ds.date)
    
    narr = np.zeros((nlat, nlon, nt))
    
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(bloomch4[i,:,:], lat_in, lon_in,
                                 lat_out, lon_out)
    return(narr)
    

def getnaeiandedgarCH4(lon_out, lat_out):
    """
    Get combined NAEI and EDGAR emissions using NAEI 2013 and EDGAR 2010
    Regrids to desired lats and lons.
    Note sure how useful this is...
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    
    Note:
        Raw NAEI data are in /dagage2/agage/metoffice/naei/naei_raw_priors
        Currently this is for CH4 and N20: 2012, 2013, 2015.
        A more useful function would be to read this in etc. 
        
    """
    
    fdir='/dagage2/agage/metoffice/naei/'
    fn = 'prior_UK_NAEI2013detailed_ch4_EU_EDGAR2010_UK40uncertain.txt'
    
    with open(fdir+fn, 'r') as f:
    #Get extent of grid defined for NAEI / EDGAR
        line = f.readline()
        line = line.rstrip()
        llextent = np.asarray(line.split()[0:4]).astype(float)
        
        #Get lon/lat grid spacing
        line = f.readline()
        line = line.rstrip()
        llspacing = np.asarray(line.split()[0:2]).astype(float)
        
        #Get emissions for defined grid, units are in g/m2/s
        while 1:
            line = f.readline()
            line = line.rstrip()
            if line.startswith('   i') == False:
                    continue
            else:
                k=0
                while True:
                    line = f.readline()
                    if not line: break
                    line = line.rstrip()
                    dat = np.asarray(line.split()[0:3])
                    if k == 0:
                        dataarr = dat.astype(float)
                        k = 1
                    else:
                        dataarr = np.vstack([dataarr, dat.astype(float)])
            break
        
    lon = np.unique(dataarr[:,0]*llspacing[0] + llextent[0] + llspacing[0]/2.)
    lat = np.unique(dataarr[:,1]*llspacing[1] + llextent[2] + llspacing[1]/2.)
    emisarr = dataarr[:,2].reshape(len(lon), len(lat))/16.04

    #Regrid the data to desired region
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon))
    
    narr, reg = regrid2d(emisarr, lat, lon,
                             lat_out, lon_out)
    return(narr)
    
def getNAEI(year, lon_out, lat_out, soi, sector):
 
    """
    Converts raw NAEI into gridded emissions data in mol/m2/s
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        year (int): 
            Year of interest
        soi (str):
            Which species you want to look at. 
            Currently only 'CH4' or 'N2O'
        sector (str):
            Which sector to look at. Options are: 
            "energyprod","domcom","indcom","indproc","offshore","roadtrans",
            "othertrans","waste","agric","nature","points","total",
            "totalexcship"
        
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon]
    
    Example: 
        emissions  = getNAEI(2013, lon, lat, 'ch4', 'total')
    
    """   
    
    soi = soi.lower()   
    naeidir = '/dagage2/agage/metoffice/naei/naei_raw_priors/'
    fn = naeidir+'LatLong1km_'+soi.upper()+'_'+str(year)+'.csv.gz'
    #Checks
    if not os.path.isfile(fn):
        print 'Raw NAEI file LatLong1km_'+soi.upper()+'_'+str(year)+'.csv.gz does not exist'
        print 'Returning None'
        return None
    sectorlist = ["energyprod","domcom","indcom","indproc","offshore","roadtrans",
    "othertrans","waste","agric","nature","points","total","totalexcship"]
    if sector not in sectorlist:
        print 'Sector not one of:'
        print sectorlist
        print 'Returning None'
        return None
        
    #Read emissions as tonnes/km^2/yr and naei lat and lons
    df = pd.DataFrame.from_csv(fn)   
    lat = np.asarray(df.Latitude)
    lon = np.asarray(df.Longitude)    
    emissions = np.asarray(df[sector])
    
    #Make a square grid for the emissions
    latarr = np.arange(min(lat), max(lat), 0.01)
    lonarr = np.arange(min(lon), max(lon), 0.01)
    grdemis = np.zeros((len(latarr), len(lonarr)))

    for i in range(len(emissions)):
        ilat = np.where(abs(latarr-lat[i]) == np.min(abs(latarr - lat[i])) )
        ilon = np.where(abs(lonarr-lon[i]) == np.min(abs(lonarr - lon[i])) )
        grdemis[ilat, ilon] = emissions[i]

    #Convert to mol/m2/s        
    if soi == 'ch4':
        speciesmm = 16.0425
    if soi == 'n20':
        speciesmm = 44.013        
    if year % 4 == 0:
        diy = 365
    else:
        diy = 366    
    grdemis = grdemis*1e3 / (diy * 3600*24) / speciesmm
    
    #Regrid to desired lats and lons
    narr, reg = regrid2d(grdemis, latarr, lonarr,
                                 lat_out, lon_out)

    return(narr)
