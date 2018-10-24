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
from __future__ import print_function

import numpy as np
import glob 
import h5py
import os
from acrg_grid.regrid import regrid2d
from acrg_grid import areagrid
import datetime
import pandas as pd
import xarray as xr
import datetime as dt
from dateutil.relativedelta import relativedelta
from acrg_tdmcmc.tdmcmc_post_process import molar_mass
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

data_path = os.getenv("DATA_PATH")

def getGFED(year, lon_out, lat_out, timeframe='monthly', months = [1,2,3,4,5,6,7,8,9,10,11,12], soi='CH4', incagr=False):
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
        months (list):
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
    months_str     = '01','02','03','04','05','06','07','08','09','10','11','12'
    if incagr == False:
       sources     = 'SAVA','BORF','TEMF','DEFO','PEAT'
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
        print("No %s data before 2003. Switching to monthly emissions." % timeframe)
        timeframe = 'monthly'
    
    #Check to see if desired year has a dataset
    possyears = np.empty(shape=[0,0],dtype=int)
    for f in glob.glob(directory+'/GFED4.1s_*'+'.hdf5'):
        fname = f.split('/')[-1]
        fyear = fname[9:13]      #Extract year from filename
        possyears = np.append(possyears, int(fyear))
    if year > max(possyears):
        print("%s is later than latest year in GFED v4.1 database" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if year < min(possyears):
        print("%s is earlier than earliest year in GFED v4.1 database" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
    
    #Get emissions
    string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    f = h5py.File(string, 'r')   
    lat = np.flipud(np.unique(np.array(f['lat']) - 0.125))
    lon = np.unique(np.array(f['lon']) + 0.125)
    leapyears = np.arange(1900, 2100, 4)    
    
    dim = np.array([31,28,31,30,31,30,31,31,30,31,30,31]) #rt17603: Updated to align with days in a month 
    if (np.min(abs(year-leapyears)) == 0):
        dim[1] = 29
    
    months = [m-1 for m in months] # rt17603: Reset to zero-indexed list (jan=0,feb=1 etc)
    dim = dim[months]    #Use only desired months
    
    if timeframe == 'monthly':
        emissions = np.zeros((len(dim), len(lat), len(lon)))
        # rt17603: Updated to run over months rather than range(len(dim)), using i for dim dimension
        for i,month in enumerate(months):
           convert2secs = dim[i]*24*3600
            # read in DM emissions
           string = '/emissions/'+months_str[month]+'/DM'
           DM_emissions = f[string][:]
           for source in range(len(sources)):
               # read in the fractional contribution of each source
               string = '/emissions/'+months_str[month]+'/partitioning/DM_'+sources[source]
               contribution = f[string][:]
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), and the emission factor (g per kg DM burned)
               #emissions += DM_emissions * contribution * EF[sourceindex[source]]
               #Then convert from g/m2/month to mol/m2/s
               emissions[i, :,:] += (DM_emissions * contribution * EF[sourceindex[source]]) / (EF[0] * convert2secs)
    elif timeframe == 'daily':
        emissions = np.zeros((np.sum(dim), len(lat), len(lon)))
        convert2secs = 24*3600
        d = 0
        
        for i,month in enumerate(months):
           days = dim[i]
           convert2secs = days*24*3600
            # read in DM emissions
           string = '/emissions/'+months_str[month]+'/DM'
           DM_emissions = f[string][:]
           contribution = f[string][:]
           # rt17603: Updated to run over range(days) as defined by dim[i], rather than dim[month]
           for day in range(days):
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), the emission factor (g per kg DM burned), and fractional
               #daily contribution.
               #Then convert from g/m2/day to mol/m2/s
               daystr = '/emissions/'+months_str[month]+'/daily_fraction/day_'+str(day+1)
               dayfrac = f[daystr][:]
               for source in range(len(sources)):
                   # read in the fractional contribution of each source
                   string = '/emissions/'+months_str[month]+'/partitioning/DM_'+sources[source]
                   emissions[d, :,:] += (DM_emissions * contribution * dayfrac * EF[sourceindex[source]]) / (EF[0] * convert2secs)
               d += 1
    elif timeframe == '3hourly':
        emissions = np.zeros((np.sum(dim)*8, len(lat), len(lon)))
        convert2secs = 3*3600
        h = 0
        #rt17603: Made into a list rather than a list containing one tuple.
        diurnalcyclenames = ["UTC_0-3h", "UTC_3-6h", "UTC_6-9h", "UTC_9-12h", "UTC_12-15h", "UTC_15-18h", "UTC_18-21h", "UTC_21-24h"]
        
        for i,month in enumerate(months):
           days = dim[i]
           # read in DM emissions
           string = '/emissions/'+months_str[month]+'/DM'
           DM_emissions = f[string][:]
           contribution = f[string][:]
           
           for day in range(days):
               # calculate emissions as the product of DM emissions (kg DM per 
               # m2 per month), the fraction the specific source contributes to 
               # this (unitless), the emission factor (g per kg DM burned), fractional
               #daily contribution and 3 hourly fractional contribution.
               #Then convert from g/m2/3 hour to mol/m2/s
               daystr = '/emissions/'+months_str[month]+'/daily_fraction/day_'+str(day+1)
               dayfrac = f[daystr][:]
               for dc in range(len(diurnalcyclenames)):
                   dcstr = '/emissions/'+months_str[month]+'/diurnal_cycle/'+diurnalcyclenames[dc]
                   dcfrac = f[dcstr][:]
                   for source in range(len(sources)):
                       # read in the fractional contribution of each source
                       string = '/emissions/'+months_str[month]+'/partitioning/DM_'+sources[source]
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
    At the time of making this the only species available are CH4 and N2O.
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
        print("%s is later than max year in EDGAR database" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if year < min(possyears):
        print("%s is earlier than min year in EDGAR database" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
        
    
    #Read in EDGAR data of annual mean CH4 emissions
    #units are in kg/m2/s
    edgar = edpath+'v432_'+soi+'_'+str(year)+'.0.1x0.1.nc'
    
    #Species molar mass
    speciesmm = molar_mass(soi)
#    if soi == 'CH4':
#        speciesmm = 16.0425
#    elif soi == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % soi
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
        
    
    ds = xr.open_dataset(edgar)
    soiname = 'emi_'+soi.lower()    
    tot = ds[soiname].values*1e3/speciesmm
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    mtohe = lon_in > 180
    lon_in[mtohe] = lon_in[mtohe] - 360 
    ordinds = np.argsort(lon_in)
    lon_in = lon_in[ordinds]
    tot = tot[:, ordinds] 
    
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
    
def getbloomwetlandsCH4(year, lon_out, lat_out, timeframe="monthly"):
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
        timeframe (str, optional) :
            Should the data be averaged by month or by day? 
            "monthly" = monthly averaged
            "daily" = daily averaged
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    """
    path = '/data/shared/GAUGE/CH4/'
    bloomwetlands = 'CH4_flux_wetlands_rice_global_2003_2009.nc'
    
    possyears = np.arange(7) + 2003
    if year > max(possyears):
        print("%s is later than max year in Bloom et al wetlands data" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if year < min(possyears):
        print("%s is earlier than min year in Bloom et al data" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
    
    ds = xr.open_dataset(path+bloomwetlands)
    
    d = ds.date.values.astype(str)
    ddt = np.empty_like(d)
    for i in range(len(ddt)):
       ddt[i] = datetime.datetime.strptime(d[i], '%Y%m%d')
    ds.date.values = pd.to_datetime(ddt)
    ds = ds.sel(date=str(year))
    
    #if monthly == True:
    if timeframe == "monthly":
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
        Currently this is for CH4 and N2O: 2012, 2013, 2015.
        Use getNAEI for this. 
        
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
    
def getNAEI(year, lon_out, lat_out, soi, naei_sector):
 
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
        naei_sector (str):
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
        print('Raw NAEI file LatLong1km_'+soi.upper()+'_'+str(year)+'.csv.gz does not exist')
        print('Returning None')
        return None
    sectorlist = ["energyprod","domcom","indcom","indproc","offshore","roadtrans",
    "othertrans","waste","agric","nature","points","total","totalexcship"]
    if naei_sector not in sectorlist:
        print('Sector not one of:')
        print(sectorlist)
        print('Returning None')
        return None
        
    #Read emissions as tonnes/km^2/yr and naei lat and lons
    df = pd.DataFrame.from_csv(fn)   
    lat = np.asarray(df.Latitude)
    lon = np.asarray(df.Longitude)    
    emissions = np.asarray(df[naei_sector])
    
    #Make a square grid for the emissions
    latarr = np.arange(min(lat), max(lat), 0.01)
    lonarr = np.arange(min(lon), max(lon), 0.01)
    grdemis = np.zeros((len(latarr), len(lonarr)))

    for i in range(len(emissions)):
        ilat = np.where(abs(latarr-lat[i]) == np.min(abs(latarr - lat[i])) )
        ilon = np.where(abs(lonarr-lon[i]) == np.min(abs(lonarr - lon[i])) )
        grdemis[ilat, ilon] = emissions[i]

    #Convert to mol/m2/s
    speciesmm = molar_mass(soi)     
#    if soi == 'ch4':
#        speciesmm = 16.0425
#    if soi == 'n2o':
#        speciesmm = 44.013        
    if year % 4 == 0:
        diy = 365
    else:
        diy = 366    
    grdemis = grdemis / (diy * 3600*24) / speciesmm
    
    #Regrid to desired lats and lons
    narr, reg = regrid2d(grdemis, latarr, lonarr,
                                 lat_out, lon_out)

    return(narr)

def getedgarannualsectors(year, lon_out, lat_out, edgar_sectors, soi='CH4'):
    """
    Get annual emission totals for species of interest from EDGAR v4.3.2 data
    for sector or sectors.
    Regrids to the desired lats and lons.
    If there is no data for the species you are looking at you may have to 
    download it from: 
    http://edgar.jrc.ec.europa.eu/overview.php?v=432_GHG&SECURE=123
    and placed in: 
    /data/shared/Gridded_fluxes/<species>/EDGAR_v4.3.2/<species>_sector_yearly/
    
    Args:
        year (int): 
            Year of interest
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        edgar_sectors (list):
            List of strings of EDGAR sectors to get emissions for.
            These will be combined to make one array.
            See 'Notes' for names of sectors
        soi (str):
            Which species you want to look at. 
            e.g. soi = 'CH4'
            Default = 'CH4'
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon]
        
    Notes:
        Names of EDGAR sectors:
           'powerindustry'; 
           'oilrefineriesandtransformationindustry'; 
           'combustionformanufacturing'; 
           'aviationclimbinganddescent';  
           'aviationcruise'; 
           'aviationlandingandtakeoff';  
           'aviationsupersonic'; 
           'roadtransport'; 
           'railwayspipelinesandoffroadtransport'; 
           'shipping';  
           'energyforbuildings';  
           'fuelexploitation'; 
           'nonmetallicmineralsproduction';  
           'chemicalprocesses';
           'ironandsteelproduction'; 
           'nonferrousmetalsproduction'; 
           'nonenergyuseoffuels'; 
           'solventsandproductsuse'; 
           'entericfermentation'; 
           'manuremanagement';  
           'agriculturalsoils';  
           'indirectN2Oemissionsfromagriculture'; 
           'agriculturalwasteburning';  
           'solidwastelandfills';  
           'wastewaterhandling';  
           'Solid waste incineration';  
           'fossilfuelfires'; 
           'indirectemissionsfromNOxandNH3';  
    """
    
    
    soi = soi.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = '/data/shared/Gridded_fluxes/'+soi+'/EDGAR_v4.3.2/'+soi+'_sector_yearly/'
    
    #Dictionary of codes for sectors
    secdict = {'powerindustry' : '1A1a', 
               'oilrefineriesandtransformationindustry' : '1A1b_1A1c_1A5b1_1B1b_1B2a5_1B2a6_1B2b5_2C1b',
               'combustionformanufacturing' : '1A2',
               'aviationclimbinganddescent' : '1A3a_CDS',
               'aviationcruise' : '1A3a_CRS',
               'aviationlandingandtakeoff' : '1A3a_LTO',
               'aviationsupersonic' : '1A3a_SPS',
               'roadtransport' : '1A3b',
               'railwayspipelinesandoffroadtransport' : '1A3c_1A3e',
               'shipping' : '1A3d_1C2',
               'energyforbuildings' : '1A4',
               'fuelexploitation' : '1B1a_1B2a1_1B2a2_1B2a3_1B2a4_1B2c',
               'nonmetallicmineralsproduction' : '2A',
               'chemicalprocesses': '2B',
               'ironandsteelproduction' : '2C1a_2C1c_2C1d_2C1e_2C1f_2C2',
               'nonferrousmetalsproduction' : '2C3_2C4_2C5',
               'nonenergyuseoffuels' : '2G',
               'solventsandproductsuse' :  '3',
               'entericfermentation' : '4A',
               'manuremanagement' : '4B',
               'agriculturalsoils' : '4C_4D',
               'indirectN2Oemissionsfromagriculture' : '4D3',
               'agriculturalwasteburning' : '4F',
               'solidwastelandfills' : '6A_6D',
               'wastewaterhandling' : '6B',
               'Solid waste incineration' : '6C',
               'fossilfuelfires' : '7A',
               'indirectemissionsfromNOxandNH3' : '7B_7C'           
    } 
    
    #Check to see range of years. If desired year falls outside of this range 
    #then take closest year
    possyears = np.empty(shape=[0,0],dtype=int)
    for f in glob.glob(edpath+'v432_'+soi+'_*'):
        fname = f.split('/')[-1]
        fyear = fname[9:13]      #Extract year from filename
        possyears = np.append(possyears, int(fyear))
    if year > max(possyears):
        print("%s is later than max year in EDGAR database" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if year < min(possyears):
        print("%s is earlier than min year in EDGAR database" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
    
        
    #Species molar mass
    speciesmm = molar_mass(soi)
#    if soi == 'CH4':
#        #speciesmm = 16.0425
#        speciesmm = molar_mass(soi)
#    elif soi == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % soi
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
    
    
    #Read in EDGAR data of annual mean CH4 emissions for each sector
    #These are summed together
    #units are in kg/m2/s
    tot = None
    for sec in edgar_sectors:
        edgar = edpath+'v432_'+soi+'_'+str(year)+'_IPCC_'+secdict[sec]+'.0.1x0.1.nc'    
        if os.path.isfile(edgar):
            ds = xr.open_dataset(edgar)
            soiname = 'emi_'+soi.lower()
            if tot is None:
                tot = ds[soiname].values*1e3/speciesmm
            else:
                tot += ds[soiname].values*1e3/speciesmm
        else:
            print('No annual file for sector %s and %s' % (sec, soi))
        
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon))    
    narr, reg = regrid2d(tot, lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)   

def getedgarmonthlysectors(lon_out, lat_out, edgar_sectors, months=[1,2,3,4,5,6,7,8,9,10,11,12],
                           soi='CH4'):
    """
    Get 2010 monthly emissions for species of interest from EDGAR v4.3.2 data
    for sector or sectors.
    Regrids to the desired lats and lons.
    If there is no data for the species you are looking at you may have to 
    download it from: 
    http://edgar.jrc.ec.europa.eu/overview.php?v=432_GHG&SECURE=123
    and place it in: 
    /data/shared/Gridded_fluxes/<species>/EDGAR_v4.3.2/<species>_sector_monthly/
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        edgar_sectors (list):
            List of strings of EDGAR sectors to get emissions for.
            These will be combined to make one array.
            See 'Notes' for names of sectors
        months (list of int; optional): 
            Desired months.
        soi (str, optional):
            Which species you want to look at. 
            e.g. soi = 'CH4'
            Default = 'CH4'
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [no of months, lat, lon]
        
    Notes:
        Names of EDGAR sectors:
           'powerindustry'; 
           'oilrefineriesandtransformationindustry'; 
           'combustionformanufacturing'; 
           'aviationclimbinganddescent';  
           'aviationcruise'; 
           'aviationlandingandtakeoff';  
           'aviationsupersonic'; 
           'roadtransport'; 
           'railwayspipelinesandoffroadtransport'; 
           'shipping';  
           'energyforbuildings';  
           'fuelexploitation'; 
           'nonmetallicmineralsproduction';  
           'chemicalprocesses';
           'ironandsteelproduction'; 
           'nonferrousmetalsproduction'; 
           'nonenergyuseoffuels'; 
           'solventsandproductsuse'; 
           'entericfermentation'; 
           'manuremanagement';  
           'agriculturalsoils';  
           'indirectN2Oemissionsfromagriculture'; 
           'agriculturalwasteburning';  
           'solidwastelandfills';  
           'wastewaterhandling';  
           'Solid waste incineration';  
           'fossilfuelfires'; 
           'indirectemissionsfromNOxandNH3';  
    """
    soi = soi.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = '/data/shared/Gridded_fluxes/'+soi+'/EDGAR_v4.3.2/'+soi+'_sector_monthly/'
    
    #Dictionary of codes for sectors
    secdict = {'powerindustry' : '1A1a', 
               'oilrefineriesandtransformationindustry' : '1A1b_1A1c_1A5b1_1B1b_1B2a5_1B2a6_1B2b5_2C1b',
               'combustionformanufacturing' : '1A2',
               'aviationclimbinganddescent' : '1A3a_CDS',
               'aviationcruise' : '1A3a_CRS',
               'aviationlandingandtakeoff' : '1A3a_LTO',
               'aviationsupersonic' : '1A3a_SPS',
               'roadtransport' : '1A3b',
               'railwayspipelinesandoffroadtransport' : '1A3c_1A3e',
               'shipping' : '1A3d_1C2',
               'energyforbuildings' : '1A4',
               'fuelexploitation' : '1B1a_1B2a1_1B2a2_1B2a3_1B2a4_1B2c',
               'nonmetallicmineralsproduction' : '2A',
               'chemicalprocesses': '2B',
               'ironandsteelproduction' : '2C1a_2C1c_2C1d_2C1e_2C1f_2C2',
               'nonferrousmetalsproduction' : '2C3_2C4_2C5',
               'nonenergyuseoffuels' : '2G',
               'solventsandproductsuse' :  '3',
               'entericfermentation' : '4A',
               'manuremanagement' : '4B',
               'agriculturalsoils' : '4C_4D',
               'indirectN2Oemissionsfromagriculture' : '4D3',
               'agriculturalwasteburning' : '4F',
               'solidwastelandfills' : '6A_6D',
               'wastewaterhandling' : '6B',
               'Solid waste incineration' : '6C',
               'fossilfuelfires' : '7A',
               'indirectemissionsfromNOxandNH3' : '7B_7C'           
    }
    
    print('Note that the only year for monthly emissions is 2010 so using that.')
        
    #Species molar mass
    speciesmm = molar_mass(soi)
#    if soi == 'CH4':
#        speciesmm = 16.0425
#    elif soi == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % soi
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
    
    
    #Read in EDGAR data of annual mean CH4 emissions for each sector
    #These are summed together
    #units are in kg/m2/s
    first = 0
    for month in months:
        tot = None
        for sec in edgar_sectors:
            edgar = edpath+'v432_'+soi+'_2010_'+str(month)+'_IPCC_'+secdict[sec]+'.0.1x0.1.nc'    
            if os.path.isfile(edgar):
                ds = xr.open_dataset(edgar)
                soiname = 'emi_'+soi.lower()
                if tot == None:
                    tot = ds[soiname].values*1e3/speciesmm
                else:
                    tot += ds[soiname].values*1e3/speciesmm
            else:
                print('No monthly file for sector %s' % sec)
        
            if first == 0:
                emissions = np.zeros((len(months), tot.shape[0], tot.shape[1]))
                emissions[0,:,:] = tot
            else:
                first += 1
                emissions[first,:,:] = tot
                
            
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon, len(months)))
       
    for i in range(len(months)):
       narr[:,:,i], reg = regrid2d(emissions[i,:,:], lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)

def _JULESfile(year):
    '''
    The _JULESfile function opens the correct JULES wetland file for the given
    year (int) as an xarray.Dataset object.
    '''
    path = os.path.join(data_path,"Gridded_fluxes/CH4/JULES")
    filename_jules = os.path.join(path,"u-ax751_ch4_{}.nc.gz".format(year))

    return filename_jules

def _SWAMPSfile():
    '''
    The _SWAMPSfile function open the correct SWAMPS wetland fraction file for
    the given year (int) as an xarray.Dataset object.
    '''
    path = os.path.join(data_path,"Gridded_fluxes/CH4/JULES")
    #filename_swamps = os.path.join(path,"fw_swamps-glwd_2000-2012.nc") # Previous file
    filename_swamps = os.path.join(path,"gcp-ch4_wetlands_2000-2017_05deg.nc")
    
    return filename_swamps

def _readJULESwetland(year,species="CH4"):
    '''
    The _readJULESwetland function reads and interprets the JULES wetland maps for the 
    specified year.
    
    Note: converts input "month" dimension into "time" co-ordinates containing datetime objects.
    
    Args:
        year (int) :
            Year of interest.
        species (str, optional) :
            Species of interest. At the moment this should just be "CH4"
            Default = "CH4".
    
    Returns:
        xarray.Dataset :
            Dataset of JULES flux in expected format.
    '''
    
    if year < 2009 or year > 2016:
        print("No JULES wetlands data outside the range 2009-2016 for now.")
        return None

    filename_jules = _JULESfile(year)
    flux_jules = xr.open_dataset(filename_jules)

    ## For flux_jules, reassign month unit as datetime unit
    flux_base_time = dt.datetime(year=year,month=1,day=1)
    flux_time = [flux_base_time+relativedelta(months=int(num)) for num in flux_jules.month.values]
    flux_time = np.array([np.datetime64(date) for date in flux_time])
    
    flux_jules = flux_jules.assign_coords(**{"time":("month",flux_time)})
    flux_jules = flux_jules.swap_dims({"month":"time"})

    return flux_jules

def _SWAMPSwetlandArea(year,lon_out,lat_out,month=1):
    '''
    The _SWAMPSwetlandArea function calculates the area of the wetland extent within the specified
    latitude and longitude grid and the fraction of the total global wetlands area.
    
    This is based on the wetlands fraction input for SWAMPS across a global grid.
    
    Args:
        year (int) :
            Year of interest. Should be between 2000-2012 (at the moment).
        lat_out (numpy.array) :
            Latitude grid.
        lon_out (numpy.array) :
            Longitude grid.
        month (int, optional) :
            Month to extract this area for. Should be between 1 and 12.
            Default = 1 (i.e. January)
    
    Returns:
        tuple (float,float) :
            wetland area in m2 within specified latitude and longitude grid, fraction of global
            wetland area.
    '''
    if month >= 1 and month <= 12:
        month_id = month-1
    else:
        raise ValueError("Did not recognise month input: {}. Expect value between 1 and 12.".format(month))
    
    frac_swamps = xr.open_dataset(_SWAMPSfile())
    fw = "Fw" # Name of variable within file

    area_swamps = frac_swamps[fw][month_id,:,:].values
    area_swamps = np.nan_to_num(area_swamps)

    grid_wetl = areagrid(frac_swamps.lat.values,frac_swamps.lon.values)
    wetl_area = np.sum(area_swamps*grid_wetl)
    
    frac_wetl_domain = regrid2d(area_swamps, frac_swamps.lat, frac_swamps.lon, lat_out, lon_out)[0]
    grid_domain = areagrid(lat_out,lon_out)
    wetl_area_domain = np.sum(frac_wetl_domain*grid_domain)

    return wetl_area_domain,wetl_area_domain/wetl_area

def getJULESwetlands(year,lon_out,lat_out,soi="CH4",scale_wetlands=True,total_emission=185e12):
    '''
    The getJULESwetlands function creates an emissions grid for wetlands based on JULES wetlands maps.
    Rather than using the modelled  JULES wetland fraction, this is scaled against the observed SWAMPS 
    wetland fraction instead.
    
    Args:
        year (int) :
           Year of interest. Should be between 2009-2012 (at the moment - overlap between available JULES and SWAMPS data). 
        lon_out (numpy.array) :
            Longitude grid.
        lat_out (numpy.array) :
            Latitude grid.
        soi (str, optional) :
            Species of interest. At the moment, this should be "CH4".
        scale_wetlands (bool, optional) :
            Whether to scale emissions by the wetland fraction of a total emissions value.
            Default = True.
        total_emissions (float, optional) :
            If scale_emissions=True, this is the global emissions assumed for wetlands methane emissions.
            The wetlands fraction within the output area will then be used to find the emissions based on
            total emissions.
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
            Default = 185e12 (value for total global wetlands emissions in 2012 from Saunois et al, 2016)
    
    Returns:
        numpy.array :
            Re-gridded emissions map with dimensions (lat,lon,time)
    '''
    
    if soi.upper() != "CH4":
        print("Unable to extract JULES wetland values for any species except CH4 (at the moment)")
        return None
    
    flux_jules = _readJULESwetland(year,soi)
    frac_swamps = xr.open_dataset(_SWAMPSfile())
    
    fch4_name = "fch4_wetl_npp"
    fwetl_name = "fwetl"
    fw_name = "Fw"
    
    # -1.0e30 used as fill value for JULES data - essentially zero pixels/nan
    # want to set these to 0.0
    fill_value = np.min(flux_jules[fch4_name].values) # This may incorrect if input file is changed and different fill value is specified
    fch4_fill_indices = np.where(flux_jules[fch4_name]==fill_value)
    
    flux_jules.fch4_wetl_npp.values[fch4_fill_indices] = 0.0

    flux_jules_frac = np.abs(flux_jules[fch4_name] / flux_jules[fwetl_name])
    flux_jules_frac.values = np.nan_to_num(flux_jules_frac) # Any number divided by 0.0 will be nan, so change these back to 0.0
    
    ## Multiply by fractions from SWAMPS to rescale to measured rather than simulated inundation area
    frac_swamps[fw_name].values = np.nan_to_num(frac_swamps[fw_name])
    frac_reindex = frac_swamps.reindex_like(flux_jules_frac,method="ffill")

    fch4_wetl_npp_new = flux_jules_frac*frac_reindex[fw_name]

    lat = fch4_wetl_npp_new.lat.values
    lon = fch4_wetl_npp_new.lon.values
    emissions = np.moveaxis(fch4_wetl_npp_new.values,0,2) # re-arrange from [time,lat,lon] to [lat,lon,time]

    ## Regrid
    nlat_out = len(lat_out)
    nlon_out = len(lon_out)
    nt = len(fch4_wetl_npp_new.time)
    
    narr = np.zeros((nlat_out, nlon_out, nt))
   
    for i in range(nt):
        narr[:,:,i], reg = regrid2d(emissions[:,:,i], lat, lon, lat_out, lon_out)

    # May also want to rescale wetlands as JULES output will be too low.
    # E.g. 185 Tg/yr for total bottom up wetlands emissions from Saunois et al, 2016
    # Work out fraction of wetland area within new lat-lon grid and multiply by total.
    if scale_wetlands:
        for i in range(nt):
            frac = _SWAMPSwetlandArea(year,lon_out,lat_out,month=i+1)[1]
            scale = round((total_emission*frac)/1e12,1)*1e12
            print("{:02}) Scaling total JULES wetlands emissions within domain to: {} g/yr".format(i+1,scale))
            narr[:,:,i] = scale_emissions(narr[:,:,i],soi,lat_out,lon_out,total_emission=scale)

    return narr


def scale_emissions(emissions_t,soi,lat,lon,total_emission):
    '''
    The scale_emissions function scales emissions values for one time grid to a 
    total_emission value.
    
    Args:
        emissions_t (numpy.array) :
            Emissions values for one time point. Expect array to have dimensions nlat x nlon
        soi (str) :
            Species of interest. Used to extract molar mass. e.g. "CH4"
        lat (numpy.array) :
            Latitude grid for emissions.
        lon (numpy.array) :
            Longitude grid for emissions.
        total_emission (float) :
            Total emissions of the output array. 
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
    
    Returns:
        numpy.array :
            Scaled emissions array (nlat x nlon x nt)
    '''

    # Calculate number of moles
    molmass = molar_mass(soi)
    
    total_time = 365.*3600.*24. # seconds in a year

    areas = areagrid(lat,lon)
    
    current_emission = np.sum(emissions_t*areas)*total_time*molmass
    scaling = total_emission/current_emission
    print("Current emissions total: {} g/yr (scaling needed to {} g/yr: {})".format(current_emission,total_emission,scaling))
    
    emissions_new = np.copy(emissions_t)*scaling
    #print "New total emissions: {} g/yr".format(np.sum(emissions_new*areas)*total_time*molmass)
    
    return emissions_new

def scale_emissions_all(emissions,soi,lat,lon,total_emission):
    '''
    The scale_emissions_all function scales emissions values within all time grids to the same 
    total_emission value.
    
    Args:
        emissions (numpy.array) :
            Emissions values across multiple times. Expects array to have dimensions nlat x nlon x nt
        soi (str) :
            Species of interest. Used to extract molar mass. e.g. "CH4"
        lat (numpy.array) :
            Latitude grid for emissions.
        lon (numpy.array) :
            Longitude grid for emissions.
        total_emission (float) :
            Total emissions of the output array. 
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
    
    Returns:
        numpy.array :
            Scaled emissions array (nlat x nlon x nt)
    '''

    nt = emissions.shape[2]
    
    emissions_new = np.copy(emissions)
    
    for i in range(nt):
        emissions_i = emissions[:,:,i]
        scale_emissions(emissions_i,soi,lat,lon,total_emission)        
    
    return emissions_new

def plot_emissions(emissions_t,lat,lon,fignum=None,cmap="inferno_r",show=True):
    '''
    The plot_emissions function plots emissions values as a contour plot on a map specified by the
    latitude and longitude extent.
    
    Args:
        emissions_t (numpy.array) :
            Emissions values for one time point. Expect array to have dimensions nlat x nlon
        lat (numpy.array) :
            Latitude extent of emissions in degrees. Can be just lower and upper bounds or whole 
            array of latitude values.
        lon (numpy.array) :
            Longitude extent of emissions in degrees. Can be just lower and upper bounds or whole 
            array of longitude values.
        fignum (int, optional) :
            Figure number for plot.
            Default = None
        cmap (str, optional) :
            Colour map to use for plotting. See matplotlib.org/users/colormaps.html for some options.
            Default = "inferno_r"
        show (bool, optional) :
            Whether to immediately flush the buffer and produce the plot.
            Default = True
    
    Returns:
        cartopy axis object:
            Axis object for plot so far
        
        If show is True:
            Plots emissions countour map.
    '''
    
    if fignum:
        plt.figure(fignum)
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent((lon[0],lon[-1],lat[0],lat[-1]),crs=ccrs.PlateCarree())
    ax.coastlines()
    
    plt.contourf(lon,lat,emissions_t, 60,transform=ccrs.PlateCarree(),cmap=cmap)
    plt.colorbar(orientation="horizontal",pad=0.05)
    
    if show:
        plt.show()
    
    return ax