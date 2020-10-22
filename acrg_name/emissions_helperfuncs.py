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
from __future__ import division

from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import xarray as xray
import glob 
import h5py
import os
import sys
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
from acrg_countrymask import domain_volume
from acrg_name import flux
import pyproj
import acrg_grid.regrid_xesmf
import getpass

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    data_path = os.getenv("DATA_PATH") 
else:
    from acrg_config.paths import paths
    data_path = paths.data

output_directory = os.path.join(data_path,"LPDM/emissions/")

def loadAscii(filename):
    '''
    Load an ASCII format grid of emissions from file

    Parameters
    ----------
    filename : Name of ascii grid to read in

    Returns
    -------
    output : xarray Dataset containing the data read from the file

    '''
    data = np.loadtxt(filename, skiprows = 6)
    data = np.flip(data, axis=0)
    header = np.loadtxt(filename,usecols=(1))[:6]
    
    #grab 'center' or 'corner' from header information
    cell_type = np.loadtxt(filename,usecols=(0),dtype=str)[2][3:]
    
    data[data == float(header[5])] = np.nan
    
    xValues = np.arange(header[0]) * header[4] + header[2]
    yValues = np.arange(header[1]) * header[4] + header[3]
    
    output = xr.Dataset({})
    output["lat"] = xr.DataArray(yValues,dims="lat")
    output["lon"] = xr.DataArray(xValues,dims="lon")
    output["data"] = xr.DataArray(data,dims=["lat", "lon"])
    output.attrs["celltype"] = cell_type
    output.attrs["cellsize"] = header[4]
    
    return output

def getGFED(year, lon_out, lat_out, timeframe='monthly', months = [1,2,3,4,5,6,7,8,9,10,11,12], species='CH4', incagr=False):
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
        species (str):
            Which species you want to look at as defined in 
            GFED4_Emission_Factors.csv.
            e.g. species = 'CO2'
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
       sourceindex = [7,5,3,1,11,9] #/data/
        
    
    directory    = os.path.join(data_path,'Gridded_fluxes/GFED4/fire_emissions_v4_R1_1293/data')
    species = species.upper()
    
    ###
    ###  Read in emission factors
    ###
    species_mult = [] # names of the different gas and aerosol species
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
            species_mult.append(contents[0])
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
                species_mult.append(contents[0])
                xvals = [i for i,x in enumerate(contents) if (x=='x') or (x=='')]
                for j in xvals:
                    contents[j] = 0
                EFs[k,:] = np.asarray(contents[1:])
                k += 1
        break
          
    f.close()
    
    # Find species of interest in file
    kc = 0
    for k in species_mult:
        if k.find(species) > -1:
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
    possyears.sort()
    #print "Possible years to use for GFED data: {}".format(possyears)
    if int(year) > max(possyears):
        print("%s is later than latest year in GFED v4.1 database" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if int(year) < min(possyears):
        print("%s is earlier than earliest year in GFED v4.1 database" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
    
    #Get emissions
    if year < 2017:
        string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    else:
        print("WARNING: GFED data for 2017+ is in beta.")
        string = directory+'/GFED4.1s_'+str(year)+'_beta.hdf5'
    #string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
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
            # read in DM emissions
           string = '/emissions/'+months_str[month]+'/DM'
           DM_emissions = f[string][:]
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
                   contribution = f[string][:] #rt17603: Moved to loop as contribution wasn't being assigned per source.
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
    
def getedgarannualtotals(year, lon_out, lat_out, species='CH4'):
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
        species (str):
            Which species you want to look at. 
            e.g. species = 'CH4'
            Default = 'CH4'
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon]
    
    """
    
    species = species.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = os.path.join(data_path,'Gridded_fluxes/'+species+'/EDGAR_v4.3.2/v432_'+species+'_TOTALS_nc/')
    
    #Check to see range of years. If desired year falls outside of this range 
    #then take closest year
    possyears = np.empty(shape=[0,0],dtype=int)
    for f in glob.glob(edpath+'v432_'+species+'_*'):
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
    edgar = edpath+'v432_'+species+'_'+str(year)+'.0.1x0.1.nc'
    
    #Species molar mass
    speciesmm = molar_mass(species)
#    if species == 'CH4':
#        speciesmm = 16.0425
#    elif species == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % species
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
        
    
    ds = xr.open_dataset(edgar)
    soiname = 'emi_'+species.lower()    
    tot = old_div(ds[soiname].values*1e3,speciesmm)
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
    path = os.path.join(data_path,'GAUGE/CH4/')
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
    path = os.path.join(data_path,'GAUGE/CH4/')
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
    print('There is a newer (much better) Bloom wetlands model!!!')
    print("Use the function getBloom2017() instead")
    
    path = os.path.join(data_path,'GAUGE/CH4/')
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
       ddt[i] = datetime.datetime.strptime(str(d[i]), '%Y%m%d')
    #ds.date.values = pd.to_datetime(ddt) # Unable to assign back to values in this way in later xarray versions
    ds = ds.assign_coords(date=pd.to_datetime(ddt))
    ds = ds.sel(date=str(year))
    
    #if monthly == True:
    if timeframe == "monthly":
        ds = ds.resample(date='M').mean()
    
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
    
def getBloom2017(year, lon_out, lat_out, modeltype='extended'):
    """
    Global wetlands CH4 emissions from Bloom et al 2017, regridded to desired 
    lats and lons for year of interest.
    This function takes the mean from an ensemle, but individuals models are 
    available in the Gridded_fluxes folder.
    
    Args:
        year (int):
            Year of interest.
        lon_out (array): 
            Longitudes to output the data on.
        lat_out (array):
            Latitudes to output the data on.
        modeltype (str):
            Which ensembel model to use. Options are 'extended' or 'full'.
            'full' uses more ensemble members
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    """    
    path = os.path.join(data_path,'Gridded_fluxes/CH4/Bloom2017/')
    if modeltype == 'extended':
        bloomwetlands = 'WetCHARTs_extended_ensemble_mean.nc4'
        possyears = np.arange(15) + 2001
    elif modeltype == 'full':
        bloomwetlands = 'WetCHARTs_full_ensemble_mean.nc4'
        possyears = np.arange(2) + 2009
    else:
        raise Exception('There is no model type %s. Using "extended" instead' % modeltype)
        bloomwetlands = 'WetCHARTs_extended_ensemble_mean.nc4'
        possyears = np.arange(15) + 2001
    
    
    if year > max(possyears):
        print("%s is later than max year in Bloom 2017 wetlands data" % str(year))
        print("Using %s as the closest year" % str(max((possyears))))
        year = max(possyears)
    if year < min(possyears):
        print("%s is earlier than min year in Bloom 2017 data" % str(year))
        print("Using %s as the closest year" % str(min((possyears))))
        year = min(possyears)
    
    yearloc = np.where(possyears == year)[0][0]
    
    with xr.open_dataset(path+bloomwetlands) as load:
        ds = load.load()   
    
    bloomch4 = ds.ch4.values[:,:,yearloc*12:yearloc*12+12]
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = 12
    
    narr = np.zeros((nlat, nlon, nt))
    
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(bloomch4[:,:,i], lat_in, lon_in,
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
    
def get_naei_public(filename, lon_out, lat_out):
    '''
    Regrid the publically available NAEI .asc files that can be downloaded from https://naei.beis.gov.uk/data/map-uk-das

    Parameters
    ----------
    filename : string
        Name of the .asc file to read in for naei data
    lon_out : array like
        longitude values for the output grid
    lat_out : array like
        latitude values for the output grid

    Returns
    -------
    output : xarray dataset
        Regridded emissions data

    '''
    #OSGB projection for NAEI
    OS_proj = pyproj.Proj(init='epsg:7405')
    #lat-lon projection
    ll_proj = pyproj.Proj(init='epsg:4326')
    
    raw_data = loadAscii(filename)
    
    #ensure the grid is cell-centered:
    if raw_data.attrs["celltype"].lower() == "corner":
        raw_data["lat"] = raw_data.lat + float(raw_data.attrs["cellsize"])/2
        raw_data["lon"] = raw_data.lon + float(raw_data.attrs["cellsize"])/2
        
    #convert the OSGB grid to lat-lon
    XX, YY = np.meshgrid(raw_data["lon"].values,raw_data["lat"].values)
    XX, YY = pyproj.transform(OS_proj,ll_proj, XX, YY)
    #create boundary mesh
    XX_b, YY_b = acrg_grid.regrid_xesmf.getGridCC(raw_data["lon"].values,raw_data["lat"].values)
    XX_b, YY_b = pyproj.transform(OS_proj,ll_proj, XX_b, YY_b)
    
    input_grid = xr.Dataset({'lon': (['x', 'y'], XX),
                             'lat': (['x', 'y'], YY),
                             'lon_b': (['x_b', 'y_b'], XX_b),
                             'lat_b': (['x_b', 'y_b'], YY_b)})
    
    output_grid = acrg_grid.regrid_xesmf.create_xesmf_grid_uniform_cc(lat_out, lon_out)
    
    #conda version of xesmf does not handle nans yet
    data_sans_nans = np.nan_to_num(raw_data.data, copy=True)
    regridded_data = acrg_grid.regrid_xesmf.regrid_betweenGrids(data_sans_nans, input_grid, output_grid)

    
    output = xr.Dataset({"data":(["lat", "lon"],regridded_data)},
                        coords={"lat":(["lat"], lat_out),
                                "lon":(["lon"], lon_out)})
    
    return output

def process_naei_public(naei_dir, species, year, sector, lon_out, lat_out, output_path = None):
    """

    Parameters
    ----------
    naei_dir : str
        Directory containing unzupped naei .asc files
    species : str
        species as given in naei filenames
    year : str
        Two digit year (such as 17 for 2017)
    sector : str
        Sector name as given in naei filenames
    lon_out : array like
        longitude values for the output grid
    lat_out : array like
        latitude values for the output grid
    output_path : str, optional
        The filename to write the dataset to. The default is None, which does not output to a file

    Returns
    -------
    output : xarray dataset
        Regridded emissions data

    """
    
    #run soome checks
    sectorlist = ["energyprod","domcom","indcom","indproc","offshore","roadtrans",
    "othertrans","waste","agric","nature","points","total","totarea"]
    if sector not in sectorlist:
        print('Sector {} not one of:'.format(sector))
        print(sectorlist)
        print('Returning None')
        return None
    filepath = os.path.join(naei_dir, sector+species+year+".asc")
    
    if not os.path.isfile(filepath):
        print('Raw NAEI file {} does not exist'.format(filepath))
        print('Returning None')
        return None
    
    naei_ds = get_naei_public(filepath, lon_out, lat_out)
    
    #Convert to mol/m2/s
    speciesmm = molar_mass(species)     
       
    if pd.to_datetime("20"+year).is_leap_year:
        days = 366
    else:
        days = 365   
        
    #add attributes    
    naei_ds["data"] = naei_ds["data"] / ((days * 3600*24) * speciesmm)
    naei_ds["data"].attrs["units"] = "mol/m2/s"
    naei_ds.attrs["Processed by"] = f"{getpass.getuser()}@bristol.ac.uk"
    naei_ds.attrs["Processed on"] = str(pd.Timestamp.now(tz="UTC"))
    naei_ds.attrs["Raw file used"] = filepath
    
    if output_path is not None:
        naei_ds.to_netcdf(output_path)
        
    return naei_ds
    
def processUKGHG(ukghg_dir,species,year,sector,lon_out,lat_out,output_path=None):
    """
    Processes UKGHG .nc files to resolution provided by lon_out and lat_out.
    Adds timestamp, author and raw file info to the output dataset.
    
    Args:
        ukghg_dir (str):
            Directory containing ukghg LonLat .nc files at 0.01 km resolution
        species (str):
            Species as given in ukghg filenames
        year (str):
            Year to use (e.g."2015")
        sector (str):
            Sector name as given in ukghg filenames. 
            Options: "agric","domcom","energyprod",'"indcom","indproc","natural",
                     "offshore","othertrans","roadtrans","solvents","total",'waste"
        lon_out (array):
            Longitude values for the output grid
        lat_out (array):
            Latitude values for the output grid
        output_path (str), optional:
            Filepath to write output dataset to. Default is None, which does not save output.

    Returns:
        ukghg_ds (xarray dataset):
            Regridded emissions data

    """
    
    #run soome checks
    sectorlist = ["agric","energyprod","domcom","indcom","indproc","natural","offshore",
                  "othertrans","roadtrans","solvents","total","waste"]
    if sector not in sectorlist:
        print('Sector {} not one of:'.format(sector))
        print(sectorlist)
        print('Returning None')
        return None
    
    filepath = os.path.join(ukghg_dir,
                            'uk_flux_'+sector+'_'+species+'_LonLat_0.01km_'+year+".nc")
    
    if not os.path.isfile(filepath):
        print('Raw UKGHG file {} does not exist'.format(filepath))
        print('Returning None')
        return None
    
    with xr.open_dataset(filepath) as file:
        flux = file[species+'_flux'].values[0,:,:]
        lat_in = file.latitude.values
        lon_in = file.longitude.values

    #regrid to NAME resolution
    regrid_data,reg = regrid2d(flux,lat_in,lon_in,lat_out,lon_out)
    
    ukghg_ds = xr.Dataset({"flux":(["lat", "lon"],regrid_data.data)},
                            coords={"lat":(["lat"], lat_out),
                                    "lon":(["lon"], lon_out)})
    
    ukghg_ds["flux"].attrs["units"] = "mol/m2/s"
    ukghg_ds.attrs["Processed by"] = f"{getpass.getuser()}@bristol.ac.uk"
    ukghg_ds.attrs["Processed on"] = str(pd.Timestamp.now(tz="UTC"))
    ukghg_ds.attrs["Raw file used"] = filepath
    
    if output_path is not None:
        ukghg_ds.to_netcdf(output_path+'.nc')
        print('Regridded UKGHG saved to {}'.format(output_path))
    
    return ukghg_ds

def get_US_EPA(epa_sectors=None,output_path=None):
    """
    Extracts yearly CH4 gridded emissions from EPA .nc files. 
    Sums emissions from given list of sectors.
    Regrids to USA domain lat/lons.
    
    Uses data files provided by Maasakkers, 2016 "A Gridded National Inventory 
    of U.S.Methane Emissions". Only 2012 available.
    
    Args:
        epa_sectors (list of str) (optional):
            If list, sum of sectors returned. 
            If None, total of all sectors returned.
            See below for list of possible sector names.
        output_path (str):
            Full output path and name to write output datatset to. Default is None, which
            does not save output.
        
    Returns:
        flux_ds (xarray dataset):
            Containing array of EPA fluxes, regridded to specified lat/lon.
            Dimensions of [lat,lon,time]
            
    Example:
        get_US_EPA("2012",["Enteric_Fermentation","Manure_Management"])
            
    """    

    EPAsectorlist = ["1A_Combustion_Mobile","1A_Combustion_Stationary",
                     "1B1a_Coal_Mining_Underground","1B1a_Coal_Mining_Surface",
                     "1B1a_Abandoned_Coal","1B2a_Petroleum","1B2b_Natural_Gas_Production",
                     "1B2b_Natural_Gas_Processing","1B2b_Natural_Gas_Transmission",
                     "1B2b_Natural_Gas_Distribution","2B5_Petrochemical_Production",
                     "2C2_Ferroalloy_Production","4A_Enteric_Fermentation",
                     "4B_Manure_Management","4C_Rice_Cultivation",
                     "4F_Field_Burning","5_Forest_Fires","6A_Landfills_Municipal",
                     "6A_Landfills_Industrial","6B_Wastewater_Treatment_Domestic",
                     "6B_Wastewater_Treatment_Industrial","6D_Composting"]
        
    #extract output lat/lons from fp file
    domain_vol = domain_volume("USA",os.path.join(data_path,"LPDM/fp_NAME"))
    lat_out = domain_vol[0]
    lon_out = domain_vol[1]
    
    #check for correct sector names
    if epa_sectors is not None:
    
        for EPAsector in epa_sectors:
            if EPAsector not in EPAsectorlist:
                print("EPA sector {0} is not one of: \n {1}".format(EPAsector,EPA_sectorlist))
                print("Returning None")
                return None
        
    epa_fp = os.path.join(data_path,"Gridded_fluxes/CH4/US_EPA/GEPA_Annual.nc")
    
    #epa flux in molec/cm2/s
    with xr.open_dataset(epa_fp) as epa_file:
        reference = epa_file.reference
        inv_ver = epa_file.inventory_version
        lat_in = epa_file.lat.values
        lon_in = epa_file.lon.values

        if epa_sectors is None:
            print("Calculating total of all sectors")
            sectors = list(epa_file.keys())
        else:
            print("Calculating total of given sectors")
            sectors = ["emissions_" + s for s in epa_sectors]
        
        for i,s in enumerate(sectors):
            if i == 0:
                total_flux = np.nan_to_num(x=epa_file[s].values,nan=0.)
            else:
                total_flux = np.add(total_flux,np.nan_to_num(x=epa_file[s].values,nan=0.))
    
    #epa flux in mol/m2/s
    print("Converting from molec/cm2/s to mol/m2/s")
    n_mol = 6.02214076e23
    
    total = (total_flux / n_mol ) * 1e4
    
    print("Regridding to USA domain lat lon")
    total_regrid,arr = regrid2d(total,lat_in,lon_in,lat_out,lon_out)
        
    flux_ds = xr.Dataset({"flux":(["lat", "lon","time"],np.expand_dims(total_regrid.data,2))},
                            coords={"lat":(["lat"], lat_out),
                                    "lon":(["lon"], lon_out),
                                    "time":np.array([np.datetime64("2012-01-01T00")])})
                                    
    flux_ds["flux"].attrs["units"] = "mol/m2/s"
    flux_ds.attrs["Processed by"] = f"{getpass.getuser()}@bristol.ac.uk"
    flux_ds.attrs["Processed on"] = str(pd.Timestamp.now(tz="UTC"))
    flux_ds.attrs["EPA sectors"] = epa_sectors
    flux_ds.attrs["reference"] = reference
    flux_ds.attrs["inventory_ver"] = inv_ver
    flux_ds.attrs["description"] = "Sectoral U.S. EPA CH4 emissions regridded to 'USA' domain"
    
    if output_path is not None:
        flux_ds.to_netcdf(output_path+'.nc')
        print('Regridded emissions data saved to {}'.format(output_path))
    
    return flux_ds
    

def getUKGHGandEDGAR(species,year,edgar_sectors=None,ukghg_sectors=None,output_path=None):
    """
    Extracts EDGAR and UKGHG gridded emissions data for the given year and sector.
    Regrids both datasets to the EUROPE domain and converts both the mol/m2/s.
    Embeds UKGHG emissions within the EDGAR map to create a grid of fluxes for the whole
    domain. For UK lat/lons (as defined by country-ukmo_EUROPE.nc) fluxes = UKGHG. Outside the UK fluxes = EDGAR.
    
    Option to include sources for one inventory. In this case, the sum of those sectors 
    only will be returned.
    
    CURRENTLY ONLY TESTED WITH EDGAR v5.0 CH4 2015 and UKGHG CH4 2015
        - DATA FILES FOR OTHER GASES/YEARS MAY NOT EXIST
    
    Args:
        year (str):
            Year of interest, e.g. '2015'
        edgar_sectors (list of str) (optional):
            EDGAR sectors to include. If list of values, the sum of these will be used.
            See below for list of possible sectors and full names.
        ukghg_sectors (list of str) (optional):
            UKGHG sectors to include. If list of values, the sum of these will be used.
            See below for list of possible sectors.
        output_path (str):
            Full output path and name to write output datatset to. Default is None, which
            does not save output.
            
    Returns:
        total_ds (xarray dataset):
            Containing array of regridded EDGAR fluxes with UKGHG embedded over UK lat/lons.
            Dimensions of [lat,lon,time]
            
    Example:
        getUKGHGandEDGAR("2015",["PRO_GAS","PRO_OIL"],["indcom","offshore"])
        
    Note:
        EDGAR sector names:
        "AGS" = Agricultural soils
        "AWB" = Agricultural waste burning
        "CHE" = Chemical processes
        "ENE" = Power industry
        "ENF" = Enteric fermentation
        "FFF" = Fossil fuel fires
        "IND" = Combustion for manufacturing
        "IRO" = Iron and steel production
        "MNM" = Maure management
        "PRO_COAL" = Fuel exploitation - coal
        "PRO_GAS" = Fuel exploitation - gas
        "PRO_OIL" = Fuel expoitation - oil
        "PRO" = Fuel exploitation - contains coal, oil, gas
        "RCO" = Energy for buildings
        "REF_TRF" = Oil refineries and transformational industries
        "SWD_INC" = Solid waste disposal - incineration
        "SWD_LDF" = Solid waste disposal - landfill
        "TNR_Aviation_CDS" = Aviation - climbing and descent
        "TNR_Aviation_CRS" = Aviation - cruise
        "TNR_Aviation_LTO" = Aviation - landing and takeoff 
        "TNR_Other" = Railways, pipelines and off-road transport
        "TNR_Ship" = Shipping
        "TRO" = Road transportation
        "WWT" = Waste water treatment
        
    """
        
    data_path = '/home/cv18710/work_shared/'
    
    edgarfp = os.path.join(data_path,"Gridded_fluxes",species.upper(),"EDGAR_v5.0/yearly_sectoral")
    ukghgfp = os.path.join(data_path,"Gridded_fluxes",species.upper(),"UKGHG")
    
    UKGHGsectorlist = ["agric","domcom","energyprod","indcom","indproc","natural","offshore",
                       "othertrans","roadtrans","total","waste"]
    
    EDGARsectorlist = ["AGS","AWB","CHE","ENE","ENF","FFF","IND","IRO","MNM",
                       "PRO_COAL","PRO_GAS","PRO_OIL","PRO","RCO","REF_TRF","SWD_INC",
                       "SWD_LDF","TNR_Aviation_CDS","TNR_Aviation_CRS",
                       "TNR_Aviation_LTO","TNR_Other","TNR_Ship","TRO","WWT"]
    
    #extract output lat/lons from fp file
    domain_vol = domain_volume("EUROPE",os.path.join(data_path,"LPDM/fp_NAME"))
    lat_out = domain_vol[0]
    lon_out = domain_vol[1]
    
    if edgar_sectors is not None:
        print('Including EDGAR sectors.')
    
        for EDGARsector in edgar_sectors:
            if EDGARsector not in EDGARsectorlist:
                print('EDGAR sector {0} not one of: \n {1}'.format(EDGARsector,EDGARsectorlist))
                print('Returning None')
                return None
            
        #edgar flux in kg/m2/s
        for i,sector in enumerate(edgar_sectors):
        
            edgarfn = "v50_" + species.upper() + "_" + year + "_" + sector + ".0.1x0.1.nc"

            with xr.open_dataset(os.path.join(edgarfp,edgarfn)) as edgar_file:
                edgar_flux = np.nan_to_num(edgar_file['emi_'+species.lower()].values,0.)
                edgar_lat = edgar_file.lat.values
                edgar_lon = edgar_file.lon.values

            if i == 0:
                edgar_total = edgar_flux
            else:
                edgar_total = np.add(edgar_total,edgar_flux)
            
        edgar_regrid_kg,arr = regrid2d(edgar_total,edgar_lat,edgar_lon,lat_out,lon_out)
    
        #edgar flux in mol/m2/s
        speciesmm = molar_mass(species)
        edgar_regrid = (edgar_regrid_kg.data*1e3) / speciesmm
        
        if ukghg_sectors is None:
            total_flux = edgar_regrid
            output_title = "EDGAR sectors regridded to EUROPE domain"
    
    #extract ukghg
    if ukghg_sectors is not None:
        print('Including UKGHG sectors.')
        
        for UKGHGsector in ukghg_sectors:
            if UKGHGsector not in UKGHGsectorlist:
                print('UKGHG sector {0} not one of: \n {1}'.format(UKGHGsector,UKGHGsectorlist))
                print('Returning None')
                return None    

        #ukghg flux in mol/m2/s
        for j,sector in enumerate(ukghg_sectors):

            ukghgfn = "uk_flux_" + sector + "_" + species.lower() + "_LonLat_0.01km_" + year + ".nc"

            with xr.open_dataset(os.path.join(ukghgfp,ukghgfn)) as ukghg_file:
                ukghg_flux = np.nan_to_num(ukghg_file[species.lower()+'_flux'].values[0,:,:],0.)
                ukghg_lat = ukghg_file.latitude.values
                ukghg_lon = ukghg_file.longitude.values

            if j == 0:
                ukghg_total = ukghg_flux
            else:
                ukghg_total = np.add(ukghg_total,ukghg_flux)

        ukghg_regrid,arr = regrid2d(ukghg_total,ukghg_lat,ukghg_lon,lat_out,lon_out)
   
        #if edgar, input ukghg over uk lat lons
        if edgar_sectors is not None:  
            print('Inserting UKGHG into EDGAR over UK lat/lons')
            
            with xr.open_dataset(os.path.join(data_path,'LPDM/countries/country-ukmo_EUROPE.nc')) as c_file:
                country = c_file['country'].values       

            for lat in range(edgar_regrid.shape[0]):
                for lon in range(edgar_regrid.shape[1]):
                    if country[lat,lon] == 19.0:
                        edgar_regrid[lat,lon] = ukghg_regrid.data[lat,lon]
            output_title = "EDGAR with UK lat/lons replaced with UKGHG values, regridded to EUROPE domain"
            total_flux = edgar_regrid

        else:
            total_flux = ukghg_regrid.data
            output_title = "UKGHG sectors regridded to EUROPE domain"
            
    #output to ds
    flux_ds = xr.Dataset({"flux":(["lat", "lon","time"],np.expand_dims(total_flux,2))},
                            coords={"lat":(["lat"], lat_out),
                                    "lon":(["lon"], lon_out),
                                    "time":np.array([np.datetime64(year+'-01-01T00')])})
                                    
    flux_ds["flux"].attrs["units"] = "mol/m2/s"
    flux_ds.attrs["Processed by"] = f"{getpass.getuser()}@bristol.ac.uk"
    flux_ds.attrs["Processed on"] = str(pd.Timestamp.now(tz="UTC"))
    flux_ds.attrs["EDGAR sectors"] = edgar_sectors
    flux_ds.attrs["UKGHG sectors"] = ukghg_sectors
    flux_ds.attrs["title"] = output_title
    
    if output_path is not None:
        flux_ds.to_netcdf(output_path+'.nc')
        print('Regridded emissions data saved to {}'.format(output_path))
    
    return flux_ds
            
def getNAEI(year, lon_out, lat_out, species, naei_sector):
 
    """
    Converts raw NAEI into gridded emissions data in mol/m2/s
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        year (int): 
            Year of interest
        species (str):
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
    
    species = species.lower()   
    naeidir = '/dagage2/agage/metoffice/naei/naei_raw_priors/'
    fn = naeidir+'LatLong1km_'+species.upper()+'_'+str(year)+'.csv.gz'
    #Checks
    if not os.path.isfile(fn):
        print('Raw NAEI file LatLong1km_'+species.upper()+'_'+str(year)+'.csv.gz does not exist')
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
    #df = pd.DataFrame.from_csv(fn) # This syntax is deprecated 
    df = pd.read_csv(fn)  
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
    speciesmm = molar_mass(species)     
#    if species == 'ch4':
#        speciesmm = 16.0425
#    if species == 'n2o':
#        speciesmm = 44.013        
    if year % 4 == 0:
        diy = 365
    else:
        diy = 366    
    grdemis = old_div(grdemis, (diy * 3600*24) * speciesmm)
    
    #Regrid to desired lats and lons
    narr, reg = regrid2d(grdemis, latarr, lonarr,
                                 lat_out, lon_out)

    return(narr)

def getedgarannualsectors(year, lon_out, lat_out, edgar_sectors, species='CH4'):
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
        species (str):
            Which species you want to look at. 
            e.g. species = 'CH4'
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
    
    
    species = species.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = os.path.join(data_path,'Gridded_fluxes/'+species+'/EDGAR_v4.3.2/'+species+'_sector_yearly/')
    
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
    for f in glob.glob(edpath+'v432_'+species+'_*'):
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
    speciesmm = molar_mass(species)
#    if species == 'CH4':
#        #speciesmm = 16.0425
#        speciesmm = molar_mass(species)
#    elif species == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % species
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
    
    
    #Read in EDGAR data of annual mean CH4 emissions for each sector
    #These are summed together
    #units are in kg/m2/s
    tot = None
    for sec in edgar_sectors:
        edgar = edpath+'v432_'+species+'_'+str(year)+'_IPCC_'+secdict[sec]+'.0.1x0.1.nc'    
        if os.path.isfile(edgar):
            ds = xr.open_dataset(edgar)
            soiname = 'emi_'+species.lower()
            if tot is None:
                tot = old_div(ds[soiname].values*1e3,speciesmm)
            else:
                tot += old_div(ds[soiname].values*1e3,speciesmm)
        else:
            print('No annual file for sector %s and %s' % (sec, species))
        
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon))    
    narr, reg = regrid2d(tot, lat_in, lon_in,
                             lat_out, lon_out)
    
    return(narr)   

def getedgarmonthlysectors(lon_out, lat_out, edgar_sectors, months=[1,2,3,4,5,6,7,8,9,10,11,12],
                           species='CH4'):
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
        species (str, optional):
            Which species you want to look at. 
            e.g. species = 'CH4'
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
    species = species.upper() #Make sure species is uppercase
    #Path to EDGAR files
    edpath = os.path.join(data_path,'Gridded_fluxes/'+species+'/EDGAR_v4.3.2/'+species+'_sector_monthly/')
    
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
    speciesmm = molar_mass(species)
#    if species == 'CH4':
#        speciesmm = 16.0425
#    elif species == 'N2O':
#        speciesmm = 44.013
#    else:
#        print "No molar mass for species %s." % species
#        print "Please add this and rerun the script"
#        print "Returning None"
#        return(None)
    
    
    #Read in EDGAR data of annual mean CH4 emissions for each sector
    #These are summed together
    #units are in kg/m2/s
    warnings = []
    first = 0
    for month in months:
        tot = np.array(None)
        for sec in edgar_sectors:
            edgar = edpath+'v432_'+species+'_2010_'+str(month)+'_IPCC_'+secdict[sec]+'.0.1x0.1.nc'    
            if os.path.isfile(edgar):
                ds = xr.open_dataset(edgar)
                soiname = 'emi_'+species.lower()
                if tot.any() == None:
                    tot = old_div(ds[soiname].values*1e3,speciesmm)
                else:
                    tot += old_div(ds[soiname].values*1e3,speciesmm)
            else:
                warnings.append('No monthly file for sector %s' % sec)
                #print 'No monthly file for sector %s' % sec
        
            if first == 0:
                emissions = np.zeros((len(months), tot.shape[0], tot.shape[1]))
                emissions[0,:,:] = tot
            else:
                first += 1
                emissions[first,:,:] = tot
            
    for warning in np.unique(warnings):
        print(warning)
                           
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    
    narr = np.zeros((nlat, nlon, len(months)))
       
    for i in range(len(months)):
       narr[:,:,i], reg = regrid2d(emissions[i,:,:], lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)

def getScarpelliFossilFuelsCH4(lon_out, lat_out, scarpelli_sector='all'):
    """
    NOTE THAT THIS IS NOT YET PUBLISHED. SPEAK TO ME (LUKE) BEFORE USING IT.
    The inventory is currently for 2012 methane emissions from fuel 
    exploitation (including oil/gas/coal) - this pertains to fugitive emissions 
    from these activities (so not combustion emissions). The inventory was 
    built by downscaling national emissions reported by countries to the UN to 
    the grid scale (0.1 degree resolution). Emissions were allocated to grid 
    scale based on infrastructure locations from various sources, including the 
    Global Oil and Gas Infrastructure Inventory and Geodatabase (includes 
    refineries, pipelines, storage, stations, etc.) and wells locations from 
    DrillingInfo.
    
    Args:
        lon_out (array): 
            Longitudes to output the data on
        lat_out (array):
            Latitudes to output the data on
        sector (string):
            Source of emissions. Options are: 'coal', 'oil', 'gas' or 'all'.
            Default = 'all'
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon]
    """
    
    print('NOTE THAT THE SCARPELLI DATASET IS NOT YET PUBLISHED. SPEAK TO ME (LUKE) BEFORE USING IT.')
    
    sectors = {'coal' : 'Coal', 'gas' : 'Gas_All', 'oil' : 'Oil_All', 'all' : 'Total_Fuel_Exploitation'}

    path = os.path.join(data_path,'Gridded_fluxes/CH4/Scarpelli_FossilFuel_CH4/')
    fnroot = 'Global_Fuel_Exploitation_Inventory_'
    sourcefn = path+fnroot+sectors[scarpelli_sector]+'.nc'
    
    with xr.open_dataset(sourcefn) as load:
        ds = load.load()         #Units are Mg / km2 / year
    
    ffemis = ds.emis_ch4.values/(365*24*3600)/molar_mass('ch4') #Convert to mol/m2/s
    lat_in = ds.lat.values
    lon_in = ds.lon.values
    
    narr, reg = regrid2d(ffemis, lat_in, lon_in,
                             lat_out, lon_out)
    return(narr)
    
def getYanetalRiceCH4(lon_out, lat_out):
    """
    Global rice emissions from Yan et al. 2009 (for year 2000).
    https://doi.org/10.1029/2008GB003299
    Data is broken down by month.
    
    Args:
        lon_out (array): 
            Longitudes to output the data on.
        lat_out (array):
            Latitudes to output the data on.
    
    Returns:
        narr (array): 
            Array of regridded emissions in mol/m2/s.
            Dimensions are [lat, lon, time]
    """    
    species = "ch4"
    #Base year is 2000   
    dim = np.array([31,29,31,30,31,30,31,31,30,31,30,31]) 
    speciesmm = molar_mass(species)
    months = [str(x).zfill(2) for x in np.arange(12)+1]
    nt = len(months)
    lat_in = np.arange(360)*0.5 - 89.75
    lon_in = np.arange(720)*0.5 + 180.25
    area = areagrid(lat_in,lon_in)
    emissions = np.zeros((len(lat_in),len(lon_in),nt))
    
    for m,month in enumerate(months):
        #Units Gg / 0.5x0.5 degree / month 
        datain = np.genfromtxt(data_path+'/Gridded_fluxes/CH4/YanetalRice/rice_CH4_half_dg_m'+month+'.txt', skip_header=2)
        data = datain.reshape((360,720))
        data[data < 0] = 0.
        emis = data*1e9/area/speciesmm/(dim[m]*24.*3600.)
        emissions[:,:,m] = emis
    
    nlat = len(lat_out)
    nlon = len(lon_out) 
    nt = 12
    narr = np.zeros((nlat, nlon, nt))   
    for i in range(nt):
       narr[:,:,i], reg = regrid2d(emissions[:,:,i], lat_in, lon_in,
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
    The _SWAMPSfile function opens the correct SWAMPS wetland fraction file as an 
    xarray.Dataset object.
    '''
    path = os.path.join(data_path,"Gridded_fluxes/CH4/SWAMPS")
    #filename_swamps = os.path.join(path,"fw_swamps-glwd_2000-2012.nc") # Previous version
    #filename_swamps = os.path.join(path,"gcp-ch4_wetlands_2000-2017_05deg.nc") # Without inland water
    filename_swamps = os.path.join(path,"gcp-ch4_wetlands-and-inland-water_2000-2017_025deg_norm.nc")
    
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

def _wetlandArea(frac,lon_in,lat_in,lon_out,lat_out):
    '''
    The _wetlandArea function calculates the area of the wetland extent within the specified
    latitude and longitude grid and the fraction of the total global wetlands area.
    
    Args:
        frac (np.array) :
            nlat x nlon array of fractional wetland extent.
        lat_in, lon_in (np.array, np.array) :
            Latitude and longitude grid of input fractional extent
        lat_out, lon_out (np.array, np.array)
            Output Latitude and longitude grid
    
    Returns:
        tuple (float,float) :
            wetland area in m2 within specified latitude and longitude output grid, 
            fraction of global wetland area.
    '''
    grid_wetl = areagrid(lat_in,lon_in)
    wetl_area = np.sum(frac*grid_wetl)
    
    frac_wetl_domain = regrid2d(frac, lat_in, lon_in, lat_out, lon_out)[0]
    grid_domain = areagrid(lat_out,lon_out)
    wetl_area_domain = np.sum(frac_wetl_domain*grid_domain)
    
    return wetl_area_domain,old_div(wetl_area_domain,wetl_area)

def _SWAMPSwetlandArea(year,lon_out,lat_out,month=1):
    '''
    The _SWAMPSwetlandArea function calculates the area of the wetland extent 
    within a grid and the fraction of the total global wetlands area based on 
    wetlands fraction input for SWAMPS across a global grid.
    
    Args:
        year (int) :
            Year of interest.
        lat_out (numpy.array) :
            Latitude grid.
        lon_out (numpy.array) :
            Longitude grid.
        month (int, optional) :
            Month to extract this area for. Should be between 1 and 12.
            Default = 1 (i.e. January)
    
    Returns:
        tuple (float,float) :
            wetland area in m2 within specified latitude and longitude grid, 
            fraction of global wetland area.
    '''
    if month >= 1 and month <= 12:
        month_id = month-1
    else:
        raise ValueError("Did not recognise month input: {}. Expect value between 1 and 12.".format(month))
    
    frac_swamps = xr.open_dataset(_SWAMPSfile())
    fw = "Fw" # Name of variable within file
    
    lat_in = frac_swamps.lat.values
    lon_in = frac_swamps.lon.values

    area_swamps = frac_swamps[fw][month_id,:,:].values
    area_swamps = np.nan_to_num(area_swamps)

    wetl_area_domain,wetl_frac = _wetlandArea(area_swamps,lon_in,lat_in,lon_out,lat_out)

    return wetl_area_domain,wetl_frac

    
def getJULESwetlands(year,lon_out,lat_out,species="CH4",extent_db="SWAMPS",scale_wetlands=True,
                     total_w_emission=185e12):
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
        species (str, optional) :
            Species of interest. At the moment, this should be "CH4".
        extent_db (str, optional) :
            Which database to use for the wetlands extent while using JULES for the emissions.
            Default = "SWAMPS"
        scale_wetlands (bool, optional) :
            Whether to scale emissions by the wetland fraction of a total emissions value.
            Default = True.
        total_w_emission (float, optional) :
            If scale_emissions=True, this is the global emissions assumed for wetlands methane emissions.
            The wetlands fraction within the output area will then be used to find the emissions based on
            total emissions.
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
            Default = 185e12 (value for total global wetlands emissions in 2012 from Saunois et al, 2016)
    
    Returns:
        numpy.array :
            Re-gridded emissions map with dimensions (lat,lon,time)
    '''
    
    if species.upper() != "CH4":
        print("Unable to extract JULES wetland values for any species except CH4 (at the moment)")
        return None
    
    flux_jules = _readJULESwetland(year,species)
    
    fch4_name = "fch4_wetl_npp"
    fwetl_name = "fwetl"
    
    # -1.0e30 used as fill value for JULES data - essentially zero pixels/nan
    # want to set these to 0.0
    fill_value = np.min(flux_jules[fch4_name].values) # This may incorrect if input file is changed and different fill value is specified
    fch4_fill_indices = np.where(flux_jules[fch4_name]==fill_value)
    
    #flux_jules.fch4_wetl_npp.values[fch4_fill_indices] = 0.0
    flux_jules[fch4_name].values[fch4_fill_indices] = 0.0

    ## Want to convert from units of kg/m2/s to mol/m2/s
    molmass = molar_mass(species)
    flux_jules[fch4_name].values = flux_jules[fch4_name].values*1000./molmass

    flux_jules_frac = np.abs(old_div(flux_jules[fch4_name], flux_jules[fwetl_name]))
    flux_jules_frac.values = np.nan_to_num(flux_jules_frac) # Any number divided by 0.0 will be nan, so change these back to 0.0

    if extent_db == "SWAMPS":
        ## Multiply by fractions from SWAMPS to rescale to measured rather than simulated inundation area
        frac_swamps = xr.open_dataset(_SWAMPSfile())
        fw_name = "Fw"
        
        frac_swamps[fw_name].values = np.nan_to_num(frac_swamps[fw_name])
        frac_reindex = frac_swamps.reindex_like(flux_jules_frac,method="ffill")
    elif extent_db == None:
        frac_reindex = None
    else:
        raise Exception("Input for wetland extent database not understood: {}".format(extent_db))
    
    if frac_reindex:
        fch4_wetl_npp_new = flux_jules_frac*frac_reindex[fw_name]
    else:
        fch4_wetl_npp_new = flux_jules[fch4_name]

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
            if extent_db == "SWAMPS":
                frac = _SWAMPSwetlandArea(year,lon_out,lat_out,month=i+1)[1]
            scale = round(old_div((total_w_emission*frac),1e12),1)*1e12
            print("{:02}) Scaling total JULES wetlands emissions within domain to: {} g/yr".format(i+1,scale))
            narr[:,:,i] = scale_emissions(narr[:,:,i],species,lat_out,lon_out,total_emission=scale)

    return narr


def scale_emissions(emissions_t,species,lat,lon,total_emission,return_scaling=False):
    '''
    The scale_emissions function scales emissions values for one time grid to a 
    total_emission value.
    
    Args:
        emissions_t (numpy.array) :
            Emissions values for one time point. Expect array to have dimensions nlat x nlon
            Assumes emissions are in units of mol/m2/s
        species (str) :
            Species of interest. Used to extract molar mass. e.g. "CH4"
        lat (numpy.array) :
            Latitude grid for emissions.
        lon (numpy.array) :
            Longitude grid for emissions.
        total_emission (float) :
            Total emissions of the output array. 
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
        return_scaliing (bool, optional) :
            Return the scaling factor associated with updating the emissions field.
            Default = False
    
    Returns:
        numpy.array :
            Scaled emissions array (nlat x nlon x nt)
    '''

    # Calculate number of moles
    molmass = molar_mass(species)
    
    total_time = 365.*3600.*24. # seconds in a year

    areas = areagrid(lat,lon)
    
    current_emission = np.sum(emissions_t*areas)*total_time*molmass
    scaling = total_emission/current_emission
    print("Current emissions total: {} g/yr (scaling needed to {} g/yr: {})".format(current_emission,total_emission,scaling))
    
    emissions_new = np.copy(emissions_t)*scaling
    #print "New total emissions: {} g/yr".format(np.sum(emissions_new*areas)*total_time*molmass)
    
    if return_scaling:
        return emissions_new,scaling
    else:
        return emissions_new

def scale_emissions_all(emissions,species,lat,lon,total_emission,return_scaling=False):
    '''
    The scale_emissions_all function scales emissions values within all time grids to the same 
    total_emission value.
    
    Args:
        emissions (numpy.array) :
            Emissions values across multiple times. Expects array to have dimensions nlat x nlon x nt
        species (str) :
            Species of interest. Used to extract molar mass. e.g. "CH4"
        lat (numpy.array) :
            Latitude grid for emissions.
        lon (numpy.array) :
            Longitude grid for emissions.
        total_emission (float) :
            Total emissions of the output array. 
            Value should be specified in g/yr e.g. 185Tg/yr should be 185e12.
        return_scaliing (bool, optional) :
            Return the scaling factor associated with updating the emissions field.
            Default = False
    
    Returns:
        numpy.array :
            Scaled emissions array (nlat x nlon x nt)
    '''

    nt = emissions.shape[2]
    
    emissions_new = np.copy(emissions)
    
    scaling = np.zeros(nt)
    for i in range(nt):
        emissions_i = emissions_new[:,:,i]
        #scale_emissions(emissions_i,species,lat,lon,total_emission)
        if return_scaling:
            emissions_new[:,:,i],scaling[i] = scale_emissions(emissions_i,species,lat,lon,total_emission,return_scaling)
        else:
            emissions_new[:,:,i] = scale_emissions(emissions_i,species,lat,lon,total_emission,return_scaling)
    
    if return_scaling:
        return emissions_new,scaling
    else:
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

def _define_time(year,timeframe=None,periods=None,months=[]):
    '''
    The _define_time function creates an array of datetimes across a year based on the 
    options provided.
    
    The following options and combinations can be specified:
        timeframe - description of timeframe e.g. "yearly"
        periods - number of time periods within a year
        periods & timeframe - number of time periods for that timeframe
        months - specific months within the year
        months & timeframe - specific months over that timeframe
    
    Args:
        year (int) :
            Year of interest.
        timeframe (str, optional) :
            Timeframe to output.
            Options are: "yearly", "monthly", "daily" or "3hourly"
        periods (int, optional) :
            Number of periods within a year to include.
            If no timeframe is specified, expected periods are based on being split 
            evenly across a year based on the 4 possible timeframes listed above.
        months (list, optional) :
            Specific months to include within the year. Should be numbers between 1-12.
    
    Returns:
        np.array :
            Array of datetime objects
    '''
    freq_dict = {"yearly":"AS","monthly":"MS","daily":"D","3hourly":"3H"}
    
    if (timeframe not in list(freq_dict.keys())) and (not periods) and (not months):
        raise KeyError("Did not recognise input '{}'. Should be one of: {}".format(timeframe,list(freq_dict.keys())))
    elif periods and not timeframe:
        if periods == 1:
            timeframe="yearly"
        elif periods == 12:
            timeframe="monthly"
        elif periods == 365 or periods == 366:
            timeframe="daily"
        elif periods == 365*(24/3.) or periods == 366*(24/3.):
            timeframe == "3hourly"
    elif not periods and not months and not timeframe:
        raise ValueError("One of timeframe, periods or months must be specified.")
    
    if periods:
        start = "{}-01-01".format(year)
        #datetimeindex = pd.DatetimeIndex(start=start,periods=periods,freq=freq_dict[timeframe],closed="left")
        datetimeindex = pd.date_range(start=start,periods=periods,freq=freq_dict[timeframe],closed="left")
    elif months and timeframe:
        start_datetime = ["{}-{:02}-01".format(year,int(month)) for month in months]
        end_datetime = ["{}-{:02}-01".format(year,int(month)+1) if month != 12 else "{}-{:02}-01".format(year+1,1) for month in months]
        for i,start,end in zip(list(range(len(start_datetime))),start_datetime,end_datetime):
            if i == 0:   
                #datetimeindex = pd.DatetimeIndex(start=start,end=end,freq=freq_dict[timeframe],closed="left")
                datetimeindex = pd.date_range(start=start,end=end,freq=freq_dict[timeframe],closed="left")
            else:
                #datetimeindex = datetimeindex.append(pd.DatetimeIndex(start=start,end=end,freq=freq_dict[timeframe],closed="left"))
                datetimeindex = datetimeindex.append(pd.date_range(start=start,end=end,freq=freq_dict[timeframe],closed="left"))
    elif months:
        datetime = ["{}-{:02}-01".format(year,int(month)) for month in months]
        datetimeindex = pd.to_datetime(datetime,format="%Y-%m-%d")
    else:
        start = "{}-01-01".format(year)
        end = "{}-01-01".format(year+1)
        #datetimeindex = pd.DatetimeIndex(start=start,end=end,freq=freq_dict[timeframe],closed="left")
        datetimeindex = pd.date_range(start=start,end=end,freq=freq_dict[timeframe],closed="left")
    
    time = np.array(datetimeindex)
    #time = time.astype("datetime64[ns]")
    
    return time

def _define_prior_dict(databases):
    '''
    The _define_prior_dict function creates a dictionary of details for the specified priors.
    This can be passed to the flux.write() function to include prior details within the emissions
    file attributes.
    
    Args:
        databases (list) :
            List of database names. Full listed of accepted names within create_emissions() function.
    
    Returns:
        dict (list) :
            Contains a three-item list of [version,resolution,reference] for each database.
    '''
    # Note: "natural" has global resolution of 4 x 5 degrees in Fung et al 1987 but 1 x 1 degrees in file.
    
    resolution = {"EDGAR":"0.1 x 0.1 degrees",
                  "GFED":"0.25 degrees x 0.25 degrees",
                  "JULES_wetlands":"0.5 x 0.5 degrees",
                  "SWAMPS":"0.5 x 0.5 degrees",
                  "natural":"1.0 x 1.0 degrees",
                  "soilsink":"1.0 x 1.0 degrees",
                  "Bloom":"3.0 x 3.0 degrees",
                  "NAEI_and_EDGAR":"0.234 x 0.352 degrees",
                  "NAEI":"0.009 x 0.014 degrees (1.0 x 1.0 km)",
                  "Scarpelli":"0.1 x 0.1 degrees",
                  "Bloom2017":"0.5 x 0.5 degrees"}

    prior_info = {"EDGAR":["v4.3.2",resolution["EDGAR"],"http://edgar.jrc.ec.europa.eu/overview.php?v=432&SECURE=123"],
                  "GFED":["v4.1",resolution["GFED"],"https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1293"],
                  "JULES_wetlands":["v4.1",resolution["JULES_wetlands"],"Created by: nicola.gedney@metoffice.gov.uk"],
                  "SWAMPS":["Global Carbon Project CH4 v2",resolution["SWAMPS"],"Schroeder et al. 2015\nBased on v3.2 SWAMPS\nCreated by: benjamin.poulter@nasa.gov"],
                  "natural":["v1",resolution["natural"],"Fung et al. 1987"],
                  "soilsink":["v1",resolution["soilsink"],"Bousquet et al. 2006"],
                  "Bloom":["v1",resolution["Bloom"],"Bloom et al. 2012"],
                  "NAEI_and_EDGAR":["v1",resolution["NAEI_and_EDGAR"],"Provided by: UK Met Office"],
                  "NAEI":["unknown",resolution["NAEI"],"http://naei.beis.gov.uk/data/"],
                  "Scarpelli":["v0",resolution["Scarpelli"],"Tia Scarpelli (personal communication - Luke Western)"],
                  "Bloom2017":["v1",resolution["Bloom2017"], "Bloom et al. 2017 https://doi.org/10.5194/gmd-10-2141-2017"]}

    prior_info_dict = {}
    
    for db in databases:
        prior_info_dict[db] = prior_info[db]
        if db == "JULES_wetlands":
            prior_info_dict["SWAMPS"] = prior_info["SWAMPS"]
    
    return prior_info_dict

def database_options(print_options=False):
    '''
    Defines emissions database options and associated functions.
    Set print_options=True to print out the available databases to use (and which functions
    they are associated with).
    
    Returns:
        (dict, dict, dict, dict, list):
            - Dictionary of database names and associated functions. Note: can include nested dictionaries if more than one function can be used for a database e.g. EDGAR.
            - Dictionary of species for each database.
            - Dictionary of any timeframes for each database.
            - Dictionary of sector for each database.
            - List of databases which contain climatological values rather than specific years.
    '''
    db_functions = {"GFED":getGFED,
                    "EDGAR":{"yearly":getedgarannualtotals,
                             "sector_yearly":getedgarannualsectors,
                             "sector_monthly":getedgarmonthlysectors},
                    "EDGAR_yearly":getedgarannualtotals,
                    "EDGAR_sector_yearly":getedgarannualsectors,
                    "EDGAR_sector_monthly":getedgarmonthlysectors,
                    "natural":getothernaturalCH4,
                    "soilsink":getsoilsinkCH4,
                    "Bloom":getbloomwetlandsCH4,
                    "NAEI":getNAEI,
                    "NAEI_and_EDGAR":getnaeiandedgarCH4,
                    "JULES_wetlands":getJULESwetlands,
                    "Scarpelli":getScarpelliFossilFuelsCH4,
                    "Bloom2017":getBloom2017}

    # Note GFED species not defined here yet.
    db_species = {"EDGAR":["CH4","N2O"],
                  "natural":["CH4"],
                  "soilsink":["CH4"],
                  "Bloom":["CH4"],
                  "NAEI_and_EDGAR":["CH4","N2O"],
                  "NAEI":["CH4","N2O"],
                  "JULES_wetlands":["CH4"],
                  "Scarpelli":["CH4"],
                  "Bloom2017":["CH4"]}
    
    db_timeframes = {"GFED":["monthly","daily","3hourly"],"Bloom":["daily","monthly"]}
    
    db_sector = {"GFED":"fire",
                 "EDGAR":"anthro",
                 "natural":"natural",
                 "soilsink":"soilsink",
                 "Bloom":"wetlands",
                 "NAEI_and_EDGAR":"anthro",
                 "NAEI":"anthro",
                 "JULES_wetlands":"wetlands",
                 "Scarpelli":"fossilfuels",
                 "Bloom2017":"wetlands"}
    
    db_climatology = ["natural","soilsink","Scarpelli"]

    if print_options:
        for db,fn in list(db_functions.items()):
            if not isinstance(fn,dict):
                print('"{}"'.format(db))
                print("    - {}() function".format(fn.__name__))
            else:
                print('"{}"  (uses one of multiple functions depending on inputs)'.format(db))
                for identifier,option in list(fn.items()):
                    print("   - {}() function ({} data)".format(option.__name__,identifier))
        #return None
    
    return db_functions,db_species,db_timeframes,db_sector,db_climatology

def create_emissions(databases,species,domain,year=None,lon_out=[],lat_out=[],
                     scale=1.0,write=True,output_directory=output_directory,
                     **kwargs):
    '''
    The create_emissions function can create a combined emissions file based on the list of databases
    provided and options specified.
    
    See individual functions for more details on databases and inputs.
    
    Either a domain or longitude and latitude grid values must be specified.
    
    Note: For "*EDGAR*", three databases are available based on annual data, annual sector 
    data and monthly sector data (for 2010 only). If "EDGAR" database is specified, this 
    function will attempt to discern which of these functions to use based on the inputs.
    It is recommended that you use the "*EDGAR_yearly*", "*EDGAR_sector_yearly*" or 
    "*EDGAR_sector_monthly*" keywords to ensure the expected database is used.
    
    WARNING: At the moment this function is unable to interpret "months" as an input (included
    as part of the keyword arguments) and so this input will be ignored if specified.
    This would be applicable to "GFED" and "EDGAR_sector_monthly" databases.
        
    Args:
        databases (list) :
            List of databases to use to create emissions file.
            The following inputs can be included:
                - "GFED"                 - GFED v4.1 biomass burning database
                - "EDGAR"                - EDGAR v4.3.2 anthropogenic database
                - "EDGAR_yearly"         - explicitly use yearly EDGAR database
                - "EDGAR_sector_yearly"  - explicitly use yearly sector EDGAR database (sector inputs must be included)
                - "EDGAR_sector_monthly" - explicitly use monthly sector EDGAR database (sector inputs must be included)
                - "natural"              - other natural CH4 emissions (volcanoes, termites, hydrates) from Fung et al 1987
                - "soilsink"             - soil sink CH4 emissions from Bousquet et al 2006
                - "Bloom"                - CH4 wetland emissions from Bloom et al
                - "Bloom2017"            - CH4 wetland emissions for updated WetCHARTS map from Bloom et al 2017
                - "NAEI"                 - NAEI anthropogenic database
                - "NAEI_and_EDGAR"       - Combined dataset of NAEI and EDGAR
                - "JULES_wetlands"       - CH4 wetlands emissions from combining JULES model emissions output and SWAMPS wetlands extent.
                - "Scarpelli"            - CH4 from fugitive fossil fuel emissions (NB this is not yet published)
            See database_options() function for full and correct list of options.
        species (str) :
            Species of interest. All listed databases must have data for this species.
        domain (str) :
            Name of domain e.g. "EUROPE","SOUTHAMERICA","NORTHAFRICA".
            See $DATA_PATH/LPDM/fp_NAME/* for all available domains.
            If lat_out and lon_out not specified, domain will be used to extract grid values.
        year (int) :
            Year of emissions to use.
            Unless all input databases are climatological, year must be specified.
        lon_out (list, optional) :
            Longitude grid for emissions output. Only needed if files for the domain 
            have not already been created.
        lat_out (list, optional) :
            Latitude grid for emissions output. Only needed if files for the domain 
            have not already been created.
        scale (float, optional) :
            Scale total emissions by some factor.
            By default scaling factor is 1.0, meaning no scaling is applied.
        write (bool, optional) :
            Write to netCDF file. Passes to flux.write() function to ensure correct format
            is applied.
            Default = True.
        output_directory (str, optional) :
            Base directory for output file. Domain name will be used as subdirectory name.
            Default = $DATA_PATH/LPDM/emissions/

        
        kwargs :
            Additional keyword arguments which can be taken by the underlying functions.
            See "get*" functions within this module for full list of arguments.
            
            timeframe :
                Can be one of "3hourly","daily","monthly". Optional for "GFED" and "Bloom" databases.
            edgar_sectors :
                Specific sectors for EDGAR database. MUST be specified if "EDGAR_sector_yearly"
                or "EDGAR_sector_monthly" databases are used.
            naei_sector :
                Specific sector for NAEI database. MUST be specified if "NAEI" database is used.
            scarpelli_sector :
                Specific sector for Scarpelli inventory. MUST be specified if "Scarpelli" inventory is used.
            incagr :
                Include agricultural waste burning. Optional for "GFED" database.
            extent_db :
                Use a different extent database with the JULES emissions. Optional for "JULES_wetlands"
            scale_wetlands :
                Scale wetlands to match a fraction of a total emissions value. Optional for "JULES_wetlands"
            total_w_emission :
                Total wetland emissions to scale emissions map to. Only used if scale_wetlands=True. Optional for "JULES_wetlands"
            modeltype :
                Model type to use for WetCHARTS based on options available. One of 'extended' or 'full'.
                Optional for "Bloom2017"

    Returns:
        numpy.array (nlat x nlon x nt), numpy.array (nt):
            Array of combined emissions from all databases.
            Array of time values.
        
        If write is True, writes to netcdf file.
        Output file naming convention defined within flux.write (check for final format).
    '''
    
    db_functions,db_species,db_timeframes,db_sector,db_climatology = database_options()
    
    EDGAR_options = ["EDGAR_yearly","EDGAR_sector_yearly","EDGAR_sector_monthly"]
    timeframe_options = {"monthly":12,"daily":365,"3hourly":old_div(365*24,3)}
    
    #if "EDGAR" in databases and "GFED" in databases and "sectors" not in kwargs:
    #    raise Exception("Cannot combine EDGAR and GFED ")
    
    if not os.path.exists(output_directory):
        raise IOError("Specified output directory does not exist: {}".format(output_directory))
    
    if not np.any(lat_out) and not np.any(lon_out):
        if domain:
            lat_out,lon_out,height = domain_volume(domain)
        else:
            lat_out = []
            lon_out = []
    
    if "timeframe" in kwargs:
        timeframe = kwargs["timeframe"]
        print("When applicable, using timeframe: {}".format(timeframe))
    else:
        timeframe = None
    
    if "months" in kwargs:
        ## TODO: Currently unable to work with a sub-set of months. Need to decide on sensible behaviour
        print("WARNING: Unable to create emissions using subset of months specified by: {} at the moment".format(kwargs["months"]))
        kwargs.pop("months")
    
    kwargs["species"] = species
    kwargs["lon_out"] = lon_out
    kwargs["lat_out"] = lat_out
    
    # Define list of functions to call based on inputs databases
    functions = []
    for i,database in enumerate(databases):
        if database == "EDGAR":
            # EDGAR has three functions associated functions, use inputs to work out which one to use.
            if timeframe:
                if timeframe == "yearly" and "sectors" not in kwargs:
                    functions.append(db_functions[database]["yearly"])
                elif timeframe == "yearly" and "sectors" in kwargs:
                    functions.append(db_functions[database]["sector_yearly"])
                elif timeframe == "monthly" and "sectors" in kwargs and "months" in kwargs:
                    functions.append(db_functions[database]["sector_monthly"])
                elif timeframe == "monthly" and "sectors" in kwargs and not "months" in kwargs:
                    print("Sectors specified and monthly timeframe specified but no months specified. Using EDGAR annual sectors as default.")
                    functions.append(db_functions[database]["sector_yearly"])
                elif timeframe != "monthly" and timeframe != "yearly" and "sectors" in kwargs:
                    print("Sectors specified but timeframe of monthly or yearly is not specified. Using EDGAR annual sectors as default.")
                    functions.append(db_functions[database]["sector_yearly"])
                elif timeframe == "monthly" and "sectors" not in kwargs:
                    print("Only able to extract monthly EDGAR data when sectors are also specified. Using EDGAR annual totals.")
                    functions.append(db_functions[database]["yearly"])
                elif timeframe != "yearly" and "sectors" not in kwargs:
                    print("Using EDGAR annual totals as default.")
                    functions.append(db_functions[database]["yearly"])
                else:
                    raise Exception("Did not recognise combined input for EDGAR database: {}".format(kwargs))
            elif "edgar_sectors" in kwargs:
                print("Sectors specified but timeframe of monthly or yearly is not specified. Using EDGAR annual sectors as default.")
                functions.append(db_functions[database]["sector_yearly"])
            else:
                print("No timeframe specified. Using EDGAR annual totals as default.")
                functions.append(db_functions[database]["yearly"])
        else:
            functions.append(db_functions[database])
        
        if database in EDGAR_options:
            databases = databases[:]
            databases[i] = "EDGAR"
   
    # Checks species can be resolved for all databases in list.
    for database in databases:
        if database in list(db_species.keys()):
            if species.upper() not in db_species[database]:
                raise Exception("Cannot create emissions map including '{}' database for species {}".format(database,species))
    
    # Determine if requested databases are for climatology
    climatology_db = [database for database in databases if database in db_climatology]
    if climatology_db:
        # If all requested databases contain climatological data, set the flag to True
        if databases == climatology_db:
            print("All specified databases are for climatological data rather than a specific year.") 
            print("Output will use default year for climatology data: 1900")
            climatology = True
        else:
            climatology = False
    else:
        climatology = False

    # Check if year has been defined and, if not, only continue if all databases are for climatology
    if not climatology:
        if year != None:
            kwargs["year"] = year
        else:
            raise Exception("For these input databases, year must be specified.")

    emissions_list = []
    for fn,db in zip(functions,databases):    
        
        # Extract list of parameter inputs for each function
        all_param = fn.__code__.co_varnames[:fn.__code__.co_argcount]

        # Create dictionary of relevant parameters for function from inputs
        param = {}
        for p in all_param:
            if p in kwargs:
                # For timeframe input need to check against timeframe options and choose most appropriate
                if p == "timeframe":
                    if kwargs[p] in db_timeframes[db]:
                        param[p] = kwargs[p]
                    else:
                        # If timeframe does not match any available values
                        # find the available timeframe which has closest frequency.
                        request_tf = kwargs[p]
                        nperiod = timeframe_options[request_tf]
                        
                        avail_tf = db_timeframes[db]
                        diff = [nperiod-timeframe_options[tf] for tf in avail_tf]
                        index = diff.index(min(i for i in diff if i>=0))
                        actual_tf = avail_tf[index]
                        
                        param[p] = actual_tf
                # Otherwise the parameter can just be taken from the input keyword arguments
                else:
                    param[p] = kwargs[p]
        
        print("--------------------------------")
        print("Calling function: {}(...)".format(fn.__name__))
        print("Calling with inputs: {} (not specified {}, using defaults)\n".format(list(param.keys()),[k for k in all_param if k not in param]))
        
        narr = fn(**param)
        if len(narr.shape) == 2:
            narr = np.expand_dims(narr,2)
        ntime = narr.shape[2]
        
        if climatology:
            time = _define_time(1900, periods=ntime)
        else:
            time = _define_time(year, periods=ntime)
        
        ds = xray.Dataset({"flux":(("lat","lon","time"),narr)},coords={"lat":lat_out,"lon":lon_out,"time":time})
        emissions_list.append(ds)
    
    # Re-order based on time frequency (most frequent first)
    time_freq = [emis.time.size for emis in emissions_list]
    emissions_list = [emis for tf,emis in sorted(zip(time_freq,emissions_list),reverse=True)]
    
    emissions_init = emissions_list[0] # Extract first emissions dataset
    
    # Match time grid of remaining datasets to the first and stack emissions
    for i,ds in enumerate(emissions_list[1:]):
        emissions_match = ds.reindex_like(emissions_init,method="ffill")
        emissions_init["flux"].values += emissions_match["flux"].values
    
    emissions = emissions_init["flux"].values
    time = emissions_init["time"].values
    
    if scale:
        emissions = emissions*scale
    
    if write:
        ## Constructing attribute inputs for emissions file
        if len(databases) == 1:
            db = databases[0]
            title = "{} emissions from {} database.".format(species,db)
        else:
            title = "Combined {} emissions between ".format(species)
            for i,database in enumerate(databases):
                if i == 0:
                    title += "{} ({})".format(database,db_sector[database])
                elif i < len(databases)-1:
                    title += ", {} ({})".format(database,db_sector[database])
                else:
                    title += " and {} ({}) databases.".format(database,db_sector[database])
    
        if len(databases) == 1:
            db = databases[0]
            source = db_sector[db]
        else:
            source = "all"
            
        prior_info_dict = _define_prior_dict(databases)
        
        flux_comments = ''
        if "edgar_sectors" in kwargs:
            flux_comments += "Sectors included for anthropogenic emissions (EDGAR): "
            for i,sector in enumerate(kwargs["edgar_sectors"]):
                if i == 0:
                    flux_comments += "{}".format(sector)
                elif i < len(kwargs["edgar_sectors"])-1:
                    flux_comments += ", {}".format(sector)
                else:
                    flux_comments += " and {}.\n".format(sector)
        if "naei_sector" in kwargs:
           flux_comments += "{} sector included for anthropogenic emissions (NAEI).\n".format(kwargs["naei_sector"])
        if climatology_db:
            flux_comments += "Based on climatology from database(s): {}.\n".format(databases)
        
        if not flux_comments:
            flux_comments = None
        
        print("Writing output to directory: {}".format(output_directory))
        
        flux.write(lat_out,lon_out,time,emissions,species,domain,source=source,title=title,
               prior_info_dict=prior_info_dict,flux_comments=flux_comments,climatology=climatology,
               regridder_used='acrg_grid.regrid.regrid_2D',output_directory=output_directory)
        
    return emissions,time
