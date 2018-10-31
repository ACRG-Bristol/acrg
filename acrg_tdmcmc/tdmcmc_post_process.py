# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 08:57:53 2015


Script to process Transdimensional MCMC output 

Includes:
Write netcdf and append netcdf to write output to nc file

plot_scaling - plot posterior scaling map

regions_histogram - plot histogram of number of regions

country_emissions - calculate emissions from given list of countries
                    Currently hard-wired for methane


@author: ml12574
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xray
import acrg_name as name
from acrg_grid import areagrid
from netCDF4 import Dataset
from acrg_time.convert import time2sec
import os
import acrg_agage as agage
import json
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
from matplotlib import ticker
import cartopy.crs as ccrs
from collections import OrderedDict
import datetime as dt
import getpass

acrg_path = os.getenv("ACRG_PATH")

def append_netcdf(flux_mean, flux_percentile, flux_it, country_mean, country_percentile,
                  k_mean, k_percentile,
                 lon, lat, time, country, percentile, nIt, experiment, outfile):
    #country_it,  
    """
    Function for appending flux and country output to file
    
    Only use this function once file has already been created
    Create file first using write_netcdf function below
    """         
      
    fm_name = "_".join(['flux_mean', experiment])
    fpc_name = "_".join(['flux_percentile', experiment])
    fit_name = "_".join(['flux_it', experiment])
    countrym_name = "_".join(['country_mean', experiment])
    countrypc_name = "_".join(['country_percentile', experiment])
    regionsm_name = "_".join(['regions_mean', experiment])
    regionspc_name = "_".join(['regions_percentile', experiment])
    #countryit_name = "_".join(['country_it', experiment])
    
    print "outfile",outfile
    ncF=Dataset(outfile, 'a')
    
    ncfluxm=ncF.createVariable(fm_name, 'f', ('lat', 'lon', 'time'))    
    ncfluxpc=ncF.createVariable(fpc_name, 'f', ('lat', 'lon', 'time', 'percentile'))    
    nccountrym=ncF.createVariable(countrym_name, 'f', ('country', 'time'))   
    nccountrypc=ncF.createVariable(countrypc_name, 'f', ('country', 'time', 'percentile'))
    ncregionsm=ncF.createVariable(regionsm_name, 'f', ('time'))   
    ncregionspc=ncF.createVariable(regionspc_name, 'f', ('time', 'percentile'))
    ncfluxit=ncF.createVariable(fit_name, 'f', ('lat', 'lon', 'time', 'nIt'))    
    #nccountryit=ncF.createVariable(countryit_name, 'f', ('country', 'nIt','time'))   
    
    ncfluxm[:, :, :]=flux_mean
    ncfluxm.units='mol/m2/s'
    ncfluxm.long_name='Posterior mean flux from each grid-cell'
    
    ncfluxpc[:, :, :,:]=flux_percentile
    ncfluxpc.units='mol/m2/s'
    ncfluxpc.long_name='Posterior mean flux from each grid-cell'
    
    ncfluxit[:, :, :,:]=flux_it
    ncfluxit.units='mol/m2/s'
    ncfluxit.long_name='Posterior it flux from each grid-cell'
    
    nccountrym[:, :]=country_mean
    nccountrym.units='Tg/yr'
    nccountrym.long_name='Posterior mean emissions from individual countries'
    
    nccountrypc[:, :, :]=country_percentile
    nccountrypc.units='Tg/yr'
    nccountrypc.long_name='Posterior percentile emissions from individaul countries'
    
    ncregionsm[:]=k_mean
    #ncregionsm.units=' '
    ncregionsm.long_name='Posterior mean number of regions'
    
    ncregionspc[:,:]=k_percentile
    #ncregionspc.units='
    ncregionspc.long_name='Posterior percentile number of regions'
    
    #nccountryit[:,:,:]=country_it
    #nccountryit.units='Tg/yr'
    #nccountryit.long_name='Posterior emissions distribution from individual countries'
    
    ncF.close()
    print "Appended " + experiment + " to " + outfile 

def write_netcdf(flux_mean, flux_percentile, flux_it, flux_prior, flux_ap_percentile,
                 country_mean, country_percentile, country_prior, country_ap_percentile,
                 country_index, k_mean, k_percentile,
                 lon, lat, time, country, percentile, nIt, experiment, outfile):
    #country_it,
    """
    Function for writing all flux and country output to file
    
    Only use this function for first set of variables
    Once file has been created, use append_netcdf function above
    """    
    
    time_seconds, time_reference = time2sec(time)
    
    fm_name = "_".join(['flux_mean', experiment])
    fpc_name = "_".join(['flux_percentile', experiment])
    fit_name = "_".join(['flux_it', experiment])
    countrym_name = "_".join(['country_mean', experiment])
    countrypc_name = "_".join(['country_percentile', experiment])
    regionsm_name = "_".join(['regions_mean', experiment])
    regionspc_name = "_".join(['regions_percentile', experiment])
    #countryit_name = "_".join(['country_it', experiment])
    
    #Write NetCDF file
    ncF=Dataset(outfile, 'w')
    ncF.createDimension('time', len(time))
    ncF.createDimension('lon', len(lon))
    ncF.createDimension('lat', len(lat))
    ncF.createDimension('country', len(country))
    ncF.createDimension('percentile', len(percentile))
    ncF.createDimension('nIt', nIt)
    
    nctime=ncF.createVariable('time', 'i', ('time',))
    nclon=ncF.createVariable('lon', 'f', ('lon',))
    nclat=ncF.createVariable('lat', 'f', ('lat',))
    nccountry=ncF.createVariable('country_name', str, ('country',))
    ncpercent=ncF.createVariable('percentile', 'i', ('percentile',))
    ncfluxm=ncF.createVariable(fm_name, 'f', ('lat', 'lon', 'time'))    
    ncfluxpc=ncF.createVariable(fpc_name, 'f', ('lat', 'lon', 'time', 'percentile'))    
    ncfluxit=ncF.createVariable(fit_name, 'f', ('lat', 'lon', 'time', 'nIt'))    
    nccountrym=ncF.createVariable(countrym_name, 'f', ('country', 'time'))   
    nccountrypc=ncF.createVariable(countrypc_name, 'f', ('country', 'time', 'percentile'))
    nccountryap=ncF.createVariable('country_prior', 'f', ('country', 'time'))
    nccountryappc=ncF.createVariable('country_prior_percentile', 'f', ('country', 'time', 'percentile'))
    ncfluxap=ncF.createVariable('flux_prior', 'f', ('lat', 'lon', 'time'))  
    ncfluxappc=ncF.createVariable('flux_prior_percentile', 'f', ('lat', 'lon', 'time', 'percentile'))  
    nccountryid=ncF.createVariable('country_index', 'i', ('lat', 'lon'))    
    ncregionsm=ncF.createVariable(regionsm_name, 'f', ('time'))   
    ncregionspc=ncF.createVariable(regionspc_name, 'f', ('time', 'percentile'))
    #nccountryit=ncF.createVariable(countryit_name, 'f', ('country', 'nIt','time'))   

    
    nctime[:]=time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + np.str(time_reference)
    nctime.calendar='gregorian'

    nclon[:]=lon
    nclon.units='Degrees east'
    
    nclat[:]=lat
    nclat.units='Degrees north'

    nccountry[:]=np.array(country)
    
    ncpercent[:]=np.array(percentile)
    
#    ncfluxprior[:, :, :]=flux_mean
#    ncfluxprior.units='mol/m2/s'
#    ncfluxprior.long_name='Mean flux from each grid-cell'    
    
    ncfluxm[:, :, :]=flux_mean
    ncfluxm.units='mol/m2/s'
    ncfluxm.long_name='Posterior mean flux from each grid-cell'
    
    ncfluxpc[:, :, :,:]=flux_percentile
    ncfluxpc.units='mol/m2/s'
    ncfluxpc.long_name='Posterior percentile flux from each grid-cell'
    
    ncfluxit[:, :, :,:]=flux_it
    ncfluxit.units='mol/m2/s'
    ncfluxit.long_name='Posterior it flux from each grid-cell'
    
    nccountrym[:, :]=country_mean
    nccountrym.units='Tg/yr'
    nccountrym.long_name='Posterior mean emissions from individual countries'
    
    nccountrypc[:, :, :]=country_percentile
    nccountrypc.units='Tg/yr'

    nccountrypc.long_name='Posterior percentile emissions from individual countries'
    
    ncfluxap[:,:,:]=flux_prior
    ncfluxap.units='mol/m2/s'
    ncfluxap.long_name='Prior mean flux from each grid-cell'
    
    ncfluxappc[:, :, :,:]=flux_ap_percentile
    ncfluxappc.units='mol/m2/s'
    ncfluxappc.long_name='Prior percentile flux from each grid-cell'
    
    nccountryap[:, :]=country_prior
    nccountryap.units='Tg/yr'
    nccountryap.long_name='Prior mean emissions from individual countries'
    
    nccountryappc[:, :, :]=country_ap_percentile
    nccountryappc.units='Tg/yr'
    nccountryappc.long_name='Prior percentile emissions from individual countries'
        
    nccountryid[:,:]=country_index
    #ncfluxap.units='mol/m2/s'
    nccountryid.long_name='Index corresponding to each country'
    
    ncregionsm[:]=k_mean
    #ncregionsm.units=' '
    ncregionsm.long_name='Posterior mean number of regions'
    
    ncregionspc[:,:]=k_percentile
    #ncregionspc.units='
    ncregionspc.long_name='Posterior percentile number of regions'
    
#    nccountryit[:,:,:]=country_it
#    nccountryit.units='Tg/yr'
#    nccountryit.long_name='Posterior emissions distribution from individual countries'
    
    ncF.close()
    print "Written " + experiment + " to " + outfile

def molar_mass(species):
    '''
    This function extracts the molar mass of a species from the acrg_species_info.json file.
    Returns:
        float : Molar mass of species
    '''
    species_info_file = os.path.join(acrg_path,"acrg_species_info.json")
    with open(species_info_file) as f:
            species_info=json.load(f)
    species_key = agage.synonyms(species, species_info)
    molmass = float(species_info[species_key]['mol_mass'])
    return molmass

def g2mol(value,species):
    ''' Convert a value in g to mol '''
    molmass = molar_mass(species)
    return value*molmass

def x_post_mean(ds):
    '''
    The x_post_mean function extracts the posterior values for all iterations and calculates 
    the mean. This is reshaped onto a latitude x longitude grid.
    
    Expects latitude and longitude coords within dataset to be named "lat" and "lon".
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
    
    Returns:
        numpy.array (nlat x nlon):
            Array containing mean posterior emissions values on a latitude x longitude grid.
    '''
    nlon = len(ds["lon"])
    nlat = len(ds["lat"])
    
    x_post_vit = ds.x_post_vit.values
    x_post_v_mean = np.mean(x_post_vit, axis=0)
    x_post_m = np.reshape(x_post_v_mean, (nlat,nlon))
    
    return x_post_m

def x_post_percentile(ds,percentiles=[5,16,50,84,95]):
    '''
    The x_post_percentile function extracts the posterior values for all iterations and calculates
    a set of percentiles. This is reshaped onto latitude x longitude grid.
    
    Expects latitude and longitude coords within dataset to be named "lat" and "lon".
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
    
    Returns:
        numpy.array (nlat x nlon x npercentiles):
            Array containing posterior values for each percentile on a latitude 
            x longitude grid.
    '''
    nlon = len(ds["lon"])
    nlat = len(ds["lat"])
    
    x_post_vit = ds.x_post_vit.values
    
    x_post_pl = np.zeros((nlat,nlon,len(percentiles)))
    for i,percent in enumerate(percentiles):
        x_post_pl[:,:,i] = np.reshape(np.percentile(x_post_vit, percent, axis=0),(nlat,nlon))
    
    return x_post_pl

def flux_iterations(ds):
    '''
    The flux_iterations function calculates the posterior flux values for each saved iteration 
    within the tdmcmc output.
    
    Expects latitude, longitude and nIt coords within dataset to be named "lat", "lon" and "nIt".
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
    Returns:
        numpy.array (nlat x nlon x nIt):
            Array containing posterior flux values for each saved iteration within tdmcmc.
    '''
    nIt = len(ds["nIt"])
    nlon = len(ds["lon"])
    nlat = len(ds["lat"])
    
    x_post_vit = ds.x_post_vit.values
    x_post_it = np.reshape(x_post_vit,(nIt,nlat,nlon))
    
    q_ap = ds.q_ap.values
    
    flux_it = x_post_it*q_ap
    flux_it = np.moveaxis(flux_it,0,2) # Rearranging axes to make nIt last dimension
    
    return flux_it

def flux_mean(ds):
    '''
    The flux_mean function calculates the mean flux based on the prior emissions and posterior 
    values.
    
    Expects latitude and longitude coords within dataset to be named "lat" and "lon".
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
    Returns:
        numpy.array (nlat x nlon):
                Array containing posterior flux values on a latitude x longitude grid.     
    '''
    q_ap = ds.q_ap.values
    x_post_m = x_post_mean(ds)
    
    flux_mean = x_post_m*q_ap
    
    return flux_mean

def flux_percentile(ds,percentiles=[5,16,50,84,95]):
    '''
    The flux_percentile function calculates the flux values for a set of percentiles.
    
    Expects latitude and longitude coords within dataset to be named "lat" and "lon".

    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
    
    Returns:
        numpy.array (nlat x nlon x npercentiles):
            Array containing posterior flux values for each percentile on a latitude x longitude
            grid.
    '''
    nlon = len(ds["lon"])
    nlat = len(ds["lat"])
    
    q_ap = ds.q_ap.values
    
    x_post_pl = x_post_percentile(ds,percentiles)
    flux_pl = np.zeros((nlat,nlon,len(percentiles)))
    
    for i in range(x_post_pl.shape[-1]):
        flux_pl[:,:,i] = x_post_pl[:,:,i]*q_ap
    
    return flux_pl

def flux_diff(ds):
    '''
    The flux_diff function calculates the flux difference between the prior and the posterior.
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
    Returns:
        numpy.array (nlat x nlon):
            Array containing flux difference values on a latitude x longitude grid.
    '''
    x_post_m = x_post_mean(ds)
    q_ap = ds.q_ap.values
    
    q_abs_diff = (x_post_m*q_ap-q_ap)
    q_abs_diff = q_abs_diff
    
    return q_abs_diff

def k_mean_percentile(ds,percentiles=[5,16,50,84,95]):
    '''
    The k_mean_percentile function calculates the mean number of regions (k) and a set of percentiles
    within the saved iterations.
    
    Expects nIt coord within dataset to be named "nIt".
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects data set to contain:
                k_it - number of regions for each iteration. Dimensions = nIt
    Returns:
        if percentiles are specified:
            float, np.array (npercentiles) :
                Mean and percentile array for number of regions within saved iterations.
        if percentiles are not specified:
            float :
                Mean number of regions within saved iterations
    '''
    k_it = ds.k_it.values
    k_mean = np.mean(k_it)
    
    if percentiles:
        k_percentile = np.zeros(len(percentiles))
        for i,percentile in enumerate(percentiles):
            k_percentile[i] = np.percentile(k_it, percentile, axis=0)
    else:
        k_percentile = np.array([])
    
    if percentiles:
        return k_mean,k_percentile
    else:
        return k_mean

def define_stations(ds,sites=[]):
    '''
    The define_stations function defines the latitude and longitude values for each site within
    a dataset.
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects dataset to contain:
                release_lons - Longitude values for each site. Dimension = len(sites)
                release_lats - Latitude values for each site. Dimension = len(sites)
                y_site       - Site identifier for each measurement. Dimension = nmeasure
        sites (list) :
            List of sites to look for within dataset.
    
    Returns:
        dict :
            Dictionary containing "site"_lat, "site"_lon values for each site.
    '''
    if not sites:
        sites = ds.sites.values
    
    stations={}
    for si, site in enumerate(sites):
        #if site in ds.y_site:
        stations[site+'_lon']=ds.release_lons[si].values
        stations[site+'_lat']=ds.release_lats[si].values
        if site not in ds.y_site:
            print "WARNING: Reference to site not found within dataset"
    stations['sites']=sites
    
    return stations

def unbiasedDivergingCmap(data, zero = 0, minValue = None, maxValue = None):
    """
    Calculate the normalisation of a diverging cbar around a given value
    Prevents bias due to asymetry in data affectin scale
    
    Args:
        data : numpy array of data to plot
        zero : the centre value of the cmap
        minValue : smallest value to use in calculation
        maxValue : largest value to use in calculation
    
    Returns:
        a normalization function to be fed into plot
    """
    
    if maxValue is None:
        maxValue = np.amax(data)
    if minValue is None:
        minValue = np.amin(data)
    maxRange = max(abs(maxValue-zero), abs(minValue-zero) )
    
    return Normalize(vmin = zero-maxRange, vmax = zero + maxRange, clip=True)

def plot_map(data,lon,lat, clevels=None, divergeCentre = None, cmap=plt.cm.RdBu_r, label=None,
                 smooth = False, out_filename=None, stations=None, fignum=None,
                 title=None):
    
    """
    Plot 2d scaling map of posterior x
    i.e. degree of scaling applied to prior emissions 
    
    Args:
        data     : 2D (lat,lon) array of whatever you want
        clevels  : Contour levels; defaults to np.arange(-2., 2.1, 0.1)
        cmap     : Colormap - defaults to Red Blue reverse
        label    : String - to appear at bottom
        smooth   : If true plot smooth contours, otherwise use pcolormesh
        stations : Default is None. If specified needs to be a dictionary containing:
                    sites: 'MHD', 'TAC'
                    MHD_lon: -9.02
                    MHHD_lat: 55.2
                    TAC_lon: etc...
    Returns:
        None
        
        If out_filename is None:
            Created plot is saved to file
        Else:
            Plot is displayed interactively
    """        
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=np.median(lon)))
    ax.set_extent([lon[0], lon[-1], lat[0], lat[-1]], crs=ccrs.PlateCarree())
    ax.coastlines()

    if clevels is None:
        print "Warning: using default contour levels. Include clevels keyword to change"
        clevels = np.arange(-2., 2.1, 0.1)
    else:
        clevels = np.arange(clevels[0], clevels[1], clevels[2])
        
    if smooth == True:
        if divergeCentre is None:
            plt.contourf(lon, lat, data, transform=ccrs.PlateCarree(), cmap = cmap, levels=clevels)
        else:
            norm = unbiasedDivergingCmap(data, zero=divergeCentre, minValue=clevels[0], maxValue=clevels[-1])
            plt.contourf(lon, lat, data, transform=ccrs.PlateCarree(), cmap = cmap, levels=clevels,
                    norm = norm)
        cb = plt.colorbar(orientation='horizontal', pad=0.05)
    else:
        lons, lats = np.meshgrid(lon,lat)
        if divergeCentre is None:
            norm = BoundaryNorm(clevels,ncolors=cmap.N,clip=True)
        else:
            norm = unbiasedDivergingCmap(data, zero=divergeCentre, minValue=clevels[0], maxValue=clevels[-1])
        cs = plt.pcolormesh(lons, lats, data,cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        cb = plt.colorbar(cs, orientation='horizontal', pad=0.05, extend='both')
                
    if label is not None:        
        cb.set_label(label,fontsize=14) 
    if title is not None:
        plt.title(title, fontsize=16) 
        
    if stations is not None:
       
        for si,site in enumerate(stations['sites']):
            ilon=stations[site+'_lon']
            ilat=stations[site+'_lat']
            
            plt.plot(ilon, ilat, color = 'black', marker = 'o', markersize=8,  transform=ccrs.PlateCarree())
               
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()                 
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()
    
    
def regions_histogram(k_it, out_filename=None, fignum=2):
    
    """
    Plot a histogram of the number of regions in trandimensional inversion
    
    Use this as a check for whether kmax should be larger
    """
    bin_max=np.ceil((np.max(k_it)+50.)/40.)*40
    #bin_max=100
    #bin_step = bin_max/40
    bin_step = 2
    bins=np.arange(0,bin_max,bin_step)    
    
    plt.figure(fignum)
    # the histogram of the data with histtype='step'
    n, bins, patches = plt.hist(k_it, bins, histtype='bar', color='mediumblue')
    #n, bins, patches = plt.hist(k_it, bins, normed=1, histtype='bar', rwidth=0.8)
    
    # add a 'best fit' line
#    y = mlab.normpdf( bins, np.mean(k_it), np.std(k_it))
#    l = plt.plot(bins, y, 'r--', linewidth=1)
    
    plt.xlabel('Number of regions', fontsize=18)
    plt.ylabel('Number of instances', fontsize=18)
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()
    
def country_emissions(ds_mcmc, countries, species, domain, x_post_vit=None, q_ap_abs_v=None, 
                      percentiles=[5,16,50,84,95], units=None, ocean=True, 
                      uk_split=False, fixed_map=False):
        
    """
    Generates national totals for a given list of countries
    Args: 
        ds_mcmc (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script). (post_mcmc xarray dataset)
        countries (list) :
            List of countries (see data_path/countries files for names)
            In general names are all capitalized.
        species (str) :
            CH4, CO2 etc.
        domain (str) :
            Domain name of interest.
            e.g. EUROPE, SOUTHAMERICA etc.
        x_post_vit (numpy.array/None, optional) :
            Posterior values (Dims: nIt x NGrid) (flattened on lat and lon dimensions).
            If not specified, will be extracted from ds_mcmc.
        q_ap_abs_v (numpy.array/None, optional) :
            Flux a priori values (Dims: NGrid) (flattened on lat and lon dimensions).
            If not specified, will be extracted from ds_mcmc.
        percentiles (list) :
            List of percentiles to extract.
            Default = [5,16,50,84,95]
        units (str/None, optional) :
            Need to specify to get something sensible out, otherwise returns in g/yr
            Options are Pg/yr, Tg/yr,Gg/yr Mg/yr
        ocean (bool, optional) :
            Default is True, If False only include emissions from land surface
        uk_split (bool, optional) :
            Break UK into devolved administrations
                Country names then have to be ['Eng', 'Sco', 'Wales',
                NIre', 'IRELAND', 'BENELUX']. Horribly inconsistent I know!
            Default = False
        fixed_map (bool, optional) :
            TODO: Add comment
    
    Returns: 
        (5 x numpy.array) :
            Country totals for each iteration in units specified (ncountries x nIt),
            Mean of country totals in units specified (ncountries), 
            Percentiles for each country in units specified (ncountries x npercentiles),
            Prior for each country in units specified (ncountries),
            Country index map (nlat x nlon)
       
    Output in Tg/yr
    """
    
    if x_post_vit is None:
        x_post_vit = ds_mcmc.x_post_vit.values
    if q_ap_abs_v is None:
        q_ap = ds_mcmc.q_ap.values
        q_ap_abs_v = np.ravel(q_ap)
    
    if units == 'Tg/yr':
        unit_factor=1.e12
    elif units == 'Gg/yr': 
        unit_factor=1.e9
    elif units == 'Pg/yr': 
        unit_factor=1.e15
    elif units == 'Mg/yr': 
        unit_factor=1.e6
    else:
        print 'Undefined units: outputting in g/yr - let this be a lesson to define your units'
        unit_factor=1.
        
#    acrg_path = os.getenv('ACRG_PATH')
#    with open(acrg_path + "/acrg_species_info.json") as f:
#        species_info=json.load(f)
#            
#    species_key = agage.synonyms(species, species_info)
#    
#    molmass = float(species_info[species_key]['mol_mass'])
    #units = species_info[species_key]['units']    
    
    molmass = molar_mass(species)

    lonmin=np.min(ds_mcmc.lon.values)
    lonmax=np.max(ds_mcmc.lon.values)
    latmin=np.min(ds_mcmc.lat.values)
    latmax=np.max(ds_mcmc.lat.values)
    area=areagrid(ds_mcmc.lat.values,ds_mcmc.lon.values)
    # GET COUNTRY DATA
    if uk_split == True:
        c_object=name.get_country(domain, ocean=True, ukmo=True, uk_split=uk_split)
    elif ocean == True:
        c_object=name.get_country(domain, ocean=True, ukmo=True, uk_split=False)
    else:
        c_object=name.get_country(domain, ocean=False)
    cds = xray.Dataset({'country': (['lat','lon'], c_object.country), 
                        'name' : (['ncountries'],c_object.name) },
                                        coords = {'lat': (c_object.lat),
                                        'lon': (c_object.lon)})
    country = cds.country.sel(lon=slice(lonmin,lonmax), 
                                        lat=slice(latmin,latmax))
    
    area_v=np.ravel(area)
    country_v=np.ravel(country)
    
    ncountries=len(countries)
    npercentiles = len(percentiles)
    #nIt=len(ds_mcmc.k_it.values)
    nIt=len(x_post_vit[:,0])

    nlon=len(ds_mcmc.lon)
    nlat=len(ds_mcmc.lat)
    country_v_new=np.zeros((nlon*nlat))

    if len(np.shape(q_ap_abs_v)) == 1:

        q_country_it = np.zeros((ncountries, nIt))
        q_country_mean=np.zeros((ncountries))
        q_country_percentile = np.zeros((ncountries,npercentiles))
#        q_country_50=np.zeros((ncountries))
#        q_country_05=np.zeros((ncountries))
#        q_country_95=np.zeros((ncountries))
#        q_country_16=np.zeros((ncountries))
#        q_country_84=np.zeros((ncountries))
        q_country_ap=np.zeros((ncountries))
        q_country_mode=np.zeros((ncountries))
    
        for ci,cc in enumerate(countries):
            name_country = np.where(cds.name == cc)
            c_index = np.where(country_v == name_country[0])
            country_v_new[c_index[0]]=ci+1
            for it in range(nIt):
                x_vit_temp=x_post_vit[it,:]
                if fixed_map == False:
                    q_country_it[ci,it]=np.sum(x_vit_temp[c_index[0]]*area_v[c_index[0]]
                                    *q_ap_abs_v[c_index[0]]) 
                elif fixed_map == True:
                    q_country_it[ci,it]=np.sum(x_vit_temp*area_v[c_index[0]]
                                    *q_ap_abs_v[c_index[0]]) 
            q_country_ap[ci] = np.sum(area_v[c_index[0]]*q_ap_abs_v[c_index[0]]
                                *365.*24.*3600.*molmass/unit_factor) # in Tg/yr
        
            #dum=np.histogram(q_country_it[ci,:], bins=100)                
            #q_country_mode[ci] = dum[1][dum[0]==np.max(dum[0])][0]
            
            q_country_mean[ci]=np.mean(q_country_it[ci,:])*365.*24.*3600.*molmass/unit_factor # in Tg/yr
            for pi,percentile in enumerate(percentiles):
                q_country_percentile[ci,pi] = np.percentile(q_country_it[ci,:],percentile)*365.*24.*3600.*molmass/unit_factor
            #q_country_50[ci]=np.percentile(q_country_it[ci,:],50)*365.*24.*3600.*molmass/unit_factor
            #q_country_05[ci]=np.percentile(q_country_it[ci,:],5)*365.*24.*3600.*molmass/unit_factor
            #q_country_95[ci]=np.percentile(q_country_it[ci,:],95)*365.*24.*3600.*molmass/unit_factor
            #q_country_16[ci]=np.percentile(q_country_it[ci,:],16)*365.*24.*3600.*molmass/unit_factor
            #q_country_84[ci]=np.percentile(q_country_it[ci,:],84)*365.*24.*3600.*molmass/unit_factor
        
        country_index = np.reshape(country_v_new, (nlat,nlon))   
        #return q_country_it*365.*24.*3600.*molmass/unit_factor,\
        #q_country_mean, q_country_05, q_country_16, q_country_50, q_country_84, \
        #q_country_95, q_country_ap, country_index
        return q_country_it*365.*24.*3600.*molmass/unit_factor,\
        q_country_mean, q_country_percentile, q_country_ap, country_index
        
    elif len(np.shape(q_ap_abs_v)) == 2:

        q_ap_abs_v = q_ap_abs_v.T
        ntimes = np.shape(q_ap_abs_v)[1]        
        
        q_country_it = np.zeros((ncountries, ntimes, nIt))
        q_country_mean=np.zeros((ncountries, ntimes))
#        q_country_50=np.zeros((ncountries, ntimes))
#        q_country_05=np.zeros((ncountries, ntimes))
#        q_country_95=np.zeros((ncountries, ntimes))
#        q_country_16=np.zeros((ncountries, ntimes))
#        q_country_84=np.zeros((ncountries, ntimes))
        q_country_percentile=np.zeros((ncountries, ntimes,npercentiles))
        q_country_ap=np.zeros((ncountries, ntimes))
    
        for ci,cc in enumerate(countries):
            name_country = np.where(cds.name == cc)
            c_index = np.where(country_v == name_country[0])
            country_v_new[c_index[0]]=ci+1
            for it in range(nIt):
                x_vit_temp=x_post_vit[it,:]
                if fixed_map == False:
                    q_country_it[ci,:,it]=np.sum(x_vit_temp[c_index[0],None]*area_v[c_index[0],None]
                                    *q_ap_abs_v[c_index[0],:], axis = 0) 
                elif fixed_map == True:
                    q_country_it[ci,:,it]=np.sum(x_vit_temp*area_v[c_index[0],None]
                                    *q_ap_abs_v[c_index[0],:], axis = 0)
            q_country_ap[ci,:] = np.sum(area_v[c_index[0],None]*q_ap_abs_v[c_index[0],:]
                                *365.*24.*3600.*molmass/unit_factor, axis = 0) # in Tg/yr
        
            q_country_mean[ci,:]=np.mean(q_country_it[ci,:,:], axis =1)*365.*24.*3600.*molmass/unit_factor # in Tg/yr
            for pi,percentile in enumerate(percentiles):
                q_country_percentile[ci,:,pi]=np.percentile(q_country_it[ci,:,:],percentile, axis =1)*365.*24.*3600.*molmass/unit_factor
#            q_country_50[ci,:]=np.percentile(q_country_it[ci,:,:],50, axis =1)*365.*24.*3600.*molmass/unit_factor
#            q_country_05[ci,:]=np.percentile(q_country_it[ci,:,:],5, axis =1)*365.*24.*3600.*molmass/unit_factor
#            q_country_95[ci,:]=np.percentile(q_country_it[ci,:,:],95, axis =1)*365.*24.*3600.*molmass/unit_factor
#            q_country_16[ci,:]=np.percentile(q_country_it[ci,:,:],16, axis =1)*365.*24.*3600.*molmass/unit_factor
#            q_country_84[ci,:]=np.percentile(q_country_it[ci,:,:],84, axis =1)*365.*24.*3600.*molmass/unit_factor
        

        country_index = np.reshape(country_v_new, (nlat,nlon))   
#        return q_country_it*365.*24.*3600.*molmass/unit_factor,\
#        q_country_mean, q_country_05, q_country_16, q_country_50, q_country_84, \
#        q_country_95, q_country_ap, country_index
        return q_country_it*365.*24.*3600.*molmass/unit_factor,\
        q_country_mean, q_country_percentile, q_country_ap, country_index

def country_emissions_mult(ds_list, countries, species, domain, x_post_vit=None, q_ap_abs_v=None, 
                      percentiles=[5,16,50,84,95], units=None, ocean=True, 
                      uk_split=False, fixed_map=False):
    '''
    Calculate country emissions across multiple datasets. Combine mean and percentiles into 
    arrays for all time points.
    See process.country_emissions() function for details of inputs
    Returns:
        (5 x numpy.array) :
            Country totals for each iteration in units specified (ncountries x nIt),
            Mean of country totals for each dataset in units specified (ncountries x ntime), 
            Percentiles for each country for each dataset in units specified (ncountries x npercentiles x ntime),
            Prior for each country in units specified (ncountries),
            Country index map (nlat x nlon)
    '''
    ncountries=len(countries)
    ntime = len(ds_list)
    npercentile=len(percentiles)
    
    # Constructed from all datasets
    country_mean = np.zeros((ncountries,ntime)) 
    country_percentile = np.zeros((ncountries,ntime,npercentile))

    for i,ds in enumerate(ds_list):
        country_out = country_emissions(ds, countries, species, domain, 
                                        percentiles=percentiles,units=units, 
                                        ocean=ocean, uk_split=uk_split)
        
        country_mean[:,i] = country_out[1]
        country_percentile[:,i,:] = country_out[2]

        if i == 0:
            # Should be the same for all datasets
            country_it = country_out[0]
            country_prior = country_out[3]
            country_index = country_out[4]

    return country_it,country_mean,country_percentile,country_prior,country_index

def prior_mf(ds):
    '''
    The prior_mf function calculates the y prior mole fraction
    
    Args:
        ds (xarray.Dataset) : output from the run_tdmcmc function
    
    Returns:
        np.array : y prior mole fraction   
    '''
    h_v_all = ds.h_v_all.values # nmeasure x NgridIC (Ngrid + nIC)
    
    #NgridIC = len(ds.NgridIC.values)
    #x = np.ones(NgridIC)
    #y_prior = np.dot(h_v_all,x)
    
    y_prior = np.sum(h_v_all,axis=1)
    
    return y_prior

def post_bc_inner(ds):
    '''
    The post_bc_inner function calculates the y posterior boundary conditions for the inner
    region for each iteration and the mean.
    This combines mole fraction contributions for boundary conditions on the whole domain, 
    any bias factors included and all fixed regions outside the inner region.
    
    Args:
        ds (xarray.Dataset) : output from the run_tdmcmc function
    
    Returns:
        (np.array,np.array) : y posterior inner bc iterations, y posterior inner bc mean
    '''

    nIC = ds.nIC.values         # number of initial conditions (basis functions + boundary conditions + bias)
    nIt = len(ds.nIt)           # number of iterations
    nmeasure = len(ds.nmeasure) # number of measurement points
    
    h_v_all = ds.h_v_all.values # nmeasure x NgridIC (Ngrid + nIC)
    x_it = ds.x_it.values       # nIt x kICmax (nIC + kmax  - maximum number of regions in sub-domain)
    
    y_bg_it = np.zeros((nIt,nmeasure))
    for it in range(nIt):
        y_bg_it[it,:]=np.dot(h_v_all[:,:nIC],x_it[it,:nIC]) # Find y posterior for each iteration
        #y_bg_it[it,:]=np.dot(h_v_all[:,:10],x_it[it,:10])
    
    y_bg_mean=np.mean(y_bg_it, axis=0) # Take the mean across all iterations
    
    return y_bg_it,y_bg_mean

def post_mf(ds):
    '''
    The post_mf function calculates the y posterior mole fraction ierations and mean
    
    Args:
        ds (xarray.Dataset) : output from the run_tdmcmc function
    
    Returns:
        (np.array,np.array) : y posterior mole fraction iterations, y posterior mole fraction mean
    '''

    nIC = ds.nIC.values         # number of initial conditions (basis functions + boundary conditions + bias)
    nIt = len(ds.nIt)           # number of iterations
    nmeasure = len(ds.nmeasure) # number of measurement points
    Ngrid = len(ds.Ngrid)       # number of positions in lat x lon grid
    
    h_v_all = ds.h_v_all.values       # nmeasure x NgridIC (Ngrid + nIC)
    x_post_vit = ds.x_post_vit.values # nIt x Ngrid
    x_it = ds.x_it.values             # nIt x kICmax (nIC + kmax  - maximum number of regions in sub-domain)
    
    x_post_all_it = np.zeros((nIt,Ngrid+nIC)) # Create x vector of dimesion NgridIC for every iteration
    x_post_all_it[:,:nIC] = x_it[:,:nIC] # Populate nIC entries with posterior values for each IC parameter for x_it
    x_post_all_it[:,nIC:] = x_post_vit   # Populate Ngrid entries with posterior values for each coord from x_post_vit
    
    y_post_it = np.zeros((nIt,nmeasure))
    for it in range(nIt):
        y_post_it[it,:]=np.dot(h_v_all,x_post_all_it[it,:]) # Find y posterior for each iteration
    
    y_post_mean=np.mean(y_post_it, axis=0) # Take the mean across all iterations
    
    return y_post_it,y_post_mean

def combine_timeseries(*ds_mult):
    '''
    The combine_timeseries function takes a list of output datasets from a tdmcmc run and combines
    the parameters relevant to plotting a mole fraction timeseries.

    Current parameters copied from input and combined:
        "y_time","y_site","y","sigma_y_it","sites"
    
    Posterior boundary conditions (inner) and modelled mole fractions are calculated for each run and
    combined. Both mean and full iterations are included. Additional parameters within dataset:
        "bc_inner_post","bc_inner_post_it","mf_post","mf_post_it"
    
    If any 
    
    Args:
        ds_mult (xarray.Dataset, xarray.Dataset, ...) :
            Any number of tdmcmc output datasets can be specified to be combined.
            All datasets will be combined.
    
    Returns:
        xarray.Dataset :
            Reduced, combined dataset from tdmcmc output.
    '''
    calc_data_vars = {"bc_inner_post":post_bc_inner,"mf_post":post_mf}
    #data_vars = ["y_time","y_site","y","sigma_y_it","bc_inner_post","mf_post","nIC","h_v_all","x_it","x_post_vit"]
    data_vars = ["y_time","y_site","y","sigma_y_it","bc_inner_post","mf_post"]
    match_coords = ["Ngrid","sites"]
    
    concat_dim = "nmeasure"
    iter_dim = "nIt"
    run_dim = "run"
    
    ## Check any coordinates within match_coords are the same for all the datasets.    
    match_dim = [[ds.dims[coord] for ds in ds_mult] for coord in match_coords]
    for ds in ds_mult:
        for i,dim in enumerate(match_dim):
            if len(set(dim)) != 1:
                raise Exception("Dimensions of {} do not match between input datasets.".format(match_coords[i]))
    
    ## Extract and combine data arrays for data_vars including calculating any new data variables to add
    data_arrays = OrderedDict()
    
    for dv in data_vars:
        for i,ds in enumerate(ds_mult):
            if dv in calc_data_vars:
                calc_dv_it,calc_dv = calc_data_vars[dv](ds)
                if i == 0:
                    da = xray.DataArray(calc_dv,dims=concat_dim)
                    da_it = xray.DataArray(calc_dv_it,dims=[iter_dim,concat_dim])
                else:
                    da = xray.concat([da,xray.DataArray(calc_dv,dims=concat_dim)],dim=concat_dim)
                    da_it = xray.concat([da_it,xray.DataArray(calc_dv_it,dims=[iter_dim,concat_dim])],dim=concat_dim)
            else:
                if i == 0:
                    da = ds[dv]
                else:
                    if concat_dim in da.dims:
                        da = xray.concat([da,ds[dv]],dim=concat_dim)
                    else:
                        da = xray.concat([da,ds[dv]],dim=run_dim)
        data_arrays[dv] = da
        
        # Add data arrays containing the full iteration history as well as the mean for bg and x_post
        if dv in calc_data_vars:
            data_arrays["{}_it".format(dv)] = da_it
    
    ## Create dataset from extracted data arrays
    combined_ds = xray.Dataset()
    
    for d in data_arrays.items():
        d = {d[0]:d[1]}
        combined_ds = combined_ds.assign(**d)
    
    # Add any coordinates from match_coords which are missing e.g. sites
    for coord_name in match_coords:
        if coord_name not in combined_ds.dims:
            combined_ds = combined_ds.assign_coords(**{coord_name:ds_mult[0].coords[coord_name]})
    
    # Add details the nmeasure value for each run (may be useful if we need to extract details for certain runs)
    run = [ds.attrs["Start date"] for ds in ds_mult]
    run_nmeasure = xray.DataArray([len(ds.nmeasure) for ds in ds_mult],coords={run_dim:run},dims=run_dim)
    
    combined_ds = combined_ds.assign(**{"run_nmeasure":run_nmeasure})
    combined_ds.attrs["Created by"] = getpass.getuser()
    combined_ds.attrs["File created"] = str(dt.datetime.now().replace(microsecond=0))

    return combined_ds
 
def plot_timeseries(ds, fig_text=None, ylim=None, out_filename=None, doplot=True):
    
    """
    Plot measurement timeseries of posterior and observed measurements
    Requires post_mcmc xray dataset
    For future: incorporate model & measurement uncertainty
    Plots separate subplots for each of the measurement sites - hopefully!
    
    ds = dataset output from run_tdmcmc script
    fig_text = String e.g. "CH$_{4}$ mole fraction (ppb)"
    ylim = array [ymin,ymax]
    
    Specify an out_filename to write to disk
    """
    
    if "run" in ds.coords:
        y_bg_mean = ds["bc_inner_post"].values
        y_bg_it = ds["bc_inner_post_it"].values
        y_post_mean = ds["mf_post"].values
        y_post_it = ds["mf_post_it"].values
    else:
        y_bg_it,y_bg_mean = post_bc_inner(ds)
        y_post_it,y_post_mean = post_mf(ds)

    sites = ds.coords['sites'].values
    nsites=len(sites)
    
    sigma_y_mean=np.mean(ds.sigma_y_it.values, axis=1)

    y_time = ds.y_time.values
    y_site = ds.y_site.values
    y_obs = ds.y.values
    upper = y_post_mean+sigma_y_mean
    lower = y_post_mean-sigma_y_mean
    

    if doplot is True:
        fig,ax=plt.subplots(nsites,sharex=True)
        
        if nsites > 1:
            for si,site in enumerate(sites):
                wh_site = np.where(y_site == site)
                ax[si].fill_between(y_time[wh_site[0]], upper[wh_site[0]],lower[wh_site[0]], alpha=0.6, 
                            facecolor='lightskyblue', edgecolor='lightskyblue')
                ax[si].plot(y_time[wh_site[0]],y_obs[wh_site[0]], 'ro', markersize=4, label='Observations')
                ax[si].plot(y_time[wh_site[0]],y_post_mean[wh_site[0]], color='blue', label='Modelled observations')
                ax[si].plot(y_time[wh_site[0]],y_bg_mean[wh_site[0]],color='black', 
                         label='Modelled baseline')
                if ylim is not None:
                    ax[si].set_ylim(ylim)
                start, end = ax[si].get_ylim()
                ax[si].yaxis.set_ticks(np.arange(start, end+1, (end-start)/5))
                ax[si].set_ylabel(site)
                if si == 0:
                    legend=ax[si].legend(loc='upper left')
                    for label in legend.get_texts():
                        label.set_fontsize('small')
                        
        else:
            ax.fill_between(y_time, upper,lower, alpha=0.6, 
                        facecolor='lightskyblue', edgecolor='lightskyblue')
            ax.plot(y_time,y_obs, 'ro', markersize=4, label='Observations')
            ax.plot(y_time,y_post_mean, color='blue', label='Modelled observations')
            ax.plot(y_time,y_bg_mean,color='black', 
                     label='Modelled baseline')
            start, end = ax.get_ylim()
            ax.yaxis.set_ticks(np.arange(start, end+1, (end-start)/5))
            legend=ax.legend(loc='upper left')
            for label in legend.get_texts():
                label.set_fontsize('small')
        
        if fig_text is not None:
            fig.text(0.01,0.65,fig_text, rotation=90)
        fig.autofmt_xdate()
        
        if out_filename is not None:
            plt.savefig(out_filename)
            plt.close()
        else:
            plt.show()
        
    return y_post_it, y_bg_it

def open_ds(path):
    
    """
    Function efficiently opens xray datasets.
    """
    # use a context manager, to ensure the file gets closed after use
    with xray.open_dataset(path) as ds:
        ds.load()
    return ds    

def plot_nuclei_density(ds, out_filename=None, fignum=3):
    
    """
    Plot a 2D histogram/density plot of where nuclei are located in inversion domain
    
    Requires x_ray dataset of outputs    
    Currently in rainbow colour scheme - Eurgh! #endrainbow
    
    Specify an out_filename to write to disk
    """
    ##############################
    # Density plot of nuclei locations
    lon=ds.lon.values
    lat=ds.lat.values
    plat_it=ds.plat_it.values
    plon_it=ds.plon_it.values
    
    dlon=lon[1]-lon[0]
    dlat=lat[1]-lat[0]
    nIt=len(plat_it[0,:])
    
    xedges=np.arange(np.min(lon),np.max(lon)+dlon,dlon)
    yedges=np.arange(np.min(lat),np.max(lat)+dlat,dlat)
    lon_range = (np.min(xedges), xedges[-2])
    lat_range = (np.min(yedges), yedges[-2])
    plat_it_v=np.ravel(plat_it)
    plon_it_v=np.ravel(plon_it)
    
    H, yedges, xedges = np.histogram2d(plat_it_v, plon_it_v, bins=(yedges, xedges))
    
    m2 = Basemap(projection='gall',
                llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
                llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
                resolution='l')
    
    lona, latb = np.meshgrid(xedges[:-1],yedges[:-1])
    mapa, mapb = m2(lona, latb)
    
    fig = plt.figure(5,figsize=(8,8))
    #ax = fig.add_subplot(111)
    
    m2.drawcoastlines()
    m2.drawcountries() 
    
    hlevels=np.arange(0, nIt/4, nIt/20/4) 
    
    ds=m2.contourf(mapa,mapb,H, hlevels, extend='both', cmap='YlGnBu')
    db = m2.colorbar(ds, location='bottom')
    plt.show()    
    
    
#    lon=ds.lon.values
#    lat=ds.lat.values
#    plat_it=ds.plat_it.values
#    plon_it=ds.plon_it.values
#    
#    xedges=np.arange(np.min(lon),np.max(lon)+2.,2.)
#    yedges=np.arange(np.min(lat),np.max(lat)+2.,2.)
#    lon_range = (np.min(xedges), xedges[-2])
#    lat_range = (np.min(yedges), yedges[-2])
#    plat_it_v=np.ravel(plat_it)
#    plon_it_v=np.ravel(plon_it)
#    
#    H, yedges, xedges = np.histogram2d(plat_it_v, plon_it_v, bins=(yedges, xedges))
#    
#    m2 = Basemap(projection='gall',
#                llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
#                llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
#                resolution='l')
#    
#    lona, latb = np.meshgrid(xedges[:-1],yedges[:-1])
#    mapa, mapb = m2(lona, latb)
#    
#    fig = plt.figure(fignum,figsize=(8,8))
#    ax = fig.add_subplot(111)
#    
#    m2.drawcoastlines()
#    m2.drawcountries() 
#    
#    hlevels=np.arange(np.min(H), np.max(H), (np.max(H)-np.min(H))/20) 
#    
#    ds=m2.contourf(mapa,mapb,H, hlevels)
#    db = m2.colorbar(ds, location='bottom')

    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()