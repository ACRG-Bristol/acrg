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
from mpl_toolkits.basemap import Basemap
import xray
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

def append_netcdf(flux_mean, flux_percentile, country_mean, country_percentile, 
                 lon, lat, time, country, percentile, experiment, outfile):
      
    """
    Function for appending flux and country output to file
    
    Only use this function once file has already been created
    Create file first using write_netcdf function below
    """         
      
    fm_name = "_".join(['flux_mean', experiment])
    fpc_name = "_".join(['flux_percentile', experiment])
    countrym_name = "_".join(['country_mean', experiment])
    countrypc_name = "_".join(['country_percentile', experiment])
              
    ncF=Dataset(outfile, 'a')
    
    ncfluxm=ncF.createVariable(fm_name, 'f', ('lat', 'lon', 'time'))    
    ncfluxpc=ncF.createVariable(fpc_name, 'f', ('lat', 'lon', 'time', 'percentile'))    
    nccountrym=ncF.createVariable(countrym_name, 'f', ('country', 'time'))   
    nccountrypc=ncF.createVariable(countrypc_name, 'f', ('country', 'time', 'percentile'))
    
    ncfluxm[:, :, :]=flux_mean
    ncfluxm.units='mol/m2/s'
    ncfluxm.long_name='Posterior mean flux from each grid-cell'
    
    ncfluxpc[:, :, :,:]=flux_percentile
    ncfluxpc.units='mol/m2/s'
    ncfluxpc.long_name='Posterior mean flux from each grid-cell'
    
    nccountrym[:, :]=country_mean
    nccountrym.units='Tg/yr'
    nccountrym.long_name='Posterior mean emissions from individual countries'
    
    nccountrypc[:, :, :]=country_percentile
    nccountrypc.units='Tg/yr'
    nccountrypc.long_name='Posterior percentile emissions from individaul countries'
        
    ncF.close()
    print "Appended " + experiment + " to " + outfile 

def write_netcdf(flux_mean, flux_percentile, flux_prior, flux_ap_percentile,
                 country_mean, country_percentile, country_prior, country_ap_percentile,
                 country_index,
                 lon, lat, time, country, percentile, experiment, outfile):
    
    """
    Function for writing all flux and country output to file
    
    Only use this function for first set of variables
    Once file has been created, use append_netcdf function above
    """    
    
    time_seconds, time_reference = time2sec(time)
    
    fm_name = "_".join(['flux_mean', experiment])
    fpc_name = "_".join(['flux_percentile', experiment])
    countrym_name = "_".join(['country_mean', experiment])
    countrypc_name = "_".join(['country_percentile', experiment])
    
    #Write NetCDF file
    ncF=Dataset(outfile, 'w')
    ncF.createDimension('time', len(time))
    ncF.createDimension('lon', len(lon))
    ncF.createDimension('lat', len(lat))
    ncF.createDimension('country', len(country))
    ncF.createDimension('percentile', len(percentile))
    
    nctime=ncF.createVariable('time', 'i', ('time',))
    nclon=ncF.createVariable('lon', 'f', ('lon',))
    nclat=ncF.createVariable('lat', 'f', ('lat',))
    nccountry=ncF.createVariable('country_name', 'str', ('country',))
    ncpercent=ncF.createVariable('percentile', 'i', ('percentile',))
    ncfluxm=ncF.createVariable(fm_name, 'f', ('lat', 'lon', 'time'))    
    ncfluxpc=ncF.createVariable(fpc_name, 'f', ('lat', 'lon', 'time', 'percentile'))    
    nccountrym=ncF.createVariable(countrym_name, 'f', ('country', 'time'))   
    nccountrypc=ncF.createVariable(countrypc_name, 'f', ('country', 'time', 'percentile'))
    nccountryap=ncF.createVariable('country_prior', 'f', ('country', 'time'))
    nccountryappc=ncF.createVariable('country_prior_percentile', 'f', ('country', 'time', 'percentile'))
    ncfluxap=ncF.createVariable('flux_prior', 'f', ('lat', 'lon', 'time'))  
    ncfluxappc=ncF.createVariable('flux_prior_percentile', 'f', ('lat', 'lon', 'time', 'percentile'))  
    nccountryid=ncF.createVariable('country_index', 'i', ('lat', 'lon'))    

    
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
    
    ncF.close()
    print "Written " + experiment + " to " + outfile

def plot_scaling(data,lon,lat, out_filename=None, absolute=False, 
                 uncertainty=False, uncertainty_red=False,
                 stations=None, fignum=None):
    
    """
    Plot 2d scaling map of posterior x
    i.e. degree of scaling applied to prior emissions 
    data = mean or median of x_post_vit
    """    
    lon_range = (np.min(lon), np.max(lon))
    lat_range = (np.min(lat), np.max(lat))
    m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')

    lons, lats = np.meshgrid(lon,lat)
    mapx, mapy = m(lons, lats)
    
    fig = plt.figure(fignum,figsize=(8,8))
    fig.add_axes([0.1,0.1,0.8,0.8])
    
    m.drawcoastlines()
    m.drawcountries() 
    
    def draw_screen_poly( lats, lons, m):
            x, y = m( lons, lats )
            xy = zip(x,y)
            poly = Polygon( xy, facecolor='red', alpha=0.4 )
            plt.gca().add_patch(poly)
    

    if absolute == True:
        
        clevels = np.arange(-2., 2.1, 0.1)
        #cs = m.contourf(mapx,mapy,data, clevels, extend='both', cmap='RdBu_r')
        #cb = m.colorbar(cs, location='bottom', pad="5%")
        cmap = plt.cm.RdBu_r
        norm = BoundaryNorm(clevels,
                        ncolors=cmap.N,
                        clip=True)
        cs = m.pcolormesh(mapx, mapy,
                data,cmap =cmap, norm = norm)
        cb = m.colorbar(cs, location='bottom', pad="5%", extend='both')
        #cb.set_label('Difference from prior (Tg yr$^{-1}$)')
        cb.set_label('Posterior - prior (kg/m$^{2}$/s)x$10^{9}$')

    elif uncertainty == True:
        clevels = np.arange(0., 4.1, 0.1) 
        cmap = plt.cm.YlGnBu
        norm = BoundaryNorm(clevels,
                        ncolors=cmap.N,
                        clip=True)
        cs = m.pcolormesh(mapx, mapy,
                data,cmap =cmap, norm = norm)
        cb = m.colorbar(cs, location='bottom', pad="5%", extend='both')
        cb.set_label('Normalized uncertainty')
    elif uncertainty_red == True:
        clevels = np.arange(0., 1.05, 0.05)  
        cmap = plt.cm.YlGnBu
        norm = BoundaryNorm(clevels,
                        ncolors=cmap.N,
                        clip=True)
        cs = m.pcolormesh(mapx, mapy,
                data,cmap =cmap, norm = norm)
        cb = m.colorbar(cs, location='bottom', pad="5%", extend='both')
        cb.set_label('Normalized uncertainty reduction')
    else:
        clevels = np.arange(0., 2.05, 0.05)  
        #cs = m.contourf(mapx,mapy,data, clevels, extend='both', cmap='RdBu_r')
        #cb = m.colorbar(cs, location='bottom', pad="5%") 
        cmap = plt.cm.RdBu_r
        norm = BoundaryNorm(clevels,
                        ncolors=cmap.N,
                        clip=True)
        cs = m.pcolormesh(mapx, mapy,
                data,cmap =cmap, norm = norm)
        cb = m.colorbar(cs, location='bottom', pad="5%", extend='both')
        cb.set_label('Scaling of prior') 
        
        
    if stations is not None:
        
        nsites=len(stations['sites'])
        sites=stations['sites']
        ilon=np.zeros((nsites),dtype=np.uint16)
        ilat=np.zeros((nsites),dtype=np.uint16)
        for si,site in enumerate(stations['sites']):
            ilon[si]=stations[site+'_lon']
            ilat[si]=stations[site+'_lat']
        mlon,mlat=m(lon[ilon],lat[ilat])
        site_loc=m.plot(mlon,mlat, linestyle='None',
                        marker='o', color='black', markersize=12)
        yoffset = 0.022*(m.ymax-m.ymin) 
        xoffset = 0.012*(m.xmax-m.xmin)             
        for ii in range(nsites):
            plt.text(mlon[ii]+xoffset,
                     mlat[ii]+yoffset,sites[ii], fontsize='24', color='black')   
                     
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
    
def country_emissions(ds_mcmc, x_post_vit, x_ap_abs_v, countries, species,
                      ocean=True, domain='EUROPE', al=False):
        
    """
    Generates national totals for a given list of countries
    Requires post_mcmc xray dataset
    
    Units are hard-wired at the moment, need to make non-methane specific
    Returns: mean, 5th,16th,median,84th,95th percentiles and prior for each country
    
    Output in Tg/yr
    """
    temp = os.path.split(os.path.realpath(__file__))
    acrg_path=os.path.join(temp[0],"..")
    with open(acrg_path + "/acrg_species_info.json") as f:
        species_info=json.load(f)
            
    species_key = agage.synonyms(species, species_info)
    
    molmass = float(species_info[species_key]['mol_mass'])
    units = species_info[species_key]['units']    
    

    lonmin=np.min(ds_mcmc.lon.values)
    lonmax=np.max(ds_mcmc.lon.values)
    latmin=np.min(ds_mcmc.lat.values)
    latmax=np.max(ds_mcmc.lat.values)
    area=areagrid(ds_mcmc.lat.values,ds_mcmc.lon.values)
    # GET COUNTRY DATA
    if ocean==True:
        c_object=name.get_country(domain, ocean=True, al=al)
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
    #nIt=len(ds_mcmc.k_it.values)
    nIt=len(x_post_vit[:,0])
    q_country_it = np.zeros((ncountries, nIt))
    q_country_mean=np.zeros((ncountries))
    q_country_50=np.zeros((ncountries))
    q_country_05=np.zeros((ncountries))
    q_country_95=np.zeros((ncountries))
    q_country_16=np.zeros((ncountries))
    q_country_84=np.zeros((ncountries))
    q_country_ap=np.zeros((ncountries))
    
    nlon=len(ds_mcmc.lon)
    nlat=len(ds_mcmc.lat)
    country_v_new=np.zeros((nlon*nlat))
    for ci,cc in enumerate(countries):
        name_country = np.where(cds.name == cc)
        c_index = np.where(country_v == name_country[0])
        country_v_new[c_index[0]]=ci+1
        for it in range(nIt):
            x_vit_temp=x_post_vit[it,:]
            q_country_it[ci,it]=np.sum(x_vit_temp[c_index[0]]*area_v[c_index[0]]
                                *x_ap_abs_v[c_index[0]]) 
        
        q_country_ap[ci] = np.sum(area_v[c_index[0]]*x_ap_abs_v[c_index[0]]
                                *365.*24.*3600./1.e9*molmass/1000.) # in Tg/yr
        
        
        q_country_mean[ci]=np.mean(q_country_it[ci,:])*365.*24.*3600./1.e9*molmass/1000. # in Tg/yr
        q_country_50[ci]=np.percentile(q_country_it[ci,:],50)*365.*24.*3600./1.e9*molmass/1000.
        q_country_05[ci]=np.percentile(q_country_it[ci,:],5)*365.*24.*3600./1.e9*molmass/1000.
        q_country_95[ci]=np.percentile(q_country_it[ci,:],95)*365.*24.*3600./1.e9*molmass/1000.
        q_country_16[ci]=np.percentile(q_country_it[ci,:],16)*365.*24.*3600./1.e9*molmass/1000.
        q_country_84[ci]=np.percentile(q_country_it[ci,:],84)*365.*24.*3600./1.e9*molmass/1000.
        
    country_index = np.reshape(country_v_new, (nlat,nlon))   
    return q_country_it*365.*24.*3600./1.e9*molmass/1000.,\
    q_country_mean, q_country_05, q_country_16, q_country_50, q_country_84, \
    q_country_95, q_country_ap, country_index
    
def plot_timeseries(ds, species, out_filename=None, doplot=True):
    
    """
    Plot measurement timeseries of posterior and observed measurements
    Requires post_mcmc xray dataset
    For future: incorporate model & measurement uncertainty
    Plots separate subplots for each of the measurement sites - hopefully!
    
    Specify an out_filename to write to disk
    """
    x_it=ds.x_it.values
    h_v_all=ds.h_v_all.values
    x_post_vit = ds.x_post_vit.values
    sigma_y_mean=np.mean(ds.sigma_y_it.values, axis=1)
    sites=ds.coords['sites'].values
    nsites=len(sites)
    nlon=len(ds.lon)
    nlat=len(ds.lat)
    Ngrid=nlon*nlat
    nIt=len(ds.nIt)
    nIC=ds.nIC.values
    nmeasure=len(ds.nmeasure)

    x_post_all_it=np.zeros((nIt,Ngrid+nIC))
    y_post_it = np.zeros((nIt,nmeasure))
    y_bg_it = np.zeros((nIt,nmeasure))
    x_post_all_it[:,:nIC]=x_it[:,:nIC]
    x_post_all_it[:,nIC:]=x_post_vit
    for it in range(nIt):
        y_post_it[it,:]=np.dot(h_v_all,x_post_all_it[it,:])  
        y_bg_it[it,:]=np.dot(h_v_all[:,:nIC],x_it[it,:nIC])
    y_post_mean=np.mean(y_post_it, axis=0)
    y_bg_mean=np.mean(y_bg_it, axis=0)
    y_time=ds.y_time.values
    y_site = ds.y_site.values
    y_obs = ds.y.values
    upper=y_post_mean+sigma_y_mean
    lower=y_post_mean-sigma_y_mean
    
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
                start, end = ax[si].get_ylim()
                ax[si].yaxis.set_ticks(np.arange(start, end+1, 100))
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
            
        fig.text(0.01,0.65,'CH$_{4}$ mole fraction (ppb)', rotation=90)
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