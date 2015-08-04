# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 08:57:53 2015


Script to process Transdimensional MCMC output 

Includes:
Write netcdf and append netcdf to write output to nc file

plot_scaling - plot posterior scaling map

regions_histogram - plot histogram of number of regions

country_emissions - calculate emissions from given list of countries
                    Currnetly hard-wired for methane

Requires:

Output files from tdmcmc_template.py stored in the form:

"output_" + network + "_" + species +  "_" + months[-1] + ".nc"

@author: ml12574
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xray
import acrg_name_xray as name
from acrg_grid import areagrid
import glob
import pandas
import matplotlib.dates as mdates
#import datetime
from netCDF4 import Dataset
from acrg_time.convert import time2sec
import os
import acrg_agage as agage
import json

def append_netcdf(flux_mean, flux_percentile, country_mean, country_percentile, 
                 lon, lat, time, country, percentile, experiment, outfile):
       
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
    ncfluxm.long_name='Mean flux from each grid-cell'
    
    ncfluxpc[:, :, :,:]=flux_percentile
    ncfluxpc.units='mol/m2/s'
    ncfluxpc.long_name='Mean flux from each grid-cell'
    
    nccountrym[:, :]=country_mean
    nccountrym.units='Tg/yr'
    nccountrym.long_name='Mean emissions from individual countries'
    
    nccountrypc[:, :, :]=country_percentile
    nccountrypc.units='Tg/yr'
    nccountrypc.long_name='Percentile emissions from individaul countries'
        
    ncF.close()
    print "Appended " + experiment + " to " + outfile 

def write_netcdf(flux_mean, flux_percentile, country_mean, country_percentile, 
                 lon, lat, time, country, percentile, experiment, outfile):
    
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
    ncfluxm.long_name='Mean flux from each grid-cell'
    
    ncfluxpc[:, :, :,:]=flux_percentile
    ncfluxpc.units='mol/m2/s'
    ncfluxpc.long_name='Mean flux from each grid-cell'
    
    nccountrym[:, :]=country_mean
    nccountrym.units='Tg/yr'
    nccountrym.long_name='Mean emissions from individual countries'
    
    nccountrypc[:, :, :]=country_percentile
    nccountrypc.units='Tg/yr'
    nccountrypc.long_name='Percentile emissions from individaul countries'
        
    ncF.close()
    print "Written " + experiment + " to " + outfile

def plot_scaling(data,lon,lat, out_filename=None, fignum=1):
    
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
    clevels = np.arange(0., 2.0, 0.02)  
 
    cs = m.contourf(mapx,mapy,data, clevels, extend='both', cmap='RdBu_r')
    cb = m.colorbar(cs, location='bottom', pad="5%")
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()
    
    
def regions_histogram(k_it, out_filename=None, fignum=2):
    
    bin_max=np.ceil((np.max(k_it)+50.)/40.)*40
    bin_step = bin_max/20
    bins=np.arange(0,bin_max,bin_step)    
    
    plt.figure(fignum)
    # the histogram of the data with histtype='step'
    n, bins, patches = plt.hist(k_it, bins, normed=1, histtype='bar', rwidth=0.8)
    plt.xlabel('# of regions')
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()
    
def country_emissions(ds_mcmc, x_post_vit, x_ap_abs_v, countries, species, domain='EUROPE'):
    
    acrg_path="/home/ml12574/work/programs/Python/acrg"
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
    c_object=name.get_country(domain, ocean=True)
    cds = xray.Dataset({'country': (['lat','lon'], c_object.country), 
                        'name' : (['ncountries'],c_object.name) },
                                        coords = {'lat': (c_object.lat),
                                        'lon': (c_object.lon)})
    country = cds.country.sel(lon=slice(lonmin,lonmax), 
                                        lat=slice(latmin,latmax))
    
    area_v=np.ravel(area)
    country_v=np.ravel(country)
    
    ncountries=len(countries)
    nIt=len(ds_mcmc.k_it.values)
    q_country_it = np.zeros((ncountries, nIt))
    q_country_mean=np.zeros((ncountries))
    q_country_50=np.zeros((ncountries))
    q_country_05=np.zeros((ncountries))
    q_country_95=np.zeros((ncountries))
    q_country_16=np.zeros((ncountries))
    q_country_84=np.zeros((ncountries))
    q_country_ap=np.zeros((ncountries))
    
    for ci,cc in enumerate(countries):
        name_country = np.where(cds.name == cc)
        c_index = np.where(country_v == name_country[0])
        
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
        
        
    return q_country_mean, q_country_05, q_country_16, q_country_50, q_country_84, q_country_95, q_country_ap
    
def plot_timeseries(ds, species, out_filename=None):
    
    
    x_it=ds.x_it.values
    h_v_all=ds.h_v_all.values
    x_post_vit = ds.x_post_vit.values
    nlon=len(ds.lon)
    nlat=len(ds.lat)
    Ngrid=nlon*nlat
    nIt=len(ds.nIt)
    nIC=ds.nIC.values
    nmeasure=len(ds.nmeasure)
    
    y_obs = ds.y.values
    
    x_post_all_it=np.zeros((nIt,Ngrid+nIC))
    y_post_it = np.zeros((nIt,nmeasure))
    y_bg_it = np.zeros((nIt,nmeasure))
       
    x_post_all_it[:,:nIC]=x_it[:,:nIC]
    x_post_all_it[:,nIC:]=x_post_vit
    
    for it in range(nIt):
        y_post_it[it,:]=np.dot(h_v_all,x_post_all_it[it,:])  
        y_bg_it[it,:]=np.dot(h_v_all[:,:nIC],x_it[it,:nIC])
    
    y_post_mean=np.mean(y_post_it, axis=0)
    y_post_05 = np.percentile(y_post_it, 5, axis=0)  
    y_post_50 = np.percentile(y_post_it, 50, axis=0)
    y_post_95 = np.percentile(y_post_it, 95, axis=0)
    
    y_bg_mean=np.mean(y_bg_it, axis=0)
    y_bg_05 = np.percentile(y_bg_it, 5, axis=0)  
    y_bg_50 = np.percentile(y_bg_it, 50, axis=0)
    y_bg_95 = np.percentile(y_bg_it, 95, axis=0)
    
    
    fig,ax=plt.subplots()
    #time = pandas.to_datetime(ds.time.values)
    ax.fill_between(np.arange(nmeasure), y_post_05,y_post_95, alpha=0.2, 
                facecolor='skyblue', edgecolor='skyblue')
    ax.plot(np.arange(nmeasure),y_post_mean, color='blue', label='y posterior')
    
    ax.plot(np.arange(nmeasure),y_obs, color='green', label='y_obs')
    
    ax.set_ylabel('Mole Fraction (ppb)')
    legend=ax.legend(loc='upper left')
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('small')
    
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()

    

def open_ds(path):
        # use a context manager, to ensure the file gets closed after use
        with xray.open_dataset(path) as ds:
            ds.load()
        return ds    
#####################################################################          
# Read in post_mcmc dataset

########### INPUTS ####################

#months=["2013-07-01", "2013-08-01", "2013-09-01", "2013-10-01", 
#       "2013-11-01", "2013-12-01", "2014-01-01", "2014-02-01",
#       "2014-03-01", "2014-04-01", "2014-05-01", "2014-06-01"]
months=["2014-03-01"]
output_directory = "/home/ml12574/work/programs/Python/td-mcmc/outputs/"

species="CH4"
network='DECC'
experiment="MHD_TAC_RGL_TTA"
countries=np.asarray(['UNITED KINGDOM', 'IRELAND', 'FRANCE', 'GERMANY', 
                      'DENMARK', 'BELGIUM', 'NETHERLANDS', 'LUXEMBOURG'])
percentile = [5,16,50,84,95]
outfile="flux_NAME-Bristol_ch4.nc"

# SUBROUTINE OPTIONS

write_outfile=False
append_outfile=False
calc_country=True
plot_scale_map=True
plot_regions = True
plot_y_timeseries=True



print 'Beginning post processing'


ncountries=len(countries)
ntime = len(months)
npercentile=len(percentile)

country_mean=np.zeros((ncountries,ntime))  
country_percentile=np.zeros((ncountries,ntime,npercentile)) 
country_50=np.zeros((ncountries,ntime))  
country_05=np.zeros((ncountries,ntime))  
country_95=np.zeros((ncountries,ntime))  
country_16=np.zeros((ncountries,ntime))  
country_84=np.zeros((ncountries,ntime))  

country_prior=np.zeros((ncountries,ntime))  

f=glob.glob(output_directory + \
        "output_" + network + "_" + species + "_" + months[-1] + ".nc")
if len(f) > 0:
    ds0 = open_ds(f[0]) 
    nlon=len(ds0.lon)
    nlat=len(ds0.lat)
else:
    raise LookupError("Try a different base file to get nlon and nlat")

flux_mean = np.zeros((nlat,nlon,ntime))
flux_percentile = np.zeros((nlat,nlon,ntime,npercentile))

files = []
#q_country={}
q_country={}
for tt,ti in enumerate(months):
        f=glob.glob(output_directory + \
        "output_" + network + "_" + species + "_" + ti + ".nc")
        if len(f) > 0:
            files += f
            ds = open_ds(f[0])
    # Process x_it, and y_it etc.
        
    
            nlon=len(ds.lon)
            nlat=len(ds.lat)
            Ngrid=nlon*nlat
            nIt=len(ds.nIt)
            nIC=ds.nIC.values
            nmeasure=len(ds.nmeasure)
            
            x_post_vit=np.zeros((nIt,Ngrid))
            
            sigma_y_mean=np.zeros((nmeasure))
            
            sigma_y_mean=np.mean(ds.sigma_y_it.values, axis=1)
            
            regions_it=ds.regions_it.values
            x_it=ds.x_it.values
            #q_agg_it=ds.q_agg_it.values
            k_it=ds.k_it.values
            x_post_vit=ds.x_post_vit.values
            x_ap_abs=ds.x_ap_abs.values
            x_ap_abs_v=np.ravel(x_ap_abs)
            
    
            x_post_v_mean=np.mean(x_post_vit, axis=0)
            x_post_v_05 = np.percentile(x_post_vit, 5, axis=0)
            x_post_v_16 = np.percentile(x_post_vit, 16, axis=0)
            x_post_v_50 = np.percentile(x_post_vit, 50, axis=0)
            x_post_v_84 = np.percentile(x_post_vit, 84, axis=0)
            x_post_v_95 = np.percentile(x_post_vit, 95, axis=0)
            

            flux_mean[:,:,tt] = np.reshape(x_post_v_mean*x_ap_abs_v, (nlat,nlon))
            flux_percentile[:,:,tt,0] = np.reshape(x_post_v_05*x_ap_abs_v, (nlat,nlon))
            flux_percentile[:,:,tt,1] = np.reshape(x_post_v_16*x_ap_abs_v, (nlat,nlon))
            flux_percentile[:,:,tt,2] = np.reshape(x_post_v_50*x_ap_abs_v, (nlat,nlon))
            flux_percentile[:,:,tt,3] = np.reshape(x_post_v_84*x_ap_abs_v, (nlat,nlon))
            flux_percentile[:,:,tt,4] = np.reshape(x_post_v_95*x_ap_abs_v, (nlat,nlon))
            
            if calc_country == True:
                country_mean[:,tt], country_05[:,tt], country_16[:,tt], \
                country_50[:,tt], country_84[:,tt], country_95[:,tt], country_prior[:,tt] \
                = country_emissions(ds,x_post_vit, x_ap_abs_v,countries, species)        
            
            if plot_scale_map == True:
                lon=np.asarray(ds.lon.values)
                lat=np.asarray(ds.lat.values)
                x_post_mean = np.reshape(x_post_v_mean, (nlat,nlon))
                plot_scaling(x_post_mean, lon,lat, fignum=tt)

            if plot_regions == True:
                regions_histogram(k_it, fignum=tt+20)
                
            if plot_y_timeseries == True:
                plot_timeseries(ds, species)

    
country_mean[country_mean==0.]=np.nan 
country_05[country_05==0.]=np.nan
country_16[country_16==0.]=np.nan
country_50[country_50==0.]=np.nan
country_84[country_84==0.]=np.nan
country_95[country_95==0.]=np.nan

country_percentile[:,:,0] = country_05
country_percentile[:,:,1] = country_16
country_percentile[:,:,2] = country_50
country_percentile[:,:,3] = country_84
country_percentile[:,:,4] = country_95


d0=pandas.to_datetime(months)
lon=np.asarray(ds.lon.values)
lat=np.asarray(ds.lat.values)


#%%
#outfile="/home/ml12574/work/programs/Python/td-mcmc/flux_NAME-Bristol_ch4.nc"
outfile="flux_NAME-Bristol_ch4.nc"
if write_outfile == True:
    write_netcdf(flux_mean, flux_percentile,  country_mean, country_percentile, 
                     lon, lat, d0, countries, percentile, experiment, outfile)
elif append_outfile == True:
    append_netcdf(flux_mean, flux_percentile, country_mean, country_percentile, 
                     lon, lat, d0, countries, percentile, experiment, outfile)








