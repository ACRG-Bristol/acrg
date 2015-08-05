# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:24:46 2015

Template file for creating plots with output of tdmcmc

Uses tdmcmc_post_process


@author: ml12574
"""

########### INPUTS ####################

import numpy as np
import tdmcmc_post_process as process
import glob
import pandas


dates=["2014-03-01"]
species="CH4"
output_directory = "/path/to/tdmcmmc/outputs/"
#output_directory = "/home/ml12574/work/programs/Python/td-mcmc/outputs/"
outfile="outfile_name.nc"
network='DECC'
experiment="MHD_TAC_RGL_TTA"
countries=np.asarray(['UNITED KINGDOM', 'IRELAND', 'FRANCE', 'GERMANY', 
                      'DENMARK', 'BELGIUM', 'NETHERLANDS', 'LUXEMBOURG'])
percentile = [5,16,50,84,95]

# POST-PROCESSING OPTIONS
write_outfile=False
append_outfile=False
calc_country=True
plot_scale_map=True
plot_regions = True
plot_y_timeseries=True
plot_density = True


#results = post_process(species, dates, network, output_directory, countries=countries,
#                       write_outfile=False, append_outfile=False, calc_country=True,
#                       plot_scale_map=True, plot_regions=True, plot_y_timeseries=True)

print 'Beginning post processing'


ncountries=len(countries)
ntime = len(dates)
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
        "output_" + network + "_" + species + "_" + dates[-1] + ".nc")
if len(f) > 0:
    ds0 = process.open_ds(f[0]) 
    nlon=len(ds0.lon)
    nlat=len(ds0.lat)
else:
    raise LookupError("Try a different base file to get nlon and nlat")

flux_mean = np.zeros((nlat,nlon,ntime))
flux_percentile = np.zeros((nlat,nlon,ntime,npercentile))

files = []
#q_country={}
q_country={}
for tt,ti in enumerate(dates):
        f=glob.glob(output_directory + \
        "output_" + network + "_" + species + "_" + ti + ".nc")
        if len(f) > 0:
            files += f
            ds = process.open_ds(f[0])
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
                = process.country_emissions(ds,x_post_vit, x_ap_abs_v,countries, species)        
            
            if plot_scale_map == True:
                lon=np.asarray(ds.lon.values)
                lat=np.asarray(ds.lat.values)
                x_post_mean = np.reshape(x_post_v_mean, (nlat,nlon))
                process.plot_scaling(x_post_mean, lon,lat, fignum=tt, out_filename=None)

            if plot_regions == True:
                process.regions_histogram(k_it, fignum=tt+20, out_filename=None)
                
            if plot_y_timeseries == True:
                process.plot_timeseries(ds, species, out_filename=None)
                
            if plot_density == True:
                process.plot_nuclei_density(ds, out_filename=None, fignum=tt+40)

    
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

d0=pandas.to_datetime(dates)
lon=np.asarray(ds.lon.values)
lat=np.asarray(ds.lat.values)


#%%
#outfile="/home/ml12574/work/programs/Python/td-mcmc/flux_NAME-Bristol_ch4.nc"
outfile="flux_NAME-Bristol_ch4.nc"
if write_outfile == True:
    process.write_netcdf(flux_mean, flux_percentile,  country_mean, country_percentile, 
                     lon, lat, d0, countries, percentile, experiment, outfile)
elif append_outfile == True:
    process.append_netcdf(flux_mean, flux_percentile, country_mean, country_percentile, 
                     lon, lat, d0, countries, percentile, experiment, outfile)