# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:24:46 2015

Template file for creating plots with output of tdmcmc

Uses tdmcmc_post_process


@author: ml12574 (updated by rt17603)
"""
from __future__ import print_function
from __future__ import absolute_import

########### INPUTS ####################

import numpy as np
import hbmcmc_post_process as process
import pandas
import os
import sys
import matplotlib.pyplot as plt

from acrg_config.paths import paths
acrg_path = paths.acrg


if __name__=="__main__":

    #### GENERAL INPUTS ####

    dates=[] # Can be a list of one date or many dates
    species="species"
    domain="domain"
    runname='runname'
    
    output_directory = "path_to_output" # ** UPDATE OUTPUT DIRECTORY **
    
    #### POST-PROCESSING OPTIONS ####
    
    calc_country=False
    plot_scale_map=True
    plot_abs_map=True
    plot_diff_map=True
    plot_y_timeseries=True
    plot_countryemissions=True
    
    #### PARAMETERS FOR EACH POST-PROCESSING OPTION ####
    # Only used if relevant value above is set to True
    
    # calc_country
    country_file = "path_to_country_file"
    country_unit_prefix = None # Specify prefix, e.g., 'T' for teragram, default of None is 'g'
    countries = None # If not None, specify as list e.g., ['country1', 'country2'], default is all listed in country_file
    
    # plot_scale_map
    grid_scale_map = True
    s_clevels = None # Set to None to set to defaults.
    s_cmap = plt.cm.YlGnBu
    s_smooth = True
    s_out_filename = None  # None means plot will not be written to file

    # plot_abs_map
    grid_abs_map = True
    a_clevels = None  # Set to None to set to defaults.
    a_cmap = plt.cm.YlGnBu
    a_smooth = False
    a_out_filename = None  # None means plot will not be written to file
    
    # plot_diff_map
    grid_diff_map = True
    d_clevels = None # Set to None to set to defaults.
    d_cmap = plt.cm.RdBu_r
    d_smooth = False
    d_out_filename = None  # None means plot will not be written to file
    
    # plot_y_timeseries
    # combine timeseries will currently fail as needs to be re-written (ag12733 8/5/20)
    combine_timeseries = False # Plot y timeseries on one axis for multiple input files.
    y_out_filename = None

    # plot_countryemissions
    countries_to_plot = [] # list of countries
    CI_to_plot = 95 # 95 or 68 CI
    
    #### IMPLEMENT PROCESSING OPTIONS ####
    print('Beginning post processing')
    
    if output_directory == "/path/to/output/directory/":
        raise Exception("Please set output directory.")
    if not os.path.isdir(output_directory):
        raise Exception("Output directory: {} does not exist.".format(output_directory))
    
    # Extract datasets from file
    ds_list,filenames = process.extract_hbmcmc_files(output_directory,species, domain, runname,dates,return_filenames=True)
    dates = process.check_missing_dates(filenames,dates)
    
    lat = ds_list[0]["lat"]
    lon=ds_list[0]["lon"]

##  set lonmin, lonmax, latmin, latmax below if you only want to plot a subset of domain
##  need to find a better way to set these inputs at top
#     lonmin_ds = lon[np.where(np.isclose(lon, lonmin, atol = 0.3, rtol=0))[0][0]]
#     lonmax_ds = lon[np.where(np.isclose(lon, lonmax, atol = 0.3, rtol=0))[0][0]]
#     latmin_ds = lat[np.where(np.isclose(lat, latmin, atol = 0.2, rtol=0))[0][0]]
#     latmax_ds = lat[np.where(np.isclose(lat, latmax, atol = 0.2, rtol=0))[0][0]]
     
#     lon = lon.sel(lon=slice(lonmin_ds,lonmax_ds))
#     lat = lat.sel(lat=slice(latmin_ds,latmax_ds))
    
    ## Plot scaling map
    if plot_scale_map == True:
        process.plot_scale_map(ds_list, lat = lat, lon = lon, grid=grid_scale_map,clevels=s_clevels,
                               cmap=s_cmap,labels=None,title=None,smooth=s_smooth,
                               out_filename=s_out_filename,extend="both")
    
    ## Plot absolute difference map
    if plot_diff_map == True:
        process.plot_diff_map(ds_list,species, lat = lat, lon = lon,grid=grid_diff_map,clevels=d_clevels, 
                               cmap=d_cmap,labels=None,title=None,smooth=d_smooth,
                               out_filename=d_out_filename,extend="both")

    ## Plot absolute map
    if plot_abs_map == True:
        process.plot_abs_map(ds_list,species, lat = lat, lon = lon,grid=grid_abs_map,clevels=a_clevels, 
                               cmap=a_cmap,labels=None,title=None,smooth=a_smooth,
                               out_filename=a_out_filename,extend="max")

     
    ## Plot y timeseries
    if plot_y_timeseries == True:
        # combine_timeseries will currently fail as needs to be re-written for new outputs (ag12733 1/5/20)
        if combine_timeseries and len(ds_list) > 1:
            ds_combined = process.combine_timeseries(*ds_list)
            process.plot_timeseries(ds_combined, fig_text=None, 
                                                          ylim=None, out_filename=y_out_filename)
        else:
            for i,ds in enumerate(ds_list):
                if y_out_filename:
                    stub,ext = os.path.splitext(y_out_filename)
                    y_out_filename_n = "{}_{}{}".format(stub,i+1,ext)
                else:
                    y_out_filename_n = None
                process.plot_timeseries(ds, fig_text=None, 
                                                          ylim=None, out_filename=y_out_filename_n)
 

    ## Calculate country or area totals defined in country_file or the standard country defintion
    if calc_country == True:
        cntrymean_arr, cntry68_arr, cntry95_arr, cntryprior_arr \
         = process.country_emissions_mult(ds_list, species, domain, \
                                          country_file=country_file, country_unit_prefix=country_unit_prefix, countries = countries)
        
        units = country_unit_prefix + 'g'
    else:
        countries = list(ds["countrynames"].values)
        
        cntrymean_arr = np.zeros((len(ds_list), len(countries)))
        cntry68_arr = np.zeros((len(ds_list),len(countries),2))
        cntry95_arr = np.zeros((len(ds_list),len(countries),2))
        cntryprior_arr = np.zeros((len(ds_list),len(countries)))

        for i,ds in enumerate(ds_list):
            cntrymean = ds["countrymean"]
            cntry68 = ds["country68"]
            cntry95 = ds["country95"]
            cntryprior = ds["countryprior"]
            units = ds.countrymean.attrs['units']

            cntrymean_arr[i,:] = cntrymean
            cntry68_arr[i,:,:] = cntry68
            cntry95_arr[i,:,:] = cntry95
            cntryprior_arr[i,:] = cntryprior
        
    if plot_countryemissions == True and len(dates)>1:
        for ii,cntry in enumerate(countries_to_plot):
            cntry_ind = countries.index(cntry)
            if CI_to_plot==95:
                cntryCI_arr = cntry95_arr
            else:
                cntryCI_arr = cntry68_arr
            process.plot_country_timeseries(cntrymean_arr[:,cntry_ind], cntryCI_arr[:,cntry_ind,:],
                                            cntryprior_arr[:,cntry_ind], dates, country_label = cntry, units = units, figsize = (7,3))

