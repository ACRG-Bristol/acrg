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
import os
import matplotlib.pyplot as plt

acrg_path=os.getenv('ACRG_PATH')

def country_emissions_multi(ds_list, countries, species, domain, x_post_vit=None, q_ap_abs_v=None, 
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
        country_out = process.country_emissions(ds, countries, species, domain, 
                                                percentiles=percentiles,units=units, 
                                                ocean=ocean, uk_split=uk_split)
        
        country_mean[:,tt] = country_out[1]
        country_percentile[:,tt,:] = country_out[2]
        
        if i == 0:
            # Should be the same for all datasets
            country_it = country_out[0]
            country_prior = country_out[3]
            country_index = country_out[4]

    return country_it,country_mean,country_percentile,country_prior,country_index

def create_output_param(ds_list,dates,experiment,nc_outfile,
                        country_dict=None,percentiles=[5,16,50,84,95],mode="write"):
    '''
    The create_output_param creates a dictionary containing the parameters to be passed to
    the write_outfile or append_outfile functions.
    '''
    
    param = {}
     
    param["lat"] = ds_list[0].lat.values
    param["lon"] = ds_list[0].lon.values
    
    param["time"] = pandas.to_datetime(dates)
    param["country"] = countries
    param["experiment"] = experiment
    param["outfile"] = nc_outfile
    param["percentile"] = percentiles
    
    nlat = len(param["lat"])
    nlon = len(param["lon"])
    nIt = len(ds_list[0].nIt)
    ntime = len(dates)
    npercentile = len(percentiles)
    
    param["nIt"] = nIt
    param["flux_it"] = np.zeros((nlat, nlon, ntime, nIt))
    param["flux_mean"] = np.zeros((nlat,nlon,ntime))
    param["flux_percentile"] = np.zeros((nlat,nlon,ntime,npercentile))
    
    param["k_mean"] = np.zeros(ntime)
    param["k_percentile"] = np.zeros((ntime,npercentile))
    param["flux_prior"] = np.zeros((nlat,nlon,ntime))
    
    for i,ds in enumerate(ds_list):
        
        param["flux_it"][:,:,i,:] = process.flux_iterations(ds)
        param["flux_mean"][:,:,i] = process.flux_mean(ds)
        param["flux_percentile"][:,:,i,:] = process.flux_percentile(ds,percentiles)
        param["flux_prior"][:,:,i] = ds.q_ap.values.copy()
        
        param["k_mean"][i],param["k_percentile"][i,:] = process.k_mean_percentile(ds,percentiles)

    # Define flux_ap_percentile as zeros for now
    param["flux_ap_percentile"] = np.zeros((nlat,nlon,ntime,npercentile))

    if country_dict:
        param.update(country_dict)

    if mode == "write":
        function = process.write_netcdf
    elif mode == "append":
        function = process.append_netcdf

    param_from_fn = function.__code__.co_varnames[:function.__code__.co_argcount]
    checked_param = {}
    for p in param_from_fn:
        if p in param:
            checked_param[p] = param[p]
    param = checked_param
    
    return param
    

if __name__=="__main__":

    #### GENERAL INPUTS ####

    dates=["2012-01-01"] # Can be a list of one date or many dates
    species="ch4"
    domain="SOUTHAMERICA"

    network="GOSAT-BRAZIL"
    experiment="MHD_TAC"
    countries=np.asarray(['UNITED KINGDOM', 'IRELAND', 'FRANCE', 'GERMANY', 
                          'DENMARK', 'BELGIUM', 'NETHERLANDS', 'LUXEMBOURG'])
    percentiles = [5,16,50,84,95]
    
    #output_directory = "/path/to/tdmcmc/outputs/" # ** UPDATE OUTPUT DIRECTORY **
    output_directory = "/home/rt17603/data/tdmcmc_output/SOUTHAMERICA_ch4/"
    
    #### POST-PROCESSING OPTIONS ####
    write_outfile=False
    append_outfile=True
    
    calc_country=True # At the moment calc_country must be True to write or append to nc file.

    plot_scale_map=False
    plot_abs_map=False
    plot_y_timeseries=False
    plot_regions=False
    plot_density=False
    
    #### PARAMETERS FOR EACH POST-PROCESSING OPTION ####
    # Only used if relevant value above is set to True
    
    # write_outfile / append_outfile
    nc_outfile = "flux_NAME-Bristol_ch4.nc"
    
    # calc_country
    units = 'Tg/yr'
    ocean = False
    uk_split = False
    
    # plot_scale_map
    s_clevels = np.arange(0.,2.1,0.1) # Set to None to set to defaults.
    s_cmap = plt.cm.RdBu_r
    s_smooth = True
    s_out_filename = None  # None means plot will not be written to file
    
    # plot_abs_map
    d_clevels = np.arange(-1.,1.05,0.05) # Set to None to set to defaults.
    d_cmap = plt.cm.RdBu_r
    d_smooth = False
    d_out_filename = None  # None means plot will not be written to file
    
    # plot_y_timeseries
    y_out_filename = None  # None means plot will not be written to file
    
    # plot regions
    r_out_filename = None  # None means plot will not be written to file
    
    # plot_density
    de_out_filename = None  # None means plot will not be written to file
    
    
    #### IMPLEMENT PROCESSING OPTIONS ####
    print 'Beginning post processing'
    
    #results = post_process(species, dates, network, output_directory, countries=countries,
    #                       write_outfile=False, append_outfile=False, calc_country=True,
    #                       plot_scale_map=True, plot_regions=True, plot_y_timeseries=True)
    
    if output_directory == "/path/to/output/directory/":
        raise Exception("Please set output directory.")
    if not os.path.isdir(output_directory):
        raise Exception("Output directory: {} does not exist.".format(output_directory))
    
    
    
    # Extract datasets from file
    ds_list = []
    for tt,date in enumerate(dates):
        fname_search = "output_{network}_{species}_{date}.nc".format(network=network,species=species,date=date)
        fname_search = os.path.join(output_directory,fname_search)
        filename = glob.glob(fname_search)
        if len(filename) > 0:
                ds = process.open_ds(filename[0])
                ds_list.append(ds)
    
    if not ds_list:
        raise Exception("No data found for dates: {}, species, {}, network {}".format(dates,species,network))
    
    ## Calculate country totals
    if calc_country == True:
        country_it,country_mean,country_percentile,country_prior,country_index \
         = country_emissions_multi(ds_list, countries, species, domain, percentiles,
                                units=units, ocean=ocean, uk_split=uk_split)
        country_data = {}
        country_data["country_mean"] = country_mean
        country_data["country_percentile"] = country_percentile
        country_data["country_prior"] = country_prior
        country_data["country_index"] = country_index
        country_data["country"] = countries
        # Setting a priori country values as zero for now
        country_data["country_ap_percentile"] = np.zeros((len(countries),len(ds_list),len(percentiles)))
    else:
        country_data = None
    
    ## Plot scaling map
    if plot_scale_map == True:
        for ds in ds_list:
            lon=np.asarray(ds.lon.values)
            lat=np.asarray(ds.lat.values)
            x_post_mean = process.x_post_mean(ds)
            
            stations = process.define_stations(ds)
            process.plot_map(x_post_mean,lon,lat,clevels=s_clevels, 
                                 cmap=s_cmap,label=None,
                                 smooth=s_smooth,fignum=None, 
                                 stations=stations,out_filename=s_out_filename)
    
    ## Plot absolute difference map
    if plot_abs_map == True:
        for ds in ds_list:
            lon=np.asarray(ds.lon.values)
            lat=np.asarray(ds.lat.values)
          
            x_post_mean = process.x_post_mean(ds)
            q_abs_diff = process.g2mol(process.flux_diff(ds),species)*1e6
            stations = process.define_stations(ds)
            
            process.plot_map(q_abs_diff,lon,lat,clevels=d_clevels, 
                                 cmap=d_cmap,label=None,
                                 smooth=d_smooth,fignum=None, 
                                 stations=stations,out_filename=d_out_filename)
     
    ## Plot y timeseries
    if plot_y_timeseries == True:
        for ds in ds_list:
            y_post_it,y_bg_it=process.plot_timeseries(ds, fig_text=None, 
                                                      ylim=None, out_filename=y_out_filename)

    ## Plot histogram of regions across iterations
    if plot_regions == True:
        for ds in ds_list:
            k_it=ds.k_it.values
            process.regions_histogram(k_it, fignum=None, out_filename=r_out_filename)
    
    ## Plot density of nuclei within trans-dimensional inversion
    if plot_density == True:
        for ds in ds_list:
            process.plot_nuclei_density(ds, out_filename=de_out_filename, fignum=None)
 
    if write_outfile == True or append_outfile == True:
        nc_outfile=os.path.join(output_directory,nc_outfile)
        
        if calc_country != True:
            raise Exception("At the moment, calc_country must be set to True to write or append to file.")
        
        if write_outfile:
            mode = "write"
            param = create_output_param(ds_list,dates,experiment,nc_outfile,
                                        country_dict=country_data,mode=mode,
                                        percentiles=percentiles)
            print "Writing flux to file: {}".format(nc_outfile)
            process.write_netcdf(**param)
        elif append_outfile:
            mode = "append"
            param = create_output_param(ds_list,dates,experiment,nc_outfile,
                            country_dict=country_data,mode=mode,
                            percentiles=percentiles)
            print "Appending flux to file: {}".format(nc_outfile)
            process.append_netcdf(**param)
