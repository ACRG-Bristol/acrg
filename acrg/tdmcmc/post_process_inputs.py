# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:24:46 2015

Template file for creating plots with output of tdmcmc

Uses tdmcmc_post_process

@author: ml12574 (updated by rt17603)
"""
import numpy as np
from . import tdmcmc_post_process as process
import glob
import pandas
import os
import matplotlib.pyplot as plt

from acrg.config.paths import Paths

acrg_path = Paths.acrg

def extract_tdmcmc_files(directory,species,network,dates,return_filenames=False):
    '''
    Find tdmcmc output filenames based on naming convention:
        "directory"/output_"network"_"species"_"date".nc e.g. output_AGAGE_ch4_2008-01-01.nc
    Open as xarray.Dataset objects and return as a list.
    '''
    ds_list = []
    filenames = []
    for tt,date in enumerate(dates):
        fname_search = "output_{network}_{species}_{date}.nc".format(network=network,species=species,date=date)
        fname_search = os.path.join(directory,fname_search)
        filename = glob.glob(fname_search)
        if len(filename) > 0:
            ds = process.open_ds(filename[0])
            ds_list.append(ds)
            filenames.append(filename[0])
    
    if not ds_list:
        raise Exception("No data found for dates: {}, species, {}, network {}".format(dates,species,network))
    
    if return_filenames:
        return ds_list,filenames
    else:
        return ds_list

def check_missing_dates(filenames,dates,labels=[]):
    
    if len(filenames) != len(dates):
        no_data = []
        for date in dates:
            for fname in filenames:
                if date in fname:
                    break
            else:
                no_data.append(date)
        
        new_dates = dates[:]
        for date in no_data:
            new_dates.remove(date)
        if labels == dates:
            labels = new_dates
        dates = new_dates
    
    if labels:
        return dates,labels
    else:
        return dates

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
    domain="EUROPE"
    # Name that specifies the naming of the country file. This allows one to specify
    # a region that is smaller than sub-transd to compute totals over different areas
    # Example: domain_name = 'THD_coastal' corresponds to a country file called
    # country_THD_coastal.nc
    # country_THD_coastal can define a region smaller than sub-transd but MUST 
    # be on the same lat/lon grid!
    domain_name = domain     
    # country_dir allows you to specify a different directory where a custom
    #country file can live. Default of None points to the /data/shared/LPDM/countries/
    #Example: country_dir = '/data/ag12733/countries/'
    country_dir = None

    network="AGAGE"
    experiment="MHD_TAC"
    countries=np.asarray(['UNITED KINGDOM', 'IRELAND', 'FRANCE', 'GERMANY', 
                          'DENMARK', 'BELGIUM', 'NETHERLANDS', 'LUXEMBOURG'])
    percentiles = [5,16,50,84,95]
    
    output_directory = "/path/to/output/directory" # ** UPDATE OUTPUT DIRECTORY **
    
    #### POST-PROCESSING OPTIONS ####
    write_outfile=False
    append_outfile=False
    
    calc_country=False # At the moment calc_country must be True to write or append to nc file.

    plot_scale_map=False
    plot_abs_map=False
    plot_diff_map=False
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
    grid_scale_map = True
    s_clevels = None # Set to None to set to defaults.
    s_cmap = plt.cm.RdBu_r
    s_smooth = True
    s_out_filename = None  # None means plot will not be written to file

    # plot_abs_map
    grid_abs_map = True
    a_clevels = None  # Set to None to set to defaults.
    a_cmap = plt.cm.RdBu_r
    a_smooth = False
    a_out_filename = None  # None means plot will not be written to file
    
    # plot_diff_map
    grid_diff_map = True
    d_clevels = None # Set to None to set to defaults.
    d_cmap = plt.cm.RdBu_r
    d_smooth = False
    d_out_filename = None  # None means plot will not be written to file
    
    # plot_y_timeseries
    combine_timeseries = True # Plot y timeseries on one axis for multiple input files.
    y_out_filename = None
    
    # plot regions
    r_out_filename = None  # None means plot will not be written to file
    
    # plot_density
    de_out_filename = None  # None means plot will not be written to file
    
    
    #### IMPLEMENT PROCESSING OPTIONS ####
    print('Beginning post processing')
    
    if output_directory == "/path/to/output/directory/":
        raise Exception("Please set output directory.")
    if not os.path.isdir(output_directory):
        raise Exception("Output directory: {} does not exist.".format(output_directory))
    
    # Extract datasets from file
    ds_list,filenames = extract_tdmcmc_files(output_directory,species,network,dates,return_filenames=True)
    dates = check_missing_dates(filenames,dates)
    
    ## Calculate country totals
    if calc_country == True:
        country_it,country_mean,country_percentile,country_prior,country_index \
         = process.country_emissions_mult(ds_list, countries, species, domain_name, percentiles,
                                units=units, ocean=ocean, uk_split=uk_split, country_dir = country_dir)
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
        process.plot_scale_map(ds_list,grid=grid_scale_map,clevels=s_clevels, 
                               cmap=s_cmap,labels=None,title=None,smooth=s_smooth,
                               out_filename=s_out_filename,extend="both")
    
    ## Plot absolute difference map
    if plot_diff_map == True:
        process.plot_diff_map(ds_list,species,grid=grid_diff_map,clevels=d_clevels, 
                               cmap=d_cmap,labels=None,title=None,smooth=d_smooth,
                               out_filename=d_out_filename,extend="both")

    ## Plot absolute map
    if plot_abs_map == True:
        process.plot_abs_map(ds_list,species,grid=grid_abs_map,clevels=a_clevels, 
                               cmap=a_cmap,labels=None,title=None,smooth=a_smooth,
                               out_filename=a_out_filename,extend="max")

     
    ## Plot y timeseries
    if plot_y_timeseries == True:
        if combine_timeseries and len(ds_list) > 1:
            ds_combined = process.combine_timeseries(*ds_list)
            y_post_it,y_bg_it=process.plot_timeseries(ds_combined, fig_text=None, 
                                                          ylim=None, out_filename=y_out_filename)
        else:
            for i,ds in enumerate(ds_list):
                if y_out_filename:
                    stub,ext = os.path.splitext(y_out_filename)
                    y_out_filename_n = "{}_{}{}".format(stub,i+1,ext)
                else:
                    y_out_filename_n = None
                y_post_it,y_bg_it=process.plot_timeseries(ds, fig_text=None, 
                                                          ylim=None, out_filename=y_out_filename_n)

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
            print("Writing flux to file: {}".format(nc_outfile))
            process.write_netcdf(**param)
        elif append_outfile:
            mode = "append"
            param = create_output_param(ds_list,dates,experiment,nc_outfile,
                            country_dict=country_data,mode=mode,
                            percentiles=percentiles)
            print("Appending flux to file: {}".format(nc_outfile))
            process.append_netcdf(**param)
