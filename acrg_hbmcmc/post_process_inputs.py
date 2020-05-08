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
import matplotlib.pyplot as plt

acrg_path=os.getenv('ACRG_PATH')

if __name__=="__main__":

    #### GENERAL INPUTS ####

    dates=["2014-01-01"] # Can be a list of one date or many dates
    species="chbr3"
    domain="CARIBBEAN"
    runname='test'
    # The following pertaining to domain_name will have to be edited once country_emissions
    # is re-coded in hbmcmc_post_process.py (ag12733 8/5/20)
    # Name that specifies the naming of the country file. This allows one to specify
    # a region that is smaller than sub-transd to compute totals over different areas
    # Example: domain_name = 'THD_coastal' corresponds to a country file called
    # country_THD_coastal.nc
    # country_THD_coastal can define a region smaller than sub-transd but MUST 
    # be on the same lat/lon grid!
    domain_name = domain     
    # country_dir is optional and allows you to specify a different directory where a custom
    #country file can live. Default of None points to the /data/shared/LPDM/countries/
    #Example: country_dir = '/data/ag12733/countries/'
    country_dir = None
    
    output_directory = "/work/ag12733/Python_outputs/Barbados_CHBr3/" # ** UPDATE OUTPUT DIRECTORY **
    
    #### POST-PROCESSING OPTIONS ####
    
    calc_country=False # At the moment calc_country must be True to write or append to nc file.

    plot_scale_map=False
    plot_abs_map=False
    plot_diff_map=False
    plot_y_timeseries=True
    
    #### PARAMETERS FOR EACH POST-PROCESSING OPTION ####
    # Only used if relevant value above is set to True
    
    # plot_scale_map
    grid_scale_map = True
    s_clevels = None # Set to None to set to defaults.
    s_cmap = plt.cm.RdBu_r
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
    d_cmap = plt.cm.YlGnBu
    d_smooth = False
    d_out_filename = None  # None means plot will not be written to file
    
    # plot_y_timeseries
    # combine timeseries will currently fail as needs to be re-written (ag12733 8/5/20)
    combine_timeseries = False # Plot y timeseries on one axis for multiple input files.
    y_out_filename = None
    
    
    #### IMPLEMENT PROCESSING OPTIONS ####
    print('Beginning post processing')
    
    if output_directory == "/path/to/output/directory/":
        raise Exception("Please set output directory.")
    if not os.path.isdir(output_directory):
        raise Exception("Output directory: {} does not exist.".format(output_directory))
    
    # Extract datasets from file
    ds_list,filenames = process.extract_hbmcmc_files(output_directory,species, domain, runname,dates,return_filenames=True)
    dates = process.check_missing_dates(filenames,dates)

#     NEEDS TO BE REWRITTEN FOR HBMCMC OUTPUT (ag12733 1/5/20)
#     ## Calculate country totals
#     if calc_country == True:
#         country_it,country_mean,country_percentile,country_prior,country_index \
#          = process.country_emissions_mult(ds_list, countries, species, domain_name, percentiles,
#                                 units=units, ocean=ocean, uk_split=uk_split, country_dir=country_dir)
#         country_data = {}
#         country_data["country_mean"] = country_mean
#         country_data["country_percentile"] = country_percentile
#         country_data["country_prior"] = country_prior
#         country_data["country_index"] = country_index
#         country_data["country"] = countries
#         # Setting a priori country values as zero for now
#         country_data["country_ap_percentile"] = np.zeros((len(countries),len(ds_list),len(percentiles)))
#     else:
#         country_data = None
    
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
 