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
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import xarray as xray
import acrg_name as name
from acrg_grid import areagrid
from netCDF4 import Dataset
from acrg_time.convert import time2sec
import os
import sys
import json
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
from matplotlib import ticker
import matplotlib.ticker as mticker
from cartopy.feature import BORDERS
from collections import OrderedDict
import datetime as dt
import getpass
import acrg_convert as convert
import pymc3 as pm
import matplotlib.dates as mdates
import pandas as pd

from acrg_config.paths import paths
acrg_path = paths.acrg

# Get acrg_site_info file
with open(os.path.join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f,object_pairs_hook=OrderedDict)

def check_platform(site,network=None):
    '''
    This function extracts platform (if specified) for the site from acrg_site_info.json file.
    network can be specified if site is associated with more than one. If not specified, the first
    network will be used by default.
    Returns:
        str : Platform type (e.g. "site", "satellite", "aircraft")
    '''
    site_info_file = os.path.join(acrg_path,"acrg_site_info.json")
    with open(site_info_file) as f:
        site_info=json.load(f,object_pairs_hook=OrderedDict)
    if network is None:
        network = list(site_info[site].keys())[0]
    if "platform" in site_info[site][network].keys():
        return site_info[site][network]["platform"]
    else:
        return None

def define_stations(ds,sites=None,use_site_info=False):
    '''
    The define_stations function defines the latitude and longitude values for each site within
    a dataset.
    
    If sites is not specified, if the platform for the site is listed as "aircraft" or 
    "satellite" in the acrg_site_info.json file then no values are included in the stations
    dictionary for this site.
    
    Args:
        ds (xarray.Dataset) :
            Output from run_tdmcmc() function (tdmcmc_inputs.py script).
            Expects dataset to contain:
                sitelons - Longitude values for each site. Dimension = len(sites)
                sitelats - Latitude values for each site. Dimension = len(sites)
                y_site       - Site identifier for each measurement. Dimension = nmeasure
        sites (list/None, optional) :
            List of sites to look for within dataset.
            If not specified, the sites will be extracted from the input dataset assuming a 
            data variable "sites" is included within the dataset.
        use_site_info (bool, optional) :
            Use positions from acrg_site_info.json file rather than extract them from the tdmcmc dataset.
            Default = False.
    
    Returns:
        dict :
            Dictionary containing sitelats, sitelons for each site.
    '''

    if sites is None:
        sites = list(ds.sitenames.values.astype(str))
        for site in sites:
            if check_platform(site) == "aircraft" or check_platform(site) == "satellite":
                sites.remove(site)
                
    stations={}
    
    if use_site_info:
        for site in sites:
            network = list(site_info[site].keys())[0]
            stations[site+'lons'] = [site_info[site][network]["latitude"]]
            stations[site+'lats'] = [site_info[site][network]["longitude"]]
    else:
        for site in sites:
            wh = np.where(ds.sitenames.values.astype(str) == site)[0]
            if len(wh) > 0:
                si = wh[0]
                #if site in ds.y_site:
                stations[site+'lons']=ds.sitelons[si].values
                stations[site+'lats']=ds.sitelats[si].values
            elif len(wh) == 0:
                print("WARNING: Reference to site not found within dataset")

    if sites:
        stations['sites']=sites
    else:
        stations = None
    
    return stations

def subplot_fmt(num,row_dims=[3,2,4],fill=False):
    '''
    The subplot_fmt function decides the placement of a grid of figures dependent on the number.
    The row_dims input determines which placement is preferable for the user.
    
    Args:
        num (int) :
            Number of figures to be placed
        row_dims (list, optional) : 
            Row dimensions in order of preference.
            For the default row_dims=[3,2,4] the preferences of placement is as follows:
                - equal rows of 3
                - equal rows of 2
                - equal rows of 4
            If none of the above are possible the format will be num x number of columns if fill 
            is True or the configuration suitable for num+1 if fill is False.
        fill (bool, optional) :
            All panels in subplot must be filled. If not, for uneven numbers an extra panel will
            be added which will be left blank when plotting.
            Default = False (i.e. allow an empty panel to be included within subplot)
                         
    Returns:
        List (int): [row_num,col_num]
                    2 item list containing the row number and column number for the subplots. 
    '''
    for r in row_dims:
        if not num%r:
            subplot = [r,num//r]
            break
    else:
        if fill or num == 1:
            subplot = [1,num]
        else:
            for r in row_dims:
                if not (num+1)%r:
                    subplot = [r,(num+1)//r]
                    break
    
    return subplot

def set_clevels(data,num_tick=20.,tick=None,centre_zero=False,above_zero=False,rescale=False,robust=False):
    '''
    The set_clevels function defines a set of contour levels for plotting based on the inputs 
    values.
    
    Args:
        data (iterable) :
            Data which will be plotted.
        num_tick (int) :
            Number of ticks on axis within levels.
            Either this or tick should be specified.
            Default = 20
        tick (int/None) :
            Tick interval to use between minimum and maximum data values.
            Either this or num_tick should be specified.
            Default = None i.e. use num_tick rather than set an explicit tick interval
        centre_zero (bool, optional) :
            Explictly centre levels around zero.
            Default = False.
        rescale (bool, optional) :
            Rescale according to the most appropriate unit.
            This will rescale based on 10^3 and return the scaling factor used.
            Default = False
        robust (bool, optional) :
            Based on xarray plotting. This finds the 2nd and 98th percentiles (rather
            than min and max) to account for any outliers which would cause the range to 
            be too large.
            Default = False
    
    Returns:
        np.array[,float] :
            Array of levels values based on min and max of input data.
            Also returns scaling factor if rescale=True.
    '''
    if robust:
        if above_zero:
            q_min = np.percentile(data[data>0],2)
            q_max = np.percentile(data,98)            
        else:
            q_min = np.percentile(data,2)
            q_max = np.percentile(data,98)
    else:
        if above_zero:
            q_min = np.min(data[data>0])
            q_max = np.max(data)            
        else:
            q_min = np.min(data)
            q_max = np.max(data)
    
    scale = 1

    if rescale:
        # Allow q to be rescaled according to the most appropriate unit
        while abs(q_max) <= 1e-3:
            q_max*=1e3
            q_min*=1e3
            scale*=1e-3

    if centre_zero:
        # If q_max and q_min are above and below zero, centre around zero.
        if q_min < 0 and q_max > 0:
            if abs(q_max) > abs(q_min):
                q_min = -1*q_max
            elif abs(q_max) < abs(q_min):
                q_max = -1*q_min
        else:
            print('Cannot centre on zero as min and max are not less than and greater than zero respectively')
        
    if not tick and num_tick:
        tick = (q_max-q_min)/num_tick
    elif not tick and not num_tick:
        raise Exception("Either tick or num_tick must be specified to define levels.")
    
    levels = np.arange(q_min,q_max,tick)
    
    if rescale:
        return levels,scale
    else:
        return levels

def unbiasedDivergingCmap(data, zero = 0, minValue = None, maxValue = None):
    """
    Calculate the normalisation of a diverging cbar around a given value
    Prevents bias due to asymetry in data affecting scale
    
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

def plot_map(data, lon, lat, clevels=None, divergeCentre = None, cmap=plt.cm.RdBu_r, borders=True,
             label=None, smooth = False, out_filename=None, stations=None, fignum=None,
                 title=None, extend="both", figsize=None, fig=None, ax=None, show=True):
    
    """
    Plot 2d map of data
    e.g. scaling map of posterior x i.e. degree of scaling applied to prior emissions 
    
    Args:
        data (numpy.array) : 
            2D (lat,lon) array of whatever you want
        lon (numpy.array) : 
            Longitude array matching to data grid
        lat (numpy.array) : 
            Latitude array  matching to data grid
        clevels (numpy.array, optional) : 
            Array of contour levels; defaults to np.arange(-2., 2.1, 0.1)
        divergeCentre (float/None, optional):
            Default is None, to replicate original clevels behaviour.
	        If given a float, this value is used to manually set the centre value of a diverging cmap,
	        while using the min and max values of clevels as the min and max values of the cmap.
        cmap (matplotlib.cm, optional) : 
            Colormap object; defaults to Red Blue reverse
        borders (bool, optional) :
            Add country borders as well as coastlines. Default = True.
        label (str, optional) : 
            Label to appear underneath the colorbar. Default = None
        smooth (bool, optional) : 
            If True plot smooth contours; otherwise use pcolormesh. Default = False.
        out_filename (str, optional) :
            Output filename. If this is specified the plot will be written to file. 
            Will be shown interactively otherwise (if show=True). Default = None.
        stations (dict, optional) : 
            Default is None. If specified needs to be a dictionary containing the list of sites
            and site locations for each site. For example:
                {"sites": ['MHD', 'TAC'],
                 MHD_lon: -9.02,
                 MHD_lat: 55.2,
                 TAC_lon: etc...
        fignum (int, optional) : 
            Figure number for created plot. Default = None
        title (str, optional) : 
            Title for the plot or sub-plot. Default = None
        extend (str, optional) :
            Extend colorbar for out-of-range values.
            Options are [ 'neither' | 'both' | 'min' | 'max' ]
            Default = "both". Set to "neither" to not extend.
        figsize (tuple/None, optional) :
            Figure size tuple if creating a new fig object.
            Default = None.
        fig (matplotlib.pyplot.Figure, optional) :
            Figure object for plot. If not specified this will be created. Default = None
        ax (matplotlib.pyplot.Axes, optional) :
            Axes object for plot domain. If not specified this will be created. Default = None
        show (bool, optional) :
            Whether to plot immediately upon completion of plotting within this function.
            Note that out_filename supercedes this option and plot will be written to file even
            if this is set to True.
            Default = True.
        
    Returns:
        None
        
        If out_filename is None:
            Created plot is saved to file
        Else:
            Plot is displayed interactively
    """
    if fig is None and ax is None:
        fig = plt.figure(fignum,figsize=figsize)
        ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=np.median(lon)))
    
    ax.set_extent([lon[0], lon[-1], lat[0], lat[-1]], crs=ccrs.PlateCarree())
    ax.coastlines()
    if borders:
        ax.add_feature(BORDERS,linewidth=0.5)
        
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=0, color='gray', linestyle='-')
    gl.xlabels_top = False
    gl.ylabels_right = False
#     gl.xlocator = mticker.FixedLocator([-66, -65, -64, -63])
#     gl.ylocator = mticker.FixedLocator([-15, -14, -13, -12])
    
    if clevels is None:
        #print "Warning: using default contour levels. Include clevels keyword to change"
        #clevels = np.arange(-2., 2.1, 0.1)
        print("Warning: using default contour levels which uses 2nd-98th percentile. Include clevels keyword to change.")
        clevels = set_clevels(data,robust=True)

        
    if smooth == True:
        if divergeCentre is None:
            cp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), cmap = cmap, 
                             levels=clevels, extend=extend)
        else:
            norm = unbiasedDivergingCmap(data, zero=divergeCentre, minValue=clevels[0], maxValue=clevels[-1])
            cp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), cmap = cmap, levels=clevels,
                    norm = norm, extend=extend)
        cb = plt.colorbar(cp, ax=ax, orientation='horizontal', pad=0.05)
    else:
        lons, lats = np.meshgrid(lon,lat)
        if divergeCentre is None:
            norm = BoundaryNorm(clevels,ncolors=cmap.N,clip=True)
        else:
            norm = unbiasedDivergingCmap(data, zero=divergeCentre, minValue=clevels[0], maxValue=clevels[-1])
        cs = ax.pcolormesh(lons, lats, data,cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        cb = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.1, extend=extend)
                
    if label is not None:        
        cb.set_label(label)#,fontsize=14) 
    if title is not None:
        #ax.set_title(title, fontsize=16) 
        fig.suptitle(title)#, fontsize=16)
        
    if stations is not None:
        for si,site in enumerate(stations['sites']):
            ilon=stations[site+'lons']
            ilat=stations[site+'lats']
            ax.plot(ilon, ilat, color = 'black', marker = 'o', markersize=8,  transform=ccrs.PlateCarree())
               
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()                 
    
    fig.tight_layout()
    if out_filename is not None:
        fig.savefig(out_filename)
        #plt.close(fig=fig)
    elif show:
        fig.show()
    
    return fig,ax

def plot_map_mult(data_all, lon, lat, grid=True, subplot="auto", clevels=None, divergeCentre=None, 
                 centre_zero=False,cmap=plt.cm.RdBu_r, borders=True, labels=None,
                 smooth=False, out_filename=None, stations=None, fignum=None,
                 title=None, extend="both",figsize=None):
    '''
    Uses plot_map function to plot a set of maps either on a grid or as separate figures.
    If plotting on a grid the subplots are either determined automatically based on shape of 
    input or using subplot input.
    
    Expect data_all to either be:
         - a numpy.array of the shape: nlat x nlon (x ngrid)
         - list of numpy.array objects each of shape nlat x nlon.
    Either the ngrid dimension or the len of the list is taken as the number of panels to 
    include on the plot.
    
    Args:
        data_all (numpy.array/list) :
            Multiple lat-lon grids to be plotted on one figure as a set of sub-plots or as 
            multiple figures.
            Can either be a list of grids or an array of dimension nlat x nlon (x ngrid).
        lon (numpy.array) :
            Longitude array matching to longitude points is each grid in grid_data.
        lat (numpy.array) :
            Latitude array matching to longitude points is each grid in grid_data.
        grid (bool, optional) :
            Whether to plot on a grid.
            Default = True.
        subplot (str/list, optional) :
            If grid is True, subplot grid to use. If this is set to "auto" this will be 
            automatically determined based on the size of ngrid (see subplot_fmt() function).
            Otherwise, this should be a two item list of [nrows, ncols]
            Default = "auto".
        labels (str/list, optional) :
            Can specify either one label for all plots (str) or a different label for 
            each plot as a list.
            If list is specified, it must match ngrid length.
        
        See plot_map() function for definition of remaining inputs.
    
    Returns:
        None
        
        If out_filename specified:
            Plot is written to file
        Otherwise:
            Plot is displayed interactively
    '''
    if isinstance(data_all,list):
        data_all = np.moveaxis(np.stack(data_all),0,2)
    elif isinstance(data_all,np.ndarray):
        if len(data_all.shape) == 2:
            data_all = np.expand_dims(data_all,2)
        elif len(data_all.shape) != 3:
            raise Exception("Did not understand input for data array to plot. Shape: {}".format(data_all.shape))
    
    nlat,nlon,nrun = data_all.shape
    if nlat != len(lat) or nlon != len(lon):
        raise Exception("First two dimensions of data_all ({},{}) must match length of lat ({}) and lon ({}) co-ordinates.".format(nlat,nlon,len(lat),len(lon)))
  
    if subplot == "auto" and grid:
        subplot = subplot_fmt(nrun)
    elif grid is False:
        subplot = [1,1]

    if isinstance(labels,list):
        if len(labels) == 1:
            labels = labels*nrun
        elif len(labels) != nrun:
            print("Unable to apply labels to sub-plots. Length of the list ({}) does not match the number of plots ({}).".format(len(labels),nrun))
            labels = [None]*nrun
    else:
        labels = [labels]*nrun
    
    if stations is None:
        stations = [None]*nrun
    else:
        if isinstance(stations,dict):
            stations = [stations]*nrun
        elif len(stations) != nrun:
            print("Unable to apply station positions to sub-plots. Number of station dictionaries ({}}) does not match the number of plots ({}).".format(len(stations),nrun))
            labels = [None]*nrun
    if not grid and nrun > 1:
        if out_filename:
            base,ext = os.path.splitext(out_filename)
            out_filename = []
        out_filename = None
        fignum = None
    
    if clevels is None:
        # Standarises clevels across all plots.
        clevels = set_clevels(data_all,centre_zero=centre_zero,robust=True)
    
    for i in range(nrun):
        data = data_all[...,i]
        
        if i == 0 and grid:
            fig = plt.figure(fignum,figsize=figsize)
            position = i+1
        elif grid is False:
            fig = plt.figure(figsize=figsize)
            position = 1
        else:
            position = i+1
        
        ax = fig.add_subplot(subplot[0],subplot[1],position,projection=ccrs.PlateCarree())

        if i < nrun-1 and grid:
            plot_map(data,lon,lat,clevels=clevels, divergeCentre = divergeCentre, 
                 cmap=cmap, borders=borders, label=labels[i], smooth=smooth, stations=stations[i],
                 title=None, extend=extend, out_filename=None, show=False, ax=ax, fig=fig)
        else:
            plot_map(data,lon,lat,clevels=clevels, divergeCentre = divergeCentre, 
                 cmap=cmap, borders=borders, label=labels[i], smooth=smooth, stations=stations[i],
                 title=title, extend=extend, out_filename=out_filename, show=True, ax=ax, fig=fig)

def plot_scale_map(ds_list, lat=None, lon=None, grid=True, clevels=None, divergeCentre=None, centre_zero=False,
                   cmap=plt.cm.YlGnBu, borders=True, labels=None, plot_stations=True,
                   use_site_info=False,
                   smooth=False, out_filename=None, fignum=None, title=None, extend="both",
                   figsize=None):
    '''
    The plot_scale_map function plots 2D scaling map(s) of posterior x. This is the degree of 
    scaling which has been applied to prior emissions.
    
    Args:
        ds_list (list) :
            List of xarray.Dataset objects. Each dataset is an output from run_tdmcmc() 
            function (tdmcmc_inputs.py script).
            Expects each data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
        lat (data array):
            Data array of lat values to plot over - must match values in ds exactly
        lon (data array):
            Data array of lon values to plot over - must match values in ds exactly
        grid (bool, optional) :
            Whether to plot the posterior on one figure as a grid or on individual plots.
        labels (str/list, optional) :
            Can specify either one label for all plots (str) or a different label for 
            each plot.
            If list is specified, it must match number of datasets in ds_list.
        plot_stations (bool, optional) :
            Plot site positions on the output map. Will not plot aircraft or satellite positions.
        use_site_info (bool, optional) :
            If plotting site positions, use positions from acrg_site_info.json file rather
            than extract them from the tdmcmc dataset.
            Default = False.
        
        See plot_map() function for definition of remaining inputs.
    
    Returns:
        None
        
        If out_filename specified:
            Plot is written to file
        Otherwise:
            Plot is displayed interactively
    '''

    if plot_stations:
        stations = [define_stations(ds,use_site_info=use_site_info) for ds in ds_list]
    else:
        stations = None
        
    if lat is None:
        lat = ds_list[0]["lat"]
    if lon is None:
        lon=ds_list[0]["lon"]
        
    ds_list = [ds.sel(lon=slice(np.min(lon.values),np.max(lon.values)),lat=slice(np.min(lat.values),np.max(lat.values))) for ds in ds_list]
    x_post_mean_list = [ds.scalingmean for ds in ds_list]
    
    plot_map_mult(x_post_mean_list, lon=lon, lat=lat, grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, centre_zero=centre_zero, 
                  cmap=cmap, labels=labels, smooth=smooth, out_filename=out_filename, 
                  stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return x_post_mean_list

def plot_abs_map(ds_list, species, lat=None, lon=None, grid=True, clevels=None, divergeCentre=None, 
                   cmap=plt.cm.YlGnBu, borders=True, labels=None, plot_stations=True,
                   use_site_info=False,
                   smooth=False, out_filename=None, fignum=None, title=None, extend="both",
                   figsize=None):
    '''
    The plot_abs_map function plots 2D map(s) of posterior x in g/m2/s.
    
    Args:
        ds_list (list) :
            List of xarray.Dataset objects. Each dataset is an output from run_tdmcmc() 
            function (tdmcmc_inputs.py script).
            Expects each data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
        lat (data array):
            Data array of lat values to plot over - must match values in ds exactly
        lon (data array):
            Data array of lon values to plot over - must match values in ds exactly
        species (str) :
            Species for the tdmcmc output.
        grid (bool, optional) :
            Whether to plot the posterior on one figure as a grid or on individual plots.
        labels (str/list) :
            Can specify either one label for all plots (str) or a different label for 
            each plot.
            If list is specified, it must match number of datasets in ds_list.
        plot_stations (bool, optional) :
            Plot site positions on the output map. Will not plot aircraft or satellite positions.
        use_site_info (bool, optional) :
            If plotting site positions, use positions from acrg_site_info.json file rather
            than extract them from the tdmcmc dataset.
            Default = False.
        
        See plot_map() function for definition of remaining inputs.
    
    Returns:
        None
        
        If out_filename specified:
            Plot is written to file
        Otherwise:
            Plot is displayed interactively
    '''

    if plot_stations:
        stations = [define_stations(ds,use_site_info=use_site_info) for ds in ds_list]
    else:
        stations = None
        
    if lat is None:
        lat = ds_list[0]["lat"]
    if lon is None:
        lon=ds_list[0]["lon"]
        
    ds_list = [ds.sel(lon=slice(np.min(lon.values),np.max(lon.values)),lat=slice(np.min(lat.values),np.max(lat.values))) for ds in ds_list]
    q_abs_list = [convert.mol2g(ds.fluxmean,species) for ds in ds_list]
    
        
    plot_map_mult(q_abs_list, lon=lon, lat=lat, grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, cmap=cmap, labels=labels, 
                  smooth=smooth, out_filename=out_filename, stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return q_abs_list

def plot_diff_map(ds_list, species, lat = None, lon = None, grid=True, clevels=None, divergeCentre=None, 
                   centre_zero=True,cmap=plt.cm.RdBu_r, borders=True, labels=None, plot_stations=True,
                   use_site_info=False,
                   smooth=False, out_filename=None, fignum=None, title=None, extend="both",
                   figsize=None):
    '''
    The plot_diff_map function plots 2D map(s) of the difference between the prior and 
    posterior x in g/m2/s.
    
    Args:
        ds_list (list) :
            List of xarray.Dataset objects. Each dataset is an output from run_tdmcmc() 
            function (tdmcmc_inputs.py script).
            Expects each data set to contain:
                x_post_vit - posterior values for each iteration flattened along lat-lon axis.
                             Dimensions = nIt x NGrid (nlat x nlon)
                q_ap       - a priori flux values on a latitude x longitude grid.
                             Dimensions = nlat x nlon
        lat (data array):
            Data array of lat values to plot over - must match values in ds exactly
        lon (data array):
            Data array of lon values to plot over - must match values in ds exactly
        species (str) :
            Species for the tdmcmc output.
        grid (bool, optional) :
            Whether to plot the posterior on one figure as a grid or on individual plots.
        labels (str/list) :
            Can specify either one label for all plots (str) or a different label for 
            each plot.
            If list is specified, it must match number of datasets in ds_list.
        plot_stations (bool, optional) :
            Plot site positions on the output map. Will not plot aircraft or satellite positions.
        use_site_info (bool, optional) :
            If plotting site positions, use positions from acrg_site_info.json file rather
            than extract them from the tdmcmc dataset.
            Default = False.
        
        See plot_map() function for definition of remaining inputs.
    
    Returns:
        None
        
        If out_filename specified:
            Plot is written to file
        Otherwise:
            Plot is displayed interactively
    '''

    if plot_stations:
        stations = [define_stations(ds,use_site_info=use_site_info) for ds in ds_list]
    else:
        stations = None

    if lat is None:
        lat = ds_list[0]["lat"]
    if lon is None:
        lon=ds_list[0]["lon"]
        
    ds_list = [ds.sel(lon=slice(np.min(lon.values),np.max(lon.values)),lat=slice(np.min(lat.values),np.max(lat.values))) for ds in ds_list]        
    q_diff_list = [convert.mol2g((ds.fluxmean - ds.fluxapriori),species) for ds in ds_list]
        
    plot_map_mult(q_diff_list, lon=lon, lat=lat, grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, centre_zero=centre_zero,
                  cmap=cmap, labels=labels, smooth=smooth, out_filename=out_filename, stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return q_diff_list

def country_emissions(ds, species, domain, country_file=None, country_unit_prefix=None, countries = None):

    c_object = name.get_country(domain, country_file=country_file)
    cntryds = xray.Dataset({'country': (['lat','lon'], c_object.country), 
                    'name' : (['ncountries'],c_object.name) },
                                    coords = {'lat': (c_object.lat),
                                    'lon': (c_object.lon)})
    # this allows the mcmc output to be sliced to match the size of a smaller country file
    # for example, if you want to look at a particular part of the ocean
    # the country file has to be rectangular (no odd shapes), however the contents inside can have any shape
    # (i.e. does not have to be the same size as the domain, but has to have the grid cells line up)

    lonmin_cds = np.min(cntryds.lon.values)
    lonmax_cds = np.max(cntryds.lon.values)
    latmin_cds = np.min(cntryds.lat.values)
    latmax_cds = np.max(cntryds.lat.values)

    # This step is here because of floating point differences in the two datasets
    # Convert the lat/lon bounds from the country dataset to the values in the ds
    lon_ds = ds.lon.values
    lat_ds = ds.lat.values
    lonmin_ds = lon_ds[np.where(np.isclose(lon_ds, lonmin_cds, atol = 0.01, rtol=0))[0][0]]
    lonmax_ds = lon_ds[np.where(np.isclose(lon_ds, lonmax_cds, atol = 0.01, rtol=0))[0][0]]
    latmin_ds = lat_ds[np.where(np.isclose(lat_ds, latmin_cds, atol = 0.01, rtol=0))[0][0]]
    latmax_ds = lat_ds[np.where(np.isclose(lat_ds, latmax_cds, atol = 0.01, rtol=0))[0][0]]
       
    ds = ds.sel(lon=slice(lonmin_ds,lonmax_ds),lat=slice(latmin_ds,latmax_ds))

    lon=ds["lon"].values
    lat=ds["lat"].values
    aprioriflux=ds["fluxapriori"]
    outs=ds["xtrace"]
    bfarray = ds["basisfunctions"]
    nui = ds["UInum"]
    steps = ds["stepnum"]
        
    area = areagrid(lat, lon)
    
    if countries is None:
        cntrynames = cntryds.name.values
    else:
        cntrynames = countries
    cntrygrid = cntryds.country.values
    cntrymean = np.zeros((len(cntrynames)))
    cntry68 = np.zeros((len(cntrynames), len(nui)))
    cntry95 = np.zeros((len(cntrynames), len(nui)))
    cntrysd = np.zeros(len(cntrynames))
    cntryprior = np.zeros(len(cntrynames))
    molarmass = convert.molar_mass(species)

    unit_factor = convert.prefix(country_unit_prefix)
    if country_unit_prefix is None:
        country_unit_prefix=''
    country_units = country_unit_prefix + 'g'

    for ci, cntry in enumerate(cntrynames):
        cntrytottrace = np.zeros(len(steps))
        cntrytotprior = 0
        for bf in range(int(np.max(bfarray))+1):
            bothinds = np.logical_and(cntrygrid == ci, bfarray.values==bf)
            cntrytottrace += np.sum(area[bothinds].ravel()*aprioriflux.values[bothinds].ravel()* \
                           3600*24*365*molarmass)*outs[:,bf]/unit_factor
            cntrytotprior += np.sum(area[bothinds].ravel()*aprioriflux.values[bothinds].ravel()* \
                   3600*24*365*molarmass)/unit_factor
        cntrymean[ci] = np.mean(cntrytottrace)
        cntry68[ci, :] = pm.stats.hpd(np.expand_dims(cntrytottrace,axis=1), 0.68)
        cntry95[ci, :] = pm.stats.hpd(np.expand_dims(cntrytottrace,axis=1), 0.95)
        cntryprior[ci] = cntrytotprior 
        
    return cntrymean, cntry68, cntry95, cntryprior
    
def country_emissions_mult(ds_list, species, domain, country_file=None, country_unit_prefix=None, countries = None):
    '''
    Calculate country emissions across multiple datasets.
    See process.country_emissions() function for details of inputs
    Returns:
        cntry_mean_list: list of mean country totals
        cntry_68_list: list of 68 CI country totals
        cntry_95_list: list of 95 CI country totals
    '''

    cntrymean_arr = np.zeros((len(ds_list), len(countries)))
    cntry68_arr = np.zeros((len(ds_list),len(countries),2))
    cntry95_arr = np.zeros((len(ds_list),len(countries),2))
    cntryprior_arr = np.zeros((len(ds_list),len(countries)))
    
    for i,ds in enumerate(ds_list):
        cntrymean, cntry68, cntry95, cntryprior = country_emissions(ds, species, domain, \
                                                        country_file=country_file, country_unit_prefix=country_unit_prefix, countries = countries)

        cntrymean_arr[i,:] = cntrymean
        cntry68_arr[i,:,:] = cntry68
        cntry95_arr[i,:,:] = cntry95
        cntryprior_arr[i,:] = cntryprior
        
    return cntrymean_arr, cntry68_arr, cntry95_arr, cntryprior_arr

def plot_country_timeseries(country_mean, country_CI, country_prior, d0, prior_label='Prior', posterior_label='Posterior', y_label='emissions', country_label = None, units = None, figsize=(7,3)):
    
    """
    Plot  timeseries of country emissions
    """
    fig,ax=plt.subplots(figsize=figsize)
    
    d0 = pd.to_datetime(d0)
    
    ax.plot(d0,country_prior, label=prior_label, color=(0,0.42,0.64))
    ax.plot(d0,country_mean, label=posterior_label, color = (0.78,0.32,0))
 
    ax.fill_between(d0, country_CI[:,0],country_CI[:,1], 
        facecolor=(0.78,0.32,0), edgecolor=(1,1,1), alpha=0.3) 

      
    ax.set_ylabel((country_label + " " + y_label + ", " + units + " yr$^{-1}$"), fontsize=10, fontweight = 'bold')
    ax.set_xlabel('Date', fontsize=10, fontweight = 'bold')
    legend=ax.legend(loc='upper left', labelspacing=0.1, fontsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))


# def combine_timeseries(*ds_mult):
#     '''
#     The combine_timeseries function takes a list of output datasets from a tdmcmc run and combines
#     the parameters relevant to plotting a mole fraction timeseries.

#     Current parameters copied from input and combined:
#         "y_time","y_site","y","sigma_y_it","sites"
    
#     Posterior boundary conditions (inner) and modelled mole fractions are calculated for each run and
#     combined. Both mean and full iterations are included. Additional parameters within dataset:
#         "bc_inner_post","bc_inner_post_it","mf_post","mf_post_it"
    
#     If any 
    
#     Args:
#         ds_mult (xarray.Dataset, xarray.Dataset, ...) :
#             Any number of tdmcmc output datasets can be specified to be combined.
#             All datasets will be combined.
    
#     Returns:
#         xarray.Dataset :
#             Reduced, combined dataset from tdmcmc output.
#     '''
#     calc_data_vars = {"bc_inner_post":post_bc_inner,"mf_post":post_mf,
#                       "bc_inner_prior":prior_bc_inner,"mf_prior":prior_mf}
#     #data_vars = ["y_time","y_site","y","sigma_y_it","bc_inner_post","mf_post","nIC","h_v_all","x_it","x_post_vit"]
#     data_vars = ["y_time","y_site","y","sigma_y_it","bc_inner_post","mf_post",
#                  "bc_inner_prior","mf_prior"]

# #    if prior:
# #        calc_data_vars["mf_prior"] = prior_mf
# #        data_vars.append("mf_prior")
# #    if bc_prior:
# #        calc_data_vars["bc_inner_prior"] = prior_bc_inner
# #        data_vars.append("bc_inner_prior")

#     match_coords = ["Ngrid","sites"]
    
#     concat_dim = "nmeasure"
#     iter_dim = "nIt"
#     run_dim = "run"
    
#     ## Check any coordinates within match_coords are the same for all the datasets.    
#     match_dim = [[ds.dims[coord] for ds in ds_mult] for coord in match_coords]
#     for ds in ds_mult:
#         for i,dim in enumerate(match_dim):
#             if len(set(dim)) != 1:
#                 raise Exception("Dimensions of {} do not match between input datasets.".format(match_coords[i]))
    
#     ## Extract and combine data arrays for data_vars including calculating any new data variables to add
#     data_arrays = OrderedDict()
    
#     for dv in data_vars:
#         for i,ds in enumerate(ds_mult):
#             if dv in calc_data_vars:
#                 if "post" in dv: # Posterior data variables
#                     calc_dv_it,calc_dv = calc_data_vars[dv](ds)
#                     if i == 0:
#                         da = xray.DataArray(calc_dv,dims=concat_dim)
#                         da_it = xray.DataArray(calc_dv_it,dims=[iter_dim,concat_dim])
#                     else:
#                         da = xray.concat([da,xray.DataArray(calc_dv,dims=concat_dim)],dim=concat_dim)
#                         da_it = xray.concat([da_it,xray.DataArray(calc_dv_it,dims=[iter_dim,concat_dim])],dim=concat_dim)
#                 else:
#                     calc_dv = calc_data_vars[dv](ds)
#                     if i == 0:
#                         da = xray.DataArray(calc_dv,dims=concat_dim)
#                     else:
#                         da = xray.concat([da,xray.DataArray(calc_dv,dims=concat_dim)],dim=concat_dim)
#             else:
#                 if i == 0:
#                     da = ds[dv]
#                 else:
#                     if concat_dim in da.dims:
#                         da = xray.concat([da,ds[dv]],dim=concat_dim)
#                     else:
#                         da = xray.concat([da,ds[dv]],dim=run_dim)
#         data_arrays[dv] = da
        
#         # Add data arrays containing the full iteration history as well as the mean for bg and x_post
#         if dv in calc_data_vars and "post" in dv:
#             data_arrays["{}_it".format(dv)] = da_it
    
#     ## Create dataset from extracted data arrays
#     combined_ds = xray.Dataset()
    
#     for d in data_arrays.items():
#         d = {d[0]:d[1]}
#         combined_ds = combined_ds.assign(**d)
    
#     # Add any coordinates from match_coords which are missing e.g. sites
#     for coord_name in match_coords:
#         if coord_name not in combined_ds.dims:
#             combined_ds = combined_ds.assign_coords(**{coord_name:ds_mult[0].coords[coord_name]})
    
#     # Add details the nmeasure value for each run (may be useful if we need to extract details for certain runs)
#     run = [ds.attrs["Start date"] for ds in ds_mult]
#     run_nmeasure = xray.DataArray([len(ds.nmeasure) for ds in ds_mult],coords={run_dim:run},dims=run_dim)
    
#     combined_ds = combined_ds.assign(**{"run_nmeasure":run_nmeasure})
#     combined_ds.attrs["Created by"] = getpass.getuser()
#     combined_ds.attrs["File created"] = str(dt.datetime.now().replace(microsecond=0))

#     return combined_ds
 
def plot_timeseries(ds, fig_text=None, ylim=None, out_filename=None, doplot=True, figsize=None,
                    plot_prior=False, plot_bc_prior=False,):
    
    """
    Plot measurement timeseries of posterior and observed measurements
    Requires post_mcmc xray dataset
    For future: incorporate model & measurement uncertainty
    Plots separate subplots for each of the measurement sites - hopefully!
    
    Args:
        ds (xarray dataset) : dataset output from run_tdmcmc script
        fig_text (String) : e.g. "CH$_{4}$ mole fraction (ppb)"
        ylim (array) : y-axis limits [ymin,ymax]
        out_filename (string) : Filename to save file
        doplot (bool) : Plot to console? (optional)
        figsize (tuple) : Specify size of figure as a two-item tuple.
        plot_prior (bool) : Plot mole fraction prior.
        plot_bc_prior (bool) : Plot inner boundary conditions prior.
        lower_percentile (float) : Lower percentile of predicted time series 
                                   uncertainty (default = 16)
        upper_percentile (float) : Upper percentile of predicted time series 
                                   uncertainty (default = 84)
        
    
    Specify an out_filename to write to disk
    """
    

    y_bg_mean = ds["YmodmeanBC"].values

    y_post_mean = ds["Ymodmean"].values
    
    if plot_prior:
        y_prior = ds["Yapriori"].values
    if plot_bc_prior:
        y_bc_prior = ds["YaprioriBC"].values

    sitenames = ds["sitenames"].values
    nsites=len(sitenames)
    
    y_time = ds.Ytime.values
    y_site = ds.siteindicator.values
    y_obs = ds.Yobs.values
    upper = ds.Ymod95.values[:,1]
    lower = ds.Ymod95.values[:,0]
    

    if doplot is True:
        fig,ax=plt.subplots(nsites,sharex=True,figsize=figsize)
        
        if nsites > 1:
            for si,site in enumerate(sitenames):
                wh_site = np.where(y_site == np.where(sitenames == site)[0][0])

                y_time_site = y_time[wh_site[0]]
                y_bg_site = y_bg_mean[wh_site[0]]
                y_post_site = y_post_mean[wh_site[0]]
                upper_site = upper[wh_site[0]]
                lower_site = lower[wh_site[0]]
                
                ax[si].fill_between(y_time_site, upper_site,lower_site, alpha=0.6, 
                            facecolor='lightskyblue', edgecolor='lightskyblue')
                
                ax[si].plot(y_time[wh_site[0]],y_obs[wh_site[0]], 'ro', markersize=4, label='Observations')
                
                ax[si].plot(y_time_site,y_post_site, color='blue', label='Modelled observations')
                
                ax[si].plot(y_time_site,y_bg_site,color='black', label='Modelled bounday conditions')
                
                if plot_prior:
                    ax[si].plot(y_time[wh_site[0]],y_prior[wh_site[0]], color='green', label='Prior')
                
                if plot_bc_prior:
                    ax[si].plot(y_time[wh_site[0]],y_bc_prior[wh_site[0]], color='0.6', label='Prior boundary conditions')
                    
                
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
                     label='Modelled boundary conditions')

            if plot_prior:
                ax.plot(y_time,y_prior, color='green', label='Prior')
            
            if plot_bc_prior:
                ax.plot(y_time,y_bc_prior, color='0.6', label='Prior boundary conditions')

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
        
    return

def open_ds(path):
    
    """
    Function efficiently opens xray datasets.
    """
    # use a context manager, to ensure the file gets closed after use
    with xray.open_dataset(path) as ds:
        ds.load()
    return ds    

def extract_hbmcmc_files(directory,species,domain, runname,dates,return_filenames=False):
    '''
    Find hbmcmc output filenames based on naming convention:
        "directory"/"species"+"domain"+"runname"_"date".nc"
    Open as xarray.Dataset objects and return as a list.
    '''
    species = species.upper()
    domain = domain.upper()
    ds_list = []
    filenames = []
    for tt,date in enumerate(dates):
        fname_search = "{species}_{domain}_{runname}_{date}.nc".format(domain=domain,species=species,date=date, runname=runname)
        fname_search = os.path.join(directory,fname_search)
        filename = glob.glob(fname_search)
        if len(filename) > 0:
            ds = open_ds(filename[0])
            ds_list.append(ds)
            filenames.append(filename[0])
    
    if not ds_list:
        raise Exception("No data found for dates {}, species, {}, runname {}, domain".format(dates,species,runname, domain))
    
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

def calculate_DIC(ds, silence=False):
    """
    Calculates the Deviance information criterion (DIC) for an inversion.
    It does this using two different definitions:
    1) Spiegelhalter et al. (2002) https://doi.org/10.1111/1467-9868.00353
    2) Gelman et al. (2004) http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf
    
    The DIC is similar to metrics like AIC (Akaike Information Criterion) but 
    better suited to hierarchical models and MCMC.
    What this is useful for is for testing things such as whether increasing the number of 
    basis functions improves things, or is it just fitting to the noise and making uncertainty larger.
    The lower the DIC the better.
    It should be noted that this doesn't give a hard and fast description of the suitability of a 
    statistical model but is the most useful indicator I can find.
    
    Args:
        ds (xarray dataset)      : dataset output from hbmcmc 
        silence (bool, optional) : Set to True to not print to screen  
    Returns:
        float : DIC using Spiegelhalter et al. (2002) definition
        float : DIC using Gelman et al. (2004) definition
    """
    sitenames = ds.sitenames.values
    # Calculate the log-likelihood of the means of the parameters    
    sig_arr=np.empty(len(ds.nmeasure.values))
    sig_trace_arr = np.empty((len(ds.nmeasure.values), len(ds.steps)))

    for n,site in enumerate(sitenames):
        site_idx=np.where(ds['siteindicator']==np.where(sitenames==site)[0][0])[0]
        sig_arr[site_idx]=np.mean(ds.sigtrace.values[:,n,:])
        sig_trace_arr[site_idx,:] =  ds.sigtrace.values[:,n,:].T

    y = ds.Yobs.values  
    hx_all_post=np.dot(ds.xsensitivity.values,np.mean(ds.xtrace.values,axis=0))
    hbc_all_post=np.dot(ds.bcsensitivity.values,np.mean(ds.bctrace.values,axis=0))
    mubar = hx_all_post+hbc_all_post
    D_thetabar= -0.5*(np.sum(2*np.log(sig_arr)) + np.sum((y-mubar)**2 / sig_arr**2) + \
                      np.log(2*np.pi)*len(y))

    # Calculate the mean log-likelihood
    mu_trace = np.dot(ds.xsensitivity.values,ds.xtrace.values.T) + \
                          np.dot(ds.bcsensitivity.values,ds.bctrace.values.T)
    D_theta_trace = -0.5*(np.sum(2*np.log(sig_trace_arr), axis=0) + \
                          np.sum((np.expand_dims(y, axis=1)-mu_trace)**2 / sig_trace_arr**2, axis=0) + \
                          np.log(2*np.pi)*len(y)) 
    Dbar_theta = np.mean(D_theta_trace, axis=0)
    p_DIC = 2*(D_thetabar - Dbar_theta)
    p_DIC_alt = 2*np.var(D_theta_trace)

    
    # DIC using Spiegelhalter et al. (2002) definition
    DIC_1 = -2*D_thetabar + 2*p_DIC
    # DIC using Gelman et al. (2004) definition
    DIC_2 = -2*D_thetabar + 2*p_DIC_alt

    if not silence:
        print("DIC using Spiegelhalter et al. (2002) definition")
        print(DIC_1)
        print("DIC using Gelman et al. (2004) definition")
        print(DIC_2)
    
    return DIC_1, DIC_2