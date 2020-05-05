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
import acrg_obs
import json
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
from matplotlib import ticker
from cartopy.feature import BORDERS
from collections import OrderedDict
import datetime as dt
import getpass

acrg_path = os.getenv("ACRG_PATH")

# Get acrg_site_info file
with open(os.path.join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f,object_pairs_hook=OrderedDict)

def molar_mass(species):
    '''
    This function extracts the molar mass of a species from the acrg_species_info.json file.
    Returns:
        float : Molar mass of species
    '''
    species_info_file = os.path.join(acrg_path,"acrg_species_info.json")
    with open(species_info_file) as f:
            species_info=json.load(f)
    species_key = acrg_obs.read.synonyms(species, species_info)
    molmass = float(species_info[species_key]['mol_mass'])
    return molmass

def mol2g(value,species):
    ''' Convert a value in moles to grams '''
    molmass = molar_mass(species)
    return value*molmass

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
                site_lon - Longitude values for each site. Dimension = len(sites)
                site_lat - Latitude values for each site. Dimension = len(sites)
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
            Dictionary containing "site"_lat, "site"_lon values for each site.
    '''

    if sites is None:
        sites = list(ds.sitenames.values.astype(str))
#        sites = [site.astype("unicode") for site in sites]
        for site in sites:
            if check_platform(site) == "aircraft" or check_platform(site) == "satellite":
                sites.remove(site)
                
    stations={}
    
    if use_site_info:
        for site in sites:
            network = list(site_info[site].keys())[0]
            stations[site+'_lon'] = [site_info[site][network]["latitude"]]
            stations[site+'_lat'] = [site_info[site][network]["longitude"]]
    else:
        for site in sites:
            wh = np.where(ds.sitenames.values.astype(str) == site)[0]
            if len(wh) > 0:
                si = wh[0]
                #if site in ds.y_site:
                stations[site+'_lon']=ds.site_lon[si].values
                stations[site+'_lat']=ds.site_lat[si].values
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
    
    #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=np.median(lon)))
    ax.set_extent([lon[0], lon[-1], lat[0], lat[-1]], crs=ccrs.PlateCarree())
    ax.coastlines()
    if borders:
        ax.add_feature(BORDERS,linewidth=0.5)

    if clevels is None:
        #print "Warning: using default contour levels. Include clevels keyword to change"
        #clevels = np.arange(-2., 2.1, 0.1)
        print("Warning: using default contour levels which uses 2nd-98th percentile. Include clevels keyword to change.")
        clevels = set_clevels(data,robust=True)
    #else:
    #    clevels = np.arange(clevels[0], clevels[1], clevels[2])
        
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
        cb = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, extend=extend)
                
    if label is not None:        
        cb.set_label(label)#,fontsize=14) 
    if title is not None:
        #ax.set_title(title, fontsize=16) 
        fig.suptitle(title)#, fontsize=16)
        
    if stations is not None:
        for si,site in enumerate(stations['sites']):
            ilon=stations[site+'_lon']
            ilat=stations[site+'_lat']
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

def plot_scale_map(ds_list, grid=True, clevels=None, divergeCentre=None, centre_zero=True,
                   cmap=plt.cm.RdBu_r, borders=True, labels=None, plot_stations=True,
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
    x_post_mean_list = [ds.meanscaling for ds in ds_list]
    
    plot_map_mult(x_post_mean_list, lon=ds_list[0]["lon"], lat=ds_list[0]["lat"], grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, centre_zero=centre_zero, 
                  cmap=plt.cm.RdBu_r, labels=labels, smooth=smooth, out_filename=out_filename, 
                  stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return x_post_mean_list

def plot_abs_map(ds_list, species, grid=True, clevels=None, divergeCentre=None, 
                   cmap=plt.cm.RdBu_r, borders=True, labels=None, plot_stations=True,
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
    q_abs_list = [mol2g(ds.meanflux,species) for ds in ds_list]
    
    plot_map_mult(q_abs_list, lon=ds_list[0]["lon"], lat=ds_list[0]["lat"], grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, cmap=cmap, labels=labels, 
                  smooth=smooth, out_filename=out_filename, stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return q_abs_list

def plot_diff_map(ds_list, species, grid=True, clevels=None, divergeCentre=None, 
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
    
    q_diff_list = [mol2g((ds.meanflux - ds.aprioriflux),species) for ds in ds_list]
    
    plot_map_mult(q_diff_list, lon=ds_list[0]["lon"], lat=ds_list[0]["lat"], grid=grid,
                  clevels=clevels, divergeCentre=divergeCentre, centre_zero=centre_zero,
                  cmap=plt.cm.RdBu_r, labels=labels, smooth=smooth, out_filename=out_filename, stations=stations, fignum=fignum, 
                  title=title, extend=extend, figsize=figsize)
    return q_diff_list

    
# def country_emissions(ds_mcmc, countries, species, domain, x_post_vit=None, q_ap_abs_v=None, 
#                       percentiles=[5,16,50,84,95], units=None, ocean=True, 
#                       uk_split=False, fixed_map=False, country_dir = None):
        
#     """
#     Generates national totals for a given list of countries
#     Args: 
#         ds_mcmc (xarray.Dataset) :
#             Output from run_tdmcmc() function (tdmcmc_inputs.py script). (post_mcmc xarray dataset)
#         countries (list) :
#             List of countries (see data_path/countries files for names)
#             In general names are all capitalized.
#         species (str) :
#             CH4, CO2 etc.
#         domain (str) :
#             Domain name of interest.
#             e.g. EUROPE, SOUTHAMERICA etc.
#         x_post_vit (numpy.array/None, optional) :
#             Posterior values (Dims: nIt x NGrid) (flattened on lat and lon dimensions).
#             If not specified, will be extracted from ds_mcmc.
#         q_ap_abs_v (numpy.array/None, optional) :
#             Flux a priori values (Dims: NGrid) (flattened on lat and lon dimensions).
#             If not specified, will be extracted from ds_mcmc.
#         percentiles (list) :
#             List of percentiles to extract.
#             Default = [5,16,50,84,95]
#         units (str/None, optional) :
#             Need to specify to get something sensible out, otherwise returns in g/yr
#             Options are Pg/yr, Tg/yr,Gg/yr Mg/yr
#         ocean (bool, optional) :
#             Default is True, If False only include emissions from land surface
#         uk_split (bool, optional) :
#             Break UK into devolved administrations
#                 Country names then have to be ['Eng', 'Sco', 'Wales',
#                 NIre', 'IRELAND', 'BENELUX']. Horribly inconsistent I know!
#             Default = False
#         fixed_map (bool, optional) :
#             TODO: Add comment
    
#     Returns: 
#         (5 x numpy.array) :
#             Country totals for each iteration in units specified (ncountries x nIt),
#             Mean of country totals in units specified (ncountries), 
#             Percentiles for each country in units specified (ncountries x npercentiles),
#             Prior for each country in units specified (ncountries),
#             Country index map (nlat x nlon)
       
#     Output in Tg/yr
#     """
    
#     if x_post_vit is None:
#         x_post_vit = ds_mcmc.x_post_vit.values
#     if q_ap_abs_v is None:
#         q_ap = ds_mcmc.q_ap.values
#         q_ap_abs_v = np.ravel(q_ap)
    
#     if units == 'Tg/yr':
#         unit_factor=1.e12
#     elif units == 'Gg/yr': 
#         unit_factor=1.e9
#     elif units == 'Pg/yr': 
#         unit_factor=1.e15
#     elif units == 'Mg/yr': 
#         unit_factor=1.e6
#     else:
#         print('Undefined units: outputting in g/yr - let this be a lesson to define your units')
#         unit_factor=1.
    
#     molmass = molar_mass(species)
    
#     lon=ds_mcmc.lon.values
#     lat=ds_mcmc.lat.values

#     lonmin=np.min(ds_mcmc.lon.values)
#     lonmax=np.max(ds_mcmc.lon.values)
#     latmin=np.min(ds_mcmc.lat.values)
#     latmax=np.max(ds_mcmc.lat.values)
    
#     nlon=len(ds_mcmc.lon)
#     nlat=len(ds_mcmc.lat)    
    
#     # GET COUNTRY DATA
#     if uk_split == True:
#         c_object=name.get_country(domain, ocean=True, ukmo=True, uk_split=uk_split)
#     elif ocean == True:
#         c_object=name.get_country(domain, ocean=True, ukmo=True, uk_split=False)
#     else:
#         c_object=name.get_country(domain, ocean=False)
#     cds = xray.Dataset({'country': (['lat','lon'], c_object.country), 
#                         'name' : (['ncountries'],c_object.name) },
#                                         coords = {'lat': (c_object.lat),
#                                         'lon': (c_object.lon)})
#     country = cds.country.sel(lon=slice(lonmin,lonmax), 
#                                         lat=slice(latmin,latmax))
    
#     lonmin_cds = np.min(country.lon.values.astype('float32'))
#     lonmax_cds = np.max(country.lon.values.astype('float32'))
#     latmin_cds = np.min(country.lat.values.astype('float32'))
#     latmax_cds = np.max(country.lat.values.astype('float32'))

#     x_post_vit=ds_mcmc.x_post_vit.values
#     nIt=len(x_post_vit[:,0])
#     x_post_vit = x_post_vit[:,:, np.newaxis]
#     x_post_it = np.reshape(x_post_vit, (nIt, nlat,nlon))
    
#     x_post_it_ds = xray.Dataset({"x_post_it" : 
#                         (["nIt", "lat", "lon"], x_post_it)},
#                        coords = {"lon": lon,
#                                  "lat": lat})
        
#     ds_mcmc = ds_mcmc.sel(lon=slice(lonmin_cds,lonmax_cds), 
#                                         lat=slice(latmin_cds,latmax_cds))

#     x_post_it_ds = x_post_it_ds.sel(lon=slice(lonmin_cds,lonmax_cds), 
#                                         lat=slice(latmin_cds,latmax_cds))
    
#     x_post_it = x_post_it_ds["x_post_it"].values
    
#     q_ap=ds_mcmc.q_ap.values
#     q_ap_abs_v=np.ravel(q_ap)
    
#     area=areagrid(ds_mcmc.lat.values,ds_mcmc.lon.values)
#     area_v=np.ravel(area)
#     country_v=np.ravel(country)
    
#     ncountries=len(countries)
#     npercentiles = len(percentiles)

#     country_v_new=np.zeros((nlon*nlat))

#     if len(np.shape(q_ap_abs_v)) == 1:

#         q_country_it = np.zeros((ncountries, nIt))
#         q_country_mean=np.zeros((ncountries))
#         q_country_percentile = np.zeros((ncountries,npercentiles))
#         q_country_ap=np.zeros((ncountries))
#         q_country_mode=np.zeros((ncountries))
    
#         for ci,cc in enumerate(countries):
#             name_country = np.where(cds.name == cc)
#             c_index = np.where(country_v == name_country[0])
#             country_v_new[c_index[0]]=ci+1
#             for it in range(nIt):
#                 x_vit_temp=np.ravel(x_post_it[it,:,:])
#                 if fixed_map == False:
#                     q_country_it[ci,it]=np.sum(x_vit_temp[c_index[0]]*area_v[c_index[0]]
#                                     *q_ap_abs_v[c_index[0]]) 
#                 elif fixed_map == True:
#                     q_country_it[ci,it]=np.sum(x_vit_temp*area_v[c_index[0]]
#                                     *q_ap_abs_v[c_index[0]]) 
#             q_country_ap[ci] = np.sum(area_v[c_index[0]]*q_ap_abs_v[c_index[0]]
#                                 *365.*24.*3600.*molmass/unit_factor) # in Tg/yr
        
#             #dum=np.histogram(q_country_it[ci,:], bins=100)                
#             #q_country_mode[ci] = dum[1][dum[0]==np.max(dum[0])][0]
            
#             q_country_mean[ci]=np.mean(q_country_it[ci,:])*365.*24.*3600.*molmass/unit_factor # in Tg/yr
#             for pi,percentile in enumerate(percentiles):
#                 q_country_percentile[ci,pi] = np.percentile(q_country_it[ci,:],percentile)*365.*24.*3600.*molmass/unit_factor
        
#         country_index = np.reshape(country_v_new, (nlat,nlon))   

#         return q_country_it*365.*24.*3600.*molmass/unit_factor,\
#         q_country_mean, q_country_percentile, q_country_ap, country_index
        
#     elif len(np.shape(q_ap_abs_v)) == 2:

#         q_ap_abs_v = q_ap_abs_v.T
#         ntimes = np.shape(q_ap_abs_v)[1]        
        
#         q_country_it = np.zeros((ncountries, ntimes, nIt))
#         q_country_mean=np.zeros((ncountries, ntimes))

#         q_country_percentile=np.zeros((ncountries, ntimes,npercentiles))
#         q_country_ap=np.zeros((ncountries, ntimes))
    
#         for ci,cc in enumerate(countries):
#             name_country = np.where(cds.name == cc)
#             c_index = np.where(country_v == name_country[0])
#             country_v_new[c_index[0]]=ci+1
#             for it in range(nIt):
#                 x_vit_temp=np.ravel(x_post_it[it,:,:])
# #                x_vit_temp=x_post_vit[it,:]
#                 if fixed_map == False:
#                     q_country_it[ci,:,it]=np.sum(x_vit_temp[c_index[0],None]*area_v[c_index[0],None]
#                                     *q_ap_abs_v[c_index[0],:], axis = 0) 
#                 elif fixed_map == True:
#                     q_country_it[ci,:,it]=np.sum(x_vit_temp*area_v[c_index[0],None]
#                                     *q_ap_abs_v[c_index[0],:], axis = 0)
#             q_country_ap[ci,:] = np.sum(area_v[c_index[0],None]*q_ap_abs_v[c_index[0],:]
#                                 *365.*24.*3600.*molmass/unit_factor, axis = 0) # in Tg/yr
        
#             q_country_mean[ci,:]=np.mean(q_country_it[ci,:,:], axis =1)*365.*24.*3600.*molmass/unit_factor # in Tg/yr
#             for pi,percentile in enumerate(percentiles):
#                 q_country_percentile[ci,:,pi]=np.percentile(q_country_it[ci,:,:],percentile, axis =1)*365.*24.*3600.*molmass/unit_factor

        

#         country_index = np.reshape(country_v_new, (nlat,nlon))   

#         return q_country_it*365.*24.*3600.*molmass/unit_factor,\
#         q_country_mean, q_country_percentile, q_country_ap, country_index

# def country_emissions_mult(ds_list, countries, species, domain, x_post_vit=None, q_ap_abs_v=None, 
#                       percentiles=[5,16,50,84,95], units=None, ocean=True, 
#                       uk_split=False, fixed_map=False, country_dir=None):
#     '''
#     Calculate country emissions across multiple datasets. Combine mean and percentiles into 
#     arrays for all time points.
#     See process.country_emissions() function for details of inputs
#     Returns:
#         (5 x numpy.array) :
#             Country totals for each iteration in units specified (ncountries x nIt),
#             Mean of country totals for each dataset in units specified (ncountries x ntime), 
#             Percentiles for each country for each dataset in units specified (ncountries x npercentiles x ntime),
#             Prior for each country in units specified (ncountries),
#             Country index map (nlat x nlon)
#     '''
#     ncountries=len(countries)
#     ntime = len(ds_list)
#     npercentile=len(percentiles)
    
#     # Constructed from all datasets
#     country_mean = np.zeros((ncountries,ntime)) 
#     country_percentile = np.zeros((ncountries,ntime,npercentile))
#     country_prior = np.zeros((ncountries,ntime))

#     for i,ds in enumerate(ds_list):
#         country_out = country_emissions(ds, countries, species, domain, 
#                                         percentiles=percentiles,units=units, 
#                                         ocean=ocean, uk_split=uk_split)
        
#         country_mean[:,i] = country_out[1]
#         country_percentile[:,i,:] = country_out[2]
#         country_prior[:,i] = country_out[3]

#         if i == 0:
#             # Should be the same for all datasets
#             country_it = country_out[0]
#             #country_prior = country_out[3]
#             country_index = country_out[4]

#     return country_it,country_mean,country_percentile,country_prior,country_index

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
    

    y_bg_mean = ds["YmodBC"].values

    y_post_mean = ds["Ymod"].values
    
    if plot_prior:
        y_prior = ds["Yapriori"].values
    if plot_bc_prior:
        y_bc_prior = ds["YaprioriBC"].values


    sites = ds["sitenames"].values
    nsites=len(sites)
    
#     sigma_y_mean=ds.Ymod.values

    y_time = ds.Ytime.values
    y_site = ds.siteindicator.values
    y_obs = ds.Y.values
    upper = ds.Ymod95.values[:,1]
    lower = ds.Ymod95.values[:,0]
    

    if doplot is True:
        fig,ax=plt.subplots(nsites,sharex=True,figsize=figsize)
        
        if nsites > 1:
            for si,site in enumerate(sites):
                wh_site = np.where(y_site == np.where(sitenames == site)[0][0])

                y_time_site = y_time[wh_site[0]]
                y_bg_site = y_bg_mean[wh_site[0]]
                y_post_site = y_post_mean[wh_site[0]]
                upper_site = upper[wh_site[0]]
                lower_site = lower[wh_site[0]]
                
                ax[si].fill_between(y_time_site, upper_site,lower_site, alpha=0.6, 
                            facecolor='lightskyblue', edgecolor='lightskyblue')
                ax[si].plot(y_time[wh_site[0]],y_obs[wh_site[0]], 'ro', markersize=4, label='Observations')
                
                #ax[si].plot(y_time[wh_site[0]],y_post_mean[wh_site[0]], color='blue', label='Modelled observations')
                #ax[si].plot(y_time[wh_site[0]],y_bg_mean[wh_site[0]],color='black', 
                #         label='Modelled baseline')
                ax[si].plot(y_time_site,y_post_site, color='blue', label='Modelled observations')
                ax[si].plot(y_time_site,y_bg_site,color='black', label='Modelled baseline')
                
                if plot_prior:
                    ax[si].plot(y_time[wh_site[0]],y_prior[wh_site[0]], color='green', label='Prior')
                
                if plot_bc_prior:
                    ax[si].plot(y_time[wh_site[0]],y_bc_prior[wh_site[0]], color='0.6', label='Prior baseline (bc inner)')
                    
                
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

            if plot_prior:
                ax.plot(y_time,y_prior, color='green', label='Prior')
            
            if plot_bc_prior:
                ax.plot(y_time,y_bc_prior, color='0.6', label='Prior baseline (bc inner)')

            start, end = ax.get_ylim()
            ax.yaxis.set_ticks(np.arange(start, end+1, (end-start)/5))
            #ax.tick_params(axis='both',labelsize=20)
            #legend=ax.legend(loc='upper left',fontsize=20)
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



