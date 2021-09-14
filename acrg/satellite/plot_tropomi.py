#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:48:12 2020

@author: rt17603
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from . import tropomi

def plot_tropomi_orbit(filename, parameter="methane_mixing_ratio_bias_corrected",
                       projection="Geostationary",central_longitude=0.0,
                       set_extent=None,set_nan_zero=True,
                       cmap="YlOrBr",
                       fig=None,subplot=[1,1,1]):
    '''
    Plot data from one tropomi oribit in its native resolution.
    
    Args:
        filename (str) :
            Full filename for tropomi data (including path).
        parameter (str, optional) :
            Parameter within tropomi data to plot (3D - time, scanline, ground_pixel only).
            i.e. not set up to plot parameters with additional dimensions e.g. averaging_kernel
            Default = "methane_mixing_ratio_bias_corrected"
        projection (str, optional) :
            Map projection to use when plotting. Currently limited to:
                - Geostationary (satellite view)
                - PlateCarree (equirectangular lat and lon, not good at the poles)
        central_longitude (float, optional) :
            Where to place central longitude when plotting.
            Default = 0.0
        set_extent (list/None, optional) :
            Only plot within given longitude and latitude bounds by passing a 4-item 
            list containing:
                [lon0, lon1, lat0, lat1]
            To not include this option pass a None value.
            Default = None
        set_nan_zero (bool, optional) :
            Set any NaN values within the data to 0 when plotting.
            Default = True
        cmap (str, optional) :
            Colour map to use when plotting the meshgrid. See matplotlib documentation
            for the full list of options.
            Default = "YlOrBr" (Yellow Orange Brown)
        fig (Figure object/None, optional) :
            Pass an explicit matplotlib Figure object to use when creating the plot
            Default = None          
        subplot (list, optional) 
            Only relevant if fig is specified. Allows plot to placed as a subplot
            of the figure object. This should a 3-item list containing 
            [nrows, ncols, number] e.g. [2, 1, 2]
            Default = [1, 1, 1] (one subplot)
    
    ## TODO: Add option to externally set min and max parameter values to be displayed
    # currently just uses the full value range.
    ## TODO: This has been tested on methane data files only so far.
    # Test this works with other tropomi data files.
    
    Returns:
        fig, ax - for plot
    '''
    
    
    data = tropomi.preProcessFile(filename)

    ## Create parameters associated with chosen projection
    if projection == "Geostationary":
        subplot_kw={"projection":ccrs.Geostationary(central_longitude=central_longitude)}
    elif projection == "PlateCarree":
        subplot_kw={"projection":ccrs.PlateCarree(central_longitude=central_longitude)}
    
    # Reduce from 3D to 2D by selecting on time axis
    # - FYI for each tropomi file time for the overall date (1 time value)
    # - FYI for full time values for each swath (ground_pixel) use "delta_time" data variable
    param = data[parameter].isel(time=0).values
    latitude_c = data["lat_corners"].isel(time=0)
    longitude_c = data["lon_corners"].isel(time=0)

    # Create plotting object
    if fig is None:
        fig, ax = plt.subplots(subplot_kw=subplot_kw)
    else:
        ax = fig.add_subplot(subplot[0],subplot[1],subplot[2],**subplot_kw)
    
    # Reduce lon-lat extent if specified
    if set_extent is not None:
        ax.set_extent(set_extent)
    
    # Plot as a colour mesh grid using longitude and latitude corners for each value
    ax.pcolormesh(longitude_c.values,latitude_c.values,param,
                  transform=ccrs.PlateCarree(),cmap=cmap)
    
    ax.coastlines() # Add coastlines to map

    return fig, ax       