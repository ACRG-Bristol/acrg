#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:13:25 2017

@author: as13988
"""

# Generates N HSV or RGB tuples spread over colour space 
class generatecolours:
    def __init__(self, N):
        
        import colorsys        
        
        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        
        self.RGB = RGB_tuples
        self.HSV = HSV_tuples
        

# Code to set up input for contour plotting
# data: either an xarray DataArray containing lat and lon dimensions or a class containing lat and lon 
# lon_range: manually set a lon range
# lat_range: manually set a lat range
# bottom_left: I don't know what the point of this option is!
# map_resolution: resolution of boundary database to use. Can be c (crude), l (low), i (intermediate), 
#      h (high), f (full) or None. If None, no boundary data will be read in. Resolution drops off by roughly 
#      80% between datasets. Higher res datasets are much slower to draw. Default l
    
class plot_map_setup:
    def __init__(self, data,
                 lon_range = None, lat_range = None,
                 bottom_left = False,
                 map_resolution = "l"):

        import xarray
        
        # Extract out the lats and lons
        # Check if it's an x-array 
        if isinstance(data, xarray.DataArry):
            lat = data.lat.values
            lon = data.lon.values
        else:
            lat = data.lat
            lon = data.lon

            
        if lon_range is None:
            lon_range = (min(lon), max(lon))
        if lat_range is None:
            lat_range = (min(lat), max(lat))
        
        m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')
        
        if bottom_left == False:
            lons, lats = np.meshgrid(lon,lat)
        else:
            dlon = lon[1] - lon[0]
            dlat = lat[1] -lat[0]            
            lons, lats = np.meshgrid(lon - dlon,
                                     lat - dlat)
        
        x, y = m(lons, lats)
        
        self.x = x
        self.y = y
        self.m = m
