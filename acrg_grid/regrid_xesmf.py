#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 09:35:18 2019

@author: al18242
"""

import xesmf
import xarray as xr
import numpy as np

def getGridCC(lon, lat):
    """
    Create a cell centered meshgrid for from 1D arrays of lon, lat
    This meshgrid defines the bounds of each cell
    """
    dx = lon[2]-lon[1]
    dy = lat[2]-lat[1]
    lon = np.append(lon, lon[-1] + dx)
    lat = np.append(lat, lat[-1] + dy)
    lon -= dx/2.
    lat -= dy/2.
    LON, LAT = np.meshgrid(lon, lat)
    return LON, LAT

def create_xesmf_grid_uniform_cc(lat, lon):
    """
    Creates a Dataset ready to be used by the xesmf regridder from 1d arrays of lat and lon
    """
    LON, LAT = np.meshgrid(lon, lat)
    LON_b, LAT_b = getGridCC(lon, lat)
    
    grid = xr.Dataset({'lon': (['x', 'y'], LON),
                         'lat': (['x', 'y'], LAT),
                         'lon_b': (['x_b', 'y_b'], LON_b),
                         'lat_b': (['x_b', 'y_b'], LAT_b)})
    
    return grid

def regrid_uniform_cc(data, input_lat, input_lon, output_lat, output_lon):
    """
    Regrid data between two uniform, cell centered grids using a conservative method
    
    inputs
        data - numpy array
        input_lat - 1d numpy array of cell centre latitudes for data
        input_lon - 1d numpy array of cell centre longitudes for data
        output_lat - 1d numpy array of cell centre latitudes for desired output
        output_lon - 1d numpy array of cell centre longitudes for desired output
        
    returns
        data regridded as numpy array
    """
    
    input_grid = create_xesmf_grid_uniform_cc(input_lat, input_lon)
    
    output_grid = create_xesmf_grid_uniform_cc(output_lat, output_lon)
    
    regridder = xesmf.Regridder(input_grid, output_grid, 'conservative')
    regridded = regridder( data )
    return regridded

def regrid_betweenGrids(data, input_grid, output_grid, method="conservative", 
                        reuse_weights=False):
    """
    Regrid data from predefined input_grid and output_grid
    
    Inputs
        data - numpy array
        
        input_grid, output_grid -  Dataset of the form:
            
            xr.Dataset({'lat': (['x', 'y'], LAT),
                         'lon': (['x', 'y'], LON),
                         'lat_b': (['x_b', 'y_b'], LAT_b),
                         'lon_b': (['x_b', 'y_b'], LON_b)})
    
            where lat and lon give cell centre locations, and lat_b and lon_b give the cell bounds (corners)
        
        method - string describing method to use.
            Should be one of "conservative" or "conservative_normed".
            Note that to use the conservative_normed method you need to have installed the "masking" development
            branch of the xesmf package (contact Daniel for more details).
            
            With the 'masking' branch of xESMF you can include a mask in the input_grid to ignore nan values

        reuse_weights - Reuse the same weights when regridding.
            From the xesmf documentaion:
                "the output data field is linearly dependent on the input data field. 
                Any linear transform can be viewed as a matrix-vector multiplication y=Ax,
                where A is a matrix containing regridding weights, and x, y are input and 
                output data fields flatten to 1D. Computing the weight matrix A
                is expensive, but A only depends on input and output grids, 
                not on input data. 
                That means we can use the same A on different input fields x, as long 
                as the grid structure is not changed."
    returns
        regridded numpy array
    """
    #if 'mask' in input_grid:
    #    # !!! Requires the 'masking' branch of xESMF to be manually installed !!!
    #    method = 'conservative_normed'
    #else:
    #method = 'conservative'
    regridder = xesmf.Regridder(input_grid, output_grid, method,reuse_weights=reuse_weights)
    regridded = regridder( data )

    
    return regridded