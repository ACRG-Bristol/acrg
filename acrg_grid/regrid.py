# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 16:58:34 2015

@author: chxmr and Mark Lunt
"""

import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube
from iris.analysis import AreaWeighted
import datetime as dt


def regrid2d(array_in, lat_in, lon_in,
             lat_out, lon_out):
    '''
    2D mass-conservative regrid
    
    Only handles a 2D array at the moment. We should add time.
    
    Inputs:
        array_in: 2D field to regrid
        lat_in: latitudes corresponding to array_in
        lon_in: longitude corresponding to array_in
        lat_out: latitude to regrid onto
        lon_out: longitude to regrid onto
        
    Returns a 2D array of dimensions [lat_out, lon_out]
    '''


    # Define input cube
    cube_lat_in = DimCoord(lat_in,
                           standard_name='latitude',
                           units='degrees')
    cube_lon_in = DimCoord(lon_in,
                           standard_name='longitude',
                           units='degrees')
    cube_in = Cube(array_in,
                   dim_coords_and_dims=[(cube_lat_in, 0),
                                        (cube_lon_in, 1)])                                   
    cube_in.coord('latitude').guess_bounds()
    cube_in.coord('longitude').guess_bounds()

    # Define output grid
    cube_lat_out = DimCoord(lat_out,
                            standard_name='latitude',
                            units='degrees')
    cube_lon_out = DimCoord(lon_out,
                            standard_name='longitude',
                            units='degrees')
    cube_out = Cube(np.zeros((len(lat_out), len(lon_out)),
                             np.float32),
                    dim_coords_and_dims=[(cube_lat_out, 0),
                                         (cube_lon_out, 1)])
    cube_out.coord('latitude').guess_bounds()
    cube_out.coord('longitude').guess_bounds()
    
    # Regrid
    print("Regridding...")
    cube_regridded = cube_in.regrid(cube_out,
                                    AreaWeighted(mdtol=1.))
    print(cube_regridded.summary(shorten=True))

    return cube_regridded.data,cube_regridded         


def regrid3d(array_in, lat_in, lon_in,
             lat_out, lon_out, time):
             
    '''
    3D mass-conservative regrid using a cached regridder
    
    Inputs:
        array_in: 3D field to regrid - Lat, Lon, Time
        lat_in: latitudes corresponding to array_in
        lon_in: longitude corresponding to array_in
        lat_out: latitude to regrid onto
        lon_out: longitude to regrid onto
        
    Returns a 3D array of dimensions [lat_out, lon_out, time_out]
    '''

    def get_cube_in(array_in, lat_in, lon_in, time):
        #Define input grid
        cube_lat_in = DimCoord(lat_in,
                               standard_name='latitude',
                               units='degrees')
        cube_lon_in = DimCoord(lon_in,
                               standard_name='longitude',
                               units='degrees')
        cube_time = DimCoord(time,
                               standard_name='time',
                               units='seconds')
        cube_in = Cube(array_in[:,:,:],
                       dim_coords_and_dims=[(cube_lat_in, 0),
                                        (cube_lon_in, 1),
                                        (cube_time, 2)])                                   
        cube_in.coord('latitude').guess_bounds()
        cube_in.coord('longitude').guess_bounds()
        cube_in.coord('time').guess_bounds()
        
        return cube_in

    def get_cube_out(lat_out, lon_out, time):
        # Define output grid
        cube_lat_out = DimCoord(lat_out,
                                standard_name='latitude',
                                units='degrees')
        cube_lon_out = DimCoord(lon_out,
                                standard_name='longitude',
                                units='degrees')
        cube_time = DimCoord(time,
                               standard_name='time',
                               units='seconds')
        cube_out = Cube(np.zeros((len(lat_out), len(lon_out), len(time)),
                                 np.float32),
                                 dim_coords_and_dims=[(cube_lat_out, 0),
                                                      (cube_lon_out, 1),
                                                      (cube_time, 2)])
        cube_out.coord('latitude').guess_bounds()
        cube_out.coord('longitude').guess_bounds()
        cube_out.coord('time').guess_bounds()
        
        return cube_out
    
    # Regrid
    print("Getting cube in and cube out")
    cube_in = get_cube_in(array_in, lat_in, lon_in, time)    
    cube_out = get_cube_out(lat_out, lon_out,time)     
    
    array_out = np.zeros((len(lat_out), len(lon_out), len(time)))
    
    print("Regridding...")
    cube_regridded = cube_in.regrid(cube_out,
                                    AreaWeighted(mdtol=1.))
    print(cube_regridded.summary(shorten=True))
    array_out = cube_regridded.data

    return array_out          