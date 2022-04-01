# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 16:58:34 2015

@author: chxmr and Mark Lunt
"""
import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube
from iris.analysis import AreaWeighted


def regrid2d(array_in, lat_in, lon_in,
             lat_out, lon_out, global_grid=False):
    '''2D mass-conservative regrid
    
    Args:
        array_in (arrray): 
            2D field to regrid
        lat_in (array): 
            latitudes corresponding to array_in
        lon_in (array): 
            longitude corresponding to array_in
        lat_out (array): 
            latitude to regrid onto
        lon_out (array): 
            longitude to regrid onto
        global_grid (boolean):
            True if grid wraps around Earth, otherwise it does not realise 
            that 0 and 360 deg longitude are the same, and your grid will have 
            a strange gap in. 
        
    Returns:
        array: 
            regridded 2D array of dimensions [lat_out, lon_out]
        iris 'Cube': regridded iris 'Cube' object
        
    Example:
        new_array, newcube = regrid2d(array_in, lat_in, lon_in, lat_out, lon_out)
    '''


    # Define input cube
    cube_lat_in = DimCoord(lat_in,
                           standard_name='latitude',
                           units='degrees')
    cube_lon_in = DimCoord(lon_in,
                           standard_name='longitude',
                           units='degrees',
                           circular=global_grid)
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
                            units='degrees',
                            circular=global_grid)
    cube_out = Cube(np.zeros((len(lat_out), len(lon_out)),
                             np.float32),
                    dim_coords_and_dims=[(cube_lat_out, 0),
                                         (cube_lon_out, 1)])
    if global_grid:
        # If global grid then need to have lat limits that stop at the poles. 
        # Guessing bounds does not do this.
        lat_bounds = np.zeros((len(lat_out),2))
        lat_bounds[1:,0] = (lat_out[1:]+lat_out[:-1])/2
        lat_bounds[:-1,1] = (lat_out[1:]+lat_out[:-1])/2
        lat_bounds[0,0] = lat_out[0]
        lat_bounds[-1,1] = lat_out[-1]
        cube_out.coord('latitude').bounds = lat_bounds[:,:]
    else:
        cube_out.coord('latitude').guess_bounds()
    cube_out.coord('longitude').guess_bounds()
    
    # Regrid
    print("Regridding...")
    cube_regridded = cube_in.regrid(cube_out,
                                    AreaWeighted(mdtol=1.))
    print(cube_regridded.summary(shorten=True))

    return cube_regridded.data,cube_regridded         


def regrid3d(array_in, lat_in, lon_in,
             lat_out, lon_out, time, global_grid = False):
    '''3D mass-conservative regrid using a cached regridder
    
    Args:
        array_in (array): 
            3D field to regrid -- Dimensions of [lat_in, lon_in, time_in]
        lat_in (array): 
            latitudes corresponding to array_in
        lon_in (array): 
            longitude corresponding to array_in
        lat_out (array): 
            latitude to regrid onto
        lon_out (array): 
            longitude to regrid onto
        time (array): 
            times to regrid onto
        global_grid (boolean):
            True if grid wraps around Earth, otherwise it does not realise 
            that 0 and 360 deg longitude are the same, and your grid will have 
            a strange gap in. 
        
    Returns
        array: 
            regridded 3D array of dimensions [lat_out, lon_out, time]
        
    Example:
        array_out = regrid3d(array_in, lat_in, lon_in, lat_out, lon_out, time)
    '''

    def get_cube_in(array_in, lat_in, lon_in, time, global_grid):
        #Define input grid
        cube_lat_in = DimCoord(lat_in,
                               standard_name='latitude',
                               units='degrees')
        cube_lon_in = DimCoord(lon_in,
                               standard_name='longitude',
                               units='degrees',
                               circular = global_grid)
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

    def get_cube_out(lat_out, lon_out, time, global_grid):
        # Define output grid
        cube_lat_out = DimCoord(lat_out,
                                standard_name='latitude',
                                units='degrees')
        cube_lon_out = DimCoord(lon_out,
                                standard_name='longitude',
                                units='degrees',
                                circular = global_grid)
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
    cube_in = get_cube_in(array_in, lat_in, lon_in, time, global_grid)    
    cube_out = get_cube_out(lat_out, lon_out,time, global_grid)     
    
    array_out = np.zeros((len(lat_out), len(lon_out), len(time)))
    
    print("Regridding...")
    cube_regridded = cube_in.regrid(cube_out,
                                    AreaWeighted(mdtol=1.))
    print(cube_regridded.summary(shorten=True))
    array_out = cube_regridded.data

    return array_out          