#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:47:07 2018

@author: lw13938

"""
import xarray as xr
import numpy as np
import os
from datetime import datetime as dt
from acrg_config.paths import paths
import cams_inversion

data_path = paths.data

def interpheight(nesw, fp_height, species, lonorlat=None, reverse=None, z_dim="level"):
    """
    Interpolates the CAMS data to the NAME heights
    
    Args:
        nesw (dataset)
            The N, E, S or W BC CAMS boundary data
        fp_height (array)
            NAME footprint heights
        species (string)
            The gas species of interest
        lonorlat (string)
            Whether you're interpolating along the 'longitude' 
            (N or S) or the 'latitude' (E or W).
        reverse (bool/None)
            Whether height values within is nesw input are in reverse order 
            (i.e. nesw["z"] level 1 values > nesw["z"] level 2 values).
            Default = None. If this is set to None this will be automatically determined.
        z_dim (str, optional) : Used if reverse is not defined to extract appropriate
            height values to compare.
            
    Returns:
        dataset : CAMS BC data interpolated to NAME heights
        
    """
    if lonorlat == 'longitude':     
        interp = np.zeros((len(fp_height),len(nesw.longitude) ))
    elif lonorlat == 'latitude':
        interp = np.zeros((len(fp_height),len(nesw.latitude) ))
    else:
        print("Please specify either lonorlat='longitude' or 'latitude'")
        return None
    
    if reverse is None:
        z_coords = nesw[z_dim].values
        z_0 = nesw['z'].sel(**{z_dim:z_coords[0]}).values[0]
        z_1 = nesw['z'].sel(**{z_dim:z_coords[1]}).values[0]

        if z_1 >= z_0:
            reverse=False
        elif z_0 > z_1:
            reverse=True        
    
    for j in range(len(nesw['z'][0,:])):
        if reverse == True:
            interp[:,j] = np.interp(fp_height, nesw['z'][:,j][::-1], nesw[species][:,j][::-1])
        elif reverse == False:
            interp[:,j] = np.interp(fp_height, nesw['z'][:,j], nesw[species][:,j])
            
    ds2 = xr.DataArray(interp, coords=[fp_height, nesw[lonorlat].values], dims=['height', lonorlat])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def interplonlat(nesw, fp_lonorlat, species, lonorlat=None, reverse=None, verbose=False):
    """
    Interpolates the CAMS data to the NAME longitudes and latitudes
    
    Args:
        nesw (dataset)
            The N, E, S or W BC CAMS boundary data
        fp_lonorlat (array)
            NAME footprint longitudes or latitudes
        species (string)
            The gas species of interest
        lonorlat (string)
            Whether you're interpolating along the 'longitude' 
            (N or S) or the 'latitude' (E or W).
        reverse (bool/None)
            Whether lon or lat values within nesw input are in reverse order
            (i.e. nesw[lonorlat][0] > nesw[lonorlat][1]).
            Default = None. If this is set to None this will be automatically determined.
            
    Returns:
        dataset : CAMS BC data interpolated to NAME longitudes or latitudes
        
    """

    if reverse is None:
        if nesw[lonorlat].values[1] >= nesw[lonorlat].values[0]:
            reverse=False
        elif nesw[lonorlat].values[0] > nesw[lonorlat].values[1]:
            reverse=True
    
    if verbose:
        print("reverse", lonorlat, reverse, nesw[lonorlat][0], nesw[lonorlat][1])
    
    interp = np.zeros(( len(nesw.height),len(fp_lonorlat) ))
    for j in range(len(nesw.height)):
        if reverse == True:
            interp[j, :] = np.interp(fp_lonorlat, nesw[lonorlat].values[::-1], nesw[species][j,:][::-1])
        elif reverse == False:
            interp[j, :] = np.interp(fp_lonorlat, nesw[lonorlat].values, nesw[species][j,:])
            
    ds2 = xr.DataArray(interp, coords=[nesw.height.values, fp_lonorlat], dims=['height', lonorlat[0:3]])
    ds2 = ds2.to_dataset(name=species)
    return ds2

def write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, st_date, species, domain, outdir, gridsize,
                           from_climatology=False):
    """
    Writes the CAMS BC data to a ncdf file.
    
    Args:
        vmr_n (array)
            Molar ratio at northern boundary
        vmr_e (array)
            Molar ratio at eastern boundary
        vmr_s (array)
            Molar ratio at western boundary
        vmr_w (array)
            Molar ratio at southern boundary
        st_date (str)
            Start date of form "YYYY-MM-dd"
        species (str)
            The gas species e.g. 'ch4', 'co2'
        domain (str)
            The domain which you want the boundary conditions for.
        outdir (str)
            Base level directory for writing output file (will be appended by /LPDM/bc/"DOMAIN").
        gridsize (int/float)
            Resolution of CAMS output in degrees.
            Possible are: 0.125, 0.25, 0.4, 0.5, 0.75, 1, 1.125, 1.5, 2, 2.5, 3
        from_climatology (bool)
            Whether BCs were produced from climatologies
    
    Returns
        netcdf file: Boundary conditions at domain boundaries
    """
    BC_edges = vmr_n.merge(vmr_e).merge(vmr_s).merge(vmr_w)
    BC_edges.expand_dims('time',2)
    BC_edges.coords['time'] = (dt.strptime(st_date, '%Y-%m-%d'))
    
    BC_edges.attrs['title']           = f"ECMWF CAMS {species} volume mixing ratios at domain edges"
    BC_edges.attrs['CAMS_resolution'] = gridsize
    BC_edges.attrs['author']          = os.getenv('USER')
    BC_edges.attrs['date_created']    = np.str(dt.today())
    
    if not os.path.isdir(os.path.join(outdir, 'LPDM', 'bc', domain.upper())):
        os.makedirs(os.path.join(outdir, 'LPDM', 'bc', domain.upper()))
    
    BC_filename = cams_inversion.bc_filename(domain           = domain,
                                             species          = species,
                                             start_date       = st_date,
                                             from_climatology = from_climatology)
#     year_month  = dt.strptime(st_date, '%Y-%m-%d').strftime('%Y%m')
#     BC_filename = os.path.join(outdir, 'LPDM', 'bc', domain.upper(), f'{species.lower()}_{domain}_{year_month}.nc')
    BC_edges.to_netcdf(path = os.path.join(outdir, 'LPDM', 'bc', domain.upper(), BC_filename), mode = 'w')