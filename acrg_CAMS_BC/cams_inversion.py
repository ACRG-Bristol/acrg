#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:05:12 2019

@author: rt17603
"""

import os
import numpy as np
import xarray as xr
from datetime import datetime as dt
import acrg_name.name as name
from acrg_satellite.gosat import extract_files
import acrg_CAMS_BC as CAMS_BC

data_path = os.getenv("DATA_PATH")
cams_directory = os.path.join(data_path,"ECMWF_CAMS/CAMS_inversion/")

def readCAMSInversion(start,end,cams_directory=cams_directory):
    
    species="ch4"

    search_str = "cams73_latest_{species}*.nc".format(species=species)
    
    files = extract_files(cams_directory,search_str,start=start,end=end,day=False)
    ds_list = [name.open_ds(file) for file in files]
    ds = xr.concat(ds_list,dim="time")
    
    ds = ds.rename({"CH4":species})
    
    return ds

def convertCAMSaltitude(ds):
    
    hlevels = ds.hlevel.values
    h_min = hlevels[0]
    h_max = hlevels[-1]

    low_alt = ds["altitude"].sel(**{"hlevel":slice(h_min,h_max-1)})
    high_alt = ds["altitude"].sel(**{"hlevel":slice(h_min+1,h_max)})
    
    altitude_diff = high_alt.values - low_alt.values
    z = low_alt + altitude_diff/2.
    
    z_dims = tuple([dim if dim!="hlevel" else "level" for dim in ds["altitude"].dims])
    
    ds = ds.assign(**{"z":(z_dims,z)})
    ds["z"] = ds["z"].transpose(*("time","level","latitude","longitude"))
    
    return ds

def bc_filename(domain,species,st_date,climatology=False):
    if climatology:
        clim_str = "_climatology"
    else:
        clim_str = ""
    return "{species}_{domain}_{date}{clim}.nc".format(species=species.lower(),domain=domain,date=dt.strptime(st_date, '%Y-%m-%d').strftime('%Y%m'),clim=clim_str)


def create_CAMS_BC(ds,fp_lat,fp_lon,fp_height,date,domain,species="ch4",outdir=None,cams_directory=cams_directory):

    lat_grid=np.mean(ds["latitude"][1:]-ds["latitude"][:-1])
    lon_grid=np.mean(ds["longitude"][1:]-ds["longitude"][:-1])
    gridsize="{} x {}".format(lat_grid,lon_grid)
    
    ds = convertCAMSaltitude(ds)
    ds = ds.mean('time')
    ds[species].values *= 1e-9 
    
    #Select the gridcells closest to the edges of the  domain and make sure outside of fp
    lat_n =(np.abs(ds.coords['latitude'].values - max(fp_lat))).argmin()
    if ds.coords['latitude'].values[lat_n] < np.max(fp_lat) and lat_n != 0:
        lat_n -= 1
    lat_s = (np.abs(ds.coords['latitude'].values - min(fp_lat))).argmin()
    if ds.coords['latitude'].values[lat_s] > np.min(fp_lat) and lat_s != (len(ds.coords['latitude'].values)-1):
        lat_s += 1
    lon_e =  (np.abs(ds.coords['longitude'].values - max(fp_lon))).argmin()
    if ds.coords['longitude'].values[lon_e] < max(fp_lon) and lon_e != (len(ds.coords['longitude'].values)-1):
        lon_e += 1
    lon_w =  (np.abs(ds.coords['longitude'].values - min(fp_lon))).argmin()
    if ds.coords['longitude'].values[lon_w] > min(fp_lon) and lon_w != 0:
        lon_e -= 1
    
    print("ds, z",ds["z"].dims)
    print("ds, ch4",ds["ch4"].dims)
    
    #Cut to these and then interpolate
    north = ds.sel(latitude = ds.coords['latitude'][lat_n],
                   longitude = slice(ds.coords['longitude'][lon_w],ds.coords['longitude'][lon_e])).drop(['latitude'])
    south = ds.sel(latitude = ds.coords['latitude'][lat_s],
                   longitude = slice(ds.coords['longitude'][lon_w],ds.coords['longitude'][lon_e])).drop(['latitude'])
    east = ds.sel(longitude = ds.coords['longitude'][lon_e],
                  latitude = slice(ds.coords['latitude'][lat_s],ds.coords['latitude'][lat_n])).drop(['longitude'])
    west = ds.sel(longitude = ds.coords['longitude'][lon_w],
                  latitude = slice(ds.coords['latitude'][lat_s],ds.coords['latitude'][lat_n])).drop(['longitude'])
    
    #print("z",north["z"].dims)
    #print("ch4",north["ch4"].dims)
    #print(north.dims["level"],north.dims["longitude"])
    
    #import pdb
    #pdb.set_trace()
    
    #vmr_n = CAMS_BC.interplonlat(CAMS_BC.interpheight(north, fp_height, species, reverse=False, lonorlat='longitude'), fp_lon, species, lonorlat='longitude').rename({species : 'vmr_n'})   
    #vmr_s = CAMS_BC.interplonlat(CAMS_BC.interpheight(south, fp_height, species, reverse=False, lonorlat='longitude'), fp_lon, species, lonorlat='longitude').rename({species : 'vmr_s'}) 
    #vmr_e = CAMS_BC.interplonlat(CAMS_BC.interpheight(east, fp_height, species, reverse=False, lonorlat='latitude'), fp_lat, species, reverse=False, lonorlat='latitude').rename({species : 'vmr_e'}) 
    #vmr_w = CAMS_BC.interplonlat(CAMS_BC.interpheight(west, fp_height, species, reverse=False, lonorlat='latitude'), fp_lat, species, reverse=False, lonorlat='latitude').rename({species : 'vmr_w'})      

    vmr_n = CAMS_BC.interplonlat(CAMS_BC.interpheight(north, fp_height, species, lonorlat='longitude'), fp_lon, species, lonorlat='longitude').rename({species : 'vmr_n'})   
    vmr_s = CAMS_BC.interplonlat(CAMS_BC.interpheight(south, fp_height, species, lonorlat='longitude'), fp_lon, species, lonorlat='longitude').rename({species : 'vmr_s'}) 
    vmr_e = CAMS_BC.interplonlat(CAMS_BC.interpheight(east, fp_height, species, lonorlat='latitude'), fp_lat, species, lonorlat='latitude').rename({species : 'vmr_e'}) 
    vmr_w = CAMS_BC.interplonlat(CAMS_BC.interpheight(west, fp_height, species, lonorlat='latitude'), fp_lat, species, lonorlat='latitude').rename({species : 'vmr_w'})      

    CAMS_BC.write_CAMS_BC_tonetcdf(vmr_n, vmr_e, vmr_s, vmr_w, date, species, domain, outdir, gridsize)  