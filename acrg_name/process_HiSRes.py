#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 12:09:42 2019

@author: al18242

Functions used to generate and process HiSRes (High spatial resolution) footprints for NAME
"""

import acrg_name as name
import xarray as xr
import numpy as np
import os
import pandas as pd
from acrg_name import emissions_helperfuncs as emfuncs

def process_all(domain, site, height,
                force_met_empty = False,
                processed_folder = "Processed_Fields_files",
                processed_folder_HR = "Processed_Fields_files_HR",
                base_dir = "/data/al18242/name_out/"):
    """
    Process high resolution footprints, by seperately processing the low and high resolution outputs
    These are then combined together into a single dataset ready for further use
    """
    
    #process low resolution:
    name.process.process_all(domain, site, 
                         height,
                         force_met_empty = force_met_empty,
                         processed_folder = processed_folder,
                         base_dir=base_dir)
    #process high resolution:
    name.process.process_all(domain, site, 
                             height,
                             force_met_empty = force_met_empty,
                             fields_folder = "Fields_files_HR",
                             processed_folder = processed_folder_HR,
                             base_dir=base_dir)
    
    #open processed datasets to be combined
    output_base = "{}{}_{}_{}/".format(base_dir,domain,site,height)
    filenames = os.listdir("{}{}".format(output_base, processed_folder_HR))
    for i, filename in enumerate(filenames):
        combine_date(output_base, processed_folder, processed_folder_HR, filename)
    
def getOverlapParameters(lat_low, lon_low, lat_high, lon_high):
    """
    calculates parameters from the low and high resolution grids for flattening and computing overlap
    """
    lats, lons = np.meshgrid(lat_low, lon_low)
    lowsize = len(lon_low) * len(lat_low)
    lons_low=lons.reshape(lowsize)
    lats_low=lats.reshape(lowsize)
    
    #process high
    lats, lons = np.meshgrid(lat_high, lon_high)
    highsize = len(lon_high) * len(lat_high)
    lons_high=lons.reshape(highsize)
    lats_high=lats.reshape(highsize)
    
    #get indicies where grids overlap
    lon_first_i = np.where(lon_low > lons_high[0])[0][0]
    lon_last_i = np.where(lon_low < lons_high[-1])[0][-1]
    
    lat_first_i = np.where(lat_low > lats_high[0])[0][0]
    lat_last_i = np.where(lat_low < lats_high[-1])[0][-1]
    
    lon_indicies = np.arange(lon_first_i, lon_last_i+1,1)
    lat_indicies = np.arange(lat_first_i, lat_last_i+1,1)
    
    indicies_to_remove_lon = np.where(np.in1d(lons_low, lon_low[lon_indicies]))[0]
    indicies_to_remove_lat = np.where(np.in1d(lats_low, lat_low[lat_indicies]))[0]
    indicies_to_remove = np.intersect1d(indicies_to_remove_lon,indicies_to_remove_lat)
    
    lons_out = np.delete(lons_low, indicies_to_remove)
    lats_out = np.delete(lats_low, indicies_to_remove)
    lons_out = np.append(lons_out,lons_high)
    lats_out = np.append(lats_out,lats_high)
    
    return (lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out)
    
def combine_date(output_base, processed_folder, processed_folder_HR, filename):
    """
    Combine the low and high resolution outputs into a single file
    """
    fp_low = xr.open_dataset("{}{}/{}".format(output_base, processed_folder, filename))
    fp_high = xr.open_dataset("{}{}/{}".format(output_base, processed_folder_HR, filename))
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(fp_low.lat.values, fp_low.lon.values, fp_high.lat.values, fp_high.lon.values)
    
    #flatten footprints ready to combine
    fp_low_flat = fp_low.copy(deep=True)
    fp_low_flat["fp"] = (["index", "time"] ,fp_low.fp.values.reshape(lowsize, -1, order="F"))
    fp_high_flat = fp_high.copy(deep=True)
    fp_high_flat["fp"] = (["index", "time"] ,fp_high.fp.values.reshape(highsize, -1, order="F"))
    
    lons_out = np.delete(lons_low, indicies_to_remove)
    lats_out = np.delete(lats_low, indicies_to_remove)
    fp_out = np.delete(fp_low_flat.fp, indicies_to_remove, axis=0)
    
    lons_out=np.append(lons_out,lons_high)
    lats_out=np.append(lats_out,lats_high)
    fp_out=np.append(fp_out,fp_high_flat.fp,axis=0)
    
    #create output file from existing
    output_file = fp_low.copy(deep=True)
    output_file["fp_low"] = fp_low.fp
    output_file["fp_high"] = (["lat_high", "lon_high", "time"], fp_high.fp.values)
    output_file["lat_high"] = (["lat_high"], fp_high.lat.values, {"units":"Degrees_north"})
    output_file["lon_high"] = (["lon_high"], fp_high.lon.values, {"units":"Degrees_east"})
    output_file["fp"] = (["index", "time"], fp_out)
    
    output_file["index_lons"] = (["index"], lons_out)
    output_file["index_lats"] = (["index"], lats_out)
    
    #save output
    out_directory = "{}Processed_Fields_files_combined".format(output_base)
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)
    output_file.to_netcdf("{}/{}".format(out_directory, filename))
    
def getFlux(ds, output_dir, name):
    """
    get Flux using edgar and NAEI for footprint file given
    Currently only good for static annual emissions
    """
    year = pd.to_datetime(ds.time.values[0]).year
    
    edgar_low = emfuncs.getedgarannualtotals(year, ds.lon.values, ds.lat.values, species='CH4')
    naei_low = emfuncs.getNAEI(year, ds.lon.values, ds.lat.values, species="CH4", naei_sector="total")
    combined = naei_low.data
    combined[np.where(naei_low.mask == True)] = edgar_low[np.where(naei_low.mask == True)]
    combined[np.where(naei_low.data == 0)] = edgar_low[np.where(naei_low.data == 0)]
    naei_high = emfuncs.getNAEI(year, ds.lon_high.values, ds.lat_high.values, species="CH4", naei_sector="total")
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(ds.lat.values, ds.lon.values, ds.lat_high.values, ds.lon_high.values)
        
    combined_low = combined
    combined = combined.reshape(lowsize, order="F")
    combined = np.delete(combined, indicies_to_remove)
    combined = np.append(combined, naei_high.reshape(highsize, order="F"))
    
    data_vars = {
            "low_res":(["lat", "lon", "time"], np.expand_dims(combined_low, axis=2)),
            "high_res":(["lat_high", "lon_high", "time"], np.expand_dims(naei_high,axis=2)),
            "flux":(["index", "time"], np.expand_dims(combined,axis=1)),
            "index_lat":(["index"], lats_out),
            "index_lon":(["index"], lons_out)
            }
    coords = {
            "lat":np.unique(lats_low),
            "lon":np.unique(lons_low),
            "lat_high":np.unique(lats_high),
            "lon_high":np.unique(lons_high),
            "time":[np.datetime64(str( pd.to_datetime(ds.time.values[0]).year), 'Y')]
            }
    output = xr.Dataset(data_vars=data_vars, coords=coords)
    output.to_netcdf("{}/{}_{}.nc".format(output_dir, name,year))
    
#    import acrg_name.flux as flux
#    from datetime import datetime
#    time = np.array([np.datetime64("2015-01-01")])#np.array([np.datetime64(pd.to_datetime("2015-01-01"))])
#    time=np.array(['2015-01-01']).astype('datetime64[ns]')
#    flux.write(lat,lon,time,np.expand_dims(combined,2),"CH4-NAEI-fixed","EUROPE",source=None,title="NAEI in EDGAR",
#                   prior_info_dict={"TODO":["todo","todo","todo"]},flux_comments="NAEI + EDGAR",climatology=False,
#                   regridder_used='acrg_grid.regrid.regrid_2D',output_directory="/data/al18242/flux/")

def makeBasisFromExisting():
    #flux = name.name.flux("EUROPE", "CH4-BTT-5", flux_directory="/data/al18242/flux_HR/")
    with xr.open_dataset("/data/al18242/flux_hr/ch4-BTT-5.nc") as ds:
        flux = ds.load()
    with xr.open_dataset("/data/shared/NAME/basis_functions/EUROPE/sub-country_mask_uk-split_EUROPE_2014.nc") as ds:
        existing = ds.load()
    basis_out = existing.copy(deep=True)
    
    #do the high res section as a simple grid for now
    start_region = np.amax(existing.basis.values) + 1
    high_res_basis = np.zeros_like(flux.high_res)
    high_res_basis[:,:,0] = start_region
    xVals = np.arange(22,55)
    yVals = np.arange(22,55)
    X, Y = np.meshgrid(yVals, xVals)
    Xc = np.floor( (X-np.amin(X))/1.0)
    Yc = np.floor( (Y-np.amin(Y))/1.0)
    rows = np.amax(Yc-np.amin(Yc))-1
    high_res_basis[X.astype(int),Y.astype(int),0] = 1+start_region + Yc + Xc * rows
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(flux.lat_low.values, flux.lon_low.values, flux.lat_high.values, flux.lon_high.values)
        
    combined = existing.basis.values.reshape(lowsize, -1,order="F")
    combined = np.delete(combined, indicies_to_remove, axis=0)
    combined = np.append(combined, high_res_basis.reshape(highsize, -1,order="F"),axis=0)
    basis_out["basis_old"] = basis_out["basis"].copy(deep=True)
    basis_out["basis"] = (["index", "time"], combined.astype(np.int32))
    basis_out["lat_high"] = flux.lat_high
    basis_out["lon_high"] = flux.lon_high
    
    basis_out.to_netcdf("/data/al18242/basis_hr/EUROPE/sub_btt-1-london_EUROPE_2000.nc")
    
    