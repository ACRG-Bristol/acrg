#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:09:10 2019

@author: al18242
"""
import numpy as np
import acrg_name as name
from matplotlib import pyplot as plt

def forwardModel(fp, flux):
    '''
    Calculate the predicted mf as fp*flux + bg
    
    inputs:
        fp - footprint dataset with fp_low and fp_high
        flux - flux dataset with low_res and high_res
        
    outputs:
        1. fp_low * flux_low
        2. fp_high * flux_high
        3. 1 + 2
        4. bg mf level, common to low and high resolutions
    '''
    output_model_low = np.sum( flux.low_res.values * fp.fp_low.values, axis=(0,1))
    output_model = np.sum( flux.flux.values * fp.fp.values, axis=0)
    output_model_high = np.sum( flux.high_res.values * fp.fp_high.values, axis=(0,1))
    #temp = temp + footprint[site].bc
    bc=name.boundary_conditions("EUROPE", "CH4",
                                start=fp.time.values[0],
                                end=fp.time.values[-1])
    bc=bc.reindex_like(fp, 'ffill')
    bg = (fp.particle_locations_n*bc.vmr_n).sum(["height", "lon"]) + \
        (fp.particle_locations_e*bc.vmr_e).sum(["height", "lat"]) + \
        (fp.particle_locations_s*bc.vmr_s).sum(["height", "lon"]) + \
        (fp.particle_locations_w*bc.vmr_w).sum(["height", "lat"])
    output_model = (output_model+bg) * 1e9
    output_model_low = (output_model_low+bg) * 1e9
    output_model_high = (output_model_high) * 1e9
    
    return output_model_low, output_model_high, output_model, bg*1e9

def plotMultiResMesh(data_low, lat_low, lon_low, data_high, lat_high, lon_high,
              vmin, vmax, scale_high = False, logPlot=True, cmap="viridis"):
    '''
    Plots low resolution and high resolution data on top of eachother using pcolormesh
    '''
    dlat_low = lat_low[1] - lat_low[0]
    dlon_low = lon_low[1] - lon_low[0]
    low_lat_bounds = np.append(lat_low,lat_low[-1] + dlat_low) - dlat_low/2.0
    low_lon_bounds = np.append(lon_low,lon_low[-1] + dlon_low) - dlon_low/2.0
    lats_low, lons_low = np.meshgrid(low_lat_bounds, low_lon_bounds)
    
    dlat_high = lat_high[1] - lat_high[0]
    dlon_high = lon_high[1] - lon_high[0]
    high_lat_bounds = np.append(lat_high,lat_high[-1] + dlat_high) - dlat_high/2.0
    high_lon_bounds = np.append(lon_high,lon_high[-1] + dlon_high) - dlon_high/2.0
    lats_high, lons_high = np.meshgrid(high_lat_bounds, high_lon_bounds)
    if (scale_high):
        data_high *= ((dlat_low/dlat_high)**2)
    fig, ax = plt.subplots()
    if (logPlot):
        ax.pcolormesh(lons_low, lats_low, np.log10(data_low.T), vmin=vmin, vmax=vmax, cmap=cmap)
        ax.pcolormesh(lons_high, lats_high, data_high.T*0, vmin=0, vmax=100, cmap="Greys")
        cs = ax.pcolormesh(lons_high, lats_high, np.log10(data_high.T), vmin=vmin, vmax=vmax, cmap=cmap)
    else:
        ax.pcolormesh(lons_low, lats_low, data_low.T, vmin=vmin, vmax=vmax, cmap=cmap)
        ax.pcolormesh(lons_high, lats_high, data_high.T*0, vmin=0, vmax=100, cmap="Greys")
        cs = ax.pcolormesh(lons_high, lats_high, data_high.T, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar(cs, ax=ax)
    return ax
    
def unflattenArray(stacked_array, template):
    '''
    Seperates the low res and hi res components of a flattened array
    
    inputs:
        stacked_array - flat array created in process_HiSRes
        template - dataset containing lat, lon, lat_high, lon_high
        
    output:
        1. unflattened low-res data
        2. unflattened hi-res data
    '''
    #undo the flattening+combining of high and low res grids
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
            name.process_HiSRes.getOverlapParameters(template.lat.values, template.lon.values,
                                     template.lat_high.values, template.lon_high.values)
    
    if len(stacked_array.shape) == 2:
        newaxis = stacked_array.shape[1]
    else:
        newaxis = -1
    
    out_high = stacked_array[-highsize:].reshape( (template.dims["lat_high"],template.dims["lon_high"],newaxis), order="F")
    #out_low missing values are given nans
    #insert inserts at indicies based on original array, so each insert location must be incremented to work as expected
    out_low = np.insert(stacked_array[:-highsize],
                      [ (indicies_to_remove.tolist()[i] - i) for i in range(len(indicies_to_remove))],
                      np.nan,
                      axis=0)
    out_low = out_low.reshape((template.dims["lat"],template.dims["lon"],newaxis), order="F")
    
    return out_high, out_low

def postProcess(ds, basis):
    '''
    Calculate useful output parameters from a tdmcmc output
    
    inputs:
        ds - dataset of tdmcmc output
        basis - basis file used in the tdmcmc inversion
        
    output:
        object containing various useful parameters
    '''
    output = lambda: None
    output.H_bg = ds.h_v_all[:,:4].values
    output.x_v_bg = ds.x_it[:,:4].values
    output.bg_post = np.zeros((len(ds.nIt), len(ds.nmeasure)))
    for it in range(len(ds.nIt)):
        output.bg_post[it,:] = np.dot(output.H_bg, output.x_v_bg[it,:])
        
    output.H = ds.h_v_all[:,4:ds.nIC.values].values
    output.x_v = ds.x_it[:,4:ds.nIC.values].values
    output.y_post = np.zeros((len(ds.nIt), len(ds.nmeasure)))
    for it in range(len(ds.nIt)):
        output.y_post[it,:] = np.dot(output.H, output.x_v[it,:])
    
    output.x_fine = np.zeros(len(basis.index))
    output.h_fine = np.zeros((len(basis.index),ds.h_v_all.shape[0]))
    for i in range(int(np.max(basis.basis))):
        wh_ri = np.where(basis.basis[:,0] == i+1)[0]
        output.x_fine[wh_ri] = np.mean(output.x_v[:,i],axis=0)
        output.h_fine[wh_ri,:] = np.expand_dims(output.H[:,i],axis=0)
              
    output.x_high, output.x_low = unflattenArray(output.x_fine, basis)
    output.h_high, output.h_low = unflattenArray(output.h_fine, basis)
    
    return output
