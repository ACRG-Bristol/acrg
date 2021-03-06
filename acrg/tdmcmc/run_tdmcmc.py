# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:44:30 2016

@author: ml12574

Template script for running uncorrelated trans-dimensional inversion

Uses acrg_hbtdmcmc_uncorr.f90

Please note: The default is not to perform parallel tempering and the rjmcmc 
is performed using a single chain on a single processor. If you think parallel
tempering is required this is easily changed. Speak to me (ML).

After the function definition at the top of the script the next section
is the inputs. This contains all the basic stuff you will probably want and 
have to change or tune.

Further changes to stepsizes and parameter uncertainties can be applied 
further down the script once the dimensions have been set. You will almost 
certainly need to tune your stepsizes appropriately for different parameter 
types (i.e. baseline and emissions). You will also want to adjust the 
uncertainties differently (emissions uncertainty > baseline uncertainty). 

@author: ml12574
"""
import pandas
import datetime as dt
from numba import jit
import time as run_time
import xarray as xray
import os
import re
from collections import OrderedDict
import numpy as np

import acrg.name as name
import acrg.obs as acrg_obs


@jit(nopython=True)
def closest_grid(region, lon, lat, plon, plat, pind):
    lai=0
    for la in lat:
        loi=0
        for lo in lon:
            maxdist=1e6
            for pi in pind:
                dist=(la - plat[pi])*(la - plat[pi]) + (lo - plon[pi])*(lo - plon[pi])
                if dist < maxdist:
                    region[lai, loi]=pi
                    maxdist=dist
            loi+=1
        lai+=1
    return region

def get_nsigma_y(fp_data_H,start_date, end_date, sites, 
                  nmeasure, sigma_values, bl_period=10, bl_split=False,
                  levels=None):      
    """
    Defines the indices of the measurement vector which are described by different 
    uncertainty hyperparameters (sigma_model)
    By default all sigma_models are split by site and time.
    Alternative option for subdivision by site and boundary layer depth. Can be set
    when calling the function.
    E.g.
    To split by 10 day periods:
        R_indices, ydim1,ydim2, sigma_model0 = get_nsigma_y(fp_data_H,start_date, end_date,
        sites, nmeasure, bl_split = False, bl_period=10,
        sigma_values=50.)
        
    To split by BL depth:
        R_indices, ydim1,ydim2, sigma_model0 = get_nsigma_y(fp_data_H,start_date, end_date,
        sites, nmeasure, bl_split = True, bl_period=None,
        sigma_values=50.)
        
    """  
    nsites = len(sites)      
    d0=pandas.to_datetime(start_date)
    d1=pandas.to_datetime(end_date)
    delta = d1 - d0
    ndays = delta.days
       
    y_bl=np.zeros((nmeasure))
    
    nsigma=0
    nsigma_max = np.int(np.ceil(ndays/np.float(bl_period)))
    ntime_stn=np.zeros((nsites))
    if levels is not None:
        ngroups=len(levels)-1
        
    ydim1=0
    sigma_models=[]
    
    for si in range(nsites):
        fp_data_H3 = fp_data_H[sites[si]].dropna("time", how="all")        
        nsigma_stn=0
        
        mf_time_temp=fp_data_H3.time.values
        mf_time_temp2=pandas.to_datetime(mf_time_temp)
        pblh_temp=fp_data_H3.PBLH.values
        #mf_mod_temp=fp_data_H3.mf_mod.values
        ntime_stn[si]=len(mf_time_temp)
        
        bl_start=d0
        
        if bl_split is True:
            mf_mod_temp=fp_data_H3.mf_mod.values
            for ti in range(ngroups):
                 
                wh = np.where(np.logical_and(pblh_temp>=levels[ti],
                                           pblh_temp<levels[ti+1]))
                                                                       
                if len(wh[0]) > 0:
                    y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                    
                    #sigma_models.append(sigma_values[ti])
                    if levels[ti]<499:
                        #if levels[ti]>499:
                        sigma_models.append(sigma_values)
                    else:
                        if len(wh[0]) > 1:
                            sigma_models.append(np.std(mf_mod_temp[wh[0]]))
                        else: 
                            sigma_models.append(20.)
                    nsigma_stn+=1
                    
                n_obs = len(wh[0])
                if n_obs > ydim1:
                    ydim1 = n_obs*1
                                      
    
            nsigma+=nsigma_stn
        
        else:
            for ti in range(nsigma_max):
                    bl_end=bl_start+dt.timedelta(days=bl_period)
                    
                    wh=np.where(np.logical_and(mf_time_temp2>=bl_start,
                                           mf_time_temp2<bl_end))
                    #                                
                    if len(wh[0]) > 0:
                        y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                        sigma_models.append(sigma_values)
                        nsigma_stn+=1
                        
                    bl_start=bl_start+dt.timedelta(days=bl_period)
                    n_obs = len(wh[0])
                    if n_obs > ydim1:
                        ydim1 = n_obs*1
                                          
        
            nsigma+=nsigma_stn
    
    # INDEX R
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2, np.asarray(sigma_models)

def add_local_ratio(fp_data_H,return_release=True):
    '''
    The add_local_ratio function adds an additional "local_ratio" parameter to all xarray.Dataset objects
    within fp_data_H dictionary.

    Note: footprints_data_merge() function creates a dictionary of xarray.Dataset objects with one for each 
    site. fp_sensivity() takes the output from footprints_data_merge() and adds sensivity matrix to each
    dataset in the dictionary.
    Each dataset should contain the following data variables:
        "mf"
        "fp"
        "release_lon"
        "release_lat"
        "sub_fp" (additional "sub_lon", "sub_lat" coords associated with this data variable)
 
    WARNING: Changes input object in place
 
    Args:
        fp_data_H (dict) : 
            xarray.Dataset. Expects output from name.fp_sensitivity() function.
        return_release (bool) : 
            Whether or not to also return the release_lats and release_lons
                                arrays.
                          
    Returns:
        if return_release:
             tuple(dict (xarray Datasets and identifiers), np.array, np.array) : 
                 fp_data_H with "local_ratio" data variable added, release_lons, release_lats
        else:
            dict (xarray Datasets and identifiers) : fp_data_H with "local_ratio" data variable added
    '''
    
    sites = [key for key in list(fp_data_H.keys()) if key[0] != '.']
    
    release_lons=np.zeros((len(sites)))
    release_lats=np.zeros((len(sites)))
    
    for si, site in enumerate(sites):
        release_lons[si]=fp_data_H[site].release_lon[0].values
        release_lats[si]=fp_data_H[site].release_lat[0].values
        dlon=fp_data_H[site].sub_lon[1].values-fp_data_H[site].sub_lon[0].values
        dlat=fp_data_H[site].sub_lat[1].values-fp_data_H[site].sub_lat[0].values
        local_sum=np.zeros((len(fp_data_H[site].mf)))
       
        for ti in range(len(fp_data_H[site].mf)):
            release_lon=fp_data_H[site].release_lon[ti].values
            release_lat=fp_data_H[site].release_lat[ti].values
            wh_rlon = np.where(abs(fp_data_H[site].sub_lon.values-release_lon) < dlon/2.)
            wh_rlat = np.where(abs(fp_data_H[site].sub_lat.values-release_lat) < dlat/2.)
            if np.any(wh_rlon[0]) and np.any(wh_rlat[0]):
                local_sum[ti] = np.sum(fp_data_H[site].sub_fp[
                                    wh_rlat[0][0]-2:wh_rlat[0][0]+3,wh_rlon[0][0]-2:wh_rlon[0][0]+3,ti].values)/\
                                np.sum(fp_data_H[site].fp[:,:,ti].values)
            else:
                local_sum[ti] = 0.0
            
        local_ds = xray.Dataset({'local_ratio': (['time'], local_sum)},
                                        coords = {'time' : (fp_data_H[site].coords['time'])})
    
        fp_data_H[site] = fp_data_H[site].merge(local_ds)
    
    if return_release:
        return fp_data_H,release_lons,release_lats
    else:
        return fp_data_H

def average_period(ds,av_period,dim="time"):
    '''
    The average_period function resamples the input dataset along the specified dimension to the average
    period.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing time series dimension.
        av_period (str) : 
            Averaging period to apply to timeseries dimension of dataset.
            E.g. "2H" (2 hours), "30min" or "30T" (30 minutes)
            See http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases
            for list of all frequencies available.
        dim (str, optional) : 
            Timeseries dimension within Dataset. Default = "time".
    
    Returns:
        xarray.Dataset: 
            original dataset averaged along the specified dimension
    '''
    ds_av = ds.resample(indexer={'time':av_period}).mean()
    ds_av = ds_av.dropna(dim, how="all")
    
    return ds_av

def average_period_fp(fp_data_H,av_period_site,dim="time"):
    '''
    The average_period_fp function resamples all datasets within a footprints_data_merge() dictionary over
    the averaging periods.
    Note: av_period_site should contain the same number of entries as there are sites in fp_data_H
    
    WARNING: Datasets within dictionary are changed in place.
    
    Args:
        fp_data_H (dict) : 
            Dictionary of datasets. Output from footprints_data_merge() function.
        av_period (list) : 
            Averaging periods (str objects), one for each site.
            E.g. "2H" (2 hours), "30min" or "30T" (30 minutes)        
        dim (str,optional) : 
            Timeseries dimension within Dataset. Default = "time".
        
    Returns:
        dict: 
            fp_data_H with datasets averaged along the specified dimension    
    '''
    
    sites = [key for key in list(fp_data_H.keys()) if key[0] != '.']

    #fp_data_H_av = {}
    for si, site in enumerate(sites):
        av_period = av_period_site[si]
        if av_period is not None:
            fp_data_H[site] = average_period(fp_data_H[site],av_period,dim=dim)
    
    return fp_data_H

def reorder_dims(fp_data_H,first_dims=["time"]):
    '''
    The reorder_dims function orders the dimensions to always have the first_dims before
    any other dimensions e.g. 'time'

    All other dimensions should retain their previous order as long as they are consistent
    between data variables e.g. if "height","lat" on one data variable and "lat","height" on
    another then one of these will be rearranged.
    
    Use ds.transpose() function.
    
    Args:
        fp_data_H (dict) :
            Dictionary of datasets. Output from footprints_data_merge() function.
        first_dims (list) :
            Dimensions to always include first.
            Default = ["time"]
    
    Returns:
        dict:
            fp_data_H with re-ordered dimensions
    '''

    sites = [key for key in list(fp_data_H.keys()) if key[0] != '.']
    
    for site in sites:
        fp_data_H_site = fp_data_H[site]
        
        dims = first_dims
        for dv in fp_data_H_site.data_vars:
            for d in fp_data_H_site[dv].dims:
                if d not in dims:
                    dims.append(d)
        
        all_dims = fp_data_H_site.dims
        if len(dims) != len(all_dims):
            for d in all_dims:
                if d not in dims:
                    dims.append(d)

        fp_data_H[site] = fp_data_H_site.transpose(*dims)
    
    return fp_data_H

def run_tdmcmc(sites,meas_period,av_period,species,start_date ,end_date, 
    domain,network,fp_basis_case ,bc_basis_case,rjmcmc,para_temp,
    bl_period,kmin,kmax,k_ap,nIt,burn_in,nsub,
    nbeta,beta,sigma_model_pdf,sigma_model_ap, 
    sigma_model_hparams,stepsize_sigma_y,stepsize_clon,stepsize_clat,
    stepsize_bd,stepsize_all,stepsize_pdf_p1_all,stepsize_pdf_p2_all,
    pdf_param1,pdf_param2,pdf_p1_hparam1,pdf_p1_hparam2,pdf_p2_hparam1,
    pdf_p2_hparam2,x_pdf ,pdf_param1_pdf,pdf_param2_pdf,inv_type,
    output_dir,fp_dir=None, flux_dir = None, data_dir=None, basis_dir=None, bc_basis_dir=None, bc_dir = None,
    inlet=None,instrument=None,
    tau_ap=None, tau_hparams=None, stepsize_tau=None, tau_pdf=None,
    bl_split=False, bl_levels=None, filters=None, max_level=None, site_modifier={}, prior_uncertainty=False,
    include_bias=True):
    #%%
    
    if para_temp is True:
        print ("Parallel tempering is true: have you remembered to uncomment " 
                "call OMP_SET_NUM_THREADS in the .f90 script? ")
        print ("No? Well do it and recompile with fopenmp otherwise this will take nbeta x longer")
    
    #########################################
    # READ IN DATA AND FOOTPRINTS THEN MERGE
    
    corr_type={"uncorrelated":False,
              "corr":False,
              "evencorr":True}
    data = acrg_obs.get_obs(sites, species, start_date = start_date, end_date = end_date, average = meas_period, 
                          keep_missing=corr_type[inv_type],max_level=max_level,inlet=inlet,instrument=instrument,
                          data_directory = data_dir)
    
    fp_all = name.footprints_data_merge(data, domain=domain, calc_bc=True, fp_directory = fp_dir, flux_directory = flux_dir, bc_directory = bc_dir,
                                        site_modifier=site_modifier)
    
    
    if fp_basis_case in ("INTEM"):    
        fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case='transd', basis_directory = basis_dir)
        basis_func = name.name.basis(domain = domain, basis_case = 'INTEM', basis_directory = basis_dir)
    else:                            
        fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case=fp_basis_case, basis_directory = basis_dir)
        
    fp_data_H2=name.bc_sensitivity(fp_data_H2, domain=domain,basis_case=bc_basis_case, bc_basis_directory = bc_basis_dir)
    
    ###########################################################################
    # CALCULATE DEGREE OF LOCALNESS FOR EACH FOOTPRINT
    fp_data_H2,release_lons,release_lats = add_local_ratio(fp_data_H2)
    
    for si, site in enumerate(sites): 
        fp_data_H2[site].attrs['Domain']=domain
#        fp_data_H2[site].attrs['Height']=fp_heights[site] # ** fp_heights needs to be defined **
    
    ###########################################################################
    # APPLY FILTERS TO DATASET

    if filters is not None:
        fp_data_H5 = name.filtering(fp_data_H2, filters,keep_missing=corr_type[inv_type])
    else:
        fp_data_H5 = fp_data_H2.copy()
     
    ###########################################################################
    # APPLY AVERAGING PERIOD
    fp_data_H = average_period_fp(fp_data_H5,av_period,dim="time")
    
    lat = np.asarray(fp_data_H[sites[0]].sub_lat)
    lon = np.asarray(fp_data_H[sites[0]].sub_lon)
    nlat=len(lat)
    nlon=len(lon)
    lonmin=np.min(lon)
    lonmax=np.max(lon)
    latmin=np.min(lat)
    latmax=np.max(lat)
    Ngrid = nlon*nlat  # Define underlying grid    

    #### EXPLICITLY REORDER DIMENSIONS TO ALWAYS HAVE TIME FIRST
    # Had to add this function when av_period is not applied (is None), as otherwise, the
    # dimensions are not necessarily time first and this is relied upon later when
    # constructing the inputs to pass to the FORTRAN code.
    fp_data_H = reorder_dims(fp_data_H,first_dims=["time"]) # added 09/01/2019

    ###########################################################
    # CHECK IF A BIAS VALUE NEEDS TO BE INCLUDED
    
    if include_bias:
        for site in sites:
            if "GOSAT" in site and len(sites) > 1: # Will need to update to base on platform rather than searching for "GOSAT"
                nBias = 1 
                break
        else:
            nBias = 0
    else:
        nBias = 0
    
    site_sat = [site for site in sites if "GOSAT" in site]
    
    ###########################################################
    # EVERYTHING NEEDS TO BE ARRAYS FOR MCMC
    #STACK FPs, FLUXES AND OBS
    y = []
    y_site = []
    y_time = []
    y_error=[]
    #H_bc5=[]
    #H_bc=[]
    local_ratio=[]
    pblh=[]
    wind_speed=[]
    
    sub_flux_temp = fp_all[".flux"]["all"].sel(lon=lon, lat=lat, method="ffill")
    
    for si, site in enumerate(sites):
            
        fp_data_H3 = fp_data_H[site].dropna("time", how="all")  
        attributes = [key for key in list(fp_data_H3.keys()) if key[0] != '.']  
        y.append(fp_data_H3.mf.values)     
        y_site.append([site for i in range(len(fp_data_H3.coords['time']))])
        y_time.append(fp_data_H3.coords['time'].values)
        #H_bc5.append(fp_data_H3.bc.values)
        #H_bc.append(fp_data_H3.H_bc.values)
        #sub_flux_temp = fp_data_H['.flux']['all'].sel(lon=lon, lat=lat, method="ffill") #fp_data_H3.flux.sel(lon=lon, lat=lat, method="nearest")
        local_ratio.append(fp_data_H3.local_ratio.values)
        pblh.append(fp_data_H3.PBLH.values)
        wind_speed.append(fp_data_H3.wind_speed.values)
        
        if prior_uncertainty == True:
            # To calculate prior uncertainity reduction we can set the error on the measurements to a very 
            # large value so that the measurements don't have any influence when calculating the posterior.
            # Setting error to 10000*input mol fraction values.
            y_error.append(fp_data_H3.mf.values*10000)
        elif 'dmf' in attributes:    
            y_error.append(fp_data_H3.dmf.values)
        elif 'vmf' in attributes:   
            y_error.append(fp_data_H3.vmf.values)
        else:
            print("No variability or repeatability in data file - use a default value")
            dmf_default = 0.002
            print("Default value used: {}. Only appropriate for methane".format(dmf_default))
            y_error.append(dmf_default*fp_data_H3.mf.values) # 0.002 only appropriate for methane
        
        H_bc_1=fp_data_H3.H_bc
        
        # If Bias factor included, add this as the first term for H_bc for all sites (1 for sat, 0 for other)  
        if nBias:
            if site in site_sat:
                bias = np.ones([len(H_bc_1.time),1])
            else:
                bias = np.zeros([len(H_bc_1.time),1])
            H_bias = xray.DataArray(bias,coords=[("time",H_bc_1.time.values),("region_bc",np.array([0]))])
            H_bc_1["region_bc"] = np.arange(1,len(H_bc_1["region_bc"])+1)
            H_bc_1 = xray.concat((H_bias,H_bc_1),dim="region_bc")
        
        if si ==0:
            H_fixed2=fp_data_H3.H
            H_vary2=fp_data_H3.sub_H
            q_ap2=sub_flux_temp
            #H_bc2=fp_data_H3.H_bc
            H_bc2=H_bc_1
        else:
            H_fixed2=xray.concat((H_fixed2,fp_data_H3.H), dim="time")    
            H_vary2=xray.concat((H_vary2,fp_data_H3.sub_H),dim="time" ) 
            q_ap2=xray.concat((q_ap2,sub_flux_temp), dim="time") 
            #H_bc2=xray.concat((H_bc2,fp_data_H3.H_bc), dim="time") 
            H_bc2=xray.concat((H_bc2,H_bc_1), dim="time") 

    if fp_data_H3.H.dims[0] != "time": 
            axis_insert = 0
    else:
            axis_insert = 1

    q_ap2=q_ap2["flux"].transpose("time","lat","lon")

   
    if H_fixed2.dims[0] != "time":
        H_fixed2=H_fixed2.transpose()
        H_vary2=H_vary2.transpose("time","sub_lat","sub_lon")
        #q_ap2=q_ap2.transpose("time","lat","lon")
#    if H_bc2.dims[0] !="time":
#        H_bc2=H_bc2.transpose()
        
    H_fixed=H_fixed2.values
    H_vary=H_vary2.values
    q_ap=q_ap2.values
    H_bc=H_bc2.values

    y = np.hstack(y)
    y_site = np.hstack(y_site)
    y_time = np.hstack(y_time)
    y_error=np.hstack(y_error)
    #H_bc5 = np.hstack(H_bc5)
    #H_bc = np.hstack(H_bc)
    local_ratio=np.hstack(local_ratio)
    pblh=np.hstack(pblh)
    wind_speed=np.hstack(wind_speed)
    
    q_ap0=q_ap[0,:,:].copy()
    #q_ap_v = np.ravel(q_ap0)
    nmeasure=len(y)
    h_v = np.zeros((nmeasure,Ngrid))
    #local_sum=np.zeros((nmeasure))
    for ti in range(nmeasure):                        
        # Already multiplied by q in fp_senitivity            
        h_v[ti,:] = np.ravel(H_vary[ti,:,:])   #*q_ap_v   # Create sensitivty matrix spatially vectorised

    #%%
    #################################################
    if inv_type in ('evencorr', 'corr'):
        # Define obs only where finite data exists         
        wh_temp=np.where(np.logical_and(np.isfinite(y_error),np.isfinite(y)))
        timeindex_nonzero=wh_temp[0]
        tindex_zero_temp = np.arange(nmeasure)
        timeindex_zero=np.delete(tindex_zero_temp, timeindex_nonzero)
        if len(timeindex_zero) > 0:
            y[timeindex_zero]=0.
            y_error[timeindex_zero]=1.e12
            
    ############################################################
    # Create IC  
    #nBC=len(fp_data_H[sites[0]].region_bc)
    nfixed = len(fp_data_H[sites[0]].region)
    
    numsites=len(sites)
    nmeasure_site = np.zeros((numsites))
    for ss, si in enumerate(sites):
        wh_site = np.ravel(np.where(y_site == si))
        nmeasure_site[ss]=len(wh_site)
     
    ################################################
    # CALCULATE INDICES OF Y CORRESPONDING TO DIFFERENT SIGMA_Ys   
    R_indices, ydim1, ydim2,sigma_model0 = get_nsigma_y(fp_data_H,start_date, end_date, sites, 
                  nmeasure, sigma_model_ap, bl_period=bl_period, bl_split=bl_split, 
                  levels=bl_levels)     
    
    ## Define H_bc based on time, so each month is scaled individually
    
    pd_start=pandas.to_datetime(start_date)
    pd_end=pandas.to_datetime(end_date)

    # Calculate number of months in inversion period
    if pd_end.day == 1:
        nmonths = pd_end.to_period('M') - pd_start.to_period('M')   
    else:
        nmonths = pd_end.to_period('M') - pd_start.to_period('M')+1
    
    # In pandas version for python 3.* nmonths is a MonthEnd object and so needs to be converted to an integer.
    if isinstance(nmonths,pandas.tseries.offsets.MonthEnd):
        nmonths = nmonths.n

        
    nBC_basis = len(fp_data_H[sites[0]].region_bc)
    nBC = nBC_basis*nmonths   # No. of bc_basis functions x nmonths

    nBC+=nBias # Add bias to input BC (if needed)
    #nIC=nBC+nfixed+nBias

    nIC=nBC+nfixed
    h_agg0 = np.zeros((nmeasure,k_ap+nIC)) # Will be stacked in order of Bias term (if present), H_bc, H_fixed, H_vary

    if nBias:
	    h_agg0[:,0] = H_bc[:,0] # zeroth region has been added to H_bc to be appropriate for bias term (1-sat,0-other)

    pdy_time = pandas.to_datetime(y_time)
    months = np.arange(pd_start.to_period('M').month, pd_start.to_period('M').month +nmonths)
    months2 = months.copy()
    months2[months>12]=months2[months>12]-12    # Make sure all months in range 1-12
    for mn,month in enumerate(months2):
        wh_month = np.where(pdy_time.to_period('M').month == month)[0]
        if len(wh_month > 0):
            if nBias: # 0th term is for the bias, set 1:nBC to H_bc region terms
	            h_agg0[wh_month,nBias+mn*nBC_basis:nBias+(mn+1)*nBC_basis] = H_bc[wh_month,nBias:]  # Assign H_agg separately for each month
            else:
                h_agg0[wh_month,mn*nBC_basis:(mn+1)*nBC_basis] = H_bc[wh_month,:]  # Assign H_agg separately for each month
    
    #rt17603: Added on 26/07/2018 - Bug fix: H_fixed wasn't being assigned, equivalent values were just 0.0 in h_agg0 (and h_agg)
    h_agg0[:,nBC:nIC] = H_fixed
    
    x_agg=np.zeros((k_ap+nIC))+1.
    
    if nBias:
    	x_agg[0] = 0 # Set zero scaling on bias initially and allow to vary from here

    #%%
    # Define prior model uncertainty
    ####################################################
    #sigma_model0 = np.zeros((ydim2)) 
    model_error = np.zeros(nmeasure)
    #sigma_model0[:]=sigma_model_ap
    sigma_measure=y_error.copy()
    model_error[:] = sigma_model0[0]	
    
    if inv_type == 'uncorrelated':
        sigma_model_hparams= sigma_model_hparams*sigma_model_ap
    elif inv_type == 'evencorr':
        sigma_model_hparam1=sigma_model0*sigma_model_hparams[0]
        sigma_model_hparam2=sigma_model0*sigma_model_hparams[1]
        # DEFINE TAU PARAMS AND DELTATIME
        deltatime=[float(s) for s in re.findall(r'\d+', av_period[0])]
        nsite_max=np.max(nmeasure_site)
        error_structure=np.zeros((nmeasure))+1.
    elif inv_type == 'corr':
        sigma_model_hparam1=sigma_model0*sigma_model_hparams[0]
        sigma_model_hparam2=sigma_model0*sigma_model_hparams[1]
        deltatime=np.zeros((nmeasure,nmeasure))+1.e12

        for ss, si in enumerate(sites):
            wh_site = np.ravel(np.where(y_site == si))
            nmeasure_site[ss]=len(wh_site)
            
            for whi in wh_site:
                tdelta = np.absolute(y_time[wh_site]-y_time[whi]).astype('timedelta64[m]')
                deltatime_site=tdelta/np.timedelta64(1, 'h')
                deltatime[whi,wh_site] = deltatime_site

        nsite_max=np.max(nmeasure_site)
        error_structure=np.zeros((nmeasure))+1.
        
    stepsize_sigma_y_all=np.zeros((ydim2))
    stepsize_sigma_y_all[:]=stepsize_sigma_y
    #%%
    
    # Define prior model and regions with uniform distribution
    #######################################
    kICmax=kmax+nIC              # nIC and kmax already defined at top of file
    
    # Set up different starting nuclei locations for each chain 
    plon=np.zeros((kmax,nbeta))
    plat=np.zeros((kmax,nbeta))
    regions_v=np.zeros((Ngrid,nbeta),dtype=np.uint16)
    h_agg=np.zeros((nmeasure, kICmax,nbeta))
    n0=np.zeros((nmeasure,nbeta))    
    
    
    #plon0 = np.random.uniform(lonmin, lonmax, k_ap) # Lon locs of nuclei
    #plat0 = np.random.uniform(latmin, latmax, k_ap) # Lat locs of nuclei
    
    for ib in range(nbeta):
    
        plon[:k_ap,ib] = np.random.uniform(lonmin, lonmax, k_ap) # Lon locs of nuclei
        plat[:k_ap,ib] = np.random.uniform(latmin, latmax, k_ap) # Lat locs of nuclei
        #plon[:k_ap,ib] = plon0
        #plat[:k_ap,ib] = plat0
        
        if fp_basis_case in ("INTEM"):
            basis_func.coords['lon']=fp_data_H3.lon
            basis_func.coords['lat']=fp_data_H3.lat
            regions_temp = basis_func.basis.sel(lon=lon, lat=lat, method='nearest')
            regions0 = regions_temp[:,:,0].values
        #    regions0 = basis_func.basis[:,:,0].values
            regions_v0 = np.ravel(regions0)  
            regions_v[:,ib]=regions_v0.copy() 
        
        else:
            region = np.zeros((nlat, nlon), dtype=np.uint16)
            regions0=closest_grid(region, lon, lat, plon[:k_ap,ib], plat[:k_ap,ib], \
                    np.arange(0, k_ap, dtype=np.uint16))
            regions_v0 = np.ravel(regions0)
            regions_v[:,ib]=regions_v0.copy()+1
    
        for ri in range(k_ap):
            wh_ri = np.where(regions_v0 == ri)
            for ti in range(nmeasure):
                h_agg0[ti,ri+nIC]=np.sum(h_v[ti,wh_ri])

        y_model = np.dot(h_agg0,x_agg) 
        n0_ap = y_model-y
    
        h_agg[:,:k_ap+nIC,ib] = h_agg0.copy()
        n0[:,ib]=n0_ap.copy()
    
    #################################
           
    #%%
    # MCMC Parameters
    #########################################
    nit_sub=nIt//nsub
    k=np.zeros((nbeta),dtype=np.int)+k_ap
    
    x=np.zeros((kICmax,nbeta))
    sigma_model = np.zeros((ydim2,nbeta))
    for ib in range(nbeta):  
        x[:k_ap+nIC,ib]=x_agg.copy()  
        sigma_model[:,ib]=sigma_model0.copy()
    
    # nIC1 dimension Stuff
    ############################################################
    nIC1=nIC+1
    
    #%%
    #########################################
    sigma_clon = stepsize_clon*1.
    sigma_clat = stepsize_clat*1.
    sigma_bd=np.mean(x_agg[nIC:])*stepsize_bd
    
    if inv_type in ('evencorr', 'corr'):
        tau=np.zeros((numsites,nbeta))+tau_ap
        y_pdf=1
        y_hparam1=np.min(y)*0.8
        y_hparam2=np.max(y)*1.25
        stepsize_y = 40.
        nzero=len(timeindex_zero)
        
        if nzero > 0:
            timeindex_zero_in=timeindex_zero+1
            y[timeindex_zero]=y_model[timeindex_zero]*0.95
            sigma_measure[timeindex_zero]=np.nanmean(sigma_measure)
        else:
            timeindex_zero_in=[-999]
            nzero=1       
        
    if rjmcmc == True:
        rjmcmc_in=1
    else:
        rjmcmc_in=0
        
    if para_temp == True:
        para_temp_in=1
    else:
        para_temp_in=0
        
    #BEGIN ITERATIONS
    ##################################################
    print('Starting MCMC...')
    startt = run_time.time()
    
    if inv_type == 'uncorrelated':
        if para_temp:
            import acrg_tdmcmc.tdmcmc_uncorr_pt as tdmcmc_uncorr
        else:
            from acrg.tdmcmc import tdmcmc_uncorr as tdmcmc_uncorr

        k_it, x_out, regions_out, plon_out, plat_out, sigma_model_out,sigma_y_out, \
        n0T_out,pdf_param1_out,pdf_param2_out, accept, reject, \
        accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
        accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, \
        stepsize_x_out, stepsize_p1_out, stepsize_p2_out, \
        stepsize_sigma_y_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
        accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
        accept_sigma_y_all, reject_sigma_y_all = tdmcmc_uncorr.hbtdmcmc(
        beta,k, x, h_agg,y,n0, plon, plat, regions_v, 
        pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, 
        R_indices, sigma_model_hparams, stepsize_sigma_y_all, sigma_model_pdf, 
        sigma_clon, sigma_clat, rjmcmc_in, para_temp_in,
        lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
        pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
        pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
        nIt, nsub, nit_sub, nIC, 
        nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, nIC1)
    
    elif inv_type == 'evencorr':
        if para_temp:
            import acrg.tdmcmc.tdmcmc_evencorr_pt as tdmcmc_evencorr
        else:
            from acrg.tdmcmc import tdmcmc_evencorr as tdmcmc_evencorr
        
        k_it, x_out, regions_out, plon_out, plat_out, sigma_y_out, sigma_model_out, \
        n0T_out,pdf_param1_out,pdf_param2_out, tau_out, y_out,accept, reject, \
        accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
        accept_swap, reject_swap, accept_sigma_y, reject_sigma_y, \
        accept_tau, reject_tau, accept_y, reject_y, \
        stepsize_x_out, stepsize_p1_out, stepsize_p2_out, stepsize_sigma_y_out, \
        stepsize_tau_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
        accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
        accept_sigma_y_all, reject_sigma_y_all, accept_tau_all, \
        reject_tau_all = tdmcmc_evencorr.transd_evencorr.hbtdmcmc(beta,k, x,
        h_agg,y,n0, plon, plat, regions_v, 
        pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, error_structure,
        R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y_all, sigma_model_pdf, 
        tau, tau_hparams, stepsize_tau, tau_pdf, deltatime,
        y_hparam1, y_hparam2, stepsize_y, y_pdf, timeindex_zero_in,
        sigma_clon, sigma_clat, rjmcmc_in, para_temp_in, nmeasure_site, nsite_max, 
        lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
        pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
        pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
        nIt, nsub, nit_sub, nIC,
        nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, numsites,nIC1,nzero)
     
    
    
    elif inv_type == 'corr':
        if para_temp:
            import acrg.tdmcmc.tdmcmc_corr_pt as tdmcmc_corr
        else:
            from acrg.tdmcmc import tdmcmc_corr as tdmcmc_corr
        
        k_it, x_out, regions_out, plon_out, plat_out, sigma_y_out, sigma_model_out, \
        n0T_out,pdf_param1_out,pdf_param2_out, tau_out, y_out,accept, reject, \
        accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
        accept_swap, reject_swap, accept_sigma_y, reject_sigma_y, \
        accept_tau, reject_tau, accept_y, reject_y, \
        stepsize_x_out, stepsize_p1_out, stepsize_p2_out, stepsize_sigma_y_out, \
        stepsize_tau_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
        accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
        accept_sigma_y_all, reject_sigma_y_all, accept_tau_all, \
        reject_tau_all= tdmcmc_corr.transd_corr.hbtdmcmc(beta,k, x,
        h_agg,y,n0, plon, plat, regions_v, 
        pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, error_structure,
        R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y_all, sigma_model_pdf, 
        tau, tau_hparams, stepsize_tau, tau_pdf, deltatime,
        y_hparam1, y_hparam2, stepsize_y, y_pdf, timeindex_zero_in,
        sigma_clon, sigma_clat, rjmcmc_in, para_temp_in, nmeasure_site, nsite_max, 
        lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
        pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
        pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
        nIt, nsub, nit_sub, nIC,
        nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, numsites,nIC1,nzero)
    
    endt=run_time.time()
    
    print('Finished MCMC in ', endt-startt)
        
    print('Beginning post processing')
    x_post_vit=np.zeros((nit_sub,Ngrid))    
    regions_it=np.transpose(regions_out)-1
    x_it=np.transpose(x_out)
    
    for it in range(nit_sub):    
        for zz in range(k_it[it]):
            wh_reg = np.where(regions_it[it,:] == zz)
            x_post_vit[it,wh_reg] = x_it[it,zz+nIC]
    
    print('Everything done')
    endt2 = run_time.time()
    print(endt2-endt)
    
    
    ##########################################
    h_v_all=np.zeros((nmeasure,Ngrid+nIC))  
    h_v_all[:,:nIC]=h_agg0[:,:nIC]
    h_v_all[:,nIC:]=h_v
      
    ##################################################################################
    # SAVE MCMC output in a dataset and write to netcdf
    # Set up post-mcmc dataset
    if inv_type == 'uncorrelated':
        props_temp=["birth", "death", "move", "sigma_y", "swap"]
    else:
        props_temp=["birth", "death", "move", "sigma_y", "tau", "swap"]
        
    #props = np.zeros((nIC1+len(props_temp)),dtype=object)
    props = np.zeros((nIC1+len(props_temp)),dtype="S{}".format(max([len(s) for s in props_temp])))
    accepts = np.zeros((nIC1+len(props_temp)))
    rejects = np.zeros((nIC1+len(props_temp)))
    if nBias:
    	props[0] = "bias"
    	for ii in range(1,nBC):
            props[ii]="bc"+str(ii)
    else:
    	for ii in range(0,nBC):
            props[ii]="bc"+str(ii)
    for jj in range(nfixed):
        props[jj+nBC]="fixed"+str(jj)
    props[nIC1-1] = "vary"
    props[nIC1:] = props_temp
    
    accepts[:nIC1]=accept
    rejects[:nIC1]=reject
    
    if inv_type == 'uncorrelated':
        accepts[nIC1:] = [accept_birth,accept_death,
                                 accept_move,accept_sigma_y, accept_swap]
        rejects[nIC1:] = [reject_birth,reject_death,
                                 reject_move,reject_sigma_y, reject_swap]
    else:
         accepts[nIC1:] = [accept_birth,accept_death,
                             accept_move,accept_sigma_y, accept_tau, accept_swap]
         rejects[nIC1:] = [reject_birth,reject_death,
                             reject_move,reject_sigma_y, reject_tau, reject_swap]                        
    
    #Do I need to store both x_it and x_post_vit. Can't I just store x_post_vit_all[nit,NgridIC]?
    if inv_type == 'evencorr':    
        y[timeindex_zero]=np.nan
        sigma_measure[timeindex_zero]=np.nan
        sigma_y_out[timeindex_zero,:]=np.nan
  
    post_mcmc = xray.Dataset(OrderedDict((("x_it",(["nIt", "kICmax"],x_it)),
                            ("k_it", (["nIt"],k_it)),
                            ("x_post_vit", (["nIt", "Ngrid"],x_post_vit)),
                            ("regions_it", (["nIt", "Ngrid"],regions_it)),
                            ("plon_it", (["kmax","nIt"],plon_out)),                       
                            ("plat_it", (["kmax","nIt"],plat_out)),                     
                            ("sigma_y_it", (["nmeasure", "nIt"],sigma_y_out)),
                            ("sigma_measure", (["nmeasure"],sigma_measure)),
                            ("sigma_model_it", (["ydim2", "nIt"],sigma_model_out)),
                            ("R_indices", (["ydim1", "ydim2"],R_indices)),
                            ("pdf_p1_it", (["kICmax","nIt"],pdf_param1_out)),
                            ("pdf_p2_it", (["kICmax","nIt"], pdf_param2_out)),                     
                            ("y", (["nmeasure"], y)),
                            ("y_time", (["nmeasure"], y_time)),
                            ("y_site", (["nmeasure"], y_site)),
                            ("y_prior", (["nmeasure"], y_model)),
                            ("release_lons", (["sites"], release_lons)),
                            ("release_lats", (["sites"], release_lats)),
                            ("accepts", (["proposal"],accepts)),
                            ("rejects", (["proposal"],rejects)),                          
                            ("stepsize", (["nIC1"], stepsize_x_out)),
                            ("stepsize_pdf_p1", (["nIC1"], stepsize_p1_out)),
                            ("stepsize_pdf_p2", (["nIC1"], stepsize_p2_out)),
                            ("stepsize_sigma_y", (["ydim2"], stepsize_sigma_y_out)),
                            ("h_v_all", (["nmeasure","NgridIC"],h_v_all)), 
                            ("q_ap", (["lat", "lon"],q_ap0)),
                            ("dates", (["ndates"], [start_date,end_date])),
                            ("measure_av", (["sites"], av_period)),
                            ("PBLH",  (["nmeasure"], pblh)),
                            ("wind_speed",  (["nmeasure"], wind_speed)),
                            ("local_ratio",  (["nmeasure"], local_ratio)),
                            ("nIC", nIC),
                            ("nfixed", nfixed))),
                            coords={"lon":lon, "lat": lat})    
    
    if inv_type in ('evencorr', 'corr'):
        # could use assign?
        post_mcmc.update({'tau_it': (["sites","nIt"], tau_out),
                          'stepsize_tau': stepsize_tau_out})
        
    # Add some global attributes with all further info about the run:    
    post_mcmc.attrs["bc_basis_case"]=bc_basis_case
    post_mcmc.attrs["fp_basis_case"]=fp_basis_case
    post_mcmc.attrs["iterations"]=str(nIt)
    post_mcmc.attrs["burn-in"]=str(burn_in)
    post_mcmc.attrs['Start date'] = start_date
    post_mcmc.attrs['End date'] = end_date
    post_mcmc.attrs['Filters'] = filters
    post_mcmc.attrs['Parallel tempering?'] = str(para_temp)
    post_mcmc.attrs['Inversion type'] = inv_type
        
    post_mcmc.coords["proposal"]=props
    post_mcmc.coords["sites"]=sites
    
    # Also:
    #Attributes: Sites, av_period, basis_case, accepts and rejects names
    #output_directory="/home/ml12574/work/programs/Python/my_acrg/td_uncorr/"
    
    #Output files from tdmcmc_template.py stored in the form:
    # "output_" + network + "_" + species +  "_" + date + ".nc"
    
    network_w = network.split('/')[-1]
    
    fname=os.path.join(output_dir,
                        "output_" + network_w + "_" + species + "_" + start_date + ".nc")

    for key in list(post_mcmc.keys()):
        post_mcmc[key].encoding['zlib'] = True
    post_mcmc.to_netcdf(path=fname, mode='w')

    return post_mcmc

