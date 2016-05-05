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
import tdmcmc_uncorr
import tdmcmc_evencorr
import acrg_name as name
import numpy as np
import acrg_agage as agage
import pandas
import datetime as dt
from numba import jit
import time as run_time
import xray
import os
import re


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

def get_nsigma_y(fp_data_H,start_date, end_date, bl_period, sites, 
                  nmeasure,threshold):
    
    nsites = len(sites)    
    d0=pandas.to_datetime(start_date)
    d1=pandas.to_datetime(end_date)
    delta = d1 - d0
    ndays = delta.days
       
    y_bl=np.zeros((nmeasure))
    
    nsigma=0
    nsigma_max = np.int(np.ceil(ndays/np.float(bl_period)))
    ntime_stn=np.zeros((nsites))
    
    ydim1=0
    
    for si in range(nsites):
        fp_data_H3 = fp_data_H[sites[si]].dropna("time", how="all")        
        nsigma_stn=0
        mf_time_temp=fp_data_H3.time.values       
        mf_time_temp2=pandas.to_datetime(mf_time_temp)
        
        ntime_stn[si]=len(mf_time_temp)
        for ji in range(2):
            bl_start=d0
            if ji == 0:
        
                for ti in range(nsigma_max):
                    bl_end=bl_start+dt.timedelta(days=bl_period)
                    
                    wh=np.where(np.logical_and(np.logical_and(mf_time_temp2>=bl_start,
                                           mf_time_temp2<bl_end),
                                           fp_data_H3.local_ratio < threshold))
                    # NEED ANOTHER INDEX IF ABOVE LOCAL RATIO THRESHOLD                               
                    if len(wh[0]) > 0:
                        y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                        nsigma_stn+=1
                        
                    bl_start=bl_start+dt.timedelta(days=bl_period)
                    n_obs = len(wh[0])
                    if n_obs > ydim1:
                        ydim1 = n_obs*1
                        
            elif ji ==1:
                 for ti in range(nsigma_max):
                    bl_end=bl_start+dt.timedelta(days=bl_period)
                    
                    wh=np.where(np.logical_and(np.logical_and(mf_time_temp2>=bl_start,
                                           mf_time_temp2<bl_end),
                                           fp_data_H3.local_ratio >= threshold))
                    # NEED ANOTHER INDEX IF ABOVE LOCAL RATIO THRESHOLD                               
                    if len(wh[0]) > 0:
                        y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                        nsigma_stn+=1
                        
                    bl_start=bl_start+dt.timedelta(days=bl_period)
                    n_obs = len(wh[0])
                    if n_obs > ydim1:
                        ydim1 = n_obs*1
    
        nsigma+=nsigma_stn
    
    # INDEX R
    print ydim1, nsigma
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2



def run_tdmcmc(sites,meas_period,av_period,species,start_date ,end_date, 
    domain,network,fp_basis_case ,bc_basis_case,rjmcmc,bl_period,kmin,kmax,
    k_ap,nIt,burn_in,nsub,threshold,nbeta,beta,sigma_model_pdf,sigma_model_ap, 
    sigma_model_hparams,stepsize_sigma_y,stepsize_clon,stepsize_clat,
    stepsize_bd,stepsize_all,stepsize_pdf_p1_all,stepsize_pdf_p2_all,
    pdf_param1,pdf_param2,pdf_p1_hparam1,pdf_p1_hparam2,pdf_p2_hparam1,
    pdf_p2_hparam2,x_pdf ,pdf_param1_pdf,pdf_param2_pdf,inv_type,
    output_dir,tau=None, tau_hparams=None, stepsize_tau=None, tau_pdf=None):
    #%%
    #########################################
    # READ IN DATA AND FOOTPRINTS THEN MERGE
    corr_type={"uncorrelated":False,
              "correlated":False,
              "evencorr":True}
    data = agage.get_obs(sites, species, start = start_date, end = end_date, average = meas_period, 
                          keep_missing=corr_type[inv_type])
    
    
    fp_all = name.footprints_data_merge(data, domain=domain, species=species, calc_bc=True)
                                        
    fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case=fp_basis_case)
    fp_data_H2=name.bc_sensitivity(fp_data_H2, domain=domain,basis_case=bc_basis_case)
    
    ###########################################################################
    # CALCULATE DEGREE OF LOCALNESS FOR EACH FOOTPRINT
    release_lons=np.zeros((len(sites)))
    release_lats=np.zeros((len(sites)))
    for si, site in enumerate(sites):
        release_lons[si]=fp_data_H2[site].release_lon[0].values
        release_lats[si]=fp_data_H2[site].release_lat[0].values
        dlon=fp_data_H2[site].sub_lon[1].values-fp_data_H2[site].sub_lon[0].values
        dlat=fp_data_H2[site].sub_lat[1].values-fp_data_H2[site].sub_lat[0].values
        wh_rlon = np.where(abs(fp_data_H2[site].sub_lon.values-release_lons[si]) < dlon/2.)
        wh_rlat = np.where(abs(fp_data_H2[site].sub_lat.values-release_lats[si]) < dlat/2.)
        local_sum=np.zeros((len(fp_data_H2[site].mf)))
        sub_lon=fp_data_H2[site].sub_lon.values
        sub_lat=fp_data_H2[site].sub_lat.values
        sub_flux_temp = fp_data_H2[site].flux.sel(lon=slice(np.min(sub_lon),np.max(sub_lon)), 
                                        lat=slice(np.min(sub_lat),np.max(sub_lat)))
        for ti in range(len(fp_data_H2[site].mf)):      
            local_sum[ti] = np.sum(fp_data_H2[site].sub_fp[
            wh_rlat[0]-1:wh_rlat[0]+2,wh_rlon[0]-1:wh_rlon[0]+2,ti].values*
            sub_flux_temp[
            wh_rlat[0]-1:wh_rlat[0]+2,wh_rlon[0]-1:wh_rlon[0]+2,ti].values)/np.sum(
            fp_data_H2[site].sub_fp[:,:,ti].values*sub_flux_temp[:,:,ti])
            
        local_ds = xray.Dataset({'local_ratio': (['time'], local_sum)},
                                        coords = {'time' : (fp_data_H2[site].coords['time'])})
    
        fp_data_H2[site] = fp_data_H2[site].merge(local_ds)
    
    fp_data_H = {}
    
    for si, site in enumerate(sites):
        site_ds = fp_data_H2[site].resample(av_period[si], dim = "time")
        site_ds2= site_ds.dropna("time", how="all")
        fp_data_H[site] = site_ds2
        
    
    lat = np.asarray(fp_data_H[sites[0]].sub_lat)
    lon = np.asarray(fp_data_H[sites[0]].sub_lon)
    nlat=len(lat)
    nlon=len(lon)
    lonmin=np.min(lon)
    lonmax=np.max(lon)
    latmin=np.min(lat)
    latmax=np.max(lat)
    Ngrid = nlon*nlat  # Define underlying grid    
    
    ###########################################################
    # EVERYTHING NEEDS TO BE ARRAYS FOR MCMC
    #STACK FPs, FLUXES AND OBS
    y = []
    y_site = []
    y_time = []
    y_error=[]
    H_bc5=[]
    local_ratio=[]
    
    for si, site in enumerate(sites):
              
        fp_data_H3 = fp_data_H[site].dropna("time", how="all")  
        attributes = [key for key in fp_data_H3.keys() if key[0] != '.']  
        y.append(fp_data_H3.mf.values)     
        y_site.append([site for i in range(len(fp_data_H3.coords['time']))])
        y_time.append(fp_data_H3.coords['time'].values)
        H_bc5.append(fp_data_H3.bc.values)
        sub_flux_temp = fp_data_H3.flux.sel(lon=slice(lonmin,lonmax), 
                                        lat=slice(latmin,latmax))
        local_ratio.append(fp_data_H3.local_ratio.values)
        if 'dmf' in attributes:    
            y_error.append(fp_data_H3.dmf.values)
        elif 'vmf' in attributes:   
            y_error.append(fp_data_H3.vmf.values)
        else:
            y_error.append(0.002*fp_data_H3.mf.values) # 0.002 only appropriate for methane
        if si ==0:
            H_fixed2=fp_data_H3.H
            H_vary2=fp_data_H3.sub_H
            q_ap2=sub_flux_temp
            H_bc2=fp_data_H3.H_bc
                  
        else:
            H_fixed2=xray.concat((H_fixed2,fp_data_H3.H), dim="time")    
            H_vary2=xray.concat((H_vary2,fp_data_H3.sub_H),dim="time" ) 
            q_ap2=xray.concat((q_ap2,sub_flux_temp), dim="time") 
            H_bc2=xray.concat((H_bc2,fp_data_H3.H_bc), dim="time") 
    
    if H_fixed2.dims[0] != "time":
        H_fixed2=H_fixed2.transpose()
        H_vary2=H_vary2.transpose("time","sub_lat","sub_lon")
        q_ap2=q_ap2.transpose("time","lat","lon")
    if H_bc2.dims[0] !="time":
        H_bc2=H_bc2.transpose()
        
    H_fixed=H_fixed2.values
    H_vary=H_vary2.values
    q_ap=q_ap2.values
    H_bc=H_bc2.values
     
    y = np.hstack(y)
    y_site = np.hstack(y_site)
    y_time = np.hstack(y_time)
    y_error=np.hstack(y_error)
    H_bc5 = np.hstack(H_bc5)
    local_ratio=np.hstack(local_ratio)
    
    q_ap0=q_ap[0,:,:].copy()
    #q_ap_v = np.ravel(q_ap0)
    nmeasure=len(y)
    print nmeasure
    h_v = np.zeros((nmeasure,Ngrid))
    local_sum=np.zeros((nmeasure))
    for ti in range(nmeasure):                        
        # Already multiplied by q in fp_senitivity            
        h_v[ti,:] = np.ravel(H_vary[ti,:,:])   #*q_ap_v   # Create sensitivty matrix spatially vectorised
    
    #%%
    #################################################
    if inv_type is 'evencorr':
        # Define obs only where finite data exists         
        wh_temp=np.where(np.logical_and(np.isfinite(y_error),np.isfinite(y)))
        timeindex_nonzero=wh_temp[0]
        tindex_zero_temp = np.arange(nmeasure)
        timeindex_zero=np.delete(tindex_zero_temp, timeindex_nonzero)
        if len(timeindex_zero > 0):
            y[timeindex_zero]=0.
            y_error[timeindex_zero]=1.e12
            
    ############################################################
    # Create IC  
    nBC=len(fp_data_H[sites[0]].region_bc)
    nfixed = len(fp_data_H[sites[0]].region)
    
    numsites=len(sites)
    nmeasure_site = np.zeros((numsites))
    for ss, si in enumerate(sites):
        wh_site = np.ravel(np.where(y_site == si))
        nmeasure_site[ss]=len(wh_site)
     
    ################################################
    # CALCULATE INDICES OF Y CORRESPONDING TO DIFFERENT SIGMA_Ys   
    R_indices, ydim1, ydim2 = get_nsigma_y(fp_data_H,start_date, end_date, bl_period, sites, 
                      nmeasure,threshold)
    nIC=nBC+nfixed
    
    
    h_agg0 = np.zeros((nmeasure,k_ap+nIC))
    x_agg=np.zeros((k_ap+nIC))+1.  
    h_agg0[:,:nBC]=H_bc.copy()
    h_agg0[:,nBC:nIC]=H_fixed.copy()
    
    #%%
    # Define prior model uncertainty
    ####################################################
    sigma_model0 = np.zeros((ydim2)) 
    model_error = np.zeros(nmeasure)
    sigma_model0[:]=sigma_model_ap
    sigma_measure=y_error.copy()
    model_error[:] = sigma_model0[0]	
    
    if inv_type is 'evencorr':
        sigma_model_hparam1=sigma_model0*sigma_model_hparams[0]
        sigma_model_hparam2=sigma_model0*sigma_model_hparams[1]
        # DEFINE TAU PARAMS AND DELTATIME
        deltatime=[float(s) for s in re.findall(r'\d+', av_period[0])]
        nsite_max=np.max(nmeasure_site)
        
        
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
    
    
    plon0 = np.random.uniform(lonmin, lonmax, k_ap) # Lon locs of nuclei
    plat0 = np.random.uniform(latmin, latmax, k_ap) # Lat locs of nuclei
    
    for ib in range(nbeta):
    
        #plon[:k_ap,ib] = np.random.uniform(lonmin, lonmax, k_ap) # Lon locs of nuclei
        #plat[:k_ap,ib] = np.random.uniform(latmin, latmax, k_ap) # Lat locs of nuclei
        plon[:k_ap,ib] = plon0
        plat[:k_ap,ib] = plat0
        
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
    nit_sub=nIt/nsub
    
    
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
    
    
    #BEGIN ITERATIONS
    
    ##################################################
    print 'Starting MCMC...'
    startt = run_time.time()
    
    if inv_type is 'uncorrelated':
        k_it, x_out, regions_out, plon_out, plat_out, sigma_model_out,sigma_y_out, \
        n0T_out,pdf_param1_out,pdf_param2_out, accept, reject, \
        accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
        accept_swap, accept_sigma_y, reject_sigma_y, \
        tot_acc_x, tot_acc_p1, tot_acc_p2, tot_acc_sigma_y = tdmcmc_uncorr.hbtdmcmc(
        beta,k, x, h_agg,y,n0, plon, plat, regions_v, 
        pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, 
        R_indices, sigma_model_hparams, stepsize_sigma_y_all, sigma_model_pdf, 
        sigma_clon, sigma_clat, rjmcmc, 
        lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
        pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
        pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
        nIt, nsub, nit_sub, nIC, 
        nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, nIC1)
    
    elif inv_type is 'evencorr':
        k_it, x_out, regions_out, plon_out, plat_out, sigma_y_out, sigma_model_out, \
        n0T_out,pdf_param1_out,pdf_param2_out, tau_out, accept, reject, \
        accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
        accept_swap, reject_swap, accept_sigma_y, reject_sigma_y, \
        accept_tau, reject_tau = tdmcmc_evencorr.transd_evencorr.hbtdmcmc(beta,k, x,
        h_agg,y,n0, plon, plat, regions_v, 
        pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, 
        R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y_all, sigma_model_pdf, 
        tau, tau_hparams, stepsize_tau, tau_pdf, deltatime,
        sigma_clon, sigma_clat, rjmcmc, nmeasure_site, nsite_max, 
        lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
        pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
        pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
        nIt, nsub, nit_sub, nIC,
        nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, numsites,nIC1)
    
    endt=run_time.time()
    
    print 'Finished MCMC in ', endt-startt
    
    
    print 'Beginning post processing'
    x_post_vit=np.zeros((nit_sub,Ngrid))
    x_post_v_mean=np.zeros((Ngrid))
    x_post_v_95=np.zeros((Ngrid))
    x_post_v_05=np.zeros((Ngrid))
    
    # Having run MCMC then need to map regions back onto grid
    sigma_y_mean=np.zeros((nmeasure))
    
    for yi in range(nmeasure):
        sigma_y_mean[yi] = np.mean(sigma_y_out[yi,:])
    
    regions_it=np.transpose(regions_out)-1
    x_it=np.transpose(x_out)
    
    for it in range(nit_sub):    
        for zz in range(k_it[it]):
            wh_reg = np.where(regions_it[it,:] == zz)
            x_post_vit[it,wh_reg] = x_it[it,zz+nIC]
    
    endt3 = run_time.time()
    print endt3-endt
    for ii in range(Ngrid):
        x_post_v_mean[ii] = np.mean(x_post_vit[:,ii]) 
        x_post_v_95[ii] = np.percentile(x_post_vit[:,ii],95)
        x_post_v_05[ii] = np.percentile(x_post_vit[:,ii],5)
    
    x_post_mean = np.reshape(x_post_v_mean, (nlat,nlon))
    xrange_90v=x_post_v_95-x_post_v_05
    xrange_90=np.reshape(xrange_90v, (nlat,nlon))
    
    print 'Everything done'
    endt2 = run_time.time()
    print endt2-endt
    
    ##########################################
    # y_post
    h_v_all=np.zeros((nmeasure,Ngrid+nIC))
    x_post_all = np.zeros(Ngrid+nIC)
    
    for xi in range(nIC):
        x_post_all[xi]=np.mean(x_it[:,xi])
    
    x_post_all[nIC:]=x_post_v_mean    
    h_v_all[:,:nIC]=h_agg0[:,:nIC]
    h_v_all[:,nIC:]=h_v
    
    x_post_all_it=np.zeros((nit_sub,Ngrid+nIC))
    y_post_it = np.zeros((nit_sub,nmeasure))
    y_post=np.zeros(nmeasure)
    y_bg_it = np.zeros((nit_sub,nmeasure))
    y_bg=np.zeros(nmeasure)
    
    x_post_all_it[:,:nIC]=x_it[:,:nIC]
    x_post_all_it[:,nIC:]=x_post_vit
    
    for it in range(nit_sub):
        y_post_it[it,:]=np.dot(h_v,x_post_all_it[it,nIC:])  
        y_bg_it[it,:]=np.dot(h_v_all[:,:nBC],x_it[it,:nBC])
    
    y_post_it=y_post_it+y_bg_it
    
    for ti in range(nmeasure):
        y_post[ti]=np.mean(y_post_it[:,ti])
        y_bg[ti]=np.mean(y_bg_it[:,ti])
    
    ##################################################################################
    # SAVE MCMC output in a dataset and write to netcdf
    # Set up post-mcmc dataset
    
    props_temp=["birth", "death", "move", "sigma_y", "swap"]
    props = np.zeros((nIC1+len(props_temp)),dtype=object)
    accepts = np.zeros((nIC1+len(props_temp)))
    rejects = np.zeros((nIC1+len(props_temp)))
    for ii in range(nBC):
        props[ii]="bc"+str(ii)
    for jj in range(nfixed):
        props[jj+nBC]="fixed"+str(jj)
    props[nIC1-1] = "vary"
    props[nIC1:] = props_temp
    
    accepts[:nIC1]=accept
    rejects[:nIC1]=reject
    
    accepts[nIC1:] = [accept_birth,accept_death,
                             accept_move,accept_sigma_y, accept_swap]
    
    reject_swap = (nIt/2-accept_swap)
    rejects[nIC1:] = [reject_birth,reject_death,
                             reject_move,reject_sigma_y, reject_swap]
    
    #Do I need to store both x_it and x_post_vit. Can't I just store x_post_vit_all[nit,NgridIC]?
    if inv_type is 'evencorr':    
        y[timeindex_zero]=np.nan
        sigma_measure[timeindex_zero]=np.nan
        sigma_y_out[timeindex_zero,:]=np.nan
    
    
    post_mcmc = xray.Dataset({"x_it": (["nIt", "kICmax"],x_it),
                            "k_it": (["nIt"],k_it),                        
                            "x_post_vit": (["nIt", "Ngrid"],x_post_vit),                       
                            "regions_it": (["nIt", "Ngrid"],regions_it),                       
                            "plon_it": (["kmax","nIt"],plon_out),                       
                            "plat_it": (["kmax","nIt"],plat_out),                     
                            "sigma_y_it": (["nmeasure", "nIt"],sigma_y_out),
                            "sigma_measure": (["nmeasure"],sigma_measure),
                            "sigma_model_it": (["ydim2", "nIt"],sigma_model_out),
                            "R_indices": (["ydim1", "ydim2"],R_indices),
                            "pdf_p1_it": (["kICmax","nIt"],pdf_param1_out),
                            "pdf_p2_it": (["kICmax","nIt"], pdf_param2_out),                     
                            "y": (["nmeasure"], y),
                            "y_time": (["nmeasure"], y_time),
                            "y_site": (["nmeasure"], y_site),
                            "release_lons": (["sites"], release_lons),
                            "release_lats": (["sites"], release_lats),
                            "accepts": (["proposal"],
                            accepts),
                            "rejects": (["proposal"],
                            rejects),                          
                            "stepsize": (["nIC1"], tot_acc_x),
                            "stepsize_pdf_p1": (["nIC1"], tot_acc_p1),
                            "stepsize_pdf_p2": (["nIC1"], tot_acc_p2),
                            "stepsize_sigma_y": (["ydim2"], tot_acc_sigma_y),
                            "h_v_all": (["nmeasure","NgridIC"],
                            h_v_all), 
                            "q_ap": (["lat", "lon"],
                            q_ap0),
                            "dates": (["ndates"], [start_date,end_date]),
                            "measure_av": (["sites"], av_period),
                            "nIC": nIC,
                            "nfixed": nfixed},
                            coords={"lon":lon, "lat": lat})
    
    if inv_type in ('evencorr', 'correlated'):
        post_mcmc.update({'tau_it': (["numsites","nIt"], tau_out)})   
    
    post_mcmc.attrs["bc_basis_case"]=bc_basis_case
    post_mcmc.attrs["fp_basis_case"]=fp_basis_case
    post_mcmc.attrs["iterations"]=str(nIt)
    post_mcmc.attrs["burn-in"]=str(burn_in)
    post_mcmc.coords["proposal"]=props
    post_mcmc.coords["sites"]=sites
    # Also:
    #Attributes: Sites, av_period, basis_case, accepts and rejects names
    #output_directory="/home/ml12574/work/programs/Python/my_acrg/td_uncorr/"
    
    #Output files from tdmcmc_template.py stored in the form:
    # "output_" + network + "_" + species +  "_" + date + ".nc"
    
    fname=os.path.join(output_dir,
                        "output_" + network + "_" + species + "_" + start_date + ".nc")
    for key in post_mcmc.keys():
        post_mcmc[key].encoding['zlib'] = True
    post_mcmc.to_netcdf(path=fname, mode='w')

    return post_mcmc