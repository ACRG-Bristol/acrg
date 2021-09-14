#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 11:17:54 2017

Make a real data 1-dimensional rjmcmc 

Idea is to break the emissions space up in time

For this first version just keep BC and fixed reigons separate from rest


@author: Mark Lunt
"""
from __future__ import print_function
import acrg_name as name
import numpy as np
import acrg_obs
import pandas
import datetime as dt
from numba import jit
import time as run_time
import xray
import os
import cams_reanalysis
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from acrg_tdmcmc import tdmcmc_post_process as process
import scipy
from acrg_grid import areagrid

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH")
    data_path = os.getenv("DATA_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg
    data_path = paths.data

#%%
# Commands for argument parsing from a shell script

#parser = argparse.ArgumentParser(description='This is a demo script by Mark.')
#parser.add_argument("start", help="Start date string yyyy-mm-dd")                  
#parser.add_argument("end", help="End date sting yyyy-mm-dd")
#parser.add_argument("out", help="Output date sting yyyy-mm")
##parser.add_argument("percent", help="Percentile for localness filter")
##parser.add_argument("diff", help="agl-asl mf difference")
#args = parser.parse_args()
#
#start_date = args.start
#end_date = args.end
#out_date = args.out

#############################################################
#%%
start_date = '2014-12-25'
end_date = '2015-02-05'
out_date = '2015-01-01'

sites=['RGL','TAC']
meas_period = ['2H','2H']  # Frequency to read in measurements 
heights=["90m","185m"]
av_period=None   # Frequency to average footprints and measuerements for inversion


fp_heights={}
fp_heights['MHD']= "10m"
fp_heights['TAC']= "185m"
fp_heights['RGL']= "90m"
fp_heights['TTA']= "222m"
fp_heights['BSD']= "250m"
fp_heights['HFD']= "100m"

inlets={}
inlets['MHD']= 100
inlets['TAC']= 185
inlets['RGL']= 100
inlets['TTA']= 222
inlets['BSD']= 250
inlets['HFD']= 100

species='CH4'
domain='EUROPE'
run_name="test"


lapse_dir = data_path + "/LPDM/vertical_profiles/"
output_directory="path/to/output/"
fp_basis_case = 'country_mask'
bc_basis_case = 'NESW'

emissions_name = 'ch4-naei-total'

sector_names=['ch4-wetlands-landuse-uk-only','ch4-naei-uk-only-other',
             'ch4-naei-uk-only-offshore','ch4-naei-uk-only-waste',
             'ch4-naei-uk-only-agric']
sector_short=['wetlands','other',
             'offshore','waste',
             'agric']             
countries = ["SW", "NW", "NE", "SE", "sub_land", "sub_sea", "FRANCE", "GERMANY", "BENELUX", 
             "NORWAY", "DENMARK", "IRELAND", "UK-wetlands","UK-other","UK-offshore",
             "UK-waste","UK-agric"]
#######################################################
# DO YOU WANT TO DO REVERSIBLE JUMP OR NOT?????
reversible_jump = True          # True = do reversible jump; False = don't
parallel_tempering = True      # True = do parallel tempering

# DO YOU WANT CORRELATED MEASUREMENTS OR NOT??
inv_type = 'corr'    # Options are 'uncorr', 'corr'

bl_period = 31      # No. of days for which each sigma_model value applies
kmin=2             # Minimum number of regions
kmax=30           # Maximum number of regions 
k_ap =5         # Starting number of regions
nIt=2000          # of iterations
burn_in=2000      # of discarded burn-in iterations 
nsub=50           # nsub=100=store every 100th iteration)
          
#####################################################
# PARALLEL TEMPERING PARAMETERS          
#nbeta=2            # Number of parallel chains - needs to be defined even if no Parallel tempering
#beta= np.array((1.,1./2.)) # Values of beta for tempered chains 
nbeta=8            # Number of parallel chains - needs to be defined even if no Parallel tempering
beta= np.array((1.,0.5, 0.25,0.125,0.0625, 1./32., 1./64., 1./128.)) # Values of beta for tempered chains 
#################################################################################
# DEFINE FORM OF PDFs FOR EMISSIONS, EMISSIONS UNCERTAINTIES AND MODEL UNCERTAINTY
x_pdf0=3        # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
pdf_param1_pdf0 = 1
pdf_param2_pdf0 = 1
sigma_model_pdf = 1   
tau_pdf=1 
##########################################################################
# Hyperparameters of emissions parameters
pdf_param10=1.         # Mean of lognormal or normal distribution
pdf_param20=0.4        # Std. of lognormal or normal distribution

pdf_p1_hparam10=0.8    # Lower bound of uniform distribution of pdf_param1
pdf_p1_hparam20=1.2    # Upper bound of uniform distribution of pdf_param1
pdf_p2_hparam10=0.05   # Lower bound of uniform distribution of pdf_param2
pdf_p2_hparam20=1.    # Upper bound of uniform distribution of pdf_param2

#pdf_param10=0.1         # Mean of lognormal or normal distribution
#pdf_param20=4.    
#pdf_p1_hparam10=0.04    # Lower bound of uniform distribution of pdf_param1
#pdf_p1_hparam20=1.    # Upper bound of uniform distribution of pdf_param1
#pdf_p2_hparam10=0.05   # Lower bound of uniform distribution of pdf_param2
#pdf_p2_hparam20=10.   

######################################################################
# Model-Measurement starting value and uncertainty parameters
sigma_model_ap = 10.   # Initial starting value of sigma_model (in same units as y, e.g ppb)
sigma_model_hparams = np.array([0.1*sigma_model_ap,5.*sigma_model_ap]) # upper and lower bounds of uniform dist. 

tau_ap = 12.
tau_hparams=np.array([1., 120.])      # upper and lower bounds of uniform dist.
stepsize_tau=5.

####################################################################
# DEFINE STEPSIZE FOR PROPOSAL DISTRIBUTIONS 
# TO TUNE INDIVIDUAL ELEMENTS SEE LINE 
stepsize=0.1   # Stepsize for proposal distirbution of x_update 
stepsize_pdf_p1=0.1  # Stepsize for proposal distirbution fof pdf_param1 
stepsize_pdf_p2=0.1   # Stepsize for proposal distirbution of pdf_param2 
stepsize_sigma_y=0.5  # Stepsize for proposal distribution of sigma_model

stepsize_chour = 5.    # Stepsize for hours for move
sigma_bd=1.        # Stepsize for change in x during birth step

# FILTERS
"""
# Include a list of filters you want to include here to filter the data
# Purpose is to remove observation times with potential biases

Options are:
 "daily_median": What it says
 "daytime": Only between 11:00 - 15:00 inclusive 
 "nighttime": Only b/w 23:00 - 03:00 inclusive
 "noon": Only 12:00 fp and obs used
 "pblh_gt_500":
 "pblh_gt_250": 
 "local_influence": Only keep times when localness is low
 "six_hr_mean":
 "ferry_loc": GAUGE-FERRY specific - Used to filter out dodgy ferry locations
 "ferry_mf": GAUGE-FERRY specific - Used to filter out dodg ferry 
 "ferry_fp_zero": GAUGE-FERRY specific 
 "lapse_rate": Based on atmopsheric stability
 "local_lapse": lBased on a combination of localness and atmospheric stability

"""
filters = ["local_lapse"]   

# END OF INPUTS SECTION
############################################################
#%%
# FUNCTIONS

def get_nsigma_y(fp_data_H, sites, nmeasure):
     
    """
    This version is based on modelled pollution peaks as means of splitting uncertainties.
    Could also use this to create a bias term for different times??
    """ 
    nsites = len(sites)      
    y_bl=np.zeros((nmeasure))
    
    pcs=[0.,25.,50.,75.,100]
    ngroups=len(pcs)
    
    nsigma=0
    ntime_stn=np.zeros((nsites))
    
    ydim1=0
    
    for si in range(nsites):
        
        fp_data_H3 = fp_data_H[sites[si]].dropna("time", how="all")        
        nsigma_stn=0  
        #mf_mod_temp=fp_data_H3.mf_mod.values
        mf_mod_temp=fp_data_H3.local_ratio.values
        #theta_slope = fp_data_H3.theta_slope.values
        ntime_stn[si]=len(mf_mod_temp)
                           
        for ti in range(ngroups-1):
            
            if ti == 0:
                wh = np.where(np.logical_and(mf_mod_temp >= np.percentile(mf_mod_temp,pcs[ti]),
                                            mf_mod_temp <= np.percentile(mf_mod_temp,pcs[ti+1])))
            else:
                wh = np.where(np.logical_and(mf_mod_temp > np.percentile(mf_mod_temp,pcs[ti]),
                                            mf_mod_temp <= np.percentile(mf_mod_temp,pcs[ti+1])))
#            if ti ==0:
#                wh = np.where(theta_slope<1.001)
#            else:
#                wh = np.where(theta_slope>=1.001)
                                           
                                                  
            if len(wh[0]) > 0:
                y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                nsigma_stn+=1
                
            n_obs = len(wh[0])
            if n_obs > ydim1:
                ydim1 = n_obs*1
                                      
    
        nsigma+=nsigma_stn
    
    # INDEX R
    print(ydim1, nsigma)
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2        

@jit(nopython=True)
def closest_point(region, hour, phour, pind):    
    hri=0
    for hr in hour:
        maxdist=1e12
        for pi in pind:
            dist= (hr - phour[pi])*(hr - phour[pi])
            if dist < maxdist:
                region[hri]=pi
                maxdist=dist
        hri+=1        
    return region
    
#%%
if reversible_jump == True:
    rjmcmc=1           # 1 = do reversible jump; any other number = don't
else:
    rjmcmc = 0
if parallel_tempering == True:    
    para_temp=1
else:
    para_temp=0
#########################################
# READ IN DATA AND FOOTPRINTS THEN MERGE

data = acrg_obs.get_obs(sites, species, start_date = start_date, end_date = end_date, average = meas_period, 
                      keep_missing=False, height=heights)                      


fp_all = name.footprints_data_merge(data, domain=domain, species=species, calc_bc=True, 
                                    height=fp_heights, emissions_name=emissions_name,
                                    interp_vmr_freq=None)                         

                        
fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case=fp_basis_case)


fp_data_H2=name.bc_sensitivity(fp_data_H2, domain=domain,basis_case=bc_basis_case)


###########################################################################
# CALCULATE DEGREE OF LOCALNESS FOR EACH FOOTPRINT
release_lons=np.zeros((len(sites)))
release_lats=np.zeros((len(sites)))
for si, site in enumerate(sites):
    release_lons[si]=fp_data_H2[site].release_lon[0].values
    release_lats[si]=fp_data_H2[site].release_lat[0].values
    dlon=fp_data_H2[site].lon[1].values-fp_data_H2[site].lon[0].values
    dlat=fp_data_H2[site].lat[1].values-fp_data_H2[site].lat[0].values
    local_sum=np.zeros((len(fp_data_H2[site].mf)))
    
    for ti in range(len(fp_data_H2[site].mf)):
        release_lon=fp_data_H2[site].release_lon[ti].values
        release_lat=fp_data_H2[site].release_lat[ti].values
        wh_rlon = np.where(abs(fp_data_H2[site].lon.values-release_lon) <= dlon/2.)        
        if len(wh_rlon[0] > 1):
            wh_rlon = [np.min(wh_rlon)]
        wh_rlat = np.where(abs(fp_data_H2[site].lat.values-release_lat) <= dlat/2.)
        if len(wh_rlat[0] > 1):
            wh_rlat = [np.min(wh_rlat)]
        #wh_rlon = np.where(abs(fp_data_H2[site].sub_lon.values-release_lon) < dlon/2.)
        #wh_rlat = np.where(abs(fp_data_H2[site].sub_lat.values-release_lat) < dlat/2.)
        local_sum[ti] = np.sum(fp_data_H2[site].fp[
        wh_rlat[0]-2:wh_rlat[0]+3,wh_rlon[0]-2:wh_rlon[0]+3,ti].values)/np.sum(
        fp_data_H2[site].fp[:,:,ti].values)#*topog_factor[site]
#        local_sum[ti] = np.sum(fp_data_H2[site].fp[
#        wh_rlat[0]-2:wh_rlat[0]+3,wh_rlon[0]-2:wh_rlon[0]+3,ti].values)
    
    local_da = xray.DataArray(local_sum, 
                              coords=[('time', fp_data_H2[site].coords['time'])])
    
    fp_data_H2[site]['local_ratio'] = local_da

    fp_data_H2[site].attrs['Domain']=domain
    fp_data_H2[site].attrs['Height']=fp_heights[site]
    fp_data_H2[site].attrs['inlet']=inlets[site]
    
    
    if site == "GAUGE-FERRY":
        slope_all=np.zeros((len(fp_data_H2[site].time)))
        std_all=np.zeros((len(fp_data_H2[site].time)))
        lapse_time = fp_data_H2[site].time
    else:
        lapse_fname = os.path.join(lapse_dir, 
                            site +"_vertical_profile.nc") 
                            
        lapse_ds = process.open_ds(lapse_fname)
        #lapse_ds['theta_slope'] = lapse_ds.theta_slope/np.mean(lapse_ds.theta_slope)
        lapse_ds['theta_slope'] = lapse_ds.theta_slope/2.6
        
    
    lapse_temp = lapse_ds.theta_slope.reindex_like(fp_data_H2[site],method="nearest")  
    fp_data_H2[site]["theta_slope"] = lapse_temp
    #fp_data_H2[site] = name.name.combine_datasets(fp_data_H2[site],lapse_ds)


if filters is not None:
       fp_data_H5 = name.filtering(fp_data_H2, filters, keep_missing=False) 
else:
        fp_data_H5 = fp_data_H2.copy()


if av_period == None:
    fp_data_H=fp_data_H5.copy()   
else:    
    fp_data_H = {}   
    for si, site in enumerate(sites):
        site_ds = fp_data_H5[site].resample(av_period[si], dim = "time")
        site_ds2= site_ds.dropna("time", how="all")
        fp_data_H[site] = site_ds2    
  


#%% 
# ADD UK sector emissions    
# 


if len(sector_names) > 0:

    for si,sec in enumerate(sector_names):
        rename_dict={}
        rename_dict['flux'] = "flux_"+sector_short[si]
        sector_ds = name.name.flux(domain, sec)
        sector_ds2=sector_ds.rename(rename_dict)
        
        if si == 0:
            all_sec_ds = sector_ds2
        else:
            
            
            # Fluxes not all of same time dimension so this doesn't work 
            all_sec_ds["flux_" + sector_short[si]] = (sector_ds2["flux_"+sector_short[si]].reindex_like(
                                                        all_sec_ds, method="nearest"))
        
    
    
    #fp_data_H[sites[0]] =  name.name.combine_datasets(fp_data_H[sites[0]], all_sec_ds)        
            
    for sec_si in sector_short: 
        #sec_temp = all_sec_ds["flux_"+sec_si].reindex_like(fp_data_H[site],method="nearest")  
        for site in sites:
            sec_temp = all_sec_ds["flux_"+sec_si].reindex_like(fp_data_H[site],method="nearest") 
            #mf_temp = fp_data_H[site].fp
            fp_data_H[site]["mf_mod_"+ sec_si]=(
            fp_data_H[site].fp*sec_temp).sum(["lat", "lon"])

#%%

# This section only works if UK is final country in basis functions

basis_ds=name.name.basis(domain, basis_case = fp_basis_case)
area=areagrid(basis_ds.lat.values,basis_ds.lon.values)
area_v=np.ravel(area)

basis_v=np.ravel(basis_ds.basis[:,:,0])

region_coords = np.arange(1,len(sector_names)+1)+len(fp_data_H[sites[0]].region)
sector_coords = np.arange(1,len(sector_names)+len(fp_data_H[sites[0]].region))
# Restructue H_fixed and mf_mod_sectors
for si, site in enumerate(sites):
    dum_ds = fp_data_H[site]
    h_sector = np.zeros((len(dum_ds.time),len(sector_coords)))
    if len(dum_ds.H[0,:]) == len(dum_ds.region):
        h_sector[:,:len(dum_ds.region)-1]=dum_ds.H[:,:-1].values
    else:
        h_sector[:,:len(dum_ds.region)-1]=dum_ds.H[:-1,:].values.T
    
    flux_site_sector = np.zeros((len(dum_ds.time),len(sector_coords)))
    
    for bb in range(len(dum_ds.region)-1):
        wh_basis = np.where(basis_v == bb+1)
        flux_temp = dum_ds.flux
        for ti in range(len(dum_ds.time)):
            if len(flux_temp[0,0,:]) == len(dum_ds.time):
                flux_area_v = np.ravel(flux_temp[:,:,ti])*area_v
            else:
                flux_area_v = np.ravel(flux_temp[ti,:,:])*area_v
            flux_site_sector[ti,bb] = np.sum(flux_area_v[wh_basis])
    
    if len(sector_names) > 0:
        for sc, sec in enumerate(sector_short): 
            h_sector[:,sc+len(dum_ds.region)-1] = dum_ds["mf_mod_"+sec].values
                 
            # NEED TO CHANGE THIS PART
            sec_temp = all_sec_ds["flux_"+sec].reindex_like(dum_ds,method="nearest")  
            for ti in range(len(dum_ds.time)):
                
                flux_temp = sec_temp[:,:,ti]
                flux_site_sector[ti,sc+len(dum_ds.region)-1] = np.sum(flux_temp*area)
    
    da_temp1 = xray.DataArray(h_sector, 
                              coords=[('time', dum_ds.time), ('sector', sector_coords)])
    da_temp2 = xray.DataArray(flux_site_sector, 
                              coords=[('time', dum_ds.time), ('sector', sector_coords)])
    fp_data_H[site]["H_sector"] = da_temp1
    fp_data_H[site]["flux_sector"] = da_temp2
    

#%%
###########################################################
# EVERYTHING NEEDS TO BE ARRAYS FOR MCMC
#STACK FPs, FLUXES AND OBS
y = []
y_site = []
y_time = []
y_error=[]
h_sector_raw=[]


for si, site in enumerate(sites):
          
    fp_data_H3 = fp_data_H[site].dropna("time", how="all")  
    attributes = [key for key in fp_data_H3.keys() if key[0] != '.']  
    y.append(fp_data_H3.mf.values)   
    h_sector_raw.append(fp_data_H3.mf_mod.values)   
    y_site.append([site for i in range(len(fp_data_H3.coords['time']))])
    y_time.append(fp_data_H3.coords['time'].values)
    
    if 'dmf' in attributes:    
        y_error.append(fp_data_H3.dmf.values)
    elif 'vmf' in attributes:   
        y_error.append(fp_data_H3.vmf.values)
    else:
        y_error.append(0.002*fp_data_H3.mf.values) # 0.002 only appropriate for methane
    if si ==0:
        H_bc2=fp_data_H3.H_bc
        H_sector2 = fp_data_H3.H_sector
        flux_sector2=fp_data_H3.flux_sector
                   
    else:
        H_bc2=xray.concat((H_bc2,fp_data_H3.H_bc), dim="time") 
        H_sector2=xray.concat((H_sector2,fp_data_H3.H_sector), dim="time")  
        flux_sector2=xray.concat((flux_sector2,fp_data_H3.flux_sector), dim="time")  
        
if H_bc2.dims[0] !="time":
    H_bc2=H_bc2.transpose()
    
H_bc=H_bc2.values
H_sector2=H_sector2.values
flux_sector=flux_sector2.values

y = np.hstack(y)
y_site = np.hstack(y_site)
y_time = np.hstack(y_time)
y_error=np.hstack(y_error)
h_sector_raw = np.hstack(h_sector_raw)

nmeasure=len(y)
y_hour = (y_time-np.datetime64(start_date))/np.timedelta64(1,'h')
max_hour=np.max(y_hour)
unique_hour = np.unique(y_hour)
hour_min = np.min(y_hour)
hour_max = np.max(y_hour)

y_day = (y_time-np.datetime64(start_date))/np.timedelta64(1,'D')
day_min = np.floor(np.min(y_day))
day_max = np.floor(np.max(y_day))

unique_time = pandas.to_datetime(np.unique(y_time))

#%%
############################################################
# Create IC  
nBC=len(fp_data_H[sites[0]].region_bc)
nfixed = len(fp_data_H[sites[0]].sector)
nIC=0


k_raw = nBC + nfixed
kmax_all = k_raw*kmax

numsites=len(sites)
nmeasure_site = np.zeros((numsites))
for ss, si in enumerate(sites):
    wh_site = np.ravel(np.where(y_site == si))
    nmeasure_site[ss]=len(wh_site)
 
print(nmeasure_site, nmeasure)

################################################
# CALCULATE INDICES OF Y CORRESPONDING TO DIFFERENT SIGMA_Ys   
#R_indices, ydim1, ydim2 = get_nsigma_y(fp_data_H,start_date, end_date, bl_period, sites, 
#                  nmeasure)
R_indices, ydim1, ydim2 = get_nsigma_y(fp_data_H, sites, nmeasure)


h_agg0 = np.zeros((nmeasure,k_raw*k_ap+nIC))
                
x_agg=np.zeros((k_raw*k_ap+nIC))+1.  

h_raw = np.zeros((nmeasure,k_raw))

h_raw[:,:nBC]=H_bc.copy()
h_raw[:,nBC:k_raw]= H_sector2.copy()


#%%

# Define prior model uncertainty
####################################################
sigma_model0 = np.zeros((ydim2)) 
sigma_model0[:]=sigma_model_ap
sigma_measure=y_error.copy()

stepsize_sigma_y_all=np.zeros((ydim2))+stepsize_sigma_y

error_structure=np.zeros((nmeasure))+1.

stepsize_sigma_y_all=sigma_model0*0.1
sigma_model_hparam1=sigma_model0*sigma_model_hparams[0]
sigma_model_hparam2=sigma_model0*sigma_model_hparams[1]
#%%
# DEFINE TAU PARAMS AND DELTATIME
if inv_type is 'corr':
    # DEFINE TAU PARAMS AND DELTATIME
    deltatime=np.zeros((nmeasure,nmeasure))+1.e12
    
    for ss, si in enumerate(sites):
        wh_site = np.ravel(np.where(y_site == si))
        nmeasure_site[ss]=len(wh_site)
        
        for whi in wh_site:
            #deltatime_site=np.zeros((nmeasure))
            #tdelta = np.absolute(y_time[wh_site]-y_time[whi]).astype('timedelta64[h]')
            tdelta = np.absolute(y_time[wh_site]-y_time[whi]).astype('timedelta64[m]')
            deltatime_site=tdelta/np.timedelta64(1, 'h')
            
            #deltatime_site[wh_site] = np.absolute(y_time[wh_site]-y_time[whi])
            deltatime[whi,wh_site] = deltatime_site


nsite_max=np.max(nmeasure_site)


#%%
# Define prior model and regions with uniform distribution
#######################################
kICmax=kmax_all+nIC              # nIC and kmax already defined at top of file

# Set up different starting nuclei locations for each chain 
phour=np.zeros((kmax,k_raw,nbeta))

regions_v=np.zeros((nmeasure,k_raw,nbeta),dtype=np.uint16)
h_agg=np.zeros((nmeasure, kICmax,nbeta))
n0=np.zeros((nmeasure,nbeta))    
n0T=np.zeros((nbeta))
sector_index0 = np.zeros((k_ap*k_raw), dtype=np.uint16)
k_sector0 = np.zeros((k_raw), dtype=np.uint16)   
sector_index = np.zeros((kmax_all,nbeta), dtype=np.uint16)
k_sector = np.zeros((k_raw,nbeta), dtype=np.uint16)   
 
k_sector0[:] = k_ap  
#phour0 = np.random.uniform(hour_min, hour_max, k_ap)

phour0 = np.floor(np.random.uniform(day_min, day_max, k_ap))*24.
  
   
region = np.zeros((nmeasure), dtype=np.uint16)
regions_v0=closest_point(region, y_hour, phour0, \
        np.arange(0, k_ap, dtype=np.uint16))

for si in range(k_raw):
    for ri in range(k_ap):
        wh_ri = np.where(regions_v0 == ri)      
        #for ti in range(nmeasure):
        h_agg0[wh_ri,nIC+(si*k_ap)+ri]=h_raw[wh_ri,si]
        sector_index0[si*k_ap+ri] = si    
  
y_model = np.dot(h_agg0,x_agg) 
n0_ap = y_model-y

for ib in range(nbeta):  
    sector_index[:k_ap*k_raw,ib]=sector_index0 + 1
    k_sector[:,ib]=k_sector0.copy()
    n0[:,ib]=n0_ap.copy()
    h_agg[:,nIC:nIC+k_raw*k_ap,ib] = h_agg0.copy()
    for si in range(k_raw):
        phour[:k_ap,si,ib] = phour0    
        regions_v[:,si,ib]=regions_v0 +1


     
   

#################################
     
#%%
# MCMC Parameters
#########################################
nit_sub=nIt/nsub

t0=run_time.time()

k=np.zeros((nbeta),dtype=np.int)+k_ap*k_raw

x=np.zeros((kICmax,nbeta))
sigma_model = np.zeros((ydim2,nbeta))
tau=np.zeros((numsites,nbeta))+tau_ap

for ib in range(nbeta):  
    x[:k_ap*k_raw+nIC,ib]=x_agg.copy()  
    sigma_model[:,ib]=sigma_model0.copy()

# nIC1 dimension Stuff
############################################################
nIC1=nIC+k_raw
pdf_param1 = np.zeros((kICmax,nbeta))
pdf_param2 = np.zeros((kICmax,nbeta))

#%%
#########################################
sigma_chour = stepsize_chour*1.

################################################
# TUNING OF INDIVIDUAL PARAMETER STEPSIZES AND UNCERTAINTIES

stepsize_all=np.zeros((nIC1))+stepsize
stepsize_pdf_p1_all=np.zeros((nIC1))+(stepsize_pdf_p1*pdf_param10)
stepsize_pdf_p2_all=np.zeros((nIC1))+(stepsize_pdf_p2*pdf_param20)

pdf_param1[:,:]=pdf_param10
pdf_param2[:,:]=pdf_param20

pdf_p1_hparam1=np.zeros((nIC1))+pdf_p1_hparam10
pdf_p1_hparam2=np.zeros((nIC1))+pdf_p1_hparam20

pdf_p2_hparam1=np.zeros((nIC1))+pdf_p2_hparam10
pdf_p2_hparam2=np.zeros((nIC1))+pdf_p2_hparam20

x_pdf = np.zeros((nIC1), dtype=np.int)+x_pdf0

pdf_param1_pdf=pdf_param1_pdf0
pdf_param2_pdf=pdf_param2_pdf0

##################################################
print('Starting MCMC...')
startt = run_time.time()
if inv_type is 'corr':
    from acrg_tdmcmc import tdmcmc_time_corr
    
    k_it, x_out, regions_out, phour_out, sigma_y_out, sigma_model_out, n0T_out, \
    pdf_param1_out,pdf_param2_out, k_sector_out, sector_index_out, tau_out, accept, reject, \
    accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
    accept_swap, reject_swap, accept_sigma_y, reject_sigma_y, \
    accept_tau, reject_tau, \
    step_out, step_p1_out, step_p2_out, step_sigma_y_out, \
    step_tau_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
    accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
    accept_sigma_y_all, reject_sigma_y_all, accept_tau_all, \
    reject_tau_all = tdmcmc_time_corr.transd_corr.hbtdmcmc(beta,k, x,
    h_agg,y,n0, phour, regions_v, 
    pdf_param1, pdf_param2, y_hour, h_raw, sigma_model, sigma_measure, error_structure,
    R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y_all, sigma_model_pdf, 
    tau, tau_hparams, stepsize_tau, tau_pdf, deltatime,
    sigma_chour, k_sector, sector_index, rjmcmc, para_temp, nmeasure_site, nsite_max, 
    hour_min, hour_max, sigma_bd, kmin, x_pdf, burn_in, 
    pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
    pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
    nIt, nsub, nit_sub, nIC, nbeta, kmax, kICmax, nmeasure,
    ydim1, ydim2, numsites,nIC1,k_raw, kmax_all)


elif inv_type is 'uncorr':
    from acrg_tdmcmc import tdmcmc_time_uncorr
    
    k_it, x_out, regions_out, phour_out, sigma_model_out,sigma_y_out, \
    n0T_out,pdf_param1_out,pdf_param2_out, k_sector_out, sector_index_out, accept, reject, \
    accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
    accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, \
    stepsize_x_out, stepsize_p1_out, stepsize_p2_out, \
    stepsize_sigma_y_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
    accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
    accept_sigma_y_all, reject_sigma_y_all = tdmcmc_time_uncorr.hbtdmcmc(
    beta,k, x, h_agg,y,n0, phour, regions_v, 
    pdf_param1, pdf_param2, y_hour, h_raw, sigma_model, sigma_measure, 
    R_indices, sigma_model_hparams, stepsize_sigma_y_all, sigma_model_pdf, 
    sigma_chour, k_sector, sector_index, rjmcmc, para_temp,
    hour_min, hour_max, sigma_bd, kmin, x_pdf, burn_in, 
    pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
    pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
    nIt, nsub, nit_sub, nIC, 
    nbeta, kmax, kICmax, nmeasure, ydim1, ydim2, nIC1, k_raw, kmax_all)   
    
    
endt=run_time.time()

print('Finished MCMC in ', endt-startt)

print('Beginning post processing')

# Need to somehow map the output back onto a regular time grid I think. 


x_post_temp=np.zeros((nit_sub,k_raw,nmeasure))

# Having run MCMC then need to map regions back onto grid
sigma_y_mean=np.zeros((nmeasure))

for yi in range(nmeasure):
    sigma_y_mean[yi] = np.mean(sigma_y_out[yi,:])

k_sector_it = np.transpose(k_sector_out)       
sector_index_it = np.transpose(sector_index_out)-1    
regions_it=np.transpose(regions_out)-1   
#regions_it=regions_out.transpose(2,0,1)-1
x_it=np.transpose(x_out)

for it in range(nit_sub):
    for si in range(k_raw):  
        wh_si = np.where(sector_index_it[it,:k_it[it]] == si)[0]
        for zz in range(k_sector_it[it,si]):
            wh_reg = np.where(regions_it[it,si,:] == zz)            
            z2 = wh_si[zz]
            x_post_temp[it,si,wh_reg] = x_it[it,z2+nIC]

endt3 = run_time.time()
print(endt3-endt)

#%%
x_post_vit = np.zeros((nit_sub, k_raw, len(unique_hour)))
flux_sector_post = np.zeros((len(unique_hour),k_raw-nBC))

for ti in range(len(unique_hour)):
    
    wh = np.where(y_hour == unique_hour[ti])[0]
    x_post_vit[:,:,ti] = x_post_temp[:,:,wh[0]]
    
    flux_sector_post[ti,:]=flux_sector[wh[0],:]


#x_post_vit = x_post_temp[:,:nmeas_max]

x_post_v_mean = np.mean(x_post_vit, axis=0)
x_post_v_95 = np.percentile(x_post_vit,95,axis=0)
x_post_v_05 = np.percentile(x_post_vit,5,axis=0)

#x_post_mean = np.reshape(x_post_v_mean, (nlat,nlon))
xrange_90v=x_post_v_95-x_post_v_05
#xrange_90=np.reshape(xrange_90v, (nlat,nlon))

print('Everything done')
endt2 = run_time.time()
print(endt2-endt)

##########################################
# y_post
y_post_it = np.zeros((nit_sub,nmeasure))
for it in range(nit_sub):
    #if nBC > 0:
    #    y_post_it[it,:] = (np.dot(h_agg0[:,:nBC],x_it[it,:nBC]) + 
    #                        h_raw*x_post_temp[it,:])
    #else:
    for si in range(k_raw):
        y_post_it[it,:] = y_post_it[it,:] + h_raw[:,si]*x_post_temp[it,si,:]

y_post=np.mean(y_post_it,axis=0)



rmean_len = len(unique_hour)-20
running_mean= np.zeros((len(unique_hour)-20))


for ti in range(rmean_len):
    running_mean[ti] = np.mean(x_post_v_mean[-1,ti:ti+20])

props_temp=["birth", "death", "move", "sigma_y", "swap"]
#props_temp=["birth", "death", "move", "sigma_y", "tau", "swap"]
#strs = ["" for ii in range(nIC1)]
props = np.zeros((nIC1+len(props_temp)),dtype=object)
accepts = np.zeros((nIC1+len(props_temp)))
rejects = np.zeros((nIC1+len(props_temp)))
for ii in range(nBC):
    props[ii]="bc"+str(ii)
    
#countries = ["SW", "NW", "NE", "SE", "sub_land", "sub_sea", "FRANCE", "GERMANY", "BENELUX", 
#             "NORWAY", "DENMARK", "IRELAND", "UNITED KINGDOM"]

    
for jj in range(nfixed):
    #props[jj+nBC]="fixed"+str(jj)
    props[jj+nBC]=countries[jj]
#props[nIC1-1] = "vary"
props[nIC1:] = props_temp

accepts[:nIC1]=accept
rejects[:nIC1]=reject

accepts[nIC1:] = [accept_birth,accept_death,
                         accept_move,accept_sigma_y, accept_swap]

rejects[nIC1:] = [reject_birth,reject_death,
                         reject_move,reject_sigma_y, reject_swap]
                         

fp_heights_out=[]                         
for si,site in enumerate(sites):
    fp_heights_out.append(fp_heights[site])

#%%  
post_mcmc = xray.Dataset({"x_it": (["nIt", "kICmax"],x_it),
                        "k_it": (["nIt"],k_it),                        
                        "x_post_vit": (["nIt", "k_raw", "time"],x_post_vit),                       
                        "regions_it": (["nIt", "k_raw", "nmeasure"],regions_it),                       
                        "phour_it": (["kmax", "k_raw", "nIt"],phour_out),                                                               
                        "sigma_y_it": (["nmeasure", "nIt"],sigma_y_out),
                        "sigma_measure": (["nmeasure"],sigma_measure),
                        "sigma_model_it": (["ydim2", "nIt"],sigma_model_out),
                        "R_indices": (["ydim1", "ydim2"],R_indices), 
                        "pdf_p1_it": (["kICmax","nIt"],pdf_param1_out),
                        "pdf_p2_it": (["kICmax","nIt"], pdf_param2_out), 
                        "y": (["nmeasure"], y),
                        "y_time": (["nmeasure"], y_time),
                        "y_site": (["nmeasure"], y_site),
                        "y_post_it": (["nIt", "nmeasure"], y_post_it),
                        "release_lons": (["sites"], release_lons),
                        "release_lats": (["sites"], release_lats),
                        "accepts": (["proposal"],accepts),
                        "rejects": (["proposal"],rejects),
                        "h_raw": (["nmeasure","k_raw"],h_raw), 
                        "dates": (["ndates"], [start_date,end_date]),
                        "heights": (["sites"],heights),
                        "fp_heights": (["sites"],fp_heights_out),
                        "flux_sector": (["time", "sector"], flux_sector_post),
                        #"PBLH":  (["nmeasure"], pblh),
                        #"wind_speed":  (["nmeasure"], wind_speed),
                        #"lapse_rate":  (["nmeasure"], lapse_rate),
                        #"local_ratio":  (["nmeasure"], local_ratio),
                        "nIC": nIC,
                        "nBC": nBC,
                        "nfixed": nfixed},
                        coords={"time":unique_time})
                        #coords={"lon":lon, "lat": lat, "time":unique_time})

if av_period != None:
    post_mcmc["measure_av"]= (("sites"), av_period)
    
if inv_type == 'corr':
    post_mcmc['tau_it'] = (("sites", "nIt"), tau_out)


post_mcmc.attrs["bc_basis_case"]='NESW'
post_mcmc.attrs["fp_basis_case"]=fp_basis_case
post_mcmc.attrs["iterations"]=str(nIt)
post_mcmc.attrs["burn-in"]=str(burn_in)
post_mcmc.coords["proposal"]=props
post_mcmc.coords["sites"]=sites
post_mcmc.coords["sector"]=countries
post_mcmc.attrs['Filters'] = filters
post_mcmc.attrs['Parallel_tempering'] = str(parallel_tempering)
post_mcmc.attrs['Inversion_type'] = inv_type

#Output files from tdmcmc_template.py stored in the form:
# "output_" + network + "_" + species +  "_" + date + ".nc"

fname=os.path.join(output_directory,
                    "output_" + run_name + "_" + species + "_" + out_date + ".nc")
for key in post_mcmc.keys():
        post_mcmc[key].encoding['zlib'] = True                    
post_mcmc.to_netcdf(path=fname, mode='w')
