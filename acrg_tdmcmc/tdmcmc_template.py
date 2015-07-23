# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:32:41 2015

Template script for running Trans-dimensional inversion

Uses hb_td_mcmc_acrg.f90

Relies on y=sum(h_v)*mean(q_v)##(x)
Solves for x

@author: ml12574
"""

import hbtdmcmc
import acrg_name_xray as name
import numpy as np
import acrg_agage as agage
import matplotlib.pyplot as plt
import pandas
import datetime as dt
from numba import jit
import time as run_time
import xray
from mpl_toolkits.basemap import Basemap
from acrg_grid import areagrid
import os
import argparse

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

def get_nsigma_y(fp_data_H,start_date, end_date, bl_period):
    
    sites = [key for key in fp_data_H.keys() if key[0] != '.']
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
        bl_start=d0
        mf_time_temp=fp_data_H3.time.values
        #mf_time_temp=fp_data_H[sites[si]].time.values
        mf_time_temp2=pandas.to_datetime(mf_time_temp)
        
        ntime_stn[si]=len(mf_time_temp)
        for ti in range(nsigma_max):
            bl_end=bl_start+dt.timedelta(days=bl_period)
            
            wh=np.where(np.logical_and(mf_time_temp2>=bl_start,
                                       mf_time_temp2<bl_end))
                                          
            if np.sum(wh) > 0:
                y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
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
    
    return R_indices, ydim1, ydim2
    

#%%
#############################################################
#parser = argparse.ArgumentParser(description='This is a demo script by Mark.')
#parser.add_argument("start", help="Start date string yyyy-mm-dd")                  
#parser.add_argument("end", help="End date sting yyyy-mm-dd")
#args = parser.parse_args()

#start_date = args.start
#end_date = args.end
    
sites=['MHD', 'TAC', 'RGL', 'TTA']
av_period=['2H', '2H', '2H', '2H']
species='CH4'
start_date = '2014-04-01'
end_date = '2014-05-01'
domain='EUROPE'


bl_period=10           # No. of days for which each sigma_y value applies
kmin=4       
kmax=250    
k_ap = 100   
nIt=50000        # of iterations
burn_in=50000     # of iterations to discard
nsub=100       # nsub=100=store every 100th iteration)
#sigma_y=0.02     # Model-measurement uncertainty, in % (e.g. 0.05)
                 
nbeta=8       # Number of parallel chains
beta= np.array((1.,1./2.,1./4.,1./8.,1./16.,1./32.,1./64.,1./128)) 


x_pdf=np.array([3,3])  # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
pdf_param1_pdf = 1
pdf_param2_pdf = 1

# Uniform
#pdf_param10=np.array([0.8,0.5])
#pdf_param20=np.array([1.25,2.])

# Gaussian/Lognormal
pdf_param10=np.array([1.,1.])
pdf_param20=np.array([0.2,0.2])

pdf_hyperparam1=pdf_param10*2
pdf_hyperparam2=pdf_param20*2

stepsize=0.9    # Stepsize for x_update in MCMC - defaults otherwise
stepsize_pdf_p1=0.1*0.
stepsize_pdf_p2=0.1*0.
############################################################
#%%
data = agage.get_obs(sites, species, start = start_date, end = end_date, average = av_period)


fp_all = name.footprints_data_merge(data, domain=domain, species=species)

fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case='transd')

fp_data_H2=name.bc_sensitivity(fp_data_H2, domain=domain)

#fp_data_H = name.filtering(fp_data_H2, ["pblh_gt_500"])
fp_data_H = name.filtering(fp_data_H2, ["pblh_gt_500", "daily_median"])

lat = np.asarray(fp_data_H[sites[0]].sub_lat)
lon = np.asarray(fp_data_H[sites[0]].sub_lon)
nlat=len(lat)
nlon=len(lon)
lonmin=np.min(lon)
lonmax=np.max(lon)
latmin=np.min(lat)
latmax=np.max(lat)
Ngrid = nlon*nlat  # Define underlying grid    

########### Stack fps, fluxes and obs #######################
y = []
y_site = []
y_time = []
H_bc=[]

for si, site in enumerate(sites):
    
    fp_data_H3 = fp_data_H[site].dropna("time", how="all")
    y.append(fp_data_H3.mf.values)     
    y_site.append([site for i in range(len(fp_data_H3.coords['time']))])
    y_time.append(fp_data_H3.coords['time'].values)
    H_bc.append(fp_data_H3.bc.values)
    sub_flux_temp = fp_data_H3.flux.sel(lon=slice(lonmin,lonmax), 
                                    lat=slice(latmin,latmax))
    
    if si ==0:
        H_fixed2=fp_data_H3.H
        H_vary2=fp_data_H3.sub_fp
        q_ap2=sub_flux_temp
              
    else:
        H_fixed2=xray.concat((H_fixed2,fp_data_H3.H), dim="time")    
        H_vary2=xray.concat((H_vary2,fp_data_H3.sub_fp),dim="time" )  
        q_ap2=xray.concat((q_ap2,sub_flux_temp), dim="time") 

if H_fixed2.dims[0] != "time":
    H_fixed2=H_fixed2.transpose()
    H_vary2=H_vary2.transpose("time","sub_lat","sub_lon")
    q_ap2=q_ap2.transpose("time","lat","lon")
    
H_fixed=H_fixed2.values
H_vary=H_vary2.values
q_ap=q_ap2.values
 
y = np.hstack(y)
y_site = np.hstack(y_site)
y_time = np.hstack(y_time)
H_bc = np.hstack(H_bc)

nmeasure=len(y)
h_v = np.zeros((nmeasure,Ngrid))
for ti in range(nmeasure):                                    
    h_v[ti,:] = np.ravel(H_vary[ti,:,:])   # Create sensitivty matrix spatially vectorised
    

#%%
############################################################
# Create IC  
#nBC=len(fp_data_H[sites[0]].region_bc)
nBC=1
nfixed = len(fp_data_H[sites[0]].region)

R_indices, ydim1, ydim2 = get_nsigma_y(fp_data_H,start_date, end_date, bl_period)
nIC=nBC+nfixed


h_agg0 = np.zeros((nmeasure,k_ap+nIC))
q_agg0 = np.zeros((k_ap+nIC))
x_agg=np.zeros((k_ap+nIC))+1.  

q_agg0[:nIC]=1.
h_agg0[:,0]=H_bc
h_agg0[:,nBC:nIC]=H_fixed

#%%
# Define prior model uncertainty
####################################################
    
sigma_model0 = np.zeros((ydim2)) 
y_error = np.zeros(nmeasure)

sigma_model0[:] = np.mean(y)*0.025   
sigma_measure = 0.002*y   # ONLY APPROPRIATE FOR METHANE
y_error[:] = sigma_model0[0]	

sigma_y0 = np.sqrt(y_error**2 + sigma_measure**2)

sigma_model_hparams = np.array([0.,2.*sigma_model0[0]])
stepsize_sigma_y=0.5
sigma_model_pdf = 1   # UNIFORM

#%%

# Define prior model and regions with uniform distribution
#######################################
kICmax=kmax+nIC              # nIC and kmax already defined at top of file

q_ap0=q_ap[0,:,:]
q_ap_v = np.ravel(q_ap0)
# Set up different starting nuclei locations for each chain 
plon=np.zeros((kmax,nbeta))
plat=np.zeros((kmax,nbeta))
regions_v=np.zeros((Ngrid,nbeta),dtype=np.uint16)
h_agg=np.zeros((nmeasure, kICmax,nbeta))
q_agg=np.zeros((kICmax,nbeta))
n0=np.zeros((nmeasure,nbeta))    
n0T=np.zeros((nbeta))
for ib in range(nbeta):

    plon[:k_ap,ib] = np.random.uniform(lonmin, lonmax, k_ap) # Lon locs of nuclei
    plat[:k_ap,ib] = np.random.uniform(latmin, latmax, k_ap) # Lat locs of nuclei


    region = np.zeros((nlat, nlon), dtype=np.uint16)
    regions0=closest_grid(region, lon, lat, plon[:k_ap,ib], plat[:k_ap,ib], \
            np.arange(0, k_ap, dtype=np.uint16))
    regions_v0 = np.ravel(regions0)
    regions_v[:,ib]=regions_v0.copy()+1

    for ri in range(k_ap):
        wh_ri = np.where(regions_v0 == ri)
        q_agg0[ri+nIC]=np.mean(q_ap_v[wh_ri])          
        for ti in range(nmeasure):
            h_agg0[ti,ri+nIC]=np.sum(h_v[ti,wh_ri])*q_agg0[ri+nIC]
  
    y_model = np.dot(h_agg0,x_agg) 
    n0_ap = y_model-y
    n0T_ap = np.sum((n0_ap/sigma_y0)**2)
  
  
    h_agg[:,:k_ap+nIC,ib] = h_agg0.copy()
    q_agg[:k_ap+nIC,ib]=q_agg0.copy()
    n0[:,ib]=n0_ap.copy()
    n0T[ib]=n0T_ap*1.
    
    
#################################


#%%
# MCMC Parameters
#########################################
nit_sub=nIt/nsub

t0=run_time.time()

k=np.zeros((nbeta),dtype=np.int)+k_ap

x=np.zeros((kICmax,nbeta))
pdf_param1 = np.zeros((2,nbeta))
pdf_param2 = np.zeros((2,nbeta))
sigma_y = np.zeros((nmeasure,nbeta))
sigma_model = np.zeros((ydim2,nbeta))
for ib in range(nbeta):  
    x[:k_ap+nIC,ib]=x_agg.copy()
    pdf_param1[:,ib]=pdf_param10
    pdf_param2[:,ib]=pdf_param20
    sigma_y[:,ib]=sigma_y0.copy()
    sigma_model[:,ib]=sigma_model0.copy()


#stepsize=mu_x*stepsize

sigma_clon = (lon[1]-lon[0])*8.
sigma_clat = (lat[1]-lat[0])*5.
sigma_bd=np.mean(x_agg[nIC:])*2.

#BEGIN ITERATIONS
##################################################
print 'Starting MCMC...'
startt = run_time.time()

## x_pdf needs to be array(2), 1st term for IC, 2nd for emissions
## kICmax = kmax+nIC

k_it, x_out, regions_out, plon_out, plat_out, sigma_y_out, n0T_out, \
pdf_param1_out, pdf_param2_out, q_agg_out, accept, reject, \
accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
accept_swap, accept_sigma_y, reject_sigma_y = hbtdmcmc.hbtdmcmc(
beta,k, x, h_agg,q_agg,y,n0,n0T, plon, plat, regions_v, 
pdf_param1, pdf_param2, lon,lat, h_v, q_ap_v, sigma_y, sigma_model, sigma_measure, 
R_indices, sigma_model_hparams, stepsize_sigma_y, sigma_model_pdf, 
sigma_clon, sigma_clat,  
lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
pdf_hyperparam1, pdf_hyperparam2, pdf_param1_pdf, pdf_param2_pdf, 
stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
nIt, nsub, nit_sub, nIC, 
nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2)

 

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
q_agg_it=np.transpose(q_agg_out)

for it in range(nit_sub):    
    for zz in range(k_it[it]):
        wh_reg = np.where(regions_it[it,:] == zz)
        x_post_vit[it,wh_reg] = x_it[it,zz+nIC]*q_agg_it[it,zz+nIC]

# THe nested loop above takes time - how can I speed it up?

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
    y_bg_it[it,:]=np.dot(h_v_all[:,:nIC],x_it[it,:nIC])

y_post_it=y_post_it+y_bg_it

for ti in range(nmeasure):
    y_post[ti]=np.mean(y_post_it[:,ti])
    y_bg[ti]=np.mean(y_bg_it[:,ti])

#PLOTTING

data=(x_post_mean-q_ap0)  # mol/m2/s
area=areagrid(lat,lon)
data=data*area*16.04*3600.*24.*365./1e9



lon_range = (lonmin, lonmax)
lat_range = (latmin, latmax)

m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')


lons, lats = np.meshgrid(lon,lat)
mapx, mapy = m(lons, lats)
data2=np.log10(x_post_mean)

fig = plt.figure(1,figsize=(8,8))
fig.add_axes([0.1,0.1,0.8,0.8])

m.drawcoastlines()
m.drawstates()
m.drawcountries() 

clevels = np.arange(0., 4.0e-8, 2.e-9)  

cs = m.contourf(mapx,mapy,x_post_mean, clevels, extend='both', cmap='YlGnBu')
cb = m.colorbar(cs, location='bottom', pad="5%")  


##############################
# Density plot of nuclei locations 
dlon=lon[1]-lon[0]
dlat=lat[1]-lat[0]
xedges=np.arange(lonmin,lonmax+2.,2.)
yedges=np.arange(latmin,latmax+2.,2.)
lon_range = (np.min(xedges), xedges[-2])
lat_range = (np.min(yedges), yedges[-2])
plat_out_v=np.ravel(plat_out)
plon_out_v=np.ravel(plon_out)

H, yedges, xedges = np.histogram2d(plat_out_v, plon_out_v, bins=(yedges, xedges))

m2 = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')

lona, latb = np.meshgrid(xedges[:-1],yedges[:-1])
mapa, mapb = m(lona, latb)

fig = plt.figure(3,figsize=(8,8))
ax = fig.add_subplot(111)

m2.drawcoastlines()
m2.drawcountries() 

hlevels=np.arange(np.min(H), np.max(H), 10.) 

ds=m2.contourf(mapa,mapb,H, hlevels)
db = m2.colorbar(ds, location='bottom')


data3=xrange_90/x_post_mean
fig = plt.figure(2,figsize=(8,8))
fig.add_axes([0.1,0.1,0.8,0.8])

m.drawcoastlines()
m.drawcountries() 
c2levels = np.arange(0., 4.0, 0.05) 
cs = m.contourf(mapx,mapy,data3, c2levels, extend='both', cmap='RdBu_r')
cb = m.colorbar(cs, location='bottom', pad="5%")  

################################################
# create a histogram by providing the bin edges (unequally spaced)
#
plt.figure(2)

bins = [0,100,200,300,400,500,600,700,800,900,1000]
# the histogram of the data with histtype='step'
n, bins, patches = plt.hist(k_it, bins, normed=1, histtype='bar', rwidth=0.8)


fig = plt.figure(4,figsize=(8,8))
fig.add_axes([0.1,0.1,0.8,0.8])

m.drawcoastlines()
m.drawstates()
m.drawcountries() 
levels = np.arange(-10., 10., 0.1) 
cs = m.contourf(mapx,mapy,data, levels, extend='both', cmap='RdBu_r')
cb = m.colorbar(cs, location='bottom', pad="5%")  
plt.show()

# GET COUNTRY DATA
c_object=name.get_country('EUROPE', ocean=True)
cds = xray.Dataset({'country': (['lat','lon'], c_object.country), 
                    'name' : (['ncountries'],c_object.name) },
                                    coords = {'lat': (c_object.lat),
                                    'lon': (c_object.lon)})
country = cds.country.sel(lon=slice(lonmin,lonmax), 
                                    lat=slice(latmin,latmax))

country_v=np.ravel(country)

name_uk = np.where(cds.name == 'UNITED KINGDOM')
uk_index = np.where(country_v == name_uk[0])
area_v=np.ravel(area)

x_uk = np.sum(x_post_v_mean[uk_index]*area_v[uk_index]) # in mol/s
x_uk = x_uk*16.04/1000.   # in kg/s
x_uk = x_uk*365.*24.*3600./1.e9    # in Tg/yr

uk_95 = np.sum(x_post_v_95[uk_index]*area_v[uk_index])
uk_05 = np.sum(x_post_v_05[uk_index]*area_v[uk_index])
uk_ap = np.sum(q_ap_v[uk_index]*area_v[uk_index])

uk_95 = uk_95*16.04/1000.*365.*24.*3600./1.e9   # in kg/s
uk_05 = uk_05*16.04/1000.*365.*24.*3600./1.e9 
uk_ap = uk_ap*16.04/1000.*365.*24.*3600./1.e9 
# SAVE MCMC output in a dataset and write to netcdf
# Set up post-mcmc dataset

props=["x_update", "birth", "death", "move", "sigma_y", "swap"]


post_mcmc = xray.Dataset({"x_it": (["nIt", "kICmax"],
                        x_it),
                        "k_it": (["nIt"],
                        k_it),
                        "q_agg_it": (["nIt", "kICmax"],q_agg_it),
                        "x_post_vit": (["nIt", "Ngrid"],
                        x_post_vit),
                        "regions_it": (["nIt", "Ngrid"],
                        regions_it),
                        "plon_it": (["kmax","nIt"],
                        plon_out),
                        "plat_it": (["kmax","nIt"],
                        plat_out),
                        "sigma_y_it": (["nmeasure", "nIt"],
                        sigma_y_out),
                        "n0T_it": (["nIt"],
                        n0T_out),
                        "y": (["nmeasure"], y),
                        "accepts": (["proposal"],
                        [accept,accept_birth,accept_death,
                         accept_move,accept_sigma_y, accept_swap]),
                        "rejects": (["proposal"],
                        [reject,reject_birth,reject_death,
                         reject_move,reject_sigma_y, (nIt/2-accept_swap)]),
                        "h_v_all": (["nmeasure","NgridIC"],
                        h_v_all), 
                        "q_ap": (["lat", "lon"],
                        q_ap0), 
                        "nIC": nIC,
                        "nfixed": nfixed},
                        coords={"proposal": props, "lon":lon, "lat": lat})


output_directory="/PATH/TO/OUTPUTS/"
#fname=os.path.join(output_directory,
#                    "output_" + start_date + ".nc")
#post_mcmc.to_netcdf(path=fname, mode='w')
