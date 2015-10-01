# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:32:41 2015

Template script for running fully correlated Trans-dimensional inversion

Uses acrg_full_hbtdmcmc.f90


@author: ml12574
"""

import full_hbtdmcmc
import acrg_name as name
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
import scipy.special as special

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

def get_nsigma_y(fp_data_H,start_date, end_date, bl_period, nmeasuremax):
    
    sites = [key for key in fp_data_H.keys() if key[0] != '.']
    nsites = len(sites)    
    
    d0=pandas.to_datetime(start_date)
    d1=pandas.to_datetime(end_date)
    delta = d1 - d0
    ndays = delta.days
       
    y_bl=np.zeros((nmeasuremax))
    
    nsigma=0
    nsigma_max = np.int(np.ceil(ndays/np.float(bl_period)))
    
    av_period=2
    obs_per_day=24/av_period
    ydim1 = bl_period*obs_per_day    
    
    R_indices = np.zeros((ydim1,nsigma_max), dtype=np.uint16) 
    for ii in range(nsigma_max):
        R_indices[:,ii]=np.arange(ydim1)+ydim1*ii+1
        
    ydim2=nsigma_max*1        
    
    return R_indices, ydim1, ydim2
    
   
    

#%%
#############################################################
#parser = argparse.ArgumentParser(description='This is a demo script by Mark.')
#parser.add_argument("start", help="Start date string yyyy-mm-dd")                  
#parser.add_argument("end", help="End date sting yyyy-mm-dd")
#args = parser.parse_args()

#start_date = args.start
#end_date = args.end
    
#sites=['MHD', 'TAC', 'RGL']
#av_period=['2H', '2H', '2H']
sites=['MHD']
av_period=['2H']
#sites=['MHD', 'TAC', 'RGL', 'TTA']
#av_period=['2H', '2H', '2H', '2H']
species='CH4'
start_date = '2014-04-01'
end_date = '2014-05-01'
domain='EUROPE'

"""
KEEP PLUGGING AWAY AT THIS TILL I DISCOVER WHAT IT IS THAT MAKES IT NOT WORK!!
"""

rjmcmc=1

bl_period=30           # No. of days for which each sigma_y value applies
kmin=4       
kmax=250    
k_ap = 100   
nIt=1000       # of iterations
burn_in=1000     # of iterations to discard
nsub=100      # nsub=100=store every 100th iteration)
#sigma_y=0.02     # Model-measurement uncertainty, in % (e.g. 0.05) DEFINED LATER
                 
#nbeta=8       # Number of parallel chains
#beta= np.array((1.,1./2.,1./4.,1./8.,1./16.,1./32.,1./64.,1./128)) 

#beta=np.array((1., np.exp(-0.125),np.exp(-0.25),np.exp(-0.5),
#               np.exp(-1.),np.exp(-2.), np.exp(-4.), np.exp(-8.)))

nbeta=4       # Number of parallel chains
beta= np.array((1.,1./2.,1./512.,1./64.)) 


#x_pdf=np.array([3,3])  # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
#pdf_param1_pdf = 1
#pdf_param2_pdf = 1

# Uniform
#pdf_param10=np.array([0.8,0.5])
#pdf_param20=np.array([1.25,2.])

# Gaussian/Lognormal
#pdf_param10=np.array([1.,1.])
#pdf_param20=np.array([0.2,0.5])

#pdf_p1_hparam1=np.array([0.5,0.5])
#pdf_p1_hparam2=np.array([2.,2.])
#pdf_p2_hparam1=np.array([0.05,0.05])
#pdf_p2_hparam2=np.array([0.5,1.])

pdf_param10=1.
pdf_param20=0.2
pdf_p1_hparam10=0.5
pdf_p1_hparam20=2.
pdf_p2_hparam10=0.05
pdf_p2_hparam20=0.5

x_pdf0=3        # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
pdf_param1_pdf0 = 1
pdf_param2_pdf0 = 1


stepsize=0.5   # Stepsize for x_update in MCMC - defaults otherwise
stepsize_pdf_p1=0.1*0.
stepsize_pdf_p2=0.05

deltatime=2.
tau0 = 8.

############################################################
#%%
data = agage.get_obs(sites, species, start = start_date, end = end_date, average = av_period)


fp_all = name.footprints_data_merge(data, domain=domain, species=species, calc_bc=True,
                                    full_corr=True)

fp_data_H2 = name.fp_sensitivity(fp_all, domain=domain, basis_case='transd')

fp_data_H2=name.bc_sensitivity(fp_data_H2, domain=domain)

fp_data_H = name.filtering(fp_data_H2, ["pblh_gt_500"], full_corr=True)
#fp_data_H = name.filtering(fp_data_H2, ["pblh_gt_500", "daily_median"])

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
y_all = []
y_error=[]
y_site = []
y_time = []
#H_bc=[]

for si, site in enumerate(sites):
    
    #fp_data_H3 = fp_data_H[site].dropna("time", how="all")
    fp_data_H3 = fp_data_H[site]
    
    attributes = [key for key in fp_data_H3.keys() if key[0] != '.']  
    
    y_all.append(fp_data_H3.mf.values)     
    y_site.append([site for i in range(len(fp_data_H3.coords['time']))])
    y_time.append(fp_data_H3.coords['time'].values)


    if 'dmf' in attributes:    
        y_error.append(fp_data_H3.dmf.values)
    elif 'vmf' in attributes:   
        y_error.append(fp_data_H3.vmf.values)
    else:
        y_error.append(0.002*fp_data_H3.mf.values)
    
    #H_bc.append(fp_data_H3.bc.values)
    sub_flux_temp = fp_data_H3.flux.sel(lon=slice(lonmin,lonmax), 
                                    lat=slice(latmin,latmax))
    
    if si ==0:
        H_fixed2=fp_data_H3.H
        H_vary2=fp_data_H3.sub_fp
        q_ap2=sub_flux_temp
        H_bc2=fp_data_H3.H_bc
              
    else:
        H_fixed2=xray.concat((H_fixed2,fp_data_H3.H), dim="time")    
        H_vary2=xray.concat((H_vary2,fp_data_H3.sub_fp),dim="time" )  
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
 
y_all = np.hstack(y_all)

y_error=np.hstack(y_error)
y_site = np.hstack(y_site)
y_time = np.hstack(y_time)
#H_bc = np.hstack(H_bc)

q_ap0=q_ap[0,:,:].copy()
q_ap_v = np.ravel(q_ap0)

#################################################
# Define obs only where finite data exists
  
wh_temp=np.where(np.isfinite(y_all))
timeindex_nonzero=wh_temp[0]
z=y_all[timeindex_nonzero]

sigma_measure = y_error[timeindex_nonzero]
##############################################
# Define data dimensions
numsites=len(sites)
nmeasuremax= len(y_all)
nmeasuretotal= nmeasuremax/numsites
nmeasure = len(z)


tindex_zero_temp = np.arange(nmeasuremax)
timeindex_zero=np.delete(tindex_zero_temp, timeindex_nonzero)

#############################################

h_v = np.zeros((nmeasuremax,Ngrid))
for ti in range(nmeasuremax):                                    
    h_v[ti,:] = np.ravel(H_vary[ti,:,:])*q_ap_v   # Create sensitivty matrix spatially vectorised
    

# If these lines are set then correlations revert to small - i.e tau=2 rho=40.
# Solution map looks odd, very jagged and not smooth
# Think that means it's not an acceptable approach.    
#for ti in timeindex_zero:
#    h_v[ti,:] = 0.
#    H_bc[ti,:] = 0.
#    H_fixed[ti,:] = 0.

#%%
############################################################
# Create IC  
nBC=len(fp_data_H[sites[0]].region_bc)
#nBC=1
#nBC=numsites
nfixed = len(fp_data_H[sites[0]].region)

R_indices, ydim1, ydim2 = get_nsigma_y(fp_data_H,start_date, end_date, bl_period, nmeasuremax)
nIC=nBC+nfixed


h_agg0 = np.zeros((nmeasuremax,k_ap+nIC))
#q_agg0 = np.zeros((k_ap+nIC))
x_agg=np.zeros((k_ap+nIC))+1.  

#h_bc_temp=np.asarray(H_bc)
#for ii in range(nBC):
#    h_agg0[nmeasuretotal*ii:nmeasuretotal*(ii+1),ii]=h_bc_temp[ii,:]



#q_agg0[:nIC]=1.
#h_agg0[:,0]=H_bc*0.91
h_agg0[:,:nBC]=H_bc*0.95
h_agg0[:,nBC:nIC]=H_fixed


#%%
# Define prior model uncertainty
####################################################
    
sigma_yt0 = np.zeros((ydim2)) 
sigma_ys0 = np.zeros((numsites))

sigma_yt0[:] = np.mean(z)*0.002   
sigma_ys0[:] = np.mean(z)*0.005   

sigma_yt_hparams = np.array([0.1*sigma_yt0[0],5.*sigma_yt0[0]])
sigma_ys_hparams = np.array([0.1*sigma_ys0[0],5.*sigma_ys0[0]])
stepsize_sigma_yt=0.1
stepsize_sigma_ys=0.1
sigma_yt_pdf = 1   # UNIFORM
sigma_ys_pdf = 1   # UNIFORM

tau_pdf=1 # UNIFORM
stepsize_tau=2.
tau_hparams=np.array([0., 24.])

stepsize_y=0.01

#%%

# Define prior model and regions with uniform distribution
#######################################
kICmax=kmax+nIC              # nIC and kmax already defined at top of file

# Set up different starting nuclei locations for each chain 
plon=np.zeros((kmax,nbeta))
plat=np.zeros((kmax,nbeta))
regions_v=np.zeros((Ngrid,nbeta),dtype=np.uint16)
h_agg=np.zeros((nmeasuremax, kICmax,nbeta))
y=np.zeros((nmeasuremax,nbeta))    
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
        for ti in range(nmeasuremax):
            h_agg0[ti,ri+nIC]=np.sum(h_v[ti,wh_ri])
  
    y[:,ib] = np.dot(h_agg0,x_agg) 
    
    h_agg[:,:k_ap+nIC,ib] = h_agg0.copy()
 

#################################
#%%
# Define spatial correlation
# Set up a matern covariance function

# Create Matern covariance for spatial correlation

nu0 = 0.50
rho0 = 100.
distance = np.zeros((numsites,numsites))
if numsites==1:
    distance=0.
    space_corr=1.
elif numsites==2:
    
    distance[0,1]=250.
    distance[1,0]=250.
    
    space_corr=np.zeros((2,2))
    for ii in range(2):
        space_corr[ii,ii]=1.
        
elif numsites==3:
    mhd_dist = [0.,700.,470.]
    tac_dist = [700.,0.,250.]
    rgl_dist = [470.,250.,0.]
    
    distance[0,:] = mhd_dist
    distance[1,:] = tac_dist
    distance[2,:] = rgl_dist
elif numsites==4:

    mhd_dist = [0.,700.,470.,525.]
    tac_dist = [700.,0.,250.,515.]
    rgl_dist = [470.,250.,0.,495.]
    tta_dist = [525.,515.,495.,0.]
    
    #MHD-TAC = 700 km
    # MHD-RGL = 470 km
    # MHD-TTA = 525 km
    #TAC-RGL= 250 km
    #TAC-TTA= 515 km
    #RGL-TTA = 495 km
    
    distance[0,:] = mhd_dist
    distance[1,:] = tac_dist
    distance[2,:] = rgl_dist
    distance[3,:] = tta_dist
    
# Make spatial correlation function exponential for now!!   
    arg = np.sqrt(2*nu0)*distance/rho0
    
    #call gamma(nu_double, ga_double)
    
    ga=special.gamma(nu0)
    
    space_corr_temp=np.zeros((numsites,numsites))
    space_corr_temp2=np.zeros((numsites,numsites))
    for ii in range(numsites): 
        for jj in range(numsites):
                arg_dum = arg[ii,jj]
                if (arg_dum==0):
                    space_corr_temp[ii,jj] = 1.
                else:
                    k_arg_dum = special.kv(0,arg_dum)
                    #k_arg_dum = k_arg_dum*np.exp(-1.*arg[ii,jj])
                    space_corr_temp[ii,jj] = (1./(2.**(nu0-1)*ga)*(
                                        np.sqrt(2*nu0)*distance[ii,jj]
                                        /rho0)**nu0*k_arg_dum)
                    space_corr_temp2[ii,jj] = np.exp(-1.*distance[ii,jj]/rho0)
    
    			
    space_corr=np.linalg.inv(space_corr_temp)
    
#sign,logdet = np.linalg.slogdet(space_corr_temp)
#detval_space_corr = sign*logdet


rho_pdf=1
nu_pdf=1
rho_hparams=np.array([0., 1000.])
nu_hparams=np.array([0.499,2.5])

stepsize_rho=50.
stepsize_nu=0.2*0.

#%%
# MCMC Parameters
#########################################
nIC1 = nIC+1
nit_sub=nIt/nsub

t0=run_time.time()

k=np.zeros((nbeta),dtype=np.int)+k_ap

x=np.zeros((kICmax,nbeta))
pdf_param1 = np.zeros((nIC1,nbeta))
pdf_param2 = np.zeros((nIC1,nbeta))
sigma_ys = np.zeros((numsites,nbeta))
sigma_yt = np.zeros((ydim2,nbeta))
tau=np.zeros((nbeta))
rho=np.zeros((nbeta))
nu=np.zeros((nbeta))
for ib in range(nbeta):  
    x[:k_ap+nIC,ib]=x_agg.copy()
    #pdf_param1[:,ib]=pdf_param10
    #pdf_param2[:,ib]=pdf_param20
    sigma_ys[:,ib]=sigma_ys0.copy()
    sigma_yt[:,ib]=sigma_yt0.copy()
    
    tau[ib]=tau0
    rho[ib]=rho0
    nu[ib]=nu0


#stepsize=mu_x*stepsize

sigma_clon = (lon[1]-lon[0])*12.
sigma_clat = (lat[1]-lat[0])*8.
sigma_bd=np.mean(x_agg[nIC:])*2.


timeindex_nonzero=timeindex_nonzero+1


# nIC1 dimension Stuff
############################################################

stepsize_all=np.zeros((nIC+1))+stepsize
stepsize_pdf_p1_all=np.zeros((nIC+1))+(stepsize_pdf_p1*pdf_param10)
stepsize_pdf_p2_all=np.zeros((nIC+1))+(stepsize_pdf_p2*pdf_param20)
stepsize_all[:nBC]=stepsize_all[:nBC]/100.
stepsize_pdf_p1_all[:nBC]=stepsize_pdf_p1_all[:nBC]/10.
stepsize_pdf_p2_all[:nBC]=stepsize_pdf_p2_all[:nBC]/10.
#stepsize_all[-1]=stepsize_all[-1]*1.2

pdf_param1[:,:]=pdf_param10
pdf_param2[:,:]=pdf_param20

pdf_p1_hparam1=np.zeros((nIC1))+pdf_p1_hparam10
pdf_p1_hparam2=np.zeros((nIC1))+pdf_p1_hparam20

pdf_p2_hparam1=np.zeros((nIC1))+pdf_p2_hparam10
pdf_p2_hparam2=np.zeros((nIC1))+pdf_p2_hparam20

x_pdf = np.zeros((nIC1), dtype=np.int)+x_pdf0
pdf_param1_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param1_pdf0
pdf_param2_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param2_pdf0

#BEGIN ITERATIONS
##################################################
print 'Starting MCMC...'
startt = run_time.time()

## x_pdf needs to be array(2), 1st term for IC, 2nd for emissions
## kICmax = kmax+nIC

### This one based on hb_td_mcmc_acrg.f90

k_it, x_out, regions_out, plon_out, plat_out, n0T_out, y_out, \
sigma_yt_out, sigma_ys_out, tau_out, pdf_param1_out, pdf_param2_out, rho_out, nu_out,  \
accept, reject, accept_birth, reject_birth, accept_death, reject_death, \
accept_move, reject_move, accept_swap, reject_swap, accept_sigma_yt, reject_sigma_yt, \
accept_sigma_ys, reject_sigma_ys, accept_tau, reject_tau, accept_rho, reject_rho, \
accept_y, reject_y = full_hbtdmcmc.transd_inv.hbtdmcmc(
beta,k, x, h_agg,y, z, plon, plat, regions_v, 
pdf_param1, pdf_param2, lon,lat, h_v, sigma_yt, sigma_ys, sigma_measure,  
R_indices, timeindex_nonzero, sigma_clon, sigma_clat, tau, deltatime, distance,
pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, sigma_yt_hparams, sigma_ys_hparams, tau_hparams, 
stepsize_all, stepsize_pdf_p1_all, stepsize_pdf_p2_all,  stepsize_sigma_yt, stepsize_sigma_ys, stepsize_tau, stepsize_y, 
x_pdf, pdf_param1_pdf, pdf_param2_pdf, sigma_yt_pdf, sigma_ys_pdf, tau_pdf, 
rho, nu, rho_hparams, nu_hparams, stepsize_rho, stepsize_nu, rho_pdf, nu_pdf, rjmcmc,
lonmin, lonmax, latmin, latmax, kmin, sigma_bd, burn_in, nIt,nmeasuretotal, nsub, nIC, nit_sub, 
nmeasure, nmeasuremax, ydim1, ydim2, Ngrid, nlon, kICmax,nlat, nbeta, kmax,numsites, nIC1)

endt=run_time.time()

print 'Finished MCMC in ', endt-startt


print 'Beginning post processing'
x_post_vit=np.zeros((nit_sub,Ngrid))
x_post_v_mean=np.zeros((Ngrid))
x_post_v_95=np.zeros((Ngrid))
x_post_v_05=np.zeros((Ngrid))

# Having run MCMC then need to map regions back onto grid
sigma_y_mean=np.zeros((nmeasure))

#for yi in range(nmeasure):
#    sigma_y_mean[yi] = np.mean(sigma_y_out[yi,:])


regions_it=np.transpose(regions_out)-1
x_it=np.transpose(x_out)
#q_agg_it=np.transpose(q_agg_out)

for it in range(nit_sub):    
    for zz in range(k_it[it]):
        wh_reg = np.where(regions_it[it,:] == zz)
        x_post_vit[it,wh_reg] = x_it[it,zz+nIC]

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
h_v_all=np.zeros((nmeasuremax,Ngrid+nIC))
x_post_all = np.zeros(Ngrid+nIC)

for xi in range(nIC):
    x_post_all[xi]=np.mean(x_it[:,xi])

x_post_all[nIC:]=x_post_v_mean    
h_v_all[:,:nIC]=h_agg0[:,:nIC]
h_v_all[:,nIC:]=h_v

x_post_all_it=np.zeros((nit_sub,Ngrid+nIC))
y_post_it = np.zeros((nit_sub,nmeasuremax))
y_post=np.zeros(nmeasuremax)
y_bg_it = np.zeros((nit_sub,nmeasuremax))
y_bg=np.zeros(nmeasuremax)

x_post_all_it[:,:nIC]=x_it[:,:nIC]
x_post_all_it[:,nIC:]=x_post_vit

for it in range(nit_sub):
    y_post_it[it,:]=np.dot(h_v,x_post_all_it[it,nIC:])  
    y_bg_it[it,:]=np.dot(h_v_all[:,:nIC],x_it[it,:nIC])

y_post_it=y_post_it+y_bg_it

for ti in range(nmeasuremax):
    y_post[ti]=np.mean(y_post_it[:,ti])
    y_bg[ti]=np.mean(y_bg_it[:,ti])

#PLOTTING

data=(x_post_mean*q_ap0-q_ap0)  # mol/m2/s
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

clevels = np.arange(0., 2., 0.04)  

cs = m.contourf(mapx,mapy,x_post_mean, clevels, extend='both', cmap='RdBu_r')
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

x_uk = np.sum(x_post_v_mean[uk_index]*area_v[uk_index]*q_ap_v[uk_index]) # in mol/s
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



"""

props=["x_update", "birth", "death", "move", "sigma_y", "swap"]


post_mcmc = xray.Dataset({"x_it": (["nIt", "kICmax"],
                        x_it),
                        "k_it": (["nIt"],
                        k_it),
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

#Output files from tdmcmc_template.py stored in the form:
# "output_" + network + "_" + species +  "_" + date + ".nc"
#fname=os.path.join(output_directory,
#                    "output_" + network + "_" + species + "_" + start_date + ".nc")
#post_mcmc.to_netcdf(path=fname, mode='w')
"""
