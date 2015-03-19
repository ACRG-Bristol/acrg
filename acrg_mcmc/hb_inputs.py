# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 11:37:06 2014

@author: chxmr
"""

import numpy as np
import netCDF4 as nc
import acrg_grid as grid

def hb_inputs(ncfilename, x, x_error, y, y_model, y_error, y_site_location, H):
    
    numthreads=4   
    
    ncf=nc.Dataset(ncfilename, 'w', format="NETCDF3_CLASSIC")

    #Dimensions:
    ########################################
    statesize=len(x)
    numsites=len(set(y_site_location))
    nmeasuremax=len(y)
    nmeasuretotal=nmeasuremax/numsites

    dim1=1  #number of observations that correspond to each sigma-y
    dim2=60
    dim3=1
    dim4=1

    timeindex_nonzero=np.where(np.isfinite(np.array(y).flatten()))[0] + 1
    nmeasure=len(timeindex_nonzero)

    ncf.createDimension('statesize', statesize)
    ncf.createDimension('nmeasure', nmeasure)
    ncf.createDimension('nmeasuremax', nmeasuremax)
    ncf.createDimension('nmeasuretotal', nmeasuretotal)
    ncf.createDimension('numsites', numsites)
    ncf.createDimension('dim1', dim1)
    ncf.createDimension('dim2', dim2)
    ncf.createDimension('dim3', dim3)
    ncf.createDimension('dim4', dim4)
    ncf.createDimension('scalar', 1)


    #Variables:
    ########################################
    nit = 15000
    burn_in = 5000

    z=y[timeindex_nonzero-1]  #finite observations
    D=np.diag(y_error[timeindex_nonzero-1])  #nugget term
    
    #Define PDF types 
    #  ( 0 = LOGNORMAL, 1 = EXPONENTIAL, 2 = GAUSSIAN, 3 = UNIFORM)
    x_pdf=np.zeros(statesize) + 2
    pdf_param1_pdf=np.zeros(statesize) + 2
    pdf_param2_pdf=np.zeros(statesize) + 2
    sigma_y_pdf=np.zeros(dim2)
    sigma_ys_pdf=np.zeros(numsites)
    tau_pdf=0
    nu_pdf=3
    rho_pdf=3
    
    #Priors and hyperpriors
    x_ap=x[:]
    pdf_param1=x[:]
    pdf_param2=x_error[:]
    pdf_hyperparam1=x_error[:]
    pdf_hyperparam2=x_error[:]
    
    #Temporal component of uncertainty
    sigma_y=np.zeros(dim2) + np.mean(y_error[timeindex_nonzero-1])
    sigma_y_pdf=np.zeros(dim2)
    sigma_y_hyperparam=np.zeros(dim2) + \
        np.sqrt(np.mean(y_error[timeindex_nonzero-1]))

    #Spatial component of uncertainty
    sigma_ys=np.zeros(numsites) + np.mean(y_error[timeindex_nonzero-1])
    sigma_ys_pdf=np.zeros(numsites)
    sigma_ys_hyperparam=np.zeros(numsites) + \
        np.sqrt(np.mean(y_error[timeindex_nonzero-1]))

    tau=np.zeros(dim3) + 2.
    tau_pdf = 3
    tau_hyperparam1 = 0
    tau_hyperparam2 = tau*2.

    nu=np.zeros(dim3) + 0.5
    nu_hyperparam1=0.
    nu_hyperparam2 = nu*2

    rho=np.zeros(dim3) + 500.
    rho_hyperparam1=0.
    rho_hyperparam2 = rho*2
    
    #Stepsizes
    stepsize=x_error[:]
    stepsize_pdf_param1=x_error[:]*0.1
    stepsize_pdf_param2=x_error[:]*0.01
    stepsize_y = 0.01*y_model
    stepsize_sigma_y=np.zeros(dim2) + np.mean(y_error[timeindex_nonzero-1])*0.1
    stepsize_sigma_ys=np.zeros(numsites) + np.mean(y_error[timeindex_nonzero-1])*0.1
    stepsize_tau = 0.1*tau
    stepsize_nu = 0.1*nu
    stepsize_rho = 0.1*rho

    #Time difference array
    deltatime=np.zeros((dim2, dim2))
    for di in range(dim2):
        deltatime[:, di]=np.abs(np.arange(0, dim2, 1.)-di)

    #Distance array
    distance=np.zeros((numsites, numsites))
    locations=set(y_site_location)
    for i, loc in enumerate(locations):
        distance[i, :]=grid.haversine.distancelist(loc, locations)

    #Kron_flag
    kron_flag=0

    #T_INDICES ARRAY
    t_indices=np.reshape(range(dim1*dim2), (dim2, dim1)) + 1.
    
    #nobs
    nobs=np.zeros(numsites)
    
    #sitespresent  
    sitespresent=np.zeros(numsites)

    #distribution
    distribution=np.zeros(dim4)
    
    #datenumber
    datenumber=np.zeros(nmeasuretotal)
    
    nc_x_ap=ncf.createVariable('x_ap', 'f8', ('statesize'))
    nc_x_pdf=ncf.createVariable('x_pdf', 'f8', ('statesize'))
    nc_pdf_param1_pdf=ncf.createVariable('pdf_param1_pdf', 'f8', ('statesize'))
    nc_pdf_param2_pdf=ncf.createVariable('pdf_param2_pdf', 'f8', ('statesize'))
    nc_sigma_y_pdf=ncf.createVariable('sigma_y_pdf','f8', ('dim2',))
    nc_tau_pdf=ncf.createVariable('tau_pdf','f8', ('scalar',))
    nc_nu_pdf=ncf.createVariable('nu_pdf','f8', ('scalar',))
    nc_rho_pdf=ncf.createVariable('rho_pdf','f8', ('scalar',))
    nc_sigma_ys_pdf=ncf.createVariable('sigma_ys_pdf', 'f8', ('numsites',))
    nc_stepsize=ncf.createVariable('stepsize', 'f8', ('statesize'))
    nc_stepsize_pdf_param1=ncf.createVariable('stepsize_pdf_param1', 'f8', ('statesize'))
    nc_stepsize_pdf_param2=ncf.createVariable('stepsize_pdf_param2', 'f8', ('statesize'))
    nc_stepsize_sigma_y=ncf.createVariable('stepsize_sigma_y','f8', ('dim2',))
    nc_stepsize_tau=ncf.createVariable('stepsize_tau','f8', ('scalar',))
    nc_stepsize_nu=ncf.createVariable('stepsize_nu','f8', ('scalar',))
    nc_stepsize_rho=ncf.createVariable('stepsize_rho','f8', ('scalar',))
    nc_stepsize_y=ncf.createVariable('stepsize_y','f8', ('nmeasuremax',))
    nc_stepsize_sigma_ys=ncf.createVariable('stepsize_sigma_ys', 'f8', ('numsites',))
    nc_pdf_param1=ncf.createVariable('pdf_param1', 'f8', ('statesize',))
    nc_pdf_param2=ncf.createVariable('pdf_param2', 'f8', ('statesize',))
    nc_sigma_y=ncf.createVariable('sigma_y','f8', ('dim2',))
    nc_sigma_ys=ncf.createVariable('sigma_ys', 'f8', ('numsites',))
    nc_D=ncf.createVariable('D', 'f8', ('nmeasure', 'nmeasure'))
    nc_pdf_hyperparam1=ncf.createVariable('pdf_hyperparam1', 'f8', ('statesize',))
    nc_pdf_hyperparam2=ncf.createVariable('pdf_hyperparam2', 'f8', ('statesize',))
    nc_sigma_y_hyperparam=ncf.createVariable('sigma_y_hyperparam','f8', ('dim2',))
    nc_tau_hyperparam1=ncf.createVariable('tau_hyperparam1','f8', ('scalar',))
    nc_tau_hyperparam2=ncf.createVariable('tau_hyperparam2','f8', ('scalar',))
    nc_sigma_ys_hyperparam=ncf.createVariable('sigma_ys_hyperparam', 'f8', ('numsites',))
    nc_nu_hyperparam1=ncf.createVariable('nu_hyperparam1','f8', ('scalar',))
    nc_rho_hyperparam1=ncf.createVariable('rho_hyperparam1','f8', ('scalar',))
    nc_nu_hyperparam2=ncf.createVariable('nu_hyperparam2','f8', ('scalar',))
    nc_rho_hyperparam2=ncf.createVariable('rho_hyperparam2','f8', ('scalar',))
    nc_H=ncf.createVariable('H', 'f8', ('statesize', 'nmeasuremax'))
    nc_nit=ncf.createVariable('nit','f8', ('scalar',))
    nc_burn_in=ncf.createVariable('burn_in','f8', ('scalar',))
    nc_nmeasure=ncf.createVariable('nmeasure','f8', ('scalar',))
    nc_statesize=ncf.createVariable('statesize','f8', ('scalar',))
    nc_z=ncf.createVariable('z', 'f8', ('nmeasure',))
    nc_nobs=ncf.createVariable('nobs', 'f8', ('numsites',))
    nc_T_indices=ncf.createVariable('T_indices', 'f8', ('dim2', 'dim1'))
    nc_dim1=ncf.createVariable('dim1','f8', ('scalar',))
    nc_dim2=ncf.createVariable('dim2','f8', ('scalar',))
    nc_dim3=ncf.createVariable('dim3','f8', ('scalar',))
    nc_dim4=ncf.createVariable('dim4','f8', ('scalar',))
    nc_deltatime=ncf.createVariable('deltatime','f8', ('nmeasuretotal','nmeasuretotal'))
    nc_distance=ncf.createVariable('distance', 'f8', ('numsites','numsites'))
    nc_tau=ncf.createVariable('tau','f8', ('scalar',))
    nc_nu=ncf.createVariable('nu','f8', ('scalar',))
    nc_rho=ncf.createVariable('rho','f8', ('scalar',))
    nc_nmeasuretotal=ncf.createVariable('nmeasuretotal','f8', ('scalar',))
    nc_datenumber=ncf.createVariable('datenumber', 'f8', ('nmeasuretotal',))
    nc_timeindex_nonzero=ncf.createVariable('timeindex_nonzero', 'f8', ('nmeasure',))
    nc_numsites=ncf.createVariable('numsites','f8', ('scalar',))
    nc_sitespresent=ncf.createVariable('sitespresent', 'f8', ('numsites',))
    nc_kron_flag=ncf.createVariable('kron_flag','f8', ('scalar',))
    nc_numthreads=ncf.createVariable('numthreads','f8', ('scalar',))
    nc_distribution=ncf.createVariable('distribution','f8', ('dim4',))
    
    nc_x_ap[:]=x_ap
    nc_x_pdf[:]=x_pdf
    nc_pdf_param1_pdf[:]=pdf_param1_pdf
    nc_pdf_param2_pdf[:]=pdf_param2_pdf
    nc_sigma_y_pdf[:]=sigma_y_pdf
    nc_tau_pdf[:]=tau_pdf
    nc_nu_pdf[:]=nu_pdf
    nc_rho_pdf[:]=rho_pdf
    nc_sigma_ys_pdf[:]=sigma_ys_pdf
    nc_stepsize[:]=stepsize
    nc_stepsize_pdf_param1[:]=stepsize_pdf_param1
    nc_stepsize_pdf_param2[:]=stepsize_pdf_param2
    nc_stepsize_tau[:]=stepsize_tau
    nc_stepsize_nu[:]=stepsize_nu
    nc_stepsize_rho[:]=stepsize_rho
    nc_stepsize_sigma_y[:]=stepsize_sigma_y
    nc_stepsize_y[:]=stepsize_y
    nc_stepsize_sigma_ys[:]=stepsize_sigma_ys
    nc_pdf_param1[:]=pdf_param1
    nc_pdf_param2[:]=pdf_param2
    nc_sigma_y[:]=sigma_y
    nc_sigma_ys[:]=sigma_ys
    nc_D[:]=D
    nc_pdf_hyperparam1[:]=pdf_hyperparam1
    nc_pdf_hyperparam2[:]=pdf_hyperparam2
    nc_sigma_y_hyperparam[:]=sigma_y_hyperparam
    nc_tau_hyperparam1[:]=tau_hyperparam1
    nc_tau_hyperparam2[:]=tau_hyperparam2
    nc_sigma_ys_hyperparam[:]=sigma_ys_hyperparam
    nc_nu_hyperparam1[:]=nu_hyperparam1
    nc_rho_hyperparam1[:]=rho_hyperparam1
    nc_nu_hyperparam2[:]=nu_hyperparam2
    nc_rho_hyperparam2[:]=rho_hyperparam2
    nc_H[:,:]=H
    nc_nit[:]=nit
    nc_burn_in[:]=burn_in
    nc_nmeasure[:]=nmeasure
    nc_statesize[:]=statesize
    nc_z[:]=z
    nc_nobs[:]=nobs
    nc_T_indices[:,:]=t_indices
    nc_dim1[:]=dim1
    nc_dim2[:]=dim2
    nc_dim3[:]=dim3
    nc_dim4[:]=dim4
    nc_deltatime[:]=deltatime
    nc_distance[:,:]=distance
    nc_tau[:]=tau
    nc_nu[:]=nu
    nc_rho[:]=rho
    nc_nmeasuretotal[:]=nmeasuretotal
    nc_datenumber[:]=datenumber
    nc_timeindex_nonzero[:]=timeindex_nonzero
    nc_numsites[:]=numsites
    nc_sitespresent[:]=sitespresent
    nc_kron_flag[:]=kron_flag
    nc_numthreads[:]=numthreads
    nc_distribution[:]=distribution


    ncf.close()
    
    
    