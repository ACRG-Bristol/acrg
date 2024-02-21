#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:57:23 2020

@author: lw13938
"""
import numpy as np
import pymc3 as pm
import pandas as pd
import xarray as xr
import getpass
from pathlib import Path

import acrg.name.name as name
from acrg.grid.areagrid import areagrid
from acrg.hbmcmc.inversionsetup import opends
from .inversionsetup import offset_matrix
from acrg.hbmcmc.hbmcmc_output import define_output_filename
import acrg.convert as convert
from acrg.config.version import code_version
from acrg.config.paths import Paths

from scipy import stats

import sys

data_path = Paths.data

def parsePrior(name, prior_params, shape = ()):
    """
    Parses all continuous distributions for pymc 3.8: https://docs.pymc.io/api/distributions/continuous.html
    This format requires updating when the pymc distributions update, but is safest for code execution
    
    Args:
        name (str):
            name of variable in the pymc model
        prior_params (dict):
            dict of parameters for the distribution, including 'pdf' for the distribution to use
        shape (array):
            shape of distribution to be created. Default shape = () is the same as used by pymc3
        
    """
    functionDict = {"uniform":pm.Uniform,
                    "flat":pm.Flat,
                    "halfflat":pm.HalfFlat,
                    "normal":pm.Normal,
                    "truncatednormal":pm.TruncatedNormal,
                    "halfnormal":pm.HalfNormal,
                    "skewnormal":pm.SkewNormal,
                    "beta":pm.Beta,
                    "kumaraswamy":pm.Kumaraswamy,
                    "exponential":pm.Exponential,
                    "laplace":pm.Laplace,
                    "studentt":pm.StudentT,
                    "halfstudentt":pm.HalfStudentT,
                    "cauchy":pm.Cauchy,
                    "halfcauchy":pm.HalfCauchy,
                    "gamma":pm.Gamma,
                    "inversegamma":pm.InverseGamma,
                    "weibull":pm.Weibull,
                    "lognormal":pm.Lognormal,
                    "chisquared":pm.ChiSquared,
                    "wald":pm.Wald,
                    "pareto":pm.Pareto,
                    "exgaussian":pm.ExGaussian,
                    "vonmises":pm.VonMises,
                    "triangular":pm.Triangular,
                    "gumbel":pm.Gumbel,
                    "rice":pm.Rice,
                    "logistic":pm.Logistic,
                    "logitnormal":pm.LogitNormal,
                    "interpolated":pm.Interpolated}
    
    pdf = prior_params["pdf"]
    #Get a dictionary of the pdf arguments
    params = {x: prior_params[x] for x in prior_params if x != "pdf"}
    return functionDict[pdf.lower()](name, shape=shape, **params)

def inferpymc3(Hx, Hbc, Y, error, siteindicator, sigma_freq_index,
               xprior={"pdf":"lognormal", "mu":1, "sd":1},
               bcprior={"pdf":"lognormal", "mu":0.004, "sd":0.02},
               sigprior={"pdf":"uniform", "lower":0.5, "upper":3},
               nit=2.5e5, burn=50000, tune=1.25e5, nchain=2, 
               sigma_per_site = True, 
               offsetprior={"pdf":"normal", "mu":0, "sd":0.6},
               add_offset = False, verbose=False):       
    """
    Uses pym3 module for Bayesian inference for emissions field, boundary 
    conditions and (currently) a single model error value.
    This uses a Normal likelihood but the (hyper)prior PDFs can selected by user.
    
    Args:
        Hx (array):
            Transpose of the sensitivity matrix to map emissions to measurement.
            This is the same as what is given from fp_data[site].H.values, where
            fp_data is the output from e.g. footprint_data_merge, but where it
            has been stacked for all sites.
        Hbc (array):
            Same as above but for boundary conditions
        Y (array):
            Measurement vector containing all measurements
        error (arrray):
            Measurement error vector, containg a value for each element of Y.
        siteindicator (array):
            Array of indexing integers that relate each measurement to a site
        sigma_freq_index (array):
            Array of integer indexes that converts time into periods
        xprior (dict):
            Dictionary containing information about the prior PDF for emissions.
            The entry "pdf" is the name of the analytical PDF used, see
            https://docs.pymc.io/api/distributions/continuous.html for PDFs
            built into pymc3, although they may have to be coded into the script.
            The other entries in the dictionary should correspond to the shape
            parameters describing that PDF as the online documentation,
            e.g. N(1,1**2) would be: xprior={pdf:"normal", "mu":1, "sd":1}.
            Note that the standard deviation should be used rather than the 
            precision. Currently all variables are considered iid.
        bcprior (dict):
            Same as above but for boundary conditions.
        sigprior (dict):
            Same as above but for model error.
        offsetprior (dict):
            Same as above but for bias offset. Only used is addoffset=True.
        sigma_per_site (bool):
            Whether a model sigma value will be calculated for each site independantly (True) or all sites together (False).
            Default: True
        add_offset (bool):
            Add an offset (intercept) to all sites but the first in the site list. Default False.
        verbose:
            When True, prints progress bar

            
    Returns:
        outs (array):
            MCMC chain for emissions scaling factors for each basis function.
        bcouts (array):
            MCMC chain for boundary condition scaling factors.
        sigouts (array):
            MCMC chain for model error.
        Ytrace (array):
            MCMC chain for modelled obs.
        YBCtrace (array):
            MCMC chain for modelled boundary condition.
        convergence (str):
            Passed/Failed convergence test as to whether mutliple chains
            have a Gelman-Rubin diagnostic value <1.05
        step1 (str):
            Type of MCMC sampler for emissions and boundary condition updates.
            Currently it's hardwired to NUTS (probably wouldn't change this
            unless you're doing something obscure).
        step2 (str):
            Type of MCMC sampler for model error updates.
            Currently it's hardwired to a slice sampler. This parameter is low
            diensional and quite simple with a slice sampler, although could 
            easily be changed.
     
    TO DO:
       - Allow non-iid variables
    """
    burn = int(burn)         
    
    hx = Hx.T 
    hbc = Hbc.T
    nx = hx.shape[1]
    nbc = hbc.shape[1]
    ny = len(Y)
    
    nit = int(nit)  
    
    #convert siteindicator into a site indexer
    if sigma_per_site:
        sites = siteindicator.astype(int)
        nsites = np.amax(sites)+1
    else:
        sites = np.zeros_like(siteindicator).astype(int)
        nsites = 1
    nsigmas = np.amax(sigma_freq_index)+1
    
    if add_offset:
        B = offset_matrix(siteindicator)

    with pm.Model() as model:
        x = parsePrior("x", xprior, shape=nx)
        xbc = parsePrior("xbc", bcprior, shape=nbc)
        sig = parsePrior("sig", sigprior, shape=(nsites, nsigmas))
        if add_offset:
            offset = parsePrior("offset", offsetprior, shape=nsites-1) 
            offset_vec = pm.math.concatenate( (np.array([0]), offset), axis=0)
            mu = pm.math.dot(hx,x) + pm.math.dot(hbc,xbc) + pm.math.dot(B, offset_vec)
        else:
            mu = pm.math.dot(hx,x) + pm.math.dot(hbc,xbc) 

        model_error = pm.math.maximum(np.abs(pm.math.dot(hx,x)) * sig[sites, sigma_freq_index], 2.0)       

        epsilon = pm.math.sqrt(error**2 + model_error**2)
        y = pm.Normal('y', mu=mu, sd=epsilon, observed=Y, shape=ny)
        
        step1 = pm.NUTS(vars=[x,xbc])
#        step1 = pm.Metropolis(vars=[x,xbc])
        step2 = pm.Slice(vars=[sig])
        
        trace = pm.sample(nit, tune=int(tune), chains=nchain,
                           step=[step1,step2], progressbar=verbose, cores=nchain)#step=pm.Metropolis())#  #target_accept=0.8,
        
        outs = trace.get_values(x, burn=burn)[0:int((nit)-burn)]
        bcouts = trace.get_values(xbc, burn=burn)[0:int((nit)-burn)]
        sigouts = trace.get_values(sig, burn=burn)[0:int((nit)-burn)]
        
        #Check for convergence
        gelrub = pm.rhat(trace)['x'].max()
        if gelrub > 1.05:
            print('Failed Gelman-Rubin at 1.05')
            convergence = "Failed"
        else:
            convergence = "Passed"
        
        if add_offset:
            offset_outs = trace.get_values(offset, burn=burn)[0:int((nit)-burn)]
            offset_trace = np.hstack([np.zeros((int(nit-burn),1)), offset_outs])
            YBCtrace = np.dot(Hbc.T,bcouts.T) + np.dot(B, offset_trace.T) 
            OFFtrace = np.dot(B, offset_trace.T)             
        else:
            YBCtrace = np.dot(Hbc.T,bcouts.T)
            offset_outs = outs * 0
            OFFtrace = YBCtrace * 0
        Ytrace = np.dot(Hx.T,outs.T) + YBCtrace

        # Added option to return offset_outs          

 
        return outs, bcouts, sigouts, offset_outs, Ytrace, YBCtrace, OFFtrace, convergence, step1, step2

def inferpymc3_postprocessouts(xouts, bcouts, sigouts, offset_outs, convergence, 
                               Hx, Hbc, Y, error, Ytrace, YBCtrace, offset_trace, 
                               step1, step2, 
                               xprior, bcprior, sigprior, offsetprior, Ytime, siteindicator, sigma_freq_index,
                               domain, species, sites,
                               start_date, end_date, outputname, outputpath,
                               country_unit_prefix,
                               burn, tune, nchain, sigma_per_site,
                               fp_data=None, flux_directory=None, emissions_name=None, 
                               basis_directory=None, country_file=None,
                               add_offset=False, rerun_file=None):

        """
        Takes the output from inferpymc3 function, along with some other input
        information, and places it all in a netcdf output. This function also 
        calculates the mean posterior emissions for the countries in the 
        inversion domain and saves it to netcdf.
        Note that the uncertainties are defined by the highest posterior 
        density (HPD) region and not percentiles (as the tdMCMC code). 
        The HPD region is defined, for probability content (1-a), as:
            1) P(x \in R | y) = (1-a)
            2) for x1 \in R and x2 \notin R, P(x1|y)>=P(x2|y)
        
        Args:
            xouts (array):
                MCMC chain for emissions scaling factors for each basis function.
            bcouts (array):
                MCMC chain for boundary condition scaling factors.
            sigouts (array):
                MCMC chain for model error.
            convergence (str):
                Passed/Failed convergence test as to whether mutliple chains
                have a Gelman-Rubin diagnostic value <1.05
            Hx (array):
                Transpose of the sensitivity matrix to map emissions to measurement.
                This is the same as what is given from fp_data[site].H.values, where
                fp_data is the output from e.g. footprint_data_merge, but where it
                has been stacked for all sites.
            Hbc (array):
                Same as above but for boundary conditions
            Y (array):
                Measurement vector containing all measurements
            error (arrray):
                Measurement error vector, containg a value for each element of Y.
            Ytrace (array):
                Trace of modelled y values calculated from mcmc outputs and H matrices
            YBCtrace (array):
                Trace of modelled boundary condition values calculated from mcmc outputs and Hbc matrices
            step1 (str):
                Type of MCMC sampler for emissions and boundary condition updates.
            step2 (str):
                Type of MCMC sampler for model error updates.
            xprior (dict):
                Dictionary containing information about the prior PDF for emissions.
                The entry "pdf" is the name of the analytical PDF used, see
                https://docs.pymc.io/api/distributions/continuous.html for PDFs
                built into pymc3, although they may have to be coded into the script.
                The other entries in the dictionary should correspond to the shape
                parameters describing that PDF as the online documentation,
                e.g. N(1,1**2) would be: xprior={pdf:"normal", "mu":1, "sd":1}.
                Note that the standard deviation should be used rather than the 
                precision. Currently all variables are considered iid.
            bcprior (dict):
                Same as above but for boundary conditions.
            sigprior (dict):
                Same as above but for model error.
            offsetprior (dict):
                Same as above but for bias offset. Only used is addoffset=True.
            Ytime (pandas datetime array):
                Time stamp of measurements as used by the inversion.
            siteindicator (array):
                Numerical indicator of which site the measurements belong to,
                same length at Y.
            sigma_freq_index (array):
                Array of integer indexes that converts time into periods
            domain (str):
                Inversion spatial domain.
            species (str):
                Species of interest
            sites (list):
                List of sites in inversion
            start_date (str):
                Start time of inversion "YYYY-mm-dd"
            end_date (str):
                End time of inversion "YYYY-mm-dd"
            outputname (str):
                Unique identifier for output/run name.
            outputpath (str):
                Path to where output should be saved.
            country_unit_prefix ('str', optional)
                A prefix for scaling the country emissions. Current options are: 
                'T' will scale to Tg, 'G' to Gg, 'M' to Mg, 'P' to Pg.
                To add additional options add to acrg_convert.prefix
                Default is none and no scaling will be applied (output in g).
            burn (int):
                Number of iterations burned in MCMC
            tune (int):
                Number of iterations used to tune step size
            nchain (int):
                Number of independent chains run 
            sigma_per_site (bool):
                Whether a model sigma value was be calculated for each site independantly (True) 
                or all sites together (False).
            fp_data (dict, optional):
                Output from footprints_data_merge + sensitivies
            flux_directory (str, optional):
                Directory containing the emissions data if not default
            emissions_name (dict, optional): 
                Allows emissions files with filenames that are longer than just the species name
                to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
                with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
                multiple sources can be read in simultaneously if they are added as separate entries to
                the emissions_name dictionary.
                If using HiTRes footprints, both the high and low frequency emissions files must be specified
                in a second dictionary like so: {'anth': {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}}.
                It is not a problem to have a mixture of sources, with some that use HiTRes footprints and some
                that don't.
            basis_directory (str, optional):
                Directory containing basis function file
            country_file (str, optional):
                Path of country definition file
            add_offset (bool):
                Add an offset (intercept) to all sites but the first in the site list. Default False.
            rerun_file (xarray dataset, optional):
                An xarray dataset containing the ncdf output from a previous run of the MCMC code.
                
        Returns:
            netdf file containing results from inversion
            
        TO DO:
            - Look at compressability options for netcdf output
            - I'm sure the number of inputs can be cut down or found elsewhere.
            - Currently it can only work out the country total emissions if
              the a priori emissions are constant over the inversion period
              or else monthly (and inversion is for less than one calendar year).
        """
        print("Post-processing output")
        
        #Get parameters for output file 
        nit = xouts.shape[0]
        nx = Hx.shape[0]
        ny = len(Y)
        nbc = Hbc.shape[0]
        noff = offset_outs.shape[0]

        nui = np.arange(2)
        steps = np.arange(nit)
        nmeasure = np.arange(ny)
        nparam = np.arange(nx)
        nBC = np.arange(nbc)
        nOFF = np.arange(noff)
        #YBCtrace = np.dot(Hbc.T,bcouts.T)

        # OFFSET HYPERPARAMETER
        YmodmuOFF = np.mean(offset_trace,axis=1)    # mean
        YmodmedOFF = np.median(offset_trace,axis=1) # median
        YmodmodeOFF = np.zeros(shape=offset_trace.shape[0]) # mode

        for i in range(0, offset_trace.shape[0]):
            # If sufficient no. of iterations in MCMC chains use a KDE to calculate
            # the mode. Else, mean value is used.
            if np.nanmax(offset_trace[i,:]) > np.nanmin(offset_trace[i,:]):
                xes_off = np.linspace(np.nanmin(offset_trace[i,:]), np.nanmax(offset_trace[i,:]), 200)
                kde = stats.gaussian_kde(offset_trace[i,:]).evaluate(xes_off)
                YmodmodeOFF[i] = xes_off[kde.argmax()]
            else:
                YmodmodeOFF[i] = np.mean(offset_trace[i,:])  

        Ymod68OFF = pm.stats.hdi(offset_trace.T, 0.68)
        Ymod95OFF = pm.stats.hdi(offset_trace.T, 0.95)


        # BC
        YmodmuBC = np.mean(YBCtrace, axis=1)
        YmodmedBC = np.median(YBCtrace, axis=1)
        YmodmodeBC = np.zeros(shape=YBCtrace.shape[0])

        for i in range(0, YBCtrace.shape[0]):
            # If sufficient no. of iterations in MCMC chains use a KDE to calculate
            # the mode. Else, mean value is used.
            if np.nanmax(YBCtrace[i,:]) > np.nanmin(YBCtrace[i,:]):
                xes_bc = np.linspace(np.nanmin(YBCtrace[i,:]), np.nanmax(YBCtrace[i,:]), 200)
                kde = stats.gaussian_kde(YBCtrace[i,:]).evaluate(xes_bc)
                YmodmodeBC[i] = xes_bc[kde.argmax()]
            else:
                YmodmodeBC[i] = np.mean(YBCtrace[i,:])

        Ymod95BC = pm.stats.hdi(YBCtrace.T, 0.95)
        Ymod68BC = pm.stats.hdi(YBCtrace.T, 0.68)
        YaprioriBC = np.sum(Hbc, axis=0)


        # Y-VALUES HYPERPARAMETER (XOUTS * H) 
        Ymodmu = np.mean(Ytrace, axis=1)
        Ymodmed = np.median(Ytrace, axis=1)
        Ymodmode = np.zeros(shape=Ytrace.shape[0])

        for i in range(0, Ytrace.shape[0]):
            # If sufficient no. of iterations in MCMC chains use a KDE to calculate
            # the mode. Else, mean value is used.
            if np.nanmax(Ytrace[i,:]) > np.nanmin(Ytrace[i,:]):
                xes = np.arange(np.nanmin(Ytrace[i,:]), np.nanmax(Ytrace[i,:]), 0.01)
                kde = stats.gaussian_kde(Ytrace[i,:]).evaluate(xes)
                Ymodmode[i] = xes[kde.argmax()]
            else:
                Ymodmode[i] = np.mean(Ytrace[i,:])

        Ymod95 = pm.stats.hdi(Ytrace.T, 0.95)
        Ymod68 = pm.stats.hdi(Ytrace.T, 0.68)
        Yapriori = np.sum(Hx.T, axis=1) + np.sum(Hbc.T, axis=1)
        sitenum = np.arange(len(sites))
 
        if fp_data is None and rerun_file is not None:
            lon = rerun_file.lon.values
            lat = rerun_file.lat.values
            site_lat = rerun_file.sitelats.values
            site_lon = rerun_file.sitelons.values
            bfds = rerun_file.basisfunctions
        else:
            lon = fp_data[sites[0]].lon.values
            lat = fp_data[sites[0]].lat.values
            site_lat = np.zeros(len(sites))
            site_lon = np.zeros(len(sites))
            for si, site in enumerate(sites):
                site_lat[si] = fp_data[site].release_lat.values[0]
                site_lon[si] = fp_data[site].release_lon.values[0]
            bfds = fp_data[".basis"]


        #Calculate mean and mode posterior scale map and flux field
        scalemap_mu = np.zeros_like(bfds.values)
        scalemap_mode = np.zeros_like(bfds.values)        

        for npm in nparam:
            scalemap_mu[bfds.values == (npm+1)] = np.mean(xouts[:,npm])
            if np.nanmax(xouts[:,npm]) > np.nanmin(xouts[:,npm]):
                xes = np.arange(np.nanmin(xouts[:,npm]), np.nanmax(xouts[:,npm]), 0.01)
                kde = stats.gaussian_kde(xouts[:,npm]).evaluate(xes)           
                scalemap_mode[bfds.values == (npm+1)] = xes[kde.argmax()]
            else:
                scalemap_mode[bfds.values == (npm+1)] = np.nanmean(xouts[:,npm])
    

        if rerun_file is not None:
            # Note, at the moment fluxapriori in the output is the mean apriori flux over the 
            # inversion period and so will not be identical to the original a priori flux, if 
            # it varies over the inversion period
            emissions_flux = np.expand_dims(rerun_file.fluxapriori.values,2)
        else:
            if emissions_name == None:
                emds = name.flux(domain, species, start = start_date, end = end_date, flux_directory=flux_directory)
            else:
                emds = name.flux(domain, list(emissions_name.values())[0], start = start_date, end = end_date, flux_directory=flux_directory)
            emissions_flux = emds.flux.values

        # NB. using mode scalemap for posterior fluxes
        flux = scalemap_mode*emissions_flux[:,:,0]
        
        #Basis functions to save
        bfarray = bfds.values-1
    
        #Calculate country totals   
        area = areagrid(lat, lon)
        if not rerun_file:
            c_object = name.get_country(domain, country_file=country_file)
            cntryds = xr.Dataset({'country': (['lat','lon'], c_object.country), 
                            'name' : (['ncountries'],c_object.name) },
                                            coords = {'lat': (c_object.lat),
                                            'lon': (c_object.lon)})
            cntrynames = cntryds.name.values
            cntrygrid = cntryds.country.values
        else:
            cntrynames = rerun_file.countrynames.values
            cntrygrid = rerun_file.countrydefinition.values

        cntrymean = np.zeros((len(cntrynames)))
        cntrymedian = np.zeros((len(cntrynames)))
        cntrymode = np.zeros((len(cntrynames)))
        cntry68 = np.zeros((len(cntrynames), len(nui)))
        cntry95 = np.zeros((len(cntrynames), len(nui)))
        cntrysd = np.zeros(len(cntrynames))
        cntryprior = np.zeros(len(cntrynames))
        molarmass = convert.molar_mass(species)

        unit_factor = convert.prefix(country_unit_prefix)
        if country_unit_prefix is None:
            country_unit_prefix=''
        country_units = country_unit_prefix + 'g'
        if rerun_file is not None:
            obs_units = rerun_file.Yobs.attrs["units"].split(" ")[0]
        else:
            obs_units = str(fp_data[".units"])

        # Not sure how it's best to do this if multiple months in emissions 
        # file. Now it scales a weighted average of a priori emissions
        # If a priori emissions have frequency of more than monthly then this
        # needs chaning.
        aprioriflux = np.zeros_like(area)
        if emissions_flux.shape[2] > 1:
            print("Assuming the inversion is over a year or less and emissions file is monthly")
            allmonths = pd.date_range(start_date, end_date).month[:-1].values
            allmonths -= np.min(allmonths)
            for mi in allmonths:
                aprioriflux += emissions_flux[:,:,mi]*np.sum(allmonths == mi)/len(allmonths)
        else:
            aprioriflux = np.squeeze(emissions_flux)

        for ci, cntry in enumerate(cntrynames):
            cntrytottrace = np.zeros(len(steps))
            cntrytotprior = 0
            for bf in range(int(np.max(bfarray))+1):
                bothinds = np.logical_and(cntrygrid == ci, bfarray==bf)
                cntrytottrace += np.sum(area[bothinds].ravel()*aprioriflux[bothinds].ravel()* \
                               3600*24*365*molarmass)*xouts[:,bf]/unit_factor
                cntrytotprior += np.sum(area[bothinds].ravel()*aprioriflux[bothinds].ravel()* \
                              3600*24*365*molarmass)/unit_factor

            cntrymean[ci] = np.mean(cntrytottrace)
            cntrymedian[ci] = np.median(cntrytottrace)

            if np.nanmax(cntrytottrace) > np.nanmin(cntrytottrace):
                xes = np.linspace(np.nanmin(cntrytottrace), np.nanmax(cntrytottrace), 200)
                kde = stats.gaussian_kde(cntrytottrace).evaluate(xes)
                cntrymode[ci] = xes[kde.argmax()]
            else:
                cntrymode[ci] = np.mean(cntrytottrace)

            cntrysd[ci] = np.std(cntrytottrace)
            cntry68[ci, :] = pm.stats.hdi(cntrytottrace, 0.68)
            cntry95[ci, :] = pm.stats.hdi(cntrytottrace, 0.95)
            cntryprior[ci] = cntrytotprior
            
    
        #Make output netcdf file
        outds = xr.Dataset({'Yobs':(['nmeasure'], Y),
                            'Yerror' :(['nmeasure'], error),                          
                            'Ytime':(['nmeasure'],Ytime),
                            'Yapriori':(['nmeasure'],Yapriori),
                            'Ymodmean':(['nmeasure'], Ymodmu), 
                            'Ymodmedian':(['nmeasure'], Ymodmed),
                            'Ymodmode': (['nmeasure'], Ymodmode),
                            'Ymod95':(['nmeasure','nUI'], Ymod95),
                            'Ymod68':(['nmeasure','nUI'], Ymod68),
                            'Yoffmean':(['nmeasure'], YmodmuOFF),
                            'Yoffmedian':(['nmeasure'], YmodmedOFF),
                            'Yoffmode':(['nmeasure'], YmodmodeOFF),
                            'Yoff68':(['nmeasure','nUI'], Ymod68OFF),
                            'Yoff95':(['nmeasure','nUI'], Ymod95OFF),
                            'YaprioriBC':(['nmeasure'],YaprioriBC),                            
                            'YmodmeanBC':(['nmeasure'], YmodmuBC),
                            'YmodmedianBC':(['nmeasure'], YmodmedBC),
                            'YmodmodeBC':(['nmeasure'], YmodmodeBC),
                            'Ymod95BC':(['nmeasure','nUI'], Ymod95BC),
                            'Ymod68BC':(['nmeasure','nUI'], Ymod68BC),                        
                            'xtrace':(['steps','nparam'], xouts),
                            'bctrace':(['steps','nBC'],bcouts),
                            'sigtrace':(['steps', 'nsigma_site', 'nsigma_time'], sigouts),
                            'siteindicator':(['nmeasure'],siteindicator),
                            'sigmafreqindex':(['nmeasure'],sigma_freq_index),
                            'sitenames':(['nsite'],sites),
                            'sitelons':(['nsite'],site_lon),
                            'sitelats':(['nsite'],site_lat),
                            'fluxapriori':(['lat','lon'], aprioriflux), #NOTE this is the mean a priori flux over the inversion period
                            'fluxmode':(['lat','lon'], flux),                            
                            'scalingmean':(['lat','lon'],scalemap_mu),
                            'scalingmode':(['lat','lon'],scalemap_mode),
                            'basisfunctions':(['lat','lon'],bfarray),
                            'countrymean':(['countrynames'], cntrymean),
                            'countrymedian':(['countrynames'], cntrymedian),
                            'countrymode':(['countrynames'], cntrymode),
                            'countrysd':(['countrynames'], cntrysd),
                            'country68':(['countrynames', 'nUI'],cntry68),
                            'country95':(['countrynames', 'nUI'],cntry95),
                            'countryapriori':(['countrynames'],cntryprior),
                            'countrydefinition':(['lat','lon'], cntrygrid),
                            'xsensitivity':(['nmeasure','nparam'], Hx.T),
                            'bcsensitivity':(['nmeasure', 'nBC'],Hbc.T)},
                        coords={'stepnum' : (['steps'], steps), 
                                   'paramnum' : (['nlatent'], nparam),
                                   'numBC' : (['nBC'], nBC),

                                 #  'numOFF': (['nOFF'],nOFF),

                                   'measurenum' : (['nmeasure'], nmeasure), 
                                   'UInum' : (['nUI'], nui),
                                   'nsites': (['nsite'], sitenum),
                                   'nsigma_time': (['nsigma_time'], np.unique(sigma_freq_index)),
                                   'nsigma_site': (['nsigma_site'], np.arange(sigouts.shape[1]).astype(int)),
                                   'lat':(['lat'],lat),
                                   'lon':(['lon'],lon),
                                   'countrynames':(['countrynames'],cntrynames)})
        
        outds.fluxmode.attrs["units"] = "mol/m2/s"
        outds.fluxapriori.attrs["units"] = "mol/m2/s"
        outds.Yobs.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yapriori.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymodmean.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymodmedian.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymodmode.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymod95.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymod68.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yoffmean.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yoffmedian.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yoffmode.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yoff95.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yoff68.attrs["units"] = obs_units+" "+"mol/mol"
        outds.YmodmeanBC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.YmodmedianBC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.YmodmodeBC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymod95BC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Ymod68BC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.YaprioriBC.attrs["units"] = obs_units+" "+"mol/mol"
        outds.Yerror.attrs["units"] = obs_units+" "+"mol/mol"
        outds.countrymean.attrs["units"] = country_units
        outds.countrymedian.attrs["units"] = country_units
        outds.countrymode.attrs["units"] = country_units
        outds.country68.attrs["units"] = country_units
        outds.country95.attrs["units"] = country_units
        outds.countrysd.attrs["units"] = country_units
        outds.countryapriori.attrs["units"] = country_units
        outds.xsensitivity.attrs["units"] = obs_units+" "+"mol/mol"
        outds.bcsensitivity.attrs["units"] = obs_units+" "+"mol/mol"
        outds.sigtrace.attrs["units"] = obs_units+" "+"mol/mol"
        
        outds.Yobs.attrs["longname"] = "observations"
        outds.Yerror.attrs["longname"] = "measurement error"
        outds.Ytime.attrs["longname"] = "time of measurements"
        outds.Yapriori.attrs["longname"] = "a priori simulated measurements"
        outds.Ymodmean.attrs["longname"] = "mean of posterior simulated measurements"
        outds.Ymodmedian.attrs["longname"] = "median of posterior simulated measurements"
        outds.Ymodmode.attrs["longname"] = "mode of posterior simulated measurements"
        outds.Ymod68.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated measurements"
        outds.Ymod95.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated measurements"
        outds.Yoffmean.attrs["longname"] = "mean of posterior simulated offset between measurements"
        outds.Yoffmedian.attrs["longname"] = "median of posterior simulated offset between measurements"
        outds.Yoffmode.attrs["longname"] = "mode of posterior simulated offset between measurements"
        outds.Yoff68.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated offset between measurements"
        outds.Yoff95.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated offset between measurements"
        outds.YaprioriBC.attrs["longname"] = "a priori simulated boundary conditions"
        outds.YmodmeanBC.attrs["longname"] = "mean of posterior simulated boundary conditions"
        outds.YmodmedianBC.attrs["longname"] = "median of posterior simulated boundary conditions"
        outds.YmodmodeBC.attrs["longname"] = "mode of posterior simulated boundary conditions"
        outds.Ymod68BC.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated boundary conditions"
        outds.Ymod95BC.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated boundary conditions"
        outds.xtrace.attrs["longname"] = "trace of unitless scaling factors for emissions parameters"
        outds.bctrace.attrs["longname"] = "trace of unitless scaling factors for boundary condition parameters"
        outds.sigtrace.attrs["longname"] = "trace of model error parameters"
        outds.siteindicator.attrs["longname"] = "index of site of measurement corresponding to sitenames"
        outds.sigmafreqindex.attrs["longname"] = "perdiod over which the model error is estimated"
        outds.sitenames.attrs["longname"] = "site names"
        outds.sitelons.attrs["longname"] = "site longitudes corresponding to site names"
        outds.sitelats.attrs["longname"] = "site latitudes corresponding to site names"
        outds.fluxapriori.attrs["longname"] = "mean a priori flux over period"
        outds.fluxmode.attrs["longname"] = "mean posterior flux over period"
        outds.scalingmean.attrs["longname"] = "mean scaling factor field over period"
        outds.scalingmode.attrs["longname"] = "mode scaling factor field over period"
        outds.basisfunctions.attrs["longname"] = "basis function field"
        outds.countrymean.attrs["longname"] = "mean of ocean and country totals"
        outds.countrymedian.attrs["longname"] = "median of ocean and country totals"
        outds.countrymode.attrs["longname"] = "mode of ocean and country totals"
        outds.country68.attrs["longname"] = "0.68 Bayesian credible interval of ocean and country totals"
        outds.country95.attrs["longname"] = "0.95 Bayesian credible interval of ocean and country totals"        
        outds.countrysd.attrs["longname"] = "standard deviation of ocean and country totals" 
        outds.countryapriori.attrs["longname"] = "prior mean of ocean and country totals"
        outds.countrydefinition.attrs["longname"] = "grid definition of countries"
        outds.xsensitivity.attrs["longname"] = "emissions sensitivity timeseries"   
        outds.bcsensitivity.attrs["longname"] = "boundary conditions sensitivity timeseries"  
        
        outds.attrs['Start date'] = start_date
        outds.attrs['End date'] = end_date
        outds.attrs['Latent sampler'] = str(step1)[20:33]
        outds.attrs['Hyper sampler'] = str(step2)[20:33]
        outds.attrs['Burn in'] = str(int(burn))
        outds.attrs['Tuning steps'] = str(int(tune))
        outds.attrs['Number of chains'] = str(int(nchain))
        outds.attrs['Error for each site'] = str(sigma_per_site)
        outds.attrs['Emissions Prior'] = ''.join(['{0},{1},'.format(k, v) for k,v in xprior.items()])[:-1]
        outds.attrs['Model error Prior'] = ''.join(['{0},{1},'.format(k, v) for k,v in sigprior.items()])[:-1]
        outds.attrs['BCs Prior'] = ''.join(['{0},{1},'.format(k, v) for k,v in bcprior.items()])[:-1]
        if add_offset:
            outds.attrs['Offset Prior'] = ''.join(['{0},{1},'.format(k, v) for k,v in offsetprior.items()])[:-1]
        outds.attrs['Creator'] = getpass.getuser()
        outds.attrs['Date created'] = str(pd.Timestamp('today'))
        outds.attrs['Convergence'] = convergence
        outds.attrs['Repository version'] = code_version()
        
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in outds.data_vars}
        output_filename = define_output_filename(outputpath,species,domain,outputname,start_date,ext=".nc")
        Path(outputpath).mkdir(parents=True, exist_ok=True)
        outds.to_netcdf(output_filename, encoding=encoding, mode="w")
