#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:57:23 2020

@author: lw13938
"""
import numpy as np
import acrg_name as name
import pymc3 as pm
import pandas as pd
import xarray as xr
from acrg_grid import areagrid
import getpass
from acrg_hbmcmc.inversionsetup import opends
import os
import acrg_convert as convert

data_path = os.getenv("DATA_PATH")

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

def inferpymc3(Hx, Hbc, Y, error, 
               xprior={"pdf":"lognormal", "mu":1, "sd":1},
               bcprior={"pdf":"lognormal", "mu":0.004, "sd":0.02},
               sigprior={"pdf":"uniform", "lower":0.5, "upper":3},
                nit=2.5e5, burn=50000, tune=1.25e5, nchain=2, verbose=False):       
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
        verbose:
            When True, prints progress bar

            
    Returns:
        outs (array):
            MCMC chain for emissions scaling factors for each basis function.
        bcouts (array):
            MCMC chain for boundary condition scaling factors.
        sigouts (array):
            MCMC chain for model error.
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
       - Allow more than one model-error
    """
    burn = int(burn)         
    
    hx = Hx.T 
    hbc = Hbc.T
    nx = hx.shape[1]
    nbc = hbc.shape[1]
    ny = len(Y)
    
    nit = int(nit)  

    with pm.Model() as model:
        x = parsePrior("x", xprior, shape=nx)
        xbc = parsePrior("xbc", xprior, shape=nbc)
        sig = parsePrior("sig", xprior)
        
        mu = pm.math.dot(hx,x) + pm.math.dot(hbc,xbc)        
        epsilon = pm.math.sqrt(error**2 + sig**2)
        y = pm.Normal('y', mu = mu, sd=epsilon, observed=Y, shape = ny)
        
        step1 = pm.NUTS(vars=[x,xbc])
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
        
        return outs, bcouts, sigouts, convergence, step1, step2

def inferpymc3_postprocessouts(outs,bcouts, sigouts, convergence, 
                               Hx, Hbc, Y, error, 
                               step1, step2, 
                               xprior, bcprior, sigprior,
                               lat, lon, Ytime, siteindicator, data,
                               emissions_name, domain, species, sites,
                               site_lat, site_lon,
                               start_date, end_date, outputname, outputpath,
                               basis_directory, country_directory, fp_basis_case, country_unit_prefix):
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
            outs (array):
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
            lat (array):
                Vector of latitudes at LPDM resolution for the inversion domain.
            lon (array):
                Vector of longitudes at LPDM resolution for the inversion domain.
            Ytime (pandas datetime array):
                Time stamp of measurements as used by the inversion.
            siteindicator (array):
                Numerical indicator of which site the measurements belong to,
                same length at Y.
            site_lon (array):
                Longitude of measurement sites
            site_lat (array):
                Latitude of measurement sites
            data (data array):
                Measurement data from get_obs function.
            emissions_name (dict): 
                Allows emissions files with filenames that are longer than just the species name
                to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
                with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
                multiple sources can be read in simultaneously if they are added as separate entries to
                the emissions_name dictionary.
                If using HiTRes footprints, both the high and low frequency emissions files must be specified
                in a second dictionary like so: {'anth': {'high_freq':'co2-ff-2hr', 'low_freq':'co2-ff-mth'}}.
                It is not a problem to have a mixture of sources, with some that use HiTRes footprints and some
                that don't.
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
            basis_directory (str):
                Directory containing basis function file
            country_directory (str):
                Directory containing country definition file
            fp_basis_case (str, optional):
                Name of basis function to use for emissions.
            country_unit_prefix ('str', optional)
                A prefix for scaling the country emissions. Current options are: 'T' will scale to Tg, 'G' to Gg, 'M' to Mg, 'P' to Pg.
                To add additional options add to acrg_convert.prefix
                Default is none and no scaling will be applied (output in g).
                
                
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
        nit = outs.shape[0]
        nx = Hx.shape[0]
        ny = len(Y)
        nbc = Hbc.shape[0]
        nui = np.arange(2)
        steps = np.arange(nit)
        nmeasure = np.arange(ny)
        nparam = np.arange(nx)
        nBC = np.arange(nbc)
        Ytrace = np.dot(Hx.T,outs.T) + np.dot(Hbc.T,bcouts.T)
        YBCtrace = np.dot(Hbc.T,bcouts.T)
        YmodBC = np.mean(YBCtrace, axis=1)
        Ymod95BC = pm.stats.hpd(YBCtrace.T, 0.95)
        Ymod68BC = pm.stats.hpd(YBCtrace.T, 0.68)
        YaprioriBC = np.sum(Hbc, axis=0)
        Ymod = np.mean(Ytrace, axis=1)
        Ymod95 = pm.stats.hpd(Ytrace.T, 0.95)
        Ymod68 = pm.stats.hpd(Ytrace.T, 0.68)
        Yapriori = np.sum(Hx.T, axis=1) + np.sum(Hbc.T, axis=1)
        sitenum = np.arange(len(sites))

        #Calculate mean posterior scale map and flux field
        if basis_directory is not None:
            bfds = opends(basis_directory+domain+"/"+fp_basis_case+"_"+domain+"_"+start_date[:4]+".nc")
        else:
            bfds = opends(data_path+"/LPDM/basis_functions/"+fp_basis_case+"_"+domain+"_"+start_date[:4]+".nc")
        scalemap = np.zeros_like(np.squeeze(bfds.basis.values))

        for npm in nparam:
            scalemap[bfds.basis.values[:,:,0] == (npm+1)] = np.mean(outs[:,npm])
        if emissions_name == None:
            emds = name.name.flux(domain, species, start = start_date, end = end_date)
        else:
            emds = name.name.flux(domain, list(emissions_name.values())[0], start = start_date, end = end_date)
        flux = scalemap*emds.flux.values[:,:,0]
        
        #Basis functions to save
        bfarray = np.squeeze(bfds.basis.values)-1

        #Calculate country totals   
        area = areagrid(lat, lon)
        c_object = name.get_country(domain, country_dir=country_directory)
        cntryds = xr.Dataset({'country': (['lat','lon'], c_object.country), 
                        'name' : (['ncountries'],c_object.name) },
                                        coords = {'lat': (c_object.lat),
                                        'lon': (c_object.lon)})
        cntrynames = cntryds.name.values
        cntrygrid = cntryds.country.values
        cntrymean = np.zeros((len(cntrynames)))
        cntry68 = np.zeros((len(cntrynames), len(nui)))
        cntry95 = np.zeros((len(cntrynames), len(nui)))
        cntrysd = np.zeros(len(cntrynames))
        molarmass = convert.molar_mass(species)

        unit_factor = convert.prefix(country_unit_prefix)
        if country_unit_prefix is None:
           country_unit_prefix=''
        country_units = country_unit_prefix + 'g'

        # Not sure how it's best to do this if multiple months in emissions 
        # file. Now it scales a weighted average of a priori emissions
        # If a priori emissions have frequency of more than monthly then this
        # needs chaning.
        aprioriflux = np.zeros_like(area)
        if emds.flux.values.shape[2] > 1:
            print("Assuming the inversion is over a year or less and emissions file is monthly")
            allmonths = pd.date_range(start_date, end_date).month[:-1].values
            allmonths -= np.min(allmonths)
            for mi in allmonths:
                aprioriflux += emds.flux.values[:,:,mi]*np.sum(allmonths == mi)/len(allmonths)
        else:
            aprioriflux = np.squeeze(emds.flux.values)
        for ci, cntry in enumerate(cntrynames):
            cntrytottrace = np.zeros(len(steps))
            for bf in range(int(np.max(bfarray))):
                bothinds = np.logical_and(cntrygrid == ci, bfarray==bf)
                cntrytottrace += np.sum(area[bothinds].ravel()*aprioriflux[bothinds].ravel()* \
                               3600*24*365*molarmass)*outs[:,bf]/unit_factor
            cntrymean[ci] = np.mean(cntrytottrace)
            cntry68[ci, :] = pm.stats.hpd(cntrytottrace, 0.68)
            cntry95[ci, :] = pm.stats.hpd(cntrytottrace, 0.95)
        
    
        #Make output netcdf file
        outds = xr.Dataset({'Yobs':(['nmeasure'], Y),
                            'Yerror' :(['nmeasure'], error),                          
                            'Ytime':(['nmeasure'],Ytime),
                            'Yapriori':(['nmeasure'],Yapriori),
                            'Ymodmean':(['nmeasure'], Ymod), 
                            'Ymod95':(['nmeasure','nUI'], Ymod95),
                            'Ymod68':(['nmeasure','nUI'], Ymod68),
                            'YaprioriBC':(['nmeasure'],YaprioriBC),                            
                            'YmodmeanBC':(['nmeasure'], YmodBC),
                            'Ymod95BC':(['nmeasure','nUI'], Ymod95BC),
                            'Ymod68BC':(['nmeasure','nUI'], Ymod68BC),                        
                            'xtrace':(['steps','nparam'], outs),
                            'bctrace':(['steps','nBC'],bcouts),
                            'sigtrace':(['steps'], sigouts),
                            'siteindicator':(['nmeasure'],siteindicator),
                            'sitenames':(['nsite'],sites),
                            'sitelons':(['nsite'],site_lon),
                            'sitelats':(['nsite'],site_lat),
                            'fluxapriori':(['lat','lon'], aprioriflux), #NOTE this is the mean a priori flux over the inversion period
                            'fluxmean':(['lat','lon'], flux),                            
                            'scalingmean':(['lat','lon'],scalemap),
                            'basisfunctions':(['lat','lon'],bfarray),
                            'countrymean':(['countrynames'], cntrymean),
                            'country68':(['countrynames', 'nUI'],cntry68),
                            'country95':(['countrynames', 'nUI'],cntry95),
                            'xsensitivity':(['nmeasure','nparam'], Hx.T),
                            'bcsensitivity':(['nmeasure', 'nBC'],Hbc.T)},
                        coords={'stepnum' : (['steps'], steps), 
                                   'paramnum' : (['nlatent'], nparam),
                                   'numBC' : (['nBC'], nBC),
                                   'measurenum' : (['nmeasure'], nmeasure), 
                                   'UInum' : (['nUI'], nui),
                                   'nsites': (['nsite'], sitenum),
                                   'lat':(['lat'],lat),
                                   'lon':(['lon'],lon),
                                   'countrynames':(['countrynames'],cntrynames)})
        
        outds.fluxmean.attrs["units"] = "mol/m2/s"
        outds.fluxapriori.attrs["units"] = "mol/m2/s"
        outds.Yobs.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Yapriori.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymodmean.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod95.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod68.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.YmodmeanBC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod95BC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod68BC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.YaprioriBC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Yerror.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.countrymean.attrs["units"] = country_units
        outds.country68.attrs["units"] = country_units
        outds.country95.attrs["units"] = country_units
        outds.xsensitivity.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.bcsensitivity.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.sigtrace.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        
        outds.Yobs.attrs["longname"] = "observations"
        outds.Yerror.attrs["longname"] = "measurement error"
        outds.Ytime.attrs["longname"] = "time of measurements"
        outds.Yapriori.attrs["longname"] = "a priori simulated measurements"
        outds.Ymodmean.attrs["longname"] = "mean of posterior simulated measurements"
        outds.Ymod68.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated measurements"
        outds.Ymod95.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated measurements"
        outds.YaprioriBC.attrs["longname"] = "a priori simulated boundary conditions"
        outds.YmodmeanBC.attrs["longname"] = "mean of posterior simulated boundary conditions"
        outds.Ymod68BC.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated boundary conditions"
        outds.Ymod95BC.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated boundary conditions"
        outds.xtrace.attrs["longname"] = "trace of unitless scaling factors for emissions parameters"
        outds.bctrace.attrs["longname"] = "trace of unitless scaling factors for boundary condition parameters"
        outds.sigtrace.attrs["longname"] = "trace of model error parameters"
        outds.siteindicator.attrs["longname"] = "index of site of measurement corresponding to sitenames"
        outds.sitenames.attrs["longname"] = "site names"
        outds.sitelons.attrs["longname"] = "site longitudes corresponding to site names"
        outds.sitelats.attrs["longname"] = "site latitudes corresponding to site names"
        outds.fluxapriori.attrs["longname"] = "mean a priori flux over period"
        outds.fluxmean.attrs["longname"] = "mean posterior flux over period"
        outds.scalingmean.attrs["longname"] = "mean scaling factor field over period"
        outds.basisfunctions.attrs["longname"] = "basis function field"
        outds.countrymean.attrs["longname"] = "mean of ocean and country totals"
        outds.country68.attrs["longname"] = "0.68 Bayesian credible interval of ocean and country totals"
        outds.country95.attrs["longname"] = "0.95 Bayesian credible interval of ocean and country totals"        
        outds.xsensitivity.attrs["longname"] = "emissions sensitivity timeseries"   
        outds.bcsensitivity.attrs["longname"] = "boundary conditions sensitivity timeseries"  
        
        outds.attrs['Start date'] = start_date
        outds.attrs['End date'] = end_date
        outds.attrs['Latent sampler'] = str(step1)[20:33]
        outds.attrs['Hyper sampler'] = str(step2)[20:33]
        outds.attrs['Emissions Prior'] = " ".join([str(x) for x in list(xprior.values())])
        outds.attrs['Model error Prior'] = " ".join([str(x) for x in list(sigprior.values())])
        outds.attrs['BCs Prior'] = " ".join([str(x) for x in list(bcprior.values())])
        outds.attrs['Creator'] = getpass.getuser()
        outds.attrs['Date created'] = str(pd.Timestamp('today'))
        outds.attrs['Convergence'] = convergence
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in outds.data_vars}
        outds.to_netcdf(outputpath+"/"+species.upper()+'_'+domain+'_'+outputname+'_'+start_date+'.nc', encoding=encoding, mode="w")

