#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:13:21 2020

@author: lw13938

Modules for running an MCMC inversion using PyMC3. There are also functions
to dynamically create a basis function grid based on the a priori sensitivity,
and some other functionality for setting up the inputs to this (or any) inverse
method.

If not using on an HPC (i.e. on Snowy), in the terminal you should do:
export OPENBLAS_NUM_THREADS=XX
and/or
export OMP_NUM_THREADS=XX
where XX is the number of chains you are running. If running in Spyer do this
before launching Spyder, else you will use every available thread. Apart from
being annoying it will also slow down your run due to unnecessary forking.

"""
import os
import numpy as np
import acrg_name as name
import acrg_obs as getobs
import socket
import acrg_hbmcmc.inversionsetup as setup 
import acrg_hbmcmc.inversion_pymc3 as mcmc
import acrg_hbmcmc.quadtreebasis as quadtree
import shutil

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

def fixedbasisMCMC(species, sites, domain, meas_period, start_date, 
                   end_date, outputpath, outputname, 
                   xprior={"pdf":"lognormal", "mu":1, "sd":1},
                   bcprior={"pdf":"lognormal", "mu":0.004, "sd":0.02},
                   sigprior={"pdf":"uniform", "lower":0.5, "upper":3},
                   nit=2.5e5, burn=50000, tune=1.25e5, nchain=2,
                   emissions_name=None, inlet=None, fpheight=None, instrument=None, 
                   fp_basis_case=None, bc_basis_case="NESW", 
                   obs_directory = None, country_directory = None,
                   fp_directory = None, bc_directory = None, flux_directory = None,
                   quadtree_basis=True,nbasis=100, 
                   averagingerror=True, bc_freq=None, country_unit_prefix=None,
                   verbose = False):

    """
    Script to run hierarchical Bayesian MCMC for inference of emissions using
    pymc3 to solve the inverse problem.
    
    Args:
        species (str):
            Species of interest
        sites (list):
            List of site names
        domain (str):
            Inversion spatial domain.
        meas_period (list):
            Averaging period of measurements
        start_date (str):
            Start time of inversion "YYYY-mm-dd"
        end_date (str):
            End time of inversion "YYYY-mm-dd"
        outputname (str):
            Unique identifier for output/run name.
        outputpath (str):
            Path to where output should be saved.
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
        nit (int):
            Number of iterations for MCMC
        burn (int):
            Number of iterations to burn in MCMC
        tune (int):
            Number of iterations to use to tune step size
        nchain (int):
            Number of independent chains to run (there is no way at all of 
            knowing whether your distribution has converged by running only
            one chain)    
        emissions_name (dict, optional):
            Allows emissions files with filenames that are longer than just the species name
            to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
            with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
            multiple sources can be read in simultaneously if they are added as separate entries to
            the emissions_name dictionary.
        inlet (str/list, optional):
            Specific inlet height for the site (must match number of sites)
        fpheight (dict, optional):
            Specific release height for the sites' footprints. 
            E.g. fpheight={"TAC":"185m"}(must match number of sites).
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites).
        fp_basis_case (str, optional):
            Name of basis function to use for emissions.
        bc_basis_case (str, optional):
            Name of basis case type for boundary conditions (NOTE, I don't 
            think that currently you can do anything apart from scaling NSEW 
            boundary conditions if you want to scale these monthly.)
        obs_directory (str, optional):
            Directory containing the obs data (with site codes as subdirectories)
            if not default.
        country_directory (str, optional):
            Directory containing the country definition file
        quadtree_basis (bool, optional):
            Creates a basis function file for emissions on the fly using a 
            quadtree algorithm based on the a priori contribution to the mole
            fraction if set to True.
        nbasis (int):
            Number of basis functions that you want if using quadtree derived
            basis function. This will optimise to closest value that fits with
            quadtree splitting algorithm, i.e. nbasis % 4 = 1.
        averagingerror (bool, optional):
            Adds the variability in the averaging period to the measurement 
            error if set to True.
        bc_freq (str, optional):
            The perdiod over which the baseline is estimates. Set to "monthly"
            to estimate per calendar month; set to a number of days,
            as e.g. "30D" for 30 days; or set to None to estimate to have one
            scaling for the whole inversion period.
        country_unit_prefix ('str', optional)
            A prefix for scaling the country emissions. Current options are: 
            'T' will scale to Tg, 'G' to Gg, 'M' to Mg, 'P' to Pg.
            To add additional options add to acrg_convert.prefix
            Default is none and no scaling will be applied (output in g).

            
    Returns:
        Saves an output from the inversion code using inferpymc3_postprocessouts.
        
    TO DO:
        Add a wishlist...
    """    
    data = getobs.get_obs(sites, species, start_date = start_date, end_date = end_date, 
                         average = meas_period, data_directory=obs_directory,
                          keep_missing=False,inlet=inlet, instrument=instrument)
    fp_all = name.footprints_data_merge(data, domain=domain, calc_bc=True, 
                                        height=fpheight, 
                                        fp_directory = fp_directory,
                                        bc_directory = bc_directory,
                                        flux_directory = flux_directory,
                                        emissions_name=emissions_name)
    
    if len(data[sites[0]].mf) == 0:
        print("No observations for %s to %s" % (start_date, end_date))
        return
    if sites[0] not in fp_all.keys():
        print("No footprints for %s to %s" % (start_date, end_date))
        return
    
    print('Running for %s to %s' % (start_date, end_date))
    
    #Add measurement variability in averaging period to measurement error
    if averagingerror:
        fp_all = setup.addaveragingerror(fp_all, sites, species, start_date, end_date,
                                   meas_period, inlet=inlet, instrument=instrument,
                                   obs_directory=obs_directory)
    
    #Create basis function using quadtree algorithm if needed
    if quadtree_basis:
        if fp_basis_case != None:
            print("Basis case %s supplied but quadtree_basis set to True" % fp_basis_case)
            print("Assuming you want to use %s " % fp_basis_case)
        else:
            tempdir = quadtree.quadtreebasisfunction(emissions_name, fp_all, sites, 
                          start_date, domain, species, outputname,
                          nbasis=nbasis)
            fp_basis_case= "quadtree"+species+"-"+outputname
            basis_directory = tempdir
    else:
        basis_directory = None
            
    fp_data = name.fp_sensitivity(fp_all, domain=domain, basis_case=fp_basis_case,basis_directory=basis_directory)
    fp_data = name.bc_sensitivity(fp_data, domain=domain,basis_case=bc_basis_case)
    
    for si, site in enumerate(sites):     
        fp_data[site].attrs['Domain']=domain
    
    #Get inputs ready
    error = np.zeros(0)
    Hbc = np.zeros(0)
    Hx = np.zeros(0)
    Y = np.zeros(0)
    siteindicator = np.zeros(0)
    for si, site in enumerate(sites):
        if 'vmf' in fp_data[site]:           
            error = np.concatenate((error, fp_data[site].vmf.values))
        if 'dmf' in fp_data[site]:
            error = np.concatenate((error, fp_data[site].dmf.values))
            
        Y = np.concatenate((Y,fp_data[site].mf.values)) 
        siteindicator = np.concatenate((siteindicator, np.ones_like(fp_data[site].mf.values)*si))
        if si == 0:
            Ytime=fp_data[site].time.values
        else:
            Ytime = np.concatenate((Ytime,fp_data[site].time.values ))
        
        if bc_freq == "monthly":
            Hmbc = setup.monthly_bcs(start_date, end_date, site, fp_data)
        elif bc_freq == None:
            Hmbc = fp_data[site].H_bc.values
        else:
            Hmbc = setup.create_bc_sensitivity(start_date, end_date, site, fp_data, bc_freq)
            
        if si == 0:
            Hbc = np.copy(Hmbc) #fp_data[site].H_bc.values 
            Hx = fp_data[site].H.values
        else:
            Hbc = np.hstack((Hbc, Hmbc))
            Hx = np.hstack((Hx, fp_data[site].H.values))

        #Run Pymc3 inversion
        xouts, bcouts, sigouts, convergence, step1, step2 = mcmc.inferpymc3(Hx, Hbc, Y, error, 
               xprior,bcprior, sigprior,nit, burn, tune, nchain, verbose=verbose)
        #Process and save inversion output
        mcmc.inferpymc3_postprocessouts(xouts,bcouts, sigouts, convergence, 
                               Hx, Hbc, Y, error, 
                               step1, step2, 
                               xprior, bcprior, sigprior,Ytime, siteindicator, data, fp_data,
                               emissions_name, domain, species, sites,
                               start_date, end_date, outputname, outputpath,
                               basis_directory, country_directory, fp_basis_case, country_unit_prefix)
        
        # remove the temporary basis function directory
        shutil.rmtree(tempdir)
        
        print("All done")
    
