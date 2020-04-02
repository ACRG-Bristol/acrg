#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:13:21 2020

@author: lw13938

Modules for running an MCMC inversion using PyMC3. There are also functions
to dynamically create a basis function grid based on the a priori sensitivity,
and some other functionality for setting up the inputs to this (or any) inverse
method.

It will probably be best to tease this script apart in future, so that the 
method of inference sits in one file and the other functionality sits elsewhere.

If not using on an HPC (i.e. on Snowy), in the terminal you should do:
export OPENBLAS_NUM_THREADS=XX
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

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

def fixedbasisMCMC(species, sites, domain, meas_period, start_date, 
                   end_date, outputpath, outputname, 
                   xprior={"pdf":"lognormal", "mu":1, "sd":1},
                   bcprior={"pdf":"lognormal", "mu":0.004, "sd":0.02},
                   sigprior={"pdf":"uniform", "lower":0.5, "upper":3},
                   nit=2.5e5, burn=50000, tune=1.25e5, nchain=2,
                   emissions_name=None,  height=None, instrument=None, 
                   fp_basis_case=None, bc_basis_case="NESW",  
                   quadtree_basis=True,nbasis=100, 
                   averagingerror=True, bc_monthly=True):
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
        height (str/list, optional):
            Specific inlet height for the site (must match number of sites).
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites).
        fp_basis_case (str, optional):
            Name of basis function to use for emissions.
        bc_basis_case (str, optional):
            Name of basis case type for boundary conditions (NOTE, I don't 
            think that currently you can do anything apart from scaling NSEW 
            boundary conditions if you want to scale these monthly.)
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
        bc_monthly (bool, optional):
            Set to true for the boundary conditions to be scaled each month.
            If set to False, the each boundary has one scaling factor over the
            whole inversion period.
            
    Returns:
        Saves an output from the inversion code using inferpymc3_postprocessouts.
        
    TO DO:
        Add a wishlist...
    """    
    data = getobs.get_obs(sites, species, start_date = start_date, end_date = end_date, 
                         average = meas_period, 
                          keep_missing=False,inlet=height, instrument=instrument)
    fp_all = name.footprints_data_merge(data, domain=domain, calc_bc=True, 
                                        height=height, 
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
                                   meas_period, inlet=height, instrument=instrument)
    
    #Create basis function using quadtree algorithm if needed
    if quadtree_basis:
        if fp_basis_case != None:
            print("Basis case %s supplied but quadtree_basis set to True" % fp_basis_case)
            print("Assuming you want to use %s " % fp_basis_case)
        else:
            quadtree.quadtreebasisfunction(emissions_name, fp_all, sites, 
                          start_date, domain, species, outputname,
                          nbasis=nbasis)
            fp_basis_case= "quadtree"+species+"-"+outputname
            basis_directory = os.getcwd()+"/Temp/"
    else:
        basis_directory = None
            
    fp_data = name.fp_sensitivity(fp_all, domain=domain, basis_case=fp_basis_case,basis_directory=basis_directory)
    fp_data = name.bc_sensitivity(fp_data, domain=domain,basis_case=bc_basis_case)
    
    for si, site in enumerate(sites):     
        fp_data[site].attrs['Domain']=domain
    
    lon = fp_all[sites[0]].lon.values
    lat = fp_all[sites[0]].lat.values
    site_lat = np.zeros(len(sites))
    site_lon = np.zeros(len(sites))
    for si, site in enumerate(sites):
        site_lat[si] = fp_data[site].release_lat.values[0]
        site_lon[si] = fp_data[site].release_lon.values[0]
    
    #Get inputs ready
    error = np.zeros(0)
    Hbc = np.zeros(0)
    Hx = np.zeros(0)
    Y = np.zeros(0)
    siteindicator = np.zeros(0)
    for si, site in enumerate(sites):    
        error = np.concatenate((error, fp_data[site].dmf.values))
        Y = np.concatenate((Y,fp_data[site].mf.values)) 
        siteindicator = np.concatenate((siteindicator, np.ones_like(fp_data[site].mf.values)*si))
        if si == 0:
            Ytime=fp_data[site].time.values
        else:
            Ytime = np.concatenate((Ytime,fp_data[site].time.values ))
        
        if bc_monthly:
            Hmbc = setup.monthly_bcs(start_date, end_date, site, fp_data)
        else:
            Hmbc = fp_data[site].H_bc.values 
            
        if si == 0:
            Hbc = np.copy(Hmbc) #fp_data[site].H_bc.values 
            Hx = fp_data[site].H.values
        else:
            Hbc = np.hstack((Hbc, Hmbc))
            Hx = np.hstack((Hx, fp_data[site].H.values))

        #Run Pymc3 inversion
        xouts, bcouts, sigouts, convergence, step1, step2 = mcmc.inferpymc3(Hx, Hbc, Y, error, 
               xprior,bcprior, sigprior,nit, burn, tune, nchain)
        #Process and save inversion output
        mcmc.inferpymc3_postprocessouts(xouts,bcouts, sigouts, convergence, 
                               Hx, Hbc, Y, error, 
                               step1, step2, 
                               xprior, bcprior, sigprior,
                               lat, lon, Ytime, siteindicator, data,
                               emissions_name, domain, species, sites,
                               site_lat, site_lon,
                               start_date, end_date, outputname, outputpath,
                               basis_directory, fp_basis_case)

        print("All done")
        
if __name__ == "__main__":
    # Pymc3 will use all available threads. Therefore if using Snowy
    # (or non HPC) limit the threads to nchain. This will actually probably
    # speed up run time due to unneccessary forking.
    if socket.gethostname() == 'snowy.chm.bris.ac.uk':
        print("==============================================================")
        print("Before running this I strongly recommend first doing:")
        print("export OPENBLAS_NUM_THREADS=XX")
        print("where XX is the number of chains you are running.")
        print("If running with Spyder, do this before launching Spyder.") 
        print("Otherwise it will use every available thread!")
        print("==============================================================")
    # Run it through to see if it works.
    outputname = "TEST"
    outputpath = "~/Documents/Python/fixed_MCMC/"
    sites = ["MHD"]
    species = "ch4"
    domain = "EUROPE"
    meas_period = ["12H"]
    bc_basis_case = "NESW"
    start_date = "2016-01-01"
    end_date = "2016-03-01"
    fixedbasisMCMC(species, sites, domain, meas_period, start_date, 
                   end_date, outputpath, outputname, nit=5000, burn=1000, 
                   tune=1000)
    
