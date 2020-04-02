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
import pymc3 as pm
import pandas as pd
import xarray as xr
from acrg_grid import areagrid
from acrg_tdmcmc.tdmcmc_post_process import molar_mass
import scipy.optimize
import getpass
import socket

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

def opends(fn):
    '''
    Open a netcdf dataset with xarray
    Args:
        fn (string)
            netcdf file to be opened
    Returns:
        xarray.Dataset: 
            netcdf file as  dataset 
    '''
    with xr.open_dataset(fn) as load:
        ds = load.load()
        return ds

class quadTreeNode:    
    
    def __init__(self, xStart, xEnd, yStart, yEnd):
        self.xStart = xStart
        self.xEnd = xEnd
        self.yStart = yStart
        self.yEnd = yEnd
        
        self.child1 = None #top left
        self.child2 = None #top right
        self.child3 = None #bottom left
        self.child4 = None #bottom right
    
    def isLeaf(self):
        if self.child1 or self.child2 or self.child3 or self.child4:
            return False
        else:
            return True
        
    def createChildren(self, grid, limit):
        value = np.sum(grid[self.xStart:self.xEnd, self.yStart:self.yEnd])#.values

        #stop subdividing if finest resolution or bucket level reached
        if (value < limit or
            (self.xEnd-self.xStart < 2) or (self.yEnd-self.yStart <2)):
            return

        dx = (self.xEnd-self.xStart)
        dy = (self.yEnd-self.yStart)

        #create 4 children for subdivison
        self.child1 = quadTreeNode(self.xStart, self.xStart + dx//2, self.yStart, self.yStart + dy//2)
        self.child2 = quadTreeNode(self.xStart + dx//2, self.xStart + dx, self.yStart, self.yStart + dy//2)
        self.child3 = quadTreeNode(self.xStart, self.xStart + dx//2, self.yStart + dy//2, self.yStart + dy)
        self.child4 = quadTreeNode(self.xStart + dx//2, self.xStart + dx, self.yStart + dy//2, self.yStart + dy)
        
        #apply recursion on all child nodes
        self.child1.createChildren(grid, limit)
        self.child2.createChildren(grid, limit)
        self.child3.createChildren(grid, limit)
        self.child4.createChildren(grid, limit)
        
    def appendLeaves(self, leafList):
        #recursively append all leaves/end nodes to leafList
        if (self.isLeaf()):
            leafList.append(self)
        else:
            self.child1.appendLeaves(leafList)
            self.child2.appendLeaves(leafList)
            self.child3.appendLeaves(leafList)
            self.child4.appendLeaves(leafList)
           
def quadTreeGrid(grid, limit):
    '''
    inputs:
        grid: 2d numpy array to apply quadtree division to
        limit: float to use as bucket level for defining maximum subdivision
    outputs:
        outputGrid: 2d numpy grid, same shape as grid, whose values indicate the box from boxList each index corresponds to
        boxList: list of lists, where each sublist describes the corners of a quadtree leaf
    '''
    #start with a single node the size of the entire input grid:
    parentNode = quadTreeNode(0, grid.shape[0], 0, grid.shape[1])
    parentNode.createChildren(grid, limit)

    leafList = []
    boxList = []
    parentNode.appendLeaves(leafList)
    
    outputGrid = np.zeros_like(grid)

    for i, leaf in enumerate(leafList):
        outputGrid[leaf.xStart:leaf.xEnd, leaf.yStart:leaf.yEnd] = i
        boxList.append([leaf.xStart, leaf.xEnd, leaf.yStart, leaf.yEnd])
    
    return outputGrid, boxList

def quadtreebasisfunction(emissions_name, fp_all, sites, 
                          start_date, domain, species, outputname,
                          nbasis=100):
    """
    Creates a basis function with nbasis grid cells using a quadtree algorithm.
    The domain is split with smaller grid cells for regions which contribute
    more to the a priori (above basline) mole fraction. This is based on the
    average footprint over the inversion period and the a priori emissions field.
    Output is a netcdf file saved to /Temp/<domain> in the current directory.
    The number of basis functions is optimised using dual annealing. Probably
    not the best or fastest method as there should only be one minima, but doesn't
    require the Jacobian or Hessian for optimisation.
    
    Args:
        emissions_name (dict): 
            Allows emissions files with filenames that are longer than just the species name
            to be read in (e.g. co2-ff-mth_EUROPE_2014.nc). This should be a dictionary
            with {source_name: emissions_file_identifier} (e.g. {'anth':'co2-ff-mth'}). This way
            multiple sources can be read in simultaneously if they are added as separate entries to
            the emissions_name dictionary.
        fp_and_data (dict): 
            Output from footprints_data_merge() function. Dictionary of datasets.
        sites (list): 
            List of site names (This could probably be found elsewhere)
        start_date (str): 
            String of start date of inversion
        domain (str): 
            The inversion domain
        species (str): 
            Species of interest
        outputname (str): 
            Identifier or run name
        nbasis (int): 
            Number of basis functions that you want. This will optimise to 
            closest value that fits with quadtree splitting algorithm, 
            i.e. nbasis % 4 = 1.
    
    Returns:
        Nothing. The new basis function is saved to a Temp directory.
    """
    if emissions_name == None:
        meanflux = np.squeeze(fp_all['.flux']['all'].flux.values)
    else:
        meanflux = np.squeeze(fp_all[".flux"][list(emissions_name.keys())[0]].flux.values)
    meanfp = np.zeros((fp_all[sites[0]].fp.shape[0],fp_all[sites[0]].fp.shape[1]))
    div=0
    for site in sites:
        meanfp += np.sum(fp_all[site].fp.values,axis=2)
        div += fp_all[site].fp.shape[2]
    if meanflux.shape != meanfp.shape:
        meanflux = np.mean(meanflux, axis=2)
    fps = meanfp*meanflux

    def qtoptim(x):
        basisQuad, boxes = quadTreeGrid(fps, x)
        return (nbasis - np.max(basisQuad)-1)**2
    optim = scipy.optimize.dual_annealing(qtoptim, np.expand_dims([0,100], axis=0))
    basisQuad, boxes = quadTreeGrid(fps, optim.x[0])
    
    lon = fp_all[sites[0]].lon.values
    lat = fp_all[sites[0]].lat.values    
    
    base = np.expand_dims(basisQuad+1,axis=2)
    
    time = [pd.to_datetime(start_date)]
    newds = xr.Dataset({'basis' : ([ 'lat','lon', 'time'], base)}, 
                        coords={'time':(['time'], time), 
                    'lat' : (['lat'],  lat), 'lon' : (['lon'],  lon)})    
    newds.lat.attrs['long_name'] = 'latitude' 
    newds.lon.attrs['long_name'] = 'longitude' 
    newds.lat.attrs['units'] = 'degrees_north'
    newds.lon.attrs['units'] = 'degrees_east'     
    newds.attrs['creator'] = getpass.getuser()
    newds.attrs['date created'] = str(pd.Timestamp.today())
    cwd = os.getcwd()
    if not os.path.isdir(cwd+"/Temp"):
        os.mkdir(cwd+"/Temp")
    if not os.path.isdir(cwd+"/Temp/"+domain):
        os.mkdir(cwd+"/Temp/"+domain)
    newds.to_netcdf(cwd+"/Temp/"+domain+"/quadtree"+species+"-"+outputname+"_"+domain+"_"+start_date.split("-")[0]+'.nc', mode='w')

def addaveragingerror(fp_all, sites, species, start_date, end_date, meas_period,  
                      inlet=None, instrument=None):
    """
    Adds the variablility within the averaging period to the mole fraction error.
    
    Args:
        fp_all (dict):
            Output from footprint_data_merge
        sites (list):
            List of site names
        species (str):
            Species of interest
        start_date (str):
            Start time of inversion "YYYY-mm-dd"
        end_date (str):
            End time of inversion "YYYY-mm-dd"
        meas_period (list):
            Averaging period of measurements
        inlet (str/list, optional):
            Specific inlet height for the site (must match number of sites).
        instrument (str/list, optional):
            Specific instrument for the site (must match number of sites).
            
    Returns:
        fp_all (dict):
            fp_all from input with averaging error added to the mole fraction
            error
    """
    #Add variability in measurement averaging period to repeatability 
    dataerr = getobs.get_obs(sites, species, start_date = start_date, end_date = end_date,  
                          keep_missing=False,inlet=inlet, instrument=instrument)
    for si, site in enumerate(sites):
        if min(dataerr[site].index) > pd.to_datetime(start_date):
            dataerr[site].loc[pd.to_datetime(start_date)] = \
                [np.nan for col in dataerr[site].columns]           
        # Pad with an empty entry at the end date
        if max(dataerr[site].index) < pd.to_datetime(end_date):
            dataerr[site].loc[pd.to_datetime(end_date)] = \
                [np.nan for col in dataerr[site].columns]
        # Now sort to get everything in the right order
        dataerr[site] = dataerr[site].sort_index()
        fp_all[site].dmf.values = fp_all[site].dmf.values + dataerr[site].mf.resample(meas_period[si]).std(ddof=0).dropna().values
    return fp_all

def monthly_bcs(start_date, end_date, site, fp_data):
    """
    Creates a sensitivity matrix (H-matrix) for the boundary conditions, which
    will map monthly boundary condition scalings to the observations. This is 
    for a single site.
    
    Args:
        start_date (str):
            Start time of inversion "YYYY-mm-dd"
        end_date (str):
            End time of inversion "YYYY-mm-dd"
        site (str):
            Site that you're creating it for
        fp_data (dict):
            Output from acrg_name.bc_sensitivity
            
    Returns:
        Hmbc (array):
            Sensitivity matrix by month for observations
            
    """
    allmonth= pd.date_range(start_date, end_date, freq="MS")[:-1] 
    nmonth = len(allmonth)
    curtime = pd.to_datetime(fp_data[site].time.values).to_period("M")
    pmonth = pd.to_datetime(fp_data[site].resample(time="MS").mean().time.values)
    Hmbc = np.zeros((4*nmonth, len(fp_data[site].time.values)))
    cnt=0
    for cord in range(4):
        for m in range(0,nmonth):
            if allmonth[m] not in pmonth:
                cnt += 1
                continue
            mnth = allmonth[m].month
            yr = allmonth[m].year
            mnthloc = np.where(np.logical_and(curtime.month == mnth,curtime.year == yr))[0]
            Hmbc[cnt,mnthloc] = fp_data[site].H_bc.values[cord,mnthloc] 
            cnt += 1
    return Hmbc

def inferpymc3(Hx, Hbc, Y, error, 
               xprior={"pdf":"lognormal", "mu":1, "sd":1},
               bcprior={"pdf":"lognormal", "mu":0.004, "sd":0.02},
               sigprior={"pdf":"uniform", "lower":0.5, "upper":3},
                nit=2.5e5, burn=50000, tune=1.25e5, nchain=2):       
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
       - Allow any type of sampler to sample variables
       - Allow non-iid variables
       - Allow more than one model-error
       - Neaten up the way the priors are assigned.
    """
    burn = int(burn)         
    
    hx = Hx.T 
    hbc = Hbc.T
    nx = hx.shape[1]
    nbc = hbc.shape[1]
    ny = len(Y)
    
    nit = int(nit)  

    with pm.Model() as model:
        if xprior["pdf"].lower() == "uniform":
            x = pm.Uniform('x',lower=xprior["lower"],upper=xprior["upper"], shape=nx)
        elif xprior["pdf"].lower() == "lognormal":
            x = pm.Lognormal('x', mu=xprior["mu"], sd=xprior["sd"], shape=nx)
        elif xprior["pdf"].lower() == "halfflat":
            x = pm.HalfFlat('x', shape=nx)
        else:
            print("Haven't coded in %s yet for emissions" % xprior["pdf"])

        if bcprior["pdf"].lower() == "uniform":
            xbc = pm.Uniform('xbc',lower=bcprior["lower"],upper=bcprior["upper"], shape=nbc)
        elif bcprior["pdf"].lower() == "lognormal":
            xbc = pm.Lognormal('xbc', mu=bcprior["mu"], sd=bcprior["sd"], shape=nbc)
        elif bcprior["pdf"].lower() == "halfflat":
            xbc = pm.HalfFlat('xbc', shape=nbc)
        else:
            print("Haven't coded in %s yet for boundary conditions" % bcprior["pdf"])  
        
        if sigprior["pdf"].lower() == "uniform":
            sig = pm.Uniform('sig',lower=sigprior["lower"],upper=sigprior["upper"])
        elif bcprior["pdf"].lower() == "lognormal":
            sig = pm.Lognormal('sig', mu=sigprior["mu"], sd=sigprior["sd"])
        elif sigprior["pdf"].lower() == "halfflat":
            sig = pm.HalfFlat('sig')
        else:
            print("Haven't coded in %s yet for model error" % sigprior["pdf"])  
        
        mu = pm.math.dot(hx,x) + pm.math.dot(hbc,xbc)        
        epsilon = pm.math.sqrt(error**2 + sig**2)
        y = pm.Normal('y', mu = mu, sd=epsilon, observed=Y, shape = ny)
        
        step1 = pm.NUTS(vars=[x,xbc])
        step2 = pm.Slice(vars=[sig])
        
        trace = pm.sample(nit, tune=int(tune), chains=nchain,
                           step=[step1,step2]) #, progressbar=False)#step=pm.Metropolis())#  #target_accept=0.8,
        
        outs = trace.get_values(x, burn=burn)[0:int((nit)-burn)]
        bcouts = trace.get_values(xbc, burn=burn)[0:int((nit)-burn)]
        sigouts = trace.get_values(sig, burn=burn)[0:int((nit)-burn)]
        
        #Check for convergence
        gelrub = pm.gelman_rubin(trace)['x'].max()
        if gelrub > 1.05:
            print('Failed Gelman-Rubin at 1.05 (%s)' % (round(gelrub,2)))
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
                               basis_directory, fp_basis_case):
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
            fp_basis_case (str, optional):
                Name of basis function to use for emissions.
                
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
        YmodBC = np.mean(np.dot(Hbc.T,bcouts.T), axis=1)
        YaprioriBC = np.sum(Hbc, axis=0)
        Ymod = np.mean(Ytrace, axis=1)
        Ymod95 = pm.stats.hpd(Ytrace.T, 0.05)
        Ymod68 = pm.stats.hpd(Ytrace.T, 0.32)
        Yapriori = np.sum(Hx.T, axis=1) + np.sum(Hbc.T, axis=1)
        sitenum = np.arange(len(sites))

        #Calculate mean posterior scale map and flux field'
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
        c_object = name.get_country(domain)
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
        molarmass = molar_mass(species)
        # Not sure how it's best to do this if multiple months in emissions 
        # file. Now it scales a weighted average of a priori emissions
        # If a priori emissions have frequency of more than monthly then this
        # needs chaning.
        aprioriflux = np.zeros_like(area)
        if emds.flux.values.shape[2] > 1:
            print("Assuming the inversion is over a year or less and emissions file is monthly")
            allmonths = pd.date_range(start_date, end_date).month[:-1]
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
                               3600*24*365*molarmass)*outs[:,bf]/1e9
            cntrymean[ci] = np.mean(cntrytottrace)
            cntry68[ci, :] = pm.stats.hpd(cntrytottrace, 0.32)
            cntry95[ci, :] = pm.stats.hpd(cntrytottrace, 0.05)
            cntrysd[ci] = np.std(cntrytottrace)
        
    
        #Make output netcdf file
        outds = xr.Dataset({'Y':(['nmeasure'], Y),
                            'Ytime':(['nmeasure'],Ytime),
                            'Yapriori':(['nmeasure'],Yapriori),
                            'Ymod':(['nmeasure'], Ymod), 
                            'Ymod95':(['nmeasure','nUI'], Ymod95),
                            'Ymod68':(['nmeasure','nUI'], Ymod68),
                            'YmodBC':(['nmeasure'], YmodBC),
                            'YaprioriBC':(['nmeasure'],YaprioriBC),
                            'xtrace':(['steps','nparam'], outs),
                            'bctrace':(['steps','nBC'],bcouts),
                            'sigtrace':(['steps'], sigouts),
                            'Yerror' :(['nmeasure'], error),
                            'siteindicator':(['nmeasure'],siteindicator),
                            'sitenames':(['nsite'],sites),
                            'site_lon':(['nsite'],site_lon),
                            'site_lat':(['nsite'],site_lat),
                            'aprioriflux':(['lat','lon'], aprioriflux), #NOTE this is the mean a priori flux over the inversion period 
                            'meanscaling':(['lat','lon'],scalemap),
                            'meanflux':(['lat','lon'], flux),
                            'basis_functions':(['lat','lon'],bfarray),
                            'countrytotals':(['countrynames'], cntrymean),
                            'country68':(['countrynames', 'nUI'],cntry68),
                            'country95':(['countrynames', 'nUI'],cntry95),
                            'countrysd':(['countrynames'],cntrysd)},
                        coords={'stepnum' : (['steps'], steps), 
                                   'paramnum' : (['nlatent'], nparam),
                                   'numBC' : (['nBC'], nBC),
                                   'measurenum' : (['nmeasure'], nmeasure), 
                                   'UInum' : (['nUI'], nui),
                                   'nsites': (['nsite'], sitenum),
                                   'lat':(['lat'],lat),
                                   'lon':(['lon'],lon),
                                   'countrynames':(['countrynames'],cntrynames)})
        outds.meanflux.attrs["units"] = "mol/m2/s"
        outds.Y.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Yapriori.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod95.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Ymod68.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.YmodBC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.YaprioriBC.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.Yerror.attrs["units"] = str(data[".units"])+" "+"mol/mol"
        outds.countrytotals.attrs["units"] = "Gg/yr"
        outds.country68.attrs["units"] = "Gg/yr"
        outds.country95.attrs["units"] = "Gg/yr"
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
        fp_all = addaveragingerror(fp_all, sites, species, start_date, end_date,
                                   meas_period, inlet=height, instrument=instrument)
    
    #Create basis function using quadtree algorithm if needed
    if quadtree_basis:
        if fp_basis_case != None:
            print("Basis case %s supplied but quadtree_basis set to True" % fp_basis_case)
            print("Assuming you want to use %s " % fp_basis_case)
        else:
            quadtreebasisfunction(emissions_name, fp_all, sites, 
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
            Hmbc = monthly_bcs(start_date, end_date, site, fp_data)
        else:
            Hmbc = fp_data[site].H_bc.values 
            
        if si == 0:
            Hbc = np.copy(Hmbc) #fp_data[site].H_bc.values 
            Hx = fp_data[site].H.values
        else:
            Hbc = np.hstack((Hbc, Hmbc))
            Hx = np.hstack((Hx, fp_data[site].H.values))

        #Run Pymc3 inversion
        xouts, bcouts, sigouts, convergence, step1, step2 = inferpymc3(Hx, Hbc, Y, error, 
               xprior,bcprior, sigprior,nit, burn, tune, nchain)
        #Process and save inversion output
        inferpymc3_postprocessouts(xouts,bcouts, sigouts, convergence, 
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
    
