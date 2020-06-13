#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:52:03 2020

@author: lw13938
"""
import numpy as np
import scipy.optimize
import xarray as xr
import os 
import getpass
import pandas as pd
import uuid

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
    meanfp /= div
    
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
    tempdir = cwd + "/Temp" + str(uuid.uuid4()) + "/"
    os.mkdir(tempdir)    
    os.mkdir(tempdir + domain)  
#     if not os.path.isdir(cwd+"/Temp"):
#         os.mkdir(cwd+"/Temp")
#     if not os.path.isdir(cwd+"/Temp/"+domain):
#         os.mkdir(cwd+"/Temp/"+domain)
    newds.to_netcdf(tempdir+domain+"/quadtree"+species+"-"+outputname+"_"+domain+"_"+start_date.split("-")[0]+'.nc', mode='w')
    
    return tempdir