# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 15:46:31 2015

@author: ag12733

Creat basis functions for domain
"""
from __future__ import print_function

import numpy as np
import xarray as xray
import glob
import pandas as pd
import acrg_name as name
import os
from os.path import join

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

if acrg_path is None:
    acrg_path = os.getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")
if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

# These are the default directories if no optional arguments are specified in footprints_data_merge,
# bc_sensitivity or fp_sensitivity
fields_file_path = join(data_path, 'NAME/fp/')
basis_dir = join(data_path, 'NAME/basis_functions/')
bc_basis_dir = join(data_path,'NAME/bc_basis_functions/')


def basis_blocks(domain, time, blocksize, basis_case=None): 
    """Creates basis functions in square blocks (e.g., 5x5 grid cells per region)
    
    Args:
        domain (str): 
            String of domain area
        time (str): 
            Timestamp for basis functions
        blocksize (int): 
            Number of grid cells to split into (i.e. blocksize x blocksize)
        basis_case (str, optional): 
            Labelling for basis case. 
            Default is e.g. '5x5' for blocksize of 5
    
    Returns:
        None
        Writes blocked basis functions to netcdf
        
    Example:
        basis_blocks(domain = "SOUTHASIA", time = "2012-01-01", blocksize=25)
    """
    if basis_case == None:
        basis_case = str(blocksize)+"x"+str(blocksize)
    
    # load a Fields file to get the domain

    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    basis_grid = np.zeros((len(lat),len(lon)))
    
    arr = basis_grid
    
    x,y = cut_array2d(arr,(len(lat)/blocksize,len(lon)/blocksize))
    
    for ii in range(len(x)):
        basis_grid[x[ii][0]:x[ii][1],y[ii][0]:y[ii][1]]=ii+1
    
    basis_grid = basis_grid.reshape(len(lat),len(lon),1)
    
    basis_ds = xray.Dataset({'basis':(['lat','lon','time'], basis_grid)}, \
            coords = {'lon':(['lon'], lon), \
                      'lat':(['lat'], lat), \
                      'time':(['time'],[time])})
                
    basis_ds.to_netcdf(basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc")

def basis_transd(domain, time, basis_case = "sub-transd", sub_lon_min = None,
                 sub_lon_max = None, sub_lat_min = None, sub_lat_max = None): 
    """Creates basis regions for the "transdimensional inversion".
    Creates variable central block based on sub_lat and sub_lon domain
    and four fixed regions surrounding the central box. 
    
    Args:
        domain (str): 
            String of domain area
        time (str): 
            Timestamp for basis functions
        basis_case (str, optional):
            Labelling for basis case. 
            Default is "transd" for transdimensional
        sub_lon_min (float): 
            Minimum longitude of sub domain
        sub_lon_max (float): 
            Maximum longitude of sub domain
        sub_lat_min (float): 
            Minimum latitude of sub domain
        sub_lat_max (float): 
            Maximum latitude of sub domain
    
    Returns:
        None
        Writes basis regions to netcdf
        
    Example:
        basis_transd(domain = "SOUTHASIA", time = "2012-01-01", basis_case = "transd",
            sub_lon_min = 65, sub_lon_max = 100, sub_lat_min = 5, sub_lat_max = 40)
    """  
    
    if None in [sub_lon_min, sub_lon_max, sub_lat_min, sub_lat_max]:
        print("At least one sub domain lon or lat is missing on input")
        return        
    
#     load a Fields file to get the domain

    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    basis_grid = np.zeros((len(lat),len(lon)))    
    
    sub_lon_min_arg = np.abs(lon - sub_lon_min).argmin()
    sub_lon_max_arg = np.abs(lon - sub_lon_max).argmin()
    sub_lat_min_arg = np.abs(lat - sub_lat_min).argmin()
    sub_lat_max_arg = np.abs(lat - sub_lat_max).argmin()
    
    sub_lon = lon[sub_lon_min_arg:sub_lon_max_arg]
    sub_lat = lat[sub_lat_min_arg:sub_lat_max_arg]   
     
    basis_grid[0:sub_lat_min_arg,:] = 1
    basis_grid[sub_lat_min_arg:,0:sub_lon_min_arg] = 2
    basis_grid[sub_lat_max_arg:,sub_lon_min_arg:sub_lon_max_arg] = 3
    basis_grid[sub_lat_min_arg:,sub_lon_max_arg:] = 4
   
    basis_grid = basis_grid.reshape(len(lat),len(lon),1)
    
    basis_ds = xray.Dataset({'basis':(['lat','lon','time'], basis_grid)}, \
            coords = {'lon':(['lon'], lon), \
                      'lat':(['lat'], lat), \
                      'time':(['time'],[time]) ,\
                      'sub_lat': (['sub_lat'],sub_lat),\
                      'sub_lon': (['sub_lon'],sub_lon)})
                
    basis_ds.to_netcdf(basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc")
    

def basis_bc_blocks(domain, time, basis_case = "NESW", vertical=1):
    """Creates uniform blocks for boundary conditions in each direction (NESW).
    Each direction is split into vertical slabs with number specified by input 'vertical'
    
    Args:
        domain (str):
            String of domain area
        basis_case (str): 
            Labelling for basis case. 
        time (str): 
            Timestamp for basis functions
        vertical (int, optioanl): 
            Number of vertical slabs to split boundary conditions into
            Default = 1
            
    Returns:
        None
        Writes basis boundary condition blocks to netcdf
        
    Example:
        basis_bc_blocks(domain='SOUTHASIA', basis_case='NESW', time='2012-01-01', vertical = 4)
    """
        
    files = glob.glob(fields_file_path + domain + "/*")

    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    regions = range(vertical*4)

    heights = range(500,20500,1000)
    
    basis_grid_e = np.zeros((len(heights),len(lat),len(regions),1))
    basis_grid_w = np.zeros((len(heights),len(lat),len(regions),1))
    basis_grid_n = np.zeros((len(heights),len(lon),len(regions),1))
    basis_grid_s = np.zeros((len(heights),len(lon),len(regions),1))
    
    vertical_index = np.round(np.linspace(0,len(heights), vertical+1)).astype('int')
    
    for ii in range(vertical):
        
        basis_grid_n[vertical_index[ii]:vertical_index[ii+1],:,0+ii,:] = 1    
        basis_grid_e[vertical_index[ii]:vertical_index[ii+1],:,vertical+ii,:] = 1  
        basis_grid_s[vertical_index[ii]:vertical_index[ii+1],:,2*vertical+ii,:] = 1
        basis_grid_w[vertical_index[ii]:vertical_index[ii+1],:,3*vertical+ii,:] = 1
        
    basis_ds = xray.Dataset({'bc_basis_n':(['height','lon','region','time'], basis_grid_n), \
                            'bc_basis_s':(['height','lon','region','time'], basis_grid_s), \
                            'bc_basis_e':(['height','lat','region','time'], basis_grid_e), \
                            'bc_basis_w':(['height','lat','region','time'], basis_grid_w)}, \
                coords = {'height': (['height'], heights), \
                          'lon':(['lon'], lon), \
                          'lat':(['lat'], lat), \
                          'time':(['time'],[time])})
    
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc")                
    
def basis_bc_uniform(domain, time, basis_case = "uniform"):
    """Creates creates one uniform scaling for boundary conditions for all directions
    
    Args:
        domain (str): 
            String of domain area
        basis_case (str): 
            Labelling for basis case. 
        time (str): 
            Timestamp for basis functions
            
    Returns:
        None
        Writes uniform basis boundary condition to netcdf
        
    Example:
        basis_bc_uniform(domain='SOUTHASIA', basis_case='uniform', time='2012-01-01')
    """
    
 

    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    regions = range(1)
    heights = range(500,20500,1000)
    
    basis_grid_e = np.zeros((len(heights),len(lat),len(regions),1))
    basis_grid_w = np.zeros((len(heights),len(lat),len(regions),1))
    basis_grid_n = np.zeros((len(heights),len(lon),len(regions),1))
    basis_grid_s = np.zeros((len(heights),len(lon),len(regions),1))
    
            
    basis_grid_n[:,:,:,:] = 1    
    basis_grid_e[:,:,:,:] = 1  
    basis_grid_s[:,:,:,:] = 1
    basis_grid_w[:,:,:,:] = 1
        
    basis_ds = xray.Dataset({'bc_basis_n':(['height','lon','region','time'], basis_grid_n), \
                            'bc_basis_s':(['height','lon','region','time'], basis_grid_s), \
                            'bc_basis_e':(['height','lat','region','time'], basis_grid_e), \
                            'bc_basis_w':(['height','lat','region','time'], basis_grid_w)}, \
                coords = {'height': (['height'], heights), \
                          'lon':(['lon'], lon), \
                          'lat':(['lat'], lat), \
                          'time':(['time'],[time])})
                
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + time.strftime('%m')+str(time.year) + ".nc")

def basis_bc_all_gradients(domain, time, species, units='ppb', basis_case='horiz-strat-grad'):
    """Creates five terms for bc basis regions
        (1) offset to shift entire field up and down
        (2) factor to scale the entire field
        (3) scaling of N-S gradient
        (4) scaling of E-W gradient
        (5) scaling of stratospheric gradient
        Technically file is species specific as it calculates the stratosphere 
        based on the vmm but one can assume/hope that the stratosphere doesn't 
        change by species 
        
    Args:
        domain (str): 
            String of domain area
        time (str): 
            Timestamp for basis functions
        species (str): 
            Species of interest
        units (str, optional): 
            Units of boundary conditions. Default ppb
        basis_case (str, optional): 
            Name of basis case. 
            Default 'horiz-strat-grad'
            
    Returns:
        None
        Writes basis boundary condition gradients to netcdf
        
    Example:
        basis_bc_all_gradients('SOUTHASIA', '2012-01-01', 'ch4', units='ppb', basis_case='horiz-strat-grad')
        
    """
    
    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values

    heights = range(500,20500,1000)
    
    bc_ds =  name.boundary_conditions(domain, species)
             
    bc_ds = bc_ds.sel(time = time)

    if units is 'ppb':
        bc_ds = bc_ds*1e9
    elif units is 'ppm':
        bc_ds = bc_ds*1e6
        
    vmr_n = bc_ds["vmr_n"].values
    vmr_e = bc_ds["vmr_e"].values
    vmr_s = bc_ds["vmr_s"].values
    vmr_w = bc_ds["vmr_w"].values

    #flip the east and south directions to make it continuous
    
    vmr_e_fl = np.fliplr(vmr_e)
    vmr_s_fl = np.fliplr(vmr_s)
    
    # stick together
    
    ind = [None]*4
    vmr_all = np.concatenate((vmr_n, vmr_e_fl, vmr_s_fl, vmr_w), axis = 1)
    
    # add component that will scale entire field up and down by some amount
    
    scale_n = 1*np.ones(np.shape(vmr_n))
    scale_e = 1*np.ones(np.shape(vmr_e))
    scale_s = 1*np.ones(np.shape(vmr_s))
    scale_w = 1*np.ones(np.shape(vmr_w))
    scale_all = np.concatenate((scale_n, scale_e, scale_s, scale_w), axis = 1)
    
    # find indices that correspond to each direction
    ind[0] = range(len(lon)) #N
    ind[1] = range(max(ind[0])+1,max(ind[0])+1 + len(lat)) #E
    ind[2] = range(max(ind[1])+1,max(ind[1])+1 + len(lon)) #s
    ind[3] = range(max(ind[2])+1,max(ind[2])+1 + len(lat)) #w
    
    # create field of East-West gradient
    grad_e_w = np.zeros(np.shape(vmr_n))
    slope_e_w = (10.0 - -10.0)/len(lon)
    for ii in range(len(lon)):
        grad_e_w[:,ii] = slope_e_w*ii + -10
 
    longrad = np.zeros(np.shape(vmr_all))
    longrad[:,ind[0]] = grad_e_w
    longrad[:,ind[2]] = np.fliplr(grad_e_w)

    # create field of North-South gradient
    grad_n_s = np.zeros(np.shape(vmr_e))
    slope_n_s = (10.0 - -10.0)/len(lat)
    for ii in range(len(lat)):
        grad_n_s[:,ii] = slope_n_s*ii + -10
 
    latgrad = np.zeros(np.shape(vmr_all))
    latgrad[:,ind[1]] = grad_n_s
    latgrad[:,ind[3]] = np.fliplr(grad_n_s)

    # find maximum rate of change (second derivative) between levels 10-15

    secderiv = np.diff(vmr_all, n=2, axis=0)
    tropopause = np.argmin(secderiv[10:15,:], axis=0)+10+1 # add 1 as compromise. Because of finite differencing, levels will be decreased by 2

        
    stratgradient =  np.zeros(np.shape(vmr_all)) 
    
    for ii in range(2*(len(lat)+len(lon))):
        stratgradient[tropopause[ii]:,ii] = -20*np.array(range(len(heights)-tropopause[ii]))      
        
    pcs = []
    pcs.append(scale_all)
    pcs.append(vmr_all)
    pcs.append(latgrad)
    pcs.append(longrad)
    pcs.append(stratgradient)   

                
        
    numpcs = len(pcs)
    regions = range(numpcs*4)
    
    basis_grid_n = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_e = np.zeros((len(heights),len(lat),len(regions)))
    basis_grid_s = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_w = np.zeros((len(heights),len(lat),len(regions)))


    # unravel pcs into directions
    pc_n = []
    pc_e = []
    pc_s = []
    pc_w = []
    
    for ii in range(numpcs):
        pc_n.append(pcs[ii][:,ind[0]])
        pc_e.append(pcs[ii][:,ind[1]])
        pc_s.append(pcs[ii][:,ind[2]])
        pc_w.append(pcs[ii][:,ind[3]])

    # reverse pcs for east and south to make same as VMRs
    for ii in range(numpcs):
        pc_e[ii] = np.fliplr(pc_e[ii])
        pc_s[ii] = np.fliplr(pc_s[ii])
                
    # create basis grids in format needed for inversion
        
    for ii in range(numpcs):
        basis_grid_n[:,:,ii] = pc_n[ii]/vmr_n
        basis_grid_e[:,:,ii+3] = pc_e[ii]/vmr_e
        basis_grid_s[:,:,ii+6] = pc_s[ii]/vmr_s
        basis_grid_w[:,:,ii+9] = pc_w[ii]/vmr_w
        
    basis_grid_n = basis_grid_n[:,:,:,None]
    basis_grid_e = basis_grid_e[:,:,:,None]
    basis_grid_s = basis_grid_s[:,:,:,None]
    basis_grid_w = basis_grid_w[:,:,:,None]
    

    basis_ds = xray.Dataset({'bc_basis_n':(['height','lon','region','time'], basis_grid_n), \
                            'bc_basis_s':(['height','lon','region','time'], basis_grid_s), \
                            'bc_basis_e':(['height','lat','region','time'], basis_grid_e), \
                            'bc_basis_w':(['height','lat','region','time'], basis_grid_w)}, \
                coords = {'height': (['height'], heights), \
                          'lon':(['lon'], lon), \
                          'lat':(['lat'], lat), \
                          'time':(['time'],[time])})
                
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + time.strftime('%m')+str(time.year) + ".nc")           

def basis_bc_horiz_gradients(domain, time, basis_case='horiz-grad'):
    """Creates three terms for bc basis regions for each direction
        (1) offset to shift field up and down
        (2) scaling of N-S gradient on east and west curtains
        (3) scaling of E-W gradient on north and south curtains
        
    Args:
        domain (str): 
            String of domain area
        time (str): 
            Timestamp for basis functions
        basis_case (str, optional): 
            Name of basis case. 
            Default 'horiz-grad'
            
    Returns:
        None
        Writes basis boundary condition gradients to netcdf
        
    Example:
        basis_bc_all_gradients('SOUTHASIA', '2012-01-01')
        
    """
    
    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values

    heights = range(500,20500,1000)
    
    # add component that will scale entire field up and down by some amount
    
    scale_n = 1*np.ones((len(heights), len(lon)))
    scale_e = 1*np.ones((len(heights), len(lat)))
    scale_s = 1*np.ones((len(heights), len(lon)))
    scale_w = 1*np.ones((len(heights), len(lat)))
    scale_all = np.concatenate((scale_n, scale_e, scale_s, scale_w), axis = 1)
    
    # find indices that correspond to each direction
    ind = [None]*4
    ind[0] = range(len(lon)) #N
    ind[1] = range(max(ind[0])+1,max(ind[0])+1 + len(lat)) #E
    ind[2] = range(max(ind[1])+1,max(ind[1])+1 + len(lon)) #s
    ind[3] = range(max(ind[2])+1,max(ind[2])+1 + len(lat)) #w
    
    # create field of East-West gradient
    grad_e_w = np.zeros((len(heights), len(lon)))
    slope_e_w = (1.0 - -1.0)/len(lon)
    for ii in range(len(lon)):
        grad_e_w[:,ii] = slope_e_w*ii + -1
 
    longrad = np.zeros((len(heights), 2*len(lon) + 2*len(lat)))
    longrad[:,ind[0]] = grad_e_w
    longrad[:,ind[2]] = np.fliplr(grad_e_w)

    # create field of North-South gradient
    grad_n_s = np.zeros((len(heights), len(lat)))
    slope_n_s = (1.0 - -1.0)/len(lat)
    for ii in range(len(lat)):
        grad_n_s[:,ii] = slope_n_s*ii + -1
 
    latgrad = np.zeros((len(heights), 2*len(lon) + 2*len(lat)))
    latgrad[:,ind[1]] = grad_n_s
    latgrad[:,ind[3]] = np.fliplr(grad_n_s)
        
    pcs = []
    pcs.append(scale_all)
    pcs.append(latgrad)
    pcs.append(longrad)
 

    numpcs = len(pcs)
    regions = range(numpcs*4)
    
    basis_grid_n = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_e = np.zeros((len(heights),len(lat),len(regions)))
    basis_grid_s = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_w = np.zeros((len(heights),len(lat),len(regions)))


    # unravel pcs into directions
    pc_n = []
    pc_e = []
    pc_s = []
    pc_w = []
    
    for ii in range(numpcs):
        pc_n.append(pcs[ii][:,ind[0]])
        pc_e.append(pcs[ii][:,ind[1]])
        pc_s.append(pcs[ii][:,ind[2]])
        pc_w.append(pcs[ii][:,ind[3]])

    # reverse pcs for east and south to make same as VMRs
    for ii in range(numpcs):
        pc_e[ii] = np.fliplr(pc_e[ii])
        pc_s[ii] = np.fliplr(pc_s[ii])
                
    # create basis grids in format needed for inversion
        
    for ii in range(numpcs):
        basis_grid_n[:,:,ii] = pc_n[ii]
        basis_grid_e[:,:,ii+3] = pc_e[ii]
        basis_grid_s[:,:,ii+6] = pc_s[ii]
        basis_grid_w[:,:,ii+9] = pc_w[ii]
        
    basis_grid_n = basis_grid_n[:,:,:,None]
    basis_grid_e = basis_grid_e[:,:,:,None]
    basis_grid_s = basis_grid_s[:,:,:,None]
    basis_grid_w = basis_grid_w[:,:,:,None]
    

    basis_ds = xray.Dataset({'bc_basis_n':(['height','lon','region','time'], basis_grid_n), \
                            'bc_basis_s':(['height','lon','region','time'], basis_grid_s), \
                            'bc_basis_e':(['height','lat','region','time'], basis_grid_e), \
                            'bc_basis_w':(['height','lat','region','time'], basis_grid_w)}, \
                coords = {'height': (['height'], heights), \
                          'lon':(['lon'], lon), \
                          'lat':(['lat'], lat), \
                          'time':(['time'],[time])})
    
    basis_ds.bc_basis_n.encoding = {'zlib':True}       
    basis_ds.bc_basis_s.encoding = {'zlib':True}  
    basis_ds.bc_basis_e.encoding = {'zlib':True}  
    basis_ds.bc_basis_w.encoding = {'zlib':True}                                   
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc", mode='w')           
 
               
def basis_bc_pca(domain, time, species, units='ppb', basis_case='pca', numregions = 4):
    """breaks the MOZART VMR into the first 'numregions' principal components and scales each PC
    
    Args:
        domain (str): 
            String of domain area
        time (str): 
            Timestamp for basis functions
        species (str): 
            Species of interest
        units (str, optional): 
            Units of species. Default is 'ppb' 
        basis_case (str):
            Labelling for basis case. Default is 'pca'
        numregions (int):
            Number of PCs 
            
    Returns:
        None
        Writes PCA basis boundary condition to netcdf
        
    Example:
        basis_bc_pca("SOUTHASIA", "2012-01-15", 'ch4', units='ppb', basis_case='pca', numregions = 4)
    """
#   
# will be species specific
    
    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    regions = range(numregions)
    heights = range(500,20500,1000)
    
    basis_grid_n = np.zeros((len(heights),len(lon),numregions))
    basis_grid_e = np.zeros((len(heights),len(lat),numregions))
    basis_grid_s = np.zeros((len(heights),len(lon),numregions))
    basis_grid_w = np.zeros((len(heights),len(lat),numregions))

    # find principal companents of MOZART vmr

    bc_ds =  name.boundary_conditions(domain, species)
    bc_ds = bc_ds.sel(time = time)
    
    if units is 'ppb':
        bc_ds = bc_ds*1e9
    elif units is 'ppm':
        bc_ds = bc_ds*1e6
        
    vmr_n = bc_ds["vmr_n"].values
    vmr_e = bc_ds["vmr_e"].values
    vmr_s = bc_ds["vmr_s"].values
    vmr_w = bc_ds["vmr_w"].values
    
    #flip the east and south directions to make continuous
    
    vmr_e_fl = np.fliplr(vmr_e)
    vmr_s_fl = np.fliplr(vmr_s)
    
    # stick together
    
    ind = [None]*4
    vmr_all = np.concatenate((vmr_n, vmr_e_fl, vmr_s_fl, vmr_w), axis = 1)
    ind[0] = range(len(lon))
    ind[1] = range(max(ind[0])+1,max(ind[0])+1 + len(lat))
    ind[2] = range(max(ind[1])+1,max(ind[1])+1 + len(lon))
    ind[3] = range(max(ind[2])+1,max(ind[2])+1 + len(lat))
    
    pca = np.linalg.svd(vmr_all)
    
    pcs = []
    for ii in range(numregions):
        u = pca[0][:,ii]
        u = np.reshape(u,(len(u),1))
        s = pca[1][ii]
        v = pca[2][ii,:]
        v = np.reshape(v,(1,len(v)))
        pcs.append(np.dot(u*s,v))
        
    # unravel pcs into directions
    pc_n = []
    pc_e = []
    pc_s = []
    pc_w = []
    
    for ii in range(numregions):
        pc_n.append(pcs[ii][:,ind[0]])
        pc_e.append(pcs[ii][:,ind[1]])
        pc_s.append(pcs[ii][:,ind[2]])
        pc_w.append(pcs[ii][:,ind[3]])

    # "un-reverse" pcs for east and south to match VMRs

    for ii in range(numregions):
        pc_e[ii] = np.fliplr(pc_e[ii])
        pc_s[ii] = np.fliplr(pc_s[ii])
        
        
    # create basis grids in format needed for inversion
        
    for ii in range(numregions):
        basis_grid_n[:,:,ii] = pc_n[ii]/vmr_n
        basis_grid_e[:,:,ii] = pc_e[ii]/vmr_e
        basis_grid_s[:,:,ii] = pc_s[ii]/vmr_s
        basis_grid_w[:,:,ii] = pc_w[ii]/vmr_w
        
    basis_grid_n = basis_grid_n[:,:,:,None]
    basis_grid_e = basis_grid_e[:,:,:,None]
    basis_grid_s = basis_grid_s[:,:,:,None]
    basis_grid_w = basis_grid_w[:,:,:,None]
    

    basis_ds = xray.Dataset({'bc_basis_n':(['height','lon','region','time'], basis_grid_n), \
                            'bc_basis_s':(['height','lon','region','time'], basis_grid_s), \
                            'bc_basis_e':(['height','lat','region','time'], basis_grid_e), \
                            'bc_basis_w':(['height','lat','region','time'], basis_grid_w)}, \
                coords = {'height': (['height'], heights), \
                          'lon':(['lon'], lon), \
                          'lat':(['lat'], lat), \
                          'region':(['region'],regions), \
                          'time':(['time'],[time])})
                
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + species + "_" + basis_case + "_" + domain + "_" + time.strftime('%m')+str(time.year) + ".nc")           
 
 
def cut_array2d(array, shape):
    """Needs a description for what it does
    """
    arr_shape = np.shape(array)
    xcut = np.linspace(0,arr_shape[0],shape[0]+1).astype(np.int)
    ycut = np.linspace(0,arr_shape[1],shape[1]+1).astype(np.int)
    xextent = [];    yextent = []
    for i in range(shape[0]):
        for j in range(shape[1]):
#            blocks.append(array[xcut[i]:xcut[i+1],ycut[j]:ycut[j+1]])
            xextent.append([xcut[i],xcut[i+1]])
            yextent.append([ycut[j],ycut[j+1]])
    return xextent,yextent