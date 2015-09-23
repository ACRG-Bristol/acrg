# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 15:46:31 2015

@author: ag12733

Creat basis functions for SOUTHASIA domain
"""

import numpy as np
import xray
import glob
import pandas as pd
import acrg_name as name

# read in fields file to get domain
fields_file_path = "/shared_data/air/shared/NAME/fp/"
basis_dir = "/data/shared/NAME/basis_functions/"
bc_basis_dir = "/data/shared/NAME/bc_basis_functions/"


def basis_blocks(domain = "SOUTHASIA", time = "2012-01-01", blocksize = 25, basis_case = "25x25"): 

# creates basis functions in square blocks (e.g., 5x5 grid cells per region)    
    
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

def basis_transd(domain = "SOUTHASIA", time = "2012-01-01", basis_case = "transd"): 

# creates basis blocks for the "transdimensional inversion"
# creates central block for variable basis and four fixed regions surrounding the central box   
    
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
    
    sub_lon_min = 65
    sub_lon_max = 100
    sub_lat_min = 8
    sub_lat_max = 35
    
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

def basis_bc_blocks(domain = "SOUTHASIA", basis_case = "NESW", time = "2012-01-01", vertical=4):
    
# creates uniform blocks for each direction (NESW). 
# Each direction is split into vertical slabs governed by input variable
# currently splits vertical grid

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
    
    vertical_index = np.round(np.linspace(0,len(heights), vertical+1))
    
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
                          'region':(['region'],regions), \
                          'time':(['time'],[time])})
                
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc")
                
def basis_bc_pca(domain = "SOUTHASIA", time = "2012-01-15", species='ch4', units='ppb', basis_case='pca'):
    
    files = glob.glob(fields_file_path + domain + "/*")
    
    time = pd.to_datetime(time)
    
    with xray.open_dataset(files[0]) as temp:
        fields_ds = temp.load()

    if len(files) == 0:
        print("Can't find Fields files: " + domain)
        return None
    
    lat = fields_ds["lat"].values
    lon = fields_ds["lon"].values
    
    regions = range(3)
    heights = range(500,20500,1000)
    
    basis_grid_n = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_e = np.zeros((len(heights),len(lat),len(regions)))
    basis_grid_s = np.zeros((len(heights),len(lon),len(regions)))
    basis_grid_w = np.zeros((len(heights),len(lat),len(regions)))

    # find principal companents of MOZART vmr

    bc_ds =  name.boundary_conditions(domain, species)
    bc_ds = bc_ds.sel(time = time)
    
    if units is 'ppb':
        bc_ds = bc_ds*1e9
        
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
    for ii in range(3):
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
    
    for ii in range(3):
        pc_n.append(pcs[ii][:,ind[0]])
        pc_e.append(pcs[ii][:,ind[1]])
        pc_s.append(pcs[ii][:,ind[2]])
        pc_w.append(pcs[ii][:,ind[3]])

# reverse pcs for east and south to make continuous curtain
    for ii in range(3):
        pc_e[ii] = np.fliplr(pc_e[ii])
        pc_s[ii] = np.fliplr(pc_s[ii])
        
        
# create basis grids in format needed for inversion
        
    for ii in range(3):
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
                
    basis_ds.to_netcdf(bc_basis_dir + domain +"/" + basis_case + "_" + domain + "_" + str(time.year) + ".nc")           
    
def cut_array2d(array, shape):
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