# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 16:03:08 2014
 
FUNCTIONS
___________ 
 
 
distance


Calculates the distance in km between two lat/lon pairs using the haversine formula.


Example:

distance = acrg_grid.haversine.distance(origin, destination, radius=radius)


Inputs:

origin: [lat_1, lon_1] two element array containing the lat and lon of your original point
destination: [lat_2, lon_2] two element array containing the lat and lon of your destination point
radius: radius in km. Defaults to 6371km.


Output:

distance between the two points in km



CLASSES
_________


multipledistances

Calculates the distance in km between a single lat/lon pair and a grid of lat/lons using the haversine formula.
Also identifies the minimum distance and it's location.


Example:

distances = acrg_grid.haversine.multipledistances(origin, lat, lon, radius=radius)

Inputs:

origin: [lat_1, lon_1] two element array containing the lat and lon of your original point
lat: an array of latitudes of length n
lon: an array of longitudes of lenth m
radius: radius in km. Defaults to 6371km.

Outputs:

distances= an n by m array containing the distances from every point in the nxm grid to the origin
mindist = the minimum distance
mindist_index = two element array containing the index in n and m corresponding to the minium distance
e.g. if the minimum distance is at the 2nd latitude and 17th longitude then mindist_index = [2,17]
mindist_loc = two element array containing the lat and lon corresponding to the minimum distance

@author: as1398
Based on code I nicked off a website made by Wayne Dick


28/11/2014
MLR: Added distancelist. Given a list of lon/lat pairs, outputs a list of 
    distances from the "origin"
    
    e.g. distances=distancelist((origin_lat, origin_lon), zip(lats, lons))

    if lons and lats are of length N, output will also be of length N

"""
 
import math
import numpy as np
from numba import jit
import pdb

radius = 6371 #km
 
def distance(origin, destination, radius=radius):
    lat1, lon1 = origin
    lat2, lon2 = destination

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
 
    return d
    
# Turns the multipledistances into a function rather than a class
@jit(nopython=True)
def fn_multipledistances(origin, lat, lon, distances, index):
    # Calculate the distance between the point of interest and 
    #  ALL points in the grid
    lat1 = origin[0]
    lon1 = origin[1]
    
    #distances = np.zeros(len(lat))    
    
    
    for j in index:
        
        lat2 = lat[j]
        lon2 = lon[j]

        dlat = np.radians(lat2-lat1)
        dlon = np.radians(lon2-lon1)
        
        a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) \
            * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        
        d_j = radius * c        
        
        distances[j] = d_j
        
       
    return distances
    
    
class multipledistances:
    def __init__(self, origin, lat, lon, \
        radius=radius):
          
        # Calculate the distance between the point of interest and 
        #  ALL points in the grid
        for j in np.arange(len(lon)):
            d = np.asarray([distance(origin, [lat[i],lon[j]]) \
                for i in np.arange(len(lat))])
        
            if j==0:
                distances = d
            
            else:
                distances = np.concatenate((distances,d))


        # Find the minimum distance
        mindist = min(distances)
        
        # Extract the lat and lon that correspond to the minimum distance
        min_index = np.where(distances == mindist)
        
        lat_index = min_index[0] -((min_index[0]/len(lat))*len(lat)) 
        lon_index = (min_index[0]/len(lat))
        
        mindist_index=[lat_index, lon_index]        
        
        mindist_loc=[lat[lat_index], lon[lon_index]]        
        
        
        self.distances=distances
        self.mindist = mindist
        self.mindist_index = mindist_index
        self.mindist_loc = mindist_loc


def distancelist(origin, destinations, \
        radius=radius):

    distances=np.asarray([distance(origin, destination, radius=radius) \
        for destination in destinations])
    
    return distances
    