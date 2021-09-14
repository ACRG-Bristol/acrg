# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 11:23:48 2014
@author: chxmr
"""
from builtins import range
import numpy as np

def areagrid(lat, lon):
  """Calculates grid of areas (m2) given arrays of latitudes and longitudes

  Args:
      lat (array): 
          1D array of latitudes
      lon (array): 
          1D array of longitudes
        
  Returns:
      area (array): 
          2D array of areas of of size lat x lon
      
  Example:
    import acrg_grid
    lat=np.arange(50., 60., 1.)
    lon=np.arange(0., 10., 1.)
    area=acrg_grid.areagrid(lat, lon)
    
  """
        
  re=6367500.0	#radius of Earth in m
  
  dlon=abs(np.mean(lon[1:] - lon[0:-1]))*np.pi/180.
  dlat=abs(np.mean(lat[1:] - lat[0:-1]))*np.pi/180.
  theta=np.pi*(90.-lat)/180.
  
  area=np.zeros((len(lat), len(lon)))
  
  for latI in range(len(lat)):
    if theta[latI] == 0. or np.isclose(theta[latI], np.pi):
      area[latI, :]=(re**2)*abs(np.cos(dlat/2.)-np.cos(0.))*dlon
    else:
      lat1=theta[latI] - dlat/2.
      lat2=theta[latI] + dlat/2.
      area[latI, :]=((re**2)*(np.cos(lat1)-np.cos(lat2))*dlon)

  return area
