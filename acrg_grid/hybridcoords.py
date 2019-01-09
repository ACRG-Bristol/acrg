# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 10:56:48 2014
; :Purpose:
;   Calculate mid-point hybrid-coordinate pressure fields.
;   For the ith level:
;   
;     P_i=A_i*P_0 + B_i*P_surface
;   
;   If /half is set, it is assumed that A (and B) define the coordinate system at
;   grid cell boundaries. Mid-point pressure values will the be output, using:
;     
;     P_i=(P_{i-1/2} + P_{i+1/2})/2
;
; :Inputs:
;   A: (input, floating, required), hybrid A values. If P0 is set, A is dimensionless.
;     If A P0 is not set, A must be in Pa.
;   B: (input, floating, optional), hybrid B values. Dimensionless. 
;     If left out, surface pressure independent (fixed) pressure levels are output. 
    PS: This is assumed to be of dimensions lon by lat. If you've got a lon by lat by 
    time array then you'll need to make a loop to go through each time point.
    
; :Keywords:
;   P0: (input, optional, floating) Reference pressure to multipy the A coordinate
;   Half: (input, optional, boolean) Set to True if A and B define coordinates at grid cell
;     boundaries. The output will still be at grid cell centres.
;     Note that if true, the output will have one fewer vertical elements than A and B.
;
; :Outputs:
;   3D pressure field (lon, lat, level)
;
; :Example::
; 
;   IDL> Pressure=mr_hybrid_coords(A, PS, B=B, P0=P0, half=half)
;
; :History:
; 	Written by: Matt Rigby, University of Bristol, Jul 16, 2013
    Hacked into Python by Ann Stavert 19th Nov 2014

@author: as13988
"""

import numpy as np

# Set defaults
P0=1.e5
half=0
B = 0

def hybridcoords(A, PS, B=B, P0=P0, half=half):

  #Define dimensions and output array  
  LevSize = len(A)
  LonSize = np.shape(PS)[1]
  LatSize = np.shape(PS)[0]
  #P = np.empty((LonSize, LatSize, LevSize))
  P = np.empty((LevSize, LatSize, LonSize))
  P[:] = 0

  #Check default inputs
  if len(B) > 1:
    if len(B) != LevSize:
        print("A and B must have same dimensions")
    has_B=1
  else: 
    has_B=0
  
  if len(np.shape(PS)) > 2:
      print("PS variable needs to be 2D i.e. lat by lon")
  
  #if np.shape(PS) != np.shape(P)[0:2]:
    #  PS = np.transpose(PS)
  
  #Calculate pressure
  for LevI in np.arange(LevSize):
    if has_B==0:    
        P[LevI, :, :] = A[LevI]*P0
   
    if has_B==1:
        P[LevI, :, :] =  A[LevI]*P0 + B[LevI]*PS
        


  #If /half, calculate mid-point pressure levels
  if half==1:
    P_out=np.empty((LevSize-1, LatSize, LonSize))
    P_out[:] = 0
    
    for LevI in np.arange(LevSize-1):
        P_out[LevI, :, :]=0.5*(P[LevI, :, :] + P[LevI+1, :, :])
        
    P = P_out
  

  return P
