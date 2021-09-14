#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 16:02:01 2018

@author: rt17603
"""
from __future__ import print_function
from __future__ import division

# Barometric formula

import numpy as np

#Up to 11,000m
h_b = 0         # m
P_b = 101325.00 # Pa
T_b = 288.15    # K
L_b = -0.0065   # K/m
R = 8.3144598   # J/mol/K
g = 9.80665     # m/s2
M = 0.0289644   # kg/mol

def pressure_at_height(h):
    P = P_b*(T_b/(T_b + L_b*(h-h_b)))**((g*M)/(R*L_b))    
    return P

def height_at_pressure(P):
    h = T_b*((P_b/P)**((R*L_b)/(g*M)) - 1)/L_b + h_b
    return h

if __name__ == "__main__":

    print('-----------------')
    print('Checking pressure_at_height() function')
    height = np.arange(0,1000,100)
    pressure = np.zeros(len(height))
    for i,h in enumerate(height):
        P = pressure_at_height(h)
        print('Pressure at height {0} m, {1} Pa'.format(h,P))
        pressure[i] = P
    
    print('-----------------')
    print('Checking height_at_pressure() function')
    for p in pressure:
        print('Height at pressure {0} Pa, {1} m'.format(p,height_at_pressure(p)))
    
    percentage = np.arange(1.0,11.0)
    print('-----------------')
    print('Ground pressure, {0} Pa, {1} m'.format(P_b,height_at_pressure(P_b)))
    for p in percentage:
        Pb_1 = P_b - P_b*p/100.
        print('{0}% above ground pressure, {1} Pa, {2} m'.format(p,Pb_1,height_at_pressure(Pb_1)))

    
