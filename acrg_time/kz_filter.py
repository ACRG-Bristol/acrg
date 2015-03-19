# -*- coding: utf-8 -*-
"""
; :Purpose:
;   Filter data using Kolmogorov-Zurbenko filter
;   
; :Inputs:
;   x_in: array of locations of the data to be filtered (e.g. time)
;   y_in: array to be filtered
;   sigma: array of uncertainties in y_in
;   
; :Keywords:
;   iterations: (input, integer) numeber of filter passes
;   window: (input, floating) width of filter window in the same units as x_in
;   growth: (input, Boolean) set to true to filter the growth rate of y_in
;   min_elements: (input, integer), minimum number of data points that must lie within a particular 
;     window in order for average to be calculated


Example:

    from acrg_time import kz_filter
    x_smoothed, y_smoothed = kz_filter(x, y, window=12., iterations=4)


Created on Thu Nov 20 12:27:48 2014

@author: chxmr
"""

import numpy as np
import matplotlib.dates as dates

def kz_filter(x_in, y_in, \
    growth=False, iterations=4, window=1., min_elements=1):

    x=x_in[:]
    y=np.copy(y_in)

    x=dates.date2num(x)

    if growth is True:
        y=(y[1:-1] - y[0:-2])/(x[1:-1] - x[0:-2])
        x=(x[0:-2] + x[1:-1])/2.
    else:
        y=y
        x=x
    
    ys=y
    xStart=x[0] + window
    xEnd=x[-1] - window
    xiStart=np.where(x >= xStart)[0][0]
    xiEnd=np.where(x <= xEnd)[0][-1]
    
    xs=x
    
    for k in range(iterations):
        ys_prev=ys
        xs_prev=xs
        for xi in range(len(xs)):
            wh=np.where( (xs_prev > xs_prev[xi]-window/2.) * \
              (xs_prev <= xs_prev[xi]+window/2.) * \
              np.isfinite(ys_prev))
            if len(wh[0]):
                ys[xi]=np.mean(ys_prev[wh])
            else:
                ys[xi]=float("Nan")
    
    ys[0:xiStart]=float("Nan")
    ys[xiEnd:]=float("Nan")
    
    wh=np.where(np.isfinite(ys))
    ys=ys[wh]
    xs=xs[wh]
    
    out_window=window*(float(iterations)**(0.5))
    
    xs=dates.num2date(xs)
    
    return xs, ys
