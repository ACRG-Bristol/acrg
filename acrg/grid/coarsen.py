# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:09:07 2015

@author: chxmr
"""
import numpy as np
from numba import jit

@jit(nopython = True)
def coarsen(arrayFine, latFine, lonFine, factor = 2, mean = True):
    '''Coarsen a fine grid by some integer factor.
    E.g. if the input array is 1000 x 500, and is coarsened by a factor of 2
    the output array will be 500 x 250. By default, the output is the mean
    of the inputs.
    
    Args:
        arrayFine (array): 
            array of values for fine resolution grid
        latFine (array): 
            1D latitudes of fine resolution grid
        lonFine (array): 
            1D longitudes of fine resolution grid
        factor (int, optional): 
            Factor by which to coarsen array. Default is 2
        mean (bool): 
            True/False. 
            True calculates the mean of the fine array values; False calculates the sum.
            Default is True
            
    Returns:
        out (array): 
            grid of coarsened values
        outLat (array): 
            Latitudes of coarsened grid
        outLon (array): 
            Longitudes of coarsened grid

    Example:
        out, outLat, outLon = coarsen(fineArray, fineLat, fineLon, factor = 10, mean = True)
    
    '''
    
    out = np.zeros((len(latFine)//factor, len(lonFine)//factor))
    latOut = np.zeros(len(latFine)//factor)
    lonOut = np.zeros(len(lonFine)//factor)

    for lati in range(len(latFine)//factor):
        for loni in range(len(lonFine)//factor):
            count = 0
            for f in range(factor):
                out[lati, loni] += arrayFine[lati*factor + f,
                                             loni*factor + f]
                count += 1
            if mean:
                out[lati, loni] = out[lati, loni]/count

    # Output grid
    for lati in range(len(latFine)//factor):
        for f in range(factor):
            latOut[lati] += latFine[lati*factor + f]/float(factor)
    for loni in range(len(lonFine)//factor):
        for f in range(factor):
            lonOut[loni] += lonFine[loni*factor + f]/float(factor)

    return out, latOut, lonOut