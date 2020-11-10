#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 09:38:39 2020

@author: rt17603
"""

import numpy as np

def semi_corr(std,correlation_coeff=0.5):
    
    nwindow = len(std)
    
    A = np.ones((nwindow))
    cov = np.zeros((nwindow,nwindow))
    for ii in range(nwindow):
        for bb in range(nwindow):
            if ii != bb:
                cov[ii,bb] = correlation_coeff*std[ii]*std[bb]
            else:
                cov[ii,ii] = std[ii]*std[ii]
    smoothed_uncertainty = np.sqrt(np.matmul(np.matmul(A,cov),A.T)/(nwindow*nwindow))

    return smoothed_uncertainty
