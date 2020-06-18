#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 14/05/2020

Test suite for acrg_hbmcmc module.

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_acrg_hbmcmc.py

@author: al18242
"""

import pytest

#import acrg_agage as agage
import acrg_hbmcmc.inversion_pymc3 as mcmc
import pymc3 as pm
import numpy as np
#from builtins import str

@pytest.fixture(scope="module")
def example_prior():
    ''' Define an example prior to be used throughout the test suite '''
    prior = {"pdf":"normal",
             "mu": 1.0,
             "sd": 0.1}
    
    return prior

    
def test_parsePrior(example_prior) :
    '''
    Test that the parse prior function works. May fail if a referenced function
    does not exist in the current pymc3 version
    '''
    with pm.Model() as model:
        test = mcmc.parsePrior("test", example_prior)
    assert test != None

def test_priorSetup():
    '''
    Primative test to ensure the priors are passing through the function call correctly
    '''
    Hx = np.array([[1],])
    Hbc = np.array([[1],])
    Y = np.array([2])
    error = np.array([10000])
    
    #use normal priors for easy comparison - could be expanded to other pdfs with a k-s test perhaps
    mus = [1, 2, 0.8]
    sds = [0.1, 0.3, 0.4]
    xprior = {"pdf":"normal",
              "mu":mus[0],
              "sd":sds[0]}
    bcprior = {"pdf":"normal",
              "mu":mus[1],
              "sd":sds[1]}
    sigprior = {"pdf":"normal",
              "mu":mus[2],
              "sd":sds[2]}
    
    xouts, bcouts, sigouts, convergence, step1, step2 = mcmc.inferpymc3(Hx, Hbc, Y, error, np.array([0]), np.array([0]), 
               xprior,bcprior, sigprior,nit=5e3, burn=1e3, tune=1e3, nchain=2)
    
    #test outputs match inputs
    #threshold of 0.02 chosen to be large enough to use a low nit for quicker testing
    assert ( np.abs(np.mean(xouts) - mus[0] ) < 0.02 ), "Discrepency in {}: {}".format("x mean", np.abs(np.mean(xouts) - mus[0] ))
    assert ( np.abs(np.std(xouts) - sds[0] ) < 0.02 ), "Discrepency in {}: {}".format("x std", np.abs(np.std(xouts) - sds[0] ))
    
    assert ( np.abs(np.mean(bcouts) - mus[1] ) < 0.02 ), "Discrepency in {}: {}".format("bc mean", np.abs(np.mean(bcouts) - mus[1] ))
    assert ( np.abs(np.std(bcouts) - sds[1] ) < 0.02 ), "Discrepency in {}: {}".format("bc std", np.abs(np.std(bcouts) - sds[1] ))
    
    assert ( np.abs(np.mean(sigouts) - mus[2] ) < 0.02 ), "Discrepency in {}: {}".format("sig mean", np.abs(np.mean(sigouts) - mus[2] ))
    assert ( np.abs(np.std(sigouts) - sds[2] ) < 0.02 ), "Discrepency in {}: {}".format("sig std", np.abs(np.std(sigouts) - sds[2] ))

def test_toyProblem():
    '''
    Primative test to ensure the priors are passing through the function call correctly
    '''
    n = 1000
    x_true = 1.4
    bc_true = 0.9
    sigma1_true = 1.0
    sigma2_true = 3.0
    indicators = np.zeros(n).astype(int)
    indicators_sigma = np.zeros(n).astype(int)
    indicators_sigma[(n//2):] = 1
    
    Hx = np.expand_dims( (1+0.5*np.sin(np.arange(n)*2*np.pi/n ))*10., axis=0)
    Hbc = np.ones((1,n))*1000.
    Y = Hx * x_true + Hbc * bc_true
    
    Y[0, :(n//2)] += np.random.normal(0.0, sigma1_true, n//2)
    Y[0, (n//2):] += np.random.normal(0.0, sigma2_true, n//2)
    error = np.ones(n)*0.
    
    #use normal priors for easy comparison - could be expanded to other pdfs with a k-s test perhaps
    xprior = {"pdf":"lognormal",
              "mu":0.0,
              "sd":0.4}
    bcprior = {"pdf":"normal",
              "mu":1.0,
              "sd":0.1}
    sigprior = {"pdf":"halfnormal",
              "sigma":1.0}
                
    
    xouts, bcouts, sigouts, convergence, step1, step2 = mcmc.inferpymc3(Hx, Hbc, Y, error, indicators, indicators_sigma,
               xprior,bcprior, sigprior,nit=5e3, burn=1e3, tune=1e3, nchain=2)
    
    #test outputs match inputs
    #threshold chosen to be large enough to use a low nit for quicker testing
    assert ( np.abs(np.mean(xouts) - x_true ) < 0.1 ), "Discrepency in {}: {}".format("x mean", np.abs(np.mean(xouts) - x_true ))
    
    assert ( np.abs(np.mean(bcouts) - bc_true ) < 0.02 ), "Discrepency in {}: {}".format("bc mean", np.abs(np.mean(bcouts) - bc_true ))
    
    assert ( np.abs(np.mean(sigouts[:,0,0]) - sigma1_true ) < 0.1 ), "Discrepency in {}: {}".format("sig1 mean", np.abs(np.mean(sigouts[:,0,0]) - sigma1_true ))
    
    assert ( np.abs(np.mean(sigouts[:,0,1]) - sigma2_true ) < 0.1 ), "Discrepency in {}: {}".format("sig2 mean", np.abs(np.mean(sigouts[:,0,1]) - sigma2_true ))
    