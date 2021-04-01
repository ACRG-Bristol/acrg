#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:53:05 2019

Tests for GCWerks python code: fitting.py
This function can generate data, fits a linear, exponential or logarithmic function to the data, and can plit the data
The test deterimines if the fitted coefficients are the same as a benchmark case.

@author: mi19881
"""

import pytest
import acrg_GCWerks.fitting as fitting
import acrg_GCWerks.CRDS_H2OCorr as H2OCorr
import numpy as np
import os
from os.path import join
import sys
from acrg_config.paths import paths


acrg_path = paths.acrg


@pytest.mark.parametrize("x, y, fit_type, expected", [
        (np.array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99.]), 
         np.array([-10.53584414, -10.56294572, -10.04705205, -10.03561186, -9.50719678,  -8.79064983,  -8.48946774,  -8.20812915, -7.87923221,  -7.97544694,  -7.65549334,  -7.35769347, -7.41093951,  -7.05055927,  -7.10125236,  -7.37607744, -7.00496802,  -7.01174075,  -6.5623265 ,  -6.56100477, -6.6929464 ,  -6.35009707,  -6.39093434,  -6.17854788, -6.45397669,  -6.13881583,  -6.3698625 ,  -6.60000548, -6.55779822,  -6.58201325,  -6.17987073,  -6.30485273, -6.02077888,  -6.36454151,  -6.25095878,  -6.00363181, -6.31997816,  -6.07186504,  -6.07141494,  -6.2508431 , -6.1765406 ,  -5.84324917,  -6.07672551,  -6.05177869, -6.13705916,  -6.10127939,  -6.32063014,  -5.81497993, -5.81984493,  -6.21519859,  -5.7490415 ,  -6.02698534, -6.08339257,  -5.99959782,  -5.90083872,  -6.12700246, -6.21205265,  -6.28256083,  -6.13619696,  -5.97154102, -6.27831165,  -5.97006408,  -5.79368017,  -6.10657687, -6.25091804,  -6.09070201,  -6.14210467,  -6.06106649, -6.17402278,  -6.19737545,  -6.01034219,  -5.97818203, -5.71245051,  -6.25981776,  -5.91624205,  -6.27668809, -5.94180587,  -5.84442549,  -6.28564834,  -5.91734788, -5.79551107,  -6.12705831,  -5.99805696,  -6.24795536, -5.741144  ,  -5.94997759,  -5.87541831,  -6.0874021 , -6.16782691,  -6.15318589,  -5.84082495,  -6.17147502, -6.2523395 ,  -6.01514185,  -5.92535965,  -5.7241558 , -6.22794771,  -6.13828348,  -5.76138219,  -5.97006163]),
         'lin',
         np.round(np.array([ 0.0236633 , -7.68835482]),7)),
         (np.array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99.]), 
         np.array([-10.53584414, -10.56294572, -10.04705205, -10.03561186, -9.50719678,  -8.79064983,  -8.48946774,  -8.20812915, -7.87923221,  -7.97544694,  -7.65549334,  -7.35769347, -7.41093951,  -7.05055927,  -7.10125236,  -7.37607744, -7.00496802,  -7.01174075,  -6.5623265 ,  -6.56100477, -6.6929464 ,  -6.35009707,  -6.39093434,  -6.17854788, -6.45397669,  -6.13881583,  -6.3698625 ,  -6.60000548, -6.55779822,  -6.58201325,  -6.17987073,  -6.30485273, -6.02077888,  -6.36454151,  -6.25095878,  -6.00363181, -6.31997816,  -6.07186504,  -6.07141494,  -6.2508431 , -6.1765406 ,  -5.84324917,  -6.07672551,  -6.05177869, -6.13705916,  -6.10127939,  -6.32063014,  -5.81497993, -5.81984493,  -6.21519859,  -5.7490415 ,  -6.02698534, -6.08339257,  -5.99959782,  -5.90083872,  -6.12700246, -6.21205265,  -6.28256083,  -6.13619696,  -5.97154102, -6.27831165,  -5.97006408,  -5.79368017,  -6.10657687, -6.25091804,  -6.09070201,  -6.14210467,  -6.06106649, -6.17402278,  -6.19737545,  -6.01034219,  -5.97818203, -5.71245051,  -6.25981776,  -5.91624205,  -6.27668809, -5.94180587,  -5.84442549,  -6.28564834,  -5.91734788, -5.79551107,  -6.12705831,  -5.99805696,  -6.24795536, -5.741144  ,  -5.94997759,  -5.87541831,  -6.0874021 , -6.16782691,  -6.15318589,  -5.84082495,  -6.17147502, -6.2523395 ,  -6.01514185,  -5.92535965,  -5.7241558 , -6.22794771,  -6.13828348,  -5.76138219,  -5.97006163]),
         'exp',
         np.round(np.array([ -4.8862878, 0.10717707,-6.03625635]),7))])
      
def test_acrg_GCWerks_fitting(x, y, fit_type, expected):
    '''
    Test if the coefficients are the same as those calculated for a benchmark case
    '''
    fit_coef = fitting.fit_data(x, y, fit_type = fit_type)
    assert np.all(np.round(fit_coef.coeffs,7) == expected)

def test_acrg_GCWerks_CRDS_H2OCorr():
    '''
    Test CRDS_H2OCorr.py function is creating the same corrections
    '''
    basedir = join(acrg_path,'tests/files/raw_obs/H2OCorr/')
    outdir = join(acrg_path,'tests/files/raw_obs/H2OCorr/Plots/')
    out1, out2 = H2OCorr.calcread_multi(sites = "BSD", years = "2018", basedir = basedir, outdir = outdir)
    
    assert round(out2['B_0']['a_ch4'][0],6) == -0.013575
    assert round(out2['B_0']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['B_0']['a_co2'][0],6) == -0.015701
    assert round(out2['B_0']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['B_0']['b_ch4'][0],6) == 0.000386
    assert round(out2['B_0']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['B_0']['b_co2'][0],6) == 0.000171
    assert round(out2['B_0']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['B_0']['date'] == '20180814'
    assert out2['B_0']['site'] == 'BSD'
    
    assert round(out2['B_1']['a_ch4'][0],6) == -0.013575
    assert round(out2['B_1']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['B_1']['a_co2'][0],6) == -0.015701
    assert round(out2['B_1']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['B_1']['b_ch4'][0],6) == 0.000386
    assert round(out2['B_1']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['B_1']['b_co2'][0],6) == 0.000171
    assert round(out2['B_1']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['B_1']['date'] == '20180814'
    assert out2['B_1']['site'] == 'BSD'
    
    assert round(out2['B_2']['a_ch4'][0],6) == -0.013575
    assert round(out2['B_2']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['B_2']['a_co2'][0],6) == -0.015701
    assert round(out2['B_2']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['B_2']['b_ch4'][0],6) == 0.000386
    assert round(out2['B_2']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['B_2']['b_co2'][0],6) == 0.000171
    assert round(out2['B_2']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['B_2']['date'] == '20180814'
    assert out2['B_2']['site'] == 'BSD'
    
    assert round(out2['B_8']['a_ch4'][0],6) == -0.013575
    assert round(out2['B_8']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['B_8']['a_co2'][0],6) == -0.015701
    assert round(out2['B_8']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['B_8']['b_ch4'][0],6) == 0.000386
    assert round(out2['B_8']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['B_8']['b_co2'][0],6) == 0.000171
    assert round(out2['B_8']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['B_8']['date'] == '20180814'
    assert out2['B_8']['site'] == 'BSD'
    
    assert round(out2['D_0']['a_ch4'][0],6) == -0.013575
    assert round(out2['D_0']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['D_0']['a_co2'][0],6) == -0.015701
    assert round(out2['D_0']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['D_0']['b_ch4'][0],6) == 0.000386
    assert round(out2['D_0']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['D_0']['b_co2'][0],6) == 0.000171
    assert round(out2['D_0']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['D_0']['date'] == '20180814'
    assert out2['D_0']['site'] == 'BSD'
    
    assert round(out2['D_1']['a_ch4'][0],6) == -0.013575
    assert round(out2['D_1']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['D_1']['a_co2'][0],6) == -0.015701
    assert round(out2['D_1']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['D_1']['b_ch4'][0],6) == 0.000386
    assert round(out2['D_1']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['D_1']['b_co2'][0],6) == 0.000171
    assert round(out2['D_1']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['D_1']['date'] == '20180814'
    assert out2['D_1']['site'] == 'BSD'
    
    assert round(out2['D_2']['a_ch4'][0],6) == -0.013575
    assert round(out2['D_2']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['D_2']['a_co2'][0],6) == -0.015701
    assert round(out2['D_2']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['D_2']['b_ch4'][0],6) == 0.000386
    assert round(out2['D_2']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['D_2']['b_co2'][0],6) == 0.000171
    assert round(out2['D_2']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['D_2']['date'] == '20180814'
    assert out2['D_2']['site'] == 'BSD'
    
    assert round(out2['D_8']['a_ch4'][0],6) == -0.013575
    assert round(out2['D_8']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['D_8']['a_co2'][0],6) == -0.015701
    assert round(out2['D_8']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['D_8']['b_ch4'][0],6) == 0.000386
    assert round(out2['D_8']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['D_8']['b_co2'][0],6) == 0.000171
    assert round(out2['D_8']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['D_8']['date'] == '20180814'
    assert out2['D_8']['site'] == 'BSD'
    
    assert round(out2['S_0']['a_ch4'][0],6) == -0.013575
    assert round(out2['S_0']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['S_0']['a_co2'][0],6) == -0.015701
    assert round(out2['S_0']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['S_0']['b_ch4'][0],6) == 0.000386
    assert round(out2['S_0']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['S_0']['b_co2'][0],6) == 0.000171
    assert round(out2['S_0']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['S_0']['date'] == '20180814'
    assert out2['S_0']['site'] == 'BSD'
    
    assert round(out2['S_1']['a_ch4'][0],6) == -0.013575
    assert round(out2['S_1']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['S_1']['a_co2'][0],6) == -0.015701
    assert round(out2['S_1']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['S_1']['b_ch4'][0],6) == 0.000386
    assert round(out2['S_1']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['S_1']['b_co2'][0],6) == 0.000171
    assert round(out2['S_1']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['S_1']['date'] == '20180814'
    assert out2['S_1']['site'] == 'BSD'
    
    assert round(out2['S_2']['a_ch4'][0],6) == -0.013575
    assert round(out2['S_2']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['S_2']['a_co2'][0],6) == -0.015701
    assert round(out2['S_2']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['S_2']['b_ch4'][0],6) == 0.000386
    assert round(out2['S_2']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['S_2']['b_co2'][0],6) == 0.000171
    assert round(out2['S_2']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['S_2']['date'] == '20180814'
    assert out2['S_2']['site'] == 'BSD'
    
    assert round(out2['S_8']['a_ch4'][0],6) == -0.013575
    assert round(out2['S_8']['a_ch4_sd'][0],6) == 0.000391
    assert round(out2['S_8']['a_co2'][0],6) == -0.015701
    assert round(out2['S_8']['a_co2_sd'][0],6) == 0.000140
    assert round(out2['S_8']['b_ch4'][0],6) == 0.000386
    assert round(out2['S_8']['b_ch4_sd'][0],6) == 0.000245
    assert round(out2['S_8']['b_co2'][0],6) == 0.000171
    assert round(out2['S_8']['b_co2_sd'][0],6) == 8.9e-05
    assert out2['S_8']['date'] == '20180814'
    assert out2['S_8']['site'] == 'BSD'       
    
    
    
