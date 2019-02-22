# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 14:32:15 2015

@author: aw12579
"""
from __future__ import division

from past.utils import old_div
import pandas as pd
import glob
import matplotlib.pyplot as plt   
import numpy as np
import wea_FAAM as fa 
#import wea_Bag
import datetime
import acrg_time
import acrg_agage as agage
import acrg_name_xray as name
import wea_background_insitu
import xray





def read14Cdata(site,start,end):
    if site == 'TAC':
        data= agage.get_obs([site], "co2c14", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end,height = '185m')
    else:
        data= agage.get_obs([site], "co2c14", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end)
    return data

def read13CCO2data(site,start,end):
    if site == 'TAC':
        data= agage.get_obs([site], "co2c13", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end,height = '185m')
    else:
        data= agage.get_obs([site], "co2c13", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end)
    return data

def calculate14Cenhancment(site,start,end):
    huba =read14Cdata(site,start,end)
    #data = {site: huba[site]}
    huba['.units']=1.
    fp=name.footprints_data_merge(huba, domain = "EUROPE", species = "co2c14", calc_bc = False, calc_timeseries = True) # cant calucluate boundery contitions yet, throws error message 
    return fp  
    
def readCO2data(site,start,end):
    if site == 'TAC':
        data= agage.get_obs([site], "co2", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end,height = '185m')
    else:
        data= agage.get_obs([site], "co2", network = ["GAUGE"], instrument = ["GAUGE-GLASS"] , start = start, end = end)
    return data    
    
def caclulateRH(site,start,end,permil):
    #if site == 'TAC':
    #    fp =name.footprints(site, start = start, end = end, domain='EUROPE', height ='185m', species = None) #takes long
    #else:
    #    fp =name.footprints(site, start = start, end = end, domain='EUROPE', species = None) #takes long  
    huba =read14Cdata(site,start,end)
    huba['.units']=1.
    fp=name.footprints_data_merge(huba, domain = "EUROPE", species = "co2-rh", calc_bc = False, calc_timeseries = True) # cant calucluate boundery contitions yet, throws error message 
    #flux_RH =name.flux('EUROPE', 'co2-rh')
    #comb=name.combine_datasets(flux_RH,fp)
    #ts=name.timeseries(comb)
    
    ts=fp[site].mf_mod
    rh_ts=C14_bigdelta_to_molpermol(permil,ts,-19)
    
    return rh_ts
    
def C14_molpermol_to_bigdelta(molpermol_C14,molpermol_C12=0.0004, delta_C13=-8.8):
    Rref=0.000000000001176
    Delta14=(old_div((old_div(molpermol_C14,molpermol_C12))*(1-old_div(2*(25+delta_C13),1000)),Rref)-1.0)*1000.0
    return Delta14
    
def C14_bigdelta_to_molpermol(Delta_C14,molpermol_C12=0.0004, delta_C13=-8.8):
    Rref=0.000000000001176     
    molpermol_C14=old_div((Delta_C14/1000.0+1.0)*Rref,(1.0-(2.0*(25.0+delta_C13))/1000.0)*molpermol_C12) 
    return molpermol_C14
    
def C13_molpermol_to_smalldelta(molpermol_C13,molpermol_C12=0.0004):
    delta13= (old_div(molpermol_C13,(molpermol_C12)/0.0112372)-1.0)*1000.0
    return delta13
    
def C13_smalldelta_to_molpermol(delta_C13,molpermol_C12=0.0004):
    molpermol_C13= (delta_C13/1000.0+1.0)*molpermol_C12*0.0112372
    return molpermol_C13


