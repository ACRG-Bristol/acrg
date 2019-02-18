# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:56:36 2015

@author: as13988
"""
from __future__ import print_function
from __future__ import absolute_import


import numpy as np
import os
import datetime as dt
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib
import glob
import time
import csv 
from itertools import chain
import pdb
from . import acrg_read_GCwerks as read_GCwerks
from . import acrg_ICP

# Class that contains the USNs and associated DNos
class USNsDNos():
   def __init__(self): 
       USNs = ['H-237', 'E-097B', 'H-241', 'H-239', 'H-240']
       DNos = ['D091962', 'Unknown', 'D091983', 'DO91977', 'D091982']  
       
       self.USNs = USNs
       self.DNos = DNos

# Code to run for multiple files based on UAN or DNo
class Calcmulti:
   def __init__(self, CylinderNo, basedir = '/Users/as13988/Documents/Work/Cylinders/Calibrations/RGL/'):
       
       # Find matching USN or vice versa
        if (CylinderNo).find('D09') == -1:
            # We've been given a USN so find the DNo
            USN_tank = CylinderNo
            index = (USNsDNos().USNs).index(USN_tank)
            DNo_tank = (USNsDNos().DNos)[index]
            
        else:
            # We've been given a DNo so find the USN
            DNo_tank = CylinderNo
            index = (USNsDNos().DNos).index(DNo_tank)
            USN_tank = (USNsDNos().USNs)[index]            


        files = glob.glob(basedir+'*'+USN_tank+'*.dat')
        
        for i in np.arange(len(files)):
            print('Reading file: ' + files[i])
            data_i = read_data(files[i])
            
            acrg_ICP.PlotRawMM(data_i, outputdir='/Users/as13988/Documents/Work/Cylinders/Calibrations/RGL/Plots/')
            
            #pdb.set_trace()            
            
            means_i = acrg_ICP.Calcmeans(data_i)

            #pdb.set_trace()
            key_i = data_i.site + '_' + data_i.instrument
            
            if i == 0:
                data = {key_i:data_i}
                means = {key_i:means_i.means}
            else:
                data[key_i] = data_i
                means[key_i] = means_i.means
        
        # NB: syntax for extracting the same vaiable from each dictionary item
        # b = [means[i].co2mean for i in means.keys()]
        
        self.data = data
        self.means = means
        self.files = files


class read_data:
    def __init__(self, datafile):
        
        # Read in the data using the general code
        indata =  read_GCwerks.read_gcexport_crds(datafile)    
            
        # Find making USN or vice versa
        if (indata.filename.split('.')[2]).find('DO') == -1:
            # We've been given a USN so find the DNo
            USN_tank = (indata.filename.split('.')[2])
            index = (USNsDNos().USNs).index(USN_tank)
            DNo_tank = (USNsDNos().DNos)[index]
            
        else:
            # We've been given a DNo so find the USN
            DNo_tank = (indata.filename.split('.')[2])
            index = (USNsDNos().DNos).index(DNo_tank)
            USN_tank = (USNsDNos().USNs)[index]            
            
        self.site = indata.filename.split('.')[0]
        self.instrument = indata.filename.split('.')[1]
        self.processeddate = indata.filename.split('.')[3]
        self.USN = USN_tank
        self.DNo = DNo_tank
           
        self.date = indata.date
        self.time = indata.time
        self.datetime = indata.dt_date
        self.sampletype = indata.sampletype
        self.samplename = indata.samplename
        self.standard = indata.standard
        self.port = indata.port.astype('int')
        self.cavity_temp = indata.cavity_temp.astype('float')
        self.cavity_press = indata.cavity_press.astype('float')
        self.h2o = indata.h2o
        
        self.co2_orig = indata.co2_orig
        self.co2 = indata.co2
        self.co2sd_orig = indata.co2sd_orig
        self.co2sd = indata.co2sd
        self.co2flags = indata.co2flags
        self.co2_n = indata.co2_n
        self.co2_scale = 'WMOx2007'

        self.ch4_orig = indata.ch4_orig
        self.ch4 = indata.ch4
        self.ch4sd_orig = indata.ch4sd_orig
        self.ch4sd = indata.ch4sd
        self.ch4flags = indata.ch4flags
        self.ch4_n = indata.ch4_n
        self.ch4_scale = 'WMOx2014A'

        
        if indata.nogases == 3:
            self.co_orig = indata.co_orig
            self.co = indata.co
            self.cosd_orig = indata.cosd_orig
            self.cosd = indata.cosd
            self.coflags = indata.coflags
            self.co_n = indata.co_n
            self.co_scale = 'CSIRO'     
        
        self.nogases = indata.nogases        
                
        self.filename = indata.filename
        self.datadir = indata.datadir
