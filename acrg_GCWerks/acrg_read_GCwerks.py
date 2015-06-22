# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:22:22 2015

@author: as13988
"""

import os
import numpy as np
import datetime as dt
import pdb

 # Class to read in the txt output of gcexport made using the standard peak.list and report config on Dagage2
def read_gcexport_crds(datafile):
        
        print 'Reading:'
        print datafile
        if type(datafile) == tuple:
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
        elif type(datafile) == str:
            data=np.genfromtxt(datafile, dtype=str, skip_header=3)
        
               
        
        if np.shape(data)[1] > 18 :
            data = read_gcexport_crds_new(data, datafile)
        else:
            data = read_gcexport_crds_old(data, datafile)
            
            
        return data

#:Created: 14 Jan 15 15:28 GMT
#     -      -         -            -            -    -    cavity    cavity         -       co2     co2       ch4     ch4        co      co 
#  date   time      type       sample     standard port      temp     press       h2o         C   stdev         C   stdev         C   stdev 
# 140710 093530      tank  USN20132963        H-239    6    45.000   139.999     0.532    418.42*  0.188*  1993.75*  0.204*   200.84*  5.913*

class read_gcexport_crds_old:
    def __init__(self, data, datafile):
        
        date = data[:,0]
        time = data[:,1]
        sampletype = data[:,2]
        samplename = data[:,3]
        standard = data[:,4]
        port = data[:,5]
        cavity_temp = data[:,6]
        cavity_press = data[:,7] 
        h2o_orig = data[:,8]
        co2_orig = data[:,9]
        co2sd_orig = data[:,10]
        co2_n = data[:,11]
        ch4_orig = data[:,12]
        ch4sd_orig = data[:,13]
        ch4_n = data[:,14]
        
        dirname, filename = os.path.split(datafile)

        
        co2flags = Makeflags(co2_orig)
        ch4flags = Makeflags(ch4_orig)
        
        
        # Strip flags from conc data
        co2_stripped = Stripflags(co2_orig)
        co2sd_stripped = Stripflags(co2sd_orig)
        ch4_stripped = Stripflags(ch4_orig)
        ch4sd_stripped = Stripflags(ch4sd_orig)
        h2o_stripped = Stripflags(h2o_orig)
        
        # CO data is only at some sites
        if np.shape(data)[1] == 18:
            co_orig = data[:,15]
            cosd_orig = data[:,16]
            co_n = data[:,17]
            coflags = Makeflags(co_orig)
            co_stripped = Stripflags(co_orig)
            cosd_stripped = Stripflags(cosd_orig)
            
        # make datetime variable
        dt_date = [dt.datetime.strptime(date[i] + ' ' + time[i] , "%y%m%d %H%M%S") for i in np.arange(len(date))]
           
        self.date = date
        self.time = time
        self.dt_date = dt_date
        self.datetime = dt_date
        self.sampletype = sampletype
        self.samplename = samplename
        self.standard = standard
        self.port = port.astype('int')
        self.cavity_temp = cavity_temp.astype('float')
        self.cavity_press = cavity_press.astype('float')
        self.h2o_orig = h2o_orig
        self.h2o =  h2o_stripped.gas
        
        self.co2_orig = co2_orig
        self.co2 = co2_stripped.gas
        self.co2sd_orig = co2sd_orig
        self.co2sd = co2sd_stripped.gas
        self.co2flags = co2flags.flags
        self.co2_n = co2_n
        
        self.ch4_orig = ch4_orig
        self.ch4 = ch4_stripped.gas
        self.ch4sd_orig = ch4sd_orig
        self.ch4sd = ch4sd_stripped.gas
        self.ch4flags = ch4flags.flags
        self.ch4_n = ch4_n
        
        self.nogases = 2        
        
        # CO data is only at some sites
        if data[0][-1] != 'nan':
            self.co_orig = co_orig
            self.co = co_stripped.gas
            self.cosd_orig = cosd_orig
            self.cosd = cosd_stripped.gas
            self.coflags = coflags.flags
            self.co_n = co_n
            
            self.nogases = 3
        
        self.filename = filename
        self.datadir = dirname


# New file format
#Created: 11 Jun 15 12:00 GMT
#     -      -     bin    -         -            -            -    -    cavity    cavity         -       co2       co2       co2     co2   co2       ch4       ch4       ch4     ch4   ch4        co        co        co      co    co 
#  date   time    date time      type       sample     standard port      temp     press       h2o         C       wet       dry   stdev     N         C       wet       dry   stdev     N         C       wet       dry   stdev     N 
#150605 104930  150605 1049       std USN-20121603 USN-20121603    5    45.000   139.963     1.567    519.74F   503.22F   516.28F216.126    14      5.92F  1950.83F  1990.64F154.924    14       nan       nan       nan     nan   nan 


class read_gcexport_crds_new:
    def __init__(self, data, datafile):
        
        date = data[:,0]
        time = data[:,1]
        sampletype = data[:,4]
        samplename = data[:,5]
        standard = data[:,6]
        port = data[:,7]
        cavity_temp = data[:,8]
        cavity_press = data[:,9] 
        h2o_orig = data[:,10]
        
        co2_orig = data[:,11]
        co2wet_orig = data[:,12]
        co2dry_orig = data[:,13]
        co2sd_orig = data[:,14]
        co2_n = data[:,15]
        
        ch4_orig = data[:,16]
        ch4wet_orig = data[:,17]
        ch4dry_orig = data[:,18]
        ch4sd_orig = data[:,19]
        ch4_n = data[:,20]
        
        co_orig = data[:,21]
        cowet_orig = data[:,22]
        codry_orig = data[:,23]
        cosd_orig = data[:,24]
        co_n = data[:,25]
                
                
        dirname, filename = os.path.split(datafile)

        
        coflags = Makeflags(co_orig)
        co2flags = Makeflags(co2_orig)
        ch4flags = Makeflags(ch4_orig)
        
        # Strip flags from conc data
        co2_stripped = Stripflags(co2_orig)
        co2sd_stripped = Stripflags(co2sd_orig)
        ch4_stripped = Stripflags(ch4_orig)
        ch4sd_stripped = Stripflags(ch4sd_orig)
        h2o_stripped = Stripflags(h2o_orig)
        co_stripped = Stripflags(co_orig)
        cosd_stripped = Stripflags(cosd_orig)
        
        co2wet_stripped = Stripflags(co2wet_orig)
        co2dry_stripped = Stripflags(co2dry_orig)
        ch4wet_stripped = Stripflags(ch4wet_orig)
        ch4dry_stripped = Stripflags(ch4dry_orig)
        cowet_stripped = Stripflags(cowet_orig)
        codry_stripped = Stripflags(codry_orig)

            
        # make datetime variable
        dt_date = [dt.datetime.strptime(date[i] + ' ' + time[i] , "%y%m%d %H%M%S") for i in np.arange(len(date))]
           
        self.date = date
        self.time = time
        self.dt_date = dt_date
        self.datetime = dt_date
        self.sampletype = sampletype
        self.samplename = samplename
        self.standard = standard
        self.port = port.astype('int')
        self.cavity_temp = cavity_temp.astype('float')
        self.cavity_press = cavity_press.astype('float')
        self.h2o_orig = h2o_orig
        self.h2o =  h2o_stripped.gas
        
        self.co2_orig = co2_orig
        self.co2 = co2_stripped.gas
        self.co2sd_orig = co2sd_orig
        self.co2sd = co2sd_stripped.gas
        self.co2flags = co2flags.flags
        self.co2_n = co2_n
        
        self.ch4_orig = ch4_orig
        self.ch4 = ch4_stripped.gas
        self.ch4sd_orig = ch4sd_orig
        self.ch4sd = ch4sd_stripped.gas
        self.ch4flags = ch4flags.flags
        self.ch4_n = ch4_n
        
        self.co2_dry = co2dry_stripped.gas
        self.co2_wet = co2wet_stripped.gas
        self.ch4_dry = ch4dry_stripped.gas
        self.ch4_wet = ch4wet_stripped.gas    
        
        self.nogases = 2        
        
        # CO data is only at some sites
        if data[0][-1] != 'nan':
            self.co_orig = co_orig
            self.co = co_stripped.gas
            self.cosd_orig = cosd_orig
            self.cosd = cosd_stripped.gas
            self.co_dry = codry_stripped
            self.co_wet = cowet_stripped
            self.coflags = coflags.flags
            self.co_n = co_n
            
            self.nogases = 3
        
        self.filename = filename
        self.datadir = dirname

# Class to make the flags
# nb: an * flag will =  a 1 flag
# an F flag will = a 2 flag
class Makeflags:
    def __init__(self, gas):    
        
        
        
        flag = np.zeros(shape=len(gas))
        for i in np.arange(len(gas)):
            A = (gas[i]).find('*')
            B = (gas[i]).find('F')
            
            # a * flag
            if A > 0:
                flag[i] = 1

            # a F flag
            if B > 0:
                flag[i] = 2


        self.flags = flag.astype('int')
        

# Class to strip the * and F flags from the raw data
class Stripflags:
    def __init__(self, gas):    
        
        orig_gas = gas        
        stripped_gas = np.zeros(shape=len(gas))
        
        for i in np.arange(len(gas)):
            
            A = (orig_gas[i]).find('*')
            B = (orig_gas[i]).find('F')
            
            # a * flag
            if A > 0:
                stripped_gas[i] = orig_gas[i][:(len(orig_gas[i])-1)]


            # a F flag
            if B > 0:
                stripped_gas[i] = orig_gas[i][:(len(orig_gas[i])-1)]

            if A < 0 and B < 0:
                stripped_gas[i] = orig_gas[i]

        self.gas = stripped_gas.astype('float')
        self.orig_gas = orig_gas


# Code to extract the flagged data
def Extractgood(time, gas, flags, SD=0):
    
    index = (np.array(np.where(flags == 0)))[0,:]       
    time_flag = [time[i] for i in index]          
    gas_flag = gas[index]
    
    if type(SD) != int:
        sd_flag = SD[index]
    else:
        sd_flag = 0
        
    return time_flag, gas_flag, sd_flag
