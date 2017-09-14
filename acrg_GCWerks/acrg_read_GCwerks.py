# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:22:22 2015

@author: as13988
"""

import os
import numpy as np
import datetime as dt
import pdb
import copy 

# Class to read in the txt output of CRDS gcexport made using the standard peak.list and report config on Dagage2
def read_gcexport_crds(datafile):
        
        print 'Reading:'
        print datafile
        
        if type(datafile) == tuple:
            #data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=1)
            
        elif type(datafile) == str:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=1)
            
        elif type(datafile) == list:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            datafile = datafile[0]
            data=np.genfromtxt(datafile, dtype=str, skip_header=1)
        
        elif type(datafile) == np.string_:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=1)
            
        data_good = data[2:,:]
        header_1 = data[0,:]
                
        # Ensure that second "time" is "bin_time"
        header_1[np.where(header_1 =='bin')[0]+1] ='bin'        
        
        header_2 = data[1,:]
        headers = []
        
        for i in np.arange(len(header_1)):
            if header_1[i] == '-':    
                headers.append(header_2[i])
                
            else:
                headers.append(header_1[i] + '_' + header_2[i])
        headers = np.array(headers)
        
        #if np.shape(data)[1] > 18 :
         #   data = read_gcexport_crds_new(data_good, datafile)
        #else:
        #    data = read_gcexport_crds_old(data_good, datafile)
        
        data = read_gcexport_crds_flexi(data_good, datafile, headers)

        return data

# More flexible approach using the headers in the file
class read_gcexport_crds_flexi:
    def __init__(self, data, datafile, headers):
        
        date = data[:, np.where(headers == 'date')[0]][:,0]
        time = data[:, np.where(headers == 'time')[0]][:,0]
        sampletype = data[:,np.where(headers == 'type')[0]][:,0]
        samplename = data[:,np.where(headers == 'sample')[0]][:,0]
        standard = data[:,np.where(headers == 'standard')[0]][:,0]
        port = data[:,np.where(headers == 'port')[0]][:,0]
        cavity_temp = data[:,np.where(headers == 'cavity_temp')[0]][:,0]
        cavity_press = data[:,np.where(headers == 'cavity_press')[0]][:,0] 
        h2o_orig = data[:,np.where(headers == 'h2o')[0]][:,0]
        h2osd_orig = data[:,np.where(headers == 'h2o_stdev')[0]][:,0]
        
        co2_orig = data[:, np.where(headers == 'co2_C')[0]][:,0]
        co2drift_orig = data[:, np.where(headers == 'co2_Cdrift')[0]][:,0]        
        co2wet_orig = data[:,np.where(headers == 'co2_wet')[0]][:,0]
        co2dry_orig = data[:,np.where(headers == 'co2_dry')[0]][:,0]
        co2sd_orig = data[:,np.where(headers == 'co2_stdev')[0]][:,0]
        co2_n = data[:,np.where(headers == 'co2_N')[0]][:,0]
        
        ch4_orig = data[:, np.where(headers == 'ch4_C')[0]][:,0]
        ch4drift_orig = data[:, np.where(headers == 'ch4_Cdrift')[0]][:,0]        
        ch4wet_orig = data[:, np.where(headers == 'ch4_wet')[0]][:,0]
        ch4dry_orig = data[:, np.where(headers == 'ch4_dry')[0]][:,0]
        ch4sd_orig = data[:, np.where(headers == 'ch4_stdev')[0]][:,0]
        ch4_n = data[:, np.where(headers == 'ch4_N')[0]][:,0]
        
        dirname, filename = os.path.split(datafile)
        
        co2flags = Makeflags(co2_orig)
        ch4flags = Makeflags(ch4_orig)
        
        
        # Strip flags from conc data
        co2_stripped = Stripflags(co2_orig)
        co2drift_stripped = Stripflags(co2drift_orig)
        co2sd_stripped = Stripflags(co2sd_orig)
        
        ch4_stripped = Stripflags(ch4_orig)
        ch4drift_stripped = Stripflags(ch4drift_orig)
        ch4sd_stripped = Stripflags(ch4sd_orig)
        
        h2o_stripped = Stripflags(h2o_orig)
        h2osd_stripped = Stripflags(h2osd_orig)
        
        co2wet_stripped = Stripflags(co2wet_orig)
        co2dry_stripped = Stripflags(co2dry_orig)
        ch4wet_stripped = Stripflags(ch4wet_orig)
        ch4dry_stripped = Stripflags(ch4dry_orig)

        self.nogases = 2  
        
        if 'co_C' in headers:
            co_orig = data[:, np.where(headers == 'co_C')[0]][:,0]
            codrift_orig = data[:, np.where(headers == 'co_Cdrift')[0]][:,0]        
            cowet_orig = data[:, np.where(headers == 'co_wet')[0]][:,0]
            codry_orig = data[:, np.where(headers == 'co_dry')[0]][:,0]
            cosd_orig = data[:, np.where(headers == 'co_stdev')[0]][:,0]
            co_n = data[:, np.where(headers == 'co_N')[0]][:,0]
            
            coflags = Makeflags(co_orig)
            
            co_stripped = Stripflags(co_orig)
            codrift_stripped = Stripflags(codrift_orig)
            cosd_stripped = Stripflags(cosd_orig)

            cowet_stripped = Stripflags(cowet_orig)
            codry_stripped = Stripflags(codry_orig)

            self.co_orig = co_orig
            self.co = co_stripped.gas
            self.codrift = codrift_stripped.gas
            self.cosd_orig = cosd_orig
            self.cosd = cosd_stripped.gas
            self.co_dry = codry_stripped.gas
            self.co_wet = cowet_stripped.gas
            self.coflags = coflags.flags
            self.co_n = co_n
            
            self.nogases = 3

        if 'h2o_reported' in headers:
            h2o_rep_orig = data[:,np.where(headers == 'h2o_reported')[0]][:,0]
            h2o_rep_stripped = Stripflags(h2o_rep_orig)
            self.h2o_rep_orig = h2o_rep_orig
            self.h2o_reported =  h2o_rep_stripped.gas


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
        self.h2osd_orig = h2osd_orig
        self.h2osd = h2osd_stripped.gas
        
        self.co2_orig = co2_orig
        self.co2 = co2_stripped.gas
        self.co2drift = co2drift_stripped.gas
        self.co2sd_orig = co2sd_orig
        self.co2sd = co2sd_stripped.gas
        self.co2flags = co2flags.flags
        self.co2_n = co2_n
        
        self.ch4_orig = ch4_orig
        self.ch4 = ch4_stripped.gas
        self.ch4drift = ch4drift_stripped.gas
        self.ch4sd_orig = ch4sd_orig
        self.ch4sd = ch4sd_stripped.gas
        self.ch4flags = ch4flags.flags
        self.ch4_n = ch4_n
        
        self.co2_dry = co2dry_stripped.gas
        self.co2_wet = co2wet_stripped.gas
        self.ch4_dry = ch4dry_stripped.gas
        self.ch4_wet = ch4wet_stripped.gas    
        
              
        
        self.filename = filename
        self.datadir = dirname

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
        
        pdb.set_trace()        
        
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

# More flexible approach using the headers in the file
class read_gcexport_lgr:
    def __init__(self, datafile):
        
        print 'Reading:'
        print datafile
        if type(datafile) == tuple:
            #data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=1)
        elif type(datafile) == str:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=1)
        
              
        data_good = data[2:,:]
        header_1 = data[0,:]
              
        # Ensure that second "time" is "bin_time"
        header_1[np.where(header_1 =='bin')[0]+1] ='bin'        
        
        header_2 = data[1,:]
        headers = []
        
        for i in np.arange(len(header_1)):
            if header_1[i] == '-':    
                headers.append(header_2[i])
                
            else:
                headers.append(header_1[i] + '_' + header_2[i])
        headers = np.array(headers)
        
        
        date = data_good[:, np.where(headers == 'date')[0]][:,0]
        time = data_good[:, np.where(headers == 'time')[0]][:,0]
        sampletype = data_good[:,np.where(headers == 'type')[0]][:,0]
        samplename = data_good[:,np.where(headers == 'sample')[0]][:,0]
        standard = data_good[:,np.where(headers == 'standard')[0]][:,0]
        port = data_good[:,np.where(headers == 'port')[0]][:,0]
        cavity_temp = data_good[:,np.where(headers == 'cavity_temp')[0]][:,0]
        cavity_press = data_good[:,np.where(headers == 'cavity_press')[0]][:,0] 
        h2o_orig = data_good[:,np.where(headers == 'h2o')[0]][:,0]
        
        co_orig = data_good[:, np.where(headers == 'co_C')[0]][:,0]
        co_wet_orig = data_good[:,np.where(headers == 'co_wet')[0]][:,0]
        co_dry_orig = data_good[:,np.where(headers == 'co_dry')[0]][:,0]
        co_sd_orig = data_good[:,np.where(headers == 'co_stdev')[0]][:,0]
        co_n = data_good[:,np.where(headers == 'co_N')[0]][:,0]
        
        n2o_orig = data_good[:, np.where(headers == 'n2o_C')[0]][:,0]
        n2o_wet_orig = data_good[:, np.where(headers == 'n2o_wet')[0]][:,0]
        n2o_dry_orig = data_good[:, np.where(headers == 'n2o_dry')[0]][:,0]
        n2o_sd_orig = data_good[:, np.where(headers == 'n2o_stdev')[0]][:,0]
        n2o_n = data_good[:, np.where(headers == 'n2o_N')[0]][:,0]
        
        dirname, filename = os.path.split(datafile)
        
        coflags = Makeflags(co_orig)
        n2oflags = Makeflags(n2o_orig)
        
               
        
        # Strip flags from conc data
        co_stripped = Stripflags(co_orig)
        co_sd_stripped = Stripflags(co_sd_orig)
        n2o_stripped = Stripflags(n2o_orig)
        n2o_sd_stripped = Stripflags(n2o_sd_orig)
        h2o_stripped = Stripflags(h2o_orig)
        
        co_wet_stripped = Stripflags(co_wet_orig)
        co_dry_stripped = Stripflags(co_dry_orig)
        n2o_wet_stripped = Stripflags(n2o_wet_orig)
        n2o_dry_stripped = Stripflags(n2o_dry_orig)

        self.nogases = 2  


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
        
        self.co_orig = co_orig
        self.co = co_stripped.gas
        self.co_sd_orig = co_sd_orig
        self.co_sd = co_sd_stripped.gas
        self.co_flags = coflags.flags
        self.co_n = co_n
        
        self.n2o_orig = n2o_orig
        self.n2o = n2o_stripped.gas
        self.n2o_sd_orig = n2o_sd_orig
        self.n2o_sd = n2o_sd_stripped.gas
        self.n2o_flags = n2oflags.flags
        self.n2o_n = n2o_n
        
        self.co_dry = co_dry_stripped.gas
        self.co_wet = co_wet_stripped.gas
        self.n2o_dry = n2o_dry_stripped.gas
        self.n2o_wet = n2o_wet_stripped.gas    
        
              
        
        self.filename = filename
        self.datadir = dirname

# More flexible approach using the headers in the file
class read_gcexport_medusa:
    def __init__(self, datafile):
        
        print 'Reading:'
        print datafile
        
        if type(datafile) == tuple:
            #data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=1)
        elif type(datafile) == str:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=1)
        
              
        data_good = data[2:,:]
        header_1 = data[0,:]
              
        # Ensure that second "time" is "bin_time"
        header_1[np.where(header_1 =='bin')[0]+1] ='bin'        
        
        header_2 = data[1,:]
        headers = []
        
        for i in np.arange(len(header_1)):
            if header_1[i] == '-':    
                headers.append(header_2[i])
                
            else:
                headers.append(header_1[i] + '_' + header_2[i])
        headers = np.array(headers)
        
        date = data_good[:, np.where(headers == 'date')[0]][:,0]
        time = data_good[:, np.where(headers == 'time')[0]][:,0]
        sampletype = data_good[:,np.where(headers == 'type')[0]][:,0]
        samplename = data_good[:,np.where(headers == 'sample')[0]][:,0]
        standard = data_good[:,np.where(headers == 'standard')[0]][:,0]
        port = data_good[:,np.where(headers == 'port')[0]][:,0]
        volume = data_good[:,np.where(headers == 'volume')[0]][:,0]
        
        count = 0
        species = []
        for i in headers:
            if '_C' == i[-2:]:
                orig_i =  data_good[:, np.where(headers == i)[0]][:,0]
                flags_i = Makeflags(orig_i)
                stripped_i = Stripflags(orig_i)
                
                orig_i_tag = i.split('_')[0] + '_orig'          
                flags_i_tag = i.split('_')[0] + '_flags' 
                
                setattr(self, orig_i_tag, orig_i)
                setattr(self, i, stripped_i.gas)
                setattr(self, flags_i_tag, flags_i.flags)
                
                count = count + 1
                
                species.append(i.split('_')[0])

        self.nogases = count
        self.species = species

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
        self.volume =  volume
        
        self.filename = os.path.basename(datafile)
        self.datadir = os.path.dirname(datafile)

        
# Class to make the flags
# nb: an * flag will =  a 1 flag
# an F flag will = a 2 flag
class Makeflags:
    def __init__(self, gas):    
        
        
        
        flag = np.zeros(shape=len(gas))
        for i in np.arange(len(gas)):
            A = (gas[i]).find('*')
            B = (gas[i]).find('F')
            C = (gas[i]).find('x')
            D = (gas[i]).find('H')
            
            # a * flag
            if A > 0:
                flag[i] = 1

            # a F flag
            if B > 0:
                flag[i] = 2

            # a x flag
            if C > 0:
                flag[i] = 3
                
            # These are ok it just means that they're done by height
            if D > 0:
                flag[i] = 0 
            
        self.flags = flag.astype('int')
        

# Class to strip the *,x and F flags from the raw data
class Stripflags:
    def __init__(self, gas):    
        
        orig_gas = gas        
        stripped_gas = np.zeros(shape=len(gas))
        
        # need to check if it's all nan's
        for i in np.arange(len(gas)):
                    
            flags = ['*','F','x','T','H']
            
            # one of the flags is in the string
            if any(x in orig_gas[i] for x in flags):
                stripped_gas[i] = orig_gas[i][:(len(orig_gas[i])-1)]

            else:
                stripped_gas[i] = orig_gas[i]
                
        """ Old method
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
        """
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
    
# Class to read in the txt output of gcexport made using the standard peak.list and report config on Dagage2
def read_gcexport_means(infile):

        # Adjust if the datafile is given as a tuple
        if type(infile) == tuple:
             infile = infile[0]
             
        print 'Reading:'
        print infile
        
        # Read in the end header line and the rest of the file
        data=np.genfromtxt(infile, dtype=str, skip_header=1)

        # Split out the header              
        data_good = data[1:,:]
        
        # Need to skip header line 0 and then read it in separately!
        header_2 = data[0,:]

       
                
        # read in the first header line
        # NB: I can't read this in at the same time as the first read as it has a different number of columns!
        header_1 = np.genfromtxt(infile, dtype=str, skip_footer=len(data))
        
        headers = []
        
        for i in np.arange(len(header_1)):
            headers.append(header_1[i] + '_' + header_2[i+3])

        headers = np.concatenate((header_2[0:3], np.array(headers)))
        
        data = read_gcexport_means_flexi(data_good, infile, headers)

        return data
        
# Reads the standard ports logs formatted files
def read_portslogs(site = 'MHD', std_only = None, logdir = '/Users/as13988/Documents/Work/Cylinders/Calibrations/N2OScaleConversion/'):
    import glob
    logfiles = glob.glob(logdir+ site + '_Tanks/portslogs'+ '/????')    
    
    datetime_out = []
    port_out = []    
    type_out = []
    tank_out = []
    reg_out = []
    
    for i in logfiles:
        print 'Reading file : ' + i
        indata = np.genfromtxt(i, dtype=str)
    
        
        
        if len(np.shape(indata)) > 1:
            date_str = indata[:,0]
            time_str = indata[:,1]
            port_out.extend([int(j) for j in np.array(indata[:,2])])
            tank_out.extend(indata[:,3])
            reg_out.extend(indata[:,4])
            type_out.extend(indata[:,5])
        else:
            date_str = [indata[0]]
            time_str = [indata[1]]
            port_out.append(int(indata[2]))
            tank_out.append(indata[3])
            reg_out.append(indata[4])
            type_out.append(indata[5])

        
        
        for k in np.arange(len(date_str)):
            datetime_out.append(dt.datetime.strptime(date_str[k]+' ' + time_str[k], '%y%m%d %H%M'))


    if std_only is not None:
        index = np.where(np.array(type_out) == 'std')[0]    

        port_out = np.array(port_out)[index]
        tank_out = np.array(tank_out)[index]
        reg_out = np.array(reg_out)[index]
        type_out = np.array(type_out)[index]
        datetime_out = [datetime_out[i] for i in index]

    else:
        
        port_out = np.array(port_out)
        tank_out = np.array(tank_out)
        reg_out = np.array(reg_out)
        type_out = np.array(type_out)

    dateindex = np.argsort(datetime_out)
    
    output = {}
    output['datetime'] = np.array(datetime_out)[dateindex]
    output['tank'] = tank_out[dateindex]
    output['port'] = port_out[dateindex]
    output['reg'] = reg_out[dateindex]
    output['type'] = type_out[dateindex]
    
    return output


# More flexible approach using the headers in the file
class read_gcexport_means_flexi:
    def __init__(self, data, infile, headers):
        
        # Extract columns common to all means files
        tank = data[:, np.where(headers == 'Tank')[0]][:,0]
        std = data[:, np.where(headers == 'Standard')[0]][:,0]
        date = data[:, np.where(headers == 'Date')[0]][:,0]
        
        # Want to loop through the remainder of the column headers and store them in a dictionary
        dataout = {}
        for i in headers[3:-1]:
            dataout[i] = (data[:, np.where(headers == i)[0]][:,0]).astype('float')

        # make datetime variable
        dt_date = [dt.datetime.strptime(i, "%y%m%d") for i in date]
           
        self.date = date
        self.datetime = dt_date
        self.standard = std
        self.tank = tank
        self.data = dataout
        
        dirname, filename = os.path.split(infile)
          
        self.filename = filename
        self.datadir = dirname

# Read in a MD standards file
# NOTE: I HAD TO COMMENT OUT A HEADER LINE THAT SOME GOOSE PUT IN HALFWAY THROUGH THE FILE!!!!!
def read_MDstds(infile='/Users/as13988/Documents/Work/Cylinders/Calibrations/N2OScaleConversion/MHD_Tanks/standards'):
    raw=[]    
    out =[]
    
    # read in the data file, ignore comment lines and split based on while space
    first_line = 0
    for line in open(infile):
        
        li=line.strip()
        
        if not li.startswith("#"):
            
            raw.append(line.rstrip())
            
            if first_line == 0:
                headers = line.split()
                first_line = 1
            else:
                out.append(line.split())
                
    
    # convert list/array to a dictionary labelled using the column headers
    out_dict={}
    
    for i in np.arange(len(headers)):
        
        if i == 0:
            out_dict['tank'] = np.array(out)[:,i]
        else:
            out_dict[headers[i]] = (np.array(out)[:,i]).astype(float)
          
    return out_dict

# Read in MD standards.factors file
# This assumed that there's only 4 columns of interest
def read_MDstdsfactors(infile='/Users/as13988/Documents/Work/Cylinders/Calibrations/N2OScaleConversion/MHD_Tanks/standards.factors'):
    raw=[]    
    out =[]
    
    # read in the data file, ignore comment lines and split based on while space
    first_line = 0
    for line in open(infile):
        
        li=line.strip()
        
        if not li.startswith("#"):
            
            raw.append(line.rstrip())
            
            stripped_l = line.split('#') # strip off the comment at the end!
            l = stripped_l[0].split()
            if len(l) > 1:
                out.append(stripped_l[0].split())

# convert list/array to a dictionary labelled using the column headers
    out_dict={}
    headers = ['species', 'scale','factor','units']
    
    for i in np.arange(len(headers)):
        if headers[i] == 'factor':
            
            out_dict[headers[i]] = (np.array(out)[:,i]).astype(float)
        else:
            out_dict[headers[i]] = np.array(out)[:,i]
        
    return out_dict

# Class to read in the txt output of gcexport made using the standard peak.list and report config on Dagage2
def read_gcexport_md(datafile='/Users/as13988/Documents/Work/Cylinders/Calibrations/N2OScaleConversion/MHD_Tanks/MHD_air.dat'):
        
        print 'Reading:'
        print datafile
        if type(datafile) == tuple:
            #data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
        elif type(datafile) == str:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=3)
        
        
        header_1 = np.genfromtxt(datafile, dtype=str, skip_header=1, skip_footer = len(data)+1)
        header_2 = np.genfromtxt(datafile, dtype=str, skip_header=2, skip_footer = len(data))
                
        headers = []
        
        for i in np.arange(len(header_2)):
            if header_1[i] == '-':    
                headers.append(header_2[i])
                
            else:
                headers.append(header_1[i] + '_' + header_2[i])

        # remove flagged lines

        # remove lines where conc  = nan
        # I think this removes all the flagged points
        header_index = np.where(np.array(headers) == 'N2O_C')[0]
        index = np.where(data[:,header_index] != 'nan')[0]

        data = data[index,:]
    

        out_dict={}

        # convert to a dictionary
        for i in np.arange(len(headers)):
            if headers[i] in ['date','time', 'type', 'standard']:
                # leave as string
                out_dict[headers[i]] = data[:,i]
            
            else:
                # convert to float
                out_dict[headers[i]] = (data[:,i]).astype(float)


        # make a datetime variable
        dt_date = [dt.datetime.strptime(out_dict['date'][i] + ' ' + out_dict['time'][i] , "%y%m%d %H%M%S") for i in np.arange(len(out_dict['date']))]

        out_dict['datetime'] = dt_date
        return out_dict

# Class to read in the txt output of gcexport made using the standard peak.list and report config on Dagage2
# The tank format is different to the air format
def read_gcexport_md_tank(datafile='/Users/as13988/Documents/Work/Cylinders/Calibrations/N2OScaleConversion/MHD_Tanks/MHD_air.dat'):
        
        print 'Reading:'
        print datafile
        if type(datafile) == tuple:
            #data=np.genfromtxt(datafile[0], dtype=str, skip_header=3)
            data=np.genfromtxt(datafile[0], dtype=str, skip_header=2)
        elif type(datafile) == str:
            #data=np.genfromtxt(datafile, dtype=str, skip_header=3)
            data=np.genfromtxt(datafile, dtype=str, skip_header=2)
        
        
        header_1 = np.genfromtxt(datafile, dtype=str, skip_footer = len(data)+1)
        header_2 = np.genfromtxt(datafile, dtype=str, skip_header=1, skip_footer = len(data))
                
        headers = list(header_2[0:len(header_2)-len(header_1)])
        
        for i in np.arange(len(header_1)):
           headers = headers + [(header_1[i] + '_' + header_2[i+len(header_2)-len(header_1)])]

        # remove flagged lines
        # remove lines where conc  = nan
        # I think this removes all the flagged points
        header_index = np.where(np.array(headers) == 'N2O_C')[0]
        index = np.where(data[:,header_index] != 'nan')[0]

        data = data[index,:]
    

        out_dict={}

        # convert to a dictionary
        for i in np.arange(len(headers)):
            if headers[i] in ['date','time', 'type', 'standard','Tank','Standard', 'Date']:
                # leave as string
                out_dict[headers[i]] = data[:,i]
            
            else:
                # convert to float
                out_dict[headers[i]] = (data[:,i]).astype(float)


        # make a datetime variable
        dt_date = [dt.datetime.strptime(out_dict['Date'][i], "%y%m%d") for i in np.arange(len(out_dict['Date']))]

        out_dict['datetime'] = dt_date
        return out_dict
