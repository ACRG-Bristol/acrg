# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 16:17:55 2015

@author: as13988
"""

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
import acrg_read_GCwerks as read_GCwerks

# Class that contains the USNs and associated DNos
class USNsDNos():
   def __init__(self): 
       USNs = ['USN20133143', 'USN20132963', 'USN20133065']
       DNos = ['D091968', 'D091969', 'D091970']  
       
       self.USNs = USNs
       self.DNos = DNos

# Class that matches each site to a specific colour
class SiteColours():
    def __init__(self, sitein = 0):
        colours = ['blue', 'red', 'orange', 'green', 'dodgerblue', 'magenta', 'orchid', 'purple', 'black', 'sienna', 'gray', 'limegreen']
        sites = ['RGL', 'TAC', 'HFD', 'BSD', 'TTA', 'FAAM', 'Man', 'BAS', 'MPI', 'FERRY']
        
        if sitein != 0:
            index = np.where(np.array(sites) == sitein)[0]
            
            if len(index) == 0:
                index = -1
                
        self.colours = colours
        self.sites = sites
        self.outcolour = colours[index]

# Class to calculate the wmo precision bars
class WMOLines():
    def __init__(self, conc, gas = 'CO2'):
        
        gases = ['CO2', 'CH4', 'CO', 'N2O', 'SF6']
        compatibility = [0.1, 2, 2, 0.1, 0.02]
        units = ['ppm', 'ppb','ppb','ppb','ppt']
        
        gas_index = np.where(np.array(gases) == gas)[0]
  
        # calculate the mean conc
        conc = np.array(conc)
        notnan_conc = conc[~np.isnan(conc)]
        
        mean_conc = np.zeros(2)
        mean_conc[:] = np.mean([notnan_conc])  
  
        self.plusline = mean_conc + compatibility[gas_index]
        self.minusline = mean_conc - compatibility[gas_index]
        self.units = units[gas_index]
        self.compatibility = compatibility[gas_index]

 # Class to read in the txt output of gcexport made using ICPDataExtractScript on Dagage2
# file format is 3 header lines

#:Created: 14 Jan 15 15:28 GMT
#     -      -         -            -            -    -    cavity    cavity         -       co2     co2       ch4     ch4        co      co 
#  date   time      type       sample     standard port      temp     press       h2o         C   stdev         C   stdev         C   stdev 
# 140710 093530      tank  USN20132963        H-239    6    45.000   139.999     0.532    418.42*  0.188*  1993.75*  0.204*   200.84*  5.913*

class read_data:
    def __init__(self, datafile):
        
        # Read in the data using the general code
        indata =  read_GCwerks.read_gcexport_crds(datafile)    
            
        # Find making USN or vice versa
        if (indata.filename.split('.')[2]).find('USN') != -1:
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
  
        scales = Read_scales()   
        # Match location to scale        
        in_sitekey =  indata.filename.split('.')[0] +'_'+ (indata.filename.split('.')[1]).split('_')[0]
        scale_index = np.where(scales.sitekey == in_sitekey)

        print in_sitekey
        print scale_index
        print scales.co2_scale[scale_index]
        print scales.ch4_scale[scale_index]
        print scales.co_scale[scale_index]
        

           
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
        self.co2_scale = scales.co2_scale[scale_index]
        
        self.ch4_orig = indata.ch4_orig
        self.ch4 = indata.ch4
        self.ch4sd_orig = indata.ch4sd_orig
        self.ch4sd = indata.ch4sd
        self.ch4flags = indata.ch4flags
        self.ch4_n = indata.ch4_n
        self.ch4_scale = scales.ch4_scale[scale_index]
        
        self.co_orig = indata.co_orig
        self.co = indata.co
        self.cosd_orig = indata.cosd_orig
        self.cosd = indata.cosd
        self.coflags = indata.coflags
        self.co_n = indata.co_n
        self.co_scale = scales.co_scale[scale_index]        
        
        self.nogases = indata.nogases        
                
        self.filename = indata.filename
        self.datadir = indata.datadir

# Plotting the minute means
class PlotRawMM:
    def __init__(self, data, outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', plotflagged=0):
     
     
        if plotflagged == 0:
            
            dt_co2, co2, co2sd = read_GCwerks.Extractgood(data.datetime, data.co2, data.co2flags, data.co2sd)
            dt_ch4, ch4, ch4sd = read_GCwerks.Extractgood(data.datetime, data.ch4, data.ch4flags, data.ch4sd)
            
            if data.nogases == 3:     
                 dt_co, co, cosd = read_GCwerks.Extractgood(data.datetime, data.co, data.co2flags, data.cosd)
        else:
            dt_co2 = data.datetime            
            co2 = data.co2
            co2sd = data.co2sd
            dt_ch4 = data.datetime            
            ch4 = data.ch4
            ch4sd = data.ch4sd    
            
            if data.nogases == 3:     
                 dt_co = data.datetime                 
                 co = data.co
          
          
        if data.nogases ==3:
             # Plot of 
             plt.subplot(3, 1, 1)
             plt.errorbar(dt_co2, co2, yerr=co2sd, fmt='y-')
             plt.ylim(Calcrange(co2, co2sd))
             plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
             plt.locator_params(axis ='y', nbins = 6)
             
             plt.title('All Minute Means - Unflagged')
             plt.ylabel('[CO$_2$] (ppm)')
        
             plt.subplot(3, 1, 2)
             plt.errorbar(dt_ch4, ch4, yerr=ch4sd, fmt='r-')
             plt.ylim(Calcrange(ch4, ch4sd))
             plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
             plt.locator_params(axis ='y', nbins = 6)
             plt.ylabel('[CH$_4$] (ppb)')
        
             plt.subplot(3, 1, 3)
             plt.errorbar(dt_co, co, yerr=cosd, fmt='b-')
             plt.ylim(Calcrange(co, cosd))
             plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
             plt.locator_params(axis ='y', nbins = 6)
             plt.ylabel('[CO] (ppb)')
             plt.xlabel('time')
        
               
             if plotflagged == 0:
                 plt.savefig(outputdir+'RawMM_'+data.filename[:-3]+'png', dpi=100)
             else:
                 plt.savefig(outputdir+'RawUnflagged_MM'+data.filename[:-3]+'png', dpi=100)
                 
             plt.show()
             
             
        if data.nogases == 2:
            # Plot of 
             plt.subplot(2, 1, 1)
             plt.errorbar(dt_co2, co2, yerr=co2sd, fmt='y-')
             plt.ylim(Calcrange(co2, co2sd))
             plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
             plt.locator_params(axis ='y', nbins = 6)
             plt.title('All Minute Means - Unflagged')
             plt.ylabel('[CO$_2$] (ppm)')
        
             plt.subplot(2, 1, 2)
             plt.errorbar(dt_ch4, ch4, yerr=ch4sd, fmt='r-')
             plt.ylim(Calcrange(ch4, ch4sd))
             plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
             plt.locator_params(axis ='y', nbins = 6)
             plt.ylabel('[CH$_4$] (ppb)')
             plt.xlabel('time')
        
               
             if plotflagged == 0:
                 plt.savefig(outputdir+'RawMM_'+data.filename[:-3]+'png', dpi=100)
             else:
                 plt.savefig(outputdir+'RawUnflagged_MM'+data.filename[:-3]+'png', dpi=100)
                 
             plt.show()

    
# code to make a range for the plot
# +/- 20% of the actual range
def Calcrange(data, sd, even = 0):
    extra = ((np.nanmax(data)+np.nanmax(sd)) - (np.nanmin(data)-np.nanmax(sd)))/10

    datarange = [np.nanmin(data)- (3*extra), np.nanmax(data)+(3*extra)]
    
    if even != 0:
        biggest = np.max(abs(np.array(datarange)))
        datarange = [(-1*biggest), biggest]
    return datarange


# Code to determine the number of runs and calculate mean and sds for each run and an overall mean and sd
# Note this only uses unflagged data!!!
# use read_data to read in the raw data. It does all the flag removal
class Calcmeans:
   def __init__(self, data):

    dt_co2, co2, co2sd = read_GCwerks.Extractgood(data.datetime, data.co2, data.co2flags)
    
    co2runmeans, co2runsds, co2runcount = Calcrunmeans(dt_co2, co2)
    
    dt_ch4, ch4, ch4sd = read_GCwerks.Extractgood(data.datetime, data.ch4, data.ch4flags)
    
    ch4runmeans, ch4runsds, ch4runcount = Calcrunmeans(dt_ch4, ch4)
   
   
    if data.nogases == 3:
        dt_co, co, ch4sd = read_GCwerks.Extractgood(data.datetime, data.co, data.coflags)
    
        corunmeans, corunsds, coruncount = Calcrunmeans(dt_co, co)
        
    
    co2_scale = data.co2_scale[0]
    ch4_scale = data.ch4_scale[0]
    co_scale = data.co_scale[0]
    
    
    Means = {'DNo' : data.DNo, \
            'USN' : data.USN, \
            'nogases' : data.nogases, \
            'processeddate' : data.processeddate, \
            'measurementdate' : data.date[0] + ' ' + data.time[0], \
            'filename': data.filename, \
            'co2' : np.mean(co2), \
            'co2_sd' : np.std(co2), \
            'co2_n' : len(co2), \
            'co2runmeans' : co2runmeans, \
            'co2runsds' : co2runsds, \
            'co2runcount' : co2runcount, \
            'co2_scale' : co2_scale,\
            'ch4' : np.mean(ch4), \
            'ch4_sd' : np.std(ch4), \
            'ch4_n' : len(ch4), \
            'ch4runmeans' : ch4runmeans, \
            'ch4runsds' : ch4runsds, \
            'ch4runcount' : ch4runcount, \
            'ch4_scale' : ch4_scale }
        
    if data.nogases == 3:     
        Means['co'] = np.mean(co)
        Means['co_sd'] = np.std(co)
        Means['co_n'] = len(co)
        Means['corunmeans'] = corunmeans
        Means['corunsds'] = corunsds
        Means['coruncount'] = coruncount
        Means['co_scale'] = co_scale
    else:
        Means['co'] = np.nan
        Means['co_sd'] = np.nan
        Means['co_n'] = np.nan
        Means['corunmeans'] = np.nan
        Means['corunsds'] = np.nan
        Means['coruncount'] = np.nan
        Means['co_scale'] = np.nan
     
    self.means = Means
    
# Code to separate out each run and calculate means and stds
def Calcrunmeans(time, gas):
    
    timeinmin = np.array([((time[i]-time[0]).total_seconds())/60 for i in np.arange(len(time))]) # converts time to number of minutes since the first timestamp
    
    diff = np.diff(timeinmin) # difference between consecutive minutes
    
    gaps = np.append(np.insert(np.where(diff > 19)[0],0,0),len(diff))
    
    runmeans = np.zeros(shape=len(gaps)-1)  
    runsds  = np.zeros(shape=len(gaps)-1)   
    runcount = np.zeros(shape=len(gaps)-1)  
    
    for i in np.arange(len(gaps)-1):
        runmeans[i]=np.mean(gas[gaps[i]:gaps[i+1]])
        runsds[i]=np.std(gas[gaps[i]:gaps[i+1]])   
        runcount[i]=len(gas[gaps[i]:gaps[i+1]])   
        
    return runmeans, runsds, runcount        
        
        
# Code to run for multiple files based on UAN or DNo
class Calcmulti:
   def __init__(self, CylinderNo, basedir = '/Users/as13988/Documents/Work/Cylinders/ICP/Results/RawData/'):
       
       # Find matching USN or vice versa
        if (CylinderNo).find('USN') != -1:
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
            print 'Reading file: ' + files[i]
            data_i = read_data(files[i])
            
            PlotRawMM(data_i)
            
            #pdb.set_trace()            
            
            means_i = Calcmeans(data_i)

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


# Code to write output files containing the run means and overall means for each instrument
# I haven't written this yet


      
# Code to print out the output for multiple files in a useful format
# NB: INPUT IS THE OUTPUT OF CALCMUILTI
# Format
# USN & DNo
# Date file made
# Inst/site tag, date, co2, co2_sd, co2_n, ch4, ch4_sd, ch4_n, co, co_sd, co_n, dateprocess, filename
class Printmeans:
    def __init__(self, means, outputdir = '/Users/as13988/Documents/Work/Cylinders/ICP/Results/Processed/CO2CH4CO/'):
        
        
        if type(means) != type(dict()):
            means = means.means     
        
        keys =  means.keys()        
        
        # Make the header lines
        Header_1 = [means[keys[0]]['USN'] + ' ' +  means[keys[0]]['DNo'],' ']
        Header_2 = [time.ctime(), ' ']
        Header_3 = ['Inst/site tag', 'date', 'co2', 'co2_sd', 'co2_n', 'co2_scale', 'ch4', 'ch4_sd', 'ch4_n','ch4_scale', 'co', 'co_sd', 'co_n', 'co_scale', 'GCWerks output date', 'filename']
        
        
        # Construct the columns    
        column_names = Header_3[1:]
        column_names[0] ='measurementdate'
        column_names[13] ='processeddate'
          
        outrow = []          
          
        for j in keys: # j is each of the elements in keys
            outlist = [j]
            
            for i in column_names: # i is each of the column names
                outlist.append(str(means[j][i]))
            
            outrow.append(outlist)
                
        #pdb.set_trace()    
            
        with open(outputdir + means[keys[0]]['USN'] + '_' +  means[keys[0]]['DNo'] + '_Means.csv', 'wb') as csvfile:
            
            print 'Filename : ' + outputdir + means[keys[0]]['USN'] + '_' +  means[keys[0]]['DNo'] + '_Means.csv'
            
            linewriter = csv.writer(csvfile, delimiter=',')
                
            linewriter.writerow(Header_1)
            linewriter.writerow(Header_2)
            linewriter.writerow(Header_3)
                
            for i in np.arange(len(keys)):
                linewriter.writerow(outrow[i])
                

        self.outarray = outrow
 

# Code to read in the files printed by Printmeans
class Read_co2_meansmulti:
    def __init__(self, files = 1, basedir =  '/Users/as13988/Documents/Work/Cylinders/ICP/Results/Processed/CO2CH4CO/', CylinderNo = 0):
       
        
        if type(files) == int:
            # No files given
    
            # Check if a CylinderNo is given
            # Cylinder number not given
            if type(CylinderNo) == int:
                print 'You need to give a CylinderNo OR a list of files'
                
                self.means = ['You did not give a CylinderNo OR a list of files so I could not read in any data']
            
            # Cylinder number given
            else: 
                # Find matching USN or vice versa
                if (CylinderNo).find('USN') != -1:
                    # We've been given a USN so find the DNo
                    USN_tank = CylinderNo
                    index = (USNsDNos().USNs).index(USN_tank)
                    DNo_tank = (USNsDNos().DNos)[index]
                    
                else:
                    # We've been given a DNo so find the USN
                    DNo_tank = CylinderNo
                    index = (USNsDNos().DNos).index(DNo_tank)
                    USN_tank = (USNsDNos().USNs)[index]         
                
                
        files = glob.glob(basedir+'*'+USN_tank+'*.csv')
    
    
        sitekey = []
        datetime = []

        co2 = []
        co2sd = []
        co2_n = []
        co2_scale = []

        ch4 = []
        ch4sd = []
        ch4_n = []
        ch4_scale = []

        co = []
        cosd = []
        co_n = []
        co_scale = []
        
        raw_filename = []
        site = []
        instrument = []
        processeddate = []    

        for i in np.arange(len(files)):
            print 'Reading file: ' + files[i]
            
            data_i = Read_co2_means(files[i], DNo_tank, USN_tank)
                
            #pdb.set_trace()
            
            sitekey.append(data_i.sitekey)
            datetime.append(data_i.datetime)

            co2.append(data_i.co2)
            co2sd.append(data_i.co2sd)
            co2_n.append(data_i.co2_n)
            co2_scale.append(data_i.co2_scale)

            ch4.append(data_i.ch4)
            ch4sd.append(data_i.ch4sd)
            ch4_n.append(data_i.ch4_n)
            ch4_scale.append(data_i.ch4_scale)

            co.append(data_i.co)
            cosd.append(data_i.cosd)
            co_n.append(data_i.co_n)
            co_scale.append(data_i.co_scale)
            
            raw_filename.append(data_i.filename)
            site.append(data_i.site)
            instrument.append(data_i.instrument)
            processeddate.append(data_i.processeddate)
               
               
               
        self.sitekey = list(chain(*sitekey))
        
        
        datetimes = []        
        
        for i in np.arange(len(datetime)):
            for j in np.arange(len(datetime[i])):
                datetimes.append(datetime[i][j])        
        
        self.datetime = datetimes
        self.co2 = list(chain(*co2))
        self.co2sd = list(chain(*co2sd))
        self.co2_n = list(chain(*co2_n))
        self.co2_scale = list(chain(*co2_scale))
        self.co2_diffs = list(chain(*co2)) - np.nanmean(list(chain(*co2)))
        
        self.ch4 = list(chain(*ch4))
        self.ch4sd = list(chain(*ch4sd))
        self.ch4_n = list(chain(*ch4_n))
        self.ch4_scale = list(chain(*ch4_scale)) 
        self.ch4_diffs = list(chain(*ch4)) - np.nanmean(list(chain(*ch4)))
        
        self.co = list(chain(*co))
        self.cosd = list(chain(*cosd))
        self.co_n = list(chain(*co_n))
        self.co_scale = list(chain(*co_scale))
        self.co_diffs = list(chain(*co)) - np.nanmean(list(chain(*co)))
        
        self.raw_filename = list(chain(*raw_filename))
        self.site = list(chain(*site))
        self.instrument = list(chain(*instrument))
        self.processeddate = list(chain(*processeddate))
        self.processed_filenames = files
        self.USN = USN_tank
        self.DNo = DNo_tank             
     
     
     
 # Code to read in the files printed by Printmeans
class Read_n2o_meansmulti:
    def __init__(self, files = 1, basedir =  '/Users/as13988/Documents/Work/Cylinders/ICP/Results/Processed/N2OSF6/', CylinderNo = 0):
       
        
        if type(files) == int:
            # No files given
    
            # Check if a CylinderNo is given
            # Cylinder number not given
            if type(CylinderNo) == int:
                print 'You need to give a CylinderNo OR a list of files'
                
                self.means = ['You did not give a CylinderNo OR a list of files so I could not read in any data']
            
            # Cylinder number given
            else: 
                # Find matching USN or vice versa
                if (CylinderNo).find('USN') != -1:
                    # We've been given a USN so find the DNo
                    USN_tank = CylinderNo
                    index = (USNsDNos().USNs).index(USN_tank)
                    DNo_tank = (USNsDNos().DNos)[index]
                    
                else:
                    # We've been given a DNo so find the USN
                    DNo_tank = CylinderNo
                    index = (USNsDNos().DNos).index(DNo_tank)
                    USN_tank = (USNsDNos().USNs)[index]         
                
                
        files = glob.glob(basedir+'*'+USN_tank+'*.csv')
    
    
        sitekey = []
        datetime = []

        n2o = []
        n2osd = []
        n2o_n = []
        n2o_scale = []

        sf6 = []
        sf6sd = []
        sf6_n = []
        sf6_scale = []

        raw_filename = []
        site = []
        instrument = []
        
        for i in np.arange(len(files)):
            
            data_i = Read_n2o_means(files[i], DNo_tank, USN_tank)
            
            if type(data_i.sitekey) != np.string_:
                for j in data_i.sitekey:
                    sitekey.append(j)
            else:
                sitekey.append(data_i.sitekey)
                
            datetime.append(data_i.datetime)

            n2o.append(data_i.n2o)
            n2osd.append(data_i.n2osd)
            n2o_n.append(data_i.n2o_n)
            
                
            
            if type(data_i.n2o_scale) != np.string_ :
                for k in data_i.n2o_scale:
                    n2o_scale.append(k)
            else:
                n2o_scale.append(data_i.n2o_scale)   

            sf6.append(data_i.sf6)
            sf6sd.append(data_i.sf6sd)
            sf6_n.append(data_i.sf6_n)
            
            
            if type(data_i.sf6_scale) != np.string_ :
                for l in data_i.sf6_scale:
                    sf6_scale.append(l)
            else:
                sf6_scale.append(data_i.sf6_scale)            
            
       

            raw_filename.append(data_i.filename)
            
            
            if type(data_i.site) != str :
                for n in data_i.site:
                    site.append(n)
            else:
                site.append(data_i.site)            
            
       
            
           
            
            if type(data_i.instrument) != str :
                for m in data_i.instrument:
                    instrument.append(m)
            else:
                instrument.append(data_i.instrument)
        

                
           
        self.sitekey = sitekey
        
        
        datetimes = []        
        
        for i in np.arange(len(datetime)):
            if type(datetime[i]) == list:
                for j in np.arange(len(datetime[i])):
                    datetimes.append(datetime[i][j])        
            else:
                datetimes.append(datetime[i])  
                
        self.datetime = datetimes
        
        self.n2o = list(chain(*n2o))
        self.n2osd = list(chain(*n2osd))
        self.n2o_n = list(chain(*n2o_n))
        self.n2o_scale = n2o_scale
        self.n2o_diffs = list(chain(*n2o)) - np.nanmean(list(chain(*n2o)))
        
        self.sf6 = list(chain(*sf6))
        self.sf6sd = list(chain(*sf6sd))
        self.sf6_n = list(chain(*sf6_n))
        self.sf6_scale = sf6_scale
        self.sf6_diffs = list(chain(*sf6)) - np.nanmean(list(chain(*sf6)))
        
        self.raw_filename = list(chain(*raw_filename))
        self.site = site
        self.instrument = instrument
        self.processed_filenames = files
        self.USN = USN_tank
        self.DNo = DNo_tank   
        
        
# Code to read in the files printed by Printmeans
class Read_co2_means:
    def __init__(self, infile, DNo_tank, USN_tank):
       
        if type(infile) == tuple:
            data=np.genfromtxt(infile[0], dtype=str, skip_header=3, delimiter = ',')
        elif type(infile) == str:
            data=np.genfromtxt(infile, dtype=str, skip_header=3, delimiter = ',')     
        
        sitekey = data[:,0]
        date_time = data[:,1]
        co2 = data[:,2]
        co2sd = data[:,3]
        co2_n = data[:,4]
        co2_scale = data[:,5]
        ch4 = data[:,6]
        ch4sd = data[:,7]
        ch4_n = data[:,8]
        ch4_scale = data[:,9]
        co = data[:,10]
        cosd = data[:,11]
        co_n = data[:,12]
        co_scale = data[:,13]
        GCWerks_outputdate = data[:,14]
        filename = data[:,15]
            
        # make datetime variable
        dt_date = [dt.datetime.strptime(date_time[i] , "%y%m%d %H%M%S") for i in np.arange(len(date_time))]
        
        
        #pdb.set_trace()
        
        self.sitekey = sitekey
        self.datetime = dt_date
    
        self.co2 = [float(i) for i in co2]
        self.co2sd = [float(i) for i in co2sd]
        self.co2_n = [float(i) for i in co2_n]
        self.co2_scale = co2_scale
        
        self.ch4 = [float(i) for i in ch4]
        self.ch4sd = [float(i) for i in ch4sd]
        self.ch4_n = [float(i) for i in ch4_n]
        self.ch4_scale = ch4_scale
        
        self.co = [float(i) for i in co]
        self.cosd = [float(i) for i in cosd]
        self.co_n = [float(i) for i in co_n]
        self.co_scale = co_scale         
         
        self.filename = filename
        self.site = [i.split('_')[0] for i in sitekey]
        self.instrument = ['_'.join(i.split('_')[1:]) for i in sitekey]
        self.processeddate = GCWerks_outputdate
        self.USN = USN_tank
        self.DNo = DNo_tank
       
       
       
# Code to read in the files printed by Printmeans
class Read_n2o_means:
    def __init__(self, infile, DNo_tank, USN_tank):
       
        if type(infile) == tuple:
            data=np.genfromtxt(infile[0], dtype=str, skip_header=3, delimiter = ',')
        elif type(infile) == str:
            data=np.genfromtxt(infile, dtype=str, skip_header=3, delimiter = ',')     
        
        print 'reading file : '
        print infile        
        
        
        
        if len(data) != 11: 
            sitekey = data[:,0]
            date_time = data[:,1]
            n2o = data[:,2]
            n2osd = data[:,3]
            n2o_n = data[:,4]
            n2o_scale = data[:,5]
            sf6 = data[:,6]
            sf6sd = data[:,7]
            sf6_n = data[:,8]
            sf6_scale = data[:,9]
            filename = data[:,10]
            
            # make datetime variable
            dt_date = [dt.datetime.strptime(date_time[i] , "%y%m%d %H%M%S") for i in np.arange(len(date_time))]
        
            self.n2o = [float(i) for i in n2o]
            self.n2osd = [float(i) for i in n2osd]
            self.n2o_n = [float(i) for i in n2o_n]        
            
            self.sf6 = [float(i) for i in sf6]
            self.sf6sd = [float(i) for i in sf6sd]
            self.sf6_n = [float(i) for i in sf6_n]            
            
            self.site = [i.split('_')[0] for i in sitekey]
            self.instrument = ['_'.join(i.split('_')[1:]) for i in sitekey]
                
            
        else:
            sitekey = data[0]

            date_time = data[1]
            n2o = data[2]
            n2osd = data[3]
            n2o_n = data[4]
            n2o_scale = data[5]
            sf6 = data[6]
            sf6sd = data[7]
            sf6_n = data[8]
            sf6_scale = data[9]
            filename = data[10]
            
            # make datetime variable
            dt_date = [dt.datetime.strptime(date_time , "%y%m%d %H%M%S")] 
        
            self.n2o = [float(n2o)]
            self.n2osd = [float(n2osd)]
            self.n2o_n = [float(n2o_n)]
        
            self.sf6 = [float(sf6)]
            self.sf6sd = [float(sf6sd)]
            self.sf6_n = [float(sf6_n)]
              
           
        
            self.site = (data[0].split('_'))[0]
            self.instrument = '_'.join(data[0].split('_')[1:])
        
        
        self.sitekey = sitekey
        self.datetime = dt_date
    
        self.n2o_scale = n2o_scale
        

        self.sf6_scale = sf6_scale
                 
        self.filename = filename
        self.USN = USN_tank
        self.DNo = DNo_tank



            
# Code to plot up the output for multiple files        
class Plot_co2_means:
    def __init__(self, means, outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix=''):
        
        # The number of colours is the number of unique sites
        # one colour per site
        n_colours = len(set(means.site))
          
        symbol_fmt = ['o', 'v','s','D', 'x', '8', '*', '+']
        
        # Plot of 
        fig = plt.figure()
        fig.subplots_adjust(right = 0.6)
        fig.set_size_inches(8,4)
        
        legend_spacing = np.arange(len(means.site))   
        ls_1 = 0
        ls_2 = 0
        
        
        # Plot the WMO compatibility goals
        daterange = [min(means.datetime), max(means.datetime)]  

        
        co2_lines = WMOLines(means.co2, gas='CO2')
        
        plt1 = plt.subplot(3, 1, 1)
        plt1.plot(daterange, co2_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, co2_lines.minusline,  ':', color = 'grey')
    

        
        ch4_lines = WMOLines(means.ch4, gas = 'CH4')    
    
        plt2 = plt.subplot(3, 1, 2)
        plt2.plot(daterange, ch4_lines.plusline, ':', color = 'grey')
        plt2.plot(daterange, ch4_lines.minusline, ':', color = 'grey')
        
        
        co_lines = WMOLines(means.co, gas = 'CO')
        
        plt3 = plt.subplot(3, 1, 3)
        plt3.plot(daterange, co_lines.plusline, ':', color = 'grey')
        plt3.plot(daterange, co_lines.minusline, ':', color = 'grey')
        
        
        for i in np.arange(n_colours):      
            # extract the site name
            site = list(set(means.site))[i]
            
            # find all the results from that site
            siteindex = (np.where([m.startswith(site) for m in means.sitekey]))[0]
            
            site_sitekey = [means.sitekey[j] for j in siteindex]
            site_datetime = [means.datetime[j] for j in siteindex]            
            site_co2 = [means.co2[j] for j in siteindex]
            site_co2sd = [means.co2sd[j] for j in siteindex]
            site_ch4 = [means.ch4[j] for j in siteindex]
            site_ch4sd = [means.ch4sd[j] for j in siteindex]  
            site_co = [means.co[j] for j in siteindex]
            site_cosd = [means.cosd[j] for j in siteindex]         
            site_co_scale = [means.co_scale[j] for j in siteindex]     
            site_instrument = [means.instrument[j] for j in siteindex]   
            
            
            # check if they're from the same site            
            # one symbol type per instrument
            instruments = list(set(site_instrument))
            
            import re            
            
            #pdb.set_trace()
            
            # plot each instrument separately
            for h in np.arange(len(instruments)):
                
                B = np.asarray([re.match('\\b' + instruments[h]+ '\\b', l) for l in site_instrument])
                
                instrumentindex =[]
                
                            
                
                for k in np.arange(len(B)):
                    if B[k] != None:
                        instrumentindex.append(k)
                        
            
                plot_datetime = [site_datetime[l] for l in instrumentindex]            
                plot_co2 = [site_co2[l] for l in instrumentindex]    
                plot_co2sd = [site_co2sd[l] for l in instrumentindex]    
                
                plot_ch4 = [site_ch4[l] for l in instrumentindex]    
                plot_ch4sd = [site_ch4sd[l] for l in instrumentindex]    
                plot_co = [site_co[l] for l in instrumentindex]    
                plot_cosd = [site_cosd[l] for l in instrumentindex]    
                plot_co_scale = [site_co_scale[l] for l in instrumentindex]    
                                
                plot_sitekey = [site_sitekey[l] for l in instrumentindex]
              
                
                #pdb.set_trace()
                # check if they're all from the same instrument
                plt1 = plt.subplot(3, 1, 1)
                
                daterange = (dt.date(2014,6,1),dt.date(2015,5,1))
                if sum(np.array(means.site) == 'MPI')> 0:                
                    daterange = (dt.date(2013,11,1),dt.date(2015,5,1))
                
                
                
                plt1.errorbar(plot_datetime, plot_co2, yerr=plot_co2sd, fmt=symbol_fmt[h], color = SiteColours(sitein=site).outcolour, label = plot_sitekey, markersize = 3)
                #plt1.set_ylim(Calcrange(means.co2, means.co2sd))
                plt1.set_xlim(daterange)     
                
                y_formatter = ticker.ScalarFormatter(useOffset=False)
                plt1.yaxis.set_major_formatter(y_formatter)
                
                x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
                x_tickno_formatter = ticker.MaxNLocator(5)
                
                plt1.xaxis.set_major_formatter(x_formatter)
                plt1.xaxis.set_major_locator(x_tickno_formatter)
                
                plt1.set_title('Means - '+ means.USN)
                plt1.set_ylabel('[CO$_2$] (ppm)')
            
                plt2 = plt.subplot(3, 1, 2)
                plt2.errorbar(plot_datetime, plot_ch4, yerr=plot_ch4sd, fmt=symbol_fmt[h], color = SiteColours(sitein=site).outcolour, label = plot_sitekey, markersize = 3)
                plt2.set_ylim(Calcrange(means.ch4, means.ch4sd))
                plt2.set_ylabel('[CH$_4$] (ppb)')
                plt2.set_xlim(daterange)  
                plt2.yaxis.set_major_formatter(y_formatter)
                plt2.xaxis.set_major_formatter(x_formatter)            
                plt2.xaxis.set_major_locator(x_tickno_formatter)
                
                #
                plt3 = plt.subplot(3, 1, 3)
                plt3.errorbar(plot_datetime, plot_co, yerr=plot_cosd, fmt=symbol_fmt[h], color = SiteColours(sitein=site).outcolour, label = plot_sitekey, markersize = 3)
                plt3.set_ylim(Calcrange(means.co, means.cosd))
                plt3.set_ylabel('[CO] (ppb)')
                plt3.set_xlabel('time')
                plt3.set_xlim(daterange)  
                plt3.yaxis.set_major_formatter(y_formatter)
                plt3.xaxis.set_major_formatter(x_formatter)
                plt3.xaxis.set_major_locator(x_tickno_formatter)
                
                plt.figtext(0.64, 0.85-(0.05*legend_spacing[ls_1]), symbol_fmt[h]+' - '+ plot_sitekey[0], verticalalignment='bottom', \
                    horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=8)
                
                # Plot the CO data if it exists
                if np.isnan(plot_co[0]):
                    print 'No CO'
                else:          
                    plt.figtext(0.78, 0.85-(0.03*legend_spacing[ls_2]), plot_sitekey[0] +  ' - CO '+ plot_co_scale[0], verticalalignment='bottom',\
                        horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=6)
                    ls_2 = ls_2 + 1
                
                ls_1 =ls_1 +1   

                
        plt.figtext(0.8, 0.2, symbol_fmt[h]+  ' - CO2 WMOx2007', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.8, 0.17, symbol_fmt[h]+  ' - CH4 WMOx2004', verticalalignment='bottom', horizontalalignment='left', fontsize=6)

        plt.figtext(0.64, 0.05, 'Mean +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')

        print 'Plot saved as: ' + outputdir+'Means_'+means.USN+suffix+'.png'
        plt.savefig(outputdir+'Means_'+means.USN+suffix+'.png', dpi=200)
            
        plt.show()


        
 # Code to plot up the output for multiple files        
class Plot_n2o_means:
    def __init__(self, means, outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix=''):
        
        # The number of colours is the number of unique sites
        # one colour per site
        n_colours = len(set(means.site))
        
        symbol_fmt = ['o', 'v','s','D', 'x', '8', '*', '+']
        
        # Plot of 
        fig = plt.figure()
        fig.subplots_adjust(right = 0.7)
        fig.set_size_inches(8,4)
       
        legend_spacing = np.arange(len(means.site))   
        ls =0
        
        # Plot the WMO compatibility goals
        daterange = [min(means.datetime), max(means.datetime)]  
        
        n2o_lines = WMOLines(means.n2o, gas='N2O')
        
        plt1 = plt.subplot(2, 1, 1)
        plt1.plot(daterange, n2o_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, n2o_lines.minusline,  ':', color = 'grey')
    
    
        sf6_lines = WMOLines(means.sf6, gas = 'SF6')    
    
        plt2 = plt.subplot(2, 1, 2)
        plt2.plot(daterange, sf6_lines.plusline, ':', color = 'grey')
        plt2.plot(daterange, sf6_lines.minusline, ':', color = 'grey')
        

        for i in np.arange(n_colours):      
            # extract the site name
            site = list(set(means.site))[i]
            
            # find all the results from that site
            siteindex = (np.where([m.startswith(site) for m in means.sitekey]))[0]          
              
            #pdb.set_trace()
             
            if len(siteindex) > 1:  
                site_sitekey = [means.sitekey[j] for j in siteindex]
                site_datetime = [means.datetime[j] for j in siteindex]            
                site_n2o = [means.n2o[j] for j in siteindex]
                site_n2osd = [means.n2osd[j] for j in siteindex]
                site_n2o_scale = [means.n2o_scale[j] for [j] in siteindex]
                site_sf6 = [means.sf6[j] for j in siteindex]
                site_sf6sd = [means.sf6sd[j] for j in siteindex]  
                site_sf6_scale = [means.sf6_scale[j] for [j] in siteindex] 
                site_instrument = [means.instrument[j] for j in siteindex]   
            
            else:
                site_sitekey = means.sitekey[siteindex]
                site_datetime = means.datetime[siteindex]    
                site_n2o = means.n2o[siteindex]
                site_n2osd = means.n2osd[siteindex]
                site_n2o_scale = means.n2o_scale[siteindex]
                site_sf6 = means.sf6[siteindex]
                site_sf6sd = means.sf6sd[siteindex]  
                site_sf6_scale = means.sf6_scale[siteindex] 
                site_instrument = [means.instrument[siteindex]]   
 
            # check if they're from the same site            
            # one symbol type per instrument
            instruments = list(set(site_instrument))
            
            import re            
            
            daterange = (dt.date(2014,6,1),dt.date(2015,5,1))
            
            if sum(np.array(means.site) == 'MPI') > 0:
                daterange = (dt.date(2013,11,1),dt.date(2015,5,1))

            
            # plot each instrument separately
            for h in np.arange(len(instruments)):
                
                B = np.asarray([re.match('\\b' + instruments[h]+ '\\b', l) for l in site_instrument])
                
                instrumentindex =[]
                
                            
                
                for k in np.arange(len(B)):
                    if B[k] != None:
                        instrumentindex.append(k)
                        
               
                
                if instrumentindex == 0 : 
                    plot_datetime = [site_datetime[l] for l in instrumentindex]            
                    plot_n2o = [site_n2o[l] for l in instrumentindex]    
                    plot_n2osd = [site_n2osd[l] for l in instrumentindex]    
                    plot_sf6 = [site_sf6[l] for l in instrumentindex]    
                    plot_sf6sd = [site_sf6sd[l] for l in instrumentindex]    
                    
                    plot_sitekey = [site_sitekey[l] for l in instrumentindex]
                    plot_n2o_scale = [site_n2o_scale[l] for l in instrumentindex]
                    plot_sf6_scale = [site_sf6_scale[l] for l in instrumentindex]
                
                else:
                    plot_datetime = site_datetime    
                    plot_n2o = site_n2o
                    plot_n2osd = site_n2osd
                    plot_sf6 = site_sf6
                    plot_sf6sd = site_sf6sd
                    
                    plot_sitekey = site_sitekey
                    plot_n2o_scale = site_n2o_scale
                    plot_sf6_scale = site_sf6_scale

                
                #pdb.set_trace()
                # check if they're all from the same instrument
                plt1 = plt.subplot(2, 1, 1)
                
        
                plt1.errorbar(plot_datetime, plot_n2o, yerr=plot_n2osd, fmt=symbol_fmt[h], color = SiteColours(sitein=site).outcolour, label = plot_sitekey, markersize = 3)
                plt1.set_ylim(Calcrange(means.n2o, means.n2osd))
                plt1.set_xlim(daterange)     
                
                y_formatter = ticker.ScalarFormatter(useOffset=False)
                plt1.yaxis.set_major_formatter(y_formatter)
                
                
                x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
                x_tickno_formatter = ticker.MaxNLocator(5)
                plt1.xaxis.set_major_formatter(x_formatter)
                plt1.xaxis.set_major_locator(x_tickno_formatter)
                
                plt1.set_title('Means - '+ means.USN)
                plt1.set_ylabel('[N$_2$O] (ppb)')
            
                plt2 = plt.subplot(2, 1, 2)
                plt2.errorbar(plot_datetime, plot_sf6, yerr=plot_sf6sd, fmt=symbol_fmt[h], color =  SiteColours(sitein=site).outcolour, label = plot_sitekey, markersize = 3)
                plt2.set_ylim(Calcrange(means.sf6, means.sf6sd))
                plt2.set_ylabel('[SF$_6$] (ppt)')
                plt2.set_xlim(daterange)  
                plt2.yaxis.set_major_formatter(y_formatter)
                plt2.xaxis.set_major_formatter(x_formatter)            
                plt2.xaxis.set_major_locator(x_tickno_formatter)
                
                # Symbol and site
                plt.figtext(0.72, 0.85-(0.05*legend_spacing[ls]), plot_sitekey, verticalalignment='bottom', horizontalalignment='left', color= SiteColours(sitein=site).outcolour, fontsize=8)
                
                # Scale                
                plt.figtext(0.72, 0.4-(0.03*legend_spacing[ls]), plot_n2o_scale+ ', ' + plot_sf6_scale, verticalalignment='bottom', horizontalalignment='left', color= SiteColours(sitein=site).outcolour, fontsize=6)
                
            
                #pdb.set_trace()
                
                ls =ls +1               
                
        plt.figtext(0.72, 0.42, 'N$_2$O scale, SF$_6$ scale', verticalalignment='bottom', horizontalalignment='left', color='black', fontsize=6)
        plt.figtext(0.72, 0.1, 'Mean +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color ='grey')


        plt.savefig(outputdir+'N2O_SF6_Means_'+means.USN+suffix+'.png', dpi=200)
        
        print "Figure saved as:"
        print outputdir+'N2O_SF6_Means_'+means.USN+suffix+'.png'
        
        plt.show()
       
        
# Code to read in the files printed by Printmeans
class Read_scales:
    def __init__(self, infile='/Users/as13988/Documents/Work/Cylinders/ICP/Results/RawData/CalibrationScales_20140129'):
       
        if type(infile) == tuple:
            data=np.genfromtxt(infile[0], dtype=str, skip_header=2, delimiter = ', ')
        elif type(infile) == str:
            data=np.genfromtxt(infile, dtype=str, skip_header=2, delimiter = ', ')     
        
        sitekey = data[:,0]
        
        co2_scale = data[:,1]      
        ch4_scale = data[:,2]
        co_scale = data[:,3]
        n2o_scale = data[:,4]
        sf6_scale = data[:,5]
        
        self.sitekey = sitekey
        self.co2_scale = co2_scale
        self.ch4_scale = ch4_scale
        self.co_scale = co_scale   
        self.n2o_scale = n2o_scale
        self.sf6_scale = sf6_scale
         
        self.filename = infile
       

# Make a summary plot with mean tank vlaues removed  
class Plot_summary:
    def __init__(self, N2O_D091968, N2O_D091969, N2O_D091970, CO2_D091968, CO2_D091969, CO2_D091970, \
     outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix=''):
       
        # DNos = ['D091968', 'D091969', 'D091970'] 
       
        sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site, CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        sites = list(chain(*sites))
        
        co2_datetimes = [CO2_D091968.datetime, CO2_D091969.datetime, CO2_D091970.datetime]
        co2_datetimes = list(chain(*co2_datetimes))
        
        n2o_datetimes = [N2O_D091968.datetime, N2O_D091969.datetime, N2O_D091970.datetime]
        n2o_datetimes = list(chain(*n2o_datetimes))
        
        
        
        co2_sites = [CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        co2_sites = list(chain(*co2_sites))
        
        co2 = [CO2_D091968.co2, CO2_D091969.co2, CO2_D091970.co2]
        co2 = list(chain(*co2))
        ch4 = [CO2_D091968.ch4, CO2_D091969.ch4, CO2_D091970.ch4]
        ch4 = list(chain(*ch4))
        
        co2_diffs = [CO2_D091968.co2_diffs, CO2_D091969.co2_diffs, CO2_D091970.co2_diffs]
        co2_diffs = list(chain(*co2_diffs))
        ch4_diffs = [CO2_D091968.ch4_diffs, CO2_D091969.ch4_diffs, CO2_D091970.ch4_diffs]
        ch4_diffs = list(chain(*ch4_diffs))
        
        
        n2o_sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site]
        n2o_sites = list(chain(*n2o_sites))
        
        n2o = [N2O_D091968.n2o, N2O_D091969.n2o, N2O_D091970.n2o]
        n2o = list(chain(*n2o))
        
        n2o_diffs = [N2O_D091968.n2o_diffs, N2O_D091969.n2o_diffs, N2O_D091970.n2o_diffs]
        n2o_diffs = list(chain(*n2o_diffs))
        
        
        co2sd = [CO2_D091968.co2sd, CO2_D091969.co2sd, CO2_D091970.co2sd]
        co2sd = list(chain(*co2sd))
        ch4sd = [CO2_D091968.ch4sd, CO2_D091969.ch4sd, CO2_D091970.ch4sd]
        ch4sd = list(chain(*ch4sd))
        n2osd = [N2O_D091968.n2osd, N2O_D091969.n2osd, N2O_D091970.n2osd]
        n2osd = list(chain(*n2osd))
       
        # The number of colours is the number of unique sites
        # one colour per site
        n_colours = len(set(sites))
          
               
        # Plot of 
        fig = plt.figure()
        fig.subplots_adjust(right = 0.7)
        fig.set_size_inches(8,4)

        legend_spacing = np.arange(len(sites))   
        ls_1 = 0
        
        # Plot the WMO compatibility goals
        daterange = [min(co2_datetimes), max(co2_datetimes)]  

        
        co2_lines = WMOLines(co2_diffs, gas='CO2')
        
        plt1 = plt.subplot(3, 1, 1)
        plt1.plot(daterange, co2_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, co2_lines.minusline,  ':', color = 'grey')
    
    
        ch4_lines = WMOLines(ch4_diffs, gas = 'CH4')    
    
        plt2 = plt.subplot(3, 1, 2)
        plt2.plot(daterange, ch4_lines.plusline, ':', color = 'grey')
        plt2.plot(daterange, ch4_lines.minusline, ':', color = 'grey')
        
        
        n2o_lines = WMOLines(n2o_diffs, gas = 'N2O')
        
        plt3 = plt.subplot(3, 1, 3)
        plt3.plot(daterange, n2o_lines.plusline, ':', color = 'grey')
        plt3.plot(daterange, n2o_lines.minusline, ':', color = 'grey')
        
        
        
        
        for i in np.arange(n_colours):      
            # extract the site name
            site = list(set(sites))[i]
            
            # find all the results from that site
            co2_siteindex = (np.where(np.array(co2_sites) == site))[0]
            n2o_siteindex = (np.where(np.array(n2o_sites) == site))[0]
            
            co2_site_datetime = [co2_datetimes[j] for j in co2_siteindex]     
            n2o_site_datetime = [n2o_datetimes[j] for j in n2o_siteindex]     
            
            #pdb.set_trace()
            site_co2_diffs = [co2_diffs[j] for j in co2_siteindex]
            site_co2sd = [co2sd[j] for j in co2_siteindex]
            site_ch4_diffs = [ch4_diffs[j] for j in co2_siteindex]
            site_ch4sd = [ch4sd[j] for j in co2_siteindex]  
            site_n2o_diffs = [n2o_diffs[j] for j in n2o_siteindex]
            site_n2osd = [n2osd[j] for j in n2o_siteindex]         
        
            
            
            # check if they're all from the same instrument
            plt1 = plt.subplot(3, 1, 1)
            
            daterange = (dt.date(2014,6,1),dt.date(2015,5,1))
            
            if sum(np.array(sites) == 'MPI') > 0:            
                daterange = (dt.date(2013,11,1),dt.date(2015,5,1))
            
            plt1.errorbar(co2_site_datetime, site_co2_diffs, yerr=site_co2sd, fmt='o', color = SiteColours(sitein=site).outcolour, markersize = 3)
            plt1.set_ylim(Calcrange(co2_diffs, co2sd, even=1))
            plt1.set_xlim(daterange)     
            
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
            
            x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
            x_tickno_formatter = ticker.MaxNLocator(5)
            
            plt1.xaxis.set_major_formatter(x_formatter)
            plt1.xaxis.set_major_locator(x_tickno_formatter)
            
            plt1.set_title('Cylinder intercomparison', fontsize =7)
            plt1.set_ylabel('mean - [CO$\mathregular{_2}$] (ppm)', fontsize =7)
        
            plt2 = plt.subplot(3, 1, 2)
            plt2.errorbar(co2_site_datetime, site_ch4_diffs, yerr=site_ch4sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt2.set_ylim(Calcrange(ch4_diffs, ch4sd, even=1))
            plt2.set_ylabel('mean- [CH$\mathregular{_4}$] (ppb)', fontsize =7)
            
            plt2.set_xlim(daterange)  
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.xaxis.set_major_formatter(x_formatter)            
            plt2.xaxis.set_major_locator(x_tickno_formatter)
            
            
            plt3 = plt.subplot(3, 1, 3)
            plt3.errorbar(n2o_site_datetime, site_n2o_diffs, yerr=site_n2osd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt3.set_ylim(Calcrange(n2o_diffs, n2osd, even=1))
            plt3.set_ylabel('mean - [N$\mathregular{_2}$O] (ppb)', fontsize =7)
            plt3.set_xlabel('time')
            plt3.set_xlim(daterange)  
            plt3.yaxis.set_major_formatter(y_formatter)
            plt3.xaxis.set_major_formatter(x_formatter)
            plt3.xaxis.set_major_locator(x_tickno_formatter)
            
            plt.figtext(0.74, 0.85-(0.05*legend_spacing[ls_1]), site, verticalalignment='bottom', horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=8)
            
            ls_1 =ls_1 +1   
                
        # '$\mathregular{x^2}$'
        
        plt.figtext(0.72, 0.2, 'CO$\mathregular{_2}$ - WMOx2007', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.17, 'CH$\mathregular{_4}$ - WMOx2004', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - SIO98 (DECC/GAUGE), WMOx2006 (others) ', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.11, 'Mean +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')



        plt.savefig(outputdir+'SummmaryPlot_Timeseries' + suffix + '.png', dpi=200)
            
        plt.show()        
        
        
        n_colours = len(set(sites))
              
        # Plot of diffs against concentration
        fig = plt.figure()
        fig.subplots_adjust(right = 0.7)
        fig.set_size_inches(8,4)
        
        legend_spacing = np.arange(len(sites))   
        ls_1 = 0
        
        # Plot WMO compatibility goals    
        co2_range = [min(co2),max(co2)]
        
        plt1 = plt.subplot(3, 1, 1)
        plt1.plot(co2_range, co2_lines.plusline,  ':', color = 'grey')
        plt1.plot(co2_range, co2_lines.minusline,  ':', color = 'grey')
    
    
        ch4_range = [min(ch4),max(ch4)]
        
        plt2 = plt.subplot(3, 1, 2)
        plt2.plot(ch4_range, ch4_lines.plusline, ':', color = 'grey')
        plt2.plot(ch4_range, ch4_lines.minusline, ':', color = 'grey')
        
        
        n2o_range = [min(n2o),max(n2o)]
        
        plt3 = plt.subplot(3, 1, 3)
        plt3.plot(n2o_range, n2o_lines.plusline, ':', color = 'grey')
        plt3.plot(n2o_range, n2o_lines.minusline, ':', color = 'grey')

        
        for k in np.arange(n_colours):      
            # extract the site name
            site = list(set(sites))[k]
            
            # find all the results from that site
            co2_siteindex = (np.where(np.array(co2_sites) == site))[0]
            n2o_siteindex = (np.where(np.array(n2o_sites) == site))[0]
            
            #pdb.set_trace()
            site_co2 = [co2[j] for j in co2_siteindex]
            site_co2sd = [co2sd[j] for j in co2_siteindex]
            site_ch4 = [ch4[j] for j in co2_siteindex]
            site_ch4sd = [ch4sd[j] for j in co2_siteindex]  
            site_n2o = [n2o[j] for j in n2o_siteindex]
            site_n2osd = [n2osd[j] for j in n2o_siteindex]         
            site_co2_diffs = [co2_diffs[j] for j in co2_siteindex]
            site_ch4_diffs = [ch4_diffs[j] for j in co2_siteindex]
            site_n2o_diffs = [n2o_diffs[j] for j in n2o_siteindex]
   
            
            
            # check if they're all from the same instrument
            plt1 = plt.subplot(3, 1, 1)
            
            plt1.errorbar(site_co2, site_co2_diffs, yerr=site_co2sd, fmt='o', color = SiteColours(sitein=site).outcolour, markersize = 3)
            plt1.set_ylim(Calcrange(co2_diffs, co2sd, even=1))
            plt1.set_xlim([370,490])
            
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
            
            
            plt1.set_title('Cylinder intercomparison', fontsize = 7)
            plt1.set_ylabel('mean - [CO$\mathregular{_2}$] (ppm)', fontsize = 7)
       
        
            plt2 = plt.subplot(3, 1, 2)
            plt2.errorbar(site_ch4, site_ch4_diffs, yerr=site_ch4sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt2.set_ylim(Calcrange(ch4_diffs, ch4sd, even = 1))
            plt2.set_ylabel('mean - [CH$\mathregular{_4}$] (ppb)', fontsize = 7)
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.set_xlim([1750,2250])
            
            plt3 = plt.subplot(3, 1, 3)
            plt3.errorbar(site_n2o, site_n2o_diffs, yerr=site_n2osd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt3.set_ylim(Calcrange(n2o_diffs, n2osd, even=1))
            plt3.set_ylabel('mean - [N$\mathregular{_2}$O] (ppb)', fontsize = 7)
            plt3.set_xlim([322,336])

            plt3.set_xlabel('Concentration', fontsize = 7)
            plt3.yaxis.set_major_formatter(y_formatter)
            
            plt.figtext(0.74, 0.85-(0.05*legend_spacing[ls_1]), site, verticalalignment='bottom', horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=8)
            
            ls_1 =ls_1 +1   
                
        # '$\mathregular{x^2}$'
        
        plt.figtext(0.72, 0.2, 'CO$\mathregular{_2}$ - WMOx2007', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.17, 'CH$\mathregular{_4}$ - WMOx2004', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - SIO98 (DECC/GAUGE), WMOx2006 (others) ', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.11, 'Mean +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')


        plt.savefig(outputdir+'SummmaryPlot_Conc' + suffix + '.png', dpi=200)
            
        plt.show()        
        
        