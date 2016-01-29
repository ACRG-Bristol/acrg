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
reload(read_GCwerks)

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
        colours = ['blue', 'red', 'orange', 'green', 'dodgerblue', 'gold', 'orchid', 'purple', 'black', 'sienna', 'limegreen', 'aqua', 'coral', 'deepskyblue']
        sites = ['RGL', 'TAC', 'HFD', 'BSD', 'TTA', 'FAAM', 'Man', 'BAS', 'MPI', 'FERRY', 'BT', 'HAD', 'TIL', 'WAO']
        
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
        
        gases = ['CO2', 'CH4', 'CO', 'N2O', 'SF6', 'co2', 'ch4', 'co', 'n2o', 'sf6']
        compatibility = [0.1, 2, 2, 0.1, 0.02, 0.1, 2, 2, 0.1, 0.02]
        units = ['ppm', 'ppb','ppb','ppb','ppt', 'ppm', 'ppb','ppb','ppb','ppt']
    
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
        self.h2o_reported = indata.h2o_reported
        
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
        
        if indata.nogases == 3:
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
        
        co_scale = data.co_scale[0]
    
    co2_scale = data.co2_scale[0]
    ch4_scale = data.ch4_scale[0]
    
    
    
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
        sitekey_list = list(chain(*sitekey))
        co2_list = (list(chain(*co2)))
        ch4_list = (list(chain(*ch4)))
        co_list = (list(chain(*co)))
        
        datetimes = []        
        
        for i in np.arange(len(datetime)):
            for j in np.arange(len(datetime[i])):
                datetimes.append(datetime[i][j])        
        
        
        print co2_list[np.where([i == 'MPI_Picarro' for i in sitekey_list])[0]]
        print type(co2_list[np.where([i == 'MPI_Picarro' for i in sitekey_list])[0]])
        MPI_co2 = co2_list[np.where([i == 'MPI_Picarro' for i in sitekey_list])[0]]
        MPI_ch4 = ch4_list[np.where([i == 'MPI_Picarro' for i in sitekey_list])[0]]
        MPI_co = co_list[np.where([i == 'MPI_Aerolaser' for i in sitekey_list])[0]]
        EMPA_co2 = co2_list[np.where([i == 'EMPA_CRDS_Prelim' for i in sitekey_list])[0]]
        EMPA_ch4 = ch4_list[np.where([i == 'EMPA_CRDS_Prelim' for i in sitekey_list])[0]]
        EMPA_co = co_list[np.where([i == 'EMPA_CRDS_Prelim' for i in sitekey_list])[0]]


        self.MPI_co2 = MPI_co2
        self.MPI_ch4 = MPI_ch4
        self.MPI_co = MPI_co
        self.EMPA_co2 = EMPA_co2
        self.EMPA_ch4 = EMPA_ch4
        self.EMPA_co = EMPA_co
        
        self.datetime = datetimes
        self.co2 = list(chain(*co2))
        self.co2sd = list(chain(*co2sd))
        self.co2_n = list(chain(*co2_n))
        self.co2_scale = list(chain(*co2_scale))
        self.co2_diffs = list(chain(*co2)) - np.nanmean(list(chain(*co2)))
        self.co2_MPIdiff = [i - MPI_co2 for i in co2_list]
        self.co2_EMPAdiff = [i - EMPA_co2 for i in co2_list]

        self.ch4 = list(chain(*ch4))
        self.ch4sd = list(chain(*ch4sd))
        self.ch4_n = list(chain(*ch4_n))
        self.ch4_scale = list(chain(*ch4_scale)) 
        self.ch4_diffs = list(chain(*ch4)) - np.nanmean(list(chain(*ch4)))
        self.ch4_MPIdiff = [i - MPI_ch4 for i in ch4_list]
        self.ch4_EMPAdiff = [i - EMPA_ch4 for i in ch4_list]
       
        self.co = list(chain(*co))
        self.cosd = list(chain(*cosd))
        self.co_n = list(chain(*co_n))
        self.co_scale = list(chain(*co_scale))
        self.co_diffs = list(chain(*co)) - np.nanmean(list(chain(*co)))
        self.co_MPIdiff = [i - MPI_co for i in co_list]
        self.co_EMPAdiff = [i - EMPA_co for i in co_list]

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
        
        n2o_list = (list(chain(*n2o)))         
        sf6_list = (list(chain(*sf6)))
        
        MPI_n2o = n2o_list[np.where([i == 'MPI_GC' for i in sitekey])[0]]
        MPI_sf6 = sf6_list[np.where([i == 'MPI_GC' for i in sitekey])[0]]
        EMPA_n2o = n2o_list[np.where([i == 'EMPA_CRDS_Prelim' for i in sitekey])[0]]
        EMPA_sf6 = sf6_list[np.where([i == 'EMPA_CRDS_Prelim' for i in sitekey])[0]]
        
        self.MPI_n2o = MPI_n2o
        self.MPI_sf6 = MPI_sf6
        self.EMPA_n2o = EMPA_n2o
        self.EMPA_sf6 = EMPA_sf6
        
        self.n2o = list(chain(*n2o))
        self.n2osd = list(chain(*n2osd))
        self.n2o_n = list(chain(*n2o_n))
        self.n2o_scale = n2o_scale
        self.n2o_diffs = list(chain(*n2o)) - np.nanmean(list(chain(*n2o)))
        self.n2o_MPIdiff = [ i - MPI_n2o for i in n2o_list]
        self.n2o_EMPAdiff = [ i - EMPA_n2o for i in n2o_list]
        
        self.sf6 = list(chain(*sf6))
        self.sf6sd = list(chain(*sf6sd))
        self.sf6_n = list(chain(*sf6_n))
        self.sf6_scale = sf6_scale
        self.sf6_diffs = list(chain(*sf6)) - np.nanmean(list(chain(*sf6)))
        self.sf6_MPIdiff = [i - MPI_sf6 for i in sf6_list]
        self.sf6_EMPAdiff = [i - EMPA_sf6 for i in sf6_list]
        
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
# Set useMPI flag if you want to just compare to MPI number 
class Plot_co2_means:
    def __init__(self, means, outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix='', useMPI = 0):
        
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
        
        if useMPI == 1:
            co2_lines = WMOLines(means.MPI_co2, gas='CO2')
                        
        
        plt1 = plt.subplot(3, 1, 1)
        plt1.plot(daterange, co2_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, co2_lines.minusline,  ':', color = 'grey')
    

        
        ch4_lines = WMOLines(means.ch4, gas = 'CH4')    

        if useMPI == 1:
            ch4_lines = WMOLines(means.MPI_ch4, gas='CH4')

    
        plt2 = plt.subplot(3, 1, 2)
        plt2.plot(daterange, ch4_lines.plusline, ':', color = 'grey')
        plt2.plot(daterange, ch4_lines.minusline, ':', color = 'grey')
        
        
        
        if useMPI == 1:
            co_lines = WMOLines(means.MPI_co, gas='CO')
        else:
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
                
                daterange = (dt.date(2014,6,1),dt.date(2015,10,1))
                if sum(np.array(means.site) == 'MPI')> 0:                
                    daterange = (dt.date(2013,11,1),dt.date(2015,10,1))
                
                
                
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
                
                plt.figtext(0.64, 0.85-(0.03*legend_spacing[ls_1]), symbol_fmt[h]+' - '+ plot_sitekey[0], verticalalignment='bottom', \
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

        if useMPI == 1:
            plt.figtext(0.64, 0.05, 'MPI +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')
        else:
            plt.figtext(0.64, 0.05, 'Mean +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')

        print 'Plot saved as: ' + outputdir+'Means_'+means.USN+suffix+'.png'
        plt.savefig(outputdir+'Means_'+means.USN+suffix+'.png', dpi=200)
            
        plt.show()


        
 # Code to plot up the output for multiple files       
 # Set useMPI flag to use MPI rather than mean
class Plot_n2o_means:
    def __init__(self, means, outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix='', useMPI=0):
        
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
        
        if useMPI == 1:        
            n2o_lines = WMOLines(means.MPI_n2o, gas='N2O')
            
        plt1 = plt.subplot(2, 1, 1)
        plt1.plot(daterange, n2o_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, n2o_lines.minusline,  ':', color = 'grey')
    
    
        sf6_lines = WMOLines(means.sf6, gas = 'SF6')

        if useMPI == 1:        
            sf6_lines = WMOLines(means.MPI_sf6, gas='SF6')

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
            
            daterange = (dt.date(2014,6,1),dt.date(2015,10,1))
            
            if sum(np.array(means.site) == 'MPI') > 0:
                daterange = (dt.date(2013,11,1),dt.date(2015,10,1))

            
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

        if useMPI ==1 :
            plt.figtext(0.72, 0.1, 'MPI +/- WMO compatibility goals', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color ='grey')
        else:
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
       

# Make a summary plot with mean tank values removed  
# Set use_MPI if you want to plot relative to MPI
class Plot_summary:
    def __init__(self, N2O_D091968, N2O_D091969, N2O_D091970, CO2_D091968, CO2_D091969, CO2_D091970, \
     outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix='', use_MPI = 0):
       
        # DNos = ['D091968', 'D091969', 'D091970'] 
       
        sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site, CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        sites = list(chain(*sites))
        
        co2_datetimes = [CO2_D091968.datetime, CO2_D091969.datetime, CO2_D091970.datetime]
        co2_datetimes = list(chain(*co2_datetimes))
        
        n2o_datetimes = [N2O_D091968.datetime, N2O_D091969.datetime, N2O_D091970.datetime]
        n2o_datetimes = list(chain(*n2o_datetimes))
        
        co2_sites = [CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        co2_sites = list(chain(*co2_sites))
        
        co2_sitekey = [CO2_D091968.sitekey, CO2_D091969.sitekey, CO2_D091970.sitekey]
        co2_sitekey = list(chain(*co2_sitekey))
        
        co2 = [CO2_D091968.co2, CO2_D091969.co2, CO2_D091970.co2]
        co2 = list(chain(*co2))
        ch4 = [CO2_D091968.ch4, CO2_D091969.ch4, CO2_D091970.ch4]
        ch4 = list(chain(*ch4))
        
        
        
        co2_diffs = [CO2_D091968.co2_diffs, CO2_D091969.co2_diffs, CO2_D091970.co2_diffs]
        co2_ylabel = 'Site - mean [CO$\mathregular{_2}$] (ppm)'
        
        if use_MPI == 1:
            co2_ylabel = 'Site - MPI [CO$\mathregular{_2}$] (ppm)'
            co2_diffs = [CO2_D091968.co2_MPIdiff, CO2_D091969.co2_MPIdiff, CO2_D091970.co2_MPIdiff]
            
        co2_diffs = list(chain(*co2_diffs))
                    
        
        ch4_diffs = [CO2_D091968.ch4_diffs, CO2_D091969.ch4_diffs, CO2_D091970.ch4_diffs]
        ch4_ylabel = 'Site - mean [CH$\mathregular{_4}$] (ppb)'
        
        if use_MPI == 1:
            ch4_diffs = [CO2_D091968.ch4_MPIdiff, CO2_D091969.ch4_MPIdiff, CO2_D091970.ch4_MPIdiff]
            ch4_ylabel = 'Site - MPI [CH$\mathregular{_4}$] (ppb)'

        ch4_diffs = list(chain(*ch4_diffs))
        
        
        n2o_sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site]
        n2o_sites = list(chain(*n2o_sites))
        
        n2o_sitekey = [N2O_D091968.sitekey, N2O_D091969.sitekey, N2O_D091970.sitekey]
        n2o_sitekey = list(chain(*n2o_sitekey))

        n2o = [N2O_D091968.n2o, N2O_D091969.n2o, N2O_D091970.n2o]
        n2o = list(chain(*n2o))
        
        n2o_diffs = [N2O_D091968.n2o_diffs, N2O_D091969.n2o_diffs, N2O_D091970.n2o_diffs]
        n2o_ylabel = 'Site - mean [N$\mathregular{_2}$O] (ppb)'

        if use_MPI == 1:
            n2o_diffs = [N2O_D091968.n2o_MPIdiff, N2O_D091969.n2o_MPIdiff, N2O_D091970.n2o_MPIdiff]
            n2o_ylabel = 'Site - MPI [N$\mathregular{_2}$O] (ppb)'
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
        if use_MPI == 1:
            co2_lines = WMOLines(0, gas='CO2')
            
        
        plt1 = plt.subplot(3, 1, 1)
        plt1.plot(daterange, co2_lines.plusline,  ':', color = 'grey')
        plt1.plot(daterange, co2_lines.minusline,  ':', color = 'grey')
    
    
        ch4_lines = WMOLines(ch4_diffs, gas = 'CH4')    
        if use_MPI == 1:
            ch4_lines = WMOLines(0, gas='CH4')
            
       
       
        plt2 = plt.subplot(3, 1, 2)
        plt2.plot(daterange, ch4_lines.plusline, ':', color = 'grey')
        plt2.plot(daterange, ch4_lines.minusline, ':', color = 'grey')
        
        
        n2o_lines = WMOLines(n2o_diffs, gas = 'N2O')
        if use_MPI == 1:
            n2o_lines = WMOLines(0, gas='N2O')
            
       
        plt3 = plt.subplot(3, 1, 3)
        plt3.plot(daterange, n2o_lines.plusline, ':', color = 'grey')
        plt3.plot(daterange, n2o_lines.minusline, ':', color = 'grey')
        
        
        
        
        for i in np.arange(n_colours):      
            # extract the site name
            site = list(set(sites))[i]
            print site
            
            # find all the results from that site
            co2_siteindex = (np.where(np.array(co2_sites) == site))[0]
            n2o_siteindex = (np.where(np.array(n2o_sites) == site))[0]
            
            if len(co2_siteindex) != 0:
                co2_sitekey_tag = co2_sitekey[co2_siteindex[0]]
            else:
                co2_sitekey_tag=''
                
            if len(n2o_siteindex) != 0:
                n2o_sitekey_tag = n2o_sitekey[n2o_siteindex[0]]
            else:
                n2o_sitekey_tag=''
    
            sitekey = co2_sitekey_tag + ', ' + n2o_sitekey_tag            
            print sitekey
          
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
            
            daterange = (dt.date(2014,6,1),dt.date(2015,10,1))
            
            if sum(np.array(sites) == 'MPI') > 0:            
                daterange = (dt.date(2013,11,1),dt.date(2015,10,1))
            
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
            plt1.set_ylabel(co2_ylabel, fontsize =7)
        
            plt2 = plt.subplot(3, 1, 2)
            plt2.errorbar(co2_site_datetime, site_ch4_diffs, yerr=site_ch4sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt2.set_ylim(Calcrange(ch4_diffs, ch4sd, even=1))
            plt2.set_ylabel(ch4_ylabel, fontsize =7)
            
            plt2.set_xlim(daterange)  
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.xaxis.set_major_formatter(x_formatter)            
            plt2.xaxis.set_major_locator(x_tickno_formatter)
            
            
            plt3 = plt.subplot(3, 1, 3)
            plt3.errorbar(n2o_site_datetime, site_n2o_diffs, yerr=site_n2osd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt3.set_ylim(Calcrange(n2o_diffs, n2osd, even=1))
            plt3.set_ylabel(n2o_ylabel, fontsize =7)
            plt3.set_xlabel('time')
            plt3.set_xlim(daterange)  
            plt3.yaxis.set_major_formatter(y_formatter)
            plt3.xaxis.set_major_formatter(x_formatter)
            plt3.xaxis.set_major_locator(x_tickno_formatter)

            plt.figtext(0.74, 0.85-(0.04*legend_spacing[ls_1]), sitekey, verticalalignment='bottom', horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=8)
            
            ls_1 =ls_1 +1   
                
        # '$\mathregular{x^2}$'
        
        plt.figtext(0.72, 0.2, 'CO$\mathregular{_2}$ - WMOx2007', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.17, 'CH$\mathregular{_4}$ - WMOx2004, x2004A', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        #plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - SIO98 (DECC/GAUGE), WMOx2006A (others) ', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - WMOx2006A', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        
        text1 = 'Mean +/- WMO compatibility goals'
        if use_MPI == 1:
            text1 = 'MPI +/- WMO compatibility goals'
            
        if use_EMPA == 1:    
            text1 = 'EMPA +/- WMO compatibility goals'
        
        
        plt.figtext(0.72, 0.11, text1, verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')
        
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
            
            if len(co2_siteindex) != 0:
                co2_sitekey_tag = co2_sitekey[co2_siteindex[0]]
            else:
                co2_sitekey_tag=''
                
            if len(n2o_siteindex) != 0:
                n2o_sitekey_tag = n2o_sitekey[n2o_siteindex[0]]
            else:
                n2o_sitekey_tag=''
    
            sitekey = co2_sitekey_tag + ', ' + n2o_sitekey_tag            
            print sitekey
            
            
            
            
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
            plt1.set_ylabel(co2_ylabel, fontsize = 7)
       
        
            plt2 = plt.subplot(3, 1, 2)
            plt2.errorbar(site_ch4, site_ch4_diffs, yerr=site_ch4sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt2.set_ylim(Calcrange(ch4_diffs, ch4sd, even = 1))
            plt2.set_ylabel(ch4_ylabel, fontsize = 7)
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.set_xlim([1750,2250])
            
            plt3 = plt.subplot(3, 1, 3)
            plt3.errorbar(site_n2o, site_n2o_diffs, yerr=site_n2osd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
            plt3.set_ylim(Calcrange(n2o_diffs, n2osd, even=1))
            plt3.set_ylabel(n2o_ylabel, fontsize = 7)
            plt3.set_xlim([322,336])

            plt3.set_xlabel('Concentration', fontsize = 7)
            plt3.yaxis.set_major_formatter(y_formatter)
            
            plt.figtext(0.74, 0.85-(0.04*legend_spacing[ls_1]), sitekey, verticalalignment='bottom', horizontalalignment='left', color=SiteColours(sitein=site).outcolour, fontsize=8)
            
            ls_1 =ls_1 +1   
                
        # '$\mathregular{x^2}$'
        
        plt.figtext(0.72, 0.2, 'CO$\mathregular{_2}$ - WMOx2007', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.17, 'CH$\mathregular{_4}$ - WMOx2004, x2004A', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        #plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - SIO98 (DECC/GAUGE), WMOx2006A (others) ', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.72, 0.14, 'N$\mathregular{_2}$O - WMOx2006A (others) ', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        
        text1 = 'Mean +/- WMO compatibility goals'
        if use_MPI == 1:
            text1 = 'MPI +/- WMO compatibility goals'
            
        if use_EMPA == 1:
            text1 = 'EMPA +/- WMO compatibility goals'
        
        plt.figtext(0.72, 0.11, text1, verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')


        plt.savefig(outputdir+'SummmaryPlot_Conc' + suffix + '.png', dpi=200)
            
        plt.show()        

        
# Make a summary plot for a single site
import re
class Plot_individualsite:
    def __init__(self, site, N2O_D091968, N2O_D091969, N2O_D091970, CO2_D091968, CO2_D091969, CO2_D091970, \
     outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix='', use_MPI = 0):
       
        # DNos = ['D091968', 'D091969', 'D091970'] 
       
        sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site, CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        sites = list(chain(*sites))
        
        co2_datetimes = [CO2_D091968.datetime, CO2_D091969.datetime, CO2_D091970.datetime]
        co2_datetimes = list(chain(*co2_datetimes))
        
        n2o_datetimes = [N2O_D091968.datetime, N2O_D091969.datetime, N2O_D091970.datetime]
        n2o_datetimes = list(chain(*n2o_datetimes))
        
        # Extract the data from the CO2/CH4 files
        co2_sites = [CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        co2_sites = list(chain(*co2_sites))
        
        co2_inst_all = [CO2_D091968.instrument, CO2_D091969.instrument, CO2_D091970.instrument]
        co2_inst_all = list(chain(*co2_inst_all))        
        
        co2 = [CO2_D091968.co2, CO2_D091969.co2, CO2_D091970.co2]
        co2 = list(chain(*co2))
        ch4 = [CO2_D091968.ch4, CO2_D091969.ch4, CO2_D091970.ch4]
        ch4 = list(chain(*ch4))
        co = [CO2_D091968.co, CO2_D091969.co, CO2_D091970.co]
        co = list(chain(*co))
    
        co2sd = [CO2_D091968.co2sd, CO2_D091969.co2sd, CO2_D091970.co2sd]
        co2sd = list(chain(*co2sd))
        ch4sd = [CO2_D091968.ch4sd, CO2_D091969.ch4sd, CO2_D091970.ch4sd]
        ch4sd = list(chain(*ch4sd))
        cosd = [CO2_D091968.cosd, CO2_D091969.cosd, CO2_D091970.cosd]
        cosd = list(chain(*cosd))
       
        co2_diffs = [CO2_D091968.co2_diffs, CO2_D091969.co2_diffs, CO2_D091970.co2_diffs]
        co2_ylabel = 'Site - mean [CO$\mathregular{_2}$] (ppm)'
        if use_MPI == 1:
            co2_diffs = [CO2_D091968.co2_MPIdiff, CO2_D091969.co2_MPIdiff, CO2_D091970.co2_MPIdiff]
            co2_ylabel = 'Site - MPI [CO$\mathregular{_2}$] (ppm)'
        co2_diffs = list(chain(*co2_diffs))
        
        ch4_diffs = [CO2_D091968.ch4_diffs, CO2_D091969.ch4_diffs, CO2_D091970.ch4_diffs]        
        ch4_ylabel = 'Site - mean [CH$\mathregular{_4}$] (ppb)'
        if use_MPI == 1:
            ch4_diffs = [CO2_D091968.ch4_MPIdiff, CO2_D091969.ch4_MPIdiff, CO2_D091970.ch4_MPIdiff]
            ch4_ylabel = 'Site - MPI [CH$\mathregular{_4}$] (ppb)'
        ch4_diffs = list(chain(*ch4_diffs))
        
        co_diffs = [CO2_D091968.co_diffs, CO2_D091969.co_diffs, CO2_D091970.co_diffs]
        co_ylabel = 'Site - mean [CO] (ppb)'
        if use_MPI == 1:
            co_diffs = [CO2_D091968.co_MPIdiff, CO2_D091969.co_MPIdiff, CO2_D091970.co_MPIdiff]
            co_ylabel = 'Site - MPI [CO] (ppb)'
        co_diffs = list(chain(*co_diffs))
        
        
        
        # Extract the data from the N2O/SF6 files
        n2o_sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site]
        n2o_sites = list(chain(*n2o_sites))

        n2o_inst_all = [N2O_D091968.instrument, N2O_D091969.instrument, N2O_D091970.instrument]
        n2o_inst_all = list(chain(*n2o_inst_all))
       
        n2o = [N2O_D091968.n2o, N2O_D091969.n2o, N2O_D091970.n2o]
        n2o = list(chain(*n2o))
        sf6 = [N2O_D091968.sf6, N2O_D091969.sf6, N2O_D091970.sf6]
        sf6 = list(chain(*sf6))
        
        n2osd = [N2O_D091968.n2osd, N2O_D091969.n2osd, N2O_D091970.n2osd]
        n2osd = list(chain(*n2osd))        
        sf6sd = [N2O_D091968.sf6sd, N2O_D091969.sf6sd, N2O_D091970.sf6sd]
        sf6sd = list(chain(*sf6sd))        
        
        
        n2o_diffs = [N2O_D091968.n2o_diffs, N2O_D091969.n2o_diffs, N2O_D091970.n2o_diffs]
        n2o_ylabel = 'Site - mean [N$\mathregular{_2}$O] (ppb)'
        if use_MPI == 1:
            n2o_diffs = [N2O_D091968.n2o_MPIdiff, N2O_D091969.n2o_MPIdiff, N2O_D091970.n2o_MPIdiff]
            n2o_ylabel = 'Site - MPI [N$\mathregular{_2}$O] (ppb)'
        n2o_diffs = list(chain(*n2o_diffs))
        
        sf6_diffs = [N2O_D091968.sf6_diffs, N2O_D091969.sf6_diffs, N2O_D091970.sf6_diffs]
        sf6_ylabel = 'Site - mean [SF$\mathregular{_6}$] (ppt)'
        if use_MPI == 1:
            sf6_diffs = [N2O_D091968.sf6_MPIdiff, N2O_D091969.sf6_MPIdiff, N2O_D091970.sf6_MPIdiff]
            sf6_ylabel = 'Site - MPI [SF$\mathregular{_6}$] (ppt)'
        sf6_diffs = list(chain(*sf6_diffs))

        
             
        # Find where the site matches
        regex=re.compile(".*("+site+").*")
        
        # initialise index 
        siteindex = []
        for i, j in enumerate(sites):
            if regex.search(j):
                siteindex.append(i)
             
        # Plot of 
        fig = plt.figure()
        fig.subplots_adjust(wspace = 0.4)
        fig.set_size_inches(6,6)
        
         
        # find all the results from that site
        co2_siteindex = (np.where(np.array(co2_sites) == site))[0]
        n2o_siteindex = (np.where(np.array(n2o_sites) == site))[0]

        site_co2 = [co2[j] for j in co2_siteindex]
        site_co2sd = [co2sd[j] for j in co2_siteindex]
        site_ch4 = [ch4[j] for j in co2_siteindex]
        site_ch4sd = [ch4sd[j] for j in co2_siteindex]  
        site_n2o = [n2o[j] for j in n2o_siteindex]
        site_n2osd = [n2osd[j] for j in n2o_siteindex]         
        site_co = [co[j] for j in co2_siteindex]
        site_cosd = [cosd[j] for j in co2_siteindex]
        site_sf6 = [sf6[j] for j in n2o_siteindex]
        site_sf6sd = [sf6sd[j] for j in n2o_siteindex]  
        site_co2inst = [co2_inst_all[i] for i in co2_siteindex]        
        site_n2oinst = [n2o_inst_all[i] for i in n2o_siteindex]
        
        site_co2_diffs = [co2_diffs[j] for j in co2_siteindex]
        site_ch4_diffs = [ch4_diffs[j] for j in co2_siteindex]
        site_n2o_diffs = [n2o_diffs[j] for j in n2o_siteindex]
        site_co_diffs = [co_diffs[j] for j in co2_siteindex]
        site_sf6_diffs = [sf6_diffs[j] for j in n2o_siteindex]
  
        if len(n2o_siteindex)  <= 0:
            site_n2o = [np.nan, np.nan, np.nan]
            site_n2osd = [np.nan, np.nan, np.nan]    
            #site_co = [np.nan, np.nan, np.nan]
            #site_cosd = [np.nan, np.nan, np.nan]
            site_sf6 = [np.nan, np.nan, np.nan]
            site_sf6sd = [np.nan, np.nan, np.nan]
            site_n2o_diffs = [np.nan, np.nan, np.nan]
            #site_co_diffs = [np.nan, np.nan, np.nan]
            site_sf6_diffs = [np.nan, np.nan, np.nan]

        # Plot the WMO compatibility goals
        # For CO2
        co2_lines = WMOLines(co2_diffs, gas='CO2')
        if use_MPI == 1:
            co2_lines = WMOLines(0, gas='CO2')

        if np.isnan(site_co2).all():        
            inst_tag = ''
        else:
            inst_tag = site_co2inst[np.where(np.isfinite(site_co2_diffs))[0][0]]
             
            
        # Plot up the data
        plt1 = plt.subplot(3, 2, 1)
        plt1.errorbar(site_co2, site_co2_diffs, yerr=site_co2sd, fmt='o', color = SiteColours(sitein=site).outcolour, markersize = 3)
        plt1.plot([370,490], co2_lines.plusline,  ':', color = 'grey')
        plt1.plot([370,490], co2_lines.minusline,  ':', color = 'grey')        
        #plt1.set_ylim(Calcrange(co2_diffs, co2sd, even=1))
        lims = [plt1.get_ylim()[0],plt1.get_ylim()[1], -0.2, 0.2]
        plt1.set_ylim([min(lims),max(lims)])
        plt1.set_xlim([370,490])
            
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        plt1.yaxis.set_major_formatter(y_formatter)
        
        plt1.set_title('Cylinder intercomparison', fontsize = 7)
        plt1.set_ylabel(co2_ylabel, fontsize = 7)

        plt1.text(0.82, 0.85, inst_tag, horizontalalignment='center', verticalalignment='center', transform = plt1.transAxes, fontsize=8)

    
        # For CH4
        ch4_lines = WMOLines(ch4_diffs, gas = 'CH4')    
        if use_MPI == 1:
            ch4_lines = WMOLines(0, gas='CH4')
            
        if np.isnan(site_ch4).all():            
            inst_tag = ''
        else:
            inst_tag = site_co2inst[np.where(np.isfinite(site_ch4_diffs))[0][0]]
       
        plt2 = plt.subplot(3, 2, 2)
        plt2.errorbar(site_ch4, site_ch4_diffs, yerr=site_ch4sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
        plt2.plot([1750,2250], ch4_lines.plusline, ':', color = 'grey')
        plt2.plot([1750,2250], ch4_lines.minusline, ':', color = 'grey')
        
        #plt2.set_ylim(Calcrange(ch4_diffs, ch4sd, even = 1))
        lims = [plt2.get_ylim()[0],plt2.get_ylim()[1], -3, 3]
        plt2.set_ylim([min(lims),max(lims)])
        plt2.set_ylabel(ch4_ylabel, fontsize = 7)
        plt2.set_xlabel('Concentration', fontsize=7)
        plt2.yaxis.set_major_formatter(y_formatter)
        plt2.set_xlim([1750,2250])
        
        plt2.text(0.82, 0.85, inst_tag, horizontalalignment='center', verticalalignment='center', transform = plt2.transAxes, fontsize=8)

        
        # For N2O
        # Some sites don't have N2O
        n2o_lines = WMOLines(n2o_diffs, gas = 'N2O')
        if use_MPI == 1:
            n2o_lines = WMOLines(0, gas='N2O')

        if np.isnan(site_n2o).all():            
            inst_tag = ''
        else:
            inst_tag = site_n2oinst[np.where(np.isfinite(site_n2o_diffs))[0][0]]
            
        plt3 = plt.subplot(3, 2, 3)
        plt3.errorbar(site_n2o, site_n2o_diffs, yerr=site_n2osd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
        plt3.plot([322,336], n2o_lines.plusline, ':', color = 'grey')
        plt3.plot([322,336], n2o_lines.minusline, ':', color = 'grey')
                
        #plt3.set_ylim(Calcrange(n2o_diffs, n2osd, even=1))
        plt3.set_ylabel(n2o_ylabel, fontsize = 7)
        lims = [plt3.get_ylim()[0],plt3.get_ylim()[1], -0.3, 0.3]
        plt3.set_ylim([min(lims),max(lims)])        
        plt3.set_xlim([322,336])
        plt3.set_xlabel('Concentration', fontsize = 7)
        plt3.yaxis.set_major_formatter(y_formatter)
        
        plt3.text(0.82, 0.85, inst_tag, horizontalalignment='center', verticalalignment='center', transform = plt3.transAxes, fontsize=8)

        
        # for CO
        
        co_lines = WMOLines(co_diffs, gas = 'CO')
        if use_MPI == 1:
            co_lines = WMOLines(0, gas='CO')
  
        if np.isnan(site_co).all():    
            inst_tag = ''
        else:
            inst_tag = site_co2inst[np.where(np.isfinite(site_co_diffs))[0][0]]
            
        plt4 = plt.subplot(3, 2, 4)
        plt4.errorbar(site_co, site_co_diffs, yerr=site_cosd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
        plt4.plot([100,300], co_lines.plusline, ':', color = 'grey')
        plt4.plot([100,300], co_lines.minusline, ':', color = 'grey')
        
        #plt4.set_ylim(Calcrange(co_diffs, n2osd, even=1))
        plt4.set_ylabel(co_ylabel, fontsize = 7)
        plt4.set_xlim([100,300])
        lims = [plt4.get_ylim()[0],plt4.get_ylim()[1], -3, 3]
        plt4.set_ylim([min(lims),max(lims)])
        plt4.set_xlabel('Concentration', fontsize = 7)
        plt4.yaxis.set_major_formatter(y_formatter)
         
        plt4.text(0.82, 0.85, inst_tag, horizontalalignment='center', verticalalignment='center', transform = plt4.transAxes, fontsize=8)

        # for SF6
        sf6_lines = WMOLines(sf6_diffs, gas = 'SF6')
        if use_MPI == 1:
            sf6_lines = WMOLines(0, gas='SF6')

        if np.isnan(site_sf6).all():            
            inst_tag = ''
        else:
            inst_tag = site_n2oinst[np.where(np.isfinite(site_sf6_diffs))[0][0]]
       
        plt5 = plt.subplot(3, 2, 5)
        plt5.errorbar(site_sf6, site_sf6_diffs, yerr=site_sf6sd, fmt='o', color = SiteColours(sitein=site).outcolour,  markersize = 3)
        plt5.plot([6,14], sf6_lines.plusline, ':', color = 'grey')
        plt5.plot([6,14], sf6_lines.minusline, ':', color = 'grey')
                
        #plt5.set_ylim(Calcrange(sf6_diffs, sf6sd, even=1))
        plt5.set_ylabel(sf6_ylabel, fontsize = 7)
        plt5.set_xlim([6,14])
        lims = [plt5.get_ylim()[0],plt5.get_ylim()[1], -0.03, 0.03]
        plt5.set_ylim([min(lims),max(lims)])

        plt5.set_xlabel('Concentration', fontsize = 7)
        plt5.yaxis.set_major_formatter(y_formatter)

        plt5.text(0.82, 0.85, inst_tag, horizontalalignment='center', verticalalignment='center', transform = plt5.transAxes, fontsize=8)

        plt.figtext(0.56, 0.3, site, fontsize =8)
        plt.figtext(0.56, 0.26, 'CO$\mathregular{_2}$ - WMOx2007, through NPL (HAD & TIL), Scripps (WAO)', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.56, 0.23, 'CH$\mathregular{_4}$ - WMOx2004, x2004A (DECC/GAUGE/FERRY)', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.56, 0.2, 'N$\mathregular{_2}$O - WMOx2006A, SIO98 (DECC/GAUGE)', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.56, 0.17, 'CO- WMOx2014', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        plt.figtext(0.56, 0.14, 'SF$\mathregular{_6}$ - WMOx2006, SIO-SF6 (DECC/GAUGE)', verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        
        text1 = 'Mean +/- WMO compatibility goals'
        if use_MPI == 1:
            text1 =  'MPI +/- WMO compatibility goals'
        if use_EMPA == 1:
            text1 ='EMPA +/- WMO compatibility goals'
        
        plt.figtext(0.72, 0.11, text1, verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')


        plt.savefig(outputdir+site+'_single' + suffix + '.png', dpi=200)
            
        plt.show()   
        
        
# Make a summary plot for a single gas
class Plot_individualgas:
    def __init__(self, gas, N2O_D091968, N2O_D091969, N2O_D091970, CO2_D091968, CO2_D091969, CO2_D091970, \
     outputdir='/Users/as13988/Documents/Work/Cylinders/ICP/Results/Plots/', suffix='', use_MPI = 0, n2o_scaled =None, use_EMPA=0):
       
        # DNos = ['D091968', 'D091969', 'D091970'] 
       
        sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site, CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
        sites = list(chain(*sites))
        
        
        n2o_datetimes = [N2O_D091968.datetime, N2O_D091969.datetime, N2O_D091970.datetime]
        n2o_datetimes = list(chain(*n2o_datetimes))
        
        # Extract the data
        if any(gas in i for i in ['co2','ch4','co', 'CO2', 'CH4' ,'CO']):
            sites = [CO2_D091968.site, CO2_D091969.site, CO2_D091970.site]
            datetimes = [CO2_D091968.datetime, CO2_D091969.datetime, CO2_D091970.datetime]
            inst_all = [CO2_D091968.instrument, CO2_D091969.instrument, CO2_D091970.instrument]

        
            if any(gas in i for i in ['co2', 'CO2']):
                conc = [CO2_D091968.co2, CO2_D091969.co2, CO2_D091970.co2]
                sd = [CO2_D091968.co2sd, CO2_D091969.co2sd, CO2_D091970.co2sd]
                diffs = [CO2_D091968.co2_diffs, CO2_D091969.co2_diffs, CO2_D091970.co2_diffs]
                tags = [CO2_D091968.sitekey, CO2_D091969.sitekey, CO2_D091970.sitekey]
                
                ylabel = 'Site - mean [CO$\mathregular{_2}$] (ppm)'
                linerange = [370,490]
                scale_tag = 'CO$\mathregular{_2}$ - WMOx2007, through NPL (HAD & TIL), Scripps (WAO)'
                title = 'CO$\mathregular{_2}$'
                ylims = [-0.2,0.2]
                
                if use_MPI == 1:
                    diffs = [CO2_D091968.co2_MPIdiff, CO2_D091969.co2_MPIdiff, CO2_D091970.co2_MPIdiff]
                    ylabel = 'Site - MPI [CO$\mathregular{_2}$] (ppm)'

                if use_EMPA == 1:
                    diffs = [CO2_D091968.co2_EMPAdiff, CO2_D091969.co2_EMPAdiff, CO2_D091970.co2_EMPAdiff]
                    ylabel = 'Site - EMPA [CO$\mathregular{_2}$] (ppm)'


            if any(gas in i for i in ['ch4', 'CH4']):
                conc = [CO2_D091968.ch4, CO2_D091969.ch4, CO2_D091970.ch4]
                sd = [CO2_D091968.ch4sd, CO2_D091969.ch4sd, CO2_D091970.ch4sd]
                diffs = [CO2_D091968.ch4_diffs, CO2_D091969.ch4_diffs, CO2_D091970.ch4_diffs]
                tags = [CO2_D091968.sitekey, CO2_D091969.sitekey, CO2_D091970.sitekey]
                
                ylabel = 'Site - mean [CH$\mathregular{_4}$] (ppb)'
                linerange = [1750,2250]
                scale_tag = 'CH$\mathregular{_4}$ - WMOx2004, x2004A (DECC/GAUGE/FERRY)'
                title = 'CH$\mathregular{_4}$'
                ylims = [-3, 3]
                                
                if use_MPI == 1:
                    diffs = [CO2_D091968.ch4_MPIdiff, CO2_D091969.ch4_MPIdiff, CO2_D091970.ch4_MPIdiff]
                    ylabel = 'Site - MPI [CH$\mathregular{_4}$] (ppb)'

                if use_EMPA == 1:
                        diffs = [CO2_D091968.ch4_EMPAdiff, CO2_D091969.ch4_EMPAdiff, CO2_D091970.ch4_EMPAdiff]
                        ylabel = 'Site - EMPA [CH$\mathregular{_4}$] (ppb)'

            if any(gas in i for i in ['co', 'CO']):
                conc = [CO2_D091968.co, CO2_D091969.co, CO2_D091970.co]
                sd = [CO2_D091968.cosd, CO2_D091969.cosd, CO2_D091970.cosd]
                diffs = [CO2_D091968.co_diffs, CO2_D091969.co_diffs, CO2_D091970.co_diffs]
                tags = [CO2_D091968.sitekey, CO2_D091969.sitekey, CO2_D091970.sitekey]

                ylabel = 'Site - mean [CO] (ppb)'
                linerange = [100,300]
                scale_tag =  'CO- WMOx2014'
                title = 'CO'
                ylims = [-3, 3]
                
                if use_MPI == 1:
                    diffs = [CO2_D091968.co_MPIdiff, CO2_D091969.co_MPIdiff, CO2_D091970.co_MPIdiff]
                    ylabel = 'Site - MPI [CO] (ppb)'

                if use_EMPA == 1:
                        diffs = [CO2_D091968.co_EMPAdiff, CO2_D091969.co_EMPAdiff, CO2_D091970.co_EMPAdiff]
                        ylabel = 'Site - EMPA [CO] (ppb)'


        if any(gas in i for i in ['n2o','sf6','N2O', 'SF6']):
            sites = [N2O_D091968.site, N2O_D091969.site, N2O_D091970.site]
            datetimes = [N2O_D091968.datetime, N2O_D091969.datetime, N2O_D091970.datetime]           
            inst_all = [N2O_D091968.instrument, N2O_D091969.instrument, N2O_D091970.instrument]
                    
        
            if any(gas in i for i in ['n2o', 'N2O']):
                conc = [N2O_D091968.n2o, N2O_D091969.n2o, N2O_D091970.n2o]
                sd = [N2O_D091968.n2osd, N2O_D091969.n2osd, N2O_D091970.n2osd]
                diffs = [N2O_D091968.n2o_diffs, N2O_D091969.n2o_diffs, N2O_D091970.n2o_diffs]
                tags = [N2O_D091968.sitekey, N2O_D091969.sitekey, N2O_D091970.sitekey]

                ylabel = 'Site - mean [N$\mathregular{_2}$O] (ppb)'
                linerange = [322,336]
                scale_tag = 'N$\mathregular{_2}$O - WMOx2006A, SIO98 (DECC/GAUGE)'
                suffix ='_unscaled'
                if n2o_scaled != None:
                    scale_tag = 'N$\mathregular{_2}$O - WMOx2006A (DECC/GAUGE scaled)'
                    suffix ='_scaled'
                title = 'N$\mathregular{_2}$O'
                ylims = [-1.2, 0.6]
        
                if use_MPI == 1:
                    diffs = [N2O_D091968.n2o_MPIdiff, N2O_D091969.n2o_MPIdiff, N2O_D091970.n2o_MPIdiff]
                    ylabel = 'Site - MPI [N$\mathregular{_2}$O] (ppb)'

                if use_EMPA == 1:
                    diffs = [N2O_D091968.n2o_EMPAdiff, N2O_D091969.n2o_EMPAdiff, N2O_D091970.n2o_EMPAdiff]
                    ylabel = 'Site - EMPA [N$\mathregular{_2}$O] (ppb)'


            if any(gas in i for i in ['sf6', 'SF6']):
                conc = [N2O_D091968.sf6, N2O_D091969.sf6, N2O_D091970.sf6]
                sd = [N2O_D091968.sf6sd, N2O_D091969.sf6sd, N2O_D091970.sf6sd]
                diffs = [N2O_D091968.sf6_diffs, N2O_D091969.sf6_diffs, N2O_D091970.sf6_diffs]
                tags = [N2O_D091968.sitekey, N2O_D091969.sitekey, N2O_D091970.sitekey]
                
                ylabel = 'Site - mean [SF$\mathregular{_6}$] (ppt)'
                linerange = [6,14]
                scale_tag = 'SF$\mathregular{_6}$ - WMOx2006, SIO-SF6 (DECC/GAUGE)'
                title = 'SF$\mathregular{_6}$'
                ylims = [-0.02, 0.02]
                
                if use_MPI == 1:
                    diffs = [N2O_D091968.sf6_MPIdiff, N2O_D091969.sf6_MPIdiff, N2O_D091970.sf6_MPIdiff]
                    ylabel = 'Site - MPI [SF$\mathregular{_6}$] (ppb)'

                if use_EMPA == 1:
                    diffs = [N2O_D091968.sf6_EMPAdiff, N2O_D091969.sf6_EMPAdiff, N2O_D091970.sf6_EMPAdiff]
                    ylabel = 'Site - EMPA [SF$\mathregular{_6}$] (ppb)'
 
        # combine into one long list
        sites = list(chain(*sites))  
        datetimes = list(chain(*datetimes))
        inst_all = list(chain(*inst_all))
        conc = list(chain(*conc))
        sd = list(chain(*sd))
        diffs = list(chain(*diffs))
        tags = list(chain(*tags))
             
        # Find all the sites to cycle through
        sites_all = np.unique(sites)
             
        # Plot of 
        fig = plt.figure()
        fig.subplots_adjust(wspace = 0.6)
        #fig.tight_layout()
        fig.set_size_inches(6,6)
        fig.suptitle(title, fontsize = 12)
        sites_all = ['MPI','TTA', 'BSD','HFD','RGL','TAC','WAO','BT','FERRY','FAAM','HAD','TIL', 'Man','EMPA']        
        
        for i in np.arange(len(sites_all)):
            
            # find all the results from that site
            siteindex = (np.where(np.array(sites) == sites_all[i]))[0]
            
            site_conc = [conc[j] for j in siteindex]
            site_sd = [sd[j] for j in siteindex]
            site_inst = [inst_all[j] for j in siteindex]        
            site_diffs = [diffs[j] for j in siteindex]
            site_tags = [tags[j] for j in siteindex]
            
            # Plot the WMO compatibility goals
            lines = WMOLines(diffs, gas=gas)
            if use_MPI == 1:
                lines = WMOLines(0, gas=gas)
    
            if use_EMPA == 1:
                lines = WMOLines(0, gas=gas)
    
            if np.isnan(site_conc).all():        
                inst_tag = ''
            else:
                inst_tag = site_inst[np.where(np.isfinite(site_diffs))[0][0]]
                 
            
            # Plot up the data
            plt1 = plt.subplot(5, 3, i)
            fmts = ['o', 'v','s']
            for j in np.arange(len(np.unique(site_tags))):
                #pdb.set_trace()
                tag_j = np.unique(site_tags)[j]
                j_index = (np.where(np.array(site_tags) == tag_j))[0]
                
                plt1.errorbar(np.array(site_conc)[j_index], np.array(site_diffs)[j_index], yerr=np.array(site_sd)[j_index], fmt=fmts[j], color = SiteColours(sitein=sites_all[i]).outcolour, markersize = 3)
                print sites_all[i] + ' ' + np.array(site_inst)[j_index][0]
                yloc = 0.85 - (j * 0.08)
                plt1.text(0.15, yloc, sites_all[i] + ' ' + np.array(site_inst)[j_index][0]+ '-' + fmts[j], horizontalalignment='left', verticalalignment='center', transform = plt1.transAxes, fontsize=6)
                       
                #plt1.text(0.75, 0.85-(j*0.03), sites_all[i] + ' ' + np.array(site_inst)[j_index][0], horizontalalignment='center', verticalalignment='center', transform = plt1.transAxes, fontsize=6)

            plt1.plot(linerange, lines.plusline,  ':', color = 'grey')
            plt1.plot(linerange, lines.minusline,  ':', color = 'grey')        
            #plt1.set_ylim(Calcrange(co2_diffs, co2sd, even=1))
            lims = [plt1.get_ylim()[0],plt1.get_ylim()[1],ylims[0], ylims[1]]
            plt1.set_ylim([min(lims),max(lims)])
            plt1.set_xlim(linerange)
                
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
            plt1.set_ylabel(ylabel, fontsize = 5)
            
            lims = [plt1.get_ylim()[0],plt1.get_ylim()[1],ylims[0], ylims[1]]
            plt1.set_ylim([min(lims),max(lims)])
            plt1.set_xlim(linerange)
 
            for label in (plt1.get_xticklabels() + plt1.get_yticklabels()):
                label.set_fontname('Arial')
                label.set_fontsize(6)
            
            
        plt.figtext(0.07, 0.06, scale_tag, verticalalignment='bottom', horizontalalignment='left', fontsize=6)
        
        text1 = 'Mean +/- WMO compatibility goals'
        if use_MPI == 1:
            text1 = 'MPI +/- WMO compatibility goals'
        
        if use_EMPA == 1:
            text1 = 'EMPA +/- WMO compatibility goals'
        
        
        plt.figtext(0.07, 0.04, text1, verticalalignment='bottom', horizontalalignment='left', fontsize=6, color ='grey')


        plt.savefig(outputdir+gas+'_single' + suffix + '.png', dpi=200)
            
        plt.show()   