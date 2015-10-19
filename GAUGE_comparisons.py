# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:24:18 2015

@author: as13988
"""

import acrg_agage
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas
import datetime as dt
import pandas as pd
import netCDF4
import bisect
from acrg_GCWerks import acrg_ICP

# gases available
def gases_measured():
    gases = ['co2','CO2','ch4','CH4', 'n2o','N2O','co','CO','sf6','SF6']
    
    return gases

# list of sites
class sites_list():
    def __init__(self):
        sites_tag = ['MHD','TTA','BSD','HFD','RGL','TAC', 'ferry','WAO','HAD','TIL', 'GLA','EHL']
        sites_name = ['Mace Head','Angus', 'Bilsdale','Heathfield','Ridge Hill','Tacolneston', 'Ferry','Weybourne','Haddenham','Tilney','Glatton',"Earl's Court"]
    
    
        self.sites_tag = sites_tag
        self.sites_name = sites_name


# units information
def units(gas):
    units = ['ppm', 'ppm', 'ppb','ppb', 'ppb','ppb','ppb','ppb','ppt','ppt']

    index = np.where(np.array(gases_measured()) == gas)[0]    
    unit = units[index]
    return '['+unit+']'


# PLOT CODE!
# Make the plot and save it    
def plotagainstother(site, gas, suffix='', yscale=None, xscale=None, other_site='MHD', network=None):
    
    site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag ) == site)[0]] 
    other_site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag) == other_site)[0]]    
    
    
    outputdir = '/home/as13988/GAUGE/MHDComparisons/'
    site_data = read_site(site, gas)
    
    if (site_data) is not None:
        start = min(site_data.index)
        end = max(site_data.index)
        
        if other_site == 'MHD' and gas == 'co2':
            network='LSCE'
        
        other_site_data = read_site(other_site, gas, start = start, end = end, network=network)
             
        if other_site_data is not None:             
             
            fig = plt.figure()
            ax = plt.subplot()
            ax.plot(site_data.index, site_data['mf'], 'ob', mec='blue', ms =2)
            ax.plot(other_site_data.index, other_site_data['mf'], 'or', mec='red', ms =2)
            
            # Adjust the tick number on the x axis
            x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
            x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
            ax.xaxis.set_major_formatter(x_formatter)
            ax.xaxis.set_major_locator(x_tickno_formatter)
            ax.set_ylabel(gas + ' ' + units(gas))
        
            if yscale != None:
                ax.set_ylim(yscale)
        
            if xscale != None:
                ax.set_xlim(xscale)
        
                        
        
            plt.figtext(0.2, 0.8, gas, verticalalignment='bottom', horizontalalignment='left', fontsize=12)
        
            plt.figtext(0.7, 0.8, other_site_name, verticalalignment='bottom', horizontalalignment='left', fontsize=10,color='red')
            plt.figtext(0.7, 0.75, site_name, verticalalignment='bottom', horizontalalignment='left', fontsize=10, color='blue')
        
            if suffix  != '':
                suffix = '_' + suffix
            plt.savefig(outputdir+site+'Vs' + other_site + '_' + gas+suffix + '.png', dpi=200)
            
        else:
                 print 'Can not find data for ' + gas + ' at ' + other_site 
    else:
        print 'Can not find data for ' + gas + ' at ' + site 

# Make the plot and save it    
def plotagainstflask(site, flasksite, gas, suffix='', yscale=None, xscale=None, network=None):
    
    site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag ) == site)[0]] 
    flask_site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag) == flasksite)[0]]    
    
    outputdir = '/home/as13988/GAUGE/MHDComparisons/'
    site_data = read_site(site, gas)
    
    if (site_data) is not None:
        start = min(site_data.index)
        end = max(site_data.index)
        
        flask_data = (read_flask(flasksite, gas, start = start, end = end)).data
             
        if flask_data is not None:             
            fig = plt.figure()
            ax = plt.subplot()
            ax.plot(site_data.index, site_data['mf'], 'ob', mec='blue', ms =2)
            ax.plot(flask_data.index, flask_data['mf'], 'or', mec='red', ms =2)
            
            # Adjust the tick number on the x axis
            x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
            x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
            ax.xaxis.set_major_formatter(x_formatter)
            ax.xaxis.set_major_locator(x_tickno_formatter)
            ax.set_ylabel(gas + ' ' + units(gas))
        
            if yscale != None:
                ax.set_ylim(yscale)
        
            if xscale != None:
                ax.set_xlim(xscale)
                        
        
            plt.figtext(0.2, 0.8, gas, verticalalignment='bottom', horizontalalignment='left', fontsize=12)
        
            plt.figtext(0.65, 0.8, flask_site_name + ' flask', verticalalignment='bottom', horizontalalignment='left', fontsize=10,color='red')
            plt.figtext(0.65, 0.75, site_name+ ' in situ', verticalalignment='bottom', horizontalalignment='left', fontsize=10, color='blue')
        
            if suffix  != '':
                suffix = '_' + suffix
            plt.savefig(outputdir+site+'Vs' + flasksite + '_flask_' + gas+suffix + '.png', dpi=200)
        
        # resample the dataframe to get 15min means and then match those
        #site_data_15 = site_data.resample('15T') 
        #site_data_15_SD = site_data.resample('15T', how=np.std) 
        # Time match the flask record to the closet in situ record and plot just these
        matchedindex = []
        matcheddiff = []
        extendedindex = []
        for i in np.arange(len(flask_data.index)):
            # Find the points either side of the flask data point
            after = bisect.bisect_left(site_data.index, flask_data.index[i])
            before = bisect.bisect_right(site_data.index, flask_data.index[i])
            
            afterdiff = abs((flask_data.index[i] - site_data.index[after]).total_seconds())
            beforediff = abs((flask_data.index[i] - site_data.index[before]).total_seconds())
            
            if afterdiff <= beforediff:
                matchedindex.append(after)
                matcheddiff.append(afterdiff) 
                extendedindex.append(np.arange(30)-15+after) #get the indicies 15 on each side
            else:
                matchedindex.append(before)
                matcheddiff.append(beforediff) 
                extendedindex.append(np.arange(30)-15+before)
                
        # Extract the data at the matched indicies
        matchedtime = [site_data.index[j] for j in matchedindex]
        matchedinsitu = [site_data['mf'][j] for j in matchedindex]
        matchedinsitu = np.array(matchedinsitu)
        
        extendedtime = [site_data.index[j] for j in extendedindex]
        extendedinsitu = [site_data['mf'][j] for j in extendedindex]
        extendedinsitu = np.array(extendedinsitu)      
        
        
        # Remove points which are more than 30 mins away from each other
        toofar = np.where(np.array(matcheddiff) > 1800)
        matchedinsitu[toofar[0]] = np.nan
        
             
        fig = plt.figure()
        f, (ax0, ax, ax1) = plt.subplots(3, sharex = True)
        ax0.plot(extendedtime, extendedinsitu, 'ob', mec='blue', ms =2)
        ax0.plot(flask_data.index, flask_data['mf'], 'or', mec='red', ms =2)
        
        # Adjust the tick number on the x axis
        x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
        x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
        ax0.xaxis.set_major_formatter(x_formatter)
        ax0.xaxis.set_major_locator(x_tickno_formatter)
        ax0.set_ylabel(gas + ' ' + units(gas))

        if yscale != None:
            ax0.set_ylim(yscale)
    
        if xscale != None:
            ax0.set_xlim(xscale)
        
        
        ax.plot(matchedtime, matchedinsitu, 'ob', mec='blue', ms =2)
        ax.plot(flask_data.index, flask_data['mf'], 'or', mec='red', ms =2)
        
        # Adjust the tick number on the x axis
        x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
        x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
        ax.xaxis.set_major_formatter(x_formatter)
        ax.xaxis.set_major_locator(x_tickno_formatter)
        ax.set_ylabel(gas + ' ' + units(gas))

        if yscale != None:
            ax.set_ylim(yscale)
    
        if xscale != None:
            ax.set_xlim(xscale)
        
        ax1.plot(matchedtime, matchedinsitu -  flask_data['mf'],'ob', mec='blue', ms =2)
        # Adjust the tick number on the x axis
        x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
        x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
        ax1.xaxis.set_major_formatter(x_formatter)
        ax1.xaxis.set_major_locator(x_tickno_formatter)
        ax1.set_ylabel('insitu - flask' + units(gas))
    
        if xscale != None:
            ax1.set_xlim(xscale)
        
        # Plot WMO lines
        lines = acrg_ICP.WMOLines(0, gas)
        daterange = [extendedtime[0], extendedtime[-1]]  

        ax1.plot(daterange, lines.plusline,'--')
        ax1.plot(daterange, lines.minusline,'--')

        plt.figtext(0.5, 0.9, site + ' Vs. ' + flasksite, verticalalignment='bottom', horizontalalignment='center', fontsize=12)
        plt.figtext(0.83, 0.84, 'Flask', verticalalignment='bottom', horizontalalignment='left', fontsize=10,color='red')
        plt.figtext(0.75, 0.84, 'In situ', verticalalignment='bottom', horizontalalignment='left', fontsize=10, color='blue')
        plt.figtext(0.15, 0.57, 'Closest minute', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        plt.figtext(0.15, 0.84, 'Closest 30 minutes', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        
        if suffix  != '':
            suffix = '_' + suffix
        
        plt.savefig(outputdir+site+'Vs' + flasksite + '_flaskmatched_' + gas+suffix + '.png', dpi=200)
        plt.close()
        
        # Plot differences against concentration
        ax2 = plt.subplot()        
        ax2.plot(matchedinsitu, matchedinsitu -  flask_data['mf'],'ob', mec='blue', ms =2)
        
        # Adjust the tick number on the x axis
        ax2.set_ylabel('insitu - flask' + units(gas))
        ax2.set_xlabel('insitu'+ units(gas))
        
        # Plot WMO lines
        lines = acrg_ICP.WMOLines(0, gas)
        daterange = [extendedtime[0], extendedtime[-1]]  

        ax2.plot([extendedinsitu.min(), extendedinsitu.max()], lines.plusline,'--')
        ax2.plot([extendedinsitu.min(), extendedinsitu.max()], lines.minusline,'--')
        plt.figtext(0.15, 0.84, site + ' Vs. ' + flasksite, verticalalignment='bottom', horizontalalignment='left', fontsize=12)

        if suffix  != '':
            suffix = '_' + suffix
        plt.savefig(outputdir+site+'Vs' + flasksite + '_flaskconc_' + gas+suffix + '.png', dpi=200)
 
    else:
        print 'Can not find data for ' + gas + ' at ' + site 


# READ CODE!
# Read in the site data for the given gas
def read_site(site, gas, start = '1900-01-01 00:00:00', end = '2020-01-01 00:00:00', network=None):
   
    if site in ['WAO','wao']:
        site_data = read_wao(gas).data   
        
    if site in ['GLA','gla']:
        site_data = read_gla(gas).data   
    
    if site in ['MHD','mhd']:
        site_data = read_MHD(gas, start = start, end=end)
    
    if site not in ['WAO','wao','GLA','gla','MHD','mhd']:
        site_data = acrg_agage.get(site, gas, start=start, end = end, network=network)

    return site_data


# Read in the MHD data for the same date range as the site
def read_MHD(gas, start = '1900-01-01 00:00:00', end = '2020-01-01 00:00:00'):
        
    if gas in ['co2', 'CO2']:
        MHD = acrg_agage.get('MHD', 'CO2', instrument='CRDS', network='LSCE', start=start, end = end)
    else:
        MHD = acrg_agage.get('MHD', gas, start=start, end = end)
    
    return MHD


# Read in the WAO files
class read_flask():
    def __init__(self, site, gas, start = None, end = None):
        
        # need to set some parameters to get the data to read in properly
        network = None
        height = None
        
        if site =='MHD' or gas =='mhd': 
            network = 'GAUGE'
        
        if site == 'TAC' or site == 'tac':
            network = 'GAUGE'
            height = '185m'
                
        df = acrg_agage.get(site, gas, network=network, height=height, start=start, end=end)   
                                       
        self.data = df


# Read in the WAO files
class read_wao():
    def __init__(self, gas):
        filenames = ['wao-CO2_weybourne_20120101_compilation_hourlyaverage.na',\
                    'wao-rga3_weybourne_20080306_h2-co-compilation-hourlyaverage.na',\
                    'wao-ghg-gc-fid_weybourne_20130306_ch4-compilation.na']
        datadir = '/home/as13988/GAUGE/WAO/'
        
        # need to figure out which file we want
        file_tag = gas
        OK = None
        
        if gas =='co2' or gas =='CO2': 
            infile ='wao-CO2_weybourne_20120101_compilation_hourlyaverage.na'
            OK =1
            
        if gas == 'CO' or gas == 'co':
            infile = 'wao-rga3_weybourne_20080306_h2-co-compilation-hourlyaverage.na'
            sd_header = 'CO_ppb_sd'
            OK = 1
            
        if gas =='CH4' or gas =='ch4':
            infile = 'wao-ghg-gc-fid_weybourne_20130306_ch4-compilation.na'
            sd_header = 'CH4_WT_ppb'
            OK = 1
        
        meancolumnheaders = ['CO2_ppm','CO2_ppm','CH4_ppb','CH4_ppb','N2O_ppb','N2O_ppb','CO_ppb','CO_ppb','SF6_ppb','SF6_ppb']         
        
        if OK != None:        
                
            index = np.where(np.array(gases_measured()) == gas)[0]    
            mean_header = meancolumnheaders[index]
            
            print 'Reading : ' + infile
            
            # read in the file
            f = open(datadir+infile, 'r')
            lines = f.readlines()
            
            # read in the first line which give you the number of header lines
            noheaderlines = int((lines[0].split())[0])
            
            # read in the header for the columns
            data = pandas.read_csv(datadir+infile, header=noheaderlines-1, delim_whitespace=True)
            
            if infile != 'wao-CO2_weybourne_20120101_compilation_hourlyaverage.na': 
                time = [dt.datetime(int(data['year'][i]),int(data['month'][i]),int(data['day'][i]),int(data['hour'][i])) for i in np.arange(len(data['year']))]
                sd = np.array(data[sd_header])
            else: 
                time = [dt.datetime.strptime(data.index[i]+ ' ' + data['Date_and_Time'][i], "%d/%m/%Y %H:%M:%S") for i in np.arange(len(data['Date_and_Time']))]
                sd = np.empty(len(time))
                sd[:] = np.NAN
            
            conc = np.array(data[mean_header])
            
            gooddata = np.where(conc != 9999.0)[0]
    
            goodtime = [time[j] for j in gooddata]        
            
            df = pd.DataFrame({"mf": conc[gooddata], 'vmf': sd[gooddata]},
                                      index = goodtime)             
        else:
            df = None
                                       
        self.data = df

# Read in the GLA files
class read_gla():
    def __init__(self, gas):

               
        filenames = ['GAUGE-insitu-fts_GLA_20141022_ch4-10m.nc',\
                        'GAUGE-insitu-fts_GLA_20141022_co2-10m.nc',\
                        'GAUGE-insitu-fts_GLA_20141022_co-10m.nc', \
                        'GAUGE-insitu-fts_GLA_20141022_n2o-10m.nc']

        datadir = '/data/shared/obs/GAUGE/'
        
        # need to figure out which file we want
        file_tag = gas
        OK = None
        
        if gas =='co2' or gas =='CO2': 
            infile = filenames[1]
            OK = 1
        if gas == 'CO' or gas =='co':
            infile = filenames[2]
            OK = 1
        if gas =='CH4' or gas == 'ch4':
            infile = filenames[0]
            OK = 1
        if gas =='N2O' or gas =='n2o':
            infile = filenames[3]
            OK = 1
            
        if OK == 1:        
                
            print 'Reading : ' + infile
            
            # read in the file
            ncfile = netCDF4.Dataset(datadir+infile, 'r')
            
            startdate = ncfile.variables['datestart']
            startdate_dt = [dt.datetime.strptime(str(i), '%Y%m%d') for i in startdate[:]]
            startdatesec = ncfile.variables['datesecstart']
            time = [startdate_dt[i] + dt.timedelta(seconds=long(startdatesec[i])) for i in np.arange(len(startdatesec))]

            df = pd.DataFrame({"mf": ncfile.variables['mf'][:]}, index = time)
            df['umf'] = ncfile.variables['u_mf'][:]
        else:
            df = None
                                       
        self.data = df



# Plot all the sites for a given gas
def plotallsites(gas):
    for i in sites_list().sites_tag:
        print 'Plotting for site: ' + i
        plotagainstother(i, gas)