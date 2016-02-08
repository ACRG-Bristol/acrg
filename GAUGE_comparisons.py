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
reload(acrg_agage)

# gases available
def gases_measured():
    gases = ['co2','CO2','ch4','CH4', 'n2o','N2O','co','CO','sf6','SF6']
    
    return gases

# list of sites
class sites_list():
    def __init__(self):
        sites_tag = ['MHD','TTA','BSD','HFD','RGL','TAC', 'ferry','WAO','HAD','TIL', 'GLA','EHL', 'WAO_UCAM', 'BTT']
        sites_name = ['Mace Head','Angus', 'Bilsdale','Heathfield','Ridge Hill','Tacolneston', 'Ferry','Weybourne','Haddenham','Tilney','Glatton',"Earl's Hall", 'Weybourne UCam', 'BT Tower']
    
    
        self.sites_tag = sites_tag
        self.sites_name = sites_name


# units information
def units(gas):
    units = ['ppm', 'ppm', 'ppb','ppb', 'ppb','ppb','ppb','ppb','ppt','ppt']

    index = np.where(np.array(gases_measured()) == gas)[0]   
    if len(index) == 1:
        unit = units[index]
    else:
        unit = ''
        
    return '['+unit+']'


# PLOT CODE!
# Make the plot and save it    
def plotagainstother(site, gas, suffix='', yscale=None, xscale=None, other_site='MHD',\
                     alt_dir_site = None, alt_dir_other = None, dateformatter =None, network=None,\
                     site_name = None, other_site_name = None, outdir = None, n2o_scale_site = None,\
                     n2o_scale_other = None):

    # Check if the given site matches
    if (np.where(np.array(sites_list().sites_tag) == site)[0]).size > 0: 
        if site_name is None:
            site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag ) == site)[0]]
            if n2o_scale_site:
                site_name = site_name + '_scaled'
                suffix = 'scaled'
    else:
        print 'Site ' + site + ' does not match any listed'
        print 'Please use: ' 
        for i in sites_list().sites_tag:
            print i
    
    if (np.where(np.array(sites_list().sites_tag) == other_site)[0]).size > 0:   
        if other_site_name is None:
            other_site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag) == other_site)[0]]    
        if n2o_scale_other:
            other_site_name = other_site_name + '_scaled'
    else:
        print 'Site ' + other_site + ' does not match any listed'
        print 'Please use: ' 
        for i in sites_list().sites_tag:
            print i
            
    if outdir is None:    
        outdir = '/home/as13988/GAUGE/MHDComparisons/'

    if site in ['WAO_UCAM']:
        network = 'UCAM-EA' 
    
    site_data = read_site(site, gas, network=network, alt_dir = alt_dir_site, n2o_scale = n2o_scale_site)
        
    if (site_data) is not None:
        start = min(site_data.index)
        end = max(site_data.index)
        
        if other_site == 'MHD' and gas == 'co2':
            network='LSCE'
        else:
            network = None
            
        other_site_data = read_site(other_site, gas, start = start, end = end, network=network,\
        alt_dir = alt_dir_other, n2o_scale = n2o_scale_other)
             
        if other_site_data is not None:             
             
            fig = plt.figure()
            ax = plt.subplot()
            ax.plot(site_data.index, site_data['mf'], 'ob', mec='blue', ms =2)
            ax.plot(other_site_data.index, other_site_data['mf'], 'or', mec='red', ms =2)
            
            # Adjust the tick number on the x axis
            x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
            if dateformatter is not None:
                x_formatter = matplotlib.dates.DateFormatter(dateformatter)

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
            plt.savefig(outdir+site+'Vs' + other_site + '_' + gas+suffix + '.png', dpi=200)
            
        else:
                 print 'Can not find data for ' + gas + ' at ' + other_site 
    else:
        print 'Can not find data for ' + gas + ' at ' + site 

# Make the plot and save it    
def plotagainstflask(site, flasksite, gas, suffix=None, yscale=None, xscale=None, network=None,\
                 yscale_diff=None, xscale_diff=None, n2o_scale = None, outdir = None, no_bar = None):
    
    site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag ) == site)[0]] 
    if n2o_scale:
        site_name = site_name + '_scaled'  
        suffix = 'scaled'
    
    flask_site_name = (sites_list().sites_name)[np.where(np.array(sites_list().sites_tag) == flasksite)[0]]    
    
    if outdir is None:
        outdir = '/home/as13988/GAUGE/FlaskComparisons/'
        
    site_data = read_site(site, gas, n2o_scale = n2o_scale)
    
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
        
            if suffix is None :
                plt.savefig(outdir+site+'Vs' + flasksite + '_flask_' + gas + '.png', dpi=200)
            else:
                plt.savefig(outdir+site+'Vs' + flasksite + '_flask_' + gas+'_' + suffix + '.png', dpi=200)

        
        # Time match the flask record to the closet in situ record and plot just these
        matchedindex = []
        matcheddiff = []
        extendedindex = []
        means = []
        sds = []
        
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
                means.append(np.mean(site_data['mf'][np.arange(30)-15+after]))
                sds.append(np.std(site_data['mf'][np.arange(30)-15+after]))
            else:
                matchedindex.append(before)
                matcheddiff.append(beforediff) 
                extendedindex.append(np.arange(30)-15+before)
                means.append(np.mean(site_data['mf'][np.arange(30)-15+before]))
                sds.append(np.std(site_data['mf'][np.arange(30)-15+before]))
                
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
               
        # Calculate differences
        diffs = np.array(means) -  flask_data['mf']
        
        # Remove differences flasks more than 30mins away from the in situ data
        for k in toofar[0]:
            diffs[k] = np.nan        
            sds[k] = np.nan
          
        
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
        
        #ax1.plot(matchedtime, matchedinsitu -  flask_data['mf'],'ob', mec='blue', ms =2)
        if no_bar is None:        
            ax1.errorbar(matchedtime, diffs,yerr=sds,fmt='ob', mec='blue', ms =2)
        else:
            ax1.plot(matchedtime, diffs,'ob', mec='blue', ms =2)
      
        # Adjust the tick number on the x axis
        x_formatter = matplotlib.dates.DateFormatter('%m/%Y')
        x_tickno_formatter = matplotlib.ticker.MaxNLocator(5)
        ax1.xaxis.set_major_formatter(x_formatter)
        ax1.xaxis.set_major_locator(x_tickno_formatter)
        ax1.set_ylabel('insitu - flask' + units(gas))
        
        # Plot WMO lines
        lines = acrg_ICP.WMOLines(0, gas)
        daterange = [extendedtime[0], extendedtime[-1]]  

        ax1.plot(daterange, lines.plusline,'--')
        ax1.plot(daterange, lines.minusline,'--')
        
        if xscale_diff != None:
            ax1.set_xlim(xscale_diff)
        
        if yscale_diff != None:
            ax1.set_ylim(yscale_diff)  
            
        plt.figtext(0.5, 0.9, site + ' Vs. ' + flasksite, verticalalignment='bottom', horizontalalignment='center', fontsize=12)
        plt.figtext(0.83, 0.84, 'Flask', verticalalignment='bottom', horizontalalignment='left', fontsize=10,color='red')
        plt.figtext(0.75, 0.84, 'In situ', verticalalignment='bottom', horizontalalignment='left', fontsize=10, color='blue')
        plt.figtext(0.15, 0.57, 'Closest minute', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        plt.figtext(0.15, 0.84, 'Closest 30 minutes', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        plt.figtext(0.15, 0.3, 'Closest 30 minutes diff', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        

        if suffix is None :
            plt.savefig(outdir+site+'Vs' + flasksite + '_flaskmatched_' + gas + '.png', dpi=200)
        else:
            plt.savefig(outdir+site+'Vs' + flasksite + '_flaskmatched_' + gas+'_' + suffix + '.png', dpi=200)

        plt.close()
        
        # Plot differences against concentration
        ax2 = plt.subplot()        
        ax2.errorbar(np.array(means), np.array(means) -  flask_data['mf'],yerr = sds,fmt='ob', mec='blue', ms =2)
        
        # Adjust the tick number on the x axis
        ax2.set_ylabel('insitu - flask' + units(gas))
        ax2.set_xlabel('insitu'+ units(gas))
        
        # Plot WMO lines
        lines = acrg_ICP.WMOLines(0, gas)
        daterange = [extendedtime[0], extendedtime[-1]]  

        ax2.plot([extendedinsitu.min(), extendedinsitu.max()], lines.plusline,'--')
        ax2.plot([extendedinsitu.min(), extendedinsitu.max()], lines.minusline,'--')
        plt.figtext(0.15, 0.84, site + ' Vs. ' + flasksite, verticalalignment='bottom', horizontalalignment='left', fontsize=12)
        
        if suffix is None :
            plt.savefig(outdir+site+'Vs' + flasksite + '_flaskconc_' + gas + '.png', dpi=200)
        else:
            plt.savefig(outdir+site+'Vs' + flasksite + '_flaskconc_' + gas+'_' + suffix + '.png', dpi=200)

    else:
        print 'Can not find data for ' + gas + ' at ' + site 


# READ CODE!
# Read in the site data for the given gas
def read_site(site, gas, start = '1900-01-01 00:00:00', end = '2020-01-01 00:00:00', network=None, \
            alt_dir = None, n2o_scale = None):
   
    
    if site in ['WAO','wao']:
        site_data = read_wao(gas).data   
        
    #if site in ['GLA','gla']:
        #site_data = read_gla(gas).data   
    
    if site in ['MHD','mhd']:
        site_data = read_MHD(gas, start = start, end=end, alt_dir=alt_dir)
    
    if site in ['WAO_UCAM', 'wao_ucam']:  
        site_data = acrg_agage.get('WAO', gas, start=start, end = end, network=network)

    if site in ['FAAM', 'faam']:
        if gas in ['co2','CO2', 'ch4', 'CH4']:
            instrument = 'FGGA'
        if gas in ['n2o','N2O']:
            instrument = 'QCL'
        if gas in ['co','CO']:
            instrument = 'AL5002'
            
        site_data = acrg_agage.get('FAAM', gas, instrument = instrument, start=start, end = end, network='GAUGE')

    if site in ['ferry', 'Ferry']:  
        site_data = acrg_agage.get('ferry', gas, start=start, end = end, network=network, status_flag_unflagged = [0,1,2])

    #if site not in ['WAO','wao','GLA','gla','MHD','mhd','WAO_UCAM', 'wao_ucam']:
    if site not in ['WAO','wao','MHD','mhd','WAO_UCAM', 'wao_ucam', 'ferry', 'Ferry']:
        site_data = acrg_agage.get(site, gas, start=start, end = end, network=network, alt_dir=alt_dir)


    # Temporary scaling to adjust to x2006A using 2014 factor
    if n2o_scale and gas in ['n2o', 'N2O']:
        site_data['mf'] = site_data['mf']*1.0009

    return site_data


# Read in the MHD data for the same date range as the site
def read_MHD(gas, start = '1900-01-01 00:00:00', end = '2020-01-01 00:00:00', alt_dir=None):
        
    if gas in ['co2', 'CO2']:
        MHD = acrg_agage.get('MHD', 'CO2', instrument='CRDS', network='LSCE', start=start, end = end, alt_dir = alt_dir)
    else:
        MHD = acrg_agage.get('MHD', gas, start=start, end = end, alt_dir = alt_dir)
    
    return MHD


# Read in the flask files
class read_flask():
    def __init__(self, site, gas, start = None, end = None):
        
        # need to set some parameters to get the data to reaTACd in properly
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
        
        if gas =='N2O' or gas =='n2o':
            infile = 'WAO_N2O_SF6_forAnn_151015.csv'
            sd_header = 'N2O_flag'
            OK = 2
            
        meancolumnheaders = ['CO2_ppm','CO2_ppm','CH4_ppb','CH4_ppb','N2O_ppb','N2O_ppb','CO_ppb','CO_ppb','SF6_ppb','SF6_ppb']         
        
        if OK != None:        
                
            index = np.where(np.array(gases_measured()) == gas)[0]    
            mean_header = meancolumnheaders[index]
            
            print 'Reading : ' + infile
            
            # read in the file
            f = open(datadir+infile, 'r')
            lines = f.readlines()

            # Hack for the N2O and SF6 csv file
            if OK  == 2:  
                noheaderlines = 3
                # read in the header for the columns
                data = pandas.read_csv(datadir+infile, header=noheaderlines-1, delimiter=',')
            
            # read in the first line which give you the number of header lines
            else:
                noheaderlines = int((lines[0].split())[0])
                # read in the header for the columns
                data = pandas.read_csv(datadir+infile, header=noheaderlines-1, delim_whitespace=True)
            
            
            if infile != 'wao-CO2_weybourne_20120101_compilation_hourlyaverage.na': 
                if OK == 2:
                    time = [dt.datetime.strptime(data['Date_Time'][i], "%d/%m/%Y %H:%M") for i in np.arange(len(data['Date_Time']))]                   
                    sd = np.arange(len(data.index))
                    sd[:] = np.NaN
                else:                    
                    time = [dt.datetime(int(data['year'][i]),int(data['month'][i]),int(data['day'][i]),int(data['hour'][i])) for i in np.arange(len(data['year']))]
                    sd = np.array(data[sd_header])
            else: 
                time = [dt.datetime.strptime(data.index[i]+ ' ' + data['Date_and_Time'][i], "%d/%m/%Y %H:%M:%S") for i in np.arange(len(data['Date_and_Time']))]
                sd = np.empty(len(time))
                sd[:] = np.NAN
            
            conc = np.array(data[mean_header])
            
            gooddata = np.where(conc < 9999.0)[0] # Some files have 9999.0 some have 9999.99
    
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