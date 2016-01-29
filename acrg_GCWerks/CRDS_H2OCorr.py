# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:24:37 2015

@author: as13988
"""
from acrg_GCWerks import acrg_read_GCwerks
import matplotlib.pylab as plt
import numpy as np
import matplotlib 
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
import pdb
from os import listdir
from os.path import isfile, join

# Find the files for a given site
def find_files(site, dir='/Users/as13988/Documents/Work/Picarro/H2OCorr/', date=None): 
    
    onlyfiles = [ f for f in listdir(dir) if isfile(join(dir,f)) ]
    sitefiles = [n for n in onlyfiles if site in n]
    dates = [int(a.split('.')[-2]) for a in sitefiles]
    
    if date is None:
        date = max(dates)
        
    files = [f for f in sitefiles if str(date) in f]

    dry_files = [dir + f for f in files if 'dry' in f]
    wet_files = [dir + f for f in files if 'wet' in f]
    
    return dry_files, wet_files

# Read multiple runs in
def read_data_multi(datafilelist): 
    
    datalist = []
    
    for i in datafilelist:
        data = acrg_read_GCwerks.read_gcexport_crds(i)   
        datalist.append(data)
        
    return datalist

# Reads in the GCWerks output
def read_data(datafile): 
    
    data = acrg_read_GCwerks.read_gcexport_crds(datafile)
        
    return data
    
# Plot the raw CO2, CH4 and H2O data
def plot_raw(data, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/'):
    
    splitstr = data.filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
    
    index = np.where(data.co2flags == 0)[0]
    
    plot_dt = [data.dt_date[i] for i in index]
    
    x_formatter = matplotlib.dates.DateFormatter('%H:%M')
    
    
    fig = plt.figure()
    fig.subplots_adjust(right = 0.8)
    fig.set_size_inches(8,4)
    
    # Plot CO2
    plt1 = plt.subplot(3, 1, 1)
    plt1.plot(plot_dt, data.co2_dry[index],  '-', color = 'red')
    plt1.plot(plot_dt, data.co2_wet[index],  ':', color = 'red')
    plt1.set_title(site + ' ' + testdate+ ' H2O Test')
    plt1.set_ylabel('CO$_2$ raw')
    plt1.xaxis.set_major_formatter(x_formatter)
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt1.yaxis.set_major_formatter(y_formatter)
    
    # Plot CH4
    plt2 = plt.subplot(3, 1, 2)
    plt2.plot(plot_dt, data.ch4_dry[index], '-', color = 'green')
    plt2.plot(plot_dt, data.ch4_wet[index], ':', color = 'green')
    plt2.set_ylabel('CH$_4$ raw')
    plt2.xaxis.set_major_formatter(x_formatter)
    plt2.yaxis.set_major_formatter(y_formatter)
    
    # Plot H2O
    plt3 = plt.subplot(3, 1, 3)
    plt3.plot(plot_dt, data.h2o[index], 'o-', color = 'blue')
    plt3.set_ylabel('H$_2$O raw')
    plt3.xaxis.set_major_formatter(x_formatter)
    plt3.yaxis.set_major_formatter(y_formatter)
    
    plt.figtext(0.82, 0.2, 'Dry = solid' , verticalalignment='bottom', horizontalalignment='left', fontsize=6)
    plt.figtext(0.82, 0.17, 'Wet = dashed', verticalalignment='bottom', horizontalalignment='left', fontsize=6)       

    print 'Plot saved as: ' + outdir+'Raw_'+tag+'_'+site+'.' + testdate + '.png'
    plt.savefig(outdir+'Raw_'+tag+'_'+site+'.' + testdate + '.png', dpi=200)
        
    plt.show()
    
    return 

# Plot the raw data against the H2O
def plot_H2O(data, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', species ='co2', label=''):
    
    splitstr = data.filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
    
    index = np.where(data.co2flags == 0)[0]
      
    fig = plt.figure()
    fig.subplots_adjust(left = 0.25)
    fig.set_size_inches(4,6)

    # Set the correct gas 
    wet_key = species + '_wet'
    dry_key = species + '_dry'
    wetsd_key = species + 'sd'


    labels_1 = {'co2' : 'CO$_2$ wet', 'ch4' :'CH$_4$ wet'}
    labels_2 = {'co2' : 'CO$_2$ raw', 'ch4' :'CH$_4$ raw'}
   
    
    
    # Plot dry vs CO2 wet
    plt1 = plt.subplot(2, 1, 1)
    wetdata = data.__dict__.get(wet_key)
    wetsddata = data.__dict__.get(wetsd_key)
    
    plt1.plot(data.h2o[index], wetdata[index],  'o', color = 'red')
    plt1.set_ylabel(labels_1[species])
    plt1.set_xlabel('H$_2$O')
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt1.yaxis.set_major_formatter(y_formatter)
    plt.figtext(0.27, 0.58,'Raw data - no correction', verticalalignment='bottom', horizontalalignment='left', fontsize=9)
    
    # Plot dry vs H2O
    plt2 = plt.subplot(2, 1, 2)
    drydata = data.__dict__.get(dry_key)
    #plt2.plot(data.h2o[index], drydata[index], 'o', color = 'green')
    plt2.errorbar(data.h2o[index], drydata[index], yerr= wetsddata[index], fmt='o', color = 'green')
    plt2.set_ylabel(labels_2[species])
    plt2.set_xlabel('H$_2$O raw')
    plt2.yaxis.set_major_formatter(y_formatter)
    plt.figtext(0.27, 0.42, 'Raw data - built in correction', verticalalignment='bottom', horizontalalignment='left', fontsize=9)
    
    plt.subplots_adjust(hspace = 0.3) 
 
    print 'Plot saved as: ' + outdir+'VsH2O_'+tag+'_'+site+'.' + testdate + '_'+ species+'.png'
    plt.savefig(outdir+'VsH2O_'+tag+'_'+site+'.' + testdate + '_'+ species + '.png', dpi=200)
        
    plt.show()
    
    return 


# Plot the ratio against the H2O and fit
def plot_ratio(data, drydata, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', species ='co2'):
    
    splitstr = data.filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
    
    index = np.where(data.co2flags == 0)[0]
      
    wet_key = species + '_wet'
    dry_key = species + '_dry'
    sd_key = species + 'sd'
     
    dry_means = calc_dry(drydata)
    wet_data = ((data.__dict__.get(wet_key))[index])
    dry_data = ((data.__dict__.get(dry_key))[index])
    sd_data = ((data.__dict__.get(sd_key))[index])
    
    ratio = wet_data/dry_means[dry_key]     
    
    h2o = data.h2o[index]
    #fit_params = np.polyfit(h2o, ratio, 2)
    #fit = fit_params[2] + (fit_params[1]*h2o) + (fit_params[0]*(h2o**2))
        
    fit_params = curve_fit(responsefn, h2o, ratio, sigma = sd_data )
    #print fit_params
    fit = responsefn(h2o, fit_params[0][0], fit_params[0][1])
    
    

    labels_1 = {'co2' : 'CO$_2$ wet/CO$_2$ dry', 'ch4' :'CH$_4$ wet/CH$_4$ dry'}
     
    
    
    # Plot the ratio
    fig = plt.figure()
    fig.set_size_inches(4.2,4)
    fig.subplots_adjust(left = 0.2)
    fig.subplots_adjust(bottom = 0.15)
    
    # Plot CO2 wet/ CO2 dry vs H2O
    plt3 = plt.subplot(1,1,1)
    plt3.plot(h2o, ratio, 'o', color = 'blue')
    #plt3.errorbar(h2o, ratio, yerr = sd_data, fmt= 'o', color = 'blue')
    plt3.set_xlabel('H$_2$O raw')
    plt3.set_ylabel(labels_1[species])
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt3.yaxis.set_major_formatter(y_formatter)
    
    plt.plot(data.h2o[index], fit, '-', color = 'red')
 
    fit_str = 'y = ' + '1 + ' + str("{:.3E}".format(fit_params[0][0])) + 'x + ' + str("{:.3E}".format(fit_params[0][1])) + 'x^2'
 
    plt.figtext(0.25, 0.2, fit_str , verticalalignment='bottom', horizontalalignment='left', fontsize=6)
 
 
    print 'Plot saved as: ' + outdir+'Ratio_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png'
    plt.savefig(outdir+'Ratio_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png', dpi=200)
        
    plt.show()
    plt.close()
    
    # PLot the fit  of the new correction compared to the in build correction
    plot = plot_fit(wet_data, dry_data, h2o, dry_means[wet_key], fit_params, sd_data)

    plt.figtext(0.22, 0.3, 'Built in correction', color ='red', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
    plt.figtext(0.22, 0.25, 'New correction', color = 'blue', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
    
    
    print 'Plot saved as: ' + outdir+'Corr_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png'
    plot.savefig(outdir+'Corr_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png', dpi=200)
        
    plot.show()
    plot.close()
    
    return fit_params


def plot_fit(rawwetdata, instcorrdata, h2o, drymean, fit_params, sd, colour = 'blue'):
        
    fit = responsefn(h2o, fit_params[0][0], fit_params[0][1])    
    
    CO2_corr = rawwetdata/fit
    # err = 100 - (correct value/dry raw value * 100 )
    err_new_fit = 100-(CO2_corr/drymean*100)
    err_inst_fit = 100-(instcorrdata/drymean*100)
    
    # Plot the corrected output
    fig = plt.figure()
    fig.set_size_inches(4.2,4)
    fig.subplots_adjust(left = 0.2)
    fig.subplots_adjust(bottom = 0.2)

    # Plot CO2 wet/ CO2 dry vs H2O
    plt3 = plt.subplot(1,1,1)
    plt3.errorbar(h2o, err_new_fit, yerr = sd, fmt= 'o', color = colour)
    plt3.errorbar(h2o, err_inst_fit, yerr = sd, fmt= 'o', color = 'red')
    plt3.set_xlabel('H$_2$O raw')
    plt3.set_ylabel('% error ')
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt3.yaxis.set_major_formatter(y_formatter)
         
    plt.show()
    
    return plt
    

# data =  the set of wet data that you want the different fits to be applied to
# drydata = list of means. You need to have the same number of means as you have sets of fit_parameters. 
# These can be the same mean or can vary with the fit parameters
# fit_params = a list of outputs from plot_ratio one per fit that you want to compare
def compare_fits(data, drydata, fit_params,  outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', species ='co2'):
        
    splitstr = data.filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
    
    # Create key names for extraction
    wet_key = species + '_wet'
    dry_key = species + '_dry'
    sd_key = species + 'sd'
    
    # Make index for unflagged data
    index = np.where(data.co2flags == 0)[0]
    
    # Extract the unflagged data 
    wet_data = ((data.__dict__.get(wet_key))[index])
    dry_data = ((data.__dict__.get(dry_key))[index])
    sd_data = ((data.__dict__.get(sd_key))[index])
    h2o = ((data.h2o)[index])
    
    colours = ['blue', 'green', 'purple', 'orange']

    # PLot the fit  of the new correction compared to the in build correction

    # calculate parameters
    err_inst_fit = 100-(dry_data/drydata[0]*100)


    # Plot the corrected output
    fig = plt.figure()
    fig.set_size_inches(4.2,4)
    fig.subplots_adjust(left = 0.2)
    fig.subplots_adjust(bottom = 0.15)
    
    # Plot CO2 wet/ CO2 dry vs H2O
    plt3 = plt.subplot(1,1,1)
    
    plt3.errorbar(h2o, err_inst_fit, yerr = sd_data, fmt= 'o', color = 'red')

    for i in np.arange(len(fit_params)):
        fit_params_i =fit_params[i]
        fit = responsefn(h2o, fit_params_i[0][0], fit_params_i[0][1])    
        
        CO2_corr = wet_data/fit
        
        # err = 100 - (correct value/dry raw value * 100 )
        err_new_fit = 100-(CO2_corr/drydata[i]*100)
        
        plt3.errorbar(h2o, err_new_fit, yerr = sd_data, fmt= 'o', color = colours[i])
        
        plt.figtext(0.22, 0.2 +i*0.04, 'New correction ' +str(i+1), color = colours[i], verticalalignment='bottom', horizontalalignment='left', fontsize=10)


    plt3.set_xlabel('H$_2$O raw')
    plt3.set_ylabel('% error ')
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt3.yaxis.set_major_formatter(y_formatter)
    
    plt.figtext(0.22, 0.16, 'Built in correction', color ='red', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        
    plt.show()
    
    print 'Plot saved as: ' + outdir + 'FitComparison_' + tag+ '.' + site + testdate +'.'+ species + '.png'
    plt.savefig(outdir + 'FitComparison_' + tag+ '.' + site + testdate + species + '.png', dpi=200)
    
    plt.close()


# plots each run on the same plot along with the fit
def compare_runs(data, drydata,  outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', species ='co2'):
    
    splitstr = data[0].filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
   
    # Create key names for extraction
    wet_key = species + '_wet'
    sd_key = species + 'sd'
    
    colours = ['b', 'g', 'r']    
    
    labels_1 = {'co2' : 'CO$_2$ wet/CO$_2$ dry', 'ch4' :'CH$_4$ wet/CH$_4$ dry'}    
    
    # want to do an overall fit so I need to collate the data    
    h2o_all = []
    ratio_all = []
    sigma_all = []
    
    for i in np.arange(len(data)):
        # Make index for unflagged data
        index = np.where(data[i].co2flags == 0)[0]
        
        # Extract the unflagged data 
        wet_data = ((data[i].__dict__.get(wet_key))[index])
        sd_data = ((data[i].__dict__.get(sd_key))[index])
        h2o = ((data[i].h2o)[index])
   
        # calculate the ratio for that run     
        ratio = wet_data/drydata[i]
    
        # Append all data together to do an overall fit
        h2o_all.extend(h2o)    
        ratio_all.extend(ratio)
        sigma_all.extend(sd_data)
        
        # Calculate the fit
        fit_params = curve_fit(responsefn, h2o, ratio, sigma = sd_data)
        
        if i == 0:
            # Plot the ratio
            fig = plt.figure()
            fig.set_size_inches(4.2,4)
            fig.subplots_adjust(left = 0.2)
            fig.subplots_adjust(right = 0.9)
            fig.subplots_adjust(bottom = 0.15)
            
            # Plot CO2 wet/ CO2 dry vs H2O
            plt3 = plt.subplot(1,1,1)
            
        plt3.plot(h2o, ratio, markerfacecolor ='none', markeredgecolor = colours[i], marker='.', linestyle='none')
       
        fit_str = 'y = ' + '1 + ' + str("{:.3E}".format(fit_params[0][0])) + 'x + ' + str("{:.3E}".format(fit_params[0][1])) + 'x^2'

        plt.figtext(0.25, 0.18 + i*0.04, fit_str, color=colours[i], verticalalignment='bottom', horizontalalignment='left', fontsize=6)

    # calculate a fit to all the data
    fit_params = curve_fit(responsefn, h2o_all, ratio_all, sigma = sigma_all )
    
    fit_str = 'y = ' + '1 + ' + str("{:.3E}".format(fit_params[0][0])) + 'x + ' + str("{:.3E}".format(fit_params[0][1])) + 'x^2 Fit to all data'
    plt.figtext(0.35, 0.85, fit_str, color='black', verticalalignment='bottom', horizontalalignment='left', fontsize=6)

   
    plt3.set_xlabel('H$_2$O raw')
    plt3.set_ylabel(labels_1[species])
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt3.yaxis.set_major_formatter(y_formatter)
 
    print 'Plot saved as: ' + outdir+'RatioComparison'+site+'.' + testdate+ '.'+ species + '.png'
    plt.savefig(outdir+'RatioComparison'+site+'.' + testdate+ '.'+ species + '.png', dpi=200)
        
    plt.show()
    plt.close()


def calc_dry(data, ignoreindex=0):
    
    index = np.where(data.co2flags == 0)[0]    
    
    #pdb.set_trace()
    
    if ignoreindex == 0:
        co2wet_mean = np.mean((data.co2_wet)[index])
        co2dry_mean = np.mean((data.co2_dry)[index])
    
        ch4wet_mean = np.mean((data.ch4_wet)[index])
        ch4dry_mean = np.mean((data.ch4_dry)[index])
    
        h2o_mean = np.mean((data.h2o)[index])

        co2wet_sd = np.std((data.co2_wet)[index])
        co2dry_sd = np.std((data.co2_dry)[index])
    
        ch4wet_sd = np.std((data.ch4_wet)[index])
        ch4dry_sd = np.std((data.ch4_dry)[index])
    
        h2o_sd = np.std((data.h2o)[index])

    
    else:
        co2wet_mean = np.mean(data.co2_wet)
        co2dry_mean = np.mean(data.co2_dry)
    
        ch4wet_mean = np.mean(data.ch4_wet)
        ch4dry_mean = np.mean(data.ch4_dry)
    
        h2o_mean = np.mean(data.h2o)

        co2wet_sd = np.std(data.co2_wet)
        co2dry_sd = np.std(data.co2_dry)
    
        ch4wet_sd = np.std(data.ch4_wet)
        ch4dry_sd = np.std(data.ch4_dry)
    
        h2o_sd = np.std(data.h2o)

    means = {'co2_wet' : co2wet_mean, \
            'co2_dry' : co2dry_mean, \
            'ch4_wet' : ch4wet_mean, \
            'ch4_dry': ch4dry_mean, \
            'h2o' : h2o_mean, \
            'co2_wet_sd' : co2wet_sd, \
            'co2_dry_sd' : co2dry_sd, \
            'ch4_wet_sd' : ch4wet_sd, \
            'ch4_dry_sd': ch4dry_sd, \
            'h2o_sd' : h2o_sd}
    
    return means


# Plot the raw data against the H2O
def plot_drymeans(means, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/'):
    
    fig = plt.figure()
    fig.subplots_adjust(left = 0.25)
    fig.set_size_inches(4,6)      
    
    for i in np.arange(len(means)):
     
        data = means[i]        
        
        # Plot CO2 wet
        plt1 = plt.subplot(3, 1, 1)
        plt1.errorbar(i+1, data['co2_wet'], yerr=data['co2_wet_sd'],  fmt='o', color = 'red')
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        plt1.yaxis.set_major_formatter(y_formatter)
        
        # Plot CH4 wet
        plt2 = plt.subplot(3, 1, 2)
        plt2.errorbar(i+1, data['ch4_wet'], yerr=data['ch4_wet_sd'],  fmt='o', color = 'green')
        plt2.yaxis.set_major_formatter(y_formatter)

        # Plot H2O
        plt3 = plt.subplot(3, 1, 3)
        plt3.errorbar(i+1, data['h2o'], yerr=data['h2o_sd'],  fmt='o', color = 'blue')
        plt3.yaxis.set_major_formatter(y_formatter)

    plt.subplots_adjust(hspace = 0.3) 
    
    plt1.set_ylabel('CO$_2$')
    plt2.set_ylabel('CH$_4$')
    plt3.set_ylabel('H$_2$O')
 
    print 'Plot saved as: ' + outdir+'drymeans.png'
    plt.savefig(outdir+'drymeans.png', dpi=200)
        
    plt.show()
    
    return 


def responsefn(H2O, A, B):
    
    ratio = 1 + A*H2O + B*H2O*H2O
    
    return ratio