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
import scipy
import pdb
from os import listdir
from os.path import isfile, join
from operator import itemgetter
from itertools import groupby
from scipy.odr import ODR, Model, Data, RealData

# Read in the data for a given list of sites and list of years
# calculate the fits and do the plotting
# defaults to a dry_cutoff of 0.0001
# h2o_sd is a filter based on h2o standard deviation
# co2_sd is a filter based on co2 standard deviation
# ch4_sd is a filter based on ch4 standard deviation
# ignore first ignores the first x mins of the wet section defaults to 5
# use_drygaps = set if you want to identify the wet periods as the "non dry" bits
# h2o_jump = set to the gap between consecutive h2o measurments that's used to identify the start of wet section defaults to 0.4
# H2Orange = set the H2O x axis range on the residual plots
def calcread_multi(sites, years, basedir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/', dir_suffix = '', quick_plot = None, \
h2o_jump = 0.4, ignorefirst = 5, dry_cutoff = 0.005, h2o_sd = 0.2, co2_sd = 0.05, ch4_sd = 0.4,  H2O_Range = None, No_plots = None, sd_flag = None):
    
    outputs = {}
    means_all = {}
    
    #loop through the sites
    for i in sites:
        #loop through the years
        for j in years:
            #find the files
            print 'Processing ' + i + ' ' + str(j)
            raw_files = find_files(i, date=j, indir = basedir + 'raw/'+dir_suffix)
            instcorr_files = find_files(i, date=j, indir = basedir + 'inbuiltcorr/'+dir_suffix)

            # if there are separate wet and dry files process the old way
            if raw_files[1] is not None:
                #read in the data
                wet_raw = read_data_multi(raw_files[1])
                wet_instcorr = read_data_multi(instcorr_files[1])
                dry_data = read_data_multi(raw_files[0])

            # or process the new way
            else:

                wet_raw = []
                dry_data = []
                wet_instcorr = []
                
                for k, item_k in enumerate(raw_files[0]):
                    wet_k, dry_k = define_wetdry(basedir + 'raw/'+dir_suffix +'/'+ item_k, quick_plot=quick_plot, ignorefirst = ignorefirst, h2o_jump = h2o_jump, dry_cutoff = dry_cutoff, h2o_sd = h2o_sd, co2_sd = co2_sd)
                    instwet_k, instdry_k = define_wetdry(basedir + 'inbuiltcorr/'+dir_suffix + '/'+ item_k, ignorefirst = ignorefirst, h2o_jump = h2o_jump, dry_cutoff = dry_cutoff, h2o_sd = h2o_sd, co2_sd = co2_sd)
                
                    # check that there was dry data for the run
                    # if so keep the data
                    # if not just discard it
                    if len(dry_k.co2) != 0:
                        wet_raw.extend(wet_k)    
                        dry_data.extend([dry_k for l in np.arange(len(wet_k))])
                        wet_instcorr.extend(instwet_k)          
                                    
            # calculate the fits and plot individual runs
            # THIS ASSUMES THAT THERE ARE THE SAME NUMBER OF INSTCORR FILES AS THERE ARE NORMAL FILES!
            output, means = plot_ratio_multi(wet_raw, dry_data, wet_instcorr, dir_suffix = dir_suffix, No_plots = No_plots, sd_flag = sd_flag)
            
            # plot the change from the built in residual for ALL runs
            plot_residual(means, wet_raw, wet_instcorr, dry_data, saveas=1, dir_suffix = dir_suffix, useinbuilt=1, H2O_Range = H2O_Range, sd_flag = sd_flag)
            
            outputs[i + '_' + str(j)] = output
            means_all[i + '_' + str(j)] = means
    
    return outputs, means_all


# Read in a file that contains both wet and dry data
# Identify the dry periods based on a given H2O tolerance
# The identify the wet periods and how many "runs" there have been in the episode.
# defaults to a dry_cutoff of 0.0001
# h2o_sd is a filter based on h2o standard deviation
# co2_sd is a filter based on co2 standard deviation
# ch4_sd is a filter based on ch4 standard deviation
# ignore first ignores the first x mins of the wet section defaults to 5
# use_drygaps = set if you want to identify the wet periods as the "non dry" bits
# h2o_jump = set to the gap between consecutive h2o measurments that's used to identify the start of wet section defaults to 0.4
def define_wetdry(infile, dry_cutoff = 0.0001, h2o_sd = 0.2, h2o_jump = 0.4, co2_sd = 0.05, ch4_sd = 0.4, ignorefirst = 5,\
 use_drygaps = None, quick_plot = None):

    # read in the data
    data_all = read_data(infile)

    # Extract the good H2O, CO2 and CH4 data and datetimes
    good = (np.where((data_all.co2flags == 0) & (h2o_sd > data_all.h2osd) & (co2_sd > data_all.co2sd)& (ch4_sd > data_all.ch4sd)))[0]
    
    if len(good) == 0:
        print 'There are no good points in file: ' + infile
    h2o = data_all.h2o[good]
    
    datetime = [data_all.datetime[i] for i in good]    
    
    # use the dry_sections class to extract these
    dry_data = dry_sections(data_all, dry_cutoff)
    
    if len(dry_data.co2) == 0:
        print 'There was no data with H2O < '+ str(dry_cutoff) + ' for file ' + data_all.filename
        
    # This determines the runs based on intervening dry sections but not every test run has dry bewtween it
    if use_drygaps is not None:
        # Find the gaps
        diff_days = [i.days for i in (dry_datetime - np.roll(dry_datetime, 1))]
        diff_seconds = [i.seconds for i in (dry_datetime - np.roll(dry_datetime, 1))]
        
        diff_minutes = (np.array(diff_days)*24*60) + (np.array(diff_seconds)/60)
    
        # gaps will be where the difference is >1
        # this index is the FIRST point of each block
        gaps_start = np.concatenate(([0],np.array(np.where(diff_minutes > 1)[0])))
        gaps_end = np.concatenate((np.array(np.where(diff_minutes > 1)[0]), [len(diff_minutes)]))
    
    # else do it based on where rapid periods of increase in H2O occurs
    else:
        # find where the H2O has increased rapidly
        index = np.where( (np.roll(h2o, -1) - h2o) > h2o_jump)[0]

        # check if it's accidentally picked up the last point
        # If so discard it
        if index[-1] >= len(h2o)-1:
            index = index[:-1]

        # find where these increases are grouped together
        gaps_start = []
        for k, g in groupby(enumerate(index), lambda (i,x):i-x):
            group =  map(itemgetter(1), g)
            #print group
            gaps_start.append(group[0])
        
        # define the end of the run
        if len(gaps_start) > 1:
            gaps_end =  np.concatenate((np.array(gaps_start[1:])-1, [len(h2o)-1]))
        else:
            gaps_end = [len(h2o)-1]                

        # need to remove blocks of very short length
        gap_length = [gaps_end[i] - gaps_start[i] for i in np.arange(len(gaps_start))]
        
        big_gap = np.where(np.array(gap_length) > 2)[0]
        
        gaps_start = np.array(gaps_start)[big_gap]
        gaps_end = np.array(gaps_end)[big_gap]

        # DECIDED TO IGNORE THE FIRST 15 MINS OF EACH RUN
        gaps_start = gaps_start + ignorefirst



        print 'Have identified ' + str(len(gaps_start)) + ' runs'
        
    # Extract the data for each run
    wet_runs = []
    for i, start_i in enumerate(gaps_start):
        
        # use the wet_sections class to extract these
        run_i = wet_sections(data_all, start_i, gaps_end[i], co2_sd = co2_sd, h2o_sd = h2o_sd)
        wet_runs.append(run_i) 

    # Do a quick plot of the H2O data to show if the separate runs have been identified
    if quick_plot is not None:
        plt.plot(datetime, h2o, mec='0.5', mfc='0.5', marker = 'o', linestyle = '')

        # find the number of runs
        no_runs =  len(wet_runs)
        
        for i in np.arange(no_runs):    
            run_i = wet_runs[i]
            blue = (0,0, (float(i)+1)/(float(no_runs)+1))
            plt.plot(run_i.datetime, run_i.h2o, mec=blue, mfc=blue, marker = 'o', linestyle = '')
        
        plt.plot(dry_data.datetime, dry_data.h2o, mec='r', mfc='r', marker = 'o', linestyle = '')
        plt.show()
        
    return wet_runs, dry_data


# extract the dry sections
class dry_sections:
    def __init__(self, data_all, dry_cutoff):
        good = (np.where(data_all.co2flags == 0))[0]
        h2o = data_all.h2o[good]
    
        # Find the "dry" sections
        dry_index = (np.where(h2o < dry_cutoff))[0]
        
        # want to replicate the class structure of data_all         
        for i in data_all.__dict__.keys():
            # if it's a numpy array then extract the data using the index
            if isinstance(getattr(data_all, i), np.ndarray):
                setattr(self, i, getattr(data_all, i)[good][dry_index])
                
            # if it's a list then extract the data using the list syntax and the index
            if isinstance(getattr(data_all, i), list):
                dummy_1 = [getattr(data_all, i)[j] for j in good]             
                dummy_2 = [dummy_1[j] for j in dry_index] 
                setattr(self, i, dummy_2)
                                                
            # if it's anything else keep it all
            if not isinstance(getattr(data_all, i), np.ndarray) and not isinstance(getattr(data_all, i), list):
                setattr(self, i, getattr(data_all, i))

        
# extract the wet sections
class wet_sections:
    def __init__(self, data_all, start_gap, end_gap, h2o_sd = 0.2, co2_sd = None):
    
        good = (np.where((data_all.co2flags == 0) & (h2o_sd > data_all.h2osd) & (co2_sd > data_all.co2sd)))[0]

        # want to replicate the class structure of data_all         
        for i in data_all.__dict__.keys():
            # if it's a numpy array then extract the data using the start and ends of the gaps
            if isinstance(getattr(data_all, i), np.ndarray):
                setattr(self, i, getattr(data_all, i)[good][start_gap:end_gap])
                
            # if it's a list then extract the data using the list syntax and the  start and ends of the gaps
            if isinstance(getattr(data_all, i), list):
                dummy_1 = [getattr(data_all, i)[j] for j in good]       
                dummy_2 = [dummy_1[j] for j in (np.arange(end_gap - start_gap)+start_gap)] 
                setattr(self, i, dummy_2)
                                                
            # if it's anything else keep it all
            if not isinstance(getattr(data_all, i), np.ndarray) and not isinstance(getattr(data_all, i), list):
                setattr(self, i, getattr(data_all, i))


# Find the files for a given site
# Syntax: dryfiles, wetfiles = CRDS_H2OCorr.findfiles(site)
def find_files(site, indir='/Users/as13988/Documents/Work/Picarro/H2OCorr/Raw/', date=None): 
    
    onlyfiles = [ f for f in listdir(indir) if isfile(join(indir,f)) ]
    sitefiles = [n for n in onlyfiles if site in n]
    dates = [(i).split('.')[-3] for i in sitefiles]
    
    if date is None:
        date = ''
        
    files = [f for f,g in zip(sitefiles,dates) if str(date) in g]
    dry_files = [indir + f for f in files if 'dry' in f]
    wet_files = [indir + f for f in files if 'wet' in f]    
       
    if len(files) > 0  and len(wet_files) == 0:
        return files, None
    else:       
        return dry_files, wet_files

# Read multiple runs in
def read_data_multi(datafilelist): 
    
    datalist = []
    
    if isinstance(datafilelist, tuple):
        datafilelist = datafilelist[0]
    
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
   
    print 'Plotting: ' + data.filename
    
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
    labelx = -0.1
        
    plt1 = plt.subplot(3, 1, 1)
    plt1.plot(plot_dt, data.co2_dry[index],  '.', color = 'red')
    plt1.plot(plot_dt, data.co2_wet[index],  '*', color = 'red')
    plt1.set_title(site + ' ' + testdate+ ' H2O Test')
    plt1.set_ylabel('CO$_2$ raw')
    plt1.xaxis.set_major_formatter(x_formatter)
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt1.yaxis.set_major_formatter(y_formatter)
    plt1.yaxis.set_label_coords(labelx, 0.5)
        
    # Plot CH4
    plt2 = plt.subplot(3, 1, 2)
    plt2.plot(plot_dt, data.ch4_dry[index], '.', color = 'green')
    plt2.plot(plot_dt, data.ch4_wet[index], '*', color = 'green')
    plt2.set_ylabel('CH$_4$ raw')
    plt2.xaxis.set_major_formatter(x_formatter)
    plt2.yaxis.set_major_formatter(y_formatter)
    plt2.yaxis.set_label_coords(labelx, 0.5)
    
    # Plot H2O
    plt3 = plt.subplot(3, 1, 3)
    plt3.plot(plot_dt, data.h2o[index], 'o-', color = 'blue')
    plt3.set_ylabel('H$_2$O raw')
    plt3.xaxis.set_major_formatter(x_formatter)
    plt3.yaxis.set_major_formatter(y_formatter)
    plt3.yaxis.set_label_coords(labelx, 0.5)
   
    plt.figtext(0.82, 0.2, 'Dry = solid' , verticalalignment='bottom', horizontalalignment='left', fontsize=6)
    plt.figtext(0.82, 0.17, 'Wet = dashed', verticalalignment='bottom', horizontalalignment='left', fontsize=6)       

    print 'Plot saved as: ' + outdir+ '/' + site+ '/'+'Raw_'+tag+'_'+site+'.' + testdate + '.png'
    plt.savefig(outdir+ '/' + site+ '/'+'Raw_'+tag+'_'+site+'.' + testdate + '.png', dpi=200)
        
    plt.show()
    
    return 

# Plot the raw data against the H2O
# syntax: CRDS_H2OCorr.plot_H2O(data)
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
 
    print 'Plot saved as: ' + outdir+ '/' + site+ '/'+'VsH2O_'+tag+'_'+site+'.' + testdate + '_'+ species+'.png'
    plt.savefig(outdir+ '/' + site+ '/'+'VsH2O_'+tag+'_'+site+'.' + testdate + '_'+ species + '.png', dpi=200)
        
    plt.show()
    
    return 


# plot multiple run ratios for both CO2 and CH4
# NB: wetdata and drydata are lists containing the data as read in using read_data_multi
# nocorrdata is the data without the instrument specific correction applied
def plot_ratio_multi(wetdata, drydata, nocorrdata=None, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', dir_suffix = '', h2o_tag = 'h2o', No_plots = None, sd_flag=None):

    gases = ['CO2', 'CH4']
    species = ['co2','ch4']
    
    output = {}
                 
    for j in np.arange(len(gases)):
        params=[]
        corr=[]
        err=[]
        h2o=[]
        sd=[]
        h2o_sd=[]
        filenames=[]
        dates=[]
        
        for i in np.arange(len(wetdata)):
            if nocorrdata is not None:
                nocorrdata_i = nocorrdata[i]
            else:
                nocorrdata_i = None

            # optional to use different dry periods
            if isinstance(drydata, list):            
                a,b,c,d,e,f = plot_ratio(wetdata[i], drydata[i], nocorrdata = nocorrdata_i, outdir=outdir, dir_suffix = dir_suffix, species = species[j], h2o_tag=h2o_tag, No_plots = No_plots, sd_flag = sd_flag)
            
            # Use all dry data as a single dry period       
            else:
                a,b,c,d,e, f = plot_ratio(wetdata[i], drydata, nocorrdata =  nocorrdata_i, outdir=outdir, dir_suffix = dir_suffix, species = species[j], h2o_tag=h2o_tag, No_plots = No_plots, sd_flag = sd_flag)
   
            
            params.append(a)
            corr.append(b)            
            h2o.append(c)
            err.append(d)
            sd.append(e)
            h2o_sd.append(f)
            filenames.append(wetdata[i].filename)
            dates.append(wetdata[i].datetime[0])
            
        gas_j = {'params' : params,\
                    'corr' : corr, \
                    'residual' : err, \
                    'sd' : sd, \
                    'h2o' : h2o, \
                    'h2o_sd' : h2o_sd, \
                    'filenames' : filenames, \
                    'dates' : dates}


        output[gases[j]] = gas_j
    
    # Print out the fits
    print 'Fit Y =  1 + a*X + b*(X*X)'
    print output['CO2']['filenames'][0][0:3]
    print 'CO2 parameters'
    print 'datetime, a, a_sd, b, b_sd'
    a_co2 = []
    b_co2 = []
    
    for i, item in enumerate(output['CO2']['params']):
        print str(output['CO2']['dates'][i]) + ', ' \
                    + '{:.7f}'.format(item[0]) +', ' \
                    + '{:.7f}'.format(output['CO2']['corr'][i][0]) + ', ' \
                    + '{:.7f}'.format(item[1]) + ', ' \
                    + '{:.7f}'.format(output['CO2']['corr'][i][1])
        a_co2.append(item[0])
        b_co2.append(item[1])
    
    mean_a_co2 = np.mean(np.array(a_co2))
    stdev_a_co2 = np.std(np.array(a_co2))
    mean_b_co2 = np.mean(np.array(b_co2))
    stdev_b_co2 = np.std(np.array(b_co2))
    
    print ''
    print 'mean' + ', ' \
        + '{:.7f}'.format(mean_a_co2) + ', ' \
        + '{:.7f}'.format(stdev_a_co2) +  ', ' \
        + '{:.7f}'.format(mean_b_co2) + ', ' \
        + '{:.7f}'.format(stdev_b_co2)
       
    
    print ''
    print ''
    print 'CH4 parameters'
    print 'datetime, a, a_sd, b, b_sd'
    a_ch4 = []
    b_ch4 = []

    for i, item in enumerate(output['CH4']['params']):
        print str(output['CH4']['dates'][i]) + ', ' \
                    + '{:.6f}'.format(item[0]) +', ' \
                    + '{:.6f}'.format(output['CH4']['corr'][i][0]) + ', ' \
                    + '{:.6f}'.format(item[1]) + ', ' \
                    + '{:.6f}'.format(output['CH4']['corr'][i][1])
        
        a_ch4.append(item[0])
        b_ch4.append(item[1])
        
    mean_a_ch4 = np.mean(np.array(a_ch4))
    stdev_a_ch4 = np.std(np.array(a_ch4))
    mean_b_ch4 = np.mean(np.array(b_ch4))
    stdev_b_ch4 = np.std(np.array(b_ch4))
    
    print ''
    print 'mean' + ', ' \
            + '{:.7f}'.format(mean_a_ch4) + ', ' \
            + '{:.7f}'.format(stdev_a_ch4) + ', ' \
            + '{:.7f}'.format(mean_b_ch4) + ', ' \
            + '{:.7f}'.format(stdev_b_ch4)

    means = {"a_co2" : [mean_a_co2], \
            "a_co2_sd" : [stdev_a_co2], \
            "b_co2" : [mean_b_co2], \
            "b_co2_sd" : [stdev_b_co2], \
            "a_ch4" : [mean_a_ch4], \
            "a_ch4_sd" : [stdev_a_ch4], \
            "b_ch4" : [mean_b_ch4], \
            "b_ch4_sd" : [stdev_b_ch4], \
            "site" : wetdata[0].filename.split('_')[0],\
            'date' : wetdata[0].filename.split('.')[-3]}

    return output, means

# Plot the ratio against the H2O and fit
# syntax: fit = CRDS_H2OCorr.plot_ratio(wetdata, drydata)
def plot_ratio(data, drydata, nocorrdata=None, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', dir_suffix = '', species ='co2', h2o_tag = 'h2o', No_plots = None, sd_flag = None):
 
    splitstr = data.filename.split('.')
    site = (splitstr[0])[0:3]
    testdate = splitstr[1]
    
    tag = (splitstr[0].split('_'))[-1]
    
    index = np.where(data.co2flags == 0)[0]
    
    # optional index based on co2 standard deviation
    if sd_flag is not None:
        index = np.where((np.array(data.co2flags) == 0) & (np.array(data.co2sd) < sd_flag))[0]

    wet_key = species + '_wet'
    dry_key = species + '_dry'
    sd_key = species + 'sd'
    
    
    dry_means = calc_dry(drydata)
    wet_data = ((data.__dict__.get(wet_key))[index])
    if nocorrdata is not None:
        instcorr_data = ((nocorrdata.__dict__.get(dry_key))[index])
    else:
        instcorr_data = None
    sd_data = ((data.__dict__.get(sd_key))[index])
    h2o_sd = ((data.h2osd)[index])    
     
    # pull out the correct h2o column
    # NB: GCWerks use "h2o" and ICOS "h2o_reported"
    h2o = (data.__dict__.get(h2o_tag))[index]

    # NOTE using the "wet" valve from the dry means as this is uncorrected using the in built correction
    # and as it's dry it shouldn't need correcting
    ratio = wet_data/dry_means[wet_key]     

    
    # new fitting with x weights as well as y weights
    # copied example from https://stackoverflow.com/questions/26058792/correct-fitting-with-scipy-curve-fit-including-errors-in-x
    dataforfit = RealData(h2o, ratio, sx=h2o_sd, sy=sd_data)
    model = Model(ODRfunc)
    
    odr = ODR(dataforfit, model, [1,1])
    odr.set_job(fit_type=2)
    output = odr.run()
    
    #print output.pprint()
    params = output.beta
    corr = np.square(output.sd_beta)
    # Original fitting
    #params, corr = scipy.optimize.curve_fit(fitfunc, h2o, ratio, [1,1], sd_data)


    fit = fitfunc(h2o, params[0], params[1])
    
    
    
    labels_1 = {'co2' : 'CO$_2$ wet/CO$_2$ dry', 'ch4' :'CH$_4$ wet/CH$_4$ dry'}
    xlabel = 'H$_2$O'
    
    if h2o_tag == 'h2o_reported':
        xlabel = 'H$_2$O reported'
    
    if No_plots is None:
        # Plot the ratio
        fig = plt.figure()
        fig.set_size_inches(4.2,4)
        fig.subplots_adjust(left = 0.2)
        fig.subplots_adjust(bottom = 0.15)
        
        # Plot CO2 wet/ CO2 dry vs H2O
        plt3 = plt.subplot(1,1,1)
        plt3.plot(h2o, ratio, 'o', color = 'blue')
        #plt3.errorbar(h2o, ratio, yerr = sd_data, fmt= 'o', color = 'blue')
        plt3.set_xlabel(xlabel)
        plt3.set_ylabel(labels_1[species])
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        plt3.yaxis.set_major_formatter(y_formatter)
        
     
        
        plt.plot(h2o, fit, '-', color = 'red')
     
        fit_str = 'y = ' + '1 + ' + str("{:.3E}".format(params[0])) + 'x + ' + str("{:.3E}".format(params[1])) + 'x^2'
        print fit_str
        print np.sqrt(corr)
        
        plt.figtext(0.25, 0.2, fit_str , verticalalignment='bottom', horizontalalignment='left', fontsize=8)
     
        outname = 'Ratio_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png'
        if h2o_tag == 'h2o_reported': 
           outname = 'Ratio_'+tag+'_'+site+'.' + testdate+ '_'+ species + '_ICOS.png'

        print 'Plot saved as: ' + outdir+ '/' + site+ '/'+dir_suffix+outname
        plt.savefig(outdir+ '/' + site+ '/'+dir_suffix+outname, dpi=200)
            
        plt.show()
        plt.close()
    
    # PLot the fit of the new correction compared to the in build correction
    outname = 'Corr_'+tag+'_'+site+'.' + testdate+ '_'+ species + '.png'
    if h2o_tag == 'h2o_reported': 
       outname = 'Corr_'+tag+'_'+site+'.' + testdate+ '_'+ species + '_ICOS.png'

    h2o, err_new_fit, sd, h2o_sd = plot_fit(wet_data, h2o, dry_means[wet_key], params, sd_data, h2o_sd, instcorrdata = instcorr_data, saveas = outdir + '/' + site+ '/' + dir_suffix+ outname, No_plots = No_plots, time=data.datetime)
    
    #return fit_params
    return params, np.sqrt(corr)[np.isfinite(np.sqrt(corr))], h2o, err_new_fit, sd, h2o_sd

def fitfunc(X1, a1, b1):
    f = lambda X, a, b : 1 + a*X + b*(X*X)
    return f(X1, a1, b1)
    
def ODRfunc(param, X1):
    a1 = param[0]
    b1 = param[1]
    f = lambda X, a, b : 1 + a*X + b*(X*X)
    return f(X1, a1, b1)
    
# syntax: CRDS_H2OCorr.plot_fit(wetdata[0].co2_wet, data_nocorr[0].co2_dry, wetdata[0].h2o, np.mean(drydata[0].co2_wet), fit, wetdata[0].co2sd)
def plot_fit(rawwetdata, h2o, drymean, fit_params, sd, h2o_sd, instcorrdata=None, saveas = None, colour = 'blue', No_plots = None, time=None):
    
    #f = 1 + a*X + b*(X*X)
    fit = fitfunc(h2o, fit_params[0], fit_params[1])    
    
    CO2_corr = rawwetdata/fit
    # err = 100 - (correct value/dry raw value * 100 )
    err_new_fit = CO2_corr - drymean
    if instcorrdata is not None:    
        err_inst_fit = instcorrdata- drymean
    
    if No_plots is None:
        # Plot the corrected output
        fig = plt.figure()
        fig.set_size_inches(4.2,4)
        fig.subplots_adjust(left = 0.2)
        fig.subplots_adjust(bottom = 0.2)
    
        # Plot CO2 wet/ CO2 dry vs H2O
        plt3 = plt.subplot(1,1,1)
        if instcorrdata is not None: 
            plt3.errorbar(h2o, err_inst_fit, yerr = sd, xerr = h2o_sd, fmt= 'o', color = 'red')

        plt3.errorbar(h2o, err_new_fit, yerr = sd, xerr = h2o_sd, fmt= 'o', color = colour)
     
        plt3.set_xlabel('H$_2$O raw')
        plt3.set_ylabel('corrected - dry mean')
        plt.figtext(0.25, 0.25, 'Dry conc = ' + str(drymean), verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        
        plt.figtext(0.25, 0.30, 'Instrument specific correction', verticalalignment='bottom', horizontalalignment='left', fontsize=10, color = colour)
        plt.figtext(0.25, 0.35, 'In built correction', verticalalignment='bottom', horizontalalignment='left', fontsize=10, color = 'red')
        
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        plt3.yaxis.set_major_formatter(y_formatter)
        plt3.set_xlim(0,4)     
        # Set y ranges
        if 'co2' in saveas:
            plt3.set_ylim(-0.6, 0.6)             
             
        if 'ch4' in saveas:
            plt3.set_ylim(-6, 6)             
             
        if saveas is not None:
            print 'Plot saved as: ' + saveas
            plt.savefig(saveas, dpi=200)
            
        plt.show()
        plt.close()    
    
        """
        # quick plot of the relationship between sd and residual
        if isinstance(time, list) :
            since = [(i -time[0]).total_seconds()/60 for i in time]
            plt.scatter(since, h2o, c=err_new_fit, lw =0.1)
            plt.ylabel('h2o')
            plt.xlabel('Minutes since start of run')
            plt.xlim()
            plt.savefig('/Users/as13988/Documents/Work/Picarro/H2OCorr/temp/'+ saveas.split('/')[-1])
            plt.show()
        """
    return h2o, err_new_fit, sd, h2o_sd

import colorsys
# Make a list of N distinct colours
def get_colours(num_colours):
    colours=[]
    for i in np.arange(0., 360., 360. / num_colours):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colours.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colours


# Uses the output of plot_ratio_multi
def plot_residual_multi(data, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', saveas = 'None'):
        
    f, (plt1, plt2) = plt.subplots(2, sharex=True)
    
    # Generate colours
    colours = get_colours(len(data['CO2']['h2o']))

    CO2 = data['CO2']
    CH4 = data['CH4']    
    
    for i in np.arange(len(data['CO2']['h2o'])):
        plt1.errorbar(CO2['h2o'][i], CO2['residual'][i], yerr=CO2['sd'][i],fmt= '.', color = colours[i])
        plt2.errorbar(CH4['h2o'][i], CH4['residual'][i], yerr=CH4['sd'][i],fmt='.', color = colours[i])

    plt2.set_xlabel('H$_2$O raw')
    plt1.set_ylabel('CO2')
    plt2.set_ylabel('CH4')

    plt.figtext(0.15, 0.82, 'Residual = |corrected - dry mean|', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
    


    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt1.yaxis.set_major_formatter(y_formatter)
    plt2.yaxis.set_major_formatter(y_formatter)
         
    if saveas is not None:
       print 'Plot saved as: ' + outdir + saveas
       plt.savefig(outdir+saveas, dpi=200)
        
    plt.show()
    plt.close()    
    

# Uses the "means" outputs of plot_ratio_multi and the wet and dry data as read in using read_data_multi
# The co2_dry output of GCWerks ≠ the co2_dry parameters from the CRDS if there's a water corrections file
# need to extract the data after removing the water corrections file to get the inbuilt correction
# and read it in as a second "wet_data" file
def plot_residual(means, wet_data, wet_data_nocorr, dry_data, outdir = '/Users/as13988/Documents/Work/Picarro/H2OCorr/Plots/', dir_suffix=None, saveas = None, useinbuilt = None, H2O_Range = None, sd_flag=None):
    
    # extract the wet data and concatenate it into one large data blob
    co2 = []
    co2_inbuilt = []
    co2_sd = []
    ch4 = []
    ch4_inbuilt = []
    ch4_sd = []
    h2o = []
    h2o_sd =[]
    flags = []

    if isinstance(wet_data, list):
        for i, item in enumerate(wet_data):
            co2.extend(item.co2_wet)
            ch4.extend(item.ch4_wet)
            co2_inbuilt.extend(wet_data_nocorr[i].co2_dry)
            ch4_inbuilt.extend(wet_data_nocorr[i].ch4_dry)
            
            ch4_sd.extend(item.ch4sd)    
            co2_sd.extend(item.co2sd)
            h2o_sd.extend(item.h2osd)
            
            flags.extend(item.co2flags)
        
            h2o.extend(item.h2o)
    
    
    # extract the dry data and calculate the mean
    co2_dry = []
    ch4_dry = []
    dry_flags = []
    if isinstance(dry_data, list):
        for i, item in enumerate(dry_data):
            co2_dry.extend(item.co2_wet)
            ch4_dry.extend(item.ch4_wet)
            dry_flags.extend(item.co2flags)

    # exclude flagged data
    dry_index = np.where(np.array(dry_flags) == 0)[0]
    co2_dry_mean = np.mean(np.array(co2_dry)[dry_index])
    ch4_dry_mean = np.mean(np.array(ch4_dry)[dry_index])

    index = np.where(np.array(flags) == 0)[0]
    print 'No. good flagged data points: ' + str(len(index))
    
    # optional SD flag
    if sd_flag is not None:
        index = np.where((np.array(flags) == 0) & (np.array(co2_sd) < sd_flag))[0]
    
    print 'No. data points used: ' + str(len(index))
            
    co2 = np.array(co2)[index]
    ch4 = np.array(ch4)[index]
    h2o = np.array(h2o)[index]
    co2_inbuilt = np.array(co2_inbuilt)[index]
    ch4_inbuilt = np.array(ch4_inbuilt)[index]
    co2_sd = np.array(co2_sd)[index]
    ch4_sd = np.array(ch4_sd)[index]
    h2o_sd = np.array(h2o_sd)[index]

    # calculate the corrected data and the residual
    co2_res = co2/(fitfunc(h2o, means['a_co2'], means['b_co2'])) - co2_dry_mean
    ch4_res = ch4/(fitfunc(h2o, means['a_ch4'], means['b_ch4'])) - ch4_dry_mean

    co2_res_inbuilt = co2_inbuilt - co2_dry_mean
    ch4_res_inbuilt = ch4_inbuilt - ch4_dry_mean
        
    # plot it
    f, (plt1, plt2) = plt.subplots(2, sharex=True)
    
    f.subplots_adjust(bottom = 0.2, left=0.15)    

    if useinbuilt is not None:
        plt1.errorbar(h2o, co2_res_inbuilt, yerr=co2_sd, xerr = h2o_sd, fmt= 'o', mec = 'grey', mfc = 'grey', ecolor ='grey', markersize =0.5)
        plt2.errorbar(h2o, ch4_res_inbuilt, yerr=ch4_sd, xerr = h2o_sd, fmt='o', mec = 'grey', mfc = 'grey', ecolor = 'grey', markersize =0.5)        
    
    plt1.errorbar(h2o, co2_res, yerr=co2_sd, xerr = h2o_sd,  fmt= 'o', mec = 'red', mfc = 'red', ecolor ='red', markersize =0.5)
    plt2.errorbar(h2o, ch4_res, yerr=ch4_sd, xerr = h2o_sd, fmt='o', mec = 'blue', mfc = 'blue', ecolor = 'blue', markersize =0.5)

    plt1.plot([0,4], [0,0], 'k-')
    plt2.plot([0,4], [0,0], 'k-')

    plt1.plot([0,4], [0.1,0.1], '--', color='black')
    plt2.plot([0,4], [2,2], '--', color='black')
        
    plt1.plot([0,4], [-0.1,-0.1], '--', color='black')
    plt2.plot([0,4], [-2,-2], '--', color='black')               

  
    # option to print means for a given range
    if H2O_Range is not None:
       index = np.where( h2o < H2O_Range) [0]
       index_rest = np.where( h2o >= H2O_Range)[0]
       
       print 'Residual mean H2O < ' + str(H2O_Range) + '%'
       for i,item in enumerate([index, index_rest]):
           co2_LT_mean_N = np.mean(np.abs(co2_res[item]))
           co2_LT_sd_N = np.std(np.abs(co2_res[item]))
           co2_LT_mean_O = np.mean(np.abs(co2_res_inbuilt[item]))
           co2_LT_sd_O = np.std(np.abs(co2_res_inbuilt[item]))
           
           ch4_LT_mean_N = np.mean(np.abs(ch4_res[item]))
           ch4_LT_sd_N = np.std(np.abs(ch4_res[item]))      
           ch4_LT_mean_O = np.mean(np.abs(ch4_res_inbuilt[item]))
           ch4_LT_sd_O = np.std(np.abs(ch4_res_inbuilt[item]))   
        
           print 'CO2 new mean ± 1SD = ' +  '{:.4f}'.format(co2_LT_mean_N) + ' ± ' + '{:.4f}'.format(co2_LT_sd_N) + 'mumol/mol'
           print 'CO2 inbuilt mean ± 1SD = ' +  '{:.4f}'.format(co2_LT_mean_O) + ' ± ' + '{:.4f}'.format(co2_LT_sd_O) + 'mumol/mol'
           print 'CH4 new mean ± 1SD = ' +  '{:.4f}'.format(ch4_LT_mean_N) + ' ± ' + '{:.4f}'.format(ch4_LT_sd_N) + 'nmol/mol'
           print 'CH4 inbuilt mean ± 1SD = ' +  '{:.4f}'.format(ch4_LT_mean_O) + ' ± ' + '{:.4f}'.format(ch4_LT_sd_N) + 'nmol/mol'
           print ' '
           if i == 0:
               print 'Residual mean H2O >= ' + str(H2O_Range) + '%'
       
    labelx = -0.1

    plt1.set_ylim([-0.6, 0.6])
    plt2.set_ylim([-6, 6])
    plt1.set_xlim([0,4])
    plt2.set_xlim([0,4]) 

    plt2.set_xlabel('H$_2$O raw')
    plt1.set_ylabel(r"[CO$_2$] residual""\n"r"($\mu$mol/mol)")
    plt1.yaxis.set_label_coords(labelx, 0.5)
    
    plt2.set_ylabel('[CH$_4$] residual\n(nmol/mol)')
    plt2.yaxis.set_label_coords(labelx, 0.5)

    plt.figtext(0.6, 0.03, 'Residual = corrected - dry mean', verticalalignment='bottom', horizontalalignment='left', fontsize=8)

    plt.figtext(0.63, 0.58, r'$\bar{x}$$\pm 1 \sigma = ' +  '{:.3f}'.format(np.mean(np.abs(co2_res))) + '\pm' + '{:.3f}'.format(np.std(np.abs(co2_res))) + '\,  \mu mol/mol$', verticalalignment='bottom', horizontalalignment='left', fontsize=8)
    plt.figtext(0.63, 0.2, r'$\bar{x}$$\pm 1 \sigma = ' +  '{:.3f}'.format(np.mean(np.abs(ch4_res))) + '\pm' + '{:.3f}'.format(np.std(np.abs(ch4_res))) + '\, nmol/mol$', verticalalignment='bottom', horizontalalignment='left', fontsize=8)

    if useinbuilt is not None:
        plt.figtext(0.63, 0.61, r'$\bar{x}$$\pm 1 \sigma = ' +  '{:.3f}'.format(np.mean(np.abs(co2_res_inbuilt))) + '\pm' + '{:.3f}'.format(np.std(np.abs(co2_res_inbuilt))) + '\,  \mu mol/mol$', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'grey')
        plt.figtext(0.63, 0.23, r'$\bar{x}$$\pm 1 \sigma = ' +  '{:.3f}'.format(np.mean(np.abs(ch4_res_inbuilt))) + '\pm' + '{:.3f}'.format(np.std(np.abs(ch4_res_inbuilt))) + '\, nmol/mol$', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'grey')
    

    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt1.yaxis.set_major_formatter(y_formatter)
    plt2.yaxis.set_major_formatter(y_formatter)
         
    if saveas is not None:
        outname = '/' +means['site'] + '/'+ dir_suffix+ means['site'] + '_' + means['date']+'_RunResiduals.png'
        if useinbuilt is not None:
            outname = '/' +means['site'] + '/'+ dir_suffix+ means['site'] + '_' + means['date']+'_RunResiduals_InBuiltCorr.png'
            
        print 'Plot saved as: ' + outdir +  outname
        plt.savefig(outdir+  outname, dpi=200)
        
    plt.show()
    plt.close()
    
    """
    # Quick plot of the relationship between the residuals and the sd of h2o, co2 and ch4.
    plt.scatter(co2_sd, h2o_sd, c=co2_res, lw=0)
    plt.xlabel('CO2 sd')
    plt.ylabel('H2O sd')
    plt.ylim([0, 0.2]) 
    #plt.xlim([0, 0.06])             
    plt.savefig('/Users/as13988/Documents/Work/Picarro/H2OCorr/temp/CO2'+outname.split('/')[-1])
    plt.show()

    plt.scatter(ch4_sd, h2o_sd, c=ch4_res, lw=0)
    plt.xlabel('CH4 sd')
    plt.ylabel('H2O sd')
    plt.ylim([0, 0.2]) 
    #plt.xlim([0, 1]) 
    plt.savefig('/Users/as13988/Documents/Work/Picarro/H2OCorr/temp/CH4'+outname.split('/')[-1])
    plt.show()
    """
    
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
    
    colours = ['blue', 'green', 'purple', 'orange', 'pink', 'magenta']

    # PLot the fit  of the new correction compared to the in build correction

    # calculate parameters
    err_inst_fit = np.abs(dry_data-drydata[0])

    
    # Plot the corrected output
    fig = plt.figure()
    fig.set_size_inches(4.2,4)
    fig.subplots_adjust(left = 0.2)
    fig.subplots_adjust(bottom = 0.15)
    
    # Plot CO2 wet/ CO2 dry vs H2O
    plt3 = plt.subplot(1,1,1)
    
    plt3.errorbar(h2o, err_inst_fit, yerr = sd_data, fmt= 'o', color = 'red')

    nofits = len(fit_params)
    
    if isinstance(fit_params, list) == False:
        nofits = 1

    for i in np.arange(nofits):
        fit_params_i = fit_params[i]
        
        if isinstance(fit_params, list) == False:
            fit_params_i = fit_params
            
        fit = fitfunc(h2o, fit_params_i[0], fit_params_i[1])       
        
        CO2_corr = wet_data/fit
        
        # err = 100 - (correct value/dry raw value * 100 )
        err_new_fit = np.abs(CO2_corr-drydata[i])
        
        plt3.errorbar(h2o, err_new_fit, yerr = sd_data, fmt= 'o', color = colours[i])
        
        plt.figtext(0.22, 0.2 +i*0.04, 'New correction ' +str(i+1), color = colours[i], verticalalignment='bottom', horizontalalignment='left', fontsize=10)


    plt3.set_xlabel('H$_2$O raw')
    plt3.set_ylabel('|corrected - dry mean|')
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    plt3.yaxis.set_major_formatter(y_formatter)
    
    plt.figtext(0.22, 0.16, 'Built in correction', color ='red', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
        
    
    print 'Plot saved as: ' + outdir + '/' + site+ '/'+ 'FitComparison_' + tag+ '.' + site + testdate +'.'+ species + '.png'
    plt.savefig(outdir+ '/' + site + '/' + 'FitComparison_' + tag+ '.' + site + testdate + species + '.png', dpi=200)
    plt.show()
    
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
 
    print 'Plot saved as: ' + outdir+ '/' + site+ '/'+'RatioComparison'+site+'.' + testdate+ '.'+ species + '.png'
    plt.savefig(outdir+ '/' + site+ '/'+'RatioComparison'+site+'.' + testdate+ '.'+ species + '.png', dpi=200)
        
    plt.show()
    plt.close()


def calc_dry(data, ignoreindex=0):
    index = np.where(data.co2flags == 0)[0]    
    
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
 
    print 'Plot saved as: ' + outdir+ '/' + site+ '/'+'drymeans.png'
    plt.savefig(outdir+ '/' + site+ '/'+'drymeans.png', dpi=200)
        
    plt.show()
    
    return 


def responsefn(H2O, A, B):
    
    ratio = 1 + A*H2O + B*H2O*H2O
    
    return ratio