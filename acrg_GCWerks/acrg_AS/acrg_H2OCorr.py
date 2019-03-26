# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:47:23 2015

@author: as13988
"""

from acrg_GCWerks import acrg_read_GCwerks as read_GCwerks
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib
import pdb
import numpy as np
import datetime as dt
import os
import bisect
from acrg_GCWerks import fitting
import glob
 
# Read raw data
class read_raw:
    def __init__(self, inputfile='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/DPG.E-114.DPG=5a.20150325.dat'):
        data = read_GCwerks.read_gcexport_crds(inputfile)
        
        self.data = data


# Read timefile
"""
# Timefile dew point generator tests run at Heathfield
# Configurations types
# 1 = H2O trap bypassed, through Nafion
# 2 = H2O trap bypassed, Nafion bypassed
# 3 = through H2O trap, Nafion bypassed
# 4 = through H2O trap, through Nafion
#
# Runs
# 1 = Cylinder E-114 Dew Point = 5C
# 2 = Cylinder E-114 Dew Point = 5C
# 3 = Cylinder E-114 Dew Point = 15C
# 4 = Cylinder E-114 Dew Point = 20C
# 5 = Cylinder of compressed air Dew Point = 5
# 6 = Cylinder of compressed air Dew Point = 15
# 7 = Cylinder H-249 Dew Point = 5
# 8 = Cylinder H-249 Dew Point = 10
# 9 = Cylinder UAN20070098 Dew Point = 15
# 10 = Cylinder UAN20070098 Dew Point = 5
#---------------------------------------
Cylinder	DewPoint	Config	RunNo	StartTime			EndTime				Comment
E-114		5 			2 		1 		11/03/2015 13:32		11/03/2015 16:05		# Could be dodgy as it's at a different flow rate
"""
class read_timefile:
    def __init__(self, timefile='/Users/as13988/Documents/Work/Picarro/Drying/DPGTimefile'):
        
        config_types = ['H2O trap bypassed, through Nafion', \
                         'H2O trap bypassed, Nafion bypassed', \
                         'through H2O trap, Nafion bypassed', \
                         'through H2O trap, through bypassed']

        run_types = ['Cylinder E-114 Dew Point = 5C', \
                    'Cylinder E-114 Dew Point = 5C', \
                    'Cylinder E-114 Dew Point = 15C', \
                    'Cylinder E-114 Dew Point = 20C', \
                    'Cylinder of compressed air Dew Point = 5', \
                    'Cylinder of compressed air Dew Point = 15', \
                    'Cylinder H-249 Dew Point = 5', \
                    'Cylinder H-249 Dew Point = 10', \
                    'Cylinder UAN20070098 Dew Point = 15', \
                    'Cylinder UAN20070098 Dew Point = 5']
        
        # If you leave the delimiter blank it splits based on white space
        if type(timefile) == tuple:
            timedata=np.genfromtxt(timefile[0], dtype=str, skip_header=20)
        elif type(timefile) == str:
            timedata=np.genfromtxt(timefile, dtype=str, skip_header=20)     
        
        cylinder = timedata[:,0]
        
        dewpoint = timedata[:,1].astype(int)      
        config_no = timedata[:,2].astype(int)
        run_no = timedata[:,3]
        starttime = [timedata[i,4] + ' ' + timedata[i,5] for i in np.arange(len(timedata[:,4]))]
        endtime = [timedata[i,6] + ' ' + timedata[i,7] for i in np.arange(len(timedata[:,4]))]
        
        # Convert starttime and endtime to a datetime 
        start_dt = [dt.datetime.strptime(i, "%d/%m/%Y %H:%M") for i in starttime]
        end_dt = [dt.datetime.strptime(i, "%d/%m/%Y %H:%M") for i in endtime]
        
        self.cylinder = cylinder
        self.dewpoint = dewpoint.astype(int)
        self.config_no = config_no.astype(int)
        self.run_no = run_no.astype(int)  
        self.starttime_str = starttime
        self.sendtime_str = endtime
        self.start_dt = start_dt
        self.end_dt = end_dt
        self.config_types = config_types
        self.run_types = run_types
                
        dirname, filename = os.path.split(timefile)
        
        self.filename = filename
        self.dirname = dirname


# Split minutemeans
class split_raw:
    def __init__(self, data, timefile='/Users/as13988/Documents/Work/Picarro/Drying/DPGTimefile', run_no = 0):
     
        data = data.data
     
     
        # Extract the unflagged data
        dt_co2, co2, co2sd = read_GCwerks.Extractgood(data.datetime, data.co2, data.co2flags, data.co2sd)
        dt_ch4, ch4, ch4sd = read_GCwerks.Extractgood(data.datetime, data.ch4, data.ch4flags, data.ch4sd)
        dt_h2o, h2o, dummy = read_GCwerks.Extractgood(data.datetime, data.h2o, data.co2flags)
  
        # Read in the timefile
        timedata = read_timefile(timefile=timefile)     
        
        # If the run number isn't given default to the run that is closest to the first time stamp in the data
        if run_no == 0:
            config_no_sorted = [x for (y,x) in sorted(zip(timedata.start_dt,timedata.run_no))]
            start_dt_sorted = [y for (y,x) in sorted(zip(timedata.start_dt,timedata.run_no))]
        
            run_no = config_no_sorted[bisect.bisect(start_dt_sorted, data.datetime[0])]
        
        #pdb.set_trace()

        # Find the start and end times for that run
        index = np.where(timedata.run_no == run_no)[0]
        starttimes = [timedata.start_dt[i] for i in index]
        endtimes = [timedata.end_dt[i] for i in index]
        config_no = timedata.config_no[index]
                
        
        #pdb.set_trace()
        
        # Split the time file info for that run no into the individual config types
        # Extract the corresponding data and time stamps
        config_types = list(set(config_no))
         
        datasplit = []         
        
        # go through each config type
        for j in config_types:
            config_index = np.where(config_no == j)[0]
                     
            configsplit = []            
            
            for k in config_index:
                startindex = bisect.bisect(dt_co2, starttimes[k])
                endindex = bisect.bisect(dt_co2, endtimes[k])
                
                co2times_k = dt_co2[startindex:endindex]
                co2_k = co2[startindex:endindex]
                co2sd_k = co2sd[startindex:endindex]
                h2o_k = h2o[startindex:endindex]                
                h2otimes_k = dt_h2o[startindex:endindex]
                
                startindex = bisect.bisect(dt_ch4, starttimes[k])
                endindex = bisect.bisect(dt_ch4, endtimes[k])
                
                ch4times_k = dt_ch4[startindex:endindex]
                ch4_k = ch4[startindex:endindex]
                ch4sd_k = ch4sd[startindex:endindex]
                        
                
                dict_k = {'dt_co2' :co2times_k, \
                        'co2' : co2_k, \
                        'co2sd' : co2sd_k, \
                        'dt_ch4' :ch4times_k, \
                        'ch4' : ch4_k, \
                        'ch4sd' : ch4sd_k, \
                        'dt_h2o' : h2otimes_k, \
                        'h2o' : h2o_k}
                
                configsplit.append(dict_k)
                
                
            dict_j = {'cylinder' : timedata.cylinder[index[0]], \
                    'dewpoint' : timedata.dewpoint[index[0]], \
                    'config_no' : j, \
                    'run_no' : run_no, \
                    'starttime' : [starttimes[i] for i in config_index], \
                    'endtime' : [endtimes[i] for i in config_index], \
                    'data' : configsplit, \
                    }
            datasplit.append(dict_j)
            
        self.split_data = datasplit   
        self.unflaggeddata = {'dt_co2' :dt_co2, \
                        'co2' : co2, \
                        'co2sd' : co2sd, \
                        'dt_ch4' :dt_ch4, \
                        'ch4' : ch4, \
                        'ch4sd' : ch4sd, \
                        'dt_h2o' : dt_h2o, \
                        'h2o' : h2o}
        self.timedata = timedata
        self.configs = config_types


# do fits and calculate means
class make_means:
    def __init__(self, inputfile='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/DPG.E-114.DPG=5a.20150326.dat', \
            timefile='/Users/as13988/Documents/Work/Picarro/Drying/DPGTimefile', \
            outputdir='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/', run_no = 0, species ='co2', fit_type=0):
     
        # Read in the raw data 
        data = read_raw(inputfile=inputfile)
     
        # split the raw data
        processed_data = split_raw(data, timefile=timefile, run_no = run_no)
        
        # all data
        all_data = processed_data.unflaggeddata
        all_seconds = np.array([(l - all_data['dt_co2'][0]).total_seconds() for l in all_data['dt_co2']])
            
# ____________________________________________________________________
# FITS        
# ____________________________________________________________________
# TYPE 2 = do a fit to the "No H2O trap, No nafion"        
        type2_data = processed_data.split_data[np.where(np.array(processed_data.configs) == 2)[0]]

        # cycle through each run of type 2   
        x_2 = []
        y_2 = []
        h2o_2 = []
        
        no_runs_2 = np.shape(type2_data['data'])[0]
        for j in np.arange(no_runs_2):
            x_2.extend(type2_data['data'][j]['dt_'+species])          
            y_2.extend(type2_data['data'][j][species])
            h2o_2.extend(type2_data['data'][j]['h2o'])
        
        
        # Calculate the mean h2o value for the run
        h2o_mean_2 = np.mean(h2o_2)        
        h2o_sd_2 = np.std(h2o_2)        
        
        # convert to seconds since the start of the run
        x_2_secs = np.array([(k - all_data['dt_'+species][0]).total_seconds() for k in x_2])
        
        # do the fit to the data:
        # If there are multiple runs the fit to exp/log fit       
        # If the individual run contains > 30mins worth of data then ignore the first 10mins 
        # (this should be done by the DPGTimeFile) and then do a exp/log
        # If the individual run contains < 30mins of data then just do a linear fit
        #pdb.set_trace()        
        if no_runs_2 == 1:
            if np.shape(x_2_secs)[0] > 30:            
                fit_2 = fitting.fit_data(x_2_secs, y_2, scale=1, fit_type=fit_type)
            else:
                fit_2 = fitting.fit_data(x_2_secs, y_2, scale=1, fit_type='lin')
                
        else:
            fit_2 = fitting.fit_data(x_2_secs, y_2, scale=1, fit_type=fit_type)

        # calculate the fit at all the data points
        all_fit_2 = (fitting.func(all_seconds*fit_2.x_scale, fit_2.coeffs, fit_type=fit_2.fit_type))/fit_2.y_scale
        
        
# TYPE 3 = do a fit to the "H2O trap, No nafion"        
        type3_data = processed_data.split_data[np.where(np.array(processed_data.configs) == 3)[0]]

        # cycle through each run of type 3   
        x_3 = []
        y_3 = []
        h2o_3 = []
        
        no_runs_3 = np.shape(type3_data['data'])[0]
        for j in np.arange(no_runs_3):
            x_3.extend(type3_data['data'][j]['dt_'+species])          
            y_3.extend(type3_data['data'][j][species])
            h2o_3.extend(type3_data['data'][j]['h2o'])
            
        # Calculate the mean h2o value for the run
        h2o_mean_3 = np.mean(h2o_3)        
        h2o_sd_3 = np.std(h2o_3)        


        # convert to seconds since the start of the run
        x_3_secs = np.array([(k - all_data['dt_co2'][0]).total_seconds() for k in x_3])
        
        # do the fit to the data:
        # If there are multiple runs the fit to exp/log fit       
        # If the individual run contains > 30mins worth of data then ignore the first 10mins 
        # (this should be done by the DPGTimeFile) and then do a exp/log
        # If the individual run contains < 30mins of data then just do a linear fit
        if no_runs_3 == 1:
            if np.shape(x_3_secs)[0] > 30:            
                fit_3 = fitting.fit_data(x_3_secs, y_3, scale=1, fit_type=fit_type)
            else:
                fit_3 = fitting.fit_data(x_3_secs, y_3, scale=1, fit_type='lin')
                
        else:
            fit_3 = fitting.fit_data(x_3_secs, y_3, scale=1, fit_type=fit_type)

        # calculate the fit at all the data points
        all_fit_3 = (fitting.func(all_seconds*fit_3.x_scale, fit_3.coeffs, fit_type=fit_3.fit_type))/fit_3.y_scale


# ____________________________________________________________________
# EXTRACT OTHER DATA     
# ____________________________________________________________________
# TYPE 1 - "No H2O trap, nafion"    
        # cycle through each run of type 1 
        x_1 = []
        y_1 = [] 
        h2o_1 = []
        type1_data = processed_data.split_data[np.where(np.array(processed_data.configs) == 1)[0]]
         
        for j in np.arange(np.shape(type1_data['data'])[0]):
            x_1.extend(type1_data['data'][j]['dt_'+species])          
            y_1.extend(type1_data['data'][j][species])
            h2o_1.extend(type1_data['data'][j]['h2o'])
            
        # Calculate the mean h2o value for the run
        h2o_mean_1 = np.mean(h2o_1)        
        h2o_sd_1 = np.std(h2o_1)                                
            
        # convert to seconds since the start of the run
        x_1_secs = np.array([(k - all_data['dt_co2'][0]).total_seconds() for k in x_1])
        
        # Calculate values from the fit to type 2 at type 1's time steps
        x_1_fit = (fitting.func(x_1_secs*fit_2.x_scale, fit_2.coeffs, fit_type=fit_2.fit_type))/fit_2.y_scale        
               
            
# TYPE 4 - "H2O trap, nafion"       
        # cycle through each run of type 4
        x_4 = []
        y_4 = []
        h2o_4 = []
        type4_data = processed_data.split_data[np.where(np.array(processed_data.configs) == 4)[0]]
        
        for j in np.arange(np.shape(type4_data['data'])[0]):
            x_4.extend(type4_data['data'][j]['dt_'+species])          
            y_4.extend(type4_data['data'][j][species])
            h2o_4.extend(type4_data['data'][j]['h2o'])
        
        # Calculate the mean h2o value for the run
        h2o_mean_4 = np.mean(h2o_4)        
        h2o_sd_4 = np.std(h2o_4)                
        
        x_4_secs = np.array([(k - all_data['dt_'+species][0]).total_seconds() for k in x_4])
        
        # Calculate values from the fit to type 3 at type 4's time steps
        x_4_fit = (fitting.func(x_4_secs*fit_3.x_scale, fit_3.coeffs, fit_type=fit_3.fit_type))/fit_3.y_scale        
        

# ____________________________________________________________________
# PLOTS       
# ____________________________________________________________________
# Fits to raw data    
        colours = ['purple','blue','green','red'] 
        
    # PLot fit to config 2
        fig = plt.figure()
        ax = fig.add_subplot(111)        
        fig.subplots_adjust(bottom = 0.2)

        # Plot all data
        ax.plot(all_seconds, all_data[species], '+', color = 'grey')

        # Plot type 2 data
        ax.plot(x_2_secs,y_2, 'o', color = colours[1])

        # Plot fit to type 2 data           
        ax.plot(all_seconds, all_fit_2, '-', color = colours[1])

        # Plot type 3 data
        ax.plot(x_3_secs,y_3, 's', color = colours[2])
        
        # Plot fit to type 3 data           
        ax.plot(all_seconds, all_fit_3, '-', color = colours[2])

        # PLot type 1, 3 & 4 data
        ax.plot(x_1_secs,y_1, 'o', color = colours[0])
        ax.plot(x_4_secs,y_4, 's', color = colours[3])
     
        plt.figtext(0.6, 0.05, 'H2O trap, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'green')
        plt.figtext(0.6, 0.1, 'H2O trap, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'red')
        plt.figtext(0.2, 0.1, 'H2O trap bypassed, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'purple')
        plt.figtext(0.2, 0.05, 'H2O trap bypassed, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'blue')
        
        plt.figtext(0.25, 0.5, fit_2.fit_label, color = 'blue')
        plt.figtext(0.25, 0.45, fit_3.fit_label, color = 'green')

        ax.set_xlabel('Seconds since start of run')
        ax.set_ylabel('[' + species +'] dry minute mean')
        ax.set_title(os.path.splitext(data.data.filename)[0])
        
        if 'DPG.H-249.DPG=5' in os.path.splitext(data.data.filename)[0]:
            ax.set_ylim((399,402))

        if'DPG.UAN' in os.path.splitext(data.data.filename)[0]:
            ax.set_ylim((439,440.5))

        y_formatter = ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
         
        plt.show()
        plt.savefig(outputdir+'fits/'+os.path.splitext(data.data.filename)[0]+'_Fits.png', dpi=200)
        
        plt.close()       
        
    # PLot data with fit deducted
        fig = plt.figure()
        ax = fig.add_subplot(111)   
        fig.subplots_adjust(bottom = 0.2)

        # Plot type 2 less fit
        ax.plot(x_2_secs, y_2- fit_2.fit, 'o', mfc = 'none', mec = colours[1])

        # Plot type 3 less fit
        ax.plot(x_3_secs, y_3 - fit_3.fit, 's', mfc = 'none', mec = colours[2])   
        
        # Plot type 1 data less fit to type 2
        ax.plot(x_1_secs, y_1 - x_1_fit, 'o', mfc = 'none', mec =  colours[0])    
        
        # Plot type 4 data less fit to type 3
        ax.plot(x_4_secs, y_4 - x_4_fit, 's', mfc = 'none', mec = colours[3])   
                
        plt.figtext(0.2, 0.05, 'H2O trap bypassed, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'purple')
        plt.figtext(0.2, 0.1, 'H2O trap bypassed, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'blue')
        plt.figtext(0.6, 0.1, 'H2O trap, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'green')
        plt.figtext(0.6, 0.05, 'H2O trap, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'red')

        ax.set_xlabel('Seconds since start of run')
        ax.set_ylabel('Difference from long term trend (No trap, No Nafion)')
        ax.set_title(os.path.splitext(data.data.filename)[0])

       
        plt.show()
        plt.savefig(outputdir+'fits/'+os.path.splitext(data.data.filename)[0]+'_DiffFits.png', dpi=200)
               

# ____________________________________________________________________
# CALCULATE MEANS       
# ____________________________________________________________________

# Cycle through each config and deduct the config2 fit and calculate means for each run
        means = []
        
        configs = [1,2,3,4]
        
        for i in configs:
            #print i
            typei_data = processed_data.split_data[np.where(np.array(processed_data.configs) == i)[0]]
            
            no_means = np.shape(typei_data['data'])[0]
            seconds = np.zeros(no_means)
            conc = np.zeros(no_means)            
            conc_sd = np.zeros(no_means)
            h2o = np.zeros(no_means)
            h2o_sd = np.zeros(no_means)
            fit = np.zeros(no_means)
            fit_sd = np.zeros(no_means)
            diff = np.zeros(no_means)
            diff_sd = np.zeros(no_means)
            
            for j in np.arange(np.shape(typei_data['data'])[0]):
                # Extract the data for that block
                x_i = typei_data['data'][j]['dt_'+species] 
                y_i = typei_data['data'][j][species]
    
                h2o_i = typei_data['data'][j]['h2o']
                
                # Determine the seconds for the block
                x_i_secs = np.array([(k - all_data['dt_'+species][0]).total_seconds() for k in x_i])
        
                # Calculate values from the fit to type 2 (or type 3) at given time steps
                if i < 3:
                    x_i_fit = (fitting.func(x_i_secs*fit_2.x_scale, fit_2.coeffs, fit_type=fit_2.fit_type))/fit_2.y_scale        
       
                if i > 2 :       
                    x_i_fit = (fitting.func(x_i_secs*fit_3.x_scale, fit_3.coeffs, fit_type=fit_3.fit_type))/fit_3.y_scale        
       
               # Calculate means
                seconds[j] = np.mean(x_i_secs)
                conc[j] = np.mean(y_i)
                conc_sd[j] = np.std(y_i)
                h2o[j] = np.mean(h2o_i)
                h2o_sd[j] = np.std(h2o_i)
                fit[j] = np.mean(x_i_fit)
                fit_sd[j] =  np.std(x_i_fit)
                diff[j] = np.mean(y_i - x_i_fit)
                diff_sd[j] = np.std(y_i - x_i_fit)
       
       
            # store means
            mean_i = { 'config_no' : i, \
                      'seconds' : seconds, \
                      'conc' : conc, \
                      'conc_sd' : conc_sd, \
                      'h2o' : h2o, \
                      'h2o_sd' : h2o_sd, \
                      'fit' : fit, \
                      'fit_sd' : fit_sd, \
                      'diff' : diff, \
                      'diff_sd' : diff_sd}
            
            means.append(mean_i)


# ____________________________________________________________________
# PLOT MEANS       
# ____________________________________________________________________

# Cycle through each config and plot the means for each run
        #pdb.set_trace()
        for i in np.arange(np.shape(means)[0]):
            #print i
            #print means[i]['seconds'], means[i]['diff'], means[i]['diff_sd']
           # print colours[i]
            
            ax.errorbar(means[i]['seconds'], means[i]['diff'], means[i]['diff_sd'], color = colours[i], fmt='o')

        plt.savefig(outputdir+'fits/'+os.path.splitext(data.data.filename)[0]+'_DiffFit2_means.png', dpi=200)
        plt.close()

        self.all_seconds = all_seconds
        self.conc = all_data[species]
        self.means = means
        self.data = data
        self.processed_data = processed_data
        self.h2o_means = [h2o_mean_1, h2o_mean_2, h2o_mean_3, h2o_mean_4]
        self.h2o_sd = [h2o_sd_1, h2o_sd_2, h2o_sd_3, h2o_sd_4]
                
# Configurations types
# 1 = H2O trap bypassed, through Nafion
# 2 = H2O trap bypassed, Nafion bypassed
# 3 = through H2O trap, Nafion bypassed
# 4 = through H2O trap, through Nafion
# calc run means and H2O conc    
class calc_h2o_vs_means:
    def __init__(self, cylinderno):


        # find the means for the given cylinder number
        files = glob.glob('/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/*' + cylinderno + '*20150326*dat')

        h2o = []
        h2o_sd = []
        diffconc = []
        diffconc_sd = []

        # go through each file and extract the differences and h2o concs
        for j in files:
            means = make_means(inputfile=j)
 
            # Extract diff means for configs 1 and 4
            # The diff means for configs 2 and 3 will be 0 as they're the data that's fit to
            # Diff for 1 = nafion diff wet
            # Diff for 4 = nafion diff dry
            diffconc_1 = means.means[0]['diff']
            diffconc_sd_1 = means.means[0]['diff_sd']            

            diffconc_4 = means.means[3]['diff']
            diffconc_sd_4 = means.means[3]['diff_sd']            
            
            # Combine the dry and then wet values
            diffconc.extend(diffconc_1)    
            diffconc.extend(diffconc_4)  
            diffconc_sd.extend(diffconc_sd_1)
            diffconc_sd.extend(diffconc_sd_4)            
            
            
            # For config 1 and 2 we want to use the H2O concentration from config 2
            h2o_1 = np.zeros(len(diffconc_1))
            h2o_sd_1 = np.zeros(len(diffconc_1))

            h2o_1[:] = means.h2o_means[1]
            h2o_sd_1[:] = means.h2o_sd[1]
            
            # For config 4 we want to use the H2O concentration from config 3
            h2o_4 = np.zeros(len(diffconc_4))
            h2o_sd_4 = np.zeros(len(diffconc_4))
 
            h2o_4[:] = means.h2o_means[2]
            h2o_sd_4[:] = means.h2o_sd[2]
              
            
            
            # Combine the dry and then wet values
            h2o.extend(h2o_1)
            h2o.extend(h2o_4)  
            h2o_sd.extend(h2o_sd_1)            
            h2o_sd.extend(h2o_sd_4)

        self.cylinder_no = cylinderno
        self.h2o = h2o
        self.diffconc = diffconc
        self.diffconc_sd = diffconc_sd




# plot multiple run means and H2O conc    
class plot_h2o_vs_means_multi:
    def __init__(self, cylindernos, outputdir='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/', filesuffix='', close = 0):

        opennew = 0
        
        for i in cylindernos:
            plot_h2o_vs_means(i, close = 1, opennew = opennew)
            opennew = 1            

        colours = ['orange', 'aqua', 'magenta', 'lime']
        cylinders = np.array(['E-114', 'H-249', 'UAN20070098', 'IndustAir'])        
        
        
        plt.figtext(0.1, 0.05, cylinders[0], color = colours[0])
        plt.figtext(0.1, 0.08, cylinders[1], color = colours[1])
        plt.figtext(0.3, 0.05, cylinders[2], color = colours[2])
        plt.figtext(0.3, 0.08, cylinders[3], color = colours[3])
        
        plt.show()        
        plt.savefig(outputdir+'All_DiffVsH2O.png', dpi=200)
        plt.close()


# plot run means and H2O conc    
class plot_h2o_vs_means:
    def __init__(self, cylinderno, outputdir='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/', filesuffix='', close = 0, opennew=0):

        data = calc_h2o_vs_means(cylinderno)
        
        colours = ['orange', 'aqua', 'magenta', 'lime']
        cylinders = np.array(['E-114', 'H-249', 'UAN20070098', 'IndustAir'])
        
        colour = colours[np.where(cylinders == cylinderno)[0]]
        
        if opennew == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)   
            fig.subplots_adjust(bottom = 0.2)
            ax.set_xlabel('H2O concentration (%) at the nafion')
            ax.set_ylabel('With Nafion less Without Nafion (ppm)')
            ax.set_xlim((-0.5,2.5))
            ax.set_ylim((-0.3,0.2))
            


        # Plot data
        plt.errorbar(data.h2o, data.diffconc, data.diffconc_sd, color = colour, fmt = 'o')
       
        

        # Do a fit
        fit = fitting.fit_data(data.h2o, data.diffconc, fit_type = 'lin', sigma = data.diffconc_sd)
        
        # Generate some dummy values for x
        x =np.arange(30)
        x = x/10.0
        x =x -0.5        
        x_fit = fitting.func(x, fit.coeffs, fit_type=fit.fit_type)        
        
        # Plot fit
        plt.plot(x, x_fit, '-', color = colour)
        #pdb.set_trace()

        # PLot zero lines    
        plt.plot([0,0],[-0.3,0.2], '--', color = 'black')
        plt.plot([-0.5,2.5],[0,0], '--', color = 'black')
        
        
        if close == 0:
            plt.figtext(0.25, 0.25, fit.fit_label, color = colour)
            ax.set_title(cylinderno)            
            plt.show()        
            
            plt.savefig(outputdir+cylinderno+'_DiffVsH2O.png', dpi=200)
            plt.close()

        self.fit = fit
        self.x = x
     

# Plot minutemeans
class plot_raw:
    def __init__(self, inputfile = '/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/DPG.E-114.DPG=5a.20150325.dat',\
        outputdir='/Users/as13988/Documents/Work/Picarro/Drying/GCWerksOutput/', \
        timefile='/Users/as13988/Documents/Work/Picarro/Drying/DPGTimefile', plot_unflagged = 1, filesuffix=''):
        
        data = read_raw(inputfile=inputfile)
  
        #pdb.set_trace()  
  
        splitdata = split_raw(data, timefile=timefile, run_no = 0)
        
        run_no = splitdata.split_data[0]['run_no']
        
        # plot all the data in grey
        fig = plt.figure()
        fig.subplots_adjust(bottom = 0.2)
        
        plt1 = plt.subplot(2, 1, 1)
        
        # Config colours
        # set these to be consistent with Georgina
        # 1 = H2O trap bypassed, through Nafion PURPLE
        # 2 = H2O trap bypassed, Nafion bypassed BLUE
        # 3 = through H2O trap, Nafion bypassed GREEN
        # 4 = through H2O trap, through Nafion RED
        colours = ['purple','blue','green','red','grey']        
        
        
        if plot_unflagged == 0:
            # Plot CO2
            plt1.errorbar(data.data.datetime, data.data.co2, yerr=data.data.co2sd, fmt='+', color = 'grey', markersize = 3)
            #plt1.set_ylim(Calcrange(means.co2, means.co2sd))
            #plt1.set_xlim(daterange)     
            
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
            
            #x_formatter = matplotlib.dates.DateFormatter('%H:%M %m/%Y')
            x_tickno_formatter = ticker.MaxNLocator(4)
            x_formatter = matplotlib.dates.DateFormatter('%H:%M %d/%m/%Y')
            plt1.xaxis.set_major_formatter(x_formatter)
            plt1.xaxis.set_major_locator(x_tickno_formatter)
            
            plt1.set_title(splitdata.timedata.run_types[run_no -1])
            plt1.set_ylabel('[CO$_2$] (ppm)')
            
            # Plot CH4 
            plt2 = plt.subplot(2, 1, 2)
            plt2.errorbar(data.data.datetime, data.data.ch4, yerr=data.data.ch4sd, fmt='+', color = 'grey', markersize = 3)
            #plt2.set_ylim(Calcrange(means.ch4, means.ch4sd))
            plt2.set_ylabel('[CH$_4$] (ppb)')
            #plt2.set_xlim(daterange)  
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.xaxis.set_major_formatter(x_formatter)            
            plt2.xaxis.set_major_locator(x_tickno_formatter)
        
        else:
            #PLot CO2 data
            plt1.errorbar(splitdata.unflaggeddata['dt_co2'], splitdata.unflaggeddata['co2'], yerr=splitdata.unflaggeddata['co2sd'], fmt='+', color = 'grey', markersize = 3)
            #plt1.set_ylim(Calcrange(means.co2, means.co2sd))
            #plt1.set_xlim(daterange)     
            
            y_formatter = ticker.ScalarFormatter(useOffset=False)
            plt1.yaxis.set_major_formatter(y_formatter)
            
            x_formatter = matplotlib.dates.DateFormatter('%H:%M %m/%Y')
            x_tickno_formatter = ticker.MaxNLocator(4)
            
            plt1.xaxis.set_major_formatter(x_formatter)
            plt1.xaxis.set_major_locator(x_tickno_formatter)
            
            plt1.set_title(splitdata.timedata.run_types[run_no -1])
            plt1.set_ylabel('[CO$_2$] (ppm)')
            
            # Plot CH4
            plt2 = plt.subplot(2, 1, 2)
            plt2.errorbar(splitdata.unflaggeddata['dt_ch4'], splitdata.unflaggeddata['ch4'], yerr=splitdata.unflaggeddata['ch4sd'], fmt='+', color = 'grey', markersize = 3)
            #plt2.set_ylim(Calcrange(means.ch4, means.ch4sd))
            plt2.set_ylabel('[CH$_4$] (ppb)')
            #plt2.set_xlim(daterange)  
            plt2.yaxis.set_major_formatter(y_formatter)
            plt2.xaxis.set_major_formatter(x_formatter)            
            plt2.xaxis.set_major_locator(x_tickno_formatter)
        
      
        # plot each config over the top
        # Determine the number of configurations used
        configno = np.shape(splitdata.split_data)[0]
        
        for i in np.arange(configno):
 
            config_colour = colours[splitdata.split_data[i]['config_no'] -1]           
            configdata_i = splitdata.split_data[i]['data']

            # Determine the number of runs of configuration i
            no_config_runs = np.shape(configdata_i)[0]            
            
        
            
            for j in np.arange(no_config_runs):
                
                plt1.errorbar(configdata_i[j]['dt_co2'], configdata_i[j]['co2'], yerr=configdata_i[j]['co2sd'], fmt='+', color = config_colour)
                y_formatter = ticker.ScalarFormatter(useOffset=False)
                plt1.yaxis.set_major_formatter(y_formatter)
                  
                x_formatter = matplotlib.dates.DateFormatter('%H:%M %d/%m/%Y')
                x_tickno_formatter = ticker.MaxNLocator(4)
                
                  
                plt1.xaxis.set_major_formatter(x_formatter)
                plt1.xaxis.set_major_locator(x_tickno_formatter)
    
                plt.errorbar(configdata_i[j]['dt_ch4'], configdata_i[j]['ch4'], yerr=configdata_i[j]['ch4sd'], fmt='+', color = config_colour)
                plt2.set_ylabel('[CH$_4$] (ppb)')
                plt2.yaxis.set_major_formatter(y_formatter)
                plt2.xaxis.set_major_formatter(x_formatter)            
                plt2.xaxis.set_major_locator(x_tickno_formatter)
        
        plt.figtext(0.2, 0.05, 'H2O trap bypassed, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'purple')
        plt.figtext(0.2, 0.1, 'H2O trap bypassed, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'blue')
        plt.figtext(0.6, 0.05, 'H2O trap, Nafion bypassed', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'green')
        plt.figtext(0.6, 0.1, 'H2O trap, Nafion', verticalalignment='bottom', horizontalalignment='left', fontsize=8, color = 'red')
        
        plt.savefig(outputdir+os.path.splitext(data.data.filename)[0]+filesuffix+'.png', dpi=200)
        
           
        
        
        plt.show()