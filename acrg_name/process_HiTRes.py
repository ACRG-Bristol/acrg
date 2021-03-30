#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 11:09:53 2018

@author: ew14860
"""

"""
Functions for processing HiTRes footprint files
"""

from past.utils import old_div
import acrg_name.process as proc
import datetime as dt
from dateutil import parser
import gzip
import sys
import os
import os.path
import glob
import acrg_name.name as nm
import numpy as np
import xarray as xray
import getpass
import pandas as pd
from os.path import join

from acrg_config.paths import paths
data_path = paths.data

if data_path is None:
    data_path = "/work/chxmr/shared/"
    print("Default Data directory is assumed to be /work/chxmr/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

normal_fp_dir = join(data_path, 'NAME/fp/')


def update_release_time(fname):
    """
    Updates the time of release in the HiTRes output files to be start and end 
    of the release at that timestep rather than the start and end of the release
    for the first timestep of that day.
    
    fname: path to file to be updated
    """

    filename = os.path.split(fname)[1]
    filelines = proc.extract_file_lines(fname)
    
    filestartline = filelines[4]
    fileendline = filelines[5]
    endlinesplit = fileendline.split(' ')
    startlinesplit = filestartline.split(' ')
            
    fileS = startlinesplit[-1] + ' ' + startlinesplit[-2]
    fileE = endlinesplit[-1] + ' ' + endlinesplit[-2]
    filestart = parser.parse(fileS, dayfirst = True)
    fileend = parser.parse(fileE, dayfirst = True)
    release_hr = filelines[19].split(',')[-2].split('_')[-1]
    dt_hour = filelines[20].split(',')[-2].split(' ')[2]
    
    if filestart == fileend + dt.timedelta(hours = int(dt_hour)):
        print("File dates already changed for file: %s" %filename)
    else:
        print("Modifying dates in file: %s" %filename)
        start_release = fileend + dt.timedelta(hours = int(release_hr)+int(dt_hour))
        end_release = fileend + dt.timedelta(hours = int(release_hr))
    
        startlinesplit[-2] = start_release.strftime('%H%M%Z')
        startlinesplit[-1] = start_release.strftime('%d/%m/%Y')
        endlinesplit[-2] = end_release.strftime('%H%M%Z')
        endlinesplit[-1] = end_release.strftime('%d/%m/%Y')
                
        filelines[4] = ' '.join(startlinesplit)
        filelines[5] = ' '.join(endlinesplit)
        filetext = '\n'.join(filelines)
                    
        outfilename = fname
        output = gzip.open(outfilename, 'wb')
        if (sys.version_info < (3,0)):
            output.write(filetext)
        else:
            filetext_b = str.encode(filetext)
            output.write(filetext_b)
        output.close()


def process_HiTRes(domain, site, height, year, month, user_max_hour_back,
                   base_dir = "/dagage2/agage/metoffice/NAME_output/",
                   HiTRes_fields_folder = "24HrBk_Fields_files",
                   met_folder = "Met", force_met_empty = True, normal_fp_dir = normal_fp_dir):
    
    """
    domain: string - 'EUROPE'
    site: string - 'MHD'
    height: string in magl - '10' 
    year: string - '2014'
    month: string - '06'
    user_max_hour_back: int - 24
    base_dir: string, directory where [domain]_[site]_[height]magl folder is containing HiTRes footprints
    HiTRes_fields_folder: string, name of HiTRes fields folder - '24HrBk_Fields_files'
    met_folder: string, name of met folder - 'Met'
    force_met_empty: Boolean. Default is True. At this stage met data is not needed for HiTRes footprint processing.
    Met data available in the month integrated footprint files.
    normal_fp_dir: string - where the month integrated footprints are located
    
    """

    subfolder = base_dir + domain + "_" + site + "_" + height + "magl"+ "/"
    
    HiTRes_search_string = subfolder+HiTRes_fields_folder+"/*_"+year+month+"*"
    fnames = glob.glob(HiTRes_search_string)
    fnames.sort()

    if len(fnames) == 0:
        print("Can't find high time resolution footprint files " + HiTRes_search_string)
        return None
    
    if force_met_empty == False:
        met_search_str = subfolder + met_folder + "/*.txt*"
        met_files = sorted(glob.glob(met_search_str))
        if len(met_files) == 0:
            print("Can't find MET files: " + met_search_str)
            return None
        else:
            met = proc.read_met(met_files)
    if force_met_empty == True:
        met = None
    
    print("Getting normal integrated footprints from: %s" %(normal_fp_dir))
    
    if month == '12':
        totfp = nm.footprints(site, start="%s-%s-01" %(year,month), end = "%s-%s-01" %(int(year)+1, '01'), domain=domain, height = height+'m', fp_directory= normal_fp_dir)
    else:
        totfp = nm.footprints(site, start="%s-%s-01" %(year,month), end = "%s-%02d-01" %(year, int(month)+1), domain=domain, height = height+'m', fp_directory = normal_fp_dir)

    totfp = totfp.fp.to_dataset()

    user_max_hour_back = int(user_max_hour_back)

    for fi, f in enumerate(fnames):
        #update_release_time(f)

        filename = os.path.split(f)[1]
        splitfile = filename.split('_')
        release_hr = splitfile[-3]
        hr_back = int(splitfile[-1].split('.')[0])
        
        print(hr_back)
        print(user_max_hour_back)
        
        if hr_back <= user_max_hour_back:
            
            fp0 = proc.footprint_array(f,met=met)
            if fi == 0:
                time0 = fp0.time[0]
        
            fp_droplev = fp0.fp.mean(dim='lev')
            fp0 = fp0.update({'fp':fp_droplev})
            fp_addH = np.expand_dims(fp0.fp.values, 3)
            fp_addH_var = xray.Variable(['time','lat','lon', 'H_back'], fp_addH)

            FDS= xray.Dataset({'fp' : fp_addH_var},
                              coords = {'time': fp0.time,
                                        'lat': totfp.lat,
                                        'lon': totfp.lon,
                                        'H_back': [hr_back]})
            FDS = FDS.transpose('lat','lon','time','H_back')
    
            if int(release_hr) == 0 and fp0.time.values == time0.values:
                if FDS.H_back == 0:
                    FDS0 = FDS
                else:
                    FDS0 = xray.concat([FDS0, FDS], dim='H_back')
            else:
                if FDS.H_back == 0:
                    FDS1 = FDS
                else:
                    FDS1 = xray.concat([FDS1, FDS], dim='H_back')
            
                if FDS.H_back == np.max(FDS0.H_back):
                    FDS0 = xray.concat([FDS0, FDS1], dim='time')

    file_max_hour_back = np.max(FDS.H_back.values)
    
    if file_max_hour_back != user_max_hour_back:
        print("WARNING: The maximum number of hours back specified by the user (%s)\
        is greater than the maximum number of hours back available in the high\
        time resolution footprint file (%s). Creating a footprint file with the\
        available number of hours back (%s)." %(user_max_hour_back, np.max(FDS.H_back.values), np.max(FDS.H_back.values)))
    
    FDS = FDS0.copy()
    addfp = FDS.sum(dim='H_back')
    
    #make sure totfp is on same time resolution as HiTRes fp
    #can't get mean using reindex_like
    #need to use resample
    tot_fp_time_period = int(float(totfp.time.period.split(' ')[0]))
    HiTRes_fp_time_period = int((FDS.time[1].values-FDS.time[0].values).astype('timedelta64[h]')/np.timedelta64(1,'h'))
    
    if HiTRes_fp_time_period > tot_fp_time_period:
        print("Time resolution of total footprints (%d hourly) is higher than\
        time resolution in 'time' dimension of HiTRes footprints (%d hourly).\
        Averaging total footprints to match time resolution of HiTRes footprints." %(tot_fp_time_period, HiTRes_fp_time_period))
        totfp = totfp.resample(time='%dH' %HiTRes_fp_time_period).mean()
        totfp = totfp.transpose('lat','lon','time')
    
    elif HiTRes_fp_time_period < tot_fp_time_period:
        print("Time resolution of total footprints (%d hourly) is too low for\
        time resolution in 'time' dimension of HiTRes footprints (%d hourly).\
        Aborting." %(tot_fp_time_period, HiTRes_fp_time_period))
        exit
        
    remfp = totfp.fp - addfp
    remfpval = np.expand_dims(remfp.fp.values, 3)
    remfpvar = xray.Variable(['lat', 'lon','time','H_back'], remfpval)
    remfp = remfp.update({'fp' : remfpvar,
                      'H_back' : [np.max(FDS.H_back.values)+2]})
    FDS = xray.concat([FDS, remfp], dim = 'H_back')
    FDS = FDS.rename({'fp':'fp_HiTRes'})
            
    FDSattrs = {"title":" NAME footprints, %s, particles recorded two hourly for first %d hours" %(site, file_max_hour_back),
                      "author" : getpass.getuser(),
                        "date_created" : np.str(dt.datetime.today()),
                        "fp_units" : "mol/mol/mol/m2/s",
                        "comments" : "Timestamp is time of end of release, first hour back (0) is \
                        the record of particles in the boundary layer at the end of the release." }

    FDS.attrs = FDSattrs

    #MAKE NETCDF FILE!!
    
    out_folder = subfolder + "Processed_%sHrBk_Fields_files/" %(file_max_hour_back)
    
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    FDS.fp_HiTRes.encoding = {'zlib':True}
    FDS.to_netcdf(path = out_folder+'%s-%smagl_%s_%s%s.nc' %(site, height, domain, year, month), mode = 'w')
    print("%s-%smagl_%s_%s%s.nc Done" %(site, height, domain, year, month))


def process_HiTRes_all(domain, site_height_dict, start_date, end_date, user_max_hour_back,
                       base_dir = "/dagage2/agage/metoffice/NAME_output/",
                       HiTRes_fields_folder = "24HrBk_Fields_files",
                       met_folder = "Met", force_met_empty = True, normal_fp_dir = None):
    
    dates = pd.date_range(start=start_date, end = end_date, freq = 'MS')

    for i in list(site_height_dict.keys()):
        for d in dates:
            print("Getting HiTRes footprints for %02d %04d at %s" %(d.year, d.month, i))
            process_HiTRes(domain, i, site_height_dict[i], "%04d" %d.year, "%02d" %d.month, user_max_hour_back,
                           base_dir = base_dir, HiTRes_fields_folder = HiTRes_fields_folder,
                           met_folder = met_folder, force_met_empty = force_met_empty, normal_fp_dir=normal_fp_dir)
            
