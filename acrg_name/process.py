# -*- coding: utf-8 -*-
"""Process script
Created on Thu May  7 10:34:25 2015

To use this script, you will need to have a directory, labeled as:
DOMAIN_SITE_HEIGHT (e.g. EUROPE_MHD_10magl), with the following structure:

Met/
Fields_files/
Particle_files/
Processed_Fields_files/
Observations/ (satellite/column footprints only)

Each file in Fields_files/ and Particle_files/ must end with a date string
(format YYYYMMDD) which must be present in both folders 
e.g. fields_file_YYYYMMDD.txt.gz.

To process all files in a folder run:

process_all("DOMAIN", "SITE")

Output netCDF files are created one file per month, per domain, per site
in the Processed_Fields_files directory.

@author: chxmr
"""
from __future__ import print_function
from __future__ import division

from builtins import zip
from builtins import str
from builtins import range
#from past.utils import old_div
from netCDF4 import Dataset
import netCDF4
import glob
import gzip
import datetime
import re
import datetime as dt
import numpy as np
import pandas as pd
import scipy.constants as const
from acrg_grid import areagrid
from acrg_time.convert import time2sec, sec2time
from acrg_config.version import code_version
import acrg_time.convert
import os
import json
from os.path import split, realpath, exists
import xarray as xray
from scipy.interpolate import interp1d
import copy 
import matplotlib.pyplot as plt
import getpass
import traceback
import sys
import scipy
import pdb
from multiprocessing import Pool
import acrg_obs as obs

#Default NAME output file version
#This is changed depending on presence of "Fields:" line in files
namever=3

#Time formats
UTC_format = '%H%M%Z %d/%m/%Y'
NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'

# Define default met parameters. The keys are the column headings
#  in the Pandas dataframe that is created by read_met,
#  the values are the search strings that are used in read_met
#  to read the column headings
met_default = {"time": "             T",
               "press": "Pressure (Pa)",
               "temp": "Temperature (C)",
               "PBLH": "Boundary layer depth",
               "wind": "Wind speed",
               "wind_direction": "Wind direction",
               "release_lon": "xxxxxxx",
               "release_lat": "yyyyyyy"}

timestep_for_output = 0.

# Default to home directory, but update if proper directory is specified
directory_status_log = os.getenv("HOME")

from acrg_config.paths import paths
acrg_path = paths.acrg
data_path = paths.data
lpdm_path = paths.lpdm

def unzip(filename, out_filename=None, return_filename=False, delete_zipped_file=False, verbose=True):
    '''
    Unzip .gz files
    
    Args:
        filename (str)
            name of file to be unzipped, including path to file
        out_filename (str)
            filename for unzipped file
            if None, the original filename is used
        delete_zipped_file (bool)
            if True the zipped file will be deleted and be replaced with the unzipped version
            if False both the original and gzipped file will remain
    '''
    import gzip
    
    if not os.path.exists(filename):
        print(f'File does not exist : {filename}')
        return None
    
    out_filename = filename.split('.gz')[0] if out_filename is None else out_filename
    
    with gzip.GzipFile(filename, 'rb') as in_file, open(out_filename, 'wb') as out_file:
        in_data = in_file.read()
        out_file.write(in_data)
    
    if delete_zipped_file:
        print(f'Deleting file : {filename}')
        os.remove(filename)
    
    if return_filename:
        return out_filename

def load_NAME(file_lines, namever):
    """
    Loads the Met Office's NAME III grid output files returning
    headers, column definitions and data arrays as 3 separate lists.
    For a technical specification see 
    http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf
    
    Args:
        file_lines (list): 
            list of file text lines (extracted using the extract_file_lines routine)
        namever (int): 
            the NAME version. Options are 2 or 3
        
    Returns:
        headers (list): 
            Headers of NAME output files
        column_headings (list): 
            Colums headings of NAME output files
        data_arrays (list): 
            Data of NAME output files
    
    Example: 
        headers, column_headings, data_arrays = load_NAME(file_lines, namever)
    
    Note:
        Where I've modified the original Met Office routine to accept file_lines, 
        rather than the file itself, I've put the initials MLR
    """
    # loading a file gives a generator of lines which can be progressed using the next() method. 
    # This will come in handy as we wish to progress through the file line by line.
#    file_handle = file(filename)
    
    # define a dictionary which can hold the header metadata about this file
    headers = {}
    
    # skip the NAME header of the file which looks something like 'NAME III (version X.X.X)'
#    file_handle.next()
    file_lines=file_lines[1:]
    
    # read the next 16 lines of header information, putting the form "header name:    header value" into a dictionary
    if namever == 2:
        nlines = 16
    elif namever == 3:
        nlines = 17
    for _ in range(nlines):
#        header_name, header_value = file_handle.next().split(':',1)
        #MLR        
        header_name, header_value = file_lines[0].split(':',1)
        file_lines=file_lines[1:]
        
        # strip off any spurious space characters in the header name and value
        header_name = header_name.strip()
        header_value = header_value.strip()

        # cast some headers into floats or integers if they match a given header name
        if header_name in ['X grid origin', 'Y grid origin', 'X grid resolution', 'Y grid resolution']:
            header_value = float(header_value)
        elif header_name in ['X grid size', 'Y grid size', 'Number of fields',
                             'Number of preliminary cols','Number of field cols']:
            header_value = int(header_value)
        # Commented out date formatting as Start and End of release can be infinity
        #elif header_name in ['Run time', 'Start of release', 'End of release']:
            # convert the time (currently only works for name 2 format) to python datetimes
            #if namever==2:
               #header_value = datetime.datetime.strptime(header_value, UTC_format)

        headers[header_name] = header_value

    # skip the next blank line in the file.    
#    file_handle.next()    
    #MLR
    file_lines=file_lines[1:]

    # Read the lines of column definitions
    if namever == 2:
        column_headings = {}
        for column_header_name in ['species_category', 'species', 'cell_measure', 'quantity', 'unit', 'z_level', 'time']:
#            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][:-1]
            #MLR
            column_headings[column_header_name] = [col.strip() for col in file_lines[0].split(',')][:-1]
            file_lines=file_lines[1:]
            
    elif namever == 3:
        # skip the next line (contains the word Fields:) in the file.    
#        file_handle.next()
        #MLR
        file_lines=file_lines[1:]
        column_headings = {}
        for column_header_name in ['category', 'name', 'quantity', 'species', 'unit','source','ensemble','time_averaging',
                                   'horiz_averaging','vert_averaging','prob_percentile','prob_percentile_ensemble',
                                   'prob_percentile_time','time', 'z_level', 'D']:
#            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][:-1]  
            column_headings[column_header_name] = [col.strip() for col in file_lines[0].split(',')][:-1]
            file_lines=file_lines[1:]

    # convert the time to python datetimes
    new_time_column_header = []
    for i, t in enumerate(column_headings['time']):
        # the first 4 columns aren't time at all, so don't convert them to datetimes
        if i >= 4:
            if namever == 2:
                new_time_column_header.append(datetime.datetime.strptime(t, UTC_format))
            elif namever == 3:
                new_time_column_header.append(datetime.datetime.strptime(t, NAMEIII_short_format))
        else:
            new_time_column_header.append(t)
        
    column_headings['time'] = new_time_column_header
    
    #MLR
    #This cuts off the first line of output so removing
    file_lines=file_lines[1:]
    
    # make a list of data arrays to hold the data for each column 
    data_shape = (headers['Y grid size'], headers['X grid size'])
    if namever==2:
        data_arrays = [np.zeros(data_shape, dtype=np.float32) for i in range(headers['Number of fields'])]
    elif namever==3:
        data_arrays = [np.zeros(data_shape, dtype=np.float32) for i in range(headers['Number of field cols'])]
      
    # iterate over the remaining lines which represent the data in a column form
#    for line in file_handle:
    #MLR
    for line in file_lines:
        
        # split the line by comma, removing the last empty column caused by the trailing comma
        vals = line.split(',')[:-1]
        
        # cast the x and y grid positions to floats and convert them to zero based indices
        # (the numbers are 1 based grid positions where 0.5 represents half a grid point.)
        if namever == 2:
            x = int(float(vals[0]) - 1.5)
            y = int(float(vals[1]) - 1.5)
        elif namever == 3:
            x = int(float(vals[0]) - 1)
            y = int(float(vals[1]) - 1)
        
        # populate the data arrays (i.e. all columns but the leading 4) 
        for i, data_array in enumerate(data_arrays):
            data_array[y, x] = float(vals[int(i) + 4])
    
    if 'cell_measure' in column_headings:
        #Extract the time integration period
        dt_fp = int(column_headings['cell_measure'][4][0:3])
        #This should only apply to NAME version 2 HiTRes footprints which currently do not have the correct start and end release times in the header
        if dt_fp < 24:
            end_release_date = headers['End of release'][-10:]
            end_release_tz = headers['End of release'][4:7]
            end_release_time = column_headings['species'][4][-2:]
            end_release = dt.datetime.strptime(end_release_time + "00" + end_release_tz + " " + end_release_date, '%H%M%Z %d/%m/%Y')
            start_release = end_release + dt.timedelta(hours = dt_fp)
            endreleaseline = end_release.strftime('%H%M%Z') + end_release_tz + " " + end_release.strftime('%d/%m/%Y')
            startreleaseline = start_release.strftime('%H%M%Z') + end_release_tz + " " + start_release.strftime('%d/%m/%Y')
            headers['End of release'] = endreleaseline
            headers['Start of release'] = startreleaseline
    
    return headers, column_headings, data_arrays
    

def extract_file_lines(fname):
    '''
    Determine whether file is zipped or not, and then extract text into
    file_lines variable.
    
    Args:
        fname (str): 
            file name to extract text from
        
    Returns:
        file_lines (str array): 
            lines of text extracted from fname 
        
    Example:
        fname = 'sometextfile.txt'
        file_lines = extract_file_lines(fname)
    '''
    
    if fname[-3:].upper() == ".GZ":
        #Unzip file
        f=gzip.open(fname, 'rt')
        file_text=f.read()
        f.close()
        file_lines=file_text.splitlines()
    else:
        file_lines=[i.strip() for i in open(fname).readlines()]

    return file_lines


def read_file(fname):
    '''
    Read a given fields file and break into header, column_headings and data_arrays.
     
    Args:
        fname (str): 
            File name to extract data from.
        
    Returns:
        headers (list): 
            Headers of NAME output files
        column_headings (list): 
            Colums headings of NAME output files
        data_arrays (list): 
            Data of NAME output files
    '''
    
    #Extract line-by-line file contents from either txt.gz or txt file    
    file_lines=extract_file_lines(fname)

    #Determine NAME version
    if file_lines[19].strip().upper() == 'FIELDS:':
        namever=3
    else:
        namever=2
    
    #Load header, column headings and data using Andrew Jones script
    header, column_headings, data_arrays = \
        load_NAME(file_lines, namever)

    return header, column_headings, data_arrays, namever


def define_grid(namever, header, column_headings, satellite = False, upper_level = None):
    '''
    Define output grid using file header information.
    
    Note: header, column_headings can be extracted using read_file() function and passed to this function.
    
    Args:
        header (list): 
            header information
        column_headings (list): 
            headings for columns
        satellite (bool, optional): 
            If using satellites then True. Default = False.
            For satellite instructions, see help(process.satellite_vertical_profile).
        upper_level (int, optional):
            Only needed when satellite=True.
            Highest level number from within the NAME run for the satellite data.
            Default = None.
    
    Returns:
        lons (np.array): 
            NAME grid longitudes
        lats (np.array): 
            NAME grid latitudes
        levs (list): 
            Levels of NAME grid
        time (list): 
            Time step label at start of each time period
        timeStep (float): 
            Timestep in hours
        
    Example:
        lons, lats, levs, time, timeStep = define_grid(header, column_headings, satellite = False)
        lons, lats, levs, time, timeStep = define_grid(header, column_headings, satellite = True, upper_level = 17)
        
    Note:
        Largely taken from Met Office Python code.
    '''

    #Get file info
    z_level=column_headings['z_level'][4:]
    time=column_headings['time'][4:]
    
    if namever==2:
        species_name=column_headings['species'][4:]
    else:
        species_name=column_headings['name'][4:]

    if namever == 2:
        timeFormat=UTC_format
    else:
        timeFormat=NAMEIII_short_format


    #Longitude and latitude    
    #Co-ordinates OF CENTRE OF CELL
    Xcount = header['X grid size']
    Xstep = header['X grid resolution']
    Xstart = header['X grid origin']
    
    Ycount = header['Y grid size']
    Ystep = header['Y grid resolution']
    Ystart = header['Y grid origin']
    if namever == 2:
        lons=np.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep) + Xstep/2.
        lats=np.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep) + Ystep/2.
    elif namever == 3:
        lons=np.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep)
        lats=np.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep)
    else:
        return None
    
    status_log("Bottom-left grid point (CENTRE): " + str(lons[0]) + "E, " + \
               str(lats[0]) + "N",
               print_to_screen=False)

    if satellite == False:
        
        #Get levels and sort in increasing height
        levs=list(set(z_level))
        lev_sort=[]
        for lev in levs:
            digits = re.findall(r"[-+]?\d*\.\d+|\d+", lev)
            if not digits:
                lev_sort.append(10000000.) #PUT STUFF LIKE BOUNDARY LAYER TO END
            else:
                lev_sort.append(float(digits[0]))
    
        levs=[lev for (sort, lev) in sorted(zip(lev_sort, levs))]
        nlevs=len(levs)
    
        #Time
        #ESTIMATE NUMBER OF TIMESTEPS PER FILE: ASSUMES ONLY ONE SPECIES PER FILE
        ntime=len(time)//nlevs
        
    else:
        if upper_level is None:
            status_log("Upper Level must be specified for satellite data.",print_to_screen=True,
                       error_or_warning="error")
        # Calculate number of time steps in file
        nlevs = upper_level
        levs = list(range(nlevs))
        ntime = len(time)//nlevs

    #Time steps    
    timeRef=datetime.datetime.strptime(header['End of release'], timeFormat)
    timeEnd=datetime.datetime.strptime(header['Start of release'], timeFormat)
    
    # Timestep in hours
    timeStep=(timeEnd - timeRef).total_seconds()/3600./ntime
    
    # Labelling time steps at START of each time period!
    time = [timeRef + datetime.timedelta(hours=timeStep*i) \
            for i in range(ntime)]
    
    status_log("Timestep: %d minutes" % round(timeStep*60),
               print_to_screen=False)
    status_log("Levels: %d " % nlevs,
               print_to_screen=False)
    status_log("NAME version: %d" % namever,
               print_to_screen=False)
        
    return lons, lats, levs, time, timeStep


def met_empty():
    '''
    Create an empty Met dictionary. Not recommended.
    '''
    
    met = pd.DataFrame({key: 0. for key in list(met_default.keys()) if key != "time"},
                          index = [dt.datetime(1900, 1, 1), dt.datetime(2100, 1, 1)])
    met.index.name = "time"
    met["press"] = [100000., 100000.] #Pa
    met["temp"] = [10., 10.]    #C
    
    status_log("NO MET", error_or_warning="warning")

    return met


def read_met(fnames, met_def_dict=None,vertical_profile=False,satellite=False):
    '''
    Given list of filenames, extract site meteorology and concatenate into Pandas dataframe.
    
    Args:
        fnames (list): 
            List of filenames for met data.
        met_def_dict (dict, optional): 
            Dictionary of met parameters.
            Default=None, which takes default entries.
        vertical_profile (bool, optional): 
            TODO. Default=False
        satellite (bool, optional):
            If using satellites then True. Default = False.
    
    Returns:
        output_df (pandas.Dataframe): 
            Dataframe of met data. 
        
    Example:
        df = read_met(fnames, met_def_dict=None, vertical_profile=False)
        
    Note:
        A dictionary of output column names and met file header search strings
        is given in the met_default dictionary at the top of process.py.
        Vertical profile assumes a file structure that consists of a time column 
        14 temperature (C) columns, followed by 14 pressure (Pa), and 14 potential 
        temperature (K)
    '''


    if isinstance(fnames, list):
        fnames.sort()
    else:
        fnames=[fnames]

    if met_def_dict is not None:
        # overwrite met_default with custom list for prcessing vertical profile met
        met_default2 = met_def_dict
    else:
        met_default2 = met_default
        #column_indices = {key: -1 for key, value in met_def_dict.iteritems()}
    #else:
    column_indices = {key: -1 for key, value in met_default2.items()}
    
    output_df = []
        
    for fname in fnames:
        
        #This should now check for column headings
        met = extract_file_lines(fname)
        if fname[-3:].upper() == '.GZ':
            compression="gzip"
        else:
            compression=None

        #Get rid of top of file
        for i in range(len(met)):
            cols = [col.strip() for col in met[i].split(',')]
            if len(cols)>1:
                a = i
                break

        m = pd.read_csv(fname, skiprows = a, encoding='utf-8',
                            compression=compression)
        m = m.fillna('')
        
        m = np.array(m)

        Xcol = None
        Ycol = None
        X_file = None
        Y_file = None
        
        #Find columns with header names
        for row in m:
            for coli, col in enumerate(row):
                if isinstance(col, str):
                    
                    #Work out column indices by searching through met_default
                    for key in list(met_default2.keys()):
                        if met_default2[key] in col:
                            column_indices[key] = coli

                    # Check whether there is an X and Y column
                    if "X (Lat-Long)" in col:
                        Xcol = coli
                    if "Y (Lat-Long)" in col:
                        Ycol = coli

                    # Check whether X and Y is in col header
                    if "X = " in col:
                        X_file = float(col.strip().split(" ")[2])
                    if "Y = " in col:
                        Y_file = float(col.strip().split(" ")[2])
                    
        #Find where data starts
        for i, mi in enumerate(m[:, 0]):
            if str(mi).strip() != '':
                break
        
        m2 = m[i+1:, :]
        
        #Create arrays/lists to store release locations
        if Xcol is not None:
            X=m2[:, Xcol].astype(float)
            Y=m2[:, Ycol].astype(float)
        elif X_file is not None:
            X = [X_file for i in m2[:, column_indices["time"]]]
            Y = [Y_file for i in m2[:, column_indices["time"]]]
        
        #Construct dictionary
        met_dict = {}
        if vertical_profile == True:
            vp_met_cols = {"time": 0,"temp20": 1,
               "temp40": 2,"temp60": 3,
               "temp80": 4,"temp100": 5,
               "temp120": 6,"temp140": 7,
               "temp160": 8,"temp180": 9,
               "temp200": 10,"temp220": 11,
               "temp240": 12,"temp260": 13,
               "temp280":14,"temp300": 15,
               "press20":16,"press40": 17,
               "press60":18,"press80": 19,
               "press100": 20,"press120": 21,
               "press140": 22,"press160": 23,
               "press180": 24,"press200": 25,
               "press220": 26,"press240": 27,
               "press260": 28,"press280": 29,
               "press300": 30,
               "theta20": 31,"theta40": 32,
               "theta60": 33,"theta80": 34,
               "theta100": 35,"theta120": 36,
               "theta140": 37,"theta160": 38,
               "theta180": 39,"theta200": 40,
               "theta220": 41,"theta240": 42,
               "theta260": 43,"theta280": 44,
               "theta300": 45}
 
            for key in list(met_default2.keys()):
                if key != "time":
                    met_dict[key] = m2[:, vp_met_cols[key]].astype(float)
            met_dict["release_lon"] = X
            met_dict["release_lat"] = Y
            
        else:
            for key in list(met_default2.keys()):
                if column_indices[key] != -1 and key != "time":
                    met_dict[key] = m2[:, column_indices[key]].astype(float)
            met_dict["release_lon"] = X
            met_dict["release_lat"] = Y

        #Construct dataframe
        # calling to_datetime on a series is MUCH faster than the previous list comprehension
        times = pd.Series([d.strip() for d in m2[:,column_indices["time"]]])
        output_df_file = pd.DataFrame(met_dict,
                        index=pd.to_datetime(times, format='%d/%m/%Y %H:%M %Z'))
        
        output_df.append(output_df_file)
    
    
    # Concatenate list of data frames
    output_df = pd.concat(output_df)
    if vertical_profile == True:
        output_df = output_df.sort_index()

    
    if satellite:
        try:
            re_met_match = re.compile(r"Met_[-]?\d+[.]\d+_[-]?\d+[.]\d+_\d{2,5}")
            file_match = [re.search(re_met_match,f).group() for f in fnames]
        except AttributeError:
            status_log("Could not match Met data filename to expected format.",error_or_warning="warning")
            output_df["label"] = np.empty((len(output_df.index)),dtype=str)
        else:
            label = [match.split('_')[-1] for match in file_match]
            output_df["label"] = label
        
    if satellite:
        # Sort by label axis which includes point number and level
        output_df = output_df.sort_values(by=["label"])
    else:
        # Sort the dataframe by time
        output_df = output_df.sort_index()
    
    output_df.index.name = "time"
    
    return output_df


def met_satellite_split(met):
    '''
    Splits the met input into points for satellite data if the data has been grouped into all levels and points
    for one day.
    
    Looks for an additional column in the met dataframe called "label" which contains the label
    associated with the input point. For data which has been grouped by day the label values should in
    the format "PPPLL" e.g. '00507' which translates to Point 005, Level 07.
    If this format or column is not found then the met input will be returned unchanged.
        
    Args:
        met (pandas.DataFrame/list) : 
            Output from read_met function. Expects a Dataframe but can also handle a list 
            containing one dataframe.
    
    returns:
        list(pandas.DataFrame):
            Met data split out into separate points each one containing all the associated levels.
    '''
    min_label_digit = 5 # Expect format PPPLL
    lev_digit = 2 # Expect last 2 digits of the label to be the level
    
    if isinstance(met,list):
        if len(met) == 1:
            met = met[0]
        else:
            print("Only expect pandas.Dataframe or 1 item list as input to met_satellite_split. Len: {}, input included: {}".format(len(met),met))
            return None
    
    if "label" in met.columns.values:
        if len(met["label"][0]) >= min_label_digit:
            point_num = np.unique([label[:-lev_digit] for label in met["label"].values])
            met = [met[met["label"].str.startswith(num)] for num in point_num]
        else:
            met = [met]
    else:
        met = [met]
    
    return met
 
def particle_locations(particle_file, time, lats, lons, levs, heights, id_is_lev = False,
                       satellite = False,upper_level=17):
    '''
    Read and process a particle location file.
    
    Note: lons, lats, levs, time can be derived from corresponding fields files using define_grid() function
    and passed to this function.
    
    Args:
        particle_file (str): 
            File of file containing particle locations
        time (list): 
            List of timestamps for NAME run
        lats (np.array): 
            Array of latitudes of domain for NAME run
        lons (np.array): 
            Array of longitudes of domain for NAME run
        levs (list): 
            List of levels within input data for NAME run
        heights (np.array/list): 
            TODO: Needs better explanation here.
            Particle location heights
        id_is_lev (bool, optional):
            Indicate that the input "Id" value from the particle files relates to the level rather than the 
            timestamp.
            Default = False.
        satellite (bool, optional):
            Used for satellite data. This flag is used to know how to process the "Id" value to link to time values
            and levels. Expect all time points have the same number of levels corresponding to upper_level value.
            Default = False.
        upper_level (int):
            Needed if satllite=True. Used to separate out the input "Id" values to link to both time and level 
            values.
        
    Returns:
        hist (xarray.Dataset):
            Number of particles at each boundary.
            Contains "pl_n","pl_e","pl_s""pl_w" data variables.
            Dimensions are "time", "lev", "lon", "height"
    '''

    def particle_location_edges(xvalues, yvalues, x, y):
        '''
        The particle_location_edges function calculates the histogram along with the x and y edge values
        for input xvalues, yvalues.
        
        
        '''
        dx = x[1] - x[0]
        xedges = np.append(x - dx/2., x[-1] + dx/2.)
        dy = y[1] - y[0]
        yedges = np.append(y - dy/2., y[-1] + dy/2.)
        
        hist, xe, ye = np.histogram2d(xvalues, yvalues,
                                         bins = (xedges, yedges))
        
        return hist
    
    def mean_age_edges(xvalues, yvalues, agevalues, x, y):
        '''
        The mean_age_edges function calculates the mean age of particles along with the x and y edge values
        for input xvalues, yvalues.
        
        
        '''
        dx = x[1] - x[0]
        xedges = np.append(x - dx/2., x[-1] + dx/2.)
        dy = y[1] - y[0]
        yedges = np.append(y - dy/2., y[-1] + dy/2.)
        
        # supress the warning on this divide as it will intentionally divide by 0 where no particles
        with np.errstate(divide='ignore', invalid='ignore'):
            meanvals = np.histogram2d(xvalues, yvalues, bins = (xedges, yedges), weights = agevalues)[0] \
                            / np.histogram2d(xvalues, yvalues, bins = (xedges, yedges))[0]

        meanvals[np.isnan(meanvals)]=0
        return meanvals
    
    edge_lons = [lons.min(), lons.max()]
    edge_lats = [lats.min(), lats.max()]
    dlons = lons[1] - lons[0]
    dlats = lats[1] - lats[0]
 
    # determine whether the domain is unbounded in the x direction and use this to set particle location histograms in 
    # east and west to 0
    if np.isclose((edge_lons[1] - edge_lons[0] + dlons),360,atol=1e-6,rtol=0):
        xUnbounded = True
    else:
        xUnbounded = False    
 

    hist = []
    particles = []
    
    hist = xray.Dataset({"pl_n":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "pl_e":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights)))),
                         "pl_s":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "pl_w":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights)))),
                         "mean_age_n":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "mean_age_e":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights)))),
                         "mean_age_s":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "mean_age_w":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights))))},
                        coords={"lat": lats, "lon":lons, "lev":levs, "height":heights, "time":time})
    
    #Variables to check domain extents
    particle_extremes = {"N": -90., "E": -360. ,"S": 90.,"W": 360.}
    
    status_log("Particle locations " + particle_file, print_to_screen=False)
    
    if particle_file[-3:].upper() == '.GZ':
        compression="gzip"
    else:
        compression=None
    
    df = pd.read_csv(particle_file, compression=compression, sep=r"\s+")

    particles_record = []
    
    id_values = set(np.array(df["Id"]))
    
    if satellite:
        npoints = len(id_values)
        if np.any(np.isnan(list(id_values))):
            status_log("Unable to read all id values from particle file. Check final column in particle_location* file is separated by 1 or more spaces.",
                       error_or_warning="error",print_to_screen=True)     
        if (npoints)%upper_level != 0:
            status_log("Total number of fields {} is not exactly divisable by the number of levels {}. Same number of levels must be used for each satellite point.".format(npoints,upper_level),
                       error_or_warning="error",print_to_screen=True)
        t=0
        l=0
    
    df_n = df[(df["Lat"] > edge_lats[1] - dlats/2.)]
    df_e = df[(df["Long"] > edge_lons[1] - dlons/2.)]
    df_s = df[(df["Lat"] < edge_lats[0] + dlats/2.)]
    df_w = df[(df["Long"] < edge_lons[0] + dlons/2.)]
    
    for i in id_values:
        
        if satellite:
            # Loop through time points and levels
            if i==1: # Use initial t=0,l=0 values
                pass
            elif (i-1)%upper_level: # Loop through levels while i is not exactly divisable by upper_level
                l+=1
            else: # Set time value to next value and reset level to 0
                t+=1
                l=0
            slice_dict = {"time": [t], "lev": [l]}
        elif id_is_lev:
            # Only one time step allowed for id_is_lev option
            slice_dict = {"time": [0], "lev": [i-1]}
        else:
            slice_dict = {"time": [i-1]}
            
        #Northern edge
        dfe = df_n[df_n["Id"] == i]
        hist.pl_n[slice_dict] = \
            particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                    lons, heights)
        hist.mean_age_n[slice_dict] = \
            mean_age_edges(dfe["Long"].values, dfe["Ht"].values,dfe["Age(hr)"].values,
                                    lons, heights)        
        #Eastern edge
        # If domain is unbounded in x direction (i.e. 0-360), set E and W directions manually to 0
        # This is done manually because any small number of particles that die in domain at 
        # the end of the run can be associated with E and W directions
        # In addition, the particle location files have long values from -180 to 180, but if NAME run with
        # xUnbounded flag, the lon values can take on numbers between -180 to 360 causing issues with the
        # interpretation of 'df["Long"] > edge_lons[1] - dlons/2.'
        # This is only an issue if NAME was run with xUnbounded
        
        if xUnbounded:
            hist.pl_e[slice_dict] = 0
            hist.mean_age_e[slice_dict] = 0
        else:
            dfe = df_e[df_e["Id"] == i]
            hist.pl_e[slice_dict] = \
            particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                    lats, heights)
            hist.mean_age_e[slice_dict] = \
            mean_age_edges(dfe["Lat"].values, dfe["Ht"].values, dfe["Age(hr)"].values,
                                    lats, heights)

        #Southern edge
        dfe = df_s[df_s["Id"] == i]
        hist.pl_s[slice_dict] = \
            particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                    lons, heights)
        hist.mean_age_s[slice_dict] = \
            mean_age_edges(dfe["Long"].values, dfe["Ht"].values,dfe["Age(hr)"].values,
                                    lons, heights)   
        #Western edge
        if xUnbounded:
            hist.pl_w[slice_dict] = 0 
            hist.mean_age_w[slice_dict] = 0 
        else:
            dfe = df_w[df_w["Id"] == i]
            hist.pl_w[slice_dict] = \
                particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                        lats, heights)
            hist.mean_age_w[slice_dict] = \
                mean_age_edges(dfe["Lat"].values, dfe["Ht"].values, dfe["Age(hr)"].values,
                                        lats, heights)

        #Calculate total particles and normalise
        hist_sum = hist[slice_dict].sum()
        particles = int(sum([hist_sum[key].values for key in ["pl_n", "pl_e", "pl_s", "pl_w"]]))
        particles_record.append(str(particles))

        if particles > 0.:
            for key in ["pl_n", "pl_e", "pl_s", "pl_w"]:
                hist[key][slice_dict] = hist[key][slice_dict]/particles
        else:
            status_log("No particles have reached edge",
                       error_or_warning="warning")
            if i > 1:
                # Note: i is one-indexed, t,l are zero-indexed. 
                status_log("Copying lower level/previous time step",
                           error_or_warning="warning")
                if id_is_lev:
                    # Copying from previous level
                    slice_dict_prev = {"time": [0], "lev": [i-2]}
                elif satellite:
                    if l-1 < 0:
                        # Copying from lowest level from previous time step
                        slice_dict_prev = {"time": [t-1], "lev": [0]}
                    else:
                        # Copying from previous level
                        slice_dict_prev = {"time": [t], "lev": [l-1]}
                else:
                    # Copying from previous time step
                    slice_dict_prev = {"time": [i-2]}
                for key in hist.data_vars.keys():
                    hist[key][slice_dict] = hist[key][slice_dict_prev].values
        
    # Store extremes
    if df["Lat"].max() > particle_extremes["N"]:
        particle_extremes["N"] = df["Lat"].max()
    if df["Lat"].min() < particle_extremes["S"]:
        particle_extremes["S"] = df["Lat"].min()
    if df["Long"].max() > particle_extremes["E"]:
        particle_extremes["E"] = df["Long"].max()
    if df["Long"].min() < particle_extremes["W"]:
        particle_extremes["W"] = df["Long"].min()

    status_log("Number of particles reaching edge: " + ", ".join(particles_record),
               print_to_screen = False)
    
    #Check extremes
    if particle_extremes["N"] < edge_lats[1] - dlats/2.:
        status_log("CHECK DOMAIN EDGE TO NORTH", error_or_warning="warning",
                   print_to_screen=False)
    if particle_extremes["E"] < edge_lons[1] - dlons/2.:
        status_log("CHECK DOMAIN EDGE TO EAST", error_or_warning="warning",
                   print_to_screen=False)
    if particle_extremes["S"] > edge_lats[0] + dlats/2.:
        status_log("CHECK DOMAIN EDGE TO SOUTH", error_or_warning="warning",
                   print_to_screen=False)
    if particle_extremes["W"] > edge_lons[0] + dlons/2.:
        status_log("CHECK DOMAIN EDGE TO WEST", error_or_warning="warning",
                   print_to_screen=False)

    return hist


def footprint_array(fields_file, 
                    particle_file = None,
                    met = None,
                    satellite = False,
                    time_step = None,
                    upper_level = None,
                    obs_file = None,
                    use_surface_conditions = True,
                    species = None,
                    user_max_hour_back = 24.,
                    lifetime_hrs = None):
    '''
    Convert text output from given files into arrays in an xarray.Dataset.
    
    Args:
        fields_file (str): 
            Filename of footprint text file.
        particle_file (str, optional): 
            File containing particle locations.
            Default = None.
        met (dict): 
            Dictionary of Met data.
            Output from read_met() function.
            Default = None.
        satellite (bool, optional): 
            If using satellites then True. Default = False.
            For satellite instructions, see help(process.satellite_vertical_profile).
        time_step (float, optional): 
            Timestep of footprint. 
            Needed if timestep cannot be extracted from input fields file.
            Default = None 
        upper_level (int, optional):
            Only needed when satellite=True. Highest level number from within the NAME run for the satellite data.
            Default = None.
        obs_file (str, optional):
            Optional for satellite=True; not used otherwise. Files containing observation data (.nc file).
            Time values cannot always be accurately determined from the input fields file. If obs_file is
            specified, time values will be extracted from these files.
            Default = None.
        use_surface_conditions (bool, optional) :
            Use default expected surface conditions for meteorological values
            if converting from gs/m3 to mol/mol / mol/m2/s units.
            P/T ratio is fixed as 345 based on typical surface conditions.
            Default = True.
        species (str,optional):
            Defaults to None which will process footprints for an inert species
            Otherwise will look in json file for lifetime and process a species-
            specific footprint
        user_max_hour_back (float, required when species = 'CO2'):
            Defaults to 24 hours. This is the maximum amount of time back from release time that 
            hourly footprints are calculated.
        lifetime_hrs (float, optional):
            Defaults to None which will process footprints for an inert species.
            Otherwise will process it for the specified loss timescale.
            
    Returns:
        fp (xarray.Dataset): 
            Dataset of footprint data.       
    '''
    
    global timestep_for_output

    if satellite and not upper_level:
        status_log("Upper level must be specified for satellite data.",
                   error_or_warning="error")
        return None

    status_log("Reading... " + os.path.split(fields_file)[1])

    if lifetime_hrs is None:
        lifetime_attr = "No loss applied"
    else:
        lifetime_attr = str(lifetime_hrs)
        
    # note the code will fail early if a species is not defined or if it is long-lived   
    if 'MixR_hourly' in fields_file:
        fields_ds       = Dataset(fields_file, "r", format="NETCDF4")
        lons            = np.array(fields_ds.variables["Longitude"][:])
        lats            = np.array(fields_ds.variables["Latitude"][:])
        attributes      = fields_ds.ncattrs()
        releasetime_str = [s for s in attributes if 'ReleaseTime' in s]
        releasetime     = [f.split("ReleaseTime")[1] for f in releasetime_str]
        time            = [datetime.datetime.strptime(f, '%Y%m%d%H%M') for f in releasetime]
        levs            = ['From     0 -    40m agl'] # not in the file, not sure if needed, placeholder
        timeStep        = fields_ds.getncattr('ReleaseDurationHours')
        data_arrays     = []     
        
        for rtime in releasetime:
            rt_dt       = datetime.datetime.strptime(rtime, '%Y%m%d%H%M')
            fp_grid     = np.zeros((len(lats), len(lons)))
            fields_vars = fields_ds.get_variables_by_attributes(ReleaseTime=rtime)
            outputtime  = [fields_vars[ii].getncattr('OutputTime') for ii in range(len(fields_vars))]
            outputtime  = list(sorted(set(outputtime)))
            
            for ot in outputtime:
                data    = [f for f in fields_vars if f.getncattr('OutputTime') == ot]
                xindex  = [f for f in data if 'Xindex' in f.name][0][:]-1 # Alistair's files index from 1
                yindex  = [f for f in data if 'Yindex' in f.name][0][:]-1 # Alistair's files index from 1
                fp_vals = [f for f in data if 'NAMEdata' in f.name][0][:]
                fp_grid_temp     = np.zeros((len(lats), len(lons)))
                ot_dt            = datetime.datetime.strptime(ot, '%Y%m%d%H%M')
                fp_timedelta_hrs = (rt_dt - ot_dt).total_seconds()/3600 + timeStep/2 # average time elapsed in hours
                # turn this data into a grid
                for ii in range(len(xindex)):
                    fp_grid_temp[yindex[ii], xindex[ii]] = fp_vals[ii]
                
                if species == "CO2":
                    # add to the total for that release time    
                    fp_grid+=fp_grid_temp
                    
                    hr_back = datetime.datetime.strptime(rtime, '%Y%m%d%H%M') - datetime.datetime.strptime(ot, '%Y%m%d%H%M')
                    hr_back = hr_back.total_seconds()/3600.
                    
                    FDS_rt = xray.Dataset({"fp_HiTRes": (["time", "lev", "lat", "lon", "H_back"],
                                           np.zeros((len(time), len(levs),len(lats), len(lons), int(user_max_hour_back))))},
                                          coords={"time": time, "lev": levs, "lat": lats, "lon": lons,  
                                                  "H_back": np.arange(0,user_max_hour_back)})
                    if hr_back < user_max_hour_back:
                        FDS_rt.fp_HiTRes.loc[dict(lev='From     0 -    40m agl', time=rt_dt, H_back=hr_back)] = fp_grid_temp
                else:                 
                    # add to the total for that release time    
                    fp_grid+=fp_grid_temp*np.exp(-1*fp_timedelta_hrs/lifetime_hrs) # lifetime applied
                    
            data_arrays.append(fp_grid)    
            
    else:
        header, column_headings, data_arrays, namever = read_file(fields_file)
        # Define grid, including output heights    
        lons, lats, levs, time, timeStep = define_grid(namever, header, column_headings,
                                                       satellite = satellite,
                                                       upper_level = upper_level)

    if met is None:
        met = met_empty()
        force_met_empty = True
    else:
        force_met_empty = False

    if type(met) is not list:
        met = [met]

    if satellite and obs_file:
        time = extract_time_obs(obs_file)
    
    # If time_step is input, overwrite value from NAME output file
    if time_step is not None:
        timeStep = time_step
    
    timestep_for_output = timeStep

    dheights = 1000
    heights = np.arange(0, 19001, dheights) + dheights/2.
    
    # Get area of each grid cell
    area=areagrid(lats, lons)
    
    # Get particle locations
    if particle_file is not None:
        particle_hist = particle_locations(particle_file,
                                           time, lats, lons, levs,
                                           heights, satellite = satellite,
                                           upper_level = upper_level)
    else:
        status_log("No particle location file corresponding to " + fields_file,
                   error_or_warning="error")

    nlon=len(lons)
    nlat=len(lats)
    nlev=len(levs)
    ntime=len(time)
    status_log("Time steps in file: %d" % ntime, print_to_screen = False)
    
    if 'MixR_hourly' in fields_file: 
        units_str = fields_ds.getncattr('units')
    else:
        z_level=column_headings['z_level'][4:] # anita: can be deleted?
        time_column=column_headings['time'][4:] # anita: can be deleted?
        units_column=column_headings["unit"][4:] # Find column containing NAME output units (e.g. g s/m3 or ppm s)
        units_str = units_column[0]
    
    # Set up footprint dataset
    fp = xray.Dataset({"fp": (["time", "lev", "lat", "lon"],
                              np.zeros((ntime, nlev, nlat, nlon)))},
                        coords={"lat": lats, "lon":lons, "lev": levs, 
                                "time":time})

    fp.fp.attrs = {"units": "(mol/mol)/(mol/m2/s)", "loss_lifetime_hrs": lifetime_attr}

    # anita I don't think these attributes below are read in - they are re-written in write_netcdfs
    fp.lon.attrs = {"units": "degrees_east"}
    fp.lat.attrs = {"units": "degrees_north"}
    fp.time.attrs = {"long_name": "start of time period"}
    
    
    ## Reformat met data to separate into levels
    if satellite:
        met = met_satellite_split(met)
        if len(met) != len(time):
            status_log("Expect number of extracted met points to match number of time points. " + \
                       "Check metError within NAME output to check if met has been generated correctly",
                       error_or_warning="error")

    # Add in met data
    met_dict = {key : (["time", "lev"], np.zeros((len(time), len(levs))))
                for key in list(met[0].keys())}
    met_ds = xray.Dataset(met_dict,
                          coords = {"time": (["time"], time),
                                   "lev": (["lev"], levs)})

    for i,t in enumerate(time):
        
        if len(met[0].index) == 1:
            if len(levs) > 1:
                status_log("ONLY ONE MET LEVEL. " + \
                           "ASSUMING MET CONSTANT WITH HEIGHT!",
                           error_or_warning="error")
        if satellite:
            metr = met[i] # Same number of met_dataframes in list as time points
        else:
            # Re-index met dataframe to each time point
            met[0] = met[0].tz_localize(None)
            if force_met_empty == False:
                metr = met[0][~met[0].index.duplicated(keep='first')].reindex(index = np.array([t]))
            if force_met_empty == True:
                metr = met[0].tail(1)
            if np.isnan(metr.values).any():
                
                print(t)
                raise ValueError("No met data for given date %s" % t)
            
        for key in list(metr.keys()):
            if key != "time":
                met_ds[key][{'time':[i]}] = metr[key].values.reshape((1,len(levs)))
        
    if "label" in met_ds.data_vars:
        met_ds = met_ds.drop("label")
            

    # Merge met dataset into footprint dataset
    fp = fp.merge(met_ds)

    # Add in particle locations
    if particle_file is not None:
        fp = fp.merge(particle_hist)
    
    if species == 'CO2':
        fp = fp.merge(FDS_rt)
 
    # Extract footprint from columns assuming ppm s units
    def convert_units(fp, slice_dict, column, units, use_surface_conditions = True):
        '''
        Conversion is based on inputs units
        
        If units are 'ppm s':
            Convert from [ppm s] (i.e. mu-mol/mol) units to [(mol/mol) / (mol/m2/s)].
            
            Using conversion:
                sensitivity [mu-mol/mol s] * area [m2] * molar mass [g/mol] * 1e-6 [mu] 
                    / (time [s] * release rate [g/s])
        
        If units are 'g s / m^3' (or 'gs/m3'):
            Convert from [g s/m3] units to [(mol/mol) / (mol/m2/s)].
            
            Using conversion:
                sensitivity [g s/m3] * area [m2] * RT/P [m3/mol]
                / (time [s] * release rate [g/s])
            
        Note:
            release rate is assumed to be 1. [g/s]
            molecular weight in the NAME run itself was set to 1.0, so molar mass will
            be 1.0 in this calculation as well.
        
        Optional args:
            use_surface_conditions (bool, optional) :
                If the input units are 'gs/m3', use representation surface P/T ratio of 345 rather
                than using the input meteorological data.
                Default = True
        '''
        units_no_space = units.replace(' ','')
        if units == "g s / m^3" or units == "gs/m3" or units_no_space == "gs/m^3"  or units_no_space == "gs/m3":
            if use_surface_conditions:
                molm3=345./const.R ## Surface P/T ratio we would expect over Europe (345).
            else:
                molm3=fp["press"][slice_dict].values/const.R/\
                    const.convert_temperature(fp["temp"][slice_dict].values.squeeze(),"C","K")
            fp.fp[slice_dict] = data_arrays[column]*area/ \
                (3600.*timeStep*1.)/molm3
        elif units == "ppm s" or units_no_space == "ppms":
            fp.fp[slice_dict] = data_arrays[column]*area*1e-6*1./(3600.*timeStep*1.)
        elif units == "Bq s / m^3" or units_no_space == "Bqs/m^3" or units_no_space == "Bqs/m3" or units == "Bqs/m3":
            fp.fp[slice_dict] = data_arrays[column]*area/(3600.*timeStep*1.)           
        else:
            status_log("DO NOT RECOGNISE UNITS OF {} FROM NAME INPUT (expect 'g s / m^3' or 'ppm s')".format(units),
                       error_or_warning="error")
        
        return fp

    def convert_units_ds(fp_ds, slice_dict, units, use_surface_conditions = True):
        '''
        Conversion is based on inputs units
        Input: an xarray dataset with uncoverted HiTRes footprints. 
        The following dimensions of the fp_HiTRes is expected: [time, lev, lat, lon, H_back]
        
        For each slice dictionary,
        
        If units are 'ppm s':
            Convert from [ppm s] (i.e. mu-mol/mol) units to [(mol/mol) / (mol/m2/s)].
            
            Using conversion:
                sensitivity [mu-mol/mol s] * area [m2] * molar mass [g/mol] * 1e-6 [mu] 
                    / (time [s] * release rate [g/s])
        
        If units are 'g s / m^3' (or 'gs/m3'):
            Convert from [g s/m3] units to [(mol/mol) / (mol/m2/s)].
            
            Using conversion:
                sensitivity [g s/m3] * area [m2] * RT/P [m3/mol]
                / (time [s] * release rate [g/s])
            
        Note:
            release rate is assumed to be 1. [g/s]
            molecular weight in the NAME run itself was set to 1.0, so molar mass will
            be 1.0 in this calculation as well.
        
        Optional args:
            use_surface_conditions (bool, optional) :
                If the input units are 'gs/m3', use representation surface P/T ratio of 345 rather
                than using the input meteorological data.
                Default = True
        '''
        units_no_space = units.replace(' ','')
        if units == "g s / m^3" or units == "gs/m3" or units_no_space == "gs/m^3"  or units_no_space == "gs/m3":
            if use_surface_conditions:
                molm3=345./const.R ## Surface P/T ratio we would expect over Europe (345).
            else:
                molm3=fp["press"][slice_dict].values/const.R/\
                    const.convert_temperature(fp["temp"][slice_dict].values.squeeze(),"C","K")
            fp_ds.fp_HiTRes.loc[slice_dict] = fp_ds.fp_HiTRes.loc[slice_dict].values*area/ \
                (3600.*timeStep*1.)/molm3
        elif units == "ppm s" or units_no_space == "ppms":
            fp_ds.fp_HiTRes.loc[slice_dict] = fp_ds.fp_HiTRes.loc[slice_dict].values*area*1e-6*1./(3600.*timeStep*1.)
        elif units == "Bq s / m^3" or units_no_space == "Bqs/m^3" or units_no_space == "Bqs/m3" or units == "Bqs/m3":
            fp_ds.fp_HiTRes.loc[slice_dict] = fp_ds.fp_HiTRes.loc[slice_dict].values*area/(3600.*timeStep*1.)           
        else:
            status_log("DO NOT RECOGNISE UNITS OF {} FROM NAME INPUT (expect 'g s / m^3' or 'ppm s')".format(units),
                       error_or_warning="error")
        
        return fp_ds


    if satellite:
        for t in range(len(time)):
            for l in range(len(levs)):
                slice_dict = dict(time = [t], lev = [l])
                column = t*len(levs)+l
                if units_str == 'g s / m^3' or units_str.replace(' ','') == 'gs/m^3':
                    status_log("NOT RECOMMENDED TO CREATE SATELLITE FOOTPRINTS USING CONVERTED g s / m3 UNITS. IF POSSIBILE, NAME FOOTPRINTS SHOULD BE RE-GENERATED IN UNITS OF ppm s.",
                                error_or_warning="warning")
                column = t*len(levs)+l
                fp = convert_units(fp, slice_dict, column, units_str,use_surface_conditions=use_surface_conditions)
    else:
        if species == 'CO2':
            for i in range(len(time)):
                slice_dict = dict(time = [i], lev = [0])
                fp = convert_units(fp, slice_dict, i, units_str,use_surface_conditions=use_surface_conditions) 
                for j in range(int(user_max_hour_back)):
                    slice_dict = dict(time = time[i], lev='From     0 -    40m agl', H_back=j)
                    fp = convert_units_ds(fp, slice_dict, units_str, use_surface_conditions=use_surface_conditions)
                    
        else:
            for i in range(len(time)):
                slice_dict = dict(time = [i], lev = [0])
                fp = convert_units(fp, slice_dict, i, units_str,use_surface_conditions=use_surface_conditions)
                
    if species == 'CO2':
        addfp = fp.fp_HiTRes.sum(dim='H_back')
        remfp = fp.fp - addfp
        remfpval = np.expand_dims(remfp.values, 4)
        remfpvar = xray.DataArray(remfpval,
                                  dims   = ['time','lev','lat', 'lon','H_back'],
                                  coords = {'time': fp.time, 'lev': fp.lev,
                                            'lat': fp.lat, 'lon': fp.lon,
                                            'H_back': [user_max_hour_back]})
        remfpds = remfpvar.to_dataset(name = 'fp_HiTRes')
        fp = xray.merge([fp, remfpds])
    
    fp.attrs["fp_output_units"] = units_str
    
    return fp
    
def footprint_concatenate(fields_prefix, 
                          particle_prefix = None,
                          datestr = "*",
                          met = None,
                          satellite = False,
                          time_step = None,
                          upper_level = None,
                          use_surface_conditions = True,
                          species = None,
                          user_max_hour_back=24.,
                          lifetime_hrs = None,
                          unzip=False):
    '''Given file search string, finds all fields and particle
    files, reads them and concatenates the output arrays.
    
    Args:
        fields_prefix (str): 
            prefix for fields file search string.
        particle_prefix (str, optional): 
            prefix for particle file search string.
            Default = None
        datestr (str, optional): 
            Date for particle and fields files. Default = '*'
        met (dict, optional): 
            Dictionary of met data. Default = None (empty)
        satellite (bool, optional): 
            Is it for satellite data? Default = False
        time_step (float, optional): 
            Timestep of footprint. Default = None
        upper_level (int):
            Only needed when satellite=True. Highest level number from within the NAME run for the satellite data.
            Default = None.
        use_surface_conditions (bool, optional) :
            Use default expected surface conditions for meteorological values
            if converting from gs/m3 to mol/mol / mol/m2/s units.
            P/T ratio is fixed as 345 based on typical surface conditions.
            Default = True.
        species (str,optional):
            Defaults to None which will process footprints for an inert species
            Otherwise will look in json file for lifetime and process a species-
            specific footprint
        user_max_hour_back (float, optional):
            Defaults to 24 hours. This will calculate high-time res footprints back to the
            number of hours specied here.
        lifetime_hrs (float, optional):
            Defaults to None which will process footprints for an inert species.
            Otherwise will process it for the specified loss timescale.
        unzip (bool)
            If True, unzip .nc.gz field files, defaults to False
    
    Returns:
        fp (xarray dataset): 
            Concatenated array of all field and particle files 
    
    Example:
        fp_dataset = footprint_concatenate("/dagage2/agage/metoffice/NAME_output/MY_FOOTPRINTS_FOLDER/Fields_Files/filename_prefix")
    '''

    # Find footprint files and MATCHING particle location files
    # These files are identified by their date string. Make sure this is right!
    if satellite:
        fields_files = sorted(glob.glob(fields_prefix + "*" +
                             datestr + ".txt*"))
    elif 'MixR_hourly' in fields_prefix:
        fields_files = sorted(glob.glob(fields_prefix + "*" +
                             datestr + "*.nc*"))       
    else:
        fields_files = sorted(glob.glob(fields_prefix + "*" +
                             datestr + "*.txt*"))
    
    # Search for particle files
    if 'MixR_hourly' in fields_prefix: 
        file_datestrs = [f.split(fields_prefix)[-1].split(".nc")[0].split("_")[-1] \
            for f in fields_files]
    else:
        file_datestrs = [f.split(fields_prefix)[-1].split(".txt")[0].split("_")[-1] \
                for f in fields_files]

    particle_files = []
    if particle_prefix is not None:
        for file_datestr in file_datestrs:
            
            particle_file_search_string = \
                particle_prefix + "*" + file_datestr + "*.txt*"
            particle_file_search = \
                glob.glob(particle_file_search_string)

            if particle_file_search:
                particle_files.append(particle_file_search[0])
            else:
                print("Can't find particle file " + \
                      particle_file_search_string)
                return None
                
        if len(particle_files) != len(fields_files):
            status_log("Particle files don't match fields files. SKIPPING.",
                       error_or_warning="error")
            return None
    else:
        particle_files = [None for f in fields_files]


    # Create a list of xray datasets
    fp = []
    if len(fields_files) > 0:
        if unzip:
            # unzip .nc.gz files
            # - check if the filename has .nc.gz and if there is already an unzipped version
            fields_files_unzip = [unzip(ff, return_filename=True) if '.nc.gz' in ff and
                                  [ff.split('.gz')[0] in fil for fil in fields_files].count(True)==1
                                  else None if '.nc.gz' in ff
                                  else ff for ff in fields_files]

            # filter out any Nones and any duplicates
            fields_files_unzip = [ff for ff in fields_files_unzip if ff is not None]
            fields_files_unzip = list(dict.fromkeys(fields_files_unzip))
        
        else:
            # list of field files, excluding .nc.gz files
            fields_files_unzip = [ff for ff in fields_files if '.nc.gz' not in ff]
        
        # extract footprints
        for fields_file, particle_file in zip(fields_files_unzip, particle_files):
            fp.append(footprint_array(fields_file,
                                      particle_file = particle_file,
                                      met = met,
                                      satellite = satellite,
                                      time_step = time_step,
                                      upper_level = upper_level,
                                      use_surface_conditions = use_surface_conditions,
                                      species = species,
                                      user_max_hour_back=user_max_hour_back,
                                      lifetime_hrs = lifetime_hrs))
        if unzip:
            # remove unzipped files to save space
            [os.remove(fields_file) for ff, fields_file in
             enumerate(fields_files_unzip) if '.nc.gz' in fields_files[ff]]
        
    # Concatenate
    if len(fp) > 0:
        fp = xray.concat(fp, "time")
    else:
        fp = None

    return fp

def write_netcdf(fp, outfile,
            temperature=None, pressure=None,
            wind_speed=None, wind_direction=None,
            PBLH=None, fp_HiTRes_inc=False, varname="fp", varname2="fp_HiTRes",
            release_lon = None, release_lat = None,
            particle_locations=None, particle_mean_age = None, particle_heights=None,
            global_attributes = {}, lapse_rate=None, lapse_error=None, units = None):
    '''Writes netCDF with footprints, particle locations, meteorology and release locations.
    
    Args:
        fp (xarray dataset): 
            Array of all footprint and particle fields in (mol/mol)/(mol/m2/s)
        lons (array): 
            1D array of longitudes
        lats (array): 
            1D array of latitudes
        levs (???array or list):
            TODO: Clarify inputs
            1D array of particle levels
        time (???str or datetime object): 
            TODO: Clarify inputs
            Timestamps for footprints
        outfile (str): 
            Name of output file
        temperature (array, optional): 
            Input temperature variable in C. Default=None
        pressure (array, optional): 
            Input pressure variable in hPa. Default = None
        wind_speed (array, optional): 
            Input wind speed variable in m/s. Default = None
        wind_direction (array, optional): 
            Input wind direction variable in degrees. 
            Default = None
        PBLH (array, optional): 
            Input planetary boundary layer height in m. 
            Default = None
        fp_HiTRes_inc (str, optional):
            True if writing out fp_HiTRes
        varname (str, optional): 
            Name of output footprint variable. Default = 'fp'
        release_lon (array, optional): 
            Array of release longitude locations. 
            Default = None
        release_lat (array, optional): 
            Array of release latitude locations. 
            Default = None
        particle_locations (xarray dataset, optional): 
            Number of particles at locations at each boundary
            Default=None
        particle_heights (array, optional): 
            Heights of particles leaving boundary.
            Default=None
        global_attributes (dictionary, optional): 
            Dictionary of global attributes.
            Default={} (empty)
        lapse_rate (array, optional): 
            Potential temperature gradient from 60 - 300 m in K/km. 
            Default=None
        lapse_error (array, optional): 
            Error in potential temperature gradient in K/km
            
    Returns:
        None.
        Writes netCDF with footprints, particle locations, meteorology and release locations
        
    Note: 
        netCDF4 is required, as we make use of compression.    
    '''
    
    fp_attr_loss = fp.fp.attrs["loss_lifetime_hrs"]
    
    lons = fp.lon.values.squeeze()
    lats = fp.lat.values.squeeze()
    levs = fp.lev.values
    time = fp.time.to_pandas().index.to_pydatetime()
    fp_ds = fp.fp.transpose("lat", "lon", "time").values.squeeze()
    
    if fp_HiTRes_inc == True:
        H_back = fp.H_back.values
        fp_HiTRes = fp.fp_HiTRes.transpose("lat", "lon", "time", "H_back").values.squeeze()
    
    time_seconds, time_reference = time2sec(time)   
    
    #Write NetCDF file
    ncF=Dataset(outfile, 'w')
    ncF.createDimension('time', len(time))
    ncF.createDimension('lon', len(lons))
    ncF.createDimension('lat', len(lats))
    ncF.createDimension('lev', 1)
    if fp_HiTRes_inc == True:
        ncF.createDimension('H_back', len(H_back))
    
    # pass any global attributes in fp to the netcdf file
    for key in list(global_attributes.keys()):
        ncF.__setattr__(key,global_attributes[key])
    ncF.__setattr__("author", getpass.getuser())
    ncF.__setattr__("created", str(dt.datetime.now()))
    
    
    nctime=ncF.createVariable('time', 'd', ('time',))
    nclon=ncF.createVariable('lon', 'f', ('lon',))
    nclat=ncF.createVariable('lat', 'f', ('lat',))
    nclev=ncF.createVariable('lev', 'S1', ('lev',))
    ncfp=ncF.createVariable(varname, 'f', ('lat', 'lon', 'time'), zlib = True,
                            least_significant_digit = 5)
    
    if fp_HiTRes_inc == True:
        ncH_back=ncF.createVariable('H_back', 'f', ('H_back',))
        ncfp_HiTRes=ncF.createVariable(varname2, 'f', ('lat', 'lon', 'time', 'H_back'), zlib = True,
                            least_significant_digit = 5)
    
    nctime[:]=time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + np.str(time_reference)
    nctime.label='left'
    nctime.period = str(timestep_for_output) + " hours"
    nctime.comment = 'time stamp corresponds to the beginning of each averaging period'
    nctime.calendar='gregorian'

    nclon[:]=lons
    nclon.units='Degrees_east'
    
    nclat[:]=lats
    nclat.units='Degrees_north'

    nclev[:]=np.array(levs)
    
    if fp_HiTRes_inc == True:
        ncH_back[:]=H_back
        ncH_back.units='Hours'
        ncH_back.long_name='Hours back from release time'
        
    
    ncfp[:, :, :]=fp_ds
    if units == None:
        ncfp.units='(mol/mol)/(mol/m2/s)'
    else:
        ncfp.units = units

    if fp_HiTRes_inc == True:
        ncfp_HiTRes[:, :, :]=fp_HiTRes
        if units == None:
            ncfp_HiTRes.units='(mol/mol)/(mol/m2/s)'
        else:
            ncfp_HiTRes.units = units

        
    ncfp.loss_lifetime_hrs = fp_attr_loss

    if fp_HiTRes_inc == True:
        ncfp_HiTRes.loss_lifetime_hrs = fp_attr_loss

    if temperature is not None:
        nctemp=ncF.createVariable('temperature', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        nctemp[:]=temperature
        nctemp.units='C'

    if pressure is not None:
        ncpress=ncF.createVariable('pressure', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncpress[:]=pressure/100.
        ncpress.units='hPa'

    if wind_speed is not None:
        ncwind_speed=ncF.createVariable('wind_speed', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncwind_speed[:]=wind_speed
        ncwind_speed.units='m/s'

    if wind_direction is not None:
        ncwind_direction=ncF.createVariable('wind_direction', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncwind_direction[:]=wind_direction
        ncwind_direction.units='degrees'

    if PBLH is not None:
        ncPBLH=ncF.createVariable('PBLH', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncPBLH[:]=PBLH
        ncPBLH.units='m'

    if release_lon is not None:
        ncRlon=ncF.createVariable('release_lon', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncRlon[:]=release_lon
        ncRlon.units='Degrees_east'

    if release_lat is not None:
        ncRlat=ncF.createVariable('release_lat', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        ncRlat[:]=release_lat
        ncRlat.units='Degrees_north'

    if particle_locations is not None:
        ncF.createDimension('height', len(particle_heights))
        ncHeight=ncF.createVariable('height', 'f',
                                   ('height',),
                                    zlib = True, least_significant_digit = 4)
        ncPartN=ncF.createVariable('particle_locations_n', 'f',
                                   ('height', 'lon', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncPartE=ncF.createVariable('particle_locations_e', 'f',
                                   ('height', 'lat', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncPartS=ncF.createVariable('particle_locations_s', 'f',
                                   ('height', 'lon', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncPartW=ncF.createVariable('particle_locations_w', 'f',
                                   ('height', 'lat', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncmeanageN=ncF.createVariable('mean_age_particles_n', 'f',
                                   ('height', 'lon', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncmeanageE=ncF.createVariable('mean_age_particles_e', 'f',
                                   ('height', 'lat', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncmeanageS=ncF.createVariable('mean_age_particles_s', 'f',
                                   ('height', 'lon', 'time'),
                                    zlib = True, least_significant_digit = 7)
        ncmeanageW=ncF.createVariable('mean_age_particles_w', 'f',
                                   ('height', 'lat', 'time'),
                                    zlib = True, least_significant_digit = 7)
        
        ncHeight[:]=particle_heights
        ncPartN[:, :, :]=particle_locations["N"]
        ncPartE[:, :, :]=particle_locations["E"]
        ncPartS[:, :, :]=particle_locations["S"]
        ncPartW[:, :, :]=particle_locations["W"]
        ncmeanageN[:, :, :]=particle_mean_age["N"]
        ncmeanageE[:, :, :]=particle_mean_age["E"]
        ncmeanageS[:, :, :]=particle_mean_age["S"]
        ncmeanageW[:, :, :]=particle_mean_age["W"]
        ncPartN.units=''
        ncPartN.long_name='Fraction of total particles leaving domain (N side)'
        ncPartE.units=''
        ncPartE.long_name='Fraction of total particles leaving domain (E side)'
        ncPartS.units=''
        ncPartS.long_name='Fraction of total particles leaving domain (S side)'
        ncPartW.units=''
        ncPartW.long_name='Fraction of total particles leaving domain (W side)'
        ncmeanageN.units='hrs'
        ncmeanageN.long_name='Mean age of particles leaving domain (N side)'
        ncmeanageE.units='hrs'
        ncmeanageE.long_name='Mean age of particles leaving domain (E side)'
        ncmeanageS.units='hrs'
        ncmeanageS.long_name='Mean age of particles leaving domain (S side)'
        ncmeanageW.units='hrs'
        ncmeanageW.long_name='Mean age of particles leaving domain (W side)'    
    if lapse_rate is not None:
        nclapse=ncF.createVariable('lapse_rate', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        nclapse[:]=lapse_rate
        nclapse.long_name="Potential temperature gradient from 60 - 300 m"
        nclapse.units='K/km'
        
    if lapse_error is not None:
        nclerror=ncF.createVariable('lapse_error', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        nclerror[:]=lapse_error
        nclerror.long_name="Error in potential temperature gradient"
        nclerror.units='K/km'
    
    ncF.close()
    status_log("Written... " + os.path.split(outfile)[1])


def satellite_vertical_profile(fp, satellite_obs_file, max_level):
    '''Do weighted average by satellite averaging kernel and
    pressure weight. One time point only. Expects xray.dataset
    with one time point and N vertical levels.
    
    Args:
        fp (xarray dataset):
            footprint for ONE time point. N levels in lev dimension with 
            footprints and particle locations defined at each level
        satellite_obs_file (str):
            Filename of NetCDF-format satellite observation file. 
            One per time step.
        max_level (int):
            Maximum vertical level of the retrieval that is included.
            (maximum for GOSAT to include all levels from model is 20). 
            The remaining levels are set equal to the a priori value.
                
    Returns:
        fp (dictionary):
            Footprint for column totals of satellite observations.
    
    Note:
        Requirements for running processing satellite footprints
        
        - General file naming: 
            Files/folders need to have a date string of the form
            YYYYMMDD-II at the end of the file/folder name,
            where II is a two-digit index corresponding to the
            sequence of footprints throughout the day. This modifies the usual
            naming convention for surface sites, which are only YYYYMMDD.
        - Met: 
            The Met/ folder MUST contain a set of subfolders which are labeled
            with the YYYYMMDD-II time stamp. Each of these subfolders should then
            contain nLev Met output files, where nLev are the number of vertical
            levels in the NAME output. This is required so that we know
            what pressure/temperature each NAME footprint for each vertical level
            corresponds to. The files can be named in any way, but the file names
            must be in ascending order when sorted (i.e. blah_01.txt.gz, blah_02.txt.gz).
        - Fields files:
            There must be ONE fields file per time stamp, labeled at the end
            of the file string with YYYYMMDD-II (e.g. blah_YYYYMMDD-II.txt.gz).
            The fields file must contain exactly nLev columns, one for each vertical
            level, in ascending order, left to right.
        - Particle location files:
            There must be ONE fields file per time stamp, labeled at the end
            of the file string with YYYYMMDD-II (e.g. blah_YYYYMMDD-II.txt.gz).
            The ID column in this file must contain nLev values in ascending order
            specifying the vertical level that each particle belongs to.    
    '''
    
    if max_level is None:
        status_log("MAX LEVEL REQUIRED TO PROCESS SATELLITE FOOTPRINTS",
                   error_or_warning="error")
        return None
    
    status_log("Reading satellite obs file: " + satellite_obs_file)
    with xray.open_dataset(satellite_obs_file) as f:
        sat = f.load()
    
    ntime = len(fp.time)
    
    for t in range(ntime):
        if np.abs(sat.lon.values[t] - fp.release_lon.values[t,0]) > 1.:
            status_log("Satellite longitude doesn't match footprints",
                       error_or_warning="error")
        if np.abs(sat.lat.values[t] - fp.release_lat.values[t,0]) > 1:
            status_log("Satellite latitude doesn't match footprints",
                       error_or_warning="error")
        
        if not np.allclose((fp.pl_n[t].sum() + fp.pl_e[t].sum() + \
                            fp.pl_s[t].sum() + fp.pl_w[t].sum()), \
                            len(fp.lev)):
            status_log("Particle histograms dont add up to 1 (or nlev)",
                       error_or_warning="error")
            return None

    if ntime == 1:
        if np.abs(sat.time.values[0] - fp.time.values[0]).astype(int) > 60*1e9:
            status_log("Satellite time doesn't match footprints",
                       error_or_warning="warning")
    else:
         status_log("Number of satellite time points > 1 so unable to rely on time interpolated from fields files to be correct. Not checking time within obs file against fields file time.",
                       error_or_warning="status")
    
        

    # Change timestamp to that from obs file
    #  because NAME output only has 1 minute resolution

    fp["time"] = sat.time.values

    # Interpolate pressure levels
    variables = ["fp", "pl_n", "pl_e", "pl_s", "pl_w"]
    out = {}
    lower_levels =  list(range(0,max_level))
     
    # Weight levels using pressure weight and averaging kernel
    
    fp_combined = fp[dict(lev = [0])].copy()
    
    for t in range(ntime):
        sum_ak_pw = np.sum(sat.pressure_weights.values[t][lower_levels] * \
                           sat.xch4_averaging_kernel.values[t][lower_levels])
        sum_particle_count = 0.
    
        for variable in variables:
            # get dimensions based on level 0
            fp_lev = int(np.abs(fp.press[{"time":t}] - \
                                sat.pressure_levels[dict(time = t,lev = 0)]*100.).argmin())
            var = fp[variable][{"time":t,"lev": fp_lev}].values * \
                  sat.pressure_weights.values[t][0] * \
                  sat.xch4_averaging_kernel.values[t][0]

            if t == 0:
                out[variable] = np.zeros((ntime,1)+var.shape)

            out[variable][t] = var.reshape((1,) + var.shape)
            # run levels 1 - max_level (note it is reindexed to 0, therefore use lev+1)
            for lev, press in enumerate(sat.pressure_levels.values[t][1:max_level]):
                    fp_lev = np.abs(fp.press.values[t] - press*100.).argmin()
                    var = fp[variable][{"time":t,"lev": fp_lev}].values * \
                          sat.pressure_weights.values[t][lev+1] * \
                          sat.xch4_averaging_kernel.values[t][lev+1]
                    out[variable][t] += var.reshape((1,) + var.shape)
            if variable[0:2] == "pl":
                sum_particle_count += out[variable][t].sum()
    
        # Check whether particle sum makes sense
        if not np.allclose(sum_particle_count, sum_ak_pw):
            status_log("ERROR: Particle fractions don't match averaging_kernel * " + \
                       "pressure_weight",
                       error_or_warning = "error")
            return None
    
    for variable in variables:
        fp_combined[variable].values = out[variable]
    
    fp_combined.attrs["max_level"] = max_level    
    fp_combined["lev"] = np.array(["column"])
    
    return fp_combined


def status_log(message,
               directory = None,
               error_or_warning = "status",
               print_to_screen = True):
    '''Write a log of an error or a warning to file. Will append to a file 
    
    Args:
        message(str):
            Error message
        directory (str, optional): 
            Directory of error log. Default = None
        error_or_warning (str):
            Message type. Default = "status"
        print_to_screen (bool): 
            Print to screen? Default=True
            
    Returns:
        None.
        Writes error to file
        
    Example:
        status_log("Everything has gone wrong")
            
    '''

    if directory is None:
        directory = directory_status_log

    date = dt.datetime.now()
    year = date.year
    month = date.month
    day = date.day
    
    if error_or_warning == "error":
        fname = os.path.join(directory,
                             "PROCESS_ERROR_LOG_%04d%02d%02d.txt" % 
                             (year, month, day))
    else:
        fname = os.path.join(directory,
                             "PROCESS_STATUS_LOG_%04d%02d%02d.txt" % 
                             (year, month, day))

    with open(fname, "a") as f:
        if error_or_warning == "warning":
            f.write("WARNING: %s: %s\r\n" %(str(date), message))
        elif error_or_warning == "error":
            f.write("ERROR: %s: %s\r\n" %(str(date), message))
        else:
            f.write("%s: %s\r\n" %(str(date), message))
            
    if print_to_screen:
        if error_or_warning == "warning":
            print("WARNING: " + message)
        elif error_or_warning == "error":
            print("ERROR: " + message)
        else:
            print(message)

def process_basic(fields_folder, outfile):
    """Basic processing with no meteorology or particle files
    NOT recommended, but helpful for a quick check
    
    Args:
        field_folder (str):
            folder containing fields data
        outfile (str):
            Name of outfile
    
    Returns:
        None.
        Writes outfile.
    
    """
    
    fp = footprint_concatenate(fields_folder)
    write_netcdf(fp, outfile)

def process(domain, site, height, year, month, 
            met_model = None,
            base_dir = os.path.join(data_path,"NAME_output/"),
            process_dir = os.path.join(lpdm_path,"fp_NAME/"),
            fields_folder = "MixR_files",
            particles_folder = "Particle_files",
            met_folder = ["Met_daily", "Met"],
            processed_folder = "/work/chxmr/shared/LPDM/fp_NAME/",
            force_met_empty = False,
            use_surface_conditions = True,
            satellite = False,
            obs_folder = "Observations",
            upper_level = None,
            max_level = None,
            force_update = False,
            perturbed_folder=None,
            vertical_profile=False,
            transport_model="NAME",
            units = None,
            species = None,
            user_max_hour_back = 24.):
    
    '''Process a single month of footprints for a given domain, site, height,
    year, month. 
    
    If you want to process all files for a given domain + site
    use process_all.
    
    This routine finds all fields files and particle location files which match
    the timestamp in the file name. The date strings are stored in the datestr
    variable, and these strings must match exactly for the fields and the particle
    locations. Hourly netcdf output must be .nc and not .nc.gz!
    
    At the moment, the entire Met folder is read each time.
    
    Args:
        domain (str):
            Domain of interest
        site (str):
            Observation site
        height (str):
            Height of observation, e.g. "10magl"
        year (int):
            The year
        month (int):
            1,2,3,4,5,6,7,8,9,10,11,12
        met_model (str):
            Met model used to run NAME
            Default is None and implies the standard global met model
            Alternates include 'UKV' and this met model must be in the outputfolder NAME
        base_dir (str, optional):
            Base directory containing NAME output
            Default (i.e. on BP1) ="/work/chxmr/shared/NAME_output/",
        process_dir (str, optional):
            Base directory for where processed netcdf files go
            Default (i.e. on BP1)="/work/chxmr/shared/LPDM/fp_NAME/",        
        fields_folder (str, optional):
            Folder containing fields data
            Default="Fields_files",
        particles_folder (str, optional):
            Folder containing particles data.
            Default = "Particle_files",
        met_folder (str or list, optional):
            Folder(s) containing met data
            Default = "Met"
        force_met_empty (bool, optional):
             Force the met data to be empty?
             Default = False.
        use_surface_conditions (bool, optional) :
            Use default expected surface conditions for meteorological values
            if converting from gs/m3 to mol/mol / mol/m2/s units.
            This involves fixing the P/T ratio to 345 based on typical surface conditions 
            in Europe and not using the meteorlogical data extracted from NAME release
            points.
            Default = True.
        satellite (bool, optional):
            NAME run was performed over satellite data.
            This means satellite_vertical_profile function will be used to combined data
            from multiple levels for each point as defined by upper_level parameter.
            Default = False.
        obs_folder (str, optional) :
            Folder for netCDF format observations. Won't be contained in the original
            NAME results folder and will need to be copied into the subfolder defined
            by this parameter.
            Default = "Observations"
        upper_level (int/None):
            Only needed when satellite=True. Highest level number from within the 
            NAME run for the satellite data.
            Default = None.
        max_level (int, optional):
            Only needed when satellite=True. The max level to
            process the footprints.
            Any levels above are replaced by the prior profile.
            This can be <= upper_level.
            If max_level is not specified, this will match upper_level.
            Default = None.
        force_update (bool, optional):
            By default, any existing netCDF files are NOT overwritten.
            To explicitly over-write a file, set force_update = True
        perturbed_folder (str, optional)
            Process a subfolder which has all of the required folders
            (Fields_files, etc). This is for perturbed parameter ensembles for 
            a particular site. E.g. for 
            EUROPE_BSD_110magl/Perturbed/PARAMETERNAME_VALUE
            you'd set: perturbed_folder = "Perturbed/PARAMETERNAME_VALUE".
            Default = None.
        vertical_profile (bool, optional):
            If set to True will look for vertical potential temperature met file
            and incorporate into footprint file. 
            This is a separate file from the normal met file, and is not mandatory.
            Default=False.
        transport_model (str, optional): 
            Defaults to "NAME". If "STILT", reads footprints in the
            ncdf format created by the STILT model. Other values are invalid. 
            Notset up to read satellite column footprints from STILT format.
            Default="NAME".
        species (str,optional):
            Defaults to None which will process footprints for an inert species
            Otherwise will look in json file for lifetime and process a species-
            specific footprint with the species name in the .nc file
        user_max_hour_back (float, required when species = 'CO2'):
            Defaults to 24 hours. This is the maximum amount of time back from release time that 
            hourly footprints are calculated.
        
    Returns:
        None.
        This routine outputs a copy of the xarray dataset that is written to file.
    
    '''
    # define the path for processed files, check it exists, and check that the user has write permission
    full_out_path = os.path.join(process_dir, domain)
    if not os.path.isdir(full_out_path):
        os.makedirs(full_out_path)
    
    if not os.access(full_out_path, os.W_OK):
        print('No write access for output directory, provide a different option for process_dir or change the directory permissions')
        return None
    
    global directory_status_log
    
    if met_model is None:
        subfolder = base_dir + domain + "_" + site + "_" + height + "/"
    else:
        subfolder = base_dir + domain + "_" + met_model + "_" + site + "_" + height + "/"
    directory_status_log = subfolder
    
    # do not run hourly processing if species specified but it is not pointing to the MixR_hourly directory
    # do not run hourly processing if species not specified but pointing to MixR_hourly as will be inefficient and should be MixR_files
    # do not run hourly processing if species is long-lived
    
    if species is not None and 'MixR_hourly' not in fields_folder:
        print("A species was specified so fields_folder changed to MixR_hourly")
        fields_folder = "MixR_hourly"
    elif species is None and 'MixR_hourly' in fields_folder:
        print("No species was specified. For efficiency, fields_folder should be MixR_files or Fields_files, not MixR_hourly")
        return
    
    if species is not None:
        with open(os.path.join(acrg_path,"acrg_species_info.json")) as f:
            species_info=json.load(f)
        species = obs.read.synonyms(species, species_info)
        if 'lifetime' in species_info[species].keys():
            lifetime = species_info[species]["lifetime"]
            if type(lifetime) == list:
                lifetime_hrs = acrg_time.convert.convert_to_hours(lifetime[month-1])
            else:
                lifetime_hrs = acrg_time.convert.convert_to_hours(lifetime)
            print('Lifetime of species in hours is', lifetime_hrs)
            
            if lifetime_hrs > 1440:
                print("This is a long-lived species. For efficiency, folder has been changed to MixR_files (inert species).")
                lifetime_hrs = None
                fields_folder = "MixR_files"
        elif species == 'CO2':
            print('Preparing footprints with high temporal resolution for CO2')
            lifetime_hrs = None
        else:
            print('No lifetime has been defined in species_info.json')
            return
    else:
        print('Running for an inert substance')
        lifetime_hrs = None
    

    # Check that there are no errors from the NAME run
    input_folder = subfolder + 'Input_files/'
    if os.path.isdir(input_folder):
        error_files = 'BackRun_' + domain + '_' + site + '_' + height + '_' + str(year) + str(month).zfill(2)
        error_days = []
        for file_name in os.listdir(input_folder):
            if file_name.startswith(error_files) and \
            file_name.endswith('Error.txt') and \
            os.stat(input_folder+file_name).st_size != 0: 
                error_days.append(re.search(r"[0-9]{8}(\S+)", file_name).group(0))
                error_days.sort()
                num_days = len(error_days)

        if len(error_days) > 0:
            raise Exception('This month cannot be processed as there are '+str(num_days)+' days with with errors: '+str(error_days))
        
   # Check for existance of subfolder
    if not os.path.isdir(subfolder):
        raise Exception("Subfolder: {} does not exist.\nExpect NAME output folder of format: domain_site_height".format(subfolder))
    
    if perturbed_folder != None:
        if perturbed_folder[-1] == "/":
            subfolder += perturbed_folder
        else:
            subfolder += perturbed_folder + "/"
    
    # Check that the specified transport model is valid.
    if transport_model == "STILT":
        if satellite:
            status_log("stiltfoot_array is not set up for satellite data!" +\
                       " Levels will be wrong!", error_or_warning="error")
        if not force_met_empty:
            status_log("STILT neither provides nor requires met information" +\
                       " to interpret footprints. Met will probably be set" +\
                       " to default values. Don't rely on these values!")
    elif transport_model != "NAME" and transport_model !="NAMEUKV":
        status_log(transport_model + " is not a valid transport model!" +\
                   " Unable to read footprint information!", 
                   error_or_warning="error")
    
    # Check for manual timestep (if footprints are for < 1min,
    # which is the min resolution of the NAME output file)
    timestep_file = os.path.join(subfolder,"time_step.txt")
    
    if os.path.exists(timestep_file):
        with open(timestep_file) as f:
            timeStep = float(f.read())
            if satellite:
                timeStep = timeStep/3600. # Assume time step input from satellite is in seconds rather than hours.
    elif satellite == True:
        status_log("A time step file called 'time_step.txt' must be specified for satellite data. "
                   +"The file should just contain the number of retrieval seconds for the NAME run in seconds. "
                   +"This is necessary because the time step cannot be accurately determined from the "
                   +"input fields file details.", 
                   error_or_warning="error")
        return None
    else:
        timeStep = None

    if satellite:
        if upper_level is None:
            status_log("Upper Level must be specified for satellite data.",
                   error_or_warning="error")
            return None
        elif max_level is None:
            status_log("No max level specified, so set to match upper level of footprints.",
                   error_or_warning="status")
            max_level = upper_level

    if use_surface_conditions:
        status_log("Setting use_surface_conditions to True means that a representative Pressure/Temperature "
                   +"ratio will be used instead of meteorological data from the NAME output. "
                   +"This is only applied when converting units from gs/m3",
                   error_or_warning="warning")


    # Get date strings for file search
    if satellite:
        file_search_string = subfolder + fields_folder +"/*" + str(year) \
                        + str(month).zfill(2) + "*.txt*"
        fields_files = glob.glob(file_search_string)
        datestrs = sorted([fi.split("_")[-1].split(".")[0] for fi in fields_files])

        # Check that we've found something
        if len(datestrs) == 0:
            status_log("can't find files in " + file_search_string,
                       error_or_warning="error")
            return None
    else:
        if transport_model == "STILT":
            datestrs = ["stilt" + str(year) + "x" + str(month).zfill(2) + "x"]
        else:
            datestrs = [str(year) + str(month).zfill(2)]

    # Output filename
    if met_model is None:
        if species is None:
            outfile = os.path.join(full_out_path, site + "-" + height +  \
                    "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc")
        else:
            outfile = os.path.join(full_out_path, site + "-" + height + "_" + species.lower() + \
                        "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc")
    else:
        if species is None:
            outfile = os.path.join(full_out_path, site + "-" + height + "_" + met_model + \
                    "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc")
        else:
            outfile = os.path.join(full_out_path, site + "-" + height + "_" + met_model + "_" + species.lower() + \
                        "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc")

    
    # Check whether outfile needs updating
    if not force_update:
        status_log("Testing whether file exists or needs updating: " + outfile)
        if os.path.exists(outfile):
            ncf = Dataset(outfile)
            time_test = sec2time(ncf.variables["time"][:],
                                 ncf.variables["time"].units[14:])
            #Find maximum day (minus 0.1s because last time is usually midnight on 1st of the next month)
            maxday = (max(time_test) - dt.timedelta(seconds = 0.1)).day
            ncf.close()
            if satellite:
                days = [int(datestr.split('-')[0][-2:]) for datestr in datestrs]
                if maxday >= max(days):
                    return None
            else:
                if transport_model == "STILT":
                    stilt_files = glob.glob(subfolder + fields_folder + "/stilt" + \
                                            str(year) + "x" + str(month).zfill(2) + "*.nc")
                    days = [int(os.path.split(stilt_file)[1].split("x")[2]) \
                            for stilt_file in stilt_files]
                elif fields_folder == "MixR_hourly":
                    fields_files = glob.glob(subfolder + fields_folder + "/*" + \
                                             datestrs[0] + "*.nc*")
                    days = [int(os.path.split(fields_file)[1].split("_")[-1][6:8]) \
                            for fields_file in fields_files]
                else:
                    fields_files = glob.glob(subfolder + fields_folder + "/*" + \
                                             datestrs[0] + "*.txt*")
                    days = [int(os.path.split(fields_file)[1].split("_")[-1][6:8]) \
                            for fields_file in fields_files]

                if maxday >= max(days):
                    return None

    fp = []
     
    for datestr in datestrs:

        status_log("Looking for files with date string: " + datestr + " in " + \
                   subfolder)

        # Get Met files
        if force_met_empty is not True:
            if type(met_folder) == list:
                met_files = []
                for metf in met_folder:
                    if satellite:
                        met_search_str = subfolder + metf + "/*" + datestr + "/*.txt*"
                    else:
                        if metf == "Met":
                            met_search_str = subfolder + metf + "/*" + str(year) + "*.txt*"
                        else:    
                            met_search_str = subfolder + metf + "/*.txt*"
                    met_files = met_files + sorted(glob.glob(met_search_str))
            else:
                if satellite:
                    met_search_str = subfolder + met_folder + "/*" + datestr + "/*.txt*"
                else:
                    met_search_str = subfolder + met_folder + "/*.txt*"
                met_files = sorted(glob.glob(met_search_str))
           
            if len(met_files) == 0:        
                raise FileNotFoundError(f"Can't find MET files: {met_search_str}")
                
            else:
                met = read_met(met_files,satellite=satellite)
        else:
            met = None
        
            # Get footprints
        if transport_model == "STILT":
            fp_file = stiltfoot_array(subfolder + fields_folder + "/" + datestr, 
                                      met=met, satellite=satellite,
                                      time_step=timeStep)
        else: # for NAME
            fields_prefix = subfolder + fields_folder + "/"
            if particles_folder is not None:
                particles_prefix = subfolder + particles_folder + "/"
            else:
                particles_prefix = None

            fp_file = footprint_concatenate(fields_prefix,
                                            datestr = datestr, met = met,
                                            particle_prefix = particles_prefix,
                                            satellite = satellite,
                                            time_step = timeStep,
                                            upper_level = upper_level,
                                            use_surface_conditions=use_surface_conditions,
                                            species = species,
                                            lifetime_hrs = lifetime_hrs,
                                            user_max_hour_back = user_max_hour_back)
           
        # Do satellite process
        if satellite:
            obs_path = os.path.join(subfolder,obs_folder)
            search_str = "*{label}_*.nc".format(label=datestr)
            search_str = os.path.join(obs_path,search_str)
            satellite_obs_file = glob.glob(search_str)
            
            if len(satellite_obs_file) != 1:
                status_log("There must be exactly one matching satellite " + 
                           "file in the {}/ folder".format(obs_folder),
                           error_or_warning="error")
                status_log("Files: " + ','.join(satellite_obs_file),
                           error_or_warning="error")
                return None

            fp_file = satellite_vertical_profile(fp_file,
                                                 satellite_obs_file[0], max_level = max_level)  
            
            if fp_file is None:
                return
                 
        if fp_file is not None:
            fp.append(fp_file)
                       
        
            
    if len(fp) > 0:
        
        # Concatentate
        fp = xray.concat(fp, "time")

        # ONLY OUTPUT FIRST LEVEL FOR NOW
        if len(fp.lev) > 1:
            fp = fp[dict(lev = [0])]
            status_log("ONLY OUTPUTTING FIRST LEVEL!", error_or_warning = "warning")
        fp = fp.squeeze(dim = "lev")
        
        
        # REMOVE ANY ZERO FOOTPRINTS        
        fp_nonzero = (fp.fp.sum(["lon", "lat"])).values.squeeze() > 0.
        pl_nonzero = (fp.pl_n.sum(["lon", "height"]) + \
                      fp.pl_e.sum(["lat", "height"]) + \
                      fp.pl_s.sum(["lon", "height"]) + \
                      fp.pl_w.sum(["lat", "height"])).values.squeeze() > 0.
        indices_nonzero = np.where(fp_nonzero + pl_nonzero)[0]
        fp = fp[dict(time = indices_nonzero)]
        
        # Get vertical profile met file
        if vertical_profile is True:
            vp_search_str = subfolder + "vertical_profile" + "/*.txt*"     
            vp_files = sorted(glob.glob(vp_search_str))
        
            if len(vp_files) == 0:
                status_log("Can't file Vertical Profile files: " + vp_search_str,
                           error_or_warning="error")
                return None
            else:
                    vp_met = process_vertical_profile(vp_files[0])
                    
                    # Merge vetical profile met into footprint file with same time index
                    vp_reindex = vp_met.reindex_like(fp,"nearest", tolerance =None)
                    lapse_in = vp_reindex.theta_slope.values
                    lapse_error_in=vp_reindex.slope_error.values
                                      
        else:
            lapse_in=None
            lapse_error_in=None
        
        #Adding Global Attributes to fp file
        
        #Using the first datestring to select a file to gather information about the model.
        #This information is added to the attributes of the processed file.
        #Currently works for NAME or NAMEUKV
        if transport_model == "NAME":
            if fields_folder == "MixR_files" or fields_folder == "MixR_hourly":
                name_info_file_str = os.path.join(subfolder,'MixR_files/*'+datestrs[0]+'*')
                name_info_file = glob.glob(name_info_file_str)[0]
                model_info = extract_file_lines(name_info_file)[0]
            elif fields_folder == "Fields_files":
                name_info_file_str = os.path.join(subfolder,'Fields_files/*'+datestrs[0]+'*')
                name_info_file = glob.glob(name_info_file_str)[0]
                model_info = extract_file_lines(name_info_file)[0]
            else:
                model_info = "Version unknown"
        else:
            model_info = "Version unknown"
        
        if fields_folder == "MixR_files" or fields_folder == "MixR_hourly":
            model_output = "Mixing Ratio"
        elif fields_folder == "Fields_files":
            model_output = "Air concentration"
 
        if species is None:
            output_species = "Inert"
        else:
            output_species = species

        fp.attrs["fp_output"] = model_output 
        fp.attrs["species"] = output_species
        if species is not None and species!='CO2' and lifetime_hrs is not None:
            fp.attrs["species_lifetime_hrs"] = lifetime_hrs           
        fp.attrs["model"] = transport_model
        if met_model is None:
            fp.attrs["met_model"] = "Global"
        else:
            fp.attrs["met_model"] = met_model
        fp.attrs["output_folder"] = fields_folder
        fp.attrs["model_version"] = model_info
        fp.attrs["domain"] = domain
        fp.attrs["site"] = site
        fp.attrs["inlet_height"] = height
        fp.attrs["ACRG_repository_version"] = code_version()
       
        #Write netCDF file
        #######################################
        
        # Define particle locations and mean age dictionaries (annoying)
        if "pl_n" in list(fp.keys()):
            pl = {"N": fp.pl_n.transpose("height", "lon", "time").values.squeeze(),
                  "E": fp.pl_e.transpose("height", "lat", "time").values.squeeze(),
                  "S": fp.pl_s.transpose("height", "lon", "time").values.squeeze(),
                  "W": fp.pl_w.transpose("height", "lat", "time").values.squeeze()}
            pl_ma = {"N": fp.mean_age_n.transpose("height", "lon", "time").values.squeeze(),
                     "E": fp.mean_age_e.transpose("height", "lat", "time").values.squeeze(),
                     "S": fp.mean_age_s.transpose("height", "lon", "time").values.squeeze(),
                     "W": fp.mean_age_w.transpose("height", "lat", "time").values.squeeze()}
            height_out = fp.height.values.squeeze()
        else:
            pl = None
            pl_ma = None
            height_out = None

        status_log("Writing file: " + outfile, print_to_screen=False)
        
        # Write outputs
        if species == 'CO2':
            write_netcdf(fp,
                             outfile,
                             temperature=fp["temp"].values.squeeze(),
                             pressure=fp["press"].values.squeeze(),
                             wind_speed=fp["wind"].values.squeeze(),
                             wind_direction=fp["wind_direction"].values.squeeze(),
                             PBLH=fp["PBLH"].values.squeeze(),
                             release_lon=fp["release_lon"].values.squeeze(),
                             release_lat=fp["release_lat"].values.squeeze(),
                             particle_locations = pl,
                             particle_mean_age = pl_ma,
                             particle_heights = height_out,
                             global_attributes = fp.attrs,
                             lapse_rate = lapse_in,
                             lapse_error = lapse_error_in, units = units,
                             fp_HiTRes_inc = True)
        else:
            write_netcdf(fp,
                             outfile,
                             temperature=fp["temp"].values.squeeze(),
                             pressure=fp["press"].values.squeeze(),
                             wind_speed=fp["wind"].values.squeeze(),
                             wind_direction=fp["wind_direction"].values.squeeze(),
                             PBLH=fp["PBLH"].values.squeeze(),
                             release_lon=fp["release_lon"].values.squeeze(),
                             release_lat=fp["release_lat"].values.squeeze(),
                             particle_locations = pl,
                             particle_mean_age = pl_ma,
                             particle_heights = height_out,
                             global_attributes = fp.attrs,
                             lapse_rate = lapse_in,
                             lapse_error = lapse_error_in, units = units)
            

    else:
        status_log("FAILED. Couldn't seem to find any files, or some files are missing for %s" %
                   datestr,
                   error_or_warning="error")

    return fp
    


def process_parallel(domain, site, height, years, months, kwargs={}, nprocess=4):
    """
    This script processes multiple months in parallel (but loops through years
    in serial). It basically just calls process.process in parallel.
    
    Args:
        domain (str):
            Domain of interest
        site (str):
            Observation site
        height (str):
            Height of observation, e.g. "10magl"
        years (list):
            List of years you want to process 
        month (int):
            List of months to process (e.g. np.arange(12)+1)
        kwargs (dict):
            Dictionary of optional arguments that you would pass to process, e.g.
            kwargs = {"met_folder":"Met_daily", "fields_folder" : "MixR_files"}
        nprocess (int):
            Number of threads (and number of months to run in parallel).
            Default = 4
            
    Returns:
        None
    """
    nprocess = int(nprocess)
    for year in years:
        pool = Pool(processes = nprocess)
        try:
            [pool.apply_async(process, args=(domain, site, height, year, month), kwds=kwargs) \
             for month in months]
        except:
            pool.close()
            pool.join()
            raise
        pool.close()
        pool.join()

def copy_processed(domain):
    '''
    Routine to copy files from:
        /dagage2/agage/metoffice/NAME_output/DOMAIN_SITE_HEIGHT/Processed_Fields_files
        to:
        air.chm:/data/shared/LPDM/fp_netcdf/DOMAIN/
        
    Args:
        domain (str):
            Domain of interest
            
    Returns:
        None
    '''
    import dirsync
    
    src_folder = "/dagage2/agage/metoffice/NAME_output/"
    dst_folder = "/data/shared/LPDM/fp_NAME/" + domain + "/"
    
    files = glob.glob(src_folder + domain +
        "_*/Processed_Fields_files/*.nc")
        
    folders = set([os.path.split(f)[0] for f in files])

    for f in folders:
        print("Syncing " + f + " and " + dst_folder)
        dirsync.sync(f, dst_folder, "sync")
    print("Done sync")

def test_processed_met(domain, site, height,
                       base_dir = "/dagage2/agage/metoffice/NAME_output/"):
    """Test to check that met is OK.
    
    Plots pressure, temperature, footprints and boundary particle locations.
    Outputs dataset of met.
    
    Args:
        domain (str):
            Domain of interest
        site (str):
            Observation site
        height (str):
            Height of observation, e.g. "10magl"
        base_dir (str, optional):
            Base directory containing NAME output
            Default="/dagage2/agage/metoffice/NAME_output/".
    
    Returns:
        xarray dataset containing met data
            
    """
    
    subfolder = base_dir + domain + "_" + site + "_" + height + \
                "/Processed_Fields_files/"
    files = sorted(glob.glob(subfolder + "*.nc"))

    ds = []    
    for f in files:
        with xray.open_dataset(f) as fi: 
            dsf = fi.load()

        ds.append(dsf)
    
    ds = xray.concat(ds, "time")
    
    plt.plot(ds.time, ds.pressure)
    plt.title("Pressure")
    plt.show()

    plt.plot(ds.time, ds.temperature)
    plt.title("Temperature")
    plt.show()

    plt.plot(ds.time, ds.fp.sum(["lon", "lat"]))
    plt.title("Sum(Footprint)")
    plt.show()

    plt.plot(ds.time, ds.particle_locations_n.sum(["lon", "height"]) + \
                      ds.particle_locations_e.sum(["lat", "height"]) + \
                      ds.particle_locations_s.sum(["lon", "height"]) + \
                      ds.particle_locations_w.sum(["lat", "height"]))
    plt.title("Sum(Particle locations) - should sum to 1")
    plt.show()
    
    return ds

# # Process a list of AGAGE/DECC/GAUGE files if called as main
# if __name__ == "__main__":

#     domain = "EUROPE"
    
#     sites = ["BSD", "TTA", "RGL", "MHD", "HFD", "TAC",
#              "GAUGE-FERRY", "GAUGE-FAAM",
#              "EHL", "TIL", "GLA", "WAO", "HAD"]
#     for site in sites:
#         process_all(domain, site, force_update = True)
        
def process_vertical_profile(vp_fname):
    """
    Function to process the site specific vertical pressure and temperature gradients.
    Relies on the process_met function
    In order to fit into this function the header information needs to be edited
    in the output files NAME generates. 
    Required structure given in vp_met_dict.
    
    Args:
        vp_name (str, optional): 
            Name of vertical profiles file name
    
    Returns: 
        xarray dataset containing:
         - theta_slope - the potential temperature gradient at each time point
         - slope_error - the error in this calculated gradient
    """ 
    
    vp_met_dict = {"time": "T","temp20": "TEMP-Z = 20.00000 m agl",
               "temp40": "TEMP-Z = 40.00000 m agl","temp60": "TEMP-Z = 60.00000 m agl",
               "temp80": "TEMP-Z = 80.00000 m agl","temp100": "TEMP-Z = 100.0000 m agl",
               "temp120": "TEMP-Z = 120.0000 m agl","temp140": "TEMP-Z = 140.0000 m agl",
               "temp160": "TEMP-Z = 160.0000 m agl","temp180": "TEMP-Z = 180.0000 m agl",
               "temp200": "TEMP-Z = 200.0000 m agl","temp220": "TEMP-Z = 220.0000 m agl",
               "temp240": "TEMP-Z = 240.0000 m agl","temp260": "TEMP-Z = 260.0000 m agl",
               "temp280": "TEMP-Z = 280.0000 m agl","temp300": "TEMP-Z = 300.0000 m agl",
               "press20": "PRES-Z = 20.00000 m agl","press40": "PRES-Z = 40.00000 m agl",
               "press60": "PRES-Z = 60.00000 m agl","press80": "PRES-Z = 80.00000 m agl",
               "press100": "PRES-Z = 100.0000 m agl","press120": "PRES-Z = 120.0000 m agl",
               "press140": "PRES-Z = 140.0000 m agl","press160": "PRES-Z = 160.0000 m agl",
               "press180": "PRES-Z = 180.0000 m agl","press200": "PRES-Z = 200.0000 m agl",
               "press220": "PRES-Z = 220.0000 m agl","press240": "PRES-Z = 240.0000 m agl",
               "press260": "PRES-Z = 260.0000 m agl","press280": "PRES-Z = 280.0000 m agl",
               "press300": "PRES-Z = 300.0000 m agl",
               "theta20": "THETA-Z = 20.00000 m agl","theta40": "THETA-Z = 40.00000 m agl",
               "theta60": "THETA-Z = 60.00000 m agl","theta80": "THETA-Z = 80.00000 m agl",
               "theta100": "THETA-Z = 100.0000 m agl","theta120": "THETA-Z = 120.0000 m agl",
               "theta140": "THETA-Z = 140.0000 m agl","theta160": "THETA-Z = 160.0000 m agl",
               "theta180": "THETA-Z = 180.0000 m agl","theta200": "THETA-Z = 200.0000 m agl",
               "theta220": "THETA-Z = 220.0000 m agl","theta240": "THETA-Z = 240.0000 m agl",
               "theta260": "THETA-Z = 260.0000 m agl","theta280": "THETA-Z = 280.0000 m agl",
               "theta300": "THETA-Z = 300.0000 m agl"}
    


    vp_met_df=read_met(vp_fname, met_def_dict=vp_met_dict, vertical_profile=True)
    lapse_ds=xray.Dataset.from_dataframe(vp_met_df)
    
    v_all = np.arange(60,320,20)
    x_all = np.asarray([
     lapse_ds.theta60.values,lapse_ds.theta80.values,lapse_ds.theta100.values,
     lapse_ds.theta120.values,lapse_ds.theta140.values,lapse_ds.theta160.values,
     lapse_ds.theta180.values,lapse_ds.theta200.values,lapse_ds.theta220.values,
     lapse_ds.theta240.values,lapse_ds.theta260.values,lapse_ds.theta280.values,
     lapse_ds.theta300.values])

    slope_all=np.zeros((len(lapse_ds.time)))
    std_all=np.zeros((len(lapse_ds.time)))
    for ti in range(len(lapse_ds.time)):                                                              
        slope_all[ti], dum, dum2, dum3, std_all[ti] = scipy.stats.linregress(
                                                    v_all,x_all[:,ti])
    
    lapse_time=lapse_ds.time
            
    lapse_ds2 = xray.Dataset({'theta_slope': (['time'], slope_all*1000.),
                              'slope_error': (['time'], std_all*1000.)},
                                    coords = {'time' : (lapse_time)})
    #ost_mcmc.coords["sites"]=sites
    lapse_ds2.theta_slope.attrs['units'] = "K/km"
    lapse_ds2.slope_error.attrs['units'] = "K/km"    
    lapse_ds2.theta_slope.attrs['comment'] = "Potential temperature gradient from 60 - 300 m"
    
    return lapse_ds2


def release_point(nc):
    """
    Extract the time and location of release from STILT output based on the
    identifier string. 
    
    Args:
        nc (netCDF4.Dataset): 
            netCDF4.Dataset in STILT format. 
            
    Returns:
        release_lat (float):
            Release latitude.
        release_lon (float):
            Release longitude.
        release_ht (float):
            Release height (units??).
        time (float):
            Release time. 
        time_seconds (float):
            ???
        time_reference (float):
            ???
            
    Todo:
        Update outputs with units and what they actually are.
    """
    ident = netCDF4.chartostring(nc.variables['ident'][:])[0]
    ident = ident.split('x')
    
    if len(ident)==8:
        ident = ident[0:4] + ident[5:8]

    if ident[4][-1] == 'N':
        release_lat = float(ident[4][:-1])
    elif ident[4][-1] == 'S':
        release_lat = -float(ident[4][:-1])
    
    if ident[5][-1] == 'E':
        release_lon = float(ident[5][:-1])
    elif ident[5][-1] == 'W':
        release_lon = -float(ident[5][:-1])
        
    release_ht = float(ident[6])
        
    time = "-".join(ident[0:3]) + " " + ident[3] + ":00"
    time = [acrg_time.convert.dateutil.parser.parse(time)]
    time_seconds, time_reference = time2sec(time)
    
    return release_lat, release_lon, release_ht, time, time_seconds, time_reference

def stilt_part(nc, Id, domain_N, domain_S, domain_E, domain_W, exitfile, append):
    """
    Extract particle data from STILT output and assign it the given Id.
    
    Args:
        nc (netCDF4.Dataset): 
            netCDF4.Dataset in STILT format. 
        ID (int):
            Identifier ID
        domain_N (float):
            N latitude 
        domain_S (float):
            S latitude
        domain_E (float):
            E longitude 
        domain_W (float):
            W longitude
        exitfile (str):
            Name of out file
        append (bool):
            Append header True/False
        
        Returns:
            part (dataframe):
                Particle data
            partnames (list):
                Variable names
            
    """
    part = pd.DataFrame(nc.variables['part'][:].T)
    partnames = list(netCDF4.chartostring(nc.variables['partnames'][:]))
    part.columns = partnames
    part['Id'] = Id
            
    # clean up particle dataframe
    
    part=part.rename(columns = {'index':'partno'})
    part['partno'] = part['partno'].astype('category')
    
    # set up domain boundaries
    
    north = (part['lat']>domain_N)
    south = (part['lat']<domain_S)
    east = (part['lon']>domain_E) 
    west = (part['lon']<domain_W)
    out = north | south | east | west
    part['out'] = out
    
    # compute domain exit timestep for particles that exit
    outpart = part.loc[out,:]
    exittime = outpart.groupby(['Id', 'partno']).apply(lambda df:df.time.argmin())
    if exittime.empty:
        exittime = pd.Series(None) #otherwise it will become a DataFrame
    
    # compute last timestep for particles that do not exit
    gp = part.groupby(['Id', 'partno'])
    lasttime = gp.apply(lambda df:df.time.argmin())
    leaves = gp.agg({'out' : any})
    lasttime = lasttime[~leaves['out']]
    if lasttime.empty:
        lasttime = pd.Series(None) #otherwise it will become a DataFrame
    
    # combine particles that do and do not exit
    
    end = pd.concat([exittime, lasttime])
    part = part.iloc[end.tolist()]
    part = part.loc[:,['Id','lat','lon','agl']]
    part.columns = ['Id','Lat','Long','Ht']
    
    # write particle locations to csv
    if append:
        header = False
        mode='a'
    else:
        header=True
        mode='w'
    part.to_csv(exitfile, index=False, sep=" ", header=header, mode=mode)

    return part, partnames

def stiltfoot_array(prefix, 
                    exitfile=None,
                    met=None,
                    satellite=False,
                    time_step=None):
    """
    Read data from STILT netcdf output files into arrays in an xray dataset.
    The corresponding function for NAME output is footprint_array.
    This function reads multiple files at once because STILT output is stored
    as a separate file for each release time. To read a month of output, give
    a prefix like stilt2016x07x.
    
    Args:
        prefix (str): 
            File pattern prefix (including path) for output files to read
        exitfile (str, optional): 
            File in which to temporarily store particle data.
            Default = None
        met (dict): 
            not used to interpret STILT footprints, which are already given
            in ppm/(mol/m^2/s), but can be included if available.
            Default = None.
        satellite (bool, optional): 
            Not implemented; will produce error if True.
            Default=False
        time_step (float, optional): 
            Ignored; not needed to interpret STILT footprints.
            Default=None.
            
    Returns:
        fp (xarray datset): 
            An xarray dataset as from footprint_concatenate 
            (see help(process.footprint_concatenate))
    """
    if exitfile is None:
        exitfile = prefix + "stilt_particle_data_temp.csv"

    if satellite:
        status_log("stiltfoot_array is not set up for satellite data!" +\
                   " Levels will be wrong!", error_or_warning="error")
    
    stiltfiles = glob.glob(prefix + "*.nc")
    stiltfiles.sort()
    patternfile = stiltfiles[0]
    
    # set up lat/lon grid and release location
    
    ncin = netCDF4.Dataset(patternfile)
    
    lats = ncin.variables['footlat'][:]
    lons = ncin.variables['footlon'][:]
    nlon=len(lons)
    nlat=len(lats)
    
    domain_N = max(lats)
    domain_S = min(lats)
    domain_E = max(lons)
    domain_W = min(lons)
    
    release_lat, release_lon, _, time, _, time_reference = release_point(ncin)
    
    levs=["0"]
    nlev=len(levs)
    
    # set up footprint array and particle dataframe
    
    fparrays = [np.sum(ncin.variables['foot'][:], axis=0)]
    
    part, partnames = stilt_part(ncin, 1, domain_N, domain_S, domain_E, domain_W,
                                 exitfile, append=False)
    
    # extract data from remaining files
    
    n=1
    for f in stiltfiles[1:]:
        status_log("Loading " + f)
        ncin = netCDF4.Dataset(f)
        
        this_lat, this_lon, _, this_time, _, this_reference = release_point(ncin)
        if this_reference != time_reference:
            status_log("Warning! Inconsistent time reference in " + f, \
                               error_or_warning="warning")
        
        if (this_lat != release_lat) or (this_lon != release_lon):
            status_log("Warning! Inconsistent release location in " + f + \
                               ". Skipping.", error_or_warning="warning")
            continue

        this_grid_lats = ncin.variables.get('footlat', [])[:]
        this_grid_lons = ncin.variables.get('footlon', [])[:]
        if (len(this_grid_lats)==0) or (len(this_grid_lons)==0):
            status_log("Warning! Footprint grid missing in " + f + \
                       ". Skipping.", error_or_warning="warning")
            continue
        
        if (len(this_grid_lats)!=len(lats) or len(this_grid_lons)!=len(lons)):
            status_log("Warning! Inconsistent domain in " + f + \
                               ". Skipping.", error_or_warning="warning")
            continue
        
        if (any(lats!=this_grid_lats) or any(lons!=this_grid_lons)):
            status_log("Warning! Inconsistent domain in " + f + \
                               ". Skipping.", error_or_warning="warning")
            continue
        
        p, pnames = stilt_part(ncin, n+1, domain_N, domain_S, domain_E, domain_W,
                               exitfile, append=True)
        if pnames!=partnames:
            status_log("Warning! Inconsistent particle variables in " + \
                               f + ".", error_or_warning="warning")
            
        time.extend(this_time)
        n=n+1
#        part = part.append(p, ignore_index=True)
        fparrays.append(np.sum(ncin.variables['foot'][:], axis=0))
    
    fparrays = [x/1000000 for x in fparrays] # convert ppm to mol/mol 
    
    # compute particle domain exit statistics
    
    dheights = 1000
    heights = np.arange(0, 19001, dheights) + dheights/2.
    
    particle_hist = particle_locations(exitfile,
                                       time, lats, lons, levs,
                                       heights, id_is_lev = False)
    
    # set up footprint dataset 
    ntime=len(time)
    fp = xray.Dataset({"fp": (["time", "lev", "lat", "lon"],
                              np.zeros((ntime, nlev, nlat, nlon)))},
                        coords={"lat": lats, "lon":lons, "lev": levs, 
                                "time":time})
    
    fp.fp.attrs = {"units": "(mol/mol)/(mol/m2/s)"}
    fp.lon.attrs = {"units": "degrees_east"}
    fp.lat.attrs = {"units": "degrees_north"}
    fp.time.attrs = {"long_name": "start of time period"}
    
    # assign dummy met variables 
    # This met information is not from observations and is not used to interpret 
    # STILT output, which is already given in ppm. It is here only to maintain the
    # expected format.
    if met is None:
        met = met_empty()
    if type(met) is not list:
        met = [met]
    
    met_dict = {key : (["time", "lev"], np.zeros((len(time), len(levs))))
                for key in list(met[0].keys())}
    met_ds = xray.Dataset(met_dict,
                          coords = {"time": (["time"], time),
                                    "lev": (["lev"], levs)})
    levi = 0
    levi_met = 0
    
    metr = met[levi_met].reindex(index = time, method = "pad")
    
    metr['release_lat'] = release_lat
    metr['release_lon'] = release_lon
    
    for key in list(met[levi_met].keys()):
        met_ds[key][dict(lev = [levi])] = \
        metr[key].values.reshape((len(time), 1))
            
    # merge datasets
        
    fp = fp.merge(met_ds)
    fp = fp.merge(particle_hist)
    
    # insert footprint arrays into dataset
    
    for i in range(len(time)):
        slice_dict = dict(time = [i], lev = [0])
        fp.fp[slice_dict] = fparrays[i]
        
    return fp

def open_obs(satellite_obs_file=None,subfolder=None,datestr=None,obs_folder="Observations"):
    '''
    Opens the satellite observations file. This should match
    the points within the NAME input csv and with the output of NAME.
    Used for satellite data as NAME input is not precise enough (minute precision).
    
    Can either specify satellite_obs_file or specify subfolder,datestr,obs_folder so obs_fils can be found. 
    
    Args:
        satellite_obs_file (str) :
            Observation filename.
        subfolder (str) :
            Only include if satellite_obs_file is not specified.
            Folder containing the output from NAME. Folder containing observations should be
            within this directory.
        datestr (str) :
            Only include if satellite_obs_file is not specified.
            The date string being used to search for files. This should be of the form
            'YYYYMMDD' or 'YYYYMMDD-NNN'
        obs_folder (str, optional) :
            Only include if satellite_obs_file is not specified.
            Name of the folder containing the observations data. Will be appended to subfolder
            to get the full path information.
    
     Returns:
        xarray.dataset:
            satellite file
            
    '''
    if not satellite_obs_file:
        try:
            search_str = os.path.join(subfolder,"Observations/*{}_*.nc".format(datestr))
        except AttributeError:
            status_log("There must be exactly one matching satellite " + 
                           "file in the Observations/ folder",
                           error_or_warning="error")
            status_log("Files: " + ','.join(satellite_obs_file),
                           error_or_warning="error")
            return None
        
        satellite_obs_files = glob.glob(search_str)

        if len(satellite_obs_files) > 1:
            status_log("There must be exactly one matching satellite " + 
                               "file in the Observations/ folder",
                               error_or_warning="error")
            status_log("Files: " + ','.join(satellite_obs_files),
                               error_or_warning="error")
        
        satellite_obs_file = satellite_obs_files[0]
    
    with xray.open_dataset(satellite_obs_file) as f:
        sat = f.load()
    
    return sat

def extract_time_obs(satellite_obs_file=None,subfolder=None,datestr=None,obs_folder="Observations"):
    '''
    '''
    if satellite_obs_file:
        sat = open_obs(satellite_obs_file)
    else:
        sat = open_obs(subfolder=subfolder,datestr=datestr,obs_folder=obs_folder)
    
    return sat.time.values
    
