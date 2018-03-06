# -*- coding: utf-8 -*-
"""
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

import numpy
from netCDF4 import Dataset
import glob
import gzip
import datetime
import re
import pandas
import datetime as dt
import numpy as np
import scipy.constants as const
from acrg_grid import areagrid
from acrg_time.convert import time2sec, sec2time
import os
import json
from os.path import split, realpath, exists
import xarray as xray
import shutil
from scipy.interpolate import interp1d
import copy
import dirsync
import matplotlib.pyplot as plt
import getpass
import traceback
import sys
import scipy

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

def load_NAME(file_lines, namever):
    """
    Loads the Met Office's NAME III grid output files returning
    headers, column definitions and data arrays as 3 separate lists.
    
    For a technical specification see 
    http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf
    
    Inputs to this routine are a list of file text lines (extracted using
    the extract_file_lines routine), and the NAME version.

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
    #MLR
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
        
    # skip the blank line after the column headers
#    file_handle.next()
    #MLR
    file_lines=file_lines[1:]
    
    # make a list of data arrays to hold the data for each column 
    data_shape = (headers['Y grid size'], headers['X grid size'])
    if namever==2:
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of fields'])]
    elif namever==3:
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of field cols'])]
      
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

    return headers, column_headings, data_arrays
    

def extract_file_lines(fname):
    '''
    Determine whether file is zipped or not, and then extract text into
    file_lines variable.
    '''
    
    if fname[-3:].upper() == ".GZ":
        #Unzip file
        f=gzip.open(fname, 'r')
        file_text=f.read()
        f.close()
        file_lines=file_text.splitlines()
    else:
        file_lines=[i.strip() for i in open(fname).readlines()]

    return file_lines


def read_file(fname):
    '''
    Read a given fields file and break into header, column_headings and data_arrays
    '''
    
    global namever
    
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

    return header, column_headings, data_arrays


def define_grid(header, column_headings, satellite = False):
    '''
    Largely taken from Met Office Python code.
    
    Define output grid using file header information.
    
    For satellite instructions, see satellite_vertical_profile.
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
        lons=numpy.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep) + Xstep/2.
        lats=numpy.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep) + Ystep/2.
    elif namever == 3:
        lons=numpy.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep)
        lats=numpy.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep)
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
        ntime=len(time)/nlevs
        
    else:
        # Assume only ONE TIME STEP per file
        #Levels
        nlevs = len(z_level)
        levs = range(nlevs)
        ntime = 1


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
    
    met = pandas.DataFrame({key: 0. for key in met_default.keys() if key != "time"},
                          index = [dt.datetime(1900, 1, 1), dt.datetime(2020, 1, 1)])
    met.index.name = "time"
    met["press"] = [100000., 100000.] #Pa
    met["temp"] = [10., 10.]    #C
    
    status_log("NO MET", error_or_warning="warning")

    return met


def read_met(fnames, met_def_dict=None, vertical_profile=False):
    '''
    For a given list of filenames, extract site meteorology and concatenate
    into a Pandas dataframe.
    
    A dictionary of output column names and met file header search strings
    is given in the met_default dictionary at the top of this file.
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
    column_indices = {key: -1 for key, value in met_default2.iteritems()}
    
    output_df = []
        
    for fname in fnames:
        
        #This should now check for column headings
        #Mark's met file output is daily and backwards
        met = extract_file_lines(fname)
        if fname[-3:].upper() == '.GZ':
            compression="gzip"
        else:
            compression=None

        #Get rid of top of file
        for i in range(len(met)):
            cols = [col.strip() for col in met[i].split(',')]
            if len(cols)>1:
                a=i
                break

        m = pandas.read_csv(fname, skiprows = a, 
                            compression=compression)
        m = m.fillna('')
        
        m = numpy.array(m)

        Xcol = None
        Ycol = None
        X_file = None
        Y_file = None
        
        #Find columns with header names
        for row in m:
            for coli, col in enumerate(row):
                if type(col) is str:
                    
                    #Work out column indices by searching through met_default
                    for key in met_default2.keys():
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
        for key in met_default2.keys():
            if column_indices[key] != -1 and key != "time":
                met_dict[key] = m2[:, column_indices[key]].astype(float)
        met_dict["release_lon"] = X
        met_dict["release_lat"] = Y

        #Construct dataframe
        output_df_file = pandas.DataFrame(met_dict,
                          index=[pandas.to_datetime(d, dayfirst=True)
                              for d in m2[:,column_indices["time"]]])

        output_df.append(output_df_file)
    
        status_log("Read Met file... " + os.path.split(fname)[1])
    
    # Concatenate list of data frames
    output_df = pandas.concat(output_df)
    
    # Check for missing values
    if vertical_profile == True:
        output_df = output_df[output_df["press20"] > 0.]
    else:
        output_df = output_df[output_df["press"] > 0.]
    output_df.drop_duplicates(inplace = True)
    
    # Remove duplicate indices (if found)
    output_df = output_df.groupby(level = 0).last()

    # Sort the dataframe by time
    output_df.sort_index(inplace = True)
    
    output_df.index.name = "time"
    
    return output_df


def particle_locations(particle_file, time, lats, lons, levs, heights,
                       id_is_lev = False):
    '''
    Read and process a particle location file.
    '''

    def particle_location_edges(xvalues, yvalues, x, y):
        
        dx = x[1] - x[0]
        xedges = numpy.append(x - dx/2., x[-1] + dx/2.)
        dy = y[1] - y[0]
        yedges = numpy.append(y - dy/2., y[-1] + dy/2.)
        
        hist, xe, ye = numpy.histogram2d(xvalues, yvalues,
                                         bins = (xedges, yedges))
        
        return hist
    
    edge_lons = [min(lons), max(lons)]
    edge_lats = [min(lats), max(lats)]
    dlons = lons[1] - lons[0]
    dlats = lats[1] - lats[0]
    
    hist = []
    particles = []
    
    hist = xray.Dataset({"pl_n":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "pl_e":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights)))),
                         "pl_s":(["time", "lev", "lon", "height"],
                         np.zeros((len(time), len(levs), len(lons), len(heights)))),
                         "pl_w":(["time", "lev", "lat", "height"],
                         np.zeros((len(time), len(levs), len(lats), len(heights))))},
                        coords={"lat": lats, "lon":lons, "lev":levs, "height":heights, "time":time})

    #Variables to check domain extents
    particle_extremes = {"N": -90., "E": -360. ,"S": 90.,"W": 360.}
    
    status_log("Particle locations " + particle_file, print_to_screen=False)
    
    if particle_file[-3:].upper() == '.GZ':
        compression="gzip"
    else:
        compression=None
    
    df = pandas.read_csv(particle_file, compression=compression, sep=r"\s+")
    
    particles_record = []
    
    for i in set(numpy.array(df["Id"])):
        
        if id_is_lev:
            # Only one time step allowed at the moment
            slice_dict = {"time": [0], "lev": [i-1]}
        else:
            slice_dict = {"time": [i-1]}
            
        #Northern edge
        dfe = df[(df["Lat"] > edge_lats[1] - dlats/2.) & (df["Id"] == i)]
        hist.pl_n[slice_dict] = \
            particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                    lons, heights)

        #Eastern edge
        dfe = df[(df["Long"] > edge_lons[1] - dlons/2.) & (df["Id"] == i)]
        hist.pl_e[slice_dict] = \
            particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                    lats, heights)

        #Southern edge
        dfe = df[(df["Lat"] < edge_lats[0] + dlats/2.) & (df["Id"] == i)]
        hist.pl_s[slice_dict] = \
            particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                    lons, heights)

        #Western edge
        dfe = df[(df["Long"] < edge_lons[0] + dlons/2.) & (df["Id"] == i)]
        hist.pl_w[slice_dict] = \
            particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                    lats, heights)

        #Calculate total particles and normalise
        hist_sum = hist[slice_dict].sum()
        particles = int(sum([hist_sum[key].values for key in hist_sum.keys()]))
        particles_record.append(str(particles))

#        print("Number of particles reaching edge: %f02" %particles)

        if particles > 0.:
            for key in hist.data_vars.keys():
                hist[key][slice_dict] = hist[key][slice_dict]/\
                                                particles
        else:
            status_log("No particles have reached edge",
                       error_or_warning="warning")
            if i > 1:
                status_log("Copying lower level/previous time step",
                           error_or_warning="warning")
                if id_is_lev:
                    slice_dict_prev = {"time": [0], "lev": [i-2]}
                else:
                    slice_dict_prev = {"time": [i-2]}
                for key in hist.data_vars.keys():
                    hist[key][slice_dict] = hist[key][slice_dict_prev]
        
        # Store extremes
        if max(df["Lat"]) > particle_extremes["N"]:
            particle_extremes["N"] = max(df["Lat"])
        if min(df["Lat"]) < particle_extremes["S"]:
            particle_extremes["S"] = min(df["Lat"])
        if max(df["Long"]) > particle_extremes["E"]:
            particle_extremes["E"] = max(df["Long"])
        if min(df["Long"]) < particle_extremes["W"]:
            particle_extremes["W"] = min(df["Long"])

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
                    time_step = None):
    '''
    Convert text output from a given file files into arrays in an 
    xray dataset.
    '''

    global timestep_for_output

    status_log("Reading... " + os.path.split(fields_file)[1])
    
    header, column_headings, data_arrays = read_file(fields_file)

    if met is None:
        met = met_empty()

    if type(met) is not list:
        met = [met]

    # Define grid, including output heights    
    lons, lats, levs, time, timeStep = define_grid(header, column_headings,
                                                   satellite = satellite)
                                 
    # If time_step is input, overwrite value from NAME output file
    if time_step is not None:
        timeStep = time_step
    
    timestep_for_output = timeStep

    dheights = 1000
    heights = numpy.arange(0, 19001, dheights) + dheights/2.

    # Get area of each grid cell
    area=areagrid(lats, lons)
    
    # Get particle locations
    if particle_file is not None:
        particle_hist = particle_locations(particle_file,
                                           time, lats, lons, levs,
                                           heights, id_is_lev = satellite)
    else:
        status_log("No particle location file corresponding to " + fields_file,
                   error_or_warning="error")
    
    nlon=len(lons)
    nlat=len(lats)
    nlev=len(levs)
    ntime=len(time)
    status_log("Time steps in file: %d" % ntime, print_to_screen = False)
    
    z_level=column_headings['z_level'][4:]
    time_column=column_headings['time'][4:]
    
    # Set up footprint dataset
    fp = xray.Dataset({"fp": (["time", "lev", "lat", "lon"],
                              numpy.zeros((ntime, nlev, nlat, nlon)))},
                        coords={"lat": lats, "lon":lons, "lev": levs, 
                                "time":time})

    fp.fp.attrs = {"units": "(mol/mol)/(mol/m2/s)"}
    fp.lon.attrs = {"units": "degrees_east"}
    fp.lat.attrs = {"units": "degrees_north"}
    fp.time.attrs = {"long_name": "start of time period"}
    
    
    # Add in met data
    met_dict = {key : (["time", "lev"], np.zeros((len(time), len(levs))))
                for key in met[0].keys()}
    met_ds = xray.Dataset(met_dict,
                          coords = {"time": (["time"], time),
                                    "lev": (["lev"], levs)})
    
    for levi, lev in enumerate(levs):
        # If there is only one level in the met list
        # but more levels in the NAME output, just use the same met level
        # for all NAME levels. This is probably a bad idea.
        if len(met) == 1:
            levi_met = 0
            if len(levs) > 1:
                status_log("ONLY ONE MET LEVEL. " + \
                           "ASSUMING MET CONSTANT WITH HEIGHT!",
                           error_or_warning="error")
        else:
            levi_met = levi

        metr = met[levi_met].reindex(index = time, method = "pad")

        for key in met[levi_met].keys():
            met_ds[key][dict(lev = [levi])] = \
                metr[key].values.reshape((len(time), 1))

    # Merge met dataset into footprint dataset
    fp.merge(met_ds, inplace = True)


    # Add in particle locations
    if particle_file is not None:
        fp.merge(particle_hist, inplace = True)

    # Extract footprint from columns
    def convert_to_ppt(fp, slice_dict, column):
        molm3=fp["press"][slice_dict].values/const.R/\
            const.C2K(fp["temp"][slice_dict].values.squeeze())
        fp.fp[slice_dict] = data_arrays[i]*area/ \
            (3600.*timeStep*1.)/molm3
        #The 1 assumes 1g/s emissions rate
        return fp

    if satellite:
        for i in range(len(levs)):
            slice_dict = dict(time = [0], lev = [i])
            fp = convert_to_ppt(fp, slice_dict, i)
    else:
        for i in range(len(time)):
            slice_dict = dict(time = [i], lev = [0])
            fp = convert_to_ppt(fp, slice_dict, i)
        
    return fp
    
    
def footprint_concatenate(fields_prefix, 
                          particle_prefix = None,
                          datestr = "*",
                          met = None,
                          satellite = False, time_step = None):
    '''
    For a given file search string (fields_prefix), find all fields and particle
    files, read them and concatenate the output arrays.
    
    Returns an xray dataset.
    
    Example:
    
    fp_dataset = footprint_concatenate("/dagage2/agage/metoffice/NAME_output/MY_FOOTPRINTS_FOLDER/Fields_Files/filename_prefix")
    '''
    

    # If no meteorology, get null values
    # THIS IS NOT RECOMMENDED. ERRORS of ~ ±10%.
    if met is None:
        met = met_empty()

    # Find footprint files and MATCHING particle location files
    # These files are identified by their date string. Make sure this is right!
    #fields_files = sorted(glob.glob(fields_prefix + "*" +
    #                         datestr + "*.txt*"))
    # Modified: 06/03/2018 - problems when # files > 100, point 10 was matching multiple
    fields_files = sorted(glob.glob(fields_prefix + "*" +
                             datestr + ".txt*"))
                             

    # Search for particle files                             
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
        for fields_file, particle_file in \
            zip(fields_files, particle_files):
                fp.append(footprint_array(fields_file,
                                          particle_file = particle_file,
                                          met = met,
                                          satellite = satellite,
                                          time_step = time_step))                                  
    # Concatenate
    if len(fp) > 0:
        fp = xray.concat(fp, "time")
    else:
        fp = None

    return fp


def write_netcdf(fp, lons, lats, levs, time, outfile,
            temperature=None, pressure=None,
            wind_speed=None, wind_direction=None,
            PBLH=None, varname="fp",
            release_lon = None, release_lat = None,
            particle_locations=None, particle_heights=None,
            global_attributes = {}, lapse_rate=None, lapse_error=None):
    '''
    This routine writes a netCDF file with footprints, particle locations
    meteorology, and release locations.
    
    NOTE that netCDF4 is required, as we make use of compression.    
    '''
    
    time_seconds, time_reference = time2sec(time)
    
    #Write NetCDF file
    ncF=Dataset(outfile, 'w')
    ncF.createDimension('time', len(time))
    ncF.createDimension('lon', len(lons))
    ncF.createDimension('lat', len(lats))
    ncF.createDimension('lev', 1)
    
    # pass any global attributes in fp to the netcdf file
    for key in global_attributes.keys():
        ncF.__setattr__(key,global_attributes[key])
    ncF.__setattr__("author", getpass.getuser())
    ncF.__setattr__("created", str(dt.datetime.now()))
    
    
    nctime=ncF.createVariable('time', 'd', ('time',))
    nclon=ncF.createVariable('lon', 'f', ('lon',))
    nclat=ncF.createVariable('lat', 'f', ('lat',))
    nclev=ncF.createVariable('lev', str, ('lev',))
    ncfp=ncF.createVariable(varname, 'f', ('lat', 'lon', 'time'), zlib = True,
                            least_significant_digit = 5)
    
    nctime[:]=time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + numpy.str(time_reference)
    nctime.label='left'
    nctime.period = str(timestep_for_output) + " hours"
    nctime.comment = 'time stamp corresponds to the beginning of each averaging period'
    nctime.calendar='gregorian'

    nclon[:]=lons
    nclon.units='Degrees_east'
    
    nclat[:]=lats
    nclat.units='Degrees_north'

    nclev[:]=numpy.array(levs)
    
    ncfp[:, :, :]=fp
    ncfp.units='(mol/mol)/(mol/m2/s)'

    if temperature is not None:
        nctemp=ncF.createVariable('temperature', 'f', ('time',), zlib = True,
                            least_significant_digit = 4)
        nctemp[:]=temperature
        nctemp.units='K'

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
        ncHeight[:]=particle_heights
        ncPartN[:, :, :]=particle_locations["N"]
        ncPartE[:, :, :]=particle_locations["E"]
        ncPartS[:, :, :]=particle_locations["S"]
        ncPartW[:, :, :]=particle_locations["W"]
        ncPartN.units=''
        ncPartN.long_name='Fraction of total particles leaving domain (N side)'
        ncPartE.units=''
        ncPartE.long_name='Fraction of total particles leaving domain (E side)'
        ncPartS.units=''
        ncPartS.long_name='Fraction of total particles leaving domain (S side)'
        ncPartW.units=''
        ncPartW.long_name='Fraction of total particles leaving domain (W side)'
    
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
    '''
    Do weighted average by satellite averaging kernel and
    pressure weight. One time point only. Expects xray.dataset
    with one time point and N vertical levels.
    
    fp = footprint for ONE time point. N levels in lev dimension 
        with footprints and particle locations defined at each level
    satellite_obs_file = NetCDF-format satellite observation file,
        one per time step
    max_level = the maximum vertical level of the retrieval that is included.
                (maximum for GOSAT to include all levels from model is 20). 
                The remaining levels are set equal to the a priori value.
                
    Requirements for running processing satellite footprints
    
    General file naming: Files/folders need to have a date string of the form
        YYYYMMDD-II at the end of the file/folder name,
        where II is a two-digit index corresponding to the
        sequence of footprints throughout the day. This modifies the usual
        naming convention for surface sites, which are only YYYYMMDD.
    Met: The Met/ folder MUST contain a set of subfolders which are labeled
        with the YYYYMMDD-II time stamp. Each of these subfolders should then
        contain nLev Met output files, where nLev are the number of vertical
        levels in the NAME output. This is required so that we know
        what pressure/temperature each NAME footprint for each vertical level
        corresponds to. The files can be named in any way, but the file names
        must be in ascending order when sorted (i.e. blah_01.txt.gz, blah_02.txt.gz).
    Fields files:
        There must be ONE fields file per time stamp, labeled at the end
        of the file string with YYYYMMDD-II (e.g. blah_YYYYMMDD-II.txt.gz).
        The fields file must contain exactly nLev columns, one for each vertical
        level, in ascending order, left to right.
    Particle location files:
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
        
    if np.abs(sat.lon.values[0] - fp.release_lon.values[0,0]) > 1.:
        status_log("Satellite longitude doesn't match footprints",
                   error_or_warning="error")
    if np.abs(sat.lat.values[0] - fp.release_lat.values[0,0]) > 1:
        status_log("Satellite latitude doesn't match footprints",
                   error_or_warning="error")
    if np.abs(sat.time.values[0] - fp.time.values[0]).astype(int) > 60*1e9:
        status_log("Satellite time doesn't match footprints",
                   error_or_warning="error")
    if len(fp.time.values) > 1:
        status_log("satellite comparison only for one time step at the moment",
                   error_or_warning="error")
        return fp
    if len(sat.time.values) > 1:
        status_log("ERROR: satellite comparison only for one time step at the moment",
                   error_or_warning="error")
        return fp

    if not np.allclose((fp.pl_n.sum() + fp.pl_e.sum() + \
                        fp.pl_s.sum() + fp.pl_w.sum()), \
                        len(fp.lev)):
        status_log("Particle histograms dont add up to 1 (or nlev)",
                   error_or_warning="error")
        return None

    # Change timestamp to that from obs file
    #  because NAME output only has 1 minute resolution
    fp = fp.reindex_like(sat.time, method = "nearest")

    # Interpolate pressure levels
    variables = ["fp", "pl_n", "pl_e", "pl_s", "pl_w"]
    out = {}
    lower_levels =  range(0,max_level)
        
    # Weight levels using pressure weight and averaging kernel
    sum_ak_pw = np.sum(sat.pressure_weights.values.squeeze()[lower_levels] * \
                       sat.xch4_averaging_kernel.values.squeeze()[lower_levels])
    sum_particle_count = 0.

    for variable in variables:
        # get dimensions based on level 0
        fp_lev = int(np.abs(fp.press - \
                            sat.pressure_levels[dict(lev = 0)]*100.).argmin()) 
        var = fp[variable][{"lev": fp_lev}].values.squeeze() * \
              sat.pressure_weights.values.squeeze()[0] * \
              sat.xch4_averaging_kernel.values.squeeze()[0]
        out[variable] = var.reshape((1, 1) + var.shape)
        # run levels 1 - max_level (note it is reindexed to 0, therefore use lev+1)
        for lev, press in enumerate(sat.pressure_levels.values.squeeze()[1:max_level]):
                fp_lev = np.abs(fp.press.values.squeeze() - press*100.).argmin()
                var = fp[variable][{"lev": fp_lev}].values.squeeze() * \
                      sat.pressure_weights.values.squeeze()[lev+1] * \
                      sat.xch4_averaging_kernel.values.squeeze()[lev+1]
                out[variable] += var.reshape((1, 1) + var.shape)
        if variable[0:2] == "pl":
            sum_particle_count += out[variable].sum()

    # Check whether particle sum makes sense
    if not np.allclose(sum_particle_count, sum_ak_pw):
        status_log("ERROR: Particle fractions don't match averaging_kernel * " + \
                   "pressure_weight",
                   error_or_warning = "error")
        return None
    
    # Compress dataset to one level and store column totals
    fp = fp[dict(lev = [0])]
    for variable in variables:
        fp[variable].values = out[variable]
        
    fp.attrs["max_level"] = max_level    
    fp["lev"] = numpy.array(["column"])
    
    return fp


def status_log(message,
               directory = None,
               error_or_warning = "status",
               print_to_screen = True):
    '''
    Write a log of an error or a warning to file. Will append to a file 
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
    """
    Basic processing with no meteorology or particle files
    NOT recommended, but helpful for a quick check
    """
    
    fp = footprint_concatenate(fields_folder)
    write_netcdf(fp.fp.values.squeeze(),
                 fp.lon.values.squeeze(), 
                 fp.lat.values.squeeze(), 
                 fp.lev.values.squeeze(), 
                 fp.time.to_pandas().index.to_pydatetime(), 
                 outfile)

def process(domain, site, height, year, month,
            base_dir = "/dagage2/agage/metoffice/NAME_output/",
            fields_folder = "Fields_files",
            particles_folder = "Particle_files",
            met_folder = "Met",
            force_met_empty = False,
            processed_folder = "Processed_Fields_files",
            satellite = False,
            max_level = None,
            force_update = False,
            perturbed_folder = None,
            vertical_profile=False):
    '''
    Process a single month of footprints for a given domain, site, height,
    year, month. If you want to process all files for a given domain + site
    use process_all.
    
    This routine finds all fields files and particle location files which match
    the timestamp in the file name. The date strings are stored in the datestr
    variable, and these strings must match exactly for the fields and the particle
    locations.
    
    At the moment, the entire Met folder is read each time.
    
    Options:
    force_update: By default, any existing netCDF files are NOT overwritten.
        To explicitly over-write a file, set force_update = True
    satellite: Read a "column" of footprints. There are very particular rules
        about how the met, footprints and particle location files are stored
        and named. See extra instructions in satellite_vertical_profile. 
    perturbed_folder: Process a subfolder which has all of the required folders
        (Fields_files, etc). This is for perturbed parameter ensembles for 
        a particular site. E.g. for EUROPE_BSD_110magl/Perturbed/PARAMETERNAME_VALUE
        you'd set:
            perturbed_folder = "Perturbed/PARAMETERNAME_VALUE"
    max_level: specified only for satellite data and indicates the max level to
        process the foorprints. levels above are replaced by the prior profile
    vertical_profile: If set to true will look for vertical potential temperature met file
        and incorporate into footprint file. 
        NB. This is a separate file from the normal met file, and is not mandatory.
        
    Outputs:
    This routine outputs a copy of the xray dataset that is written to file.
    '''
    
    global directory_status_log
        
    subfolder = base_dir + domain + "_" + site + "_" + height + "/"
    
    directory_status_log = subfolder
    
    if perturbed_folder is not None:
        if perturbed_folder[-1] == "/":
            subfolder += perturbed_folder
        else:
            subfolder += perturbed_folder + "/"
    
    # Check for manual timestep (if footprints are for < 1min,
    # which is the min resolution of the NAME output file)    
    if os.path.exists(subfolder + "/time_step.txt"):
        with open(subfolder + "/time_step.txt") as f:
            timeStep = float(f.read())
    else:
        timeStep = None
    
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
        datestrs = [str(year) + str(month).zfill(2)]

    # Output filename
    outfile = subfolder + processed_folder + "/" + site + "-" + height + \
                "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc"
    
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
            if satellite:
                #met_search_str = subfolder + met_folder + "/*" + datestr + "*/*.txt*"
                # Modified: 06/03/2018 - problems when # files > 100, point 10 was matching multiple
                met_search_str = subfolder + met_folder + "/*" + datestr + "/*.txt*"
            else:
                met_search_str = subfolder + met_folder + "/*.txt*"
      
            met_files = sorted(glob.glob(met_search_str))
        
            if len(met_files) == 0:
                status_log("Can't file MET files: " + met_search_str,
                           error_or_warning="error")
                return None
            else:
                if satellite:
                    met = []
                    for met_file in met_files:
                        met.append(read_met(met_file))
                else:
                    met = read_met(met_files)
        else:
            met = None

        # Get footprints
        fields_prefix = subfolder + fields_folder + "/"
        if particles_folder is not None:
            particles_prefix = subfolder + particles_folder + "/"
        else:
            particles_prefix = None

           
        fp_file = footprint_concatenate(fields_prefix,
                                            datestr = datestr, met = met,
                                            particle_prefix = particles_prefix,
                                            satellite = satellite,
                                            time_step = timeStep)
        
        # Do satellite process
        if satellite:
            #satellite_obs_file = glob.glob(subfolder + "Observations/*" + \
            #                               datestr + "*.nc")
            # Modified: 06/03/2018 - problems when # files > 100, point 10 was matching multiple
            satellite_obs_file = glob.glob(subfolder + "Observations/*" + \
                                           datestr + "_" + "*.nc")
            if len(satellite_obs_file) != 1:
                status_log("There must be exactly one matching satellite " + 
                           "file in the Observations/ folder",
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
        
        #Write netCDF file
        #######################################
        
        # Define particle locations dictionary (annoying)
        if "pl_n" in fp.keys():
            pl = {"N": fp.pl_n.transpose("height", "lon", "time").values.squeeze(),
                  "E": fp.pl_e.transpose("height", "lat", "time").values.squeeze(),
                  "S": fp.pl_s.transpose("height", "lon", "time").values.squeeze(),
                  "W": fp.pl_w.transpose("height", "lat", "time").values.squeeze()}
            height_out = fp.height.values.squeeze()
        else:
            pl = None
            height_out = None

        status_log("Writing file: " + outfile, print_to_screen=False)
        
        # Write outputs
        #try:
        write_netcdf(fp.fp.transpose("lat", "lon", "time").values.squeeze(),
                         fp.lon.values.squeeze(),
                         fp.lat.values.squeeze(),
                         fp.lev.values,
                         fp.time.to_pandas().index.to_pydatetime(),
                         outfile,
                         temperature=fp["temp"].values.squeeze(),
                         pressure=fp["press"].values.squeeze(),
                         wind_speed=fp["wind"].values.squeeze(),
                         wind_direction=fp["wind_direction"].values.squeeze(),
                         PBLH=fp["PBLH"].values.squeeze(),
                         release_lon=fp["release_lon"].values.squeeze(),
                         release_lat=fp["release_lat"].values.squeeze(),
                         particle_locations = pl,
                         particle_heights = height_out,
                         global_attributes = fp.attrs,
                         lapse_rate = lapse_in,
                         lapse_error = lapse_error_in)

    else:
        status_log("FAILED. Couldn't seem to find any files, or some files are missing for %s" %
                   datestr,
                   error_or_warning="error")

    return fp
    

def process_all(domain, site,
                heights = None,
                years_in = None,
                months_in = None,
                base_dir = "/dagage2/agage/metoffice/NAME_output/",
                force_update = False,
                satellite = False,
                perturbed_folder = None,
                max_level = None,
                force_met_empty=False,
                vertical_profile=False):
    '''
    For a given domain and site, process all available fields files (including
    multiple heights).
    If you want to specify a subset of years/months to process, use the years_in
    or months_in kewords.
    
    Keywords:
    heights: If you only want to process a subset of heights, OR IF THE HEIGHT
        INFORMATION IS NOT CONTAINED IN acrg_site_info.json, specify a list of
        height strings here.
    years_in: If you only want to process a subset of years, specify here.
    months_in: Ditto for months
    force_update: By default, this routine will not overwrite existing netcdf
        files. Set force_update = True to reprocess an existing file
    satellite: If a column of footprints needs to be processed, set to true.
        See extra instructions in satellite_vertical_profile. 
    perturbed_folder: Process a subfolder which has all of the required folders
        (Fields_files, etc). This is for perturbed parameter ensembles for 
        a particular site. E.g. for EUROPE_BSD_110magl/Perturbed/PARAMETERNAME_VALUE
        you'd set:
            perturbed_folder = "Perturbed/PARAMETERNAME_VALUE"
    '''

    acrg_path = os.getenv("ACRG_PATH")
    with open(os.path.join(acrg_path, "acrg_site_info.json")) as f:
        site_info=json.load(f)
        
    # If no height specified, run all heights
    if heights is None:
        heights = site_info[site]["height_name"]
    elif type(heights) is not list:
        heights = [heights]
    
    for height in heights:
        
        subfolder = base_dir + domain + "_" + site + "_" + height + "/"
        if perturbed_folder is not None:
            if perturbed_folder[-1] == "/":
                subfolder += perturbed_folder
            else:
                subfolder += perturbed_folder + "/"
        
        if years_in is None:
            #Find all years and months available
            #Assumes fields files are processes with _YYYYMMDD.txt.gz at the end    
            years = []
            months = []
            
            fields_files = sorted(glob.glob(subfolder + "/Fields_files/*.txt*"))
            for fields_file in fields_files:
                f = split(fields_file)[1].split("_")[-1].split('.')[0]
                years.append(int(f[0:4]))
                months.append(int(f[4:6]))
        else:
            years = copy.copy(years_in)
            months = copy.copy(months_in)

        for year, month in set(zip(years, months)):
            #try:
            out = process(domain, site, height, year, month,
                    base_dir = base_dir, force_update = force_update,
                    satellite = satellite, perturbed_folder = perturbed_folder,
                    max_level = max_level, force_met_empty = force_met_empty,
                    vertical_profile=vertical_profile)


def copy_processed(domain):
    '''
    Routine to copy files from:
        /dagage2/agage/metoffice/NAME_output/DOMAIN_SITE_HEIGHT/Processed_Fields_files
        to:
        air.chm:/data/shared/NAME/fp_netcdf/DOMAIN/
    '''
    
    src_folder = "/dagage2/agage/metoffice/NAME_output/"
    dst_folder = "/data/shared/NAME/fp/" + domain + "/"
    
    files = glob.glob(src_folder + domain +
        "_*/Processed_Fields_files/*.nc")
        
    folders = set([os.path.split(f)[0] for f in files])

    for f in folders:
        print("Syncing " + f + " and " + dst_folder)
        dirsync.sync(f, dst_folder, "sync")
    print("Done sync")

def test_processed_met(domain, site, height,
                       base_dir = "/dagage2/agage/metoffice/NAME_output/"):
    
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

# Process a list of AGAGE/DECC/GAUGE files if called as main
if __name__ == "__main__":

    domain = "EUROPE"
    
    sites = ["BSD", "TTA", "RGL", "MHD", "HFD", "TAC",
             "GAUGE-FERRY", "GAUGE-FAAM",
             "EHL", "TIL", "GLA", "WAO", "HAD"]
    for site in sites:
        process_all(domain, site, force_update = True)
        
def process_vertical_profile(vp_fname):
    """
    Function to process the site specific vertical pressure and temperature gradients.
    Relies on the process_met function
    In order to fit into this function the header information needs to be edited
    in the output files NAME generates. 
    Required structure given in vp_met_dict
    
    Outputs: xarray dataset containing:
        theta_slope - the potential temperature gradient at each time point
        slope_error - the error in this calculated gradient
    """
    
#    vp_met_dict = {"time": "             T","temp20": "TEMP-Z = 20.00000 m agl",
#               "temp40": "TEMP-Z = 40.00000 m agl","temp60": "TEMP-Z = 60.00000 m agl",
#               "temp80": "TEMP-Z = 80.00000 m agl","temp100": "TEMP-Z = 100.0000 m agl",
#               "temp120": "TEMP-Z = 120.0000 m agl","temp140": "TEMP-Z = 140.0000 m agl",
#               "temp160": "TEMP-Z = 160.0000 m agl","temp180": "TEMP-Z = 180.0000 m agl",
#               "temp200": "TEMP-Z = 200.0000 m agl","temp220": "TEMP-Z = 220.0000 m agl",
#               "temp240": "TEMP-Z = 240.0000 m agl","temp260": "TEMP-Z = 260.0000 m agl",
#               "temp280": "TEMP-Z = 280.0000 m agl","temp300": "TEMP-Z = 300.0000 m agl",
#               "press20": "PRES-Z = 20.00000 m agl","press40": "PRES-Z = 40.00000 m agl",
#               "press60": "PRES-Z = 60.00000 m agl","press80": "PRES-Z = 80.00000 m agl",
#               "press100": "PRES-Z = 100.0000 m agl","press120": "PRES-Z = 120.0000 m agl",
#               "press140": "PRES-Z = 140.0000 m agl","press160": "PRES-Z = 160.0000 m agl",
#               "press180": "PRES-Z = 180.0000 m agl","press200": "PRES-Z = 200.0000 m agl",
#               "press220": "PRES-Z = 220.0000 m agl","press240": "PRES-Z = 240.0000 m agl",
#               "press260": "PRES-Z = 260.0000 m agl","press280": "PRES-Z = 280.0000 m agl",
#               "press300": "PRES-Z = 300.0000 m agl",
#               "theta20": "THETA-Z = 20.00000 m agl","theta40": "THETA-Z = 40.00000 m agl",
#               "theta60": "THETA-Z = 60.00000 m agl","theta80": "THETA-Z = 80.00000 m agl",
#               "theta100": "THETA-Z = 100.0000 m agl","theta120": "THETA-Z = 120.0000 m agl",
#               "theta140": "THETA-Z = 140.0000 m agl","theta160": "THETA-Z = 160.0000 m agl",
#               "theta180": "THETA-Z = 180.0000 m agl","theta200": "THETA-Z = 200.0000 m agl",
#               "theta220": "THETA-Z = 220.0000 m agl","theta240": "THETA-Z = 240.0000 m agl",
#               "theta260": "THETA-Z = 260.0000 m agl","theta280": "THETA-Z = 280.0000 m agl",
#               "theta300": "THETA-Z = 300.0000 m agl"}
    #vp_fname = "/dagage2/agage/metoffice/vertical_profiles/UKV/GLATTON_Vertical_Profile.txt.gz"  
    
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

