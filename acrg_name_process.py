# -*- coding: utf-8 -*-
"""
Created on Thu May  7 10:34:25 2015

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
import xray
import shutil
from scipy.interpolate import interp1d
import copy
import dirsync

#Default NAME output file version
#This is changed depending on presence of "Fields:" line in files
namever=3

#Time formats
UTC_format = '%H%M%Z %d/%m/%Y'
NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'

met_default = {"time": "             T",
               "press": "Pressure (Pa)",
               "temp": "Temperature (C)",
               "PBLH": "Boundary layer depth",
               "wind": "Wind speed",
               "wind_direction": "Wind direction",
               "release_lon": "xxxxxxx",
               "release_lat": "yyyyyyy"}

def load_NAME(file_lines, namever):
    """
    Loads the Met Office's NAME III grid output files returning headers, column definitions and data arrays as 3 separate lists.
    
    For a technical specification see http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf
    
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
            x = float(vals[0]) - 1.5
            y = float(vals[1]) - 1.5
        elif namever == 3:
            x = float(vals[0]) - 1
            y = float(vals[1]) - 1
        
        # populate the data arrays (i.e. all columns but the leading 4) 
        for i, data_array in enumerate(data_arrays):
            data_array[y, x] = float(vals[i + 4])

    return headers, column_headings, data_arrays
    

def extract_file_lines(fname):
    
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
    
    print("Bottom-left grid point (CENTRE): " + str(lons[0]) + "E, " + \
        str(lats[0]) + "N")

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
    
    timeStep=(timeEnd - timeRef).total_seconds()/3600./ntime
    
    # Labelling time steps at END of each time period!
    time = [timeRef + datetime.timedelta(hours=timeStep*(i+1)) \
            for i in range(ntime)]
    
    print("Timestep: %d minutes" % round(timeStep*60))
    print("Levels: %d " % nlevs)
    print("NAME version: %d" % namever)
        
    return lons, lats, levs, time, timeStep


def read_met(fnames):

    if isinstance(fnames, list):
        fnames.sort()
    else:
        fnames=[fnames]

    column_indices = {key: -1 for key, value in met_default.iteritems()}
    
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
                    for key in met_default.keys():
                        if met_default[key] in col:
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
        for key in met_default.keys():
            if column_indices[key] != -1 and key != "time":
                met_dict[key] = m2[:, column_indices[key]].astype(float)
        met_dict["release_lon"] = X
        met_dict["release_lat"] = Y

        #Construct dataframe
        output_df_file = pandas.DataFrame(met_dict,
                          index=[pandas.to_datetime(d, dayfirst=True)
                              for d in m2[:,column_indices["time"]]])

        output_df.append(output_df_file)
    
        print("Read Met file " + fname)
    
    # Concatenate list of data frames
    output_df = pandas.concat(output_df)
    
    # Check for missing values
    output_df = output_df[output_df["press"] > 0.]
    output_df.drop_duplicates(inplace = True)
    
    # Remove duplicate indices (if found)
    output_df = output_df.groupby(level = 0).last()

    # Sort the dataframe by time
    output_df.sort(inplace = True)
    
    output_df.index.name = "time"
    
    return output_df


def particle_locations(particle_file, time, lats, lons, levs, heights,
                       id_is_lev = False):
    '''
    Read and process particle location files.
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
    
    hist = xray.Dataset(variables = {"pl_n":(["time", "lev", "lon", "height"],
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
    
    print("Particle locations " + particle_file)
    
    if particle_file[-3:].upper() == '.GZ':
        compression="gzip"
    else:
        compression=None
    
    df = pandas.read_csv(particle_file, compression=compression, sep=r"\s+")
    
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
        particles = sum([hist_sum[key].values for key in hist_sum.keys()])
        if particles > 0.:
            for key in hist.data_vars.keys():
                hist[key][slice_dict] = hist[key][slice_dict]/\
                                                particles
        else:
            print("WARNING: No particles have reached edge")
            if i > 1:
                print("WARNING: Copying lower level/previous time step")
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
    
    #Check extremes
    if particle_extremes["N"] < edge_lats[1] - dlats/2.:
        print("WARNING: CHECK DOMAIN EDGE TO NORTH")
    if particle_extremes["E"] < edge_lons[1] - dlons/2.:
        print("WARNING: CHECK DOMAIN EDGE TO EAST")
    if particle_extremes["S"] > edge_lats[0] + dlats/2.:
        print("WARNING: CHECK DOMAIN EDGE TO SOUTH")
    if particle_extremes["W"] > edge_lons[0] + dlons/2.:
        print("WARNING: CHECK DOMAIN EDGE TO WEST")

    return hist


def footprint_array(fields_file, 
                    particle_file = None,
                    met = None,
                    satellite = False,
                    time_step = None):

    print("Reading ... " + fields_file)
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
        print("Warning: no particle location file corresponding to " + fields_file)

    nlon=len(lons)
    nlat=len(lats)
    nlev=len(levs)
    ntime=len(time)
    print("Time steps in file: %d" % ntime)
    
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
    fp.time.attrs = {"long_name": "end of time period"}
    
    
    # Add in met data
    met_dict = {key : (["time", "lev"], np.zeros((len(time), len(levs))))
                for key in met[0].keys()}
    met_ds = xray.Dataset(met_dict,
                          coords = {"time": (["time"], time),
                                    "lev": (["lev"], levs)})

    for levi, lev in enumerate(levs):
        metr = met[levi].reindex(index = time, method = "pad")
        for key in met[levi].keys():
            met_ds[key][dict(lev = [levi])] = \
                metr[key].values.reshape((len(time), 1))
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

    # If no meteorology, get null values
    # THIS IS NOT RECOMMENDED. ERRORS of ~ ±10%.
    if met is None:
        met = met_empty()

    # Find footprint files and MATCHING particle location files
    # These files are identified by their date string. Make sure this is right!
    fields_files = sorted(glob.glob(fields_prefix + "*" +
                             datestr + "*.txt*"))

    # Search for particle files                             
    file_datestrs = [f.split(fields_prefix)[-1].split(".txt")[0].split("_")[-1] \
            for f in fields_files]

    particle_files = []
    if particle_prefix is not None:
        for file_datestr in file_datestrs:
            if not glob.glob(particle_prefix + "*" + file_datestr + "*.txt*") is False:
                particle_files.append(
                    glob.glob(particle_prefix + "*" + file_datestr + "*.txt*")[0])
    
        if len(particle_files) != len(fields_files):
            print("Particle files don't match fields files")
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
            particle_locations=None, particle_heights=None):
    
    time_seconds, time_reference = time2sec(time)
    
    #Write NetCDF file
    ncF=Dataset(outfile, 'w')
    ncF.createDimension('time', len(time))
    ncF.createDimension('lon', len(lons))
    ncF.createDimension('lat', len(lats))
    ncF.createDimension('lev', 1)
    
    nctime=ncF.createVariable('time', 'i', ('time',))
    nclon=ncF.createVariable('lon', 'f', ('lon',))
    nclat=ncF.createVariable('lat', 'f', ('lat',))
    nclev=ncF.createVariable('lev', 'str', ('lev',))
    ncfp=ncF.createVariable(varname, 'f', ('lat', 'lon', 'time'), zlib = True,
                            least_significant_digit = 5)
    
    nctime[:]=time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + numpy.str(time_reference)
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
    
    ncF.close()
    print "Written " + outfile


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
            processed_folder = "Processed_Fields_files",
            satellite = False,
            force_update = True):
    
    subfolder = base_dir + domain + "_" + site + "_" + height + "/"
    
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
        datestrs = sorted([f.split("_")[-1].split(".")[0] for f in fields_files])

        # Check that we've found something
        if len(datestrs) == 0:
            print("Error, can't find files in " + file_search_string)
            return None

    else:
        datestrs = [str(year) + str(month).zfill(2)]

    # Output filename
    outfile = subfolder + processed_folder + "/" + site + "-" + height + \
                "_" + domain + "_" + str(year) + str(month).zfill(2) + ".nc"
    
    # Check whether outfile needs updating
    if not force_update:
        print("Testing whether file exists or needs updating: " + outfile)
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

        print("Looking for files with date string: " + datestr + " in " + \
              subfolder)

        # Get Met files
        if satellite:
            met_search_str = subfolder + met_folder + "/*" + datestr + "*/*.txt*"
        else:
            met_search_str = subfolder + met_folder + "/*.txt*"
  
        met_files = sorted(glob.glob(met_search_str))
    
        if len(met_files) == 0:
            print("Can't file MET files: " + met_search_str)
            return None
        else:
            if satellite:
                met = []
                for met_file in met_files:
                    met.append(read_met(met_file))
            else:
                met = read_met(met_files)

        # Get footprints
        fields_prefix = subfolder + fields_folder + "/"
        particles_prefix = subfolder + particles_folder + "/"

        fp_file = footprint_concatenate(fields_prefix,
                                        datestr = datestr, met = met,
                                        particle_prefix = particles_prefix,
                                        satellite = satellite,
                                        time_step = timeStep)
        
        # Do satellite process
        if satellite:
            satellite_obs_file = glob.glob(subfolder + "Observations/*" + \
                                           datestr + "*.nc")
            if len(satellite_obs_file) != 1:
                print("ERROR: There must be exactly one matching satellite " + 
                        "file in the Observations/ folder")
                print("Files: " + satellite_obs_file)
                return None

            fp_file = satellite_vertical_profile(fp_file,
                                                 satellite_obs_file[0])

        if fp_file is not None:
            fp.append(fp_file)
    
    if len(fp) > 0:
        
        # Concatentate
        fp = xray.concat(fp, "time")

        #Write netCDF file
        
        # Define particle locations dictionary (annoying)
        pl = {"N": numpy.transpose(fp.pl_n.values.squeeze(), (2, 1, 0)),
              "E": numpy.transpose(fp.pl_e.values.squeeze(), (2, 1, 0)),
              "S": numpy.transpose(fp.pl_s.values.squeeze(), (2, 1, 0)),
              "W": numpy.transpose(fp.pl_w.values.squeeze(), (2, 1, 0))}

        print("Writing file: " + outfile)
        
        # Write outputs
        write_netcdf(numpy.transpose(fp.fp.values.squeeze(), (1, 2, 0)),
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
                     particle_heights = fp.height.values.squeeze())

    else:
        print("Couldn't seem to find any files")

    return fp


def satellite_vertical_profile(fp, satellite_obs_file):
    '''
    Do weighted average by satellite averaging kernel and
    pressure weight. One time point only. Expects xray.dataset
    with one time point and N vertical levels.
    
    fp = footprint for ONE time point. N levels in lev dimension 
        with footprints and particle locations defined at each level
    satellite_obs_file = NetCDF-format satellite observation file,
        one per time step
    '''
    
    
    def satellite_interpolator(var, press_in, press_out):
        '''
        Interpolator function by pressure for 2D fields
        '''
        
        initial_sum = np.sum(var)
        
        # Check if lowest levels are the same
        # This is an error with NAME vertical profiles where it assumes
        # the pressure of the surface if particles have accidentally been
        # released under ground
        press_in_edit = press_in.copy()
        if np.allclose(press_in[0], press_in[1]):
            press_in_edit[0]-= 10000.
        
        # Add extra level at the bottom for extrapolation
        dp = 10000.
        dvardp_low = (var[1,:,:] - var[0,:,:])/\
                     (press_in_edit[1] - press_in_edit[0])
        var_low = var[0,:,:] + dp*dvardp_low
        press_low = press_in_edit[0] + dp
        press_temp = numpy.concatenate((np.array(press_low).reshape(1),
                                        press_in_edit))
        var_temp = numpy.concatenate((var_low.reshape((1,) + var.shape[1:]),
                                      var))
        
        # Assume top level missing values are zero
        interpolator = interp1d(np.log(press_temp), var_temp,
                                bounds_error = False, axis = 0,
                                fill_value = 0.)

        # Interpolate
        out = interpolator(np.log(press_out))
        
        # Set negative values to zero
        out[out<0.] = 0.

        final_sum = np.sum(out)
        
        print(initial_sum/final_sum)
        
        # Re-scale to correct any difference in total
        return out*initial_sum/final_sum
        
    print("Reading satellite obs file: " + satellite_obs_file)
    sat = xray.open_dataset(satellite_obs_file)
    if np.abs(sat.lon.values[0] - fp.release_lon.values[0,0]) > 1.:
        print("WARNING: Satellite longitude doesn't match footprints")
    if np.abs(sat.lat.values[0] - fp.release_lat.values[0,0]) > 1:
        print("WARNING: Satellite latitude doesn't match footprints")        
    if np.abs(sat.time.values[0] - fp.time.values[0]) > 60*1e9:
        print("WARNING: Satellite time doesn't match footprints")
    if len(fp.time.values) > 1:
        print("ERROR: satellite comparison only for one time step at the moment")
        return fp
    if len(sat.time.values) > 1:
        print("ERROR: satellite comparison only for one time step at the moment")
        return fp

    print(fp.pl_n.sum() + fp.pl_e.sum() + fp.pl_s.sum() + fp.pl_w.sum())
        
    if not np.allclose((fp.pl_n.sum() + fp.pl_e.sum() + \
                        fp.pl_s.sum() + fp.pl_w.sum()), \
                        len(fp.lev)):
        print("ERROR: Particle histograms dont add up to 1 (or nlev)")
        return None

    # Change timestamp to that from obs file
    #  because NAME output only has 1 minute resolution
    fp = fp.reindex_like(sat.time, method = "nearest")
    print(fp.pl_n.sum() + fp.pl_e.sum() + fp.pl_s.sum() + fp.pl_w.sum())
    
    # Interpolate pressure levels
    variables = ["fp", "pl_n", "pl_e", "pl_s", "pl_w"]
    out = {}
    out_interpolated = {}
    particle_total = 0.
    for variable in variables:
        print(variable)
#        return fp[variable].values.squeeze(), \
#                 fp.press.values.squeeze(), \
#                 sat.pressure_levels.values.squeeze()*100.

        #Interpolate fields to satellite data pressure levels
        var = satellite_interpolator(fp[variable].values.squeeze(),
                                     fp.press.values.squeeze(),
                                     sat.pressure_levels.values.squeeze()*100.)
        if variable[0:2] == "pl":
            particle_total += np.sum(var)

        # Store interpolated variable for use in for loop later
        out_interpolated[variable] = var
        
        # Create empty output variable
        out[variable] = np.zeros((1, 1,) + var.shape[1:])
        
    #Check any interpolation errors
    for variable in variables:
        if not np.allclose(particle_total,
                           float(len(fp["lev"]))):
            print("ERROR: Particle count doesn't match lev")
            return None

    # Weight levels using pressure weight and averaging kernel
    sum_ak_pw = np.sum(sat.pressure_weights.values.squeeze() * \
          sat.xch4_averaging_kernel.values.squeeze())
    for variable in variables:
        for lev in range(len(sat.pressure_levels.values.squeeze())):
            out[variable] += out_interpolated[variable][lev, :, :] * \
                             sat.pressure_weights.values.squeeze()[lev] * \
                             sat.xch4_averaging_kernel.values.squeeze()[lev]
    
    # Compress dataset to one level and store column averages
    fp = fp[dict(lev = [0])]
    for variable in variables:
        fp[variable].values = out[variable]
    
    fp["lev"] = numpy.array(["column"])
    
    print(sum_ak_pw)
    print(fp.pl_n.sum() + fp.pl_e.sum() + fp.pl_s.sum() + fp.pl_w.sum())

    return fp
    

def process_agage(domain, site,
                  heights = None,
                  years_in = None, months_in = None,
                  base_dir = "/dagage2/agage/metoffice/NAME_output/",
                  force_update = False, satellite = False):

    acrg_path=split(realpath(__file__))
    
    with open(acrg_path[0] + "/acrg_site_info.json") as f:
        site_info=json.load(f)
        
    # If no height specified, run all heights
    if heights is None:
        heights = site_info[site]["height_name"]
    elif type(heights) is not list:
        heights = [heights]
    
    for height in heights:
        
        subfolder = base_dir + domain + "_" + site + "_" + height + "/"
    
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
            out = process(domain, site, height, year, month,
                base_dir = base_dir, force_update = force_update,
                satellite = satellite)


def met_empty():
    
    met = pandas.DataFrame({key: 0. for key in met_default.keys() if key != "time"},
                          index = [dt.datetime(1900, 1, 1), dt.datetime(2020, 1, 1)])
    met.index.name = "time"
    met["press"] = [100000., 100000.]
    met["temp"] = [280., 280.]
    
    print("WARNING: NO MET")

    return met

# Routine to copy files from:
#   /dagage2/agage/metoffice/NAME_output/DOMAIN_SITE_HEIGHT/Processed_Fields_files
#   to:
#   air:/data/shared/NAME/fp_netcdf/DOMAIN/
def copy_processed(domain):

    src_folder = "/dagage2/agage/metoffice/NAME_output/"
    dst_folder = "/data/shared/NAME/fp/" + domain + "/"
    
    files = glob.glob(src_folder + domain +
        "*/Processed_Fields_files/*.nc")

    folders = set([os.path.split(f)[0] for f in files])

    for f in folders:
        print("Syncing " + f + " and " + dst_folder)
        dirsync.sync(f, dst_folder, "sync")
    print("Done sync")

# Process a list of AGAGE/DECC/GAUGE files if called as main
if __name__ == "__main__":

    domain = "EUROPE"
    sites = ["BSD", "TTA", "RGL", "MHD", "HFD", "TAC",
             "GAUGE-FERRY", "GAUGE-FAAM",
             "EHL", "TIL", "GLA", "WAO", "HAD", "GSN"]
    for site in sites:
        process_agage(domain, site, force_update = True)

