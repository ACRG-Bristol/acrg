# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 17:52:10 2014

Process Alistair Manning's NAME II-format output and create netCDF file

I've modified NAME/NAMEIII_v6_1/Code_Python/NAME2NetCDF_incload.py
You can see where I've made changes by looking for the MLR comments
    and the lines that they replace commented out just above

To run for one site:
    import acrg_name_process
    acrg_process_name.process("MHD")
    
To run for all sites:
    import acrg_name_process
    acrg_name_process.run()

NOTE: 2-hours and "small" footprints and NAMEII format
    are hard wired at the moment.

There are lots of undocumented routines

Testing, testing

@author: chxmr
"""

import numpy
from netCDF4 import Dataset
import glob
import gzip
import datetime
import re
import pandas
import scipy.constants as const
from operator import itemgetter
from acrg_grid import areagrid
from acrg_time.convert import time2sec
import os
import matplotlib.pyplot as plt

#Default NAME output file version
#This is changed depending on presence of "Fields:" line in files
namever=3

#Time formats
UTC_format = '%H%M%Z %d/%m/%Y'
NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'

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
        x = float(vals[0]) - 1.5
        y = float(vals[1]) - 1.5
        
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
    lons=numpy.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep) + Xstep/2.
    
    Ycount = header['Y grid size']
    Ystep = header['Y grid resolution']
    Ystart = header['Y grid origin']
    lats=numpy.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep) + Ystep/2

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
        # Assume only one time step per file
        #Levels
        nlevs = len(z_level)
        levs = range(nlevs)
    
        #Time
        time=[datetime.datetime.strptime(header['End of release'], timeFormat)]
        ntime = 1


    #Time steps    
    timeRef=datetime.datetime.strptime(header['End of release'], timeFormat)
    timeEnd=datetime.datetime.strptime(header['Start of release'], timeFormat)
    
    timeStep=(timeEnd - timeRef).total_seconds()/3600/ntime

    time=[]
    if ntime == 1:
        time.append(timeRef)
    else:
        for i in range(ntime):
            time.append(timeRef + datetime.timedelta(hours=timeStep*(i+1)))


    print("Timestep: %d minutes" % round(timeStep*60))
    print("Levels: %d " % nlevs)
    print("NAME version: %d" % namever)
    
    
    return lons, lats, levs, time, timeStep
    

def footprint_array(header, column_headings, data_arrays, \
    met_time=None, pressure = None, temperature = None, area=None):

    STP_warning=False

    lons, lats, levs, time, timeStep = define_grid(header, column_headings)

    if area is None:
        area=areagrid(lats, lons)

    nlon=len(lons)
    nlat=len(lats)
    nlev=len(levs)
    ntime=len(time)
    print("Time steps in file: %d" % ntime)

    z_level=column_headings['z_level'][4:]
    time_column=column_headings['time'][4:]
    
    #Declare footprint array
    fp=numpy.zeros((nlat, nlon, ntime))

    temp=numpy.zeros((ntime))
    press=numpy.zeros((ntime))
    
    #Calculate conversion factor depending on meteorology
    lev=levs[0] #JUST GET FIRST LEVEL
    lev_columns=[i for i, l in enumerate(z_level) if l == lev]
    for i, column in enumerate(lev_columns):
        if pressure is None:
            STP_warning = True
            molm3=44.6429  #mol of air per m3 at STP
            temp[i]=273.
            press[i]=100000.
        else:
            #NOTE: I've added the timeStep to the Met time below
            # because time is labeled at the beginning of the period
            # and footprint is labeled at end
            time_diff=[abs(t + datetime.timedelta(0, 3600*timeStep) - time[i]) 
                for t in met_time]
            meti=min(enumerate(time_diff), key=itemgetter(1))[0]
            if abs(time[i] - met_time[meti]) > \
                    datetime.timedelta(0, timeStep*4*3600):
                print time[i], " don't have met data, assuming STP"
                molm3=44.6429  #mol of air per m3 at STP
                temp[i]=273.
                press[i]=100000.
            else:
                molm3=pressure[meti]/const.R/temperature[meti]
                temp[i]=temperature[meti]
                press[i]=pressure[meti]

        #Convert footprint to (mol/m2/s)^-1
        fp[:, :, i]=data_arrays[column]*area/ \
            (3600.*timeStep*1.)/molm3 #The 1 assumes 1g/s emissions rate

    if STP_warning:
        print("WARNING: No met data, assuming STP")

#This code can be used if you want to output all levels
#    for levi, lev in enumerate(levs):
#        lev_columns=[i for i, l in enumerate(z_level) if l == lev]
#        for i, column in enumerate(lev_columns):
#            fp[:, :, levi, i]=data_arrays[column]
    
    return fp, lons, lats, lev, time, temp, press

def read_met(fnames):
    
    if isinstance(fnames, list):
        fnames.sort()
    else:
        fnames=[fnames]
    
    T=[]
    P=[]
    wind=[]
    wind_direction=[]
    PBLH=[]
    time=[]
    
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

        #Find columns with header names
        timecol = int(numpy.where(m=='                       T')[1])
        Tcol = int(numpy.where(m=='           Temperature (C)')[1])
        Pcol = int(numpy.where(m=='             Pressure (Pa)')[1])
        PBLHcol = int(numpy.where(m=='      Boundary layer depth')[1])
        WINDcol = int(numpy.where(m=='                Wind speed')[1])
        WINDDcol = int(numpy.where(m=='  Wind direction (degrees)')[1])
        
        #Find where data starts
        for i, mi in enumerate(m[:, 0]):
            if str(mi).strip() != '':
                break
        
        m2 = m[i+1:, :]

#        #Read rest of file in so that data is not in string format
#        m2 = pandas.read_csv(fname, skiprows = a+18, 
#                             compression=compression)
#        m2 = numpy.array(m2)
        
        time = time + \
            [pandas.to_datetime(d, dayfirst=True) for d in m2[:,timecol]]

        T=T + list(const.C2K(m2[:, Tcol].astype(float)))
        P=P + list(m2[:, Pcol].astype(float))
        wind=wind + list(m2[:, WINDcol].astype(float))
        wind_direction=wind_direction + list(m2[:, WINDDcol].astype(float))
        PBLH=PBLH + list(m2[:, PBLHcol].astype(float))
        
    T=numpy.array(T)
    P=numpy.array(P)
    PBLH=numpy.array(PBLH)
    wind=numpy.array(wind)
    wind_direction=numpy.array(wind_direction)
        
    if len(time) == 0:
        time=None
        T=None
        P=None
        wind=None
        wind_direction=None
        PBLH=None
        
    return {'time': time, "T": T, "P": P, "PBLH": PBLH, "wind":wind,
             "wind_direction":wind_direction}


def concatenate_footprints(file_list, \
    met_time=None, pressure = None, temperature = None):

    file_list=sorted(file_list)
    nfiles=len(file_list)
    
    #Get first file in order to define data arrays
    header, column_headings, data_arrays = read_file(file_list[0])
    fp_first, lons, lats, lev, time, temp_first, press_first = \
        footprint_array(header, column_headings, data_arrays, \
        met_time=met_time, pressure=pressure, temperature=temperature)

    print("... read " + file_list[0])
    
    fp = [fp_first]
    temp = [temp_first]
    press = [press_first]

    #Calculate area
    area=areagrid(lats, lons)

    levlast=lev

    if nfiles > 1:
        for fi, f in enumerate(file_list[1:]):
            header, column_headings, data_arrays = read_file(f)
            fp_file, lons, lats, lev, time_file, temp_file, press_file = \
                footprint_array(header, column_headings, data_arrays, \
                met_time=met_time, pressure=pressure, temperature=temperature, \
                area=area)
            if levlast <> lev:
                print("WARNING, looks like lowest level has changed: " + fp_file)
            fp.append(fp_file)
            temp.append(temp_file)
            press.append(press_file)
            time=time+time_file
            levlast=lev
            print("... read " + f)

    fp=numpy.dstack(fp)
    temp=numpy.hstack(temp)
    press=numpy.hstack(press)

    return fp, lons, lats, lev, time, temp, press


def write_netcdf(fp, lons, lats, levs, time, outfile, \
            temperature=None, pressure=None, \
            wind_speed=None, wind_direction=None, \
            PBLH=None, varname="fp", particle_locations=None, \
            particle_heights=None):
    
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
    ncfp=ncF.createVariable(varname, 'f', ('lat', 'lon', 'time'))
    
    nctime[:]=time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + numpy.str(time_reference)
    nctime.calendar='gregorian'

    nclon[:]=lons
    nclon.units='Degrees east'
    
    nclat[:]=lats
    nclat.units='Degrees north'

    nclev[:]=numpy.array(levs)
    
    ncfp[:, :, :]=fp
    ncfp.units='(mol/mol)/(mol/m2/s)'

    if temperature is not None:
        nctemp=ncF.createVariable('temperature', 'f', ('time',))
        nctemp[:]=temperature
        nctemp.units='K'

    if pressure is not None:
        ncpress=ncF.createVariable('pressure', 'f', ('time',))
        ncpress[:]=pressure/100.
        ncpress.units='hPa'

    if wind_speed is not None:
        ncwind_speed=ncF.createVariable('wind_speed', 'f', ('time',))
        ncwind_speed[:]=wind_speed
        ncwind_speed.units='m/s'

    if wind_direction is not None:
        ncwind_direction=ncF.createVariable('wind_direction', 'f', ('time',))
        ncwind_direction[:]=wind_direction
        ncwind_direction.units='degrees'

    if PBLH is not None:
        ncPBLH=ncF.createVariable('PBLH', 'f', ('time',))
        ncPBLH[:]=PBLH
        ncPBLH.units='m'

    if particle_locations is not None:
        ncF.createDimension('height', len(particle_locations[0]["N"][0,:]))
        ncHeight=ncF.createVariable('height', 'f',
                                   ('height',))
        ncPartN=ncF.createVariable('particle_locations_n', 'f',
                                   ('height', 'lon', 'time'))
        ncPartE=ncF.createVariable('particle_locations_e', 'f',
                                   ('height', 'lat', 'time'))
        ncPartS=ncF.createVariable('particle_locations_s', 'f',
                                   ('height', 'lon', 'time'))
        ncPartW=ncF.createVariable('particle_locations_w', 'f',
                                   ('height', 'lat', 'time'))
        ncHeight[:]=particle_heights
        ncPartN[:, :, :]=numpy.transpose(
            numpy.dstack([pl["N"] for pl in particle_locations]), (1, 0, 2))
        ncPartE[:, :, :]=numpy.transpose(
            numpy.dstack([pl["E"] for pl in particle_locations]), (1, 0, 2))
        ncPartS[:, :, :]=numpy.transpose(
            numpy.dstack([pl["S"] for pl in particle_locations]), (1, 0, 2))
        ncPartW[:, :, :]=numpy.transpose(
            numpy.dstack([pl["W"] for pl in particle_locations]), (1, 0, 2))
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


def particle_locations(input_search_string, lons = None, lats = None):

    def particle_location_edges(xvalues, yvalues, x, y):
        
        dx = x[1] - x[0]
        xedges = numpy.append(x - dx/2., x[-1] + dx/2.)
        dy = y[1] - y[0]
        yedges = numpy.append(y - dy/2., y[-1] + dy/2.)
        
        hist, xe, ye = numpy.histogram2d(xvalues, yvalues,
                                         bins = (xedges, yedges))
        
        return hist
    
    files=glob.glob(input_search_string)
    files.sort()

    edge_lons = [min(lons), max(lons)]
    edge_lats = [min(lats), max(lats)]
    dlons = lons[1] - lons[0]
    dlats = lats[1] - lats[0]
    dheights = 1000
    
    heights = numpy.arange(0, 15001, 1000) + dheights/2.
    
    # Define output grid
    hist_lats = numpy.hstack([numpy.zeros(len(lats)) + lats[-1],
                           lats[::-1],
                           numpy.zeros(len(lats)) + lats[0],
                           lats] )
    hist_lons = numpy.hstack([lons,
                           numpy.zeros(len(lons)) + lons[-1],
                           lons[::-1],
                           numpy.zeros(len(lons)) + lons[0]])
    hist_heights = heights
    
    hist = []
    particles = []
    
    #Variables to check domain extents
    particle_extremes = {"N": -90., "E": -360. ,"S": 90.,"W": 360.}   
    
    for f in files:
        
        if f[-3:].upper() == '.GZ':
            compression="gzip"
        else:
            compression=None
        
        df = pandas.read_csv(f, compression=compression, sep=r"\s+")
        for i in range(1, len(set(numpy.array(df["Id"])))+1):
        
            hist_ti = {}
            particles_ti = 0.
            
            #Northern edge
            dfe = df[(df["Lat"] > edge_lats[1] - dlats/2.) & (df["Id"] == i)]
            hist_ti["N"] = \
                particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                        lons, heights)
            particles_ti += numpy.sum(hist_ti["N"])
            #Eastern edge
            dfe = df[(df["Long"] > edge_lons[1] - dlons/2.) & (df["Id"] == i)]
            hist_ti["E"] = \
                particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                        lats, heights)
            particles_ti += numpy.sum(hist_ti["E"])
            #Southern edge
            dfe = df[(df["Lat"] < edge_lats[0] + dlats/2.) & (df["Id"] == i)]
            hist_ti["S"] = \
                particle_location_edges(dfe["Long"].values, dfe["Ht"].values,
                                        lons, heights)
            particles_ti += numpy.sum(hist_ti["S"])
            #Western edge
            dfe = df[(df["Long"] < edge_lons[0] + dlons/2.) & (df["Id"] == i)]
            hist_ti["W"] = \
                particle_location_edges(dfe["Lat"].values, dfe["Ht"].values,
                                        lats, heights)
            particles_ti += numpy.sum(hist_ti["W"])

            hist_ti = {direction: hist_direction / particles_ti
                for (direction, hist_direction) in hist_ti.iteritems()}

            hist.append(hist_ti)
            particles.append(particles_ti)

            # Store extremes
            if max(df["Lat"]) > particle_extremes["N"]:
                particle_extremes["N"] = max(df["Lat"])
            if min(df["Lat"]) < particle_extremes["S"]:
                particle_extremes["S"] = min(df["Lat"])
            if max(df["Long"]) > particle_extremes["E"]:
                particle_extremes["E"] = max(df["Long"])
            if min(df["Long"]) < particle_extremes["W"]:
                particle_extremes["W"] = min(df["Long"])
    
    print(particle_extremes)
    
    return hist_lats, hist_lons, hist_heights, hist

def process_satellite_single(input_directory, output_file):
    
    # Get list of met
    met_files = sorted(\
        glob.glob(os.path.join(input_directory, "Met*/*.txt.gz")))

    met = [read_met(f) for f in met_files]
    
    files = glob.glob(os.path.join(input_directory, "Fields_*.txt.gz"))
    files = [f for f in files if "_BL_" not in f]
    
    if len(files) > 1:
        print("More than one file in directory")
        return None
    
    # Get data
    header, column_headings, data_arrays = read_file(files[0])
    lons, lats, levs, time, timeStep = define_grid(header, column_headings, 
                                                   satellite = True)

    fp=numpy.zeros((len(lats), len(lons), len(levs)))
    press = numpy.zeros(len(levs))
    temp = numpy.zeros(len(levs))

    area = areagrid(lats, lons)

    for lev in levs:
        #Convert footprint to (mol/m2/s)^-1
        #The 1 assumes 1g/s emissions rate
        fp[:, :, lev]=data_arrays[lev]*area/ \
            (3600.*timeStep*1.)/ \
            (met[lev]["P"][0]/const.R/met[lev]['T'][0])
        press[lev] = met[lev]["P"][0]
        temp[lev] = met[lev]["T"][0]

    return fp, lons, lats, levs, time, temp, press


def process_multiple(input_search_string, output_file, 
                     met_search_string=None, particle_search_string=None):

    #Get site meteorology
    if met_search_string is not None:
        met_files=sorted(glob.glob(met_search_string))
        met = read_met(met_files)
    else:
        met={"time":None, "T":None, "P":None, 
             "wind": None, "wind_direction": None, "PBLH": None}

    #Search for files
    files=glob.glob(input_search_string)
    files.sort()

    fp, lons, lats, levs, time, temp, press = \
        concatenate_footprints(files,
            met_time=met["time"], pressure=met["P"],
            temperature=met["T"])

    if particle_search_string is not None:
        hist_lats, hist_lons, hist_heights, hist = \
            particle_locations(particle_search_string, lons = lons, lats = lats)
    else:
        hist = None
        hist_heights = None

    write_netcdf(fp, lons, lats, levs, time, output_file,
            temperature=met["T"], pressure=met["P"],
            wind_speed=met["wind"], wind_direction=met["wind_direction"],
            PBLH=met["PBLH"], particle_locations = hist,
            particle_heights = hist_heights)
