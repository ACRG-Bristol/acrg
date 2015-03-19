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

#Default NAME output file version
#This is changed depending on presence of "Fields:" line in files
namever=2

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


def define_grid(header, column_headings):
    
    #Co-ordinates OF CENTRE OF CELL
    #IS THIS SAME FOR NAME II OR NAME III FORMAT????
    Xcount = header['X grid size']
    Xstep = header['X grid resolution']
    Xstart = header['X grid origin']
    lons=numpy.arange(Xstart,Xstart+Xcount*Xstep-0.01*Xstep,Xstep) + Xstep/2.
   
    Ycount = header['Y grid size']
    Ystep = header['Y grid resolution']
    Ystart = header['Y grid origin']
    lats=numpy.arange(Ystart,Ystart+Ycount*Ystep-0.01*Ystep,Ystep) + Ystep/2

    #ESTIMATE NUMBER OF TIMESTEPS PER FILE: ASSUMES ONLY ONE SPECIES PER FILE!!!
    z_level=column_headings['z_level'][4:]
    time=column_headings['time'][4:]
    if namever==2:
        species_name=column_headings['species'][4:]
    else:
        species_name=column_headings['name'][4:]

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

    ntime=len(time)/nlevs

    time=[]
    if namever == 2:
        timeFormat=UTC_format
    else:
        timeFormat=NAMEIII_short_format
        
    timeRef=datetime.datetime.strptime(header['End of release'], timeFormat)
    timeEnd=datetime.datetime.strptime(header['Start of release'], timeFormat)
    
    timeStep=(timeEnd - timeRef).total_seconds()/3600/ntime

    print("Timestep: %d minutes" % round(timeStep*60))
    print("Levels: %d " % nlevs)
    print("NAME version: %d" % namever)
    
    if ntime == 1:
        time.append(timeRef)
    else:
        for i in range(ntime):
            time.append(timeRef + \
                datetime.timedelta(hours=timeStep*(i+1)))

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

        m = pandas.read_csv(fname, skiprows = a, nrows =19, 
                            compression=compression)
        m = numpy.array(m)

        #Find columns with header names
        timecol = int(numpy.where(m=='                       T')[1])
        Tcol = int(numpy.where(m=='           Temperature (C)')[1])
        Pcol = int(numpy.where(m=='             Pressure (Pa)')[1])
        PBLHcol = int(numpy.where(m=='      Boundary layer depth')[1])
        WINDcol = int(numpy.where(m=='                Wind speed')[1])
        WINDDcol = int(numpy.where(m=='  Wind direction (degrees)')[1])
        
        #Read rest of file in so that data is not in string format
        m2 = pandas.read_csv(fname, skiprows = a+18, 
                             compression=compression)
        m2 = numpy.array(m2)
        
        time = time + \
            [pandas.to_datetime(d, dayfirst=True) for d in m2[:,timecol]]

        T=T + list(const.C2K(m2[:, Tcol]))
        P=P + list(m2[:, Pcol])
        wind=wind + list(m2[:, WINDcol])
        wind_direction=wind_direction + list(m2[:, WINDDcol])
        PBLH=PBLH + list(m2[:, PBLHcol])
            
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
            PBLH=None, varname="fp"):
    
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
    
    ncF.close()
    print "Written " + outfile


def name_site_name(site):
    
    name={"MHD": "MACE_HEAD",
          "TTA": "ANGUS_TOWER",
          "RGL": "RIDGEHILL",
          "TAC": "TACOLNESTON",
          "BSD": "BILSDALE",
          "HFD": "HEATHFIELD",
          "THD": "TRINIDAD_HEAD"}
    
    return name[site]


def process_agage(site, suffix="_small"):

    #Get site meteorology
    met_files=glob.glob('/data/shared/NAME/site_meteorology/*' + \
        name_site_name(site) + "*.txt")

    met = read_met(met_files)

    #Search for files
    files=glob.glob("/dagage2/agage/metoffice/NAME_output_" + \
        name_site_name(site) + suffix + "/*.txt.gz")
    files.sort()
    
    year=[int(f[f.rfind('_')+1:f.rfind('_')+5]) for f in files]
    years=set(year)

    for y in years:
        
        #Output file
        outFile='/data/shared/NAME/fp_netcdf/' + \
            site + suffix + '_' + str(y) + '.nc'
        
        #Site meteorology for this year
        if len(met_files) > 0:
            met_time_year, temperature_year, pressure_year = zip(* \
                [met for met in \
                zip(met_time, temperature, pressure) if met[0].year == y])
        else:
            met_time_year=None
            temperature_year=None
            pressure_year=None
        
        #Input files for this year
        fileYear=[f for f, yr in zip(files, year) if yr == y]
        fileYear.sort()
        
        fp, lons, lats, levs, time, temp, press = \
            concatenate_footprints(fileYear, \
                met_time=met_time_year, pressure=pressure_year, \
                temperature=temperature_year)
        
        write_netcdf(fp, lons, lats, levs, time, outFile, \
            temperature=temp, pressure=press)


def run_agage():
    
    sites=['TTA', 'HFD', 'RGL', 'TAC', 'MHD', 'BSD']

    for site in sites:
        process_agage(site, suffix="_small")


#def process_single(infile, outfile, met_file=None):
#
#    if met_file is not None:
#        met = read_met(met_file)
#    else:
#        met={"time":None, "T":None, "P":None}
#    
#    header, column_headings, data_arrays = read_file(infile)
#
#    fp, lons, lats, levs, time, temp, press = \
#        footprint_array(header, column_headings, data_arrays, \
#        met_time=met["time"], pressure=met["P"], temperature=met["T"])
#
#    write_netcdf(fp, lons, lats, levs, time, outfile, \
#            temperature=temp, pressure=press)


def process_multiple(input_search_string, output_file, 
                     met_search_string=None):

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

    write_netcdf(fp, lons, lats, levs, time, output_file,
            temperature=met["T"], pressure=met["P"],
            wind_speed=met["wind"], wind_direction=met["wind_direction"],
            PBLH=met["PBLH"])

