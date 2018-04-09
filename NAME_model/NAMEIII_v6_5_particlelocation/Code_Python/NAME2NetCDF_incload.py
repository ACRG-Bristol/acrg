#!/bin/python
#
# First attempt to create a NetCDF file from a NAME fields file using Python
#
######################################################################################

import numpy
import sys
import datetime
from netCDF4 import Dataset
#from datetime import datetime, timedelta
from netCDF4 import num2date, date2num

UTC_format = '%H%M%Z %d/%m/%Y'
NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'

######################################################################################
def load_NAME_III(filename,namever):
    """
    Loads the Met Office's NAME III grid output files returning headers, column definitions and data arrays as 3 separate lists.
    
    For a technical specification see http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf
    
    """
    # loading a file gives a generator of lines which can be progressed using the next() method. 
    # This will come in handy as we wish to progress through the file line by line.
    file_handle = file(filename)
    
    # define a dictionary which can hold the header metadata about this file
    headers = {}
    
    # skip the NAME header of the file which looks something like 'NAME III (version X.X.X)'
    file_handle.next()
    
    # read the next 16 lines of header information, putting the form "header name:    header value" into a dictionary
    if namever == 2:
        nlines = 16
    elif namever == 3:
        nlines = 17
    for _ in range(nlines):
        header_name, header_value = file_handle.next().split(':',1)

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

    # Origin in namever=2 format is bottom-left hand corner so alter this to centre of a grid box
    if namever == 2:
        headers['X grid origin'] = headers['X grid origin'] + headers['X grid resolution']/2
        headers['Y grid origin'] = headers['Y grid origin'] + headers['Y grid resolution']/2

    # skip the next blank line in the file.    
    file_handle.next()
        
    # Read the lines of column definitions
    if namever == 2:
        column_headings = {}
        for column_header_name in ['species_category', 'species', 'cell_measure', 'quantity', 'unit', 'z_level', 'time']:
            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][:-1]
    elif namever == 3:
        # skip the next line (contains the word Fields:) in the file.    
        file_handle.next()
        column_headings = {}
        for column_header_name in ['category', 'name', 'quantity', 'species', 'unit','source','ensemble','time_averaging',
                                   'horiz_averaging','vert_averaging','prob_percentile','prob_percentile_ensemble',
                                   'prob_percentile_time','time', 'z_level', 'D']:
            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][:-1]  

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
    file_handle.next()
    
    # make a list of data arrays to hold the data for each column 
    data_shape = (headers['Y grid size'], headers['X grid size'])
    if namever==2:
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of fields'])]
    elif namever==3:
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of field cols'])]
      
    # iterate over the remaining lines which represent the data in a column form
    for line in file_handle:
        
        # split the line by comma, removing the last empty column caused by the trailing comma
        vals = line.split(',')[:-1]
        
        # cast the x and y grid positions to floats and convert them to zero based indices
        # (the numbers are 1 based grid positions where 0.5 represents half a grid point.)
        x = int(float(vals[0])) - 1
        y = int(float(vals[1])) - 1
        
        # populate the data arrays (i.e. all columns but the leading 4) 
        for i, data_array in enumerate(data_arrays): 
            data_array[y, x] = float(vals[i + 4])
    
    return headers, column_headings, data_arrays
    
    

######################################################################################
def AddVariableII(rootgrp,column_headings,data_arrays,fieldname,index):

   field01 =  rootgrp.createVariable(fieldname,'f4',('time','latitude','longitude',))
   field01.category = column_headings['species_category'][index+4]
   field01.species = column_headings['species'][index+4]
   field01.time_averaging_integrating = column_headings['cell_measure'][index+4]
   field01.quantity = column_headings['quantity'][index+4]
   field01.units = column_headings['unit'][index+4]
   field01.Z_plus_Z_averaging = column_headings['z_level'][index+4]
   STRtime = column_headings['time'][index+4].strftime('%H%M%z %d/%m/%y')
   field01.time = STRtime

   field01[0,:,:] = data_arrays[index]

   return

######################################################################################
def AddVariableIII(rootgrp,column_headings,data_arrays,fieldname,index):

   field01 =  rootgrp.createVariable(fieldname,'f4',('time','latitude','longitude',))
   field01.category = column_headings['category'][index+4]
   field01.name = column_headings['name'][index+4]
   field01.quantity = column_headings['quantity'][index+4]
   field01.species = column_headings['species'][index+4]
   field01.units = column_headings['unit'][index+4]
   field01.source_group = column_headings['source'][index+4]
   field01.ensemble_averaging = column_headings['ensemble'][index+4]
   field01.time_averaging_integrating = column_headings['time_averaging'][index+4]
   field01.horiz_averaging_integrating = column_headings['horiz_averaging'][index+4]
   field01.vert_averaging_integrating = column_headings['vert_averaging'][index+4]
   field01.prob_percentile = column_headings['prob_percentile'][index+4]
   field01.prob_percentile_ensemble = column_headings['prob_percentile_ensemble'][index+4]
   field01.prob_percentile_time = column_headings['prob_percentile_time'][index+4]
   STRtime = column_headings['time'][index+4].strftime('%H%M%z %d/%m/%y')
   field01.output_time = STRtime
   field01.Z = column_headings['z_level'][index+4]
   field01.D = column_headings['D'][index+4]


   field01[0,:,:] = data_arrays[index]

   return

#######################################################################################
def AddGlobalAttribII(rootgrp,header):

   #Add global attributes from header information
   rootgrp.Title = header['Title']
   #STRRunTime = header['Run time'].strftime('%H%M%z %d/%m/%y')
   rootgrp.RunTime = header['Run time'] #STRRunTime
   rootgrp.MetData = header['Met data']
   #STRStartofRelease = header['Start of release'].strftime('%H%M%z %d/%m/%y')
   rootgrp.StartofRelease = header['Start of release'] #STRStartofRelease
   #STREndofRelease = header['End of release'].strftime('%H%M%z %d/%m/%y')
   rootgrp.EndofRelease = header['End of release'] #STREndofRelease
   rootgrp.ReleaseRate = header['Release rate']
   rootgrp.ReleaseLocation = header['Release location']
   rootgrp.ReleaseHeight = header['Release height']
   rootgrp.ForecastDuration = header['Forecast duration']

   return

#######################################################################################
def AddGlobalAttribIII(rootgrp,header):

   #Add global attributes from header information
   rootgrp.RunName = header['Run name']
   rootgrp.RunTime = header['Run time']
   rootgrp.MetData = header['Met data']
   rootgrp.StartofRelease = header['Start of release']
   rootgrp.EndofRelease = header['End of release']
   rootgrp.Sourcestrngth = header['Source strength']
   rootgrp.ReleaseLocation = header['Release location']
   rootgrp.ReleaseHeight = header['Release height']
   rootgrp.RunDuration = header['Run duration']

   return

#######################################################################################
def main():
   # Read filename from argument list
   fname=sys.argv[1]
   namever=int(sys.argv[2])
   outname=sys.argv[3]

   #Use AVD's custom loader to load in NAME data
   header, column_headings, data_arrays = load_NAME_III(fname,namever)

   # Create NetCDF file
   rootgrp = Dataset(outname,'w',format='NETCDF3_CLASSIC')
   
   # Use information from NAME file to create dimensions
   rootgrp.createDimension('time', None)
   Xcount = header['X grid size']
   rootgrp.createDimension('longitude',Xcount)
   Ycount = header['Y grid size']
   rootgrp.createDimension('latitude',Ycount)
   
   #Use information from NAME file to create variables
   longitudes = rootgrp.createVariable('longitude','f8',('longitude',))
   latitudes = rootgrp.createVariable('latitude','f8',('latitude',))
   times = rootgrp.createVariable('time','f8',('time',))

   #Add attributes for latitude and longitude arrays
   longitudes.units="degrees_east"
   longitudes.long_name="longitude degrees east from the greenwich meridian"
   longitudes.point_spacing="even"
   latitudes.units="degrees_north"
   latitudes.long_name="latitude degrees north from the equator"
   latitudes.point_spacing="even"

   #Populate latitudes and longitudes (Note that the -0.1*X(Y)Step is required due to rounding issues)
   Xstep = header['X grid resolution']
   Xstart = header['X grid origin']
   lons=numpy.arange(Xstart,Xstart+Xcount*Xstep-0.1*Xstep,Xstep)
   
   Ystep = header['Y grid resolution']
   Ystart = header['Y grid origin']
   lats=numpy.arange(Ystart,Ystart+Ycount*Ystep-0.1*Ystep,Ystep)

   latitudes[:]=lats
   longitudes[:]=lons
   
   #These attributes are critical for creating the time array
   times.units = 'hours since 1970-01-01 00:00:00.0'
   times.calendar = 'gregorian'
   times[0] = date2num(column_headings['time'][4],units=times.units,calendar=times.calendar)
   
   #Add header information to global attributes
   if namever == 2:
      AddGlobalAttribII(rootgrp,header)
   else:
      AddGlobalAttribIII(rootgrp,header)
   
   #Add data and column headings to NetCDF file   
   nvar = len(data_arrays)
   for i in range(nvar):
      i1 = i + 1
      fieldname = 'field'+str('%02d' % i1)
      if namever == 2:
         AddVariableII(rootgrp,column_headings,data_arrays,fieldname,i)
      else:
         AddVariableIII(rootgrp,column_headings,data_arrays,fieldname,i)
         
   #Create data array
   nlats = len(rootgrp.dimensions['latitude'])
   nlons = len(rootgrp.dimensions['longitude'])
   
   #Close NetCDF file
   rootgrp.close()

if __name__ == '__main__':
    main()
