
import datetime

import numpy
import re
from collections import namedtuple

NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'


def load_NAMEIIIField(filename):
    """
    Loads the Met Office's NAME III grid output files returning headers, column definitions and data arrays as 3 separate lists.

    For a technical specification see http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf

    """
    # loading a file gives a generator of lines which can be progressed using the next() method.
    # This will come in handy as we wish to progress through the file line by line.
    with open(filename, 'r') as file_handle:

        # define a dictionary which can hold the header metadata about this file
        headers = {}

        # Read in the first line which has a different format to the rest of the header lines
        header_value = file_handle.next().strip()
        header_name = 'NAME Version'
        headers[header_name] = header_value

        # read the next 17 lines of header information, putting in the form 
        # "header name:    header value" into a dictionary
        nlines = 17
        
        for _ in range(nlines):
            header_name, header_value = file_handle.next().split(':',1)

            # strip off any spurious space characters in the header name and value
            header_name = header_name.strip()
            header_value = header_value.strip()

            # cast some headers into floats or integers if they match a given header name
            if header_name in ['X grid origin', 'Y grid origin', 'X grid resolution', 'Y grid resolution']:
                header_value = numpy.float64(header_value)
            elif header_name in ['X grid size', 'Y grid size',
                                 'Number of preliminary cols','Number of field cols']:
                header_value = int(header_value)

            headers[header_name] = header_value

        # skip the next blank line in the file.
        file_handle.next()

        # skip the next line (contains the word Fields:) in the file.    
        file_handle.next()

        # Read the lines of column definitions - 
        # In this version a fixed order of column headers is assumed (and first 4 columns are ignored)
        column_headings = {}
        for column_header_name in ['Species Category', 'Name', 'Quantity', 'Species', 'Unit',
                                   'Sources','Ensemble Av','Time Av or Int', 'Horizontal Av or Int',
                                   'Vertical Av or Int','Prob Perc','Prob Perc Ens',
                                   'Prob Perc Time','Time', 'Z', 'D']:
            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][4:-1]  

        # Determine the coordinates of the data and store in namedtuples
        # First define the namedtuple
        Coord = namedtuple('Coord',['name','dimension','values'])

        # convert the time to python datetimes
        new_time_column_header = []
        for i, t in enumerate(column_headings['Time']):
            new_time_column_header.append(datetime.datetime.strptime(t, NAMEIII_short_format))
        column_headings['Time'] = new_time_column_header

        # Convert averaging/integrating period to hours
        column_headings['Av or Int period'] = []
        for i, t in enumerate(column_headings['Time Av or Int']):
            p = re.compile(r'(\d{0,2})(day)?\s*(\d{1,2})(hr)?\s*(\d{1,2})(min)?\s*(\w*)').search(t)
            if p:
                days = 0
                hours = 0
                minutes = 0
                if len(p.group(1)) > 0:
                    days = float(p.group(1))
                if len(p.group(3)) > 0:
                    hours = float(p.group(3))
                if len(p.group(1)) > 0:
                    minutes = float(p.group(5))
                dt = days*24.0 + hours + minutes/60.0
            else:
                dt = 0
            column_headings['Av or Int period'].append(datetime.timedelta(hours=dt))

        # build a time coordinate
        tdim = Coord(name='time',dimension=None, values=numpy.array(column_headings['Time']))

        # build regular latitude and longitude coordinates
        start = headers['X grid origin']
        step = headers['X grid resolution']
        count = headers['X grid size']
        pts = start + numpy.arange(count, dtype=numpy.float64) * step
        xdim = Coord(name='longitude',dimension=1, values=pts)

        start = headers['Y grid origin']
        step = headers['Y grid resolution']
        count = headers['Y grid size']
        pts = start + numpy.arange(count, dtype=numpy.float64) * step
        ydim = Coord(name='latitude',dimension=0, values=pts)

        coords = [xdim, ydim, tdim]

        # skip the blank line after the column headers
        file_handle.next()

        # make a list of data arrays to hold the data for each column
        data_shape = (headers['Y grid size'], headers['X grid size'])
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of field cols'])]

        # iterate over the remaining lines which represent the data in a column form
        for line in file_handle:

            # split the line by comma, removing the last empty column caused by the trailing comma
            vals = line.split(',')[:-1]

            # cast the x and y grid positions to integers and convert them to zero based indices
            x = int(float(vals[0])) - 1
            y = int(float(vals[1])) - 1

            # populate the data arrays (i.e. all columns but the leading 4)
            for i, data_array in enumerate(data_arrays):
                data_array[y, x] = float(vals[i + 4])

    return headers, column_headings, coords, data_arrays


