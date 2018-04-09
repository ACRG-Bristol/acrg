
import datetime

import numpy
import re
from collections import namedtuple

UTC_format = '%H%M%Z %d/%m/%Y'

def load_NAMEIIField(filename):
    """
    Loads the Met Office's NAME III grid output files returning headers, column definitions and data arrays as 3 separate lists.

    For a technical specification see http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf

    """
    # loading a file gives a generator of lines which can be progressed using the next() method.
    # This will come in handy as we wish to progress through the file line by line.
    with open(filename, 'r') as file_handle:

        # define a dictionary which can hold the header metadata about this file
        headers = {}

        # Read in the first line which has a different format to the rest of the header linesn
        header_value = file_handle.next().strip()
        header_name = 'NAME Version'
        headers[header_name] = header_value

        # read the next 16 lines of header information, putting the form "header name:    header value" into a dictionary
        nlines = 16
        
        for _ in range(nlines):
            header_name, header_value = file_handle.next().split(':',1)

            # strip off any spurious space characters in the header name and value
            header_name = header_name.strip()
            header_value = header_value.strip()

            # cast some headers into floats or integers if they match a given header name
            if header_name in ['X grid origin', 'Y grid origin', 'X grid resolution', 'Y grid resolution']:
                header_value = numpy.float64(header_value)
            elif header_name in ['X grid size', 'Y grid size', 'Number of fields',
                                 'Number of preliminary cols','Number of field cols']:
                header_value = int(header_value)

            headers[header_name] = header_value

        # Origin in namever=2 format is bottom-left hand corner so alter this to centre of a grid box
        headers['X grid origin'] = headers['X grid origin'] + headers['X grid resolution']/2
        headers['Y grid origin'] = headers['Y grid origin'] + headers['Y grid resolution']/2

        # skip the next blank line in the file.
        file_handle.next()

        # Read the lines of column definitions. A fixed order and 4 blank columns is assumed.
        column_headings = {}
        for column_header_name in ['Species Category', 'Species', 'Time Av or Int', 'Quantity', 'Unit', 'Z', 'Time']:
            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][4:-1]  
        
        # Determine the coordinates of the data and store in namedtuples
        # First define the namedtuple
        Coord = namedtuple('Coord',['name','dimension','values'])

        # convert the time to python datetimes
        new_time_column_header = []
        for i, t in enumerate(column_headings['Time']):
            new_time_column_header.append(datetime.datetime.strptime(t, UTC_format))
        column_headings['Time'] = new_time_column_header

        # Convert averaging/integrating period to hours
        column_headings['Av or Int period'] = []
        for i, t in enumerate(column_headings['Time Av or Int']):
            p = re.compile(r'\s*(\d{3})\s*(hr)?\s*(time)\s*(\w*)').search(t)
            hours = 0
            if p:
                if len(p.group(1)) > 0:
                    hours = float(p.group(1))
            column_headings['Av or Int period'].append(datetime.timedelta(hours=hours))

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
        data_arrays = [numpy.zeros(data_shape, dtype=numpy.float32) for i in range(headers['Number of fields'])]

        # iterate over the remaining lines which represent the data in a column form
        for line in file_handle:

            # split the line by comma, removing the last empty column caused by the trailing comma
            vals = line.split(',')[:-1]

            # cast the x and y grid positions to floats and convert them to zero based indices
            # In Namever 2 the numbers are 1 based grid positions where 0.5 represents half a grid point.
            #  so just convert vals to zero based integers using int.
            x = int(float(vals[0])) - 1
            y = int(float(vals[1])) - 1

            # populate the data arrays (i.e. all columns but the leading 4)
            for i, data_array in enumerate(data_arrays):
                data_array[y, x] = float(vals[i + 4])

    return headers, column_headings, coords, data_arrays
