
import datetime

import numpy
import re
from collections import namedtuple

NAMEIII_short_format = '%d/%m/%Y  %H:%M %Z'

def load_NAMEIIITimeseries(filename):
    """
    Loads the Met Office's NAME III time series files returning headers, column 
    definitions and data arrays as 3 separate lists.

    For a technical specification see 
    http://www-hc.metoffice.com/~apdg/nameiii/ModelDocumentation/md2_4_v2%28Output%29.pdf

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

        # read the next 17 lines of header information, putting the form "header name:    header value" into a dictionary
        nlines = 17
        #
        for _ in range(nlines):
            header_name, header_value = file_handle.next().split(':',1)

            # strip off any spurious space characters in the header name and value
            header_name = header_name.strip()
            header_value = header_value.strip()

            # cast some headers into floats or integers if they match a given header name
            if header_name in ['X grid resolution', 'Y grid resolution']:
                header_value = numpy.float64(header_value)
            elif header_name in ['X grid size', 'Y grid size', 'Number of preliminary cols',
                                 'Number of field cols']:
                header_value = int(header_value)

            headers[header_name] = header_value

        # skip the next blank line in the file.
        file_handle.next()

        # skip the next line (contains the word Fields:) in the file.    
        file_handle.next()

        # Read the lines of column definitions - currently hardwired
        column_headings = {}
        for column_header_name in ['Species Category','Name', 'Quantity','Species', 'Unit',
                                   'Sources','Ens Av','Time Av or Int','Horizontal Av or Int',
                                   'Vertical Av or Int','Prob Perc','Prob Perc Ens','Prob Perc Time',
                                   'Location','X','Y','Z','D']:
            column_headings[column_header_name] = [col.strip() for col in file_handle.next().split(',')][1:-1]

        # Determine the coordinates of the data and store in namedtuples
        # First define the namedtuple
        Coord = namedtuple('Coord',['name','dimension','values'])

        # Extract latitude and longitude information from X, Y location headers
        new_Xlocation_column_header = []
        new_Ylocation_column_header = []
        for i, t in enumerate(column_headings['X']):
            if 'Lat-Long' in t:
                p = re.compile(r'\-?[0-9]*\.[0-9]*').search(t)
                new_Xlocation_column_header.append(float(p.group(0)))
            else:
                new_Xlocation_column_header.append(t)
        for i, t in enumerate(column_headings['Y']):
            if 'Lat-Long' in t:
                p = re.compile(r'\-?[0-9]*\.[0-9]*').search(t)
                new_Ylocation_column_header.append(float(p.group(0)))
            else:
                new_Ylocation_column_header.append(t)
        column_headings['X'] = new_Xlocation_column_header
        column_headings['Y'] = new_Ylocation_column_header
        xdim = Coord(name='longitude',dimension=None, values=column_headings['X'])
        ydim = Coord(name='latitude',dimension=None, values=column_headings['Y'])

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

        # skip the blank line after the column headers
        file_handle.next()

        # make a list of data lists to hold the data for each column
        data_lists = [[] for i in range(headers['Number of field cols'])]
        time_list = []

        # iterate over the remaining lines which represent the data in a column form
        for line in file_handle:

            # split the line by comma, removing the last empty column caused by the trailing comma
            vals = line.split(',')[:-1]

            # time is stored in the first column
            t = vals[0].strip()
            time_list.append(datetime.datetime.strptime(t,NAMEIII_short_format))

            # populate the data arrays
            for i in range(headers['Number of field cols']):
                data_lists[i].append(float(vals[i + 1]))

        data_arrays = [numpy.array(l) for l in data_lists]
        time_array = numpy.array(time_list)
        tdim = Coord(name='time',dimension=0, values=time_array)

        coords = [xdim, ydim, tdim]

    return headers, column_headings, coords, data_arrays






    
