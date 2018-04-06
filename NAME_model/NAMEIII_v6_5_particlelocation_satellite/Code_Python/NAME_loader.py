"""
Loading NAME output files into IRIS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This script loads fields files (XY data) and time series files in NAME II and NAMEIII into an IRIS cube

This is done in three steps:
First the NAME_file_type determines what style of NAME file is being loaded and selects the appropriate loader
Second the data is loaded into python (the loaders are stored in separate files which are imported when
the script runs
Finally NAME_to_cube turns the in memory representation into an Iris cube
"""

from NAMEIII_Field_Loader import load_NAMEIIIField
from NAMEII_Field_Loader import load_NAMEIIField
from NAMEIII_TimeSeries_Loader import load_NAMEIIITimeseries
from NAMEII_TimeSeries_Loader import load_NAMEIITimeseries

import iris
import iris.coords as icoords
import iris.coord_systems as icoord_systems
import iris.fileformats
import iris.io.format_picker as format_picker

import datetime
import numpy

def SIunits(NAMEunit):
    """
    Converting NAME units into standard SI units
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Note that there are several issues here:
    1) Non-SI units which are not recognised by IRIS
    2) Units which are in the wrong case (case is ignored in NAME)
    3) Units where the space between SI units is missing
    4) Units where the characters used are non-standard (i.e. 'mc' for micro instead of 'u')
    """

    unit_mapper = {'Risks/m3': '1',   #Used for Bluetongue
                   'TCID50s/m3': '1', #Used for Foot and Mouth
                   'TCID50/m3': '1',  #Used for Foot and Mouth
                   'N/A': '1',        #Used for CHEMET area at risk
                   'DU': '1',         #Dobson Units
                   'ppm': '1',        #Parts per million
                   'ppb': '1',        #Parts per million
                   'lb' : '1',        #pounds
                   'oz' : '1',        #ounces
                   'Ci' : '1',        #Curies (radioactivity)
                   't'  : '1',        #tonnes
                   'deg': 'degree',   #angular degree
                   'oktas' : '1',     #oktas
                   'deg C' : 'deg_C'  #degrees Celcius
                   }
    new_unit = unit_mapper.get(NAMEunit, NAMEunit)
    
    new_unit = new_unit.replace('Kg', 'kg')
    new_unit = new_unit.replace('gs', 'g s')
    new_unit = new_unit.replace('Bqs', 'Bq s')
    new_unit = new_unit.replace('mcBq', 'uBq')
    new_unit = new_unit.replace('mcg', 'ug')

    return new_unit

def NAME_file_type(filename):
    """ Determines whether a NAME file is in NAME II or NAME III format. """
    with open(filename, 'r') as file_handle:

        head=[file_handle.next().strip() for x in xrange(17)]

        for line in head:
            if line.lstrip().startswith("Run name"):
                for line in head:
                    if line.lstrip().startswith("X grid origin"):
                        headerinfo = line.strip().split(":")
                        if headerinfo[1] == '':
                            fileloader = load_NAMEIIITimeseries
                            break
                        else:
                            fileloader = load_NAMEIIIField
                            break
            elif line.lstrip().startswith("Title"):
                for line in head:
                    if line.lstrip().startswith("Number of series"):
                        fileloader = load_NAMEIITimeseries
                        break
                    else:
                        fileloader = load_NAMEIIField
            elif  "trajectories" in line: 
                print "NAME trajectories are not currently supported in IRIS"
                exit()

    return fileloader

def apply_NAME_loader(f, fname):
    return f(fname)

def NAME_to_cube(filenames, callback):
    """Returns a generator of cubes given a filename and a callback."""

    for filename in filenames:

        fileloader = NAME_file_type(filename)

        header, column_headings, coords, data_arrays = apply_NAME_loader(fileloader,filename)

        for i, data_array in enumerate(data_arrays):
            # turn the dictionary of column headers with a list of header information for each field into a dictionary of
            # headers for just this field.
            field_headings = dict([(k, v[i]) for k, v in column_headings.iteritems()])

            # make an empty cube
            cube = iris.cube.Cube(data_array)

            # define the standard name and unit
            name = ('{} {}'.format(field_headings['Species'], field_headings['Quantity'])).upper().replace(' ', '_')
            cube.rename(name)

            # Some units are not in SI units, are missing spaces or typed in the wrong case. SIunits returns units
            # which are recognised by Iris.
            cube.units = SIunits(field_headings['Unit'])

            # define and add the singular coordinates of the field (flight level, time etc.)
            cube.add_aux_coord(icoords.AuxCoord(field_headings['Z'], long_name='z', units='no-unit'))

            # define the time unit and use it to serialise the datetime for the time coordinate
            time_unit = iris.unit.Unit('hours since epoch', calendar=iris.unit.CALENDAR_GREGORIAN)

            # build a coordinate system which can be referenced by latitude and longitude coordinates
            lat_lon_coord_system = icoord_systems.GeogCS(6371229)

            # build time, latitude and longitude coordinates 
            for coord in coords:
                pts = coord.values
                coord_sys = None
                if  coord.name == 'latitude' or coord.name =='longitude':
                    coord_units = 'degrees'
                    coord_sys = lat_lon_coord_system
                if coord.name == 'time':
                    coord_units = time_unit
                    pts = time_unit.date2num(coord.values)
                if coord.dimension >=0:
                    icoord = icoords.DimCoord(points=pts, 
                                              standard_name=coord.name, 
                                              units=coord_units, 
                                              coord_system=coord_sys)
                    if coord.name == 'time' and 'Av or Int period' in field_headings:
                        date2=coord.values-field_headings['Av or Int period']
                        bnds = time_unit.date2num(numpy.vstack((date2, coord.values)).T)
                        icoord.bounds = bnds
                    else:
                        icoord.guess_bounds(0.5)
                    cube.add_dim_coord(icoord, coord.dimension)
                else:
                    icoord = icoords.AuxCoord(points=pts[i], 
                                              standard_name=coord.name,
                                              coord_system=coord_sys,
                                              units=coord_units)
                    if coord.name == 'time' and 'Av or Int period' in field_headings:
                        date2=coord.values-field_headings['Av or Int period']
                        bnds = time_unit.date2num(numpy.vstack((date2, coord.values)).T)
                        icoord.bounds = bnds[i,:]
                    cube.add_aux_coord(icoord)


            # Add the Main Headings as attributes
            for key, value in header.iteritems():
                cube.attributes[key] = value

            # Add the Column Headings as attributes
            # First remove headings/column headings which have been used in coordinates
            headings = ['Z', 'Unit']
            for heading in headings:
                del field_headings[heading]
            for key, value in field_headings.iteritems():
                cube.attributes[key] = value

            # implement standard iris callback capability. Although callbacks are not used in this example, the standard
            # mechanism for a custom loader to implement a callback is shown:
            cube = iris.io.run_callback(callback, cube, [header, field_headings, coords, data_array], filename)

            # yield the cube created (the loop will continue when the next() element is requested)
            yield cube


# Create a format_picker specification of the NAME file format giving it a priority below NetCDF, GRIB & PP etc.
_NAME_III_spec = format_picker.FormatSpecification('Name III', format_picker.LEADING_LINE,
                                      lambda line: line.lstrip().startswith("NAME III"), NAME_to_cube,
                                      priority=3,)

# Register the NAME loader with iris
iris.fileformats.FORMAT_AGENT.add_spec(_NAME_III_spec)


# ---------------------------------------------
# |          Using the new loader             |
# ---------------------------------------------

def main():
    """ This script is to be used for testing changes to the NAME loading scripts. It will
    load an example of each NAME file type and output the cube(s) contained within the file.
    The output can then be checked against the previous output"""

    # NAME III fields file
    print "Checking NAME III Fields file\n"
    fname = '/project/NAME/apnm/NAMEIIIout/Python_testing/Fields_grid2_C1_T1_201103150000.txt'
    species_constraint = iris.Constraint('LPAR_WET_DEPOSITION')
    deposit = iris.load_cube(fname, species_constraint)
    print deposit
    
    # NAME II fields file
    print "Checking NAME II Fields file\n"
    fname = '/project/NAME/apnm/NAMEIIIout/Python_testing/Fields_grid1_201111010425.txt'
    inert_constraint = iris.Constraint('INERT-TRACER_DOSAGE', z='Boundary layer')
    cube1 = iris.load_cube(fname, inert_constraint)
    print cube1

    # Name II Time series file
    print "Checking NAME II Time Series file\n"
    fname = '/project/NAME/apnm/NAMEIIIout/Python_testing/Time_series_grid1.txt'
    cubelist = iris.load_raw(fname)
    print cubelist[0]
    
    # Name II Time series file
    print "Checking NAME III Time Series file\n"
    fname = '/project/NAME/apnm/NAMEIIIout/Python_testing/TimeSeries1_C1.txt'
    cube = iris.load_raw(fname)
    print cube[9]


if __name__ == '__main__':
    main()
