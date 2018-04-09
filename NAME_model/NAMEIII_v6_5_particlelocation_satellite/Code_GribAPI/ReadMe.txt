Read me file for ECMWF GRIB_API software.

Components: 
----------  

The ECMWF GRIB API is an application program interface accessible from C and FORTRAN programs
developed for encoding and decoding WMO FM-92 GRIB edition 1 and edition 2 messages.

Version 1.9.5, downloaded 14/12/2010 (see http://www.ecmwf.int/products/data/software/download/grib_api.html).
Any local modifications, local tables, etc. are indicated in the changes file Changes.txt.

ReadMe.txt  } This file

Changes.txt } Changes file

LinuxIntelRelease } Standard compile script created by us to simplify compilation process

lgpl-3.0.txt } Licence files (see below)
gpl-3.0.txt  }


Documentation: 
-------------  

Documentation of the GRIB API software can be found via the links
http://www.ecmwf.int/products/data/software/grib_api.html
http://www.ecmwf.int/publications/manuals/grib_api/index.html


Licence:
-------

The ECMWF GRIB API software is licensed under the GNU Lesser General Public License (see lgpl-3.0.txt)
which incorporates the terms and conditions of version 3 of the GNU General Public License (gpl-3.0.txt).


Compilation:
-----------

The source files are in ...\grib_api-1.9.5

The standard compile script LinuxIntelRelease has been added to automate the compilation process.
Source files are compiled into the shared object libraries libgrib_api.so, libgrib_api_f77.so
and libgrib_api_f90.so and then installed, along with appropriate tables and templates, into
...\SharedLibraries_Linux (see ...\grib_api-1.9.5\README for further details on compilation, etc.).
'_64bit' is added to the names of the libraries when compiled on 64-bit systems.

The shared libraries, tables and templates are then available for linking into NAME III at run time.


Run-time linking:
----------------
The following environment variables need to be set to dynamically link the GRIB API at run time:

 - GRIB_DEFINITION_PATH=${SHAREDLIB_DIR}/share/definitions
 - GRIB_SAMPLES_PATH=${SHAREDLIB_DIR}/share/samples
 - GRIB_API_INCLUDE=${SHAREDLIB_DIR}/include
 - GRIB_API_LIB=${SHAREDLIB_DIR}/lib

where ${SHAREDLIB_DIR} is the top-level directory of the GRIB API installation, and also
${GRIB_API_LIB} then needs to be added to the 'LD_LIBRARY_PATH' environment variable.
