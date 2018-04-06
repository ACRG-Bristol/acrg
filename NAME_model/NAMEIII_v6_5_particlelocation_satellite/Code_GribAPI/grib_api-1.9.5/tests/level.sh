#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh
#set -x

outfile=out.grib

file=${data_dir}/regular_gaussian_pressure_level.grib1

cat >good<<EOF
grib1: level=850 isobaricInhPa indicatorOfTypeOfLevel=100
grib2: level=850 isobaricInhPa pl 85000
grib1: level=850 isobaricInhPa indicatorOfTypeOfLevel=100
grib1: level=50 isobaricInPa indicatorOfTypeOfLevel=210
grib2: level=50 isobaricInPa pl 50
grib1: level=50 isobaricInPa indicatorOfTypeOfLevel=210
grib1: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer indicatorOfTypeOfLevel=112
grib2: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer first= sfc 1 second= sfc 2
grib1: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer indicatorOfTypeOfLevel=112
EOF

cat >level.filter<<EOF
set typeOfLevel="isobaricInhPa";
set level=850;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] [typeOfLevel] [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=1;
set typeOfLevel="isobaricInPa";
set level=50;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] [typeOfLevel] [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set typeOfLevel="depthBelowLandLayer";
set topLevel=1;
set bottomLevel=2;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] first= [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface] second= [typeOfSecondFixedSurface:s] [scaledValueOfSecondFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
EOF

    
${tools_dir}grib_filter level.filter $file > test.dump
diff good test.dump

rm -f level.filter good test.dump


