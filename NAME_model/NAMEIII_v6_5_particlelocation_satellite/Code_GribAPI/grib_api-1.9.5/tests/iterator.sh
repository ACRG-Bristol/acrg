#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

files="reduced_latlon_surface.grib1 \
      reduced_gaussian_pressure_level.grib1 \
      regular_gaussian_pressure_level.grib1 \
      reduced_gaussian_sub_area.grib1 \
      regular_latlon_surface.grib1 \
      reduced_latlon_surface.grib2  \
      reduced_gaussian_pressure_level.grib2 \
      regular_gaussian_pressure_level.grib2 \
      reduced_gaussian_sub_area.grib2 \
      regular_latlon_surface.grib2"

for f in `echo $files`
do

file=${data_dir}/$f

iterator_count=`${tools_dir}grib_get_data -m 9999:missing -f -p centre -F "%g" -w count=1 $file | grep -v Lat |wc -l `
numberOfPoints=`${tools_dir}grib_get -w count=1 -p numberOfPoints $file`

[ $numberOfPoints = ${iterator_count} ]

done

