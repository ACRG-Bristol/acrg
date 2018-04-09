#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

files="constant_field\
 reduced_gaussian_pressure_level_constant \
 reduced_latlon_surface_constant \
 regular_gaussian_pressure_level_constant \
 regular_latlon_surface_constant \
 reduced_latlon_surface \
 reduced_gaussian_pressure_level \
 reduced_gaussian_sub_area \
 regular_gaussian_pressure_level \
 regular_latlon_surface \
 reduced_gaussian_model_level \
 regular_gaussian_model_level \
 reduced_gaussian_surface \
 regular_gaussian_surface \
 spherical_pressure_level \
 spherical_model_level "

for f in `echo $files`
do

file=${data_dir}/$f

rm -f ${file}.grib1_ || true

${tools_dir}grib_set -s editionNumber=1 ${file}.grib2 ${file}.grib1_ 2> $REDIRECT > $REDIRECT


grib1Statistics=`${tools_dir}grib_get -fp numberOfValues,numberOfPoints,max,min,average,numberOfMissing ${file}.grib1_` 
grib2Statistics=`${tools_dir}grib_get -fp numberOfValues,numberOfPoints,max,min,average,numberOfMissing ${file}.grib2` 

if [ "$grib1Statistics" != "$grib2Statistics" ]
then 
  exit 1
fi

#${tools_dir}grib_compare -A1.0e-8 -c values ${file}.grib1_ ${file}.grib2 2> /dev/null > /dev/null

rm -f ${file}.grib1_ || true

done

