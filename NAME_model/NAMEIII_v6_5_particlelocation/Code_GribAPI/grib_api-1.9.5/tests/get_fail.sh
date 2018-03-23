#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

set +e 

${tools_dir}grib_get -p gribname ${data_dir}/regular_latlon_surface.grib1 2> $REDIRECT > $REDIRECT

if [ $? -eq 0 ] 
  then
  exit 1;
fi


