#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

[ -z "$GRIB_DEFINITION_PATH" ] | GRIB_DEFINITION_PATH=`${tools_dir}grib_info -d`

for file in `find $GRIB_DEFINITION_PATH -name '*.def' -print`
do
  ${tools_dir}parser $file > $REDIRECT
done

