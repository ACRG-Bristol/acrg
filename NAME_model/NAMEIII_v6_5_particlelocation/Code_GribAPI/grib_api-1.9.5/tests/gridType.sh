#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

tmpdata=grib_api.grib

rm -f $tmpdata || true

${tools_dir}grib_set -s gridType=regular_gg ${data_dir}/reduced_latlon_surface.grib1 ${tmpdata} > $REDIRECT

gridType=`${tools_dir}grib_get -p gridType $tmpdata`
if [ $gridType != "regular_gg" ]
then
    echo "Unable to change from reduced_latlon to regular_gg"
    echo $gridType
	exit 1
fi

rm -f $tmpdata || true

${tools_dir}grib_set -s gridType=reduced_gg ${data_dir}/regular_gaussian_pressure_level.grib1 ${tmpdata} > $REDIRECT

gridType=`${tools_dir}grib_get -p gridType $tmpdata`
if [ $gridType != "reduced_gg" ]
then
    echo "Unable to change from regular_gg to reduced_gg"
    echo $gridType
	exit 1
fi

rm -f $tmpdata || true

