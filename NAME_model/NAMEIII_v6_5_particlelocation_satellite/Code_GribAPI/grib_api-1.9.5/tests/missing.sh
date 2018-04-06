#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

infile="${data_dir}/missing.grib2"
outfile="${data_dir}/missing_new.grib2"

scaleFactorOfSecondFixedSurface=`${tools_dir}grib_get -w count=1 -p scaleFactorOfSecondFixedSurface $infile`
[ "$scaleFactorOfSecondFixedSurface" = "0" ]

scaledValueOfSecondFixedSurface=`${tools_dir}grib_get -w count=1 -p scaledValueOfSecondFixedSurface $infile`
[ "$scaledValueOfSecondFixedSurface" = "0" ]

${tools_dir}grib_set -s scaleFactorOfSecondFixedSurface=missing,scaledValueOfSecondFixedSurface=missing $infile $outfile

scaleFactorOfSecondFixedSurface=`${tools_dir}grib_get -w count=1 -p scaleFactorOfSecondFixedSurface $outfile`
[ "$scaleFactorOfSecondFixedSurface" = "MISSING" ]

scaledValueOfSecondFixedSurface=`${tools_dir}grib_get -w count=1 -p scaledValueOfSecondFixedSurface $outfile`
[ "$scaledValueOfSecondFixedSurface" = "MISSING" ]

set +e

${tools_dir}grib_set -s centre=missing $infile $outfile 2> $REDIRECT > $REDIRECT

[ $? -ne 0 ]

rm -f $outfile

