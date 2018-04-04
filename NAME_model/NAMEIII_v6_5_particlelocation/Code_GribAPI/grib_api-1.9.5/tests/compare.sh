#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

infile="${data_dir}/regular_latlon_surface.grib1"
outfile=${infile}.1

rm -f $outfile || true

${tools_dir}grib_set -s shortName=2d $infile $outfile

${tools_dir}grib_compare -b indicatorOfParameter,paramId,shortName $infile $outfile > $REDIRECT

