#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

infile=${data_dir}/regular_gaussian_surface.grib1
outfile=${data_dir}/set.grib1

rm -f $outfile || true

${tools_dir}grib_set -v -p levtype,centre,levtype:l,centre:l -s levtype=pl,centre=80 $infile $outfile >$REDIRECT

levtype=`${tools_dir}grib_get -p levtype $outfile`
[ $levtype = "pl" ]

centre=`${tools_dir}grib_get -p centre $outfile`
[ $centre = "cnmc" ]

levtype=`${tools_dir}grib_get -p levtype:l $outfile`
[ $levtype -eq 100 ]

centre=`${tools_dir}grib_get -p centre:l $outfile`
[ $centre -eq 80 ]

infile=${data_dir}/regular_gaussian_surface.grib2
outfile=${data_dir}/set.grib2

rm -f $outfile || true

${tools_dir}grib_set -v -p levtype:l,centre:s  -s levtype:l=100,centre:s=cnmc $infile $outfile >$REDIRECT

levtype=`${tools_dir}grib_get -p levtype $outfile`
[ $levtype = "pl" ]

centre=`${tools_dir}grib_get -p centre $outfile`
[ $centre = "cnmc" ]

levtype=`${tools_dir}grib_get -p levtype:l $outfile`
[ $levtype -eq 100 ]

centre=`${tools_dir}grib_get -p centre:l $outfile`
[ $centre -eq 80 ]

set +e

${tools_dir}grib_set -p levtype $infile $outfile 2> $REDIRECT > $REDIRECT

[ $? -ne 0 ]


