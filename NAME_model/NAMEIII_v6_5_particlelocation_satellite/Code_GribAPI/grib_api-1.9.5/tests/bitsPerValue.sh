#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

files="regular_latlon_surface.grib2 \
       regular_latlon_surface.grib1"

for file in `echo $files`
do

infile=${data_dir}/$file
outfile1=${infile}_bitsPerValue_1
outfile2=${infile}_bitsPerValue_2

${tools_dir}grib_set -r -s bitsPerValue=10 $infile $outfile1

${tools_dir}grib_set -s setBitsPerValue=10 $infile $outfile2

${tools_dir}grib_cmp $outfile1 $outfile2

rm -f $outfile1 $outfile2 || true

done

