#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.



. ./include.sh

REDIRECT=/dev/null

files="regular_latlon_surface.grib2 \
       regular_latlon_surface.grib1"

for file in `echo $files`
do

infile=${data_dir}/$file
outfile1=${infile}_decimalPrecision_1
outfile2=${infile}_decimalPrecision_2

${tools_dir}grib_set -r -s decimalScaleFactor=1,bitsPerValue=0 $infile $outfile1

${tools_dir}grib_compare -b totalLength,section4Length,section7Length,dataLength,decimalScaleFactor,binaryScaleFactor,referenceValue,bitsPerValue $infile $outfile1 > $REDIRECT

${tools_dir}grib_set -s changeDecimalPrecision=1 $infile $outfile2

${tools_dir}grib_compare -b totalLength,section4Length,section7Length,dataLength,decimalScaleFactor,binaryScaleFactor,referenceValue,bitsPerValue $infile $outfile2 > $REDIRECT

${tools_dir}grib_compare $outfile1 $outfile2

rm -f $outfile1 $outfile2 || true

done

