#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

out1=${data_dir}/multi1.out
out2=${data_dir}/multi2.out

for infile in ${data_dir}/multi.grib2 ${data_dir}/reduced_gaussian_model_level.grib1
do

if [ ! -f ${infile} ]
then
  echo no data to test
  exit 0
fi

${test_dir}/multi_from_message -m $infile > $out1
${tools_dir}grib_get -p count,shortName,level $infile > $out2

diff $out1 $out2

${test_dir}/multi_from_message $infile > $out1
${tools_dir}grib_get -M -p count,shortName,level $infile > $out2

diff $out1 $out2

rm -f $out1 $out2 | true

done
