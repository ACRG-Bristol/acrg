#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

. ./include.sh

infile=${data_dir}/mixed.grib
outfile=${data_dir}/read_any.out

if [ ! -f ${infile} ]
then
  echo no data to test
  exit 0
fi

$test_dir/read_any $infile > ${outfile}

diff ${outfile} ${data_dir}/read_any.ok

rm -f ${outfile}

count=`${tools_dir}/grib_count ${data_dir}/pad.grib`

if [ $count != 6 ]
then
  echo grib_io problem
  exit 1
fi

