#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


#set -x
. ./include.sh
infile=${data_dir}/index.grib

if [ ! -f ${infile} ]
then
  echo no data to test
  exit 0
fi

${test_dir}/index ${infile} > index.out

diff index.out ${data_dir}/index.ok

${test_dir}/read_index ${infile} > index.out

diff index.out ${data_dir}/index.ok

rm -f index.out

