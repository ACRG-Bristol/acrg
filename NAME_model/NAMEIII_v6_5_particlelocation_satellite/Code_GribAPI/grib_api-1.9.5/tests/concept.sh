#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

. ./include.sh

data="${data_dir}/*grib1 ${data_dir}/*grib2"
tmp="tmp.out"

for d in $data
do

set -e
${tools_dir}grib_ls -p gridType,packingType -w dataRepresentationType!=192 $d > $tmp

set +e
grep unknown $tmp

if [ $? -eq 0 ] 
then 
  false
  exit
fi

rm -f $tmp | true

done

