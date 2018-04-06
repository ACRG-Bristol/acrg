#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

if [ ! -f ${data_dir}/multi.grib2 ]
then
  echo no data to test
  exit 0
fi

tmpdata=grib_api.grib

rm -f $tmpdata || true

parameterNumber=`${tools_dir}grib_get -p parameterNumber -w parameterCategory=2,parameterNumber=3 ${data_dir}/multi.grib2`
if [ -z "$parameterNumber"  ] 
then
  echo ---------- grib_get failure
  exit 1
fi

${tools_dir}grib_copy -w parameterCategory=2,parameterNumber=3 ${data_dir}/multi.grib2  $tmpdata.1
${tools_dir}grib_cmp ${data_dir}/v.grib2 $tmpdata.1

cat > $tmpdata.rules <<EOF
if ( parameterCategory==2 && parameterNumber==3) {
  print "[parameterNumber]";
}
EOF

parameterNumber=`${tools_dir}grib_filter $tmpdata.rules ${data_dir}/multi.grib2`
if [ -z "$parameterNumber"  ] 
then
  echo ---------- grib_filter failure
  exit 1
fi

rm -f $tmpdata.1 $tmpdata.rules || true

