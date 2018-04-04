#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

cat > rules<<EOF
edition=2;
EOF

${tools_dir}grib_convert -g  rules ${data_dir}/gts.grib ${data_dir}/out.grib2

cgts=`${tools_dir}grib_count ${data_dir}/gts.grib`
cout=`${tools_dir}grib_count ${data_dir}/out.grib2`

if [ $cgts != $cout ]
then
  echo "error in gts test"
  exit 1
fi

rm -f ${data_dir}/out.grib2

exit 0

if [ ! -f "${data_dir}/lfpw.grib1" ] 
then
  echo no data to convert
  exit 0
fi

tmpdata=grib_api.grib

rm -f $tmpdata || true

${tools_dir}grib_convert ${data_dir}/mf.rules ${data_dir}/lfpw.grib1 $tmpdata 2> $REDIRECT > $REDIRECT
${tools_dir}grib_compare -f ${data_dir}/lfpw.grib2 $tmpdata

rm -f $tmpdata || true

