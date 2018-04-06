#!/bin/sh

. ./include.sh

if [ ! -f "${data_dir}/sample.grib2" ] 
then
  echo SKIP: $0
  exit
fi
${examples_dir}multi_write ${data_dir}/sample.grib2 ${data_dir}/multi_sample.grib2 > /dev/null

${tools_dir}grib_get -p step ${data_dir}/multi_sample.grib2 > ${data_dir}/multi_step.test

diff ${data_dir}/multi_step.test ${data_dir}/multi_step.txt

step=`${tools_dir}grib_get -M -p step ${data_dir}/multi_sample.grib2`

[ $step -eq 12 ]

rm -f ${data_dir}/multi_sample.grib2 ${data_dir}/multi_step.test | true


