#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

rm -f log | true
workdir=`pwd`

cd ${data_dir}
infile=regular_gaussian_model_level.grib1

${tools_dir}grib_ls -P count $infile > log
${tools_dir}grib_ls -p count,step $infile >> log
${tools_dir}grib_ls $infile >> log
${tools_dir}grib_ls -l 0,0,1 $infile >> log
${tools_dir}grib_get -l 0,0,1 $infile >> log
${tools_dir}grib_get -p count,step $infile >> log
${tools_dir}grib_get -P count $infile >> log
 

diff log ls.log 
#rm -f log

cd $workdir

