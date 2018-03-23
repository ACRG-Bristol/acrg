#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

files="regular_latlon_surface.grib2 regular_latlon_surface.grib1"

for file in $files
do

cat >statistics.filter<<EOF
set Ni=2;
set Nj=2;

set decimalPrecision=4;
set values={2.0,2.0,2.0,2.0};
print "values=[values]";
print "max=[max] min=[min] average=[average]";

set values={2.0,5.0,2.0,2.0};
print "values=[values]";
print "max=[max] min=[min] average=[average]";

EOF

${tools_dir}grib_filter statistics.filter ${data_dir}/$file > statistics.out

diff statistics.out ${data_dir}/statistics.out.good

done


