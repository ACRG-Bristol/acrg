#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh

REDIRECT=/dev/null

grib1=${data_dir}/regular_latlon_surface.grib1
grib2=${data_dir}/regular_latlon_surface.grib2
infile=${data_dir}/reduced_gaussian_model_level.grib1
outfile=${data_dir}/with_bitmap.grib1
outfile1=${data_dir}/without_bitmap.grib1

rm -f $outfile || true

${tools_dir}grib_set -s bitmapPresent=1 $infile $outfile >$REDIRECT

${tools_dir}grib_dump -O $infile  | grep -v FILE > $infile.dump
${tools_dir}grib_dump -O $outfile | grep -v FILE  > $outfile.dump

diff $outfile.dump ${data_dir}/bitmap.diff
diff $infile.dump ${data_dir}/no_bitmap.diff

rm -f $infile.dump $outfile.dump || true

${tools_dir}grib_set -s bitmapPresent=0 $outfile $outfile1 >$REDIRECT

${tools_dir}grib_dump -O $outfile1 | grep -v FILE > $outfile1.dump
${tools_dir}grib_dump -O $outfile  | grep -v FILE> $outfile.dump

diff $outfile1.dump ${data_dir}/no_bitmap.diff

rm -f $outfile1 $outfile1.dump $outfile $outfile.dump || true

cat > bitmap.rules<<EOF
set bitmapPresent=1;
set missingValue=1111;
set Ni=6;
set Nj=2;
set values={1,2,3,4,5,6,7,1111,1111,8,9,10};
write ;
EOF

${tools_dir}grib_filter -o out.grib1 bitmap.rules $grib1 
${tools_dir}grib_filter -o out.grib2 bitmap.rules $grib2
exit 0

${tools_dir}grib_get_data -m missing out.grib1 > out.grib1.data
${tools_dir}grib_get_data -m missing out.grib2 > out.grib2.data

diff out.grib1.data out.grib2.data

rm -f  out.grib1.data out.grib2.data out.grib1 out.grib2 bitmap.rules


