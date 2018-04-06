#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

# --- test grib edition 2 to 1 conversions with tigge data
# ---   using the new features of grib_compare that allows
# ---   namespace comparison

. ./include.sh

REDIRECT=/dev/null

#set -eax 

dir="${data_dir}/tigge/"
temp1="temp.grib1_"
temp2="temp.grib2_"

# --- Do I want to exclude any file pattern from the comparison ?
exclusion_pattern="tcw|ssr|str|skt|cap|ci|ttr"

for file in `ls ${dir}/tigge_*.grib`
do

exclude=`echo $file | awk " /$exclusion_pattern/ {print \"found\";} "`

if [ -z "$exclude" ]; then

	rm -f ${temp1} ${temp2} || true

    # 1 to 2 conversion check
	${tools_dir}grib_set -s editionNumber=1 ${file} ${temp1} 2> $REDIRECT > $REDIRECT
	${tools_dir}grib_compare -P -c data:n,geography:n ${temp1} ${file} 2> $REDIRECT > $REDIRECT

    # 2 to 1 conversion check
    ${tools_dir}grib_set -s editionNumber=2 ${temp1} ${temp2} 2> $REDIRECT > $REDIRECT
	${tools_dir}grib_compare -P -c shortName,data:n,geography:n ${temp2} ${file} 2> $REDIRECT > $REDIRECT
fi
done

rm -f ${temp1} ${temp2} || true
