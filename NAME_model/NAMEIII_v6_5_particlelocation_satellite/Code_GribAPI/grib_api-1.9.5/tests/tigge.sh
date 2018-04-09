#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

. ./include.sh

REDIRECT=/dev/null

#set -eax 

dir="${data_dir}/tigge/"


for file in `ls ${dir}/tigge_*.grib`
do
	${tigge_dir}tigge_check ${file} 2> $REDIRECT > $REDIRECT
#	${tigge_dir}tigge_check ${file} 
done

