#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


. ./include.sh
${test_dir}/julian > ${data_dir}/julian.out

diff ${data_dir}/julian.out ${data_dir}/julian.out.good

rm -f ${data_dir}/julian.out | true
