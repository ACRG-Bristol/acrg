#!/bin/sh
# Copyright 2005-2007 ECMWF
#
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.


same=`./same_int_long`

if [ $same -eq 1 ]
then
 interface=grib_f90_int.f90
else
 interface=grib_f90_long_int.f90
fi
cat grib_f90.f90.head $interface grib_f90.f90.tail > grib_f90.f90

