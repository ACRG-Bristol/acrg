# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

set -ea

echo
echo "TEST: $0"

data_dir=""

# save current working dir
save=`pwd`

if [ -z "${data_dir}" ]
then
  cd ../
  path=`pwd`
  GRIB_DEFINITION_PATH=$path/definitions
  export GRIB_DEFINITION_PATH
  GRIB_SAMPLES_PATH=$path/samples
  export GRIB_SAMPLES_PATH
  tools_dir=$path/tools/
  tigge_dir=$path/tigge/
  data_dir=$path/data
  test_dir=$path/tests
else
  tools_dir=""
  tigge_dir=""
fi

if [ -z "${GRIB_API_INCLUDE}" ]
then 
  GRIB_API_INCLUDE=`pwd`/src
fi

if [ -z "${GRIB_API_LIB}" ]
then 
  GRIB_API_LIB=`pwd`/src
fi

#${tools_dir}grib_info

# go back to current working dir
cd $save

set -u

