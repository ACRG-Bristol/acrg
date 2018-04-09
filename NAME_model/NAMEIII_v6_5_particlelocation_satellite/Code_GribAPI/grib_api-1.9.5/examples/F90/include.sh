set -ea

echo
echo "TEST: $0"

if [ -z "${data_dir}" ]
then
  cd ../../
  path=`pwd`
  GRIB_DEFINITION_PATH=$path/definitions
  export GRIB_DEFINITION_PATH
  GRIB_SAMPLES_PATH=$path/samples
  export GRIB_SAMPLES_PATH
  tools_dir=$path/tools/
  examples_dir=$path/examples/F90/
  data_dir=$path/data
else
  echo "Skipping test $0"
  exit
fi

cd ${examples_dir}

if [ -z "${GRIB_API_INCLUDE}" ]
then 
  GRIB_API_INCLUDE=`pwd`/src
fi

if [ -z "${GRIB_API_LIB}" ]
then 
  GRIB_API_LIB=`pwd`/src
fi

#${tools_dir}grib_info

set -u

