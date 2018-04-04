set -ea

echo
echo "TEST: $0"

path=$TOPBUILDDIR
GRIB_DEFINITION_PATH=$path/definitions
export GRIB_DEFINITION_PATH
GRIB_SAMPLES_PATH=$path/samples
export GRIB_SAMPLES_PATH
tools_dir=$path/tools/
examples_dir=$path/examples/python
data_dir=$path/data

PYTHONPATH=$path/python:$PYTHONPATH
export PYTHONPATH

set -u

