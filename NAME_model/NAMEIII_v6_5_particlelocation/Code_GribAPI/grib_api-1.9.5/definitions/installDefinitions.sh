#!/bin/sh

if [ $# != 2 ] 
then
echo "
usage: $0 grib_api_installation_dir definition_files_installation_dir
"
exit 1
fi

requiredLibVersion=%LIBRARY_VERSION%

grib_api_dir=$1
grib_api_bin=$1/bin

definitions=$2

if [ ! -f $grib_api_bin/grib_info ]
then
echo "
Unable to find grib_api tools in $grib_api_bin
"
exit 1
fi

set -e
version=`$grib_api_bin/grib_info -v`
defaultDefinitions=`$grib_api_bin/grib_info -d`
set +e

if [ $version != $requiredLibVersion ]
then
echo "
#################################################################
# grib_api version $version found in 
# $grib_api_dir 
# Version $requiredLibVersion is required.
# Installation aborted.
#################################################################
"
exit 1
fi

echo checking definition files compatibility...
for file in `find . -name '*.def' -print`
do
  ${grib_api_bin}/parser $file 
done
if [ $? != 0 ]
then
	echo definition files are not compatible with library version $version
	echo installation aborted
	exit 1
fi

set -e
echo compatibility check ok

echo copying definition files to $definitions.tmp~
[ ! -d $definitions.tmp~ ] || rm -rf $definitions.tmp~
mkdir -p $definitions.tmp~
cp -r * $definitions.tmp~

if [ -d $definitions ]
then
	if [ -d ${definitions}.backup~ ] 
	then
		echo "
#################################################################
# A backup definition files directory is present:
# ${definitions}.backup~
# Please rename or remove it before installing a 
# new version of definition files.
# INSTALLATION ABORTED 
#################################################################
"
  exit 1
	fi
	echo "
#################################################################
# Definition file directory found in 
# ${definitions}
# Moving $definitions to 
# ${definitions}.backup~
#################################################################
	"
	mv $definitions ${definitions}.backup~
fi

echo moving $definitions.tmp~ to ${definitions}
mv $definitions.tmp~ ${definitions}

echo "

Definition files successfully installed in:
${definitions}
"

if [ ${definitions} != $defaultDefinitions ]
then
echo "
## Please remember to set
##    GRIB_DEFINITION_PATH=${definitions}
## to activate the new definition files.
"

fi
