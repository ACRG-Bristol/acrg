#!/bin/ksh
# ******************************************************************************
#
# Project: NAME III documentation --> User Guides --> Testing Runs
#
# File:    shell script for executing NAME test runs
#
# Author:  Andrew Jones, Atmospheric Dispersion, UK Met Office
#
# Date:    31/01/2012
#
# ******************************************************************************

# Get system information and set system variables, etc.

ulimit -s unlimited

if test `uname` = 'Linux'
then
  
  if test `uname -p` = 'x86_64'
  then
    NAMEEXE='nameiii_64bit.exe'
  else
    NAMEEXE='nameiii.exe'
  fi
  
else
  
  echo 'Error: detected system is not Linux - terminating script'
  exit 1
  
fi

# Get the test directory

TESTDIR=`pwd`

#
# Unzip the sample NWP met files
#

echo -e 'Unzipping sample NWP met files (if these are still gzipped) ...'
find ./Sample_NWP_Met_Data -follow -name 'MO*.gz' -type f -exec gunzip {} \;
echo -e '... done\n'

#
# NAME run using single-site met data
#

echo -e 'Performing NAME test run using single-site met data ...'
../../Executables_Linux/${NAMEEXE} Example_SingleSiteMet.txt
echo -e '... done\n'

# Set the NAME output directory.

WORKDIR="${TESTDIR}/Output_SingleSiteMet/"

#
# Generate IDL graphics for single-site met run
#

echo -e 'Generating IDL graphics for NAME run based on single-site met data ...'
./genIDLplots.ksh "${WORKDIR}"
echo -e '... done\n'

#
# Generate python graphics for single-site met run
#

echo -e 'Generating python graphics for NAME run based on single-site met data ...'
./genPythonPlots.ksh "${WORKDIR}"
echo -e '... done\n'

#
# NAME run using NWP met data
#

echo -e 'Performing NAME test run using NWP met data ...'
../../Executables_Linux/${NAMEEXE} Example_NWPMet.txt
echo -e '... done\n'

# Set the NAME output directory.

WORKDIR="${TESTDIR}/Output_NWPMet/"

#
# Generate IDL graphics for NWP met run
#

echo -e 'Generating IDL graphics for NAME run based on NWP met data ...'
./genIDLplots.ksh "${WORKDIR}"
echo -e '... done\n'

#
# Generate python graphics for NWP met run
#

echo -e 'Generating python graphics for NAME run based on NWP met data ...'
./genPythonPlots.ksh "${WORKDIR}"
echo -e '... done\n'

echo "Script $0 completed"

exit 0
