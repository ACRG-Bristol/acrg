#!/bin/ksh
# ******************************************************************************
#
# Project: Template files and scripts for research users on NAME-JASMIN
#
# File:    metrestore script for NAME-JASMIN
#
# Author:  Andrew Jones, Atmospheric Dispersion, UK Met Office
#
# Date:    25/09/2014
#
# ******************************************************************************

# This script is designed to run with NAME on the JASMIN platform, and
# restores an archived UM met data file from the NAME-JASMIN met archive.
#
# It first checks for the existence of the file in the local met directory
# (if a local gzipped file exists then that file is unzipped), otherwise
# it searches for the met file in the main met archive and copies/unzips
# the file to the local met directory.


# Retrieve script arguments

if [[ $# -eq 2 ]] ; then
  
  METDIR=$1
  METFILENAME=$2
  
else
  
  echo "Error in script $0: two arguments METDIR and METFILENAME are needed"
  exit 1
  
fi


# Extract SUFFIX from METFILENAME

SUFFIX=`echo ${METFILENAME} | cut -f2- -d'.'`


# Set top-level met archive on JASMIN

#ARCHIVEROOTDIR='/data/shared/NAME/Met/'
ARCHIVEROOTDIR='/met_archive/NAME/Met/'


# Set archive directory for each SUFFIX type

if   [[ ${SUFFIX} == UMG_Mk9_[IM]_L59PT+(\d).pp ]] ; then

  # UMG_Mk9
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/UMG_Mk9PT"

elif   [[ ${SUFFIX} == UMG_Mk8_[IM]_L59PT+(\d).pp ]] ; then
  
  # UMG_Mk8
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/UMG_Mk8PT"
  
elif [[ ${SUFFIX} == UMG_Mk7_[IM]_L59PT+(\d).pp ]] ; then
  
  # UMG_Mk7
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/UMG_Mk7PT"
  
elif [[ ${SUFFIX} == UMG_Mk6_L59PT+(\d).pp      ]] ; then
  
  # UMG_Mk6
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/UMG_Mk6PT"
  
elif [[ ${SUFFIX} == 'UMG_Mk5_L52.pp'           ]] ; then
  
  # UMG_Mk5
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/UMG_Mk5"
  
elif [[ ${SUFFIX} == 'GLOUM6.pp'                ]] ; then
  
  # GLOUM6
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/GLOUM6pp"
  
elif [[ ${SUFFIX} == 'GLOUM6'                   ]] ; then
  
  # GLOUM6
  ARCHIVEMETDIR="${ARCHIVEROOTDIR}/Global/GLOUM6"
  
else
  
  echo "Unknown file suffix = ${SUFFIX} in call to metrestore script"
  exit 2
  
fi


# Check ARCHIVEMETDIR directory is different to METDIR

if [ "${ARCHIVEMETDIR}/" == ${METDIR} ] ; then
 
 echo "!!! ERROR: THIS IS NOT ALLOWED - PLEASE CHANGE LOCAL MET DIRECTORY !!!"
 ARCHIVEMETDIR=MISTAKE
 exit 3
 
fi


# Switch to local met directory

cd ${METDIR}


# Test for presence of met file and attempt to restore from archive if necessary

METFILE=${METDIR}${METFILENAME}
echo "metrestore: looking for ${METFILE}"

if test -f ${METFILE}
then
  
  # file already exists
  echo "metrestore: ${METFILE} exists"
  
elif test -f ${METFILE}.gz
then
  
  # gzip file already exists
  echo "metrestore: unzipping ${METFILE}.gz"
  gunzip ${METFILE}.gz
  
elif test -f ${ARCHIVEMETDIR}/${METFILENAME}.gz
then
  
  # copy met file from archive to working met directory
  echo "metrestore: copying ${METFILENAME}.gz from met archive"
  cp ${ARCHIVEMETDIR}/${METFILENAME}.gz ${METDIR}/${METFILENAME}.gz
  gunzip ${METFILE}.gz
  
else
  
  # met file not found
  echo "metrestore: ${METFILE} not found"
  exit 4
  
fi


# set file permissions for group read-write access

chmod 664 ${METFILE}

exit 0
