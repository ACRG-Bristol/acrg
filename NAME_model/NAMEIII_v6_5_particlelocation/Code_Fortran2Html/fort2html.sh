#!/bin/ksh
#
#####################################################################
# Script for html-ising the NAME-III source code
# (creates links between modules, types, functions and subroutines).
# Andrew Jones, Atmospheric Dispersion Group, 30/01/2003.
# Based on original code by Caroline Woolcock, Local Forecasting.
# Updated on 29/03/2004 and 27/05/2004 to improve pattern matching.
# Updated on 06/05/2005 to remove absolute file links and to add in
# title page generation.
# Updated on 28/07/2005 to use the system date and to take version
# number and input/output directories from settings file. Also tests
# that output directory does not already exist and then creates it.
#####################################################################
# Run info.sh
# Creates information files for the cross-links:
#  modules.txt, types.txt, functions.txt and subroutines.txt.
# It may not produce good output therefore it might be better to
# run it separately and check its output.
# info.sh

####################################################################
#
#  Read in version number and directory locations from settings file
#
####################################################################
for rec in `cat settings.txt`
do
  var_name=`echo $rec | cut -f1 -d=`
  var_name=`echo $var_name | cut -f1 -d.`
  var_value=`echo $rec | cut -f2 -d=`
  if [[ $var_name = 'NAMEIII_VersionNumber' ]]
  then
    version=`echo $var_value`
  fi
  if [[ $var_name = 'Directory_Input' ]]
  then
    InDir=`echo $var_value`
  fi
  if [[ $var_name = 'Directory_Output' ]]
  then
    OutDir=`echo $var_value`
  fi
done

#####################
#
#  Set current date
#
#####################
date=`date '+%d/%m/%Y'`

####################################################################
#
#  Check output directory does not already exist and then create it
#
####################################################################
if [[ -d $OutDir ]]
then
  echo 'Warning: output directory already exists - terminating script to avoid overwrites'
  exit 0
fi
mkdir $OutDir

#############################################
#
#  Copy title page template and keywords file
#
#############################################
cp title_template.txt $InDir/titlepage.html
cp keywords.txt $InDir/keywords.txt

#####################
#
#  Set colour scheme
#
#####################
clr_fun='"#FF6600"' # functions
clr_sub='"#CC9900"' # subroutines
clr_mod='"#009933"' # modules
clr_typ='"#336699"' # types
clr_com='"#000066"' # comments
clr_str='"#990066"' # strings
clr_dir='"#666666"' # compiler directives

############################
#
#  Change working directory
#
############################
cd $InDir

##################################
#
#  Create title page from template
#
##################################
sed -i "s%\*\*Version\*\*%$version%g" titlepage.html
sed -i "s%\*\*Date\*\*%$date%g" titlepage.html

#########################################
#
#  Create basic HTML structure
#
#########################################
echo 'Creating HTML structure'
for f90_filename in `ls *.F90`
do
  echo '<html><head><title>NAME III ('$version'): '$f90_filename'</title></head>' > $f90_filename.html
  echo '<body bgcolor="#ffffff"><pre>' >> $f90_filename.html
  cat $f90_filename >> $f90_filename.html
  echo '</pre></body></html>' >> $f90_filename.html
done
    
#####################################################################
#
#  Insert named anchors for subroutines appearing in subroutines.txt
#
#####################################################################
echo ' - Setting anchors for SUBROUTINES'
for fileline in `cat subroutines.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  sed -i "\%^Subroutine%Is%Subroutine \+\<$procname\>%<font color=$clr_sub>Subroutine <a name=$procname>$procname</a></font>%I" $filename.html
done
#################################################################
#
#  Insert named anchors for functions appearing in functions.txt
#
#################################################################
echo ' - Setting anchors for FUNCTIONS'
for fileline in `cat functions.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  sed -i "\%^Function%Is%Function \+\<$procname\>%<font color=$clr_fun>Function <a name=$procname>$procname</a></font>%I" $filename.html
done
#########################################################
#
#  Insert named anchors for types appearing in types.txt
#
#########################################################
echo ' - Setting anchors for TYPES'
for fileline in `cat types.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  sed -i "\%^Type \+::%Is%Type \+:: \+\<$procname\>%<font color=$clr_typ>Type :: <a name=$procname>$procname</a></font>%I" $filename.html
done
#############################################################
#
#  Insert named anchors for modules appearing in modules.txt
#
#############################################################
echo ' - Setting anchors for MODULES'
for fileline in `cat modules.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  sed -i "\%^Module%Is%Module \+\<$procname\>%<font color=$clr_mod>Module <a name=$procname>$procname</a></font>%I" $filename.html
done

for f90html_filename in `ls *.F90.html`
do
  echo 'Processing '$f90html_filename
#############################################################
#
#  Insert links for subroutines appearing in subroutines.txt
#
#############################################################
  echo ' - Linking SUBROUTINES'
  for fileline in `cat subroutines.txt`
  do
    filename=` echo $fileline | cut -f1 -d:`
    procname=` echo $fileline | cut -f2 -d:`
    sed -i "\%\(name=$procname\)\|\(\!.*\<$procname\>\)\|\('.*\<$procname\>.*'\)%I!s%\([^%]*\)\<$procname\>%\1<a href=$filename.html#$procname>$procname</a>%" $f90html_filename
  done
#########################################################
#
#  Insert links for functions appearing in functions.txt
#
#########################################################
  echo ' - Linking FUNCTIONS'
  for fileline in `cat functions.txt`
  do
    filename=` echo $fileline | cut -f1 -d:`
    procname=` echo $fileline | cut -f2 -d:`
    sed -i "\%\(name=$procname\)\|\(\!.*\<$procname\>\)\|\('.*\<$procname\>.*'\)%I!s%\([^%]*\)\<$procname\>%\1<a href=$filename.html#$procname>$procname</a>%" $f90html_filename
  done    
#################################################
#
#  Insert links for types appearing in types.txt
#
#################################################
  echo ' - Linking TYPES'
  for fileline in `cat types.txt`
  do
    filename=` echo $fileline | cut -f1 -d:`
    procname=` echo $fileline | cut -f2 -d:`
    sed -i "\%\(name=$procname\)\|\(\!.*\<$procname\>\)\|\('.*\<$procname\>.*'\)%I!s%\([^%]*\)\<$procname\>%\1<a href=$filename.html#$procname>$procname</a>%I" $f90html_filename
  done    
#####################################################
#
#  Insert links for modules appearing in modules.txt
#
#####################################################
  echo ' - Linking MODULES'
  for fileline in `cat modules.txt`
  do
    filename=` echo $fileline | cut -f1 -d:`
    procname=` echo $fileline | cut -f2 -d:`
    sed -i "\%\(name=$procname\)\|\(\!.*\<$procname\>\)\|\('.*\<$procname\>.*'\)%I!s%\([^%]*\)\<$procname\>%\1<a href=$filename.html#$procname>$procname</a>%I" $f90html_filename
  done    
done

###################################
#
#  Italicise any comments
#
###################################
echo 'Putting comments in italics'
for f90html_filename in `ls *.F90.html`
do
  sed -i "\%\!%{s%\!%<font color=$clr_com><i>!%;s%$%</i></font>%;}" $f90html_filename
done
###################################
#
#  Change colour of any strings
#
###################################
echo 'Changing colour of strings'
for f90html_filename in `ls *.F90.html`
do
  sed -i "\%\!.*'.*'%!s%\('.*'\)%<font color=$clr_str>\1</font>%" $f90html_filename
done
##################################################################
#
#  Additional emphasis of any $$ comments and compiler directives
#
##################################################################
echo 'Emphasising $$ comments and compiler directives'
for f90html_filename in `ls *.F90.html`
do
  sed -i 's%\(\$\$\)%<font color="#FF0000"><b>\1</b></font>%' $f90html_filename
  sed -i "s%\(# *ifdef\)%<font color=$clr_dir>\1%" $f90html_filename
  sed -i "s%\(# *endif\)%\1</font>%" $f90html_filename
done  
###################################################
#
#  Bold Fortran keywords (defined in keywords.txt)
#
###################################################
echo 'Putting Fortran keywords in bold'
for keyword in `cat keywords.txt`
do
  for f90html_filename in `ls *.F90.html`
  do
    sed -i "\%\(\!.*\<$keyword\>\)\|\('.*\<$keyword\>.*'\)%I!s%\<$keyword\>%<b>$keyword</b>%I" $f90html_filename
  done  
done

################################
#
#  Create links for subroutines
#
################################
echo 'Linking index of SUBROUTINES'
echo '<html><head><title>NAME III ('$version'): List of subroutines</title></head>' > subroutines.html
echo '<body bgcolor="#ffffff"><pre>' >> subroutines.html
for fileline in `cat subroutines.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  echo "<i>$filename</i> <a href=$filename.html#$procname><font color=$clr_sub>$procname</font></a>" >> subroutines.html
done
echo '</pre></body></html>' >> subroutines.html
##############################
#
#  Create links for functions
#
##############################
echo 'Linking index of FUNCTIONS'
echo '<html><head><title>NAME III ('$version'): List of functions</title></head>' > functions.html
echo '<body bgcolor="#ffffff"><pre>' >> functions.html
for fileline in `cat functions.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  echo "<i>$filename</i> <a href=$filename.html#$procname><font color=$clr_fun>$procname</font></a>" >> functions.html
done
echo '</pre></body></html>' >> functions.html
##########################
#
#  Create links for types
#
##########################
echo 'Linking index of TYPES'
echo '<html><head><title>NAME III ('$version'): List of types</title></head>' > types.html
echo '<body bgcolor="#ffffff"><pre>' >> types.html
for fileline in `cat types.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  echo "<i>$filename</i> <a href=$filename.html#$procname><font color=$clr_typ>$procname</font></a>" >> types.html
done
echo '</pre></body></html>' >> types.html
############################
#
#  Create links for modules
#
############################
echo 'Linking index of MODULES'
echo '<html><head><title>NAME III ('$version'): List of modules</title></head>' > modules.html
echo '<body bgcolor="#ffffff"><pre>' >> modules.html
for fileline in `cat modules.txt`
do
  filename=` echo $fileline | cut -f1 -d:`
  procname=` echo $fileline | cut -f2 -d:`
  echo "<i>$filename</i> <a href=$filename.html#$procname><font color=$clr_mod>$procname</font></a>" >> modules.html
done
echo '</pre></body></html>' >> modules.html
######################################
#
#  Move html files to final directory
#
######################################
echo 'Moving files to output directory'
for filename in `ls *.html`
do
  mv $filename $OutDir/$filename
done
exit 0
