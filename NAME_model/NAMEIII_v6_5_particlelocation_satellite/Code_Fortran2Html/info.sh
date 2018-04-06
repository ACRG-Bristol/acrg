#!/bin/ksh
#
##################################################################
# Script for extracting structural information from NAME-III code
# Andrew Jones, Atmospheric Dispersion Group, 30/01/2003.
#
# Updated on 28/07/2005 to use input directory from settings file.
# Updated on 29/07/2005 to improve pattern matching on objects -
# it is now better at excluding any types, etc. within a comment.
# Updated on 21/03/2006 to build listing of warning/error messages
# directly in html format (so that such a listing can be produced
# without needing to run the main fort2html processing routine
# although, of course, full cross-referencing needs the latter!).
##################################################################
#
#
# Read version number and input directory location from settings file
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
done

# Check for lock file in input directory (a lock file is created
# on first execution of this script and prevents a subsequent run
# - since named anchors will have been created in the initial
# Fortran code and rerunning will mess up the html links!).
if [[ -f $InDir/lock.txt ]]
then
  echo 'Please reset the Fortran files for processing and then delete the lock file'
  exit 1
else
  touch $InDir/lock.txt
fi

#
#
# Change working directory
cd $InDir
#
#
# Extract subroutine information
echo "CREATING subroutines.txt"
grep -i '^[^!]*SUBROUTINE' *.F90 | grep -vi 'END SUBROUTINE' > subroutines.txt
cut -f1 -d'(' subroutines.txt > subroutines.txt.temp
sed "s/\<SUBROUTINE\>//I" subroutines.txt.temp > subroutines.txt.temp2
sed "s/ //g" subroutines.txt.temp2 > subroutines.txt
#
#
# Extract function information
echo "CREATING functions.txt"
grep -i '^[^!]*FUNCTION' *.F90 | grep -vi 'END FUNCTION' > functions.txt
cut -f1 -d'(' functions.txt > functions.txt.temp
sed "s/\<FUNCTION\>//I" functions.txt.temp > functions.txt.temp2
sed "s/\<REAL\>//I" functions.txt.temp2 > functions.txt.temp
sed "s/ //g" functions.txt.temp > functions.txt
#
#
# Extract type information
echo "CREATING types.txt"
grep -i '^[^!]*TYPE ::' *.F90 | grep -vi 'END TYPE' > types.txt
cut -f1 -d'!' types.txt > types.txt.temp
sed "s/TYPE :://I" types.txt.temp > types.txt.temp2
sed "s/ //g" types.txt.temp2 > types.txt
#
#
# Extract module information
echo "CREATING modules.txt"
grep -i '^[^!]*MODULE' *.F90 | grep -vi 'END MODULE' | grep -vi 'USE ' | grep -vi 'MODULE PROCEDURE' | grep -vi "'" | grep -vi 'CALL' | grep -vi 'SUBROUTINE' | grep -vi 'PRIVATE.*MODULE' | grep -vi 'PUBLIC.*MODULE' > modules.txt
sed "s/\<MODULE\>//I" modules.txt > modules.txt.temp
sed "s/ //g" modules.txt.temp > modules.txt
#
#
# Remove temporary files
rm *.temp 
rm *.temp2

# Script for producing list of error messages from NAME III code
# Added 20/03/2006
echo "CREATING errors.html"
echo "Note: any unprocessed messages will appear here"
quote='"'

# Create opening section of the listing 
echo '<html><head><title>NAME III ('$version'): List of warnings and error messages</title></head>' > errors.html
echo '<body bgcolor="#ffffff"><pre>' >> errors.html
echo '<div align="center">' >> errors.html
echo '<h1> NAME III ('$version'): List of warnings and error messages </h1>' >> errors.html
echo '</div>' >> errors.html
echo '<!-- Page Index -->' >> errors.html
echo '<hr>' >> errors.html
echo '' >> errors.html

# Loop over all .F90 files in input directory
for f90_filename in `ls *.F90`
do

# Add subtitle and index link for this file
echo '<h3><a name="'${f90_filename}'"> '${f90_filename}' </h3>' >> errors.html
sed -i "/<!-- Page Index -->/ {
i\
<a href=${quote}#${f90_filename}${quote}>${f90_filename}</a>\

}" errors.html

# Add identifiers (= line numbers) to each error message
# Note: this modifies the Fortran in the input directory!
sed -i "/^[^!]*\<Call\> *\<Message\> *(/I {
=
}" $f90_filename
sed -i "/^[^!]*\<Write\> *(6, *\*) *'/I {
=
}" $f90_filename

# Incorporate these identifiers into named anchors
sed -i '/^[0123456789]\+$/ {
s/^\([0123456789]\+\)$/<a name="ErrorMsgId\1">/
N
s/\n//
}' $f90_filename

# Create copy of the Fortran for local processing
cp $f90_filename $f90_filename.temp

# Strip out all comments
sed -i '/"!"[^!]*$/! s/![^"].*$//g' $f90_filename.temp

# Strip out any Fortran content before the 'Call Message' statement
sed -i "s/\(^<a name=${quote}.*${quote}> *\).*\( *\<Call\> *\<Message\>.*\)/\1\2/g" $f90_filename.temp

# Extract all error messages to temporary files
sed -i "/^[^!]*\<Call\> *\<Message\> *(/I {
  : repeat
  s/\(^.*[^& ] *$\)/\1/
  t break
  N
  b repeat
  : break
  s/ *& *\n/ /g
  s/ \+/ /g
  s/' *\/\/ *'//g
  s/ *$//g
  w messages.temp
}" $f90_filename.temp
sed -i "/^[^!]*\<Write\> *(6, *\*) *'/I {
  : repeat
  s/\(^.*[^& ] *$\)/\1/
  t break
  N
  b repeat
  : break
  s/ *& *\n/ /g
  s/ \+/ /g
  s/' *\/\/ *'//g
  s/ *$//g
  w messages_misc.temp
}" $f90_filename.temp

# Extract each category of error message to a separate file and set up anchors
for i in 1 2 3 4
do
sed -i "/.*${i} *)$/ {
 s/^<a name=${quote}\(.*\)${quote}> *\<Call\> *\<Message\> *( *\(.*\), *${i} *)$/<li><a href=${quote}${f90_filename}.html#\1${quote}>\2<\/a>/I
 w messages_${i}.temp
}" messages.temp
echo 'Message code: '${i} >> errors.html
echo '<ul>' >> errors.html
cat messages_${i}.temp >> errors.html
echo '</ul>' >> errors.html
done
sed -i "/.*[')] *)$/ {
 s/^<a name=${quote}\(.*\)${quote}> *\<Call\> *\<Message\> *( *\(.*\) *)$/<li><a href=${quote}${f90_filename}.html#\1${quote}>\2<\/a>/I
 w messages_0.temp
}" messages.temp
echo 'Message code: null' >> errors.html
echo '<ul>' >> errors.html
cat messages_0.temp >> errors.html
echo '</ul>' >> errors.html
sed -i "{
 s/^<a name=${quote}\(.*\)${quote}> *\<Write\> *(6, *\*) *\(.*\)$/<li><a href=${quote}${f90_filename}.html#\1${quote}>\2<\/a>/I
}" messages_misc.temp
echo 'Message code: "write (6, *)"' >> errors.html
echo '<ul>' >> errors.html
cat messages_misc.temp >> errors.html
echo '</ul>' >> errors.html
echo '<br><br>' >> errors.html
echo '' >> errors.html

# Check that all error messages in messages.temp have been processed
grep "<a name=" messages.temp

# remove temp files
rm *.temp

done

# Create closing section of the listing
echo '</pre></body></html>' >> errors.html

exit 0


