#!/bin/sh

set -e

# create tarfiles.txt with the list of parameter files
./create_def.pl > /dev/null

# create the tarball with these parameter files
version=`grep parametersVersion parameters_version.def | awk 'BEGIN {FS="=";}{ sub(/ *;$/,"");print $2; }'`
cat tarfiles.txt | xargs tar zcf grib_api_parameters-v$version.tar.gz

# get the current html download page for GRIB API
dlpage=grib_api.html
rm -f $dlpage || true
cadaver http://wedit.ecmwf.int:81/products/data/software/download << EOF
get $dlpage
EOF

./inject_download_page.pl $dlpage > temp.html
mv temp.html $dlpage

# upload the updated download page
cadaver http://wedit.ecmwf.int:81/products/data/software/download << EOF
put $dlpage
cd software_files
put grib_api_parameters-v$version.tar.gz
EOF
