#!/bin/sh
# Copyright 2005-2007 ECMWF
# 
# Licensed under the GNU Lesser General Public License which
# incorporates the terms and conditions of version 3 of the GNU
# General Public License.
# See LICENSE and gpl-3.0.txt for details.

. ./include.sh

cat > list.filter <<EOF
set productDefinitionTemplateNumber=31;
set numberOfContributingSpectralBands=4;
set satelliteSeries={61,61,61,62};
set satelliteNumber={4,4,4,5};
set instrumentType={591,591,591,592};
set scaleFactorOfCentralWaveNumber={3,3,3,4};
set scaledValueOfCentralWaveNumber={26870,9272,8377,22334};
write "test.grib";
EOF

file="${data_dir}/reduced_gaussian_model_level.grib2"

${tools_dir}grib_filter  list.filter $file

cat >list.filter<<EOF
print("x.out") "productDefinitionTemplateNumber=[productDefinitionTemplateNumber]";
print("x.out") "numberOfContributingSpectralBands=[numberOfContributingSpectralBands]";
print("x.out") "satelliteSeries={[satelliteSeries',']}";
print("x.out") "satelliteNumber={[satelliteNumber',']}";
print("x.out") "instrumentType={[instrumentType',']}";
print("x.out") "scaleFactorOfCentralWaveNumber={[scaleFactorOfCentralWaveNumber',']}";
print("x.out") "scaledValueOfCentralWaveNumber={[scaledValueOfCentralWaveNumber',']}";
EOF

${tools_dir}grib_filter list.filter test.grib

cat >list.ref<<EOF
productDefinitionTemplateNumber=31
numberOfContributingSpectralBands=4
satelliteSeries={61,61,61,62}
satelliteNumber={4,4,4,5}
instrumentType={591,591,591,592}
scaleFactorOfCentralWaveNumber={3,3,3,4}
scaledValueOfCentralWaveNumber={26870,9272,8377,22334}
EOF

diff x.out list.ref

rm -f x.out list.ref list.filter test.grib

