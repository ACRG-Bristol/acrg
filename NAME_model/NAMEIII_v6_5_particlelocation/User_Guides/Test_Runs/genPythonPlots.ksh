#!/bin/ksh
##################################################################
# Sample script for plotting NAME III fields output with python. #
# Andrew Jones, Atmospheric Dispersion, UK Met Office 31/01/2012 #
##################################################################

# Set the Python scripts directory.

PYTHONDIR="./Code_Python"
cd ${PYTHONDIR}

#
# Single argument specifies the NAME output directory
#

if [ $# -eq 1 ] ; then
  
  WORKDIR=$1
  
else
  
  echo "Error in $0 - a single argument WORKDIR is required here"
  exit 1
  
fi

#
# Plot NAME output fields
#

echo " - generating output graphics in ${WORKDIR}"

# Plot 1: 3-hr average air concentration field over 0-100 m
python PlotField2d.py --datadir=${WORKDIR}  --plotdir=${WORKDIR}  --grid="grid1"  --fieldname="0-100m Av Air Conc"  --label="AirConc_0m_100m"

convert -delay 25 ${WORKDIR}AirConc_0m_100m_T*.png ${WORKDIR}Animate_AirConc_0m_100m.gif

# Plot 2: 3-hr average air concentration field over boundary-layer
python PlotField2d.py --datadir=${WORKDIR}  --plotdir=${WORKDIR}  --grid="grid1"  --fieldname="BL Av Air Conc"      --label="AirConc_BoundaryLayer"

convert -delay 25 ${WORKDIR}AirConc_BoundaryLayer_T*.png ${WORKDIR}Animate_AirConc_BoundaryLayer.gif

# Plot 3: time-integrated boundary-layer air concentration and deposition fields
python PlotField2d.py --datadir=${WORKDIR}  --plotdir=${WORKDIR}  --grid="grid1"  --fieldname="Total Deposition"    --label="TotalDeposition"

convert -delay 25 ${WORKDIR}TotalDeposition_T*.png ${WORKDIR}Animate_TotalDeposition.gif

exit 0
