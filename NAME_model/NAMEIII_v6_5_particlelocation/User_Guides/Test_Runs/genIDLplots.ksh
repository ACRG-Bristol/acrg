#!/bin/ksh
##################################################################
# Sample script for plotting NAME III fields output with IDL.    #
# Andrew Jones, Atmospheric Dispersion, UK Met Office 31/01/2012 #
##################################################################

# Set the NAME III graphics directory.

GRAPHICSDIR="../../Code_IDLGraphics"

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

idl <<DELIMIT

cd,'${GRAPHICSDIR}'

; Plot 1: 3-hr average air concentration field over 0-100 m
selfield=['Air Concentration']
sellevel=['Z = 50.00000 m agl']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['3hr 0min average']
plotfield,'${WORKDIR}',selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,gengif=1,genjpg=0,$
 multiplot=[0,1,1],plotheader=1,plotlegend=1,genanim=1,group='field',$
 SelTimeAverage=seltav,filetext='Av_0m_100m'

; Plot 2: 3-hr average air concentration field over boundary-layer
selfield=['Air Concentration']
sellevel=['Boundary layer average']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['3hr 0min average']
plotfield,'${WORKDIR}',selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,gengif=1,genjpg=0,$
 multiplot=[0,1,1],plotheader=1,plotlegend=1,genanim=1,group='field',$
 SelTimeAverage=seltav,filetext='Av_BoundaryLayer'

; Plot 3: time-integrated boundary-layer air concentration and deposition fields
selfield=['Air Concentration','Dry deposition','Wet deposition','Deposition']
sellevel=['Boundary layer average','','','']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['1day 12hr 0min integral','1day 12hr 0min integral','1day 12hr 0min integral','1day 12hr 0min integral']
plotfield,'${WORKDIR}',selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,gengif=1,genjpg=0,$
 multiplot=[0,2,2],plotheader=1,plotlegend=1,genanim=1,group='field',$
 SelTimeAverage=seltav,filetext='DosageAndDep'

DELIMIT

exit 0
