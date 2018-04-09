#!/bin/ksh
########################################################
# Script for plotting NAME III field output.           #
# Andrew Jones, Atmospheric Dispersion Group, 6/11/06. #
########################################################

# Set the output directory.
WORKDIR='/home/hc1400/aprj/NameIII/GroupFiles/IntroductoryExamples/Output'

# Set the NAME III graphics directory.
GRAPHICSDIR='/home/hc1400/aprj/NameIII/Code/NameIII/Code_PVWaveGraphics'

cd $WORKDIR
echo 'Processing output files in: ' $WORKDIR

wave <<DELIMIT

cd,'${GRAPHICSDIR}'

; Plot 1: time-averaged concentration of tracer
selfield=['Air Concentration']
sellevel=['Z = 0.0000000E+00 m agl']
selspecies='Tracer'
selgrid='Fields_grid1'
plotfield,'${WORKDIR}',selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,gengif=1,genjpg=0,$
 multiplot=[0,1,1],plotheader=1,plotlegend=1,genanim=1,group='field',$
 SelTimeAverage='No time averaging'
; SelTimeAverage='1hr 0min average'
; SelTimeAverage='1day 0hr 0min average'

DELIMIT

echo 'Finished generating graphics'
exit 0
