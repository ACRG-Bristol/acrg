pro genidlplots, GRAPHICSDIR, WORKDIR

;##################################################################
; Sample procedure for plotting NAME III fields output with IDL.
; Andrew Jones, Atmospheric Dispersion, UK Met Office, 31/01/2012
;
; Arguments
;  GRAPHICSDIR : (string)  NAME graphics directory (i.e. directory
;                           containing the NAME IDL routines).
;  WORKDIR     : (string)  NAME output directory (i.e. directory
;                           containing the NAME output files).
;##################################################################

if (n_elements(GRAPHICSDIR) eq 0) then begin
  print, 'ERROR: argument GRAPHICSDIR is required here'
  stop
endif

if (n_elements(WORKDIR) eq 0) then begin
  print, 'ERROR: argument WORKDIR is required here'
  stop
endif

;
; Plot NAME output fields
;

print," - generating output graphics in " + WORKDIR

cd, GRAPHICSDIR

; Plot 1: 3-hr average air concentration field over 0-100 m
selfield=['Air Concentration']
sellevel=['Z = 50.00000 m agl']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['3hr 0min average']
plotfield,WORKDIR,selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,$
 multiplot=[0,1,1],plotheader=1,plotlegend=1,group='field',$
 SelTimeAverage=seltav,filetext='Av_0m_100m'

; Plot 2: 3-hr average air concentration field over boundary-layer
selfield=['Air Concentration']
sellevel=['Boundary layer average']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['3hr 0min average']
plotfield,WORKDIR,selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,$
 multiplot=[0,1,1],plotheader=1,plotlegend=1,group='field',$
 SelTimeAverage=seltav,filetext='Av_BoundaryLayer'

; Plot 3: time-integrated boundary-layer air concentration and deposition fields
selfield=['Air Concentration','Dry deposition','Wet deposition','Deposition']
sellevel=['Boundary layer average','','','']
selspecies='DepositingTracer'
selgrid='Fields_grid1'
seltav=['1day 12hr 0min integral','1day 12hr 0min integral','1day 12hr 0min integral','1day 12hr 0min integral']
plotfield,WORKDIR,selgrid,selspecies=selspecies,namever=3,exact=1,$
 selfield=selfield,sellevel=sellevel,plotpmsl=0,$
 multiplot=[0,2,2],plotheader=1,plotlegend=1,group='field',$
 SelTimeAverage=seltav,filetext='DosageAndDep'

print,'Plotting finished'

end
