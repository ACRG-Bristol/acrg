pro readtshead,datadir,modtitle,tshead,selgrid=selgrid,$
               namever=namever,timeframe=timeframe,selcase=selcase

;--------------------------------------------------------------
; procedure to read timeseries header
;
; DBR 20/02/2002
;
; Arguments
;
; required:
;  datadir      : (string)  run directory
;
; output:
;  modtitle     : (string array) model info/header
;  tshead       : (string array) time series header 
;--------------------------------------------------------------

if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(selcase) eq 0)then selcase='1'

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse


if (namever eq 3)then begin
  fieldfiles=findfile(datadir+delim+selgrid+'*')
endif else begin
  ;NAME
  fieldfiles=findfile(datadir+delim+'Time_series*')
endelse

filename=fieldfiles(0)
numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

modhead1=strarr(1)
modhead2=strarr(1)

if (namever eq 3)then begin
  readtitle,filename,modtitle,modhead1,modhead2,namever=namever
endif else begin
  ;NAME
  readtitle,filename,modtitle,modhead1,modhead2,namever=namever,nrecs=13
endelse

; number of series

if (namever eq 3)then begin
  nfields=modhead2(where(modhead1 eq 'Number of field cols:'))
endif else begin
  ;NAME
  nfields=modhead2(where(modhead1 eq 'Number of series:'))
endelse

nfields=fix(nfields(0))

; series headers

if (namever eq 3)then begin

  readdata,filename,fieldhead,nskip=19,nrows=19

endif else begin
  ;NAME
 
  readdata,filename,fieldhead,nskip=13,nrows=8

endelse

  
fieldhead=strtrim(fieldhead,2)
tshead=fieldhead(1:*,*)

end
