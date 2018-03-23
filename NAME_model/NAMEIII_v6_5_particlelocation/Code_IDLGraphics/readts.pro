pro readts,datadir,tssplit,location,field,level,species,grid, $
  date,seldata,units,namever=namever,timeframe=timeframe,$
  selcase=selcase,found=found

;--------------------------------------------------------------
; procedure to read field output from name
; reads model header, time series headers and all series
;
; DBR 20/02/2002
;
; Arguments
;
; required:
;  datadir      : (string)  run directory
;  tssplit      : (string) 'LOCATION','SPECIES' or 'NONE'
;                 from namelist variable TSSPLIT - defines 
;                 how time series files were named
; location      : (string) time series location (eg 'ASCOT')
; field         : (string) time series field (eg 'Air concentration')
; species       : (string) time series species (eg 'TRACER')
; grid          : (string) 'grid1' or 'grid2' - time series grid
;
;output:
; date         : (date array) time series date/time
; seldata      : (fltarr)  time series data
; units        : (string) time series units 
;--------------------------------------------------------------

if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(selcase) eq 0)then selcase='1'
found=0 ;Assume series is not found unless set otherwise

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

;-----------------------------------------------------------------------
; define date/time structure

date={datetime,year:0,month:0,day:0,hour:0,minute:0,second:0} 

;-----------------------------------------------------------------------
; filename

if (namever eq 3)then begin
  ;NAME III
   case tssplit of 

    'LOCATION': begin
                filename=datadir+delim+grid+'_C'+selcase+'_'+location+'.txt'	    
                end   
    'SPECIES':  begin
                filename=datadir+delim+grid+'_C'+selcase+'_'+species+'.txt'		
                end   
    'NONE':     begin
                filename=datadir+delim+grid+'_C'+selcase+'.txt'
		end   
  endcase

endif else begin
  ;NAME
  case tssplit of 

    'LOCATION': begin
                filename=datadir+delim+'Time_series_'+grid+'_'+location+'.txt'	    
                end   
    'SPECIES':  begin
                filename=datadir+delim+'Time_series_'+grid+'_'+species+'.txt'		
                end   
    'NONE':     begin
                filename=datadir+delim+'Time_series_'+grid+'.txt'
		end   
  endcase
endelse

filename=strcompress(filename,/remove_all)

; read title text

modhead1=strarr(1)
modhead2=strarr(1)

if (namever eq 3)then begin
  readtitle,filename,modtitle,modhead1,modhead2,namever=namever
endif else begin
  ;NAME
  readtitle,filename,modtitle,modhead1,modhead2,namever=namever,nrecs=14
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
fieldhead=fieldhead(1:*,*)

; read data

if (namever eq 3)then begin
  readdata,filename,data,nrecs,nskip=38
endif else begin
  ;NAME
  readdata,filename,data,nrecs,nskip=21
endelse

date=replicate({datetime},nrecs)

for i=long(0),nrecs-1 do begin
  if (namever eq 3) then begin
    tempdat=strsplit(data(0,i),' ',/Extract,/Preserve_Null)
    temptime=strsplit(tempdat(1),':',/Extract,/Preserve_Null)
    tempdat=strsplit(tempdat(0),'/',/Extract,/Preserve_Null)
  endif else begin
    ;NAME
    tempdat=strsplit(data(0,i),'/',/Extract,/Preserve_Null)
    temptime=strsplit(data(1,i),':',/Extract,/Preserve_Null)
  endelse
  date(i).year=tempdat(2)
  date(i).month=tempdat(1)
  date(i).day=tempdat(0)
  date(i).hour=temptime(0)
  date(i).minute=temptime(1)
  if n_elements(temptime) eq 3 then begin
    date(i).second=temptime(2)
  endif 
endfor

if (namever eq 3) then begin
  data=float(data(1:*,*))
endif else begin
  ;NAME
  data=float(data(2:*,*))
endelse

; select field

if (namever eq 3)then begin

 index=where(strlowcase(strcompress(fieldhead(*,13),/remove_all)) $
                          eq strlowcase(strcompress(location,/remove_all)) and $
             strlowcase(strcompress(fieldhead(*,3),/remove_all)) $
                          eq strlowcase(strcompress(species,/remove_all)) and $
	     strlowcase(strcompress(fieldhead(*,2),/remove_all)) $
                          eq strlowcase(strcompress(field,/remove_all)) and $
	     strlowcase(strcompress(fieldhead(*,16),/remove_all)) $
                          eq strlowcase(strcompress(level,/remove_all)),nindex)
                         
 endif else begin

 index=where(strcompress(fieldhead(*,2),/remove_all) eq strcompress(location,/remove_all) and $
             strcompress(fieldhead(*,4),/remove_all) eq strcompress(species,/remove_all) and $
	     strcompress(fieldhead(*,5),/remove_all) eq strcompress(field,/remove_all) and $
	     strcompress(fieldhead(*,6),/remove_all) eq strcompress(level,/remove_all),nindex)

endelse


if (nindex gt 0)then begin
  found=1
  if (namever eq 3)then begin
    seldata=reform(data(index(0),*))
    units=fieldhead(index(0),4)
  endif else begin
    seldata=reform(data(index(0),*))
    units=fieldhead(index(0),7)
  endelse
endif else begin
  print,'Selected series not found'
  print,location
  print,species
  print,field
  print,level
  found=0
endelse

end
