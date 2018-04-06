pro writefield,datadir,seltime,selgrid,selfield,sellevel,$
    selspecies,lon2d,lat2d,fielddata,fieldhead,xyrel=xyrel,$
    modtitle=modtitle,namever=namever,timeframe=timeframe,$
    TimeAveragingList=TimeAveragingList,nrecs=nrecs,$
    seldatavalue=seldatavalue,selname=selname,AAR=AAR,filename=filename


;--------------------------------------------------------------
; procedure to write field output from name
; reads model header, field header and field
; overwrites single field to same file
;
; LH 06/05/2010
; adapted from readfield.pro
; DBR 20/02/2002
;
;--------------------------------------------------------------

; inputs:
;  datadir            : REFER TO PLOTFIELD INPUTS
;  seltime            :
;  selgrid            :
;  selfield           :
;  sellevel           :
;  selspecies         :
;  namever            :
;  timeframe          :
;  TimeAveragingList  :
;  seldatavalue       :
;  selname            :
;  AAR                : 
;  filename           : give user defined filename instead of name style filename

; outputs:
;  lon2D              : array of longitudes
;  lat2D              : array of latitudes
;  fielddata          : data array
;  fieldhead          : array of column headers
;  xyrel              : source location

;--------------------------------------------------------------
if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(seldatavalue) eq 0)then seldatavalue=' '
if(n_elements(selname) eq 0)then selname=' '
if(n_elements(AAR) eq 0) then AAR=0

; generate filename

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

if (n_elements(filename) eq 0) then begin
  if (namever eq 3)then begin
    ; NAME III
    filename=datadir+delim+selgrid+'_'+seltime+'.txt'
  endif else begin
    ;NAME
    filename=datadir+delim+'Fields_'+selgrid+'_'+seltime+'.txt'
  endelse
endif

;-----------------------------------------------------------------------  
; read title text

status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'0. Cannot open file:',filename
  return
endif

readtitle,filename,modtitle,modhead1,modhead2,namever=namever

;-----------------------------------------------------------------------  
; extract grid details

xo=modhead2(where(modhead1 eq 'X grid origin:'))
yo=modhead2(where(modhead1 eq 'Y grid origin:'))
nx=modhead2(where(modhead1 eq 'X grid size:'))
ny=modhead2(where(modhead1 eq 'Y grid size:'))
dx=modhead2(where(modhead1 eq 'X grid resolution:'))
dy=modhead2(where(modhead1 eq 'Y grid resolution:'))


if (namever eq 3)then begin
  ; NAME III
  nfields=modhead2(where(modhead1 eq 'Number of field cols:'))
endif else begin
  ;NAME
  nfields=modhead2(where(modhead1 eq 'Number of fields:'))
endelse


xo=float(xo(0))
yo=float(yo(0))
nx=fix(nx(0))
ny=fix(ny(0))
dx=float(dx(0))
dy=float(dy(0))

nfields=fix(nfields(0))


;-----------------------------------------------------------------------  
; define lon/lat arrays

lon2d=fltarr(nx,ny)
lat2d=fltarr(nx,ny)

if (namever eq 3)then begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      lon2d(ix,iy)=xo+float(ix)*dx
      lat2d(ix,iy)=yo+float(iy)*dy
    endfor
  endfor
endif else begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      lon2d(ix,iy)=xo+float(ix+0.5)*dx
      lat2d(ix,iy)=yo+float(iy+0.5)*dy
    endfor
  endfor
endelse

;-----------------------------------------------------------------------  
; ensure -180 to 180
; (undo any shift done by readfield)

index1=where(lon2d gt 180.0,nindex1)
if(nindex1 gt 0)then begin
  lon2d(index1)=lon2d(index1)-360.0
  index2=sort(lon2d(*,0))
  nshift=index2(0)-2
  lat2d     =shift(lat2d,nshift,0)
  lon2d     =shift(lon2d,nshift,0)
  nshift=-nshift
  fielddata=Shift(fielddata,nshift)  
endif 

;-----------------------------------------------------------------------  
; field headers

status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'1. Cannot open file:',filename
  fielddata=lon2d*0.0
  return
endif

if (namever eq 3) then begin
  ; NAME III
  readdata,filename,fieldhead,nfields=nfields+5,nrows=17,nskip=19
endif else begin
  ; NAME
  readdata,filename,fieldhead,nfields=nfields+5,nrows=7,nskip=17
endelse

fieldhead=strtrim(fieldhead,2)
fieldhead=fieldhead(4:*,*)

;-----------------------------------------------------------------------
; select field

 if (namever eq 3)then begin
   ; NAME III
   if (selname ne ' ') then begin
     ; field specified by its name
     selname=selname(0)
     index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selname,/remove_all),nindex)
   endif else begin
     ; field specified in usual way
     index=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selspecies,/remove_all) and $
                 strlowcase(strcompress(fieldhead(*,2),/remove_all)) eq strlowcase(strcompress(selfield,/remove_all)) and $
	         strcompress(fieldhead(*,14),/remove_all) eq strcompress(sellevel,/remove_all) and $
	         strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList,/remove_all) and $
	         strcompress(fieldhead(*,15),/remove_all) eq strcompress(seldatavalue,/remove_all),nindex)
   endelse
                                    
 endif else begin
   ; NAME
   index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selspecies,/remove_all) and $
               strlowcase(strcompress(fieldhead(*,3),/remove_all)) eq strlowcase(strcompress(selfield,/remove_all)) and $
	       strcompress(fieldhead(*,5),/remove_all) eq strcompress(sellevel,/remove_all),nindex)
 endelse

;-----------------------------------------------------------------------  
; read data

if(nindex gt 0)then begin

 if (namever eq 3) then begin
   ; NAME III
   readdata,filename,data,nrecs,nfields=nfields+5,nskip=36
 endif else begin
   ; NAME
   readdata,filename,data,nrecs,nfields=nfields+5,nskip=24
 endelse

 if (nrecs eq 0) then begin
   data=0
 endif else begin
    data(index(0)+4,*)=string(fielddata,format='(E24.8)')
 endelse
 
 fieldhead=reform(fieldhead(index(0),*))

endif else begin

  print,'No field found'

  ; Check for which field

  if (namever eq 3)then begin
    ; NAME III
    if (selname ne ' ') then begin
     ; field specified by its name
     index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selname,/remove_all),count)
     if (count eq 0) then print,'Cannot match name: ',selname
    endif else begin
      ; field specified in usual way
      index=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selspecies,/remove_all),count1)
      if (count1 eq 0) then print,'Cannot match species: ',selspecies
      index=where(strlowcase(strcompress(fieldhead(*,2),/remove_all)) eq strlowcase(strcompress(selfield,/remove_all)),count2)
      if (count2 eq 0) then print,'Cannot match field: ',selfield
      index=where(strcompress(fieldhead(*,14),/remove_all) eq strcompress(sellevel,/remove_all),count3)
      if (count3 eq 0) then print,'Cannot match level: ',sellevel
      index=where(strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList,/remove_all),count4)
      if (count4 eq 0) then print,'Cannot match time average: ',TimeAveragingList
      index=where(strcompress(fieldhead(*,15),/remove_all) eq strcompress(seldatavalue,/remove_all),count5)
      if (count4 eq 0) then print,'Cannot match data value: ',seldatavalue
    endelse

  endif else begin
    ; NAME
   index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selspecies,/remove_all),count1)
   if (count1 eq 0) then print,'Cannot match species: ',selspecies
   index=where(strlowcase(strcompress(fieldhead(*,3),/remove_all)) eq strlowcase(strcompress(selfield,/remove_all)),count2)
   if (count2 eq 0) then print,'Cannot match field: ',selfield
   index=where(strcompress(fieldhead(*,5),/remove_all) eq strcompress(sellevel,/remove_all),count3)
   if (count3 eq 0) then print,'Cannot match level: ',sellevel
  endelse

  status=-1

endelse

ncols=n_elements(data(*,0))
if (n_elements(data) eq 1) then begin
  print, 'No data in file ',filename
  return
endif 
if(data(ncols-1,0) eq '') then begin
  data=data(0:ncols-2,*)
  ncols=ncols-1
endif

get_lun,u
openu,u,filename

if (namever eq 3) then begin
; NAME III
  skip_lun,u,38,/Lines
endif else begin
; NAME
  skip_lun,u,26,/Lines
endelse

truncate_lun,u
for iline=0L,nrecs-1 do begin
  for item=0,3 do begin
    printf,u,data(item,iline),',',format='(a19,a1,$)'
  endfor
  if (ncols-4) gt 1 then begin
    for item=4,ncols-2 do begin
      printf,u,data(item,iline),',',format='(a24,a1,$)'
    endfor
  endif
  printf,u,data(ncols-1,iline),',',format='(a24,a1)'
endfor
close,u
free_lun,u

;-----------------------------------------------------------------------  

end
