pro readmet,datadir,seltime,selgrid,selfield,lon2d,lat2d,fielddata,$
 modtitle,fieldhead,namever=namever

;--------------------------------------------------------------
; procedure to read met data
;
; DBR 09/04/2002

;--------------------------------------------------------------

if(n_elements(namever) eq 0)then namever=2

;--------------------------------------------------------------

; generate filename

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

if (namever eq 3)then begin
  ; NAME III
  filename=datadir+delim+selgrid+'_'+seltime+'.txt'
endif else begin
  ;NAME
  filename=datadir+delim+'Fields_'+selgrid+'_'+seltime+'.txt'
endelse
print,'Reading: ',filename,' ',selfield

;-----------------------------------------------------------------------  
; read title text

status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'0. Cannot open file:',filename
  return
endif

readtitle,filename,modtitle,modhead1,modhead2,namever=namever

if (namever eq 3)then begin
  ; NAME III
  nfields=modhead2(where(modhead1 eq 'Number of field cols:'))
endif else begin
  ;NAME
  nfields=modhead2(where(modhead1 eq 'Number of fields:'))
endelse

nfields=fix(nfields(0))

;-----------------------------------------------------------------------  
; extract grid details

xo=modhead2(where(modhead1 eq 'X grid origin:'))
yo=modhead2(where(modhead1 eq 'Y grid origin:'))
nx=modhead2(where(modhead1 eq 'X grid size:'))
ny=modhead2(where(modhead1 eq 'Y grid size:'))
dx=modhead2(where(modhead1 eq 'X grid resolution:'))
dy=modhead2(where(modhead1 eq 'Y grid resolution:'))

xo=float(xo(0))
yo=float(yo(0))
nx=fix(nx(0))
ny=fix(ny(0))
dx=float(dx(0))
dy=float(dy(0))


;-----------------------------------------------------------------------  
; define lon/lat arrays

lon2d=fltarr(nx,ny)
lat2d=fltarr(nx,ny)

for ix=0,nx-1 do begin
  for iy=0,ny-1 do begin
    lon2d(ix,iy)=xo+float(ix+0.5)*dx
    lat2d(ix,iy)=yo+float(iy+0.5)*dy
  endfor
endfor

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
  index=where(strcompress(fieldhead(*,2),/remove_all) eq strcompress(selfield,/remove_all),nindex)
endif else begin
  ; NAME
  index=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selfield,/remove_all),nindex)
endelse

;-----------------------------------------------------------------------  
; read data

if(nindex ge 0)then begin

 if (namever eq 3) then begin
   ; NAME III
   readdata,filename,data,nrecs,nfields=nfields+5,nskip=36
 endif else begin
   ; NAME
   readdata,filename,data,nrecs,nfields=nfields+5,nskip=24
 endelse
 
 data=strtrim(data,2)
 xgrid=float(data(0,*))
 ygrid=float(data(1,*))
 xlon=float(data(2,*))
 ylat=float(data(3,*))
 temp=fltarr(nrecs)
 temp(*)=float(data(index(0)+4,*))
 data=temp
 
 fieldhead=reform(fieldhead(index(0),*))

endif else begin

  print,'No field found'
 
  status=-1

endelse

;-----------------------------------------------------------------------  
; generate 2d array

fielddata=fltarr(nx,ny)

if(status eq 1)then begin

  for i=0l,n_elements(xgrid)-1 do begin
    ix=fix(xgrid(i))
    iy=fix(ygrid(i))

    fielddata(ix-1,iy-1)=float(data(i))

  endfor    

endif


;-----------------------------------------------------------------------  
; ensure -180 to 180

index1=where(lon2d gt 180.0,nindex)
if(nindex gt 0)then begin
  lon2d(index1)=lon2d(index1)-360.0
  index2=sort(lon2d(*,0))
  nshift=index2(0)-2
  fielddata=Shift(fielddata,nshift,0)
  lat2d     =shift(lat2d,nshift,0)
  lon2d     =shift(lon2d,nshift,0)  
endif 


end
