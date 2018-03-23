pro readvertfield,datadir,seltime,selgrid,selfield,sellevel,$
    selspecies,lon2d,lat2d,fielddata,fieldhead,xyrel=xyrel,$
    modtitle=modtitle,namever=namever,timeframe=timeframe,$
    TimeAveragingList=TimeAveragingList,nrecs=nrecs,$
    seldatavalue=seldatavalue,selname=selname,zrel=zrel,$
    vgrid=vgrid


;--------------------------------------------------------------
; procedure to read field output from name on z-x or z-y grid
; reads model header, field header and field
;
; SJL Aug 2009
; (created by modifying readfield.pro)
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
;  vgrid              :

; outputs:
;  lon2D              : array of longitudes or latitudes
;  lat2D              : array of heights
;  fielddata          : data array
;  fieldhead          : array of column headers
;  xyrel              : source location

;--------------------------------------------------------------
if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(seldatavalue) eq 0)then seldatavalue=' '
if(n_elements(selname) eq 0)then selname=' '
if(n_elements(AAR) eq 0) then AAR=0
if(n_elements(vgrid) eq 0) then vgrid=[0.0,0.0,0.0]

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

;-----------------------------------------------------------------------  
; read title text

status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'0. Cannot open file:',filename
  return
endif

readtitle,filename,modtitle,modhead1,modhead2,namever=namever

;-----------------------------------------------------------------------  
; extract lat or lon grid details

if strmatch(strmid(sellevel,0,1),'Y') then begin
    xo=modhead2(where(modhead1 eq 'X grid origin:'))
    nx=modhead2(where(modhead1 eq 'X grid size:'))
    dx=modhead2(where(modhead1 eq 'X grid resolution:'))
endif else if strmatch(strmid(sellevel,0,1),'X') then begin
    xo=modhead2(where(modhead1 eq 'Y grid origin:'))
    nx=modhead2(where(modhead1 eq 'Y grid size:'))
    dx=modhead2(where(modhead1 eq 'Y grid resolution:'))
endif else begin
    print, 'X-dimension is unclear.'
endelse

xo=float(xo(0))
nx=fix(nx(0))
dx=float(dx(0))

;----------------------------------------------------------------------
;determine number of fields

if (namever eq 3)then begin
  ; NAME III
  nfields=modhead2(where(modhead1 eq 'Number of field cols:'))
endif else begin
  ;NAME
  nfields=modhead2(where(modhead1 eq 'Number of fields:'))
endelse

nfields=fix(nfields(0))

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
        index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selname,/remove_all),nindex)
    endif else begin
     ; field specified in usual way
        index=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selspecies,/remove_all) and $
                    strlowcase(strcompress(fieldhead(*,2),/remove_all)) eq strlowcase(strcompress(selfield,/remove_all)) and $
                    strcompress(fieldhead(*,15),/remove_all) eq strcompress(sellevel,/remove_all) and $
                    strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList,/remove_all) and $
                    strcompress(fieldhead(*,16),/remove_all) eq strcompress(seldatavalue,/remove_all),nindex)
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
        readdata,filename,data,nrecs,nfields=nfields+5,nskip=37
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
            index=where(strcompress(fieldhead(*,15),/remove_all) eq strcompress(sellevel,/remove_all),count3)
            if (count3 eq 0) then print,'Cannot match level: ',sellevel
            index=where(strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList,/remove_all),count4)
            if (count4 eq 0) then print,'Cannot match time average: ',TimeAveragingList
            index=where(strcompress(fieldhead(*,16),/remove_all) eq strcompress(seldatavalue,/remove_all),count5)
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

;-----------------------------------------------------------------------  
; define vertical array

ny = fix(max(ygrid))

maxylat = max(ylat,smax)
ylat = ylat + (ylat eq -999.0)*(maxylat+999.0)
minylat = min(ylat,smin)
;
;
if (vgrid[0] lt 1) then begin
;Assume constant intervals to rebuild hgt array
    b2 = (maxylat-minylat)/(ygrid(smax)-ygrid(smin))
    a2 = maxylat - b2*ygrid(smax)
    zz = (findgen(ny)+1.0)*b2 + a2
endif else begin
    ny = vgrid[0]
    zz = findgen(ny)*vgrid[1] + vgrid[2]
endelse
;
;Expand lat and lon arrays to 2-dimensions
lat2d = fltarr(nx,ny)

for ii = 0, nx-1 do begin
    lat2d[ii,*] = zz
endfor

;-----------------------------------------------------------------------  
; define lon/lat array

lon2d=fltarr(nx,ny)

if (namever eq 3)then begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      lon2d(ix,iy)=xo+float(ix)*dx
    endfor
  endfor
endif else begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      lon2d(ix,iy)=xo+float(ix+0.5)*dx
    endfor
  endfor
endelse

;
;-----------------------------------------------------------------------
;Extract data

fielddata=fltarr(nx,ny)
;
if(status eq 1)then begin
    for i=0L,n_elements(ygrid)-1 do begin
       ix=fix(xgrid(i))
       iy=fix(ygrid(i))
       fielddata(ix-1,iy-1)=float(data(i))
    endfor
endif

;-----------------------------------------------------------------------  
; release location

cloc=strsplit(modhead2(where(modhead1 eq 'Release location:')),' ',/extract)
len=strlen(cloc)

dir=strarr(2)
dir(0)=strmid(cloc(0),len(0)-1)
dir(1)=strmid(cloc(1),len(1)-1)

if(dir(0) eq 'W' OR dir(0) eq 'E')then begin

  xyrel=dblarr(2)
  xyrel(0)=double(strmid(cloc(0),0,len(0)-2))
  xyrel(1)=double(strmid(cloc(1),0,len(1)-2))

  if(dir(0) eq 'W')then begin
    xyrel(0)=-xyrel(0)
  endif

  if(dir(1) eq 'S')then begin
    xyrel(1)=-xyrel(1)
  endif

endif else begin
  xyrel=[-999.0,-999.0]
endelse

;release height
chgt=strsplit(modhead2(where(modhead1 eq 'Release height:')),' ',/extract)
zrel=double(strmid(chgt(0),0,len(0)-2))

;-----------------------------------------------------------------------  

end
