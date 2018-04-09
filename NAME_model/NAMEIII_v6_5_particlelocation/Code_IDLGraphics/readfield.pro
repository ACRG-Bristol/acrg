pro readfield,datadir,seltime,selgrid,selfield,sellevel,$
    selspecies,lon2d,lat2d,fielddata,fieldhead,xyrel=xyrel,$
    modtitle=modtitle,namever=namever,timeframe=timeframe,$
    TimeAveragingList=TimeAveragingList,nrecs=nrecs,$
    seldatavalue=seldatavalue,selname=selname,AAR=AAR,filename=filename, $
    nolongitudeshift=nolongitudeshift


;--------------------------------------------------------------
; procedure to read field output from name
; reads model header, field header and field
;
; DBR 20/02/2002
;
; Eike Mueller, 11/11/2010
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
;  nolongitudeshift   : if this (optional) parameter is set to 1, the 
;                       data-, lat- and longitude arrays are not shifted so
;                       that the longitude lies in the range [-180,180].

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
if(n_elements(nolongitudeshift) eq 0) then nolongitudeshift=0

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

index1=where(lon2d gt 180.0,nindex1)
if ( (nindex1 gt 0) and (nolongitudeshift eq 0) )then begin
  lon2d(index1)=lon2d(index1)-360.0
  index2=sort(lon2d(*,0))
  nshift=index2(0)-2
  lat2d     =shift(lat2d,nshift,0)
  lon2d     =shift(lon2d,nshift,0)  
endif else begin
  nshift=0
endelse


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
   ; --- NAME III ---
   if (selname ne ' ') then begin
     ; field specified by its name
     nfieldspecs = n_elements(selname)
     index = [-1]
     for i=0,nfieldspecs-1 do begin
       tmpindex=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selname(i),/remove_all),nindex)
       if (nindex ne 0) then index = [index,tmpindex(0)]
     endfor
     ; count number of elements and remove first element if anything has been found
     nindex = n_elements(index)-1
     if (nindex gt 0) then index = index(1:nindex)
   endif else begin
     ; Check that all arrays have the same size
     nfieldspecs = n_elements(selspecies)
     if ( not ( (n_elements(selfield)          eq nfieldspecs ) and $
                (n_elements(sellevel)          eq nfieldspecs ) and $
                (n_elements(TimeAveragingList) eq nfieldspecs ) and $
                (n_elements(seldatavalue)      eq nfieldspecs ) ) ) then begin
       print,'length of field specification lists selspecies, selfield, sellevel, TimeAveragingList and seldatavalue have to be identical'
       status = -1
     endif
     index = [-1]
     ; field specified in usual way
     ; loop over all array elements
     for i=0,nfieldspecs-1 do begin
       tmpindex=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selspecies(i),/remove_all) and $
                      strlowcase(strcompress(fieldhead(*,2),/remove_all)) eq strlowcase(strcompress(selfield(i),/remove_all)) and $
   	                  strcompress(fieldhead(*,14),/remove_all) eq strcompress(sellevel(i),/remove_all) and $
	                    strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList(i),/remove_all) and $
	                    strcompress(fieldhead(*,15),/remove_all) eq strcompress(seldatavalue(i),/remove_all),nindex)
       if (nindex ne 0) then index = [index,tmpindex(0)]
     endfor     
     ; count number of elements and remove first element if anything has been found
     nindex = n_elements(index)-1
     if (nindex gt 0) then index = index(1:nindex)
   endelse
                                    
 endif else begin
   ; --- NAME ---
   ; Verify that all field specification arrays have the same size
   nfieldspecs = n_elements(selspecies)
   if ( not ( (n_elements(selfield) eq nfieldspecs ) and $
              (n_elements(sellevel) eq nfieldspecs ) ) ) then begin
     print,'length of field specification lists selspecies, selfield and sellevel have to be identical'
     status = -1
   endif
   index = [-1]
   ; loop over all array elements
   for i=0,nfieldspecs-1 do begin
     tmpindex=where(strcompress(fieldhead(*,1),/remove_all)  eq strcompress(selspecies(i),/remove_all)           and $
         strlowcase(strcompress(fieldhead(*,3),/remove_all)) eq strlowcase(strcompress(selfield(i),/remove_all)) and $
	                  strcompress(fieldhead(*,5),/remove_all)  eq strcompress(sellevel(i),/remove_all),nindex)
     if (nindex ne 0) then index = [index,tmpindex(0)]
   endfor
   ; count number of elements and remove first element again if anything has been found
   nindex = n_elements(index)-1
   if (nindex gt 0) then index = index(1:nindex)
 endelse

;-----------------------------------------------------------------------  
; read data

if (nindex eq nfieldspecs) then begin

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

 if (nrecs eq 0) then begin
   data=0
 endif else begin
   temp=fltarr(nindex,nrecs)
   for i=0,nindex-1 do begin
     temp(i,*)=float(data(index(i)+4,*))
   endfor
   data=temp
 endelse
 
 fieldhead=reform(fieldhead(index,*))

endif else begin

  print,'Field(s) not found'

  ; Check for which field
  if (namever eq 3)then begin
    ; NAME III
    if (selname ne ' ') then begin
     for i=0,nfieldspecs-1 do begin
       ; field specified by its name
       index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selname(i),/remove_all),count)
       if (count eq 0) then print,'Cannot match name: ',selname(i)
     endfor
    endif else begin
      for i=0,nfieldspecs-1 do begin
        ; field specified in usual way
        index=where(strcompress(fieldhead(*,3),/remove_all) eq strcompress(selspecies(i),/remove_all),count1)
        if (count1 eq 0) then print,'Cannot match species: ',selspecies(i)
        index=where(strlowcase(strcompress(fieldhead(*,2),/remove_all)) eq strlowcase(strcompress(selfield(i),/remove_all)),count2)
        if (count2 eq 0) then print,'Cannot match field: ',selfield(i)
        index=where(strcompress(fieldhead(*,14),/remove_all) eq strcompress(sellevel(i),/remove_all),count3)
        if (count3 eq 0) then print,'Cannot match level: ',sellevel(i)
        index=where(strcompress(fieldhead(*,7),/remove_all) eq strcompress(TimeAveragingList(i),/remove_all),count4)
        if (count4 eq 0) then print,'Cannot match time average: ',TimeAveragingList(i)
        index=where(strcompress(fieldhead(*,15),/remove_all) eq strcompress(seldatavalue(i),/remove_all),count5)
        if (count4 eq 0) then print,'Cannot match data value: ',seldatavalue(i)
      endfor
    endelse

  endif else begin  
   ; NAME
   for i=0,nfieldspecs-1 do begin
     index=where(strcompress(fieldhead(*,1),/remove_all) eq strcompress(selspecies(i),/remove_all),count1)
     if (count1 eq 0) then print,'Cannot match species: ',selspecies(i)
     index=where(strlowcase(strcompress(fieldhead(*,3),/remove_all)) eq strlowcase(strcompress(selfield(i),/remove_all)),count2)
     if (count2 eq 0) then print,'Cannot match field: ',selfield(i)
     index=where(strcompress(fieldhead(*,5),/remove_all) eq strcompress(sellevel(i),/remove_all),count3)
     if (count3 eq 0) then print,'Cannot match level: ',sellevel(i)
   endfor
  endelse

  status=-1

endelse

;-----------------------------------------------------------------------  
; generate 2d array or array for Area At Risk

fielddata=fltarr(nfieldspecs,nx,ny)

if(status eq 1)then begin
  if(AAR eq 0)then begin
    ;NAME
    for i=0l,n_elements(xgrid)-1 do begin
      ix=fix(xgrid(i))
      iy=fix(ygrid(i))
      fielddata(*,ix-1,iy-1)=float(data(*,i))
    endfor
  endif else begin
  ;  Area At Risk
    fielddata=fltarr(nfieldspecs,n_elements(data),6)
    for i=0,nfieldspecs-1 do begin
      fielddata(i,*,0)=xgrid(*)
      fielddata(i,*,1)=ygrid(*)
      fielddata(i,*,2)=xlon(*)
      fielddata(i,*,3)=ylat(*)
      fielddata(i,*,4)=data(i,*)
    endfor
  endelse
endif

;-----------------------------------------------------------------------
; ensure -180 to 180

if(nindex1 gt 0)then begin  
  fielddata=Shift(fielddata,0,nshift,0)
endif 

; If first dimension has size 1, remove it for backward compatibility
fielddata_dims = size(fielddata,/dimensions)
if (fielddata_dims[0] eq 1) then begin
  fielddata = reform(fielddata,fielddata_dims[1],fielddata_dims[2])
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
  xyrel(0)=double(strmid(cloc(0),0,len(0)-1))
  xyrel(1)=double(strmid(cloc(1),0,len(1)-1))

  if(dir(0) eq 'W')then begin
    xyrel(0)=-xyrel(0)
  endif

  if(dir(1) eq 'S')then begin
    xyrel(1)=-xyrel(1)
  endif

endif else begin
  xyrel=[-999.0,-999.0]
endelse

;-----------------------------------------------------------------------  

end
