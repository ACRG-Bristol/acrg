pro readplume,filename,modtitle,xlon,ylat,zqfe,zqnh,age,$
  slon,slat,sqfe,blflag,xyrel=xyrel


;--------------------------------------------------------------
; procedure to read plume
; reads model header and all fields
;
; DBR 22/03/2002
;--------------------------------------------------------------


;-----------------------------------------------------------------------  
; read title text

modhead1=strarr(1)
modhead2=strarr(1)

status=dc_read_fixed(filename,modhead1,modhead2,/column,nrecs=11,nskip=0,$
  resize=[1,2],Format="(A21,A40)")
  
if(status ne 0)then begin
 print,'Cannot open file:',filename
 return
endif  

modtitle=modhead1+modhead2 

modhead1=strtrim(modhead1,2)
modhead2=strtrim(modhead2,2)
modtitle=strtrim(modtitle,2)


;-----------------------------------------------------------------------  
; read data

status=dc_read_free(filename,xlon,ylat,zqfe,zqnh,age,slon,slat,sqfe,$
 /column,resize=[1,2,3,4,5,6,7,8],nskip=14)

if(status ne 0)then begin
 print,'Cannot open file:',filename
 return
endif

; Negative zqfe heights indicate that the particle is in the boundary layer
; Create Boundary layer flag and negate negative zqfe heights

blflag=zqfe*0.0
indzero=where(zqfe lt 0.0,nindzero)
if nindzero gt 0 then begin
  blflag(indzero)=1.0
  zqfe(indzero)=-zqfe(indzero)
endif

;-----------------------------------------------------------------------  
; release location


cloc=modtitle(where(modhead1 eq 'Release location:'))
cloc=cloc(0)

d1=strpos(cloc,'.')
if(d1 ge 0)then begin
  d2=strpos(cloc,'.',d1+1)
endif else begin
  d2=-1
endelse

if(d1 ge 0 and d2 ge 0)then begin

  if(strmid(cloc,d1+5,1) eq 'N' or strmid(cloc,d1+5,1) eq 'S')then begin
    dn=d1
    de=d2
  endif else begin
    dn=d2
    de=d1
  endelse
  
  xyrel=fltarr(2)

  xyrel(0)=float(strmid(cloc,de-3,8))
  if(strmid(cloc,de+5,1) eq 'W')then begin
    xyrel(0)=-xyrel(0)
  endif 
  
  xyrel(1)=float(strmid(cloc,dn-2,7))
  if(strmid(cloc,dn+5,1) eq 'S')then begin
    xyrel(1)=-xyrel(1)
  endif   
  
endif else begin
  xyrel=[-999.0,-999.0]
endelse

end
