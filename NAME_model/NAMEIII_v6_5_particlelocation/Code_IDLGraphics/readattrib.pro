pro readattrib,datadir,grid,location,species,field,level, $
  matrixdates,lon2d,lat2d,matrixarr,head,modtitle,ntime

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
;  grid         : (string) 'grid1' or 'grid2' - time series grid
;  location     : (string) time series location (eg 'ASCOT')
;  field        : (string) time series field (eg 'Air concentration')
;  level        : (string) time series level (eg 'Boundary layer')
;  species      : (string) time series species (eg 'TRACER')
;
; output:
;  matrixdates  : (date array) dates 
;  lon2d        : (2-D fltarr) - longitude for each point
;  lat2d        : (2-D fltarr) - latitude for each point
;  matrixarr    : (3-D fltarr) - attribution matrix
;                 1st index is time  
;  head         : (string) header
;  modtitle     : (string array) model header
;  ntime        : (integer) number of timesteps    
;--------------------------------------------------------------

;-----------------------------------------------------------------------
; define date/time structure

date={datetime,year:-4713,month:1,day:1,hour:12,minute:0,second:0,julian:0.0D}

;--------------------------------------------------------------

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

; filename

filename=datadir+delim+'Attrib_'+grid+'_'+location+'.txt'	    
print,'Reading attribution: ',grid,' ',location,' ',$
 species,' ',field,' ',level

;-----------------------------------------------------------------------  
; read title text

readtitle,filename,modtitle,modhead1,modhead2

;-----------------------------------------------------------------------  
; extract grid details

xo=modhead2(where(modhead1 eq 'X grid origin:'))
yo=modhead2(where(modhead1 eq 'Y grid origin:'))
nx=modhead2(where(modhead1 eq 'X grid size:'))
ny=modhead2(where(modhead1 eq 'Y grid size:'))
dx=modhead2(where(modhead1 eq 'X grid resolution:'))
dy=modhead2(where(modhead1 eq 'Y grid resolution:'))

nfields=modhead2(where(modhead1 eq 'Number of series:'))

xo=float(xo(0))
yo=float(yo(0))
nx=int(nx(0))
ny=int(ny(0))
dx=float(dx(0))
dy=float(dy(0))

nfields=int(nfields(0))


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
; series headers

readdata,filename,fieldhead,nrows=6,nskip=17
fieldhead=strtrim(fieldhead,2)
fieldhead=fieldhead(5:*,*)

;-----------------------------------------------------------------------
; read data

readdata,filename,data,nrecs,nskip=24
data=strtrim(data,2)
  
date=replicate({datetime},nrecs)

for i=long(0),nrecs-1 do begin
  tempdat=strsplit(data(0,i),' ',/Extract,/Preserve_Null)
  temptime=strsplit(tempdat(1),':',/Extract,/Preserve_Null)
  tempdat=strsplit(tempdat(0),'/',/Extract,/Preserve_Null)
  
  date(i).year=tempdat(2)
  date(i).month=tempdat(1)
  date(i).day=tempdat(0)
  date(i).hour=temptime(0)
  date(i).minute=temptime(1)
  date(i).second=temptime(2)
endfor
  
date.julian=julday(date.month,date.day,date.year,date.hour,date.minute,date.second)
  
;-----------------------------------------------------------------------  
; select field

index=where(fieldhead(*,0) eq location and $
            fieldhead(*,1) eq species and $
	    fieldhead(*,3) eq field and $
	    fieldhead(*,5) eq level,nindex)


if (nindex gt 0)then begin
  seldata=reform(data(index(0),*))
  units=fieldhead(index(0),4)
endif else begin
  print,'Selected series not found'
endelse

head=reform(fieldhead(index(0),*))


;-----------------------------------------------------------------------  
; make up matrix array

times=date.julian[uniq(date.julian,sort(date.julian))]
ntime=n_elements(times)

matrixarr=fltarr(ntime,nx,ny)

for it=0l,ntime-1 do begin

  index=where(date.julian eq times(it),nit)
   
  for i=0,nit-1 do begin
    ix=int(xgrid(index(i)))
    iy=int(ygrid(index(i)))
    matrixarr(it,ix-1,iy-1)=seldata(index(i))
   endfor    
  
endfor

matrixdates=replicate({datetime},ntime)
matrixdates=date[uniq(date.julian,sort(date.julian))]

;-----------------------------------------------------------------------  
; ensure -180 to 180

index1=where(lon2d gt 180.0,nindex)
if(nindex gt 0)then begin
  lon2d(index1)=lon2d(index1)-360.0
  index2=sort(lon2d(*,0))
  nshift=index2(0)-2
  
  for i=0l,ntime-1 do begin
    fielddata=reform(matrixarr(i,*,*))
    fielddata=Shift(fielddata,nshift,0)
    matrixarr(i,*,*)=fielddata
  endfor  
  
  lat2d     =shift(lat2d,nshift,0)
  lon2d     =shift(lon2d,nshift,0)  
endif 
end
