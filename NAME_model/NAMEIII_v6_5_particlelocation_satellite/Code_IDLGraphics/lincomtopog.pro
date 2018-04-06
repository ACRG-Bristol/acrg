pro lincomtopog,datadir,seltime,selgrid,selfield,lon2d,lat2d,$
topog,topog_long,topog_lat,$
 modtitle,fieldhead,timeframe=timeframe

;--------------------------------------------------------------
; procedure to read lincom topog files data
;
; MCH 09/02/2004
;--------------------------------------------------------------

print,'start lincom topog read'

filename='/data/local/apmh/NAMEIII/NAMEIII_1_3e/NAMEIII/NameIII/SingleSite/TopogC1.txt'

;-----------------------------------------------------------------------
; define date/time structure

date={datetime,year:-4713,month:1,day:1,hour:12,minute:0,second:0,julian:0.0D}

;-----------------------------------------------------------------------  
; read title text

modhead1=strarr(1)
modhead2=strarr(1)

readtitle,filename,modtitle,modhead1,modhead2

;-----------------------------------------------------------------------  
; extract grid details

xo=modhead2(where(modhead1 eq 'X grid origin:'))
yo=modhead2(where(modhead1 eq 'Y grid origin:'))
nx=modhead2(where(modhead1 eq 'X grid size:'))
ny=modhead2(where(modhead1 eq 'Y grid size:'))
dx=modhead2(where(modhead1 eq 'X grid resolution:'))
dy=modhead2(where(modhead1 eq 'Y grid resolution:'))

nfields=modhead2(where(modhead1 eq 'Number of fields:'))

xo=float(xo(0))
yo=float(yo(0))
nx=fix(nx(0))
ny=fix(ny(0))
dx=float(dx(0))
dy=float(dy(0))

nfields=fix(nfields(0))

topog=fltarr(nx,ny)

;-----------------------------------------------------------------------  
; define lon/lat arrays

lon2d=fltarr(nx,ny)
lat2d=fltarr(nx,ny)

  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      lon2d(ix,iy)=xo+float(ix)*dx
      lat2d(ix,iy)=yo+float(iy)*dy
    endfor
  endfor

;-----------------------------------------------------------------------  
; ensure -180 to 180

index1=where(lon2d gt 180.0,nindex1)
if(nindex1 gt 0)then begin
  lon2d(index1)=lon2d(index1)-360.0
  index2=sort(lon2d(*,0))
  nshift=index2(0)-2
  lat2d     =shift(lat2d,nshift,0)
  lon2d     =shift(lon2d,nshift,0)  
endif 

;-----------------------------------------------------------------------  
; read data

 data=fltarr(2)
;
readdata,filename,data,nrecs,nfields=nfields+5,nskip=36

   if (timeframe eq 'absolute')then begin
     ; NAME III absolute time
     ;status=dc_read_free(filename,date,xlon,ylat,zheight,data,$
     ;/column, Delim=[','],get_columns=[1,2,3,4,index(0)+5],resize=[2,3,4,5],nskip=24)
   endif else if (timeframe eq 'relative')then begin
     ;date1 = REPLICATE({!DT}, 500000)
     ;status=dc_read_free(filename,date1,xlon,ylat,zheight,data,  $
     ;/column, Dt_Template=[-1], get_columns=[1,2,3,4,5], $
     ;resize=[2,3,4,5],nskip=24)
   endif else begin
     print,'Error; No time frame specified for NAME III data'
   endelse

data=strtrim(data,2)
xlon=float(data(1,*))
ylat=float(data(2,*))
zheight=float(data(3,*))
temp=fltarr(nrecs)
temp(*)=float(data(index(0)+4,*))
data=temp

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
; generate 2d array

;fielddata=fltarr(nx,ny)
;
    ; NAME III this works for files where all data points have been output for
    datacount=0L
    ; comment:  for i=0l,n_elements(xgrid)-1 do begin
    for i=0l,nx-1 do begin
      for j=0,ny-1 do begin
        topog(i,j)=float(data(datacount))
        datacount=datacount+1
      endfor
    endfor

end
