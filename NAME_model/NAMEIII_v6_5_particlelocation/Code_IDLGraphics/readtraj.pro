pro readtraj,filename,modtitle,date,x,y,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
    zqfe,zqnh,zp,zfl,traveltime,namever=namever

;-----------------------------------------------------------------------

; procedure to read trajectory model output
;
; DBR 02/2002
;
; Changes:
; 01/02/2006 ARJ  updated for NAME III trajectory output
; 02/07/2006 CSW  updated to accept variable numbers of columns 
;
;-----------------------------------------------------------------------

if (n_elements(namever) eq 0) then namever=2

;-----------------------------------------------------------------------
; define date/time structure

date={datetime,year:-4713,month:1,day:1,hour:12,minute:0,second:0,julian:0.0D} 

;-----------------------------------------------------------------------

; read title text

modtitle=strarr(1)

if (namever eq 3) then begin
  ;NAME III
  readdata,filename,modtitle,nrows=8
endif else begin
  ;NAME
 readdata,filename,modtitle,nrows=10
endelse
  
modtitle=strtrim(modtitle,2)

; read and check coord systems in the column headers (NAME III only)
; also extract particle release time from the first line of data

if (namever eq 3) then begin
  ;NAME III
  readdata,filename,headers,nrows=1,nskip=8
 
  xcoord=strtrim(headers(6),2)
  ycoord=strtrim(headers(7),2)
  zcoord=strtrim(headers(8),2)
  if (strpos(xcoord(0),'Lat-Long') eq -1 or strpos(ycoord(0),'Lat-Long') eq -1) $
    then begin
    print, 'Error in readtraj routine: coord system is not Lat-Long'
    stop
  endif 
  
  il=strpos(zcoord(0),'(')
  ir=strpos(zcoord(0),')')
  if (il eq -1 or ir eq -1) then begin
    print, 'Error in readtraj routine: invalid header for vertical coord'
    stop
  endif
  zcoord=strmid(zcoord(0),il+1,ir-il-1)
  
  ; Check for presence of met data in NAMEIII file

  mettest=n_elements(headers)
  if (mettest le 16) then begin
    nfields=15
  endif else begin
    nfields=29
  endelse
  
  readdata,filename,data,nrows=1,nskip=9
  rel_datetime=replicate({datetime},1)
  tempdat=strtrim(data(3),2)
  tempdat=strsplit(tempdat,' ',/Extract,/Preserve_Null)
  temptime=strsplit(tempdat(1),':',/Extract,/Preserve_Null)
  tempdat=strsplit(tempdat(0),'/',/Extract,/Preserve_Null)  
  
  rel_datetime.year=tempdat(2)
  rel_datetime.month=tempdat(1)
  rel_datetime.day=tempdat(0)
  rel_datetime.hour=temptime(0)
  rel_datetime.minute=temptime(1)
  rel_datetime.second=temptime(2)

  rel_datetime.julian=julday(rel_datetime.month,rel_datetime.day,rel_datetime.year,$
                        rel_datetime.hour,rel_datetime.minute,rel_datetime.second)
  
endif

; read data 

; define arrays for variables which may not be present
temp=fltarr(1)
press=fltarr(1)
ptemp=fltarr(1)
bldepth=fltarr(1)
cloud=fltarr(1)
rh=fltarr(1)
wspd=fltarr(1)
wdir=fltarr(1)


if (namever eq 3) then begin
  ;NAME III
  
  readdata,filename,data,nrecs,nskip=9
  data=strtrim(data,2)
  
  date=replicate({datetime},nrecs)

  for i=long(0),nrecs-1 do begin
    tempdat=strsplit(data(4,i),' ',/Extract,/Preserve_Null)
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

  traveltime=double(float(data(5,*))/(3600.0D*24.0D))
  x=float(data(6,*))
  y=float(data(7,*))
  z=float(data(8,*))
  
  if (mettest gt 16) then begin
    temp=float(data(21,*))
    press=float(data(22,*))
    ptemp=float(data(23,*))
    bldepth=float(data(24,*))
    cloud=float(data(25,*))
    rh=float(data(26,*))
    wspd=float(data(27,*))
    wdir=float(data(28,*))
  endif
  
  if (zcoord eq 'm agl') then begin
    zqfe=z
    zqnh=0.0
    zfl=0.0
    zp=0.0
  endif else if (zcoord eq 'm asl') then begin
    zqfe=0.0
    zqnh=z
    zfl=0.0
    zp=0.0
  endif else if (zcoord eq 'FL') then begin
    zqfe=0.0
    zqnh=0.0
    zfl=z
    zp=0.0
  endif else begin
    print, 'Error in readtraj routine: invalid vertical coord system'
    stop 
  endelse
  
endif else begin
  ;NAME
  
  readdata,filename,data,nrecs,nskip=12
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

  x=float(data(1,*))
  y=float(data(2,*))
  zqfe=float(data(3,*))
  zqnh=float(data(4,*))
  zfl=float(data(5,*))
  zp=float(data(6,*))
  
endelse

; fix up title text information for NAME III (so that NAME III returns the
; same details in modtitle as a NAME II run)

if (namever eq 3) then begin
  ;NAME III
  ; modtitle(5) = 'Forward trajectories'/'Backward trajectories'
  modtitle(5)=modtitle(4)

endif

if (namever ne 3) then begin
  modtitle(1)=modtitle(6)
endif

End
