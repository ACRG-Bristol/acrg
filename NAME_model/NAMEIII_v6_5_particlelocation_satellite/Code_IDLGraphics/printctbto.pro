pro printctbto,datadir,selgrid,selspecies=selspecies,$
 selfield=selfield,sellevel=sellevel,filetext=filetext,$
 template=template
 
;--------------------------------------------------------------------
; Procedure to generate sequence of selected plots
; 
;
; DBR Feb 2002
;
; Arguments
; required:
;  datadir      : (string)  run directory
;  selgrid      : (string)  grid ('grid1' or 'grid2')
;
; optional:
;  selspecies   : (string or string array) species(s) to plot  
;  selfield     : (string or string array) field(s) to plot
;  sellevel     : (string or string array) level(s) to plot
;
;  filetext     : (string)  text added to filenames
;  template     : string string to include when listing Fields* files
;                 eg ='0000' to search for midnight data only

;    
; note:
; (1) selfield and sellevel must contain same number of requests
;
;-----------------------------------------------------------------------

print,'filetext',filetext
;-----------------------------------------------------------------------
; arguments

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse
if(n_elements(template) eq 0)then template='*'

;-----------------------------------------------------------------------
; define date/time structure

date={datetime,year:-4713,month:1,day:1,hour:12,minute:0,second:0,julian:0.0D} 

;-----------------------------------------------------------------------
; graphics initialisation

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

;-----------------------------------------------------------------------
; find all files and read header 

fieldfiles=findfile(datadir+delim+'Fields_'+selgrid+'*'+template+'*')
fieldfiles=reverse(fieldfiles)

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

readfieldhead,fieldfiles(0),modtitle,fieldhead,xyrel

lon=xyrel(0)
lat=xyrel(1)

;---------------------------------------------------------------------------
; release start and end times

cstartdt=var_to_dt(strmid(modtitle(4),36,4),strmid(modtitle(4),33,2),$
                  strmid(modtitle(4),30,2),strmid(modtitle(4),22,2))
cenddt=var_to_dt(strmid(modtitle(5),36,4),strmid(modtitle(5),33,2),$
                  strmid(modtitle(5),30,2),strmid(modtitle(5),22,2))

cstartstr1=string(cstartdt.year,format='(i4)')+string(cstartdt.month,format='(i2.2)')+$
          string(cstartdt.day,format='(i2.2)')
cstartstr2=string(cstartdt.hour,format='(i2.2)')

cendstr1=string(cenddt.year,format='(i4)')+string(cenddt.month,format='(i2.2)')+$
          string(cenddt.day,format='(i2.2)')
cendstr2=string(cenddt.hour,format='(i2.2)')

relduration=(cstartdt.julian-cenddt.julian)*24

;Old NAME II line
;rate=real(strmid(modtitle(6),20,9))

; Change made for slight difference in NAME III file format
; even when outputting as NAME II which it must do for this script.
rate=real(strmid(modtitle(6),21,13))
mass=rate*relduration*3600
;---------------------------------------------------------------------------
;directories

runname=strarr(numfiles)
dirfile=strarr(numfiles)
for it=0,numfiles-1 do begin
  filesplit=strsplit(fieldfiles(it),delim,/extract,/preserve_null)
  for is=1,n_elements(filesplit)-2 do begin
    dirfile(it)=dirfile(it)+delim+filesplit(is)
  endfor
  runname(it)=filesplit(n_elements(filesplit)-2)
endfor
array=dirfile(where(dirfile))
array=array[uniq(array,sort(array))]
numdir=n_elements(array)

;---------------------------------------------------------------------------
; times

timelist=strarr(numfiles)
timestr=strarr(numfiles)
timedate=replicate({!datetime},numfiles)
for it=0,numfiles-1 do begin

  timelist(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12) 
  timestr(it)=strmid(timelist(it),8,4)+'Z '+strmid(timelist(it),6,2)+$
   '/'+strmid(timelist(it),4,2)+'/'+strmid(timelist(it),0,4)

endfor

timedate.year=strmid(timelist(it),0,4)
timedate.month=strmid(timelist(it),4,2)
timedate.day=strmid(timelist(it),6,2)
timedate.hour=strmid(timelist(it),8,2)
timedate.julian=julday(timedate.month,timedate.day,timedate.year,$
                       timedate.hour,timedate.minute,timedate.second)

dt=24*(timedate(0).julian-timedate(1).julian)
fctime=dt*numfiles

;---------------------------------------------------------------------------
; all species

if(n_elements(selspecies) eq 0)then begin
  array=fieldhead(where(fieldhead(*,1)),1)
  specieslist=array[uniq(array,sort(array))]
  numspecies=n_elements(specieslist)
endif else begin
  numspecies=n_elements(selspecies)
  specieslist=selspecies
endelse

;---------------------------------------------------------------------------
; all fields

if(n_elements(selfield) eq 0 or $
   n_elements(sellevel) eq 0) then begin
  numfields=n_elements(fieldhead(*,0))
  fieldlist=reform(fieldhead(*,3))
  levellist=reform(fieldhead(*,5))
endif else begin
  numfields=n_elements(selfield)
  fieldlist=selfield
  levellist=sellevel
endelse 

;---------------------------------------------------------------------------
; read first field

filename=fieldfiles(0)
filetime=strmid(fieldfiles(0),(strpos(fieldfiles(0),'Field'))+13,12) 

readfield,dirfile(0),filetime,selgrid,fieldlist(0),levellist(0),$
 specieslist(0),lon1,lat1,dat1,head1

dlon=lon1(1,0)-lon1(0,0)
dlat=lat1(0,1)-lat1(0,0)

;-----------------------------------------------------------------------     
; open file for printing

loc=strmid(runname(0),0,5)
filetem=loc+'.'+cstartstr1+cstartstr2+'.EGRR'+'.txt'
openw,1,dirfile(numfiles-1)+'/'+filetem

loc='"'+loc+'"' 

printf,1,lon,lat,cendstr1,cendstr2,cstartstr1,cstartstr2,$
 mass,fctime,dt,dt,dlon,dlat,loc,$
format='(f7.2,2x,f6.2,1x,a8,1x,a2,1x,a8,1x,a2,1x,e8.2,2x,i4,1x,i2,1x,i2,1x,f4.2,1x,f4.2,1x,a7)'

;-----------------------------------------------------------------------
; loop through fields and species

ip=-1

for ifield=0,numfields-1 do begin 
 
  field=fieldlist(ifield)
  level=levellist(ifield)
  
  for ispec=0,numspecies-1 do begin

    species=specieslist(ispec)

     
;-----------------------------------------------------------------------        
; loop through times

    
    for it=0,numfiles-1 do begin

      filename=fieldfiles(it)
      filetime=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12) 

      readfield,dirfile(it),filetime,selgrid,field,level,$
	  species,lon1,lat1,dat1,head1

      xysize=size(reform(lon1))
      nx=xysize(1)
      ny=xysize(2)
      
      for ix=0,nx-1 do begin
        for iy=0,ny-1 do begin

          if(dat1(ix,iy) gt 0)then begin
            printf,1,lat1(ix,iy),lon1(ix,iy),it+1,dat1(ix,iy),$
              format='(f7.2,2x,f7.2,2x,i4,2x,e13.7)'
          endif
        endfor
      endfor
      
    endfor


  endfor
  
endfor

close,/all

end
