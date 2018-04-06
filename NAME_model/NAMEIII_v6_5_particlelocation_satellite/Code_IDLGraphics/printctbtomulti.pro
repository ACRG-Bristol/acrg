pro printctbtomulti,datadir,selgrid,count,selspecies=selspecies,$
 selfield=selfield,sellevel=sellevel,SensorID=SensorID,filetext=filetext,$
 template=template
 
;--------------------------------------------------------------------
; Procedure to generate sequence of selected plots
; 
;
; DBR Feb 2002 - Original single detection run data processing code
; MCH Jan 2007 - major rewrite to deal with multiple detections
;                in single NAME III run.
;
; Arguments
; required:
;  datadir      : (string)  run directory
;  selgrid      : (string)  grid ('grid1' or 'grid2')
;  count        : the loop counter used to identify NAME run
;
; optional:
;  selspecies   : (string or string array) species(s) to plot  
;  selfield     : (string or string array) field(s) to plot
;  sellevel     : (string or string array) level(s) to plot
;
;  filetext     : (string)  text added to filenames
;  template     : string string to include when listing Fields* files
;                 eg ='0000' to search for midnight data only
;  SensorID     : CTBTO sensor name - used to name CTBTO upload
;                 files
;    
; note:
; (1) selfield and sellevel must contain same number of requests
;
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
; arguments

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse
if(n_elements(template) eq 0)then template='*'

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

fieldfiles=findfile(datadir+delim+'Fields_'+selgrid+'_*'+template+'*')

fieldfiles=reverse(fieldfiles)

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

readfieldhead,fieldfiles(0),modtitle,fieldhead,xyrel


;---------------------------------------------------------------------------
;directories

dirfile=strarr(numfiles)

for it=0,numfiles-1 do begin
  filesplit=strsplit(fieldfiles(it),delim,/extract,/preserve_null)
  for is=1,n_elements(filesplit)-2 do begin
    dirfile(it)=dirfile(it)+delim+filesplit(is)
  endfor
endfor

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
  numfields=n_elements(where(fieldhead(*,0)))
  fieldlist=reform(fieldhead(*,3))
  levellist=reform(fieldhead(*,5))
endif else begin
  numfields=n_elements(selfield)
  fieldlist=selfield
  levellist=sellevel
endelse 


;-----------------------------------------------------------------------     
; open file for printing

ctbto_header_file=datadir+'/ctbto_data.txt'

;Setting up strings
SensorStartDate='a'
SensorStartTime='b'
SensorStopDate='c'
SensorStopTime='d'

readdata,ctbto_header_file,data,nskip=count,Delim=[' '],nrows=1,null=0
SensorStopDate=data(4)
SensorStopTime=data(5)

loc=strmid(SensorID,0,5)
filetem=loc+'.'+SensorStopDate+SensorStopTime+'.EGRR'+'.txt'
openw,1,dirfile(numfiles-1)+'/'+filetem

; read in and write out ctbto file header information.
ctbto_header=strarr(1)

readdata,ctbto_header_file,ctbto_header,nrows=1,nskip=count,null=0
   
printf,1,ctbto_header


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
 
      filetime=strmid(fieldfiles(it),(strpos(fieldfiles(it),'.txt'))-12,12) 
print,'field',level
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
