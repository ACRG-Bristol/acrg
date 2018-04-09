pro writefile,datadir,seltime,selgrid,title,fieldhead,fielddata,lon2d,lat2d, $
              filename=filename,writezeros=writezeros, $
              boundarythreshold=boundarythreshold

; ----------------------------------------------------------------------------------------
;
; Lois Huggett
; 11/09/09
;
; Eike Mueller
; 10/11/10
;		 
; Produce output file in NAME II-style
; output files which can be plotted using plotfield.
; the fielddata can be (1+2) dimensional, where the first index labels different columns,
; the same applies for the fieldhead information, the first dimension of the fieldhead 
; array has to be identical to that of the data array.
;
; ----------------------------------------------------------------------------------------
;
; inputs:
;   datadir          : String with filepath to output files directory
;   seltime          : time string in output filename
;   selgrid          : name of grid in output filename
;   title            : String array of file title information for output file
;   fieldhead           : String array of fieldhead information for output file, can be 
;                      (1+1) dimensional
;   fielddata        : Array containing data to output, can be (1+2) dimensional
;   lon2d            : array with longitudes
;   lat2d            : array with latitudes
;   [filename]       : can be used to specify filename explictly
;   [writezeros]     : if this is set to 1, empty rows will be written to
;                      the output file (default = 0)
;   
; Note: 
; * Unless filename is given explicitly, the name of the output file is
;   datadir+'/'+'Fields_'+selgrid+'_'+seltime+'.txt'
;   where the delimiter '/' might depend on the operating system
; * The longitude and latitude arrays have to be compatible with the data
;
;

if n_elements(writezeros) eq 0 then writezeros=0
if n_elements(boundarythreshold) eq 0 then boundarythreshold=0

; ----------------------------------------------------------------------------------------

; error handling. A value of '1' indicates everything is ok
status = 1

; ----------------------------------------------------------------------------------------
; generate filename

if (strpos(!version.os,'NT') gt 0) then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

if (n_elements(filename) eq 0) then begin 
  filename=datadir+delim+'Fields_'+selgrid+'_'+seltime+'.txt'
endif


; ----------------------------------------------------------------------------------------
; if data is 2d and fieldhead is 1d, add one array dimension to store column index

if (( size(fielddata,/n_dimensions) eq 2) and (size(fieldhead,/n_dimensions) eq 1)) then begin
  nx = n_elements(fielddata(*,0))
  ny = n_elements(fielddata(0,*))
  fielddata = reform(fielddata,1,nx,ny)
  nh = n_elements(fieldhead)
  fieldhead = reform(fieldhead,1,nh)
endif

; ----------------------------------------------------------------------------------------
; check that dimensions of data and fieldhead arrays match

if ( n_elements(fielddata(*,0,0)) ne n_elements(fieldhead(*,0)) ) then begin
  print,'ERROR: size of fielddata and fieldhead arrays do not match'
  status = -1
endif

; ----------------------------------------------------------------------------------------
; set number of output columns and modify title information
ncolumns = n_elements(fielddata(*,0,0))
title(16) = ' Number of fields: '+string(ncolumns)

if (status eq 1) then begin 
   
; ----------------------------------------------------------------------------------------
; Generate output file

; split 2-D data array (also xlon and ylat) back into 1-D array for output

  xgridsize=strsplit(title(12),':',/extract)
  xgridsize=long(xgridsize(1))

  ygridsize=strsplit(title(13),':',/extract)
  ygridsize=long(ygridsize(1))

  date=strmid(seltime,6,2)+'/'+strmid(seltime,4,2)+'/'+strmid(seltime,0,4)
  time=strmid(seltime,8,4)
  datetime=strarr(ncolumns)
  datetime(*)='    '+time+'UTC '+date

; open file and write data
  get_lun,OUT
  openw,OUT,filename,error=iostat

  if (iostat eq 0) then begin
    for i=0,n_elements(title)-1 do begin
      if (i eq 0) then begin
        printf,OUT,title(i)
     endif else begin
         subtitle=strsplit(title(i),':',/extract)
         printf,OUT,subtitle(0)+':',subtitle(1),format='(a-22,a)'
     endelse
   endfor
   printf,OUT
  
   fieldheadformat='(4a20,'+strcompress(string(ncolumns),/REMOVE_ALL)+'(a24,","))'
   for iout=0,5 do begin
     printf,OUT,',',',',',',',',fieldhead(*,iout),format=fieldheadformat
   endfor
   printf,OUT,'X grid,','Y grid,','Longitude,','Latitude,',datetime,format=fieldheadformat
   printf,OUT,' '
  
  
   dataformat='(4(f19.6,a1),'+strcompress(string(ncolumns))+'(E24.8,","))'
   for igrid=0L,xgridsize-1 do begin
     for jgrid=0L,ygridsize-1 do begin
       if ( (writezeros eq 1) or (where(fielddata(*,igrid,jgrid) ne 0)) ne [-1]) then begin
          printf,OUT, igrid+1.5,',',            $
                      jgrid+1.5,',',            $
                      lon2d(igrid,jgrid),',',   $
	                    lat2d(igrid,jgrid),',',   $
                      fielddata(*,igrid,jgrid), $
                      format=dataformat
        endif
      endfor
    endfor

; open file and write data
    close,OUT
    free_lun,OUT
    
   if (boundarythreshold ge 0) then begin

     get_lun,OUT_BOUNDARIES
     openw,OUT_BOUNDARIES,filename+'.boundaries',error=iostat
     printf,OUT_BOUNDARIES,"# min_lon, max_lon, min_lat, max_lat"
     if ( size(fielddata,/n_dimensions) eq 3) then begin
       nfieldspecs = n_elements(fielddata(*,0,0))     
     endif else begin
       nfieldspecs = 1
     endelse
     count = 0
     max_lat = -90.
     min_lat = 90.
     max_lon = -180.
     min_lon = 180.
     ; Shift the longitudes into the range [-180,180]. This will give
     ; sensible plot boundaries as long as the plume does not cross the 
     ; date line:
     ;
     ; -180                         0                          180
     ; [---------##################-+----------------------------]
     ;         min_lon          max_lon
     ;
     for i = 0, nfieldspecs - 1 do begin
       for igrid=0L,xgridsize-1 do begin
         for jgrid=0L,ygridsize-1 do begin
           if (fielddata(i,igrid,jgrid) ge boundarythreshold) then begin
             max_lat = max([lat2d(igrid,jgrid),max_lat])
             min_lat = min([lat2d(igrid,jgrid),min_lat])
             lon_shifted = lon2d(igrid,jgrid)
             if (lon_shifted gt 180.0) then lon_shifted = lon_shifted - 360.
             max_lon = max([lon_shifted,max_lon])
             min_lon = min([lon_shifted,min_lon])
             count = count+1
           endif
         endfor
       endfor
     endfor
     ; If no value exceeds the threshold, plot the entire world
     if (count eq 0) then begin
       min_lat = -90.
       max_lat = 90.
       min_lon = -180.
       max_lon = 180.
     endif else begin
       if (max_lon gt 180.) then max_lon = max_lon - 360.
       if (min_lon gt 180.) then min_lon = min_lon - 360.
     endelse
     ; If max_lon > 175 and min_lon < -175 and count > 0, 
     ; the plume has probably crossed the dateline. 
     ;
     ; -180                         0                          180
     ; [###########-----------------+------------------------####]
     ;  min_lon                                           max_lon
     ;
     ; Repeat the process, but without shifting the longitudes first:
     ;
     ; 0                           180                         360
     ; [------------------------##################---------------]
     ;  min_lon                                           max_lon
     ;     
     if ((min_lon lt -175.) and (max_lon gt 175.) and (count gt 0)) then begin
       max_lon = 0.
       min_lon = 360.
       for i = 0, nfieldspecs - 1 do begin
         for igrid=0L,xgridsize-1 do begin
           for jgrid=0L,ygridsize-1 do begin
             if (fielddata(i,igrid,jgrid) ge boundarythreshold) then begin
               max_lon = max([lon2d(igrid,jgrid),max_lon])
               min_lon = min([lon2d(igrid,jgrid),min_lon])
             endif
           endfor
         endfor
       endfor
       ; shift back into the range [-180, 180]
       if (max_lon gt 180.) then max_lon = max_lon - 360.
       if (min_lon gt 180.) then min_lon = min_lon - 360.
       ; If still max_lon > 175 and min_lon < -175, the plume covers more
       ; than 180 degrees in this case plot the entire world
       ; 0                           180                         360
       ; [##-----------------------################################]
       ;  min_lon                                           max_lon       
       if ( (min_lon lt -175.) and (max_lon gt 175.) )then begin
         min_lat = -90.
         max_lat = 90.
         min_lon = -180.
         max_lon = 180.
       endif       
     endif
     
     
     dataformat='(E24.8," ",E24.8," ",E24.8," ",E24.8)'
     printf, OUT_BOUNDARIES, min_lon, max_lon, min_lat, max_lat
    close,OUT_BOUNDARIES
    free_lun,OUT_BOUNDARIES
   endif
  endif else status = -1
 endif
  
end
