pro areaatrisk,datadir,selgrid,selspecies=selspecies,$
 selfield=selfield,sellevel=sellevel,filetext=filetext,$
 template=template,last=last,namever=namever,datamap=datamap,$
 SelTimeAverage=SelTimeAverage,threshold=threshold,$
 selcase=selcase,selname=selname,sourceloc=sourceloc,$
 xdim=xdim,ydim=ydim,aspect=aspect,margin=margin,pacific=pacific,$
 AspectRatio=AspectRatio
 
;--------------------------------------------------------------------
; Procedure to generate lat/long limits in file for selected field(s)
; 
;
; Lois Huggett 23/03/09
;
; Arguments
; required:
;  datadir      : (string)  run directory
;  selgrid      : (string)  grid ('grid1' or 'grid2')
;
; optional:
;  selspecies   : (string or string array) species to calculate AAR for  
;  selfield     : (string or string array) 'Area At Risk', or exactly how
;                 it is written in the NAME output file
;  sellevel     : (string or string array) level to calculate AAR for
;  SelTimeAverage  : (string) time averaged field to plot
;  filetext     : (string)  text added to filenames
;  template     : string - string to include when listing Fields* files
;                 eg ='0000' to search for midnight data only
;  last         : integer (0 or 1) =1 plots last timestep only
;  namever      : integer indicating version of model e.g. 3 - NAME III(PPM)
;  datamap      : set to 1 for map size based on plume area above threshold level
;  threshold    : (float) AAR min level to take into account
;  selcase      : case number (name iii only - defaults to 1 if not given)
;  selname      : name of field requirement (can be used to specify
;                 a field which is not uniquely defined by the other
;                 parameters; name iii only - defaults to ' ' if not given)
;  sourceloc    : effective source location in the case of multiple sources
;  xdim         : (float) x dimension of map in km.  Defaults to 11km (roughly
;                 l:50000 on A4 paper).
;  ydim         : (float) y dimension of map in km.  Defaults to xdim/(root 2).
;  aspect       : (integer) 0=test for best map, portrait/landscape,
;                 (swapping x and y dimensions) 
;                 1=use xdim as x size and ydim as y size of map only.
;  margin       : (float) fraction of map dimension to be used as margin
;  pacific      : (integer) flag to indicate if map area is likely to cross
;                  the international date line
;  AspectRatio  : flag to crontal aspect ratio of plot area.
;                 only used if only xdim provided
;                 option are: 'A4', 'Square'
;
; note:
; (1) selfield and sellevel must contain same number of requests,
;     otherwise all attribution fields plotted
;
; example calls:

;  a NAME III call:
;  areaatrisk,'/home/fr1100/apnm/namev66/testv66','Fields_grid2',$
;   selspecies='TRACER',sourceloc=[-3.0,55.0],xdim=5.0,ydim=8.0,aspect=1,$
;   selfield='Area At Risk',namever=3,sellevel='Z = 50.00000 m agl'
;
;  a NAME II-style call:
;  areaatrisk,'/home/fr1100/apnm/namev66/testv66','grid2',pacific=1
;
;
; subroutine calls:
;  readfieldhead
;  readfield
;  mapcalc
;  readtitle
;  readdata
;  distance
;
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
; arguments

if(n_elements(selfield) eq 0) then selfield='Area At Risk'
if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse
if(n_elements(template) eq 0)then template='*'
if(n_elements(last) eq 0)then last=0
if(n_elements(namever) eq 0)then namever=2
if(n_elements(datamap) eq 0)then datamap=0
if(n_elements(threshold) eq 0)then threshold=0.5
if(n_elements(selcase) eq 0)then begin
  selcase=1
  printcasenum=0
endif else begin
  printcasenum=1
endelse
if(n_elements(xdim) eq 0)then xdim=11.0
;Moved further down script
;if(n_elements(ydim) eq 0)then ydim=xdim/1.4142
if(n_elements(aspect) eq 0)then aspect=0
if(n_elements(margin) eq 0)then margin=0.05
if(n_elements(pacific) eq 0) then pacific=0
if(n_elements(AspectRatio) eq 0) then AspectRatio='A4'

;-----------------------------------------------------------------------
; graphics initialisation

if(strpos(!version.os,'NT') ge 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

if(strpos(!version.os,'Linux') ge 0)then begin
  linux=1
endif else begin
  linux=0
endelse


;-----------------------------------------------------------------------
; find all files and readheader 

if (namever eq 3)then begin
  ; NAME III
  casenum=STRTRIM(selcase,2)
  fieldfiles=findfile(datadir+delim+selgrid+'_C'+casenum+'_T*')

  fieldfiles1=fieldfiles
  numfiles=n_elements(fieldfiles)

  FileSplit = STRSPLIT(fieldfiles(0),'_')
  NFileSplit = n_elements(FileSplit)

  TSection=strarr(numfiles)

; split file names into array elements. '_' used as dividing points
; this enables us to sort and order files correctly.
  for i=0,numfiles-1 do begin
    FileSplit = STRSPLIT(fieldfiles(i),'_',/Extract)
    TSection(i) = FileSplit(NFileSplit-2)  
  endfor

; search through TSection array to correctly order files.

  for i=1,numfiles do begin
    Tname=STRCOMPRESS('T'+string(i), /Remove_All)
    for j=0,numfiles-1 do begin
      if (TSection(j) eq Tname)then begin
        fieldfiles1(i-1) = fieldfiles(j)
      endif
    endfor
  endfor

  fieldfiles = fieldfiles1

endif else begin
  ;NAME
  fieldfiles=findfile(datadir+delim+'Fields_'+selgrid+'_*'+template+'*')
endelse

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif
  
readfieldhead,fieldfiles(0),modtitle,fieldhead,xyrel,$
  grid,namever=namever,timeframe=timeframe

; all times

timelist=strarr(numfiles)
for it=0,numfiles-1 do begin
  if (namever eq 3)then begin
      ; NAME III
      if (it eq 0) then begin
        FileTime=strarr(2)
      endif
      FileSplit = STRSPLIT(fieldfiles(it),'_',/Extract)
      NFileSplit = n_elements(FileSplit)
      TimeBit = FileSplit(NFileSplit-1)
      FileTime = STRSPLIT(TimeBit,'.t',/Extract)
      timelist(it)=FileTime(0)
  endif else begin
    ; NAME
;    timelist(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12)
    timelist(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'.txt'))-12,12) 
  endelse
endfor

if(last ne 0)then begin
  timelist=timelist(numfiles-1)
  fieldfiles=fieldfiles(numfiles-1)
  numfiles=1
endif

; all species

if (namever eq 3)then begin
  ; NAME III
  if(n_elements(selspecies) eq 0)then begin
    array=fieldhead(where(fieldhead(*,3)),3)
    specieslist=array[uniq(array,sort(array))]
    numspecies=n_elements(specieslist)
  endif else begin
    numspecies=n_elements(selspecies)
    specieslist=selspecies
  endelse
endif else begin
  ;NAME
  if(n_elements(selspecies) eq 0)then begin
    array=fieldhead(where(fieldhead(*,1)),1)
    specieslist=array[uniq(array,sort(array))]
    numspecies=n_elements(specieslist)
  endif else begin
    numspecies=n_elements(selspecies)
    specieslist=selspecies
  endelse
endelse

; all fields

if n_elements(numfields) eq 0 then begin
  numfields=n_elements(where(fieldhead(*,0)))
endif

if (namever eq 3)then begin
  ; NAME III
  
  if ((n_elements(selfield) eq 0) OR (n_elements(sellevel) eq 0)) then begin
    fieldlist=reform(fieldhead(*,2))
    fieldlist=fieldlist(where(fieldlist))
    levellist=reform(fieldhead(*,14))
    levellist=levellist(where(levellist))
  endif else begin
    numfields=n_elements(selfield)
    fieldlist=selfield
    if (n_elements(sellevel) ne numfields) then begin
      print, 'ERROR - incorrect number of entries in level list'
      exit
    endif else begin
      levellist=sellevel
    endelse
  endelse
  
  if (n_elements(SelTimeAverage) eq 0) then begin
    timeaveraginglist=reform(fieldhead(*,7))
    timeaveraginglist=timeaveraginglist(where(timeaveraginglist))
  endif else begin
    timeaveraginglist=SelTimeAverage
  endelse
  
  if (n_elements(timeaveraginglist) ne numfields) then begin
    print,'ERROR - incorrect number of entries in time-averaging list'
    exit
  endif
  
  if (n_elements(selname) eq 0) then begin
    FieldNameList=replicate(' ',numfields)
  endif else if (n_elements(selname) ne numfields) then begin
    print,'ERROR - incorrect number of entries in field name list'
    exit
  endif else begin
    FieldNameList=selname
  endelse 
   
endif else begin
  ;NAME
  
  if ((n_elements(selfield) eq 0) OR (n_elements(sellevel) eq 0)) then begin
    numfields=n_elements(where(fieldhead(*,0)))
    fieldlist=reform(fieldhead(*,3))
    fieldlist=fieldlist(where(fieldlist))
    levellist=reform(fieldhead(*,5))
    levellist=levellist(where(levellist))
  endif else begin
    numfields=n_elements(selfield)
    fieldlist=selfield
    if (n_elements(sellevel) ne numfields) then begin
      print, 'ERROR - incorrect number of entries in level list'
      exit
    endif else begin
      levellist=sellevel
    endelse
  endelse 
  
endelse


;-----------------------------------------------------------------------

; loop through fields and species


for ifield=0,numfields-1 do begin 
 
  field=fieldlist(ifield)
  level=levellist(ifield)
 
  if (namever eq 3)then begin
    ; NAME III
    TimeAveraging=TimeAveragingList(ifield)
    FieldName=FieldNameList(ifield)
  endif  

  for ispec=0,numspecies-1 do begin

    species=specieslist(ispec)

;-----------------------------------------------------------------------     
      
; read data

    for it=0,numfiles-1 do begin

      filename=fieldfiles(it)
      filetime = timelist(it)

      if (namever eq 3)then begin
        ; NAME III

        selgrid1 = STRING(selgrid,'_C',casenum,'_T',it+1)
        selgrid1 = STRCOMPRESS(selgrid1, /Remove_All)

        readfield,datadir,filetime,selgrid1,field,level,$
	     species,lon1,lat1,dat1,head1,modtitle=modtitle,nrecs=nrecs,$
             namever=namever,timeframe=timeframe,AAR=1,xyrel=xyrel,$
             timeaveraging=timeaveraging,seldatavalue=seldatavalue,selname=FieldName

      endif else begin 
        filetime = timelist(it)

        readfield,datadir,filetime,selgrid,field,level,$
	species,lon1,lat1,dat1,head1,modtitle=modtitle,nrecs=nrecs,$
        namever=namever,timeframe=timeframe,AAR=1,xyrel=xyrel
      endelse
     
     
;-----------------------------------------------------------------------   
;calculate distance from source and eliminate points outside AAR 
;swap number sign for non area at risk fields

      n=long(0)
      datr=fltarr(1,6)

      for i=long(0),nrecs-1 do begin

        if selfield(ifield) eq 'Area At Risk' then begin
	  distance,d,lat1=xyrel(1),lon1=xyrel(0),lat2=dat1(i,3),lon2=dat1(i,2)
          dat1(i,5)=d
	endif else begin
	  dat1(i,5)=0.0-dat1(i,4)
	endelse
  
        if dat1(i,4) ge threshold then begin
          n=n+1
	  if (i eq 0) then begin
	    datr=dat1(i,*)
	  endif else begin
            datr=[datr,dat1(i,*)]
	  endelse
        endif

      endfor

      nrecs=n

;-----------------------------------------------------------------------
; sort data

      sdat=sort(datr(*,5))

;-----------------------------------------------------------------------
; adjust margin based on longitude west (if ydim not specified in
; input), then adjust xdim to return requested map size

      if(n_elements(ydim) eq 0)then begin
          xrescale = xdim*(1.0+(2.0*margin))
          margin = margin + ((xyrel(0)+2.0)/12.0)^2.0
          xdim =  xrescale/(1.0+(2.0*margin))
          if(AspectRatio eq 'Square')then begin
            ydim=xdim/1.0
          endif else begin
            ydim=xdim/1.4142
          endelse
      endif
     
;-----------------------------------------------------------------------
; do map calculation

      if n gt 0 then begin

        if(n_elements(sourceloc) ne 0)then xyrel=sourceloc
	
        if datamap eq 1 then begin
          lon1=min(datr(1:*,2))
          lat1=min(datr(1:*,3))
	  lon2=max(datr(1:*,2))
	  lat2=max(datr(1:*,3))
	  distance,xdim,lat1=lat1,lon1=lon1,lat2=lat1,lon2=lon2
	  distance,ydim,lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon1
        endif
      
        mapcalc,datr,sdat,xyrel,nrecs,xdim,ydim,margin,pacific,ll,tr,ie

        chemet_shape = "landscape"
      
        if aspect eq 0 then begin
	  mapcalc,datr,sdat,xyrel,nrecs,ydim,xdim,margin,pacific,ll2,tr2,ie2
	  if ie2 gt ie then begin
	    ll=ll2
	    tr=tr2
            chemet_shape = "portrait"
  	  endif
        endif  
       
;-----------------------------------------------------------------------
; output max and min longitudes and latitudes

        ; generate filename

        if(strpos(!version.os,'NT') gt 0)then begin
          pc=1
          delim='\'
        endif else begin
          pc=0
          delim='/'
        endelse

        aarfile=datadir+delim+'AaR_'+selgrid+'_'+specieslist(ispec)
        s_ifield = strcompress(string(ifield),/remove_all)
	aarfile=aarfile+'_'+filetime+'_'+s_ifield+'.txt'
    
        get_lun, unit
        openw, unit, aarfile
        printf,unit,'LONG_MIN=',strtrim(ll(0),2)
        printf,unit,'LONG_MAX=',strtrim(tr(0),2)
        printf,unit,'LAT_MIN=',strtrim(ll(1),2)
        printf,unit,'LAT_MAX=',strtrim(tr(1),2)
        close,unit
        free_lun, unit

        shape_file = datadir+delim+'AaR_shape.txt'

        get_lun, unit
        openw, unit, shape_file
        printf,unit,'shape=',chemet_shape
        close,unit
        free_lun, unit
	
      endif else begin
        print,'No data above threshold level for file:',filename
      endelse
 
    endfor
  endfor
endfor
;-----------------------------------------------------------------------     


end
