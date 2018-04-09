pro datathreshold,datadir,selgrid,selspecies=selspecies,$
    selfield=selfield,sellevel=sellevel,filetext=filetext,template=template,$
    scalefactor=scalefactor,namever=namever,SelTimeAverage=SelTimeAverage,$
    selcase=selcase,seldatavalue=seldatavalue,selname=selname,plotvaacobs=plotvaacobs,$
    selmask=selmask,threshold=threshold,missing=missing
    
  ;--------------------------------------------------------------------
  ; Procedure to set data points to zero when masking field below threshold level 
  ;
  ; LH 06/05/2010 
  ;
  ; Changes by HJC 13/01/2011:
  ; 1. Add an input variable 'missing', used for missing fields in the NAME
  ; output file, which defaults to zero, but if set to 1, it makes the scale
  ; factor zero and then sets the field to -1.0 everywhere. 
  ; 2. Bug fix - numlevels is not defined for NAME-III, only for NAME-II, so
  ; I'm adding it to the NAME-III loop as it gives errors otherwise. 
  ;
  ; Adapted from plotfield.pro
  ; DBR Feb 2002
  ;
  ; Arguments
  ; required:
  ;  datadir      : (string)  run directory
  ;  selgrid      : (string)  grid ('grid1' or 'grid2')
  ;  selfield     : (string or string array) field(s) to overwrite
  ;  selmask      : (string or srting array) field(s) to use as mask
  ;                 must contain same number of elements as selfield
  ;
  ; optional:
  ;  selspecies   : (string or string array) species(s) to overwrite
  ;  sellevel     : (string or string array) level(s) to overwrite
  ;
  ;  SelTimeAverage  : (string) time averaged field to overwrite
  ;
  ;  filetext     : (string)  text added to filenames
  ;  template     : string - string to include when listing Fields* files
  ;                 eg ='0000' to search for midnight data only
  ;  scalefactor  : real - scaling factors to apply to data
  ;  namever      : integer indicating version of model e.g. 3 - NAME III(PPM)
  ;  selcase      : case number (name iii only - defaults to 1 if not given)
  ;  seldatavalue : data value (i.e. percentile or probability)
  ;                 (name iii only - defaults to ' ' if not given)
  ;  selname      : name of field requirement (can be used to specify
  ;                 a field which is not uniquely defined by the other
  ;                 parameters; name iii only - defaults to ' ' if not given)
  ;  plotvaacobs      : string filename conatining observations to overplot,
  ;                 see plotvaacobs.pro for file format
  ;  threshold    : level to use as threshold for masking field
  ;                 defaults to 10^-17 (for VAAC)
  ;  missing      : default is zero; set to 1 for all values to be set to -1.0
  ;
  ;-----------------------------------------------------------------------
  ; For a more full description of NAME III output format and the appropriate
  ; NAME III options for standard idl plotting please refer to the NAME III
  ; document: md2_4_v2(Output).
    
    
  ;-----------------------------------------------------------------------
  ; arguments
    
  if(n_elements(plotvaacobs) eq 0)then plotvaacobs=''
  if(n_elements(titletext) eq 0)then titletext=''
  if(n_elements(filetext) eq 0)then begin
    filetext=''
  endif else begin
    filetext=filetext+'_'
  endelse
  if(n_elements(template) eq 0)then template='*'
  if(n_elements(scalefactor) eq 0)then scalefactor=1.0
  if(n_elements(namever) eq 0)then namever=2
  if(n_elements(selcase) eq 0)then begin
    selcase=1
    printcasenum=0
  endif else begin
    printcasenum=1
  endelse
  ;if(n_elements(seldatavalue) eq 0)then seldatavalue=' '
  if(n_elements(threshold) eq 0) then threshold=1.0E-17
  if(n_elements(missing) eq 0)then missing=0
; to flag field as missing, first set it all to zero
  if (missing eq 1) then scalefactor=0.0
  
  ;-----------------------------------------------------------------------
  ; find all files and readheader
  
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
    namever=namever
    
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
  
  fieldlist=selfield
  mask=selmask
  numfields=n_elements(selfield)
  
  if (namever eq 3)then begin
    ; NAME III
  
    if (n_elements(sellevel) eq 0) then begin
      levellist=reform(fieldhead(*,14))
      levellist=levellist(0:n_elements(levellist)-2)
    endif else begin
      levellist=sellevel
    endelse
    numlevels=n_elements(levellist)
    
    if (n_elements(SelTimeAverage) eq 0) then begin
      timeaveraginglist=reform(fieldhead(*,7))
      timeaveraginglist=timeaveraginglist(0:n_elements(timeaveraginglist)-2)
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
    
    if (n_elements(seldatavalue) eq 0) then begin
      datavaluelist=replicate(' ',numfields)
    endif else if (n_elements(seldatavalue) ne numfields) then begin
      print,'ERROR - incorrect number of entries in seldatavalue list'
      exit
    endif else begin
      datavaluelist=seldatavalue
    endelse
    
  endif else begin
    ;NAME
  
    numfields=n_elements(selfield)
    fieldlist=selfield
  
    if (n_elements(sellevel) eq 0) then begin
      levellist=reform(fieldhead(*,5))
      levellist=levellist(0:n_elements(levellist)-2)
    endif else begin
      levellist=sellevel
    endelse
    numlevels=n_elements(levellist)
    
  endelse
  

  ;-----------------------------------------------------------------------
  ; loop through times

  for it=0,numfiles-1 do begin

    filename=fieldfiles(it)
    filetime = timelist(it)
  
    if (namever eq 3)then begin
      ; NAME III    
      selgrid1 = STRING(selgrid,'_C',casenum,'_T',it+1)
      selgrid1 = STRCOMPRESS(selgrid1, /Remove_All)  
    endif 
    
    ;-----------------------------------------------------------------------
    ; loop through species, fields and levels
        
    for ispec=0,numspecies-1 do begin
      
      species=specieslist(ispec)
  
      for ifield=0,numfields-1 do begin
    
        field=fieldlist(ifield)
        mask=selmask(ifield)
        
        if (namever eq 3)then begin
          ; NAME III
          TimeAveraging=TimeAveragingList(ifield)
          FieldName=FieldNameList(ifield)
          datavalue=datavaluelist(ifield)
        endif
        
        for ilevel=0,numlevels-1 do begin
          
          level=levellist(ilevel)
                  
          if (namever eq 3)then begin
            ; NAME III
            
            readfield,datadir,filetime,selgrid1,field,level,$
              species,lon1,lat1,dat1,head1,modtitle=modtitle,$
              namever=namever,timeaveraging=timeaveraging,$
              seldatavalue=datavalue,selname=FieldName,AaR=1
              
            readfield,datadir,filetime,selgrid1,mask,level,$
              species,lon1,lat1,dat2,head1,modtitle=modtitle,$
              namever=namever,timeaveraging=timeaveraging,$
              seldatavalue=datavalue,selname=FieldName,AaR=1
              
          endif else begin
            ; NAME
            
            readfield,datadir,filetime,selgrid,field,level,$
              species,lon1,lat1,dat1,head1,modtitle=modtitle,$
              namever=namever,AaR=1
            
            readfield,datadir,filetime,selgrid,mask,level,$
              species,lon1,lat1,dat2,head1,modtitle=modtitle,$
              namever=namever,AaR=1
              
          endelse
          
          dat1=dat1(*,4)
          dat2=dat2(*,4)
          
          index=where(dat2 lt threshold,ntot)
          if(ntot gt 0) then dat1(index)=0.0
          
          dat1=dat1*scalefactor
          
          index=where(dat1 lt 1e-9*max(dat1),ntot)
          if(ntot gt 0)then dat1(index)=0.0
          
          ; set all values to -1.0 if missing flag is set
          
          if (missing eq 1) then dat1=dat1-1.0
          
          if (namever eq 3)then begin
            ; NAME III
            
            writefield,datadir,filetime,selgrid1,field,level,$
              species,lon1,lat1,dat1,head1,modtitle=modtitle,$
              namever=namever,timeaveraging=timeaveraging,$
              seldatavalue=datavalue,selname=FieldName
              
          endif else begin
            ; NAME
            
            writefield,datadir,filetime,selgrid,field,level,$
              species,lon1,lat1,dat1,head1,modtitle=modtitle,$
              namever=namever
              
          endelse
        
        endfor
      endfor
    endfor 
  endfor 

end
