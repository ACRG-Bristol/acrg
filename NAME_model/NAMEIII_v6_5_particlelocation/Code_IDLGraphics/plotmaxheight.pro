pro plotmaxheight,datadir,selgrid,selfield,selspecies,$
    sellevel=sellevel,SelTimeAverage=seltimeaverage,$
    threshold=threshold,contourcols=contourcols,$
    namever=namever,zoom=zoom,filetext=filetext,projection=projection,$
    multiplot=multiplot,plotheader=plotheader,plotlegend=plotlegend,$
    type=type,gengif=gengif,genanim=genanim

    
  ;--------------------------------------------------------------------
  ; Procedure to generate sequence of plots of the maximum height
  ; at which the selected threshold is exceeded.
  ; Based on plotfield.pro and developed for VAAC volcanic ash 
  ;
  ;
  ; Claire Witham Feb 2011
  ;
  ; Arguments
  ; required:
  ;  datadir      : (string)  run directory
  ;  selgrid      : (string)  grid ('grid1' or 'grid2')
  ;  selfield     : (string or string array) field(s) to plot
  ;  selspecies   : (string or string array) species(s) to plot
  ;  sellevel     : (string or string array) level(s) to plot
  ;  selTimeAverage : (string or string array) time average to plot - NAME III
  ;  threshold    : real number specifying the threshold for air concentration or
  ;                 other required field in the units of the NAME Field files
  ;  contourcols  : (integer array) colours for the levels, based on ct=4
  ;  namever      : integer indicating version of model e.g. 3 - NAME III(PPM)
  ;  zoom         : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
  ;  filetext     : (string) optional text added to filenames
  ;  projection   : integer (1 - 16) - projection for mapping (for options see
  ;                 IDL map option) #8 (cylindrical) is standard
  ;  multiplot    : intarr(3) number of plots,
  ;                 (eg multiplot=[0,2,2] for 2 by 2 plots)
  ;                 first element can be any value
  ;  plotheader   : integer (0 or 1) =1 adds run information and title to plots
  ;  plotlegend   : integer (0 or 1) =1 plots legend below each plot
  ;  type         : string ('pixel' or 'contour') choose pixels or
  ;                  filled contours
  ;  gengif       : integer (0 or 1) =1 generates gif images (on unix)
  ;  genanim      : integer (0 or 1) =1 generates gifs animation
    
  ;-----------------------------------------------------------------------
  ; arguments
 
  if(n_elements(type) eq 0)then type='pixel'
  if(n_elements(titletext) eq 0)then titletext=''    
  if(n_elements(plotlegend) eq 0)then plotlegend=1
  if(n_elements(plotheader) eq 0)then plotheader=1
  if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
  if(n_elements(genanim) eq 0)then genanim=0
  if(n_elements(gengif) eq 0) then gengif=0
  if(n_elements(filetext) eq 0)then begin
    filetext=''
  endif else begin
    filetext=filetext+'_'
  endelse
  if(n_elements(projection) eq 0) then projection=8
  if(n_elements(template) eq 0)then template='*'
  if(n_elements(printtimeaverage) eq 0) then printtimeaverage=0
  if(n_elements(selcase) eq 0)then begin
    selcase=1
    printcasenum=0
  endif else begin
    printcasenum=1
  endelse   
  if(n_elements(metoffice) eq 0) then metoffice=1

  ;-----------------------------------------------------------------------
  ; Check correct arguments
  
  if (n_elements(sellevel) ne n_elements(contourcols)) then begin
      print,'ERROR - number of levels does not match number of contour colours'
      exit  
  endif
    
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
  
  xysize=1.0
  axiscol=0
  !p.font=0
  !x.thick=1
  !y.thick=1
  
  nxplot=multiplot(1)
  nyplot=multiplot(2)
  
  if(plotheader ne 0)then begin
    position=[0.00,0.15,1.00,0.87]
  endif else begin
    position=[0.0,0.0,1.0,1.0]
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
    namever=namever,timeframe=timeframe
    
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
  
  if n_elements(numfields) eq 0 then begin
    numfields=n_elements(where(fieldhead(*,0)))
  endif
 
  
  if (namever eq 3)then begin
    ; NAME III
  
    if ((n_elements(selfield) eq 0) OR (n_elements(sellevel) eq 0)) then begin
      fieldlist=reform(fieldhead(*,2))
      fieldlist=fieldlist(0:n_elements(fieldlist)-2)
      levellist=reform(fieldhead(*,14))
      levellist=levellist(0:n_elements(levellist)-2)
    endif else begin
      numfields=n_elements(selfield)
      fieldlist=selfield
      numlevels=n_elements(sellevel)
      levellist=sellevel
    endelse
    
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
    
  endif else begin
    ;NAME
  
    if ((n_elements(selfield) eq 0) OR (n_elements(sellevel) eq 0)) then begin
      numfields=n_elements(where(fieldhead(*,0)))
      fieldlist=reform(fieldhead(*,3))
      fieldlist=fieldlist(0:n_elements(fieldlist)-2)
      levellist=reform(fieldhead(*,5))
      levellist=levellist(0:n_elements(levellist)-2)
    endif else begin
      numfields=n_elements(selfield)
      fieldlist=selfield
      numlevels=n_elements(sellevel)      
      levellist=sellevel
    endelse
    
  endelse
    
  
  ;-----------------------------------------------------------------------
  ; map limits
  
  ; this extra read is needed when ContourGrid  is being used
  ; it ensures that the correct map extant is used.
  
  field=fieldlist(0)
  level=levellist(0)
  species=specieslist(0)
  
  if (namever eq 3)then begin
    ; NAME III
  
    TimeAveraging=TimeAveragingList(0)
;    FieldName=FieldNameList(0)
;    datavalue=datavaluelist(0)
    
    selgrid1 = selgrid+'_C'+casenum+'_T1'
    
    readfield,datadir,timelist(0),selgrid1,field,level,$
      species,lon1,lat1,dat1,head1,namever=namever,timeframe=timeframe,$
      TimeAveraging=TimeAveraging,seldatavalue=datavalue,selname=FieldName
      
  endif else begin
    ;NAME
  
    readfield,datadir,timelist(0),selgrid,field,level,$
      species,lon1,lat1,dat1,head1,namever=namever,timeframe=timeframe
    
  endelse
  
  if(n_elements(zoom) eq 0) then begin
    zoom=fltarr(4)
    maplimits,lon1,lat1,dat1,xmin,xmax,ymin,ymax,exact=exact
    zoom(0)=xmin
    zoom(1)=xmax
    zoom(2)=ymin
    zoom(3)=ymax
  endif else if (n_elements(zoom) ne 0) then begin
    xmin=zoom(0)
    xmax=zoom(1)
    ymin=zoom(2)
    ymax=zoom(3)
  endif
  

  ;-----------------------------------------------------------------------
  ; plot time sequences

      ;-----------------------------------------------------------------------
      ; loop through fields and species
    
    
      for ifield=0,numfields-1 do begin
      
        field=fieldlist(ifield)
        
        if (namever eq 3)then begin
          ; NAME III
          TimeAveraging=TimeAveragingList(ifield)
;          FieldName=FieldNameList(ifield)
;          datavalue=datavaluelist(ifield)
        endif
        
        for ispec=0,numspecies-1 do begin
        
          species=specieslist(ispec)
          
          ;-----------------------------------------------------------------------
          ; open postscript file for printing
          
          if (namever eq 3)then begin
            ; NAME III
            filetem='Plot_'+filetext+$
              'MaxHeight_'+$
              'FromFL000-FL400_'+$
              selgrid+'_'+species
          endif else begin
            filetem='Plot_'+filetext+$
              'MaxHeight_'+$
              'FromFL000-FL400_'+$
              selgrid+'_'+species
          endelse
          
          
          file=datadir+delim+filetem+'.ps'
          set_plot,'ps'
          
          device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
          device,/color
          device,filename=file
          loadct,4
          tek_color
          
          
          ;-----------------------------------------------------------------------
          ; loop through times
          
          for it=0,numfiles-1 do begin
          
            filename=fieldfiles(it)
 
            ;-----------------------------------------------------------------------
            ; loop through levels
          
            for il=0,numlevels-1 do begin

                level=levellist(il)

                if (namever eq 3) then begin
                  ; NAME III
                  levelnum=float(strmid(level,3,5))
                endif else begin
                   ; NAME                 
                  levelnum=float(strmid(level,15,3))
                endelse
             
                if (namever eq 3)then begin
                  ; NAME III
                  filetime = timelist(it)
              
                  selgrid1 = STRING(selgrid,'_C',casenum,'_T',it+1)
                  selgrid1 = STRCOMPRESS(selgrid1, /Remove_All)
              
                  readfield,datadir,filetime,selgrid1,field,level,$
                    species,lon1,lat1,dat1,head1,modtitle=modtitle,$
                    namever=namever,timeframe=timeframe,$
                    timeaveraging=timeaveraging,seldatavalue=datavalue,selname=FieldName
                
                endif else begin
                   ; NAME
                   filetime = timelist(it)
              
                   readfield,datadir,filetime,selgrid,field,level,$
                     species,lon1,lat1,dat1,head1,modtitle=modtitle,$
                     namever=namever,timeframe=timeframe
                
                endelse
                   
                ; NEW processing here for calculating max height
                ; This could probably be made more efficient

                ; If first level, read in key values           
                if (il eq 0) then begin
                  datout=dat1
                  datout(*)=-1.0
                  nx=n_elements(dat1[*,0])
                  ny=n_elements(dat1[0,*])
                endif 

                ; loop over all x and y values to determine where grid cells are above threshold
                For y=0,ny-1 do begin
                  For x=0,nx-1 do begin
                     ; if value above threshold set to that levels 
                     if (dat1[x,y] ge threshold) then begin
                        datout[x,y]=levelnum-1.0   ; subtract 1.0 to take it below the value used for contouring
                     endif                   
                  Endfor
                Endfor

             endfor
                    
             pos=(it) mod (nxplot*nyplot)
             multiplot=[pos,nxplot,nyplot]
            
            
            ;-----------------------------------------------------------------------
            ; plot fields
         
            ; Determine the number of contours and the corresponding contour text 
            ; from the number of levels requested
                        
            ncont=n_elements(levellist)
            contours1=fltarr(ncont+1)           

            if (namever eq 3) then begin
              ; NAME III
              contours1[0]=0.00
              for l=0,ncont-1 do begin
                contours1[l+1]=float(strmid(levellist[l],3,5))
              endfor                        
            endif else begin
              ;NAME          
              contours1[0]=float(strmid(levellist[0],7,3))
              for l=0,ncont-1 do begin
                contours1[l+1]=float(strmid(levellist[l],15,3))
              endfor
            endelse
            
            ; For the 8 standard levels contours1=[0., 50., 100., 150., 200., 250., 300.,350., 400.]

              if (namever eq 3)then begin
                ; NAME III
                if (printtimeaverage ne 0) then begin
                  timeavstring,head1(13),head1(7),postitle
                endif else begin
                  postitle='Valid at '+head1(13)
                endelse
              endif else begin
                ; NAME
                if (printtimeaverage ne 0) then begin
                  timeavstring,head1(6),head1(2),postitle
                endif else begin
                  postitle='Valid at '+head1(6)
                endelse
              endelse
            
            
            plotfieldpos_maxheight,lon1,lat1,datout, position=position,type=type,$
              zoom=zoom,contourlevels=contours1,contourcolors=contourcols,$
              xyrel=xyrel,fieldunits=head1(4),multiplot=multiplot,$
              modtitle=postitle,plotloc=plotloc,$
              plotlegend=plotlegend,pmsl=pmsl,u=u,v=v,topog=topog,$
              highres=highres,projection=projection,tlon=tlon,tlat=tlat,$
              plotlincomtopog=plotlincomtopog,polar=polar,$
              countries=countries,coasts=coasts,states=states,$
              placenames=placenames,origin=origin,automap=automap,plotvaacobs=plotvaacobs
              
              
            ;-----------------------------------------------------------------------
            ; annotation
              
            if(plotheader ne 0 and pos eq 0)then begin
            
              xyouts,0.5,0.98,strtrim(modtitle(0),2),$
                /normal,alignment=0.5,charsize=1.3,color=axiscol
                
              xyouts,0.5,0.94,strcompress(modtitle(1))+titletext,$
                /normal,alignment=0.5,charsize=1.3,color=axiscol
               
              if (namever eq 3)then begin
                ; NAME III
                xyouts,0.5,0.90,head1(7),$
                  /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol
                  
                xyouts,0.5,0.873,'Maximum Height of Specified Threshold: '+string(threshold),$
                  /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol
                  
              endif else begin
                ; NAME
              
                   xyouts,0.5,0.90,head1(2),$
                    /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol                   
                    
                  xyouts,0.5,0.873,'Maximum Height of Specified Threshold: '+string(threshold),$
                    /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol
                
              endelse
              
              
              xyouts,0.05,0.12,modtitle(4),/normal,charsize=0.9,color=axiscol
              xyouts,0.05,0.10,modtitle(5),/normal,charsize=0.9,color=axiscol
              xyouts,0.05,0.08,modtitle(6),/normal,charsize=0.9,color=axiscol
              xyouts,0.05,0.06,modtitle(7),/normal,charsize=0.9,color=axiscol
              xyouts,0.05,0.04,modtitle(8),/normal,charsize=0.9,color=axiscol
              
              if (namever eq 3)then begin
                ; NAME III
                xyouts,0.60,0.12,'Pollutant: '+head1(3),$
                  /normal,charsize=0.9,color=axiscol
                ;
;                if (datavaluelist(ifield) ne ' ')then begin
;                  ; NAME III data value
;                  tokens=strsplit(datavaluelist(ifield), '=', /extract)
;                  dvalue=strtrim(tokens(1), 2)
;                  xyouts,0.60,0.02,'Data value: '+dvalue,$
;                    /normal,charsize=0.9,color=axiscol
;                endif
              endif else begin
                ; NAME
                xyouts,0.60,0.12,'Pollutant: '+head1(1),$
                  /normal,charsize=0.9,color=axiscol
              endelse
              
              xyouts,0.6,0.10,modtitle(3),/normal,charsize=0.9,color=axiscol
              xyouts,0.6,0.08,modtitle(2),/normal,charsize=0.9,color=axiscol
              
              
              if (printcasenum eq 1)then begin
                ; NAME III case number
                xyouts,0.60,0.04,'Case number: '+casenum,$
                  /normal,charsize=0.9,color=axiscol
              endif

              if (metoffice eq 1) then begin
                xyouts,0.5,0.00,'Met Office Crown copyright',$
                  /normal,alignment=0.5,charsize=0.8,color=axiscol
                  
                image=read_image('MO_Master_B.jpg')
                loadct,0
                tv,image,0.82,0.91,true=1,xsize=0.15,/normal
                loadct,4
                tek_color

              endif else if (metoffice eq 2) then begin
                xyouts,0.5,0.00,'Met Office Crown copyright',$
                  /normal,alignment=0.5,charsize=0.8,color=axiscol
                  
                image=read_image('MO_Master_Mono_B.jpg')
                loadct,0
                tv,image,0.82,0.91,true=1,xsize=0.15,/normal
                loadct,4
                tek_color

              endif
              
            endif
            
            ;-----------------------------------------------------------------------
            if(pos eq ((nxplot*nyplot)-1))then begin
              erase
            endif
            
          endfor
          
          device,/close_file
          
          
          ;-----------------------------------------------------------------------
          ; generate gifs
          
          
          if((gengif eq 1 or genanim eq 1) and pc eq 0)then begin
          
            genanim,file,datadir,filetem,gengif=gengif,$
              genanim=genanim,method=method,xpixel=xpixel
              
          endif
          
        endfor
      endfor
      

  
end
