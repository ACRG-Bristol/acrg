pro plotmet,datadir,selgrid,selfield=selfield,zoom=zoom,filetext=filetext,$
 multiplot=multiplot,plotheader=plotheader,plotlegend=plotlegend,$
 genanim=genanim,contours=contours,group=group,$
 plotpmsl=plotpmsl,plotuv=plotuv,plottopog=plottopog,selcase=selcase,$
 exact=exact,gengif=gengif,genjpg=genjpg,last=last,namever=namever,polar=polar,$
 projection=projection,origin=origin,countries=countries,states=states
 
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
;  selfield     : (string or string array) field(s) to plot
;
;  zoom         : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  filetext     : (string)  text added to filenames
;  multiplot    : intarr(3) number of plots, 
;                 (eg multiplot=[0,2,2] for 2 by 2 plots)
;                 first element can be any value
;  plotheader   : integer (0 or 1) =1 adds run information and title to plots
;  plotlegend   : integer (0 or 1) =1 plots legend below each plot
;  genanim      : integer (0 or 1) =1 generates gifs and gif animation
;  contours     : fltarr(number of fields,*)
;                 array of contour values for each field  
;                 if not passed levels determined automatically
;  group        : string ('fields','time') determines what 
;                 order to plot fields
;  plotpmsl     : integer (0 or 1) =1 plots pressure contours
;  plotuv       : integer (0 or 1) =1 plots wind vectors
;  plottopog    : integer (0 or 1) =1 plots topography
;  exact        : integer (0 or 1) =1 map area based on lon/lat limits
;  gengif       : integer (0 or 1) =1 generates gif images (on unix)
;  genjpg       : integer (0 or 1) =1 generates jpg images (on unix)
;  last         : integer (0 or 1) =1 plots last timestep only
;  namever      : integer (2 or 3) defaults to 2
;  selcase      : case number (name iii only - defaults to 1 if not given)
;  projection   : integer (1 - 16) - projection for mapping (for options see
;                 IDL map option) #8 (cylindrical) is standard
;  polar        : integer, 0 for non-polar projection, -1 for south pole, 
;                 1 for north pole, projection modes 1 (stereographic)
;                 and 2 (orthographic) are recommended, but it will work with other
;                 projection modes.  The polar keyword overrides the origin keyword.
;  countries    : set to 0 to not plot national borders
;  states       : set to 0 to not plot individual states on maps of the USA
;  origin       : 2-element array [lat,long] to be used for origin of map projection.  
;                 If this keyword is not set, the automap option will be enabled and the
;                 map origin will be the centre of the map itself, (unless the polar 
;                 keyword is set).

;
;    
; example calls:

;  plotmet,'/home/fr1100/apnm/namev66/testv66',$
;   'grid2',group='time',selfield='Pmsl',$
;    multiplot=[0,2,2],plotheader=1,plotlegend=1,/genanim
;
;  plotmet,'/home/fr1100/apnm/namev66/testv66',$
;   'grid2',multiplot=[0,1,1],plotheader=1,plotlegend=1,/genanim
;
;
; subroutine calls:
;  readmet
;  readmethead
;  maplimits
;  genposition
;  plotfieldpos
;  genanim
;
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
; arguments

if(n_elements(exact) eq 0)then exact=1
if(n_elements(plottopog) eq 0)then plottopog=0
if(n_elements(plotpmsl) eq 0)then plotpmsl=0
if(n_elements(plotuv) eq 0)then plotuv=0
if(n_elements(group) eq 0)then group='time'
if(n_elements(plotlegend) eq 0)then plotlegend=1
if(n_elements(plotheader) eq 0)then plotheader=1
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(genanim) eq 0)then genanim=1
if(n_elements(gengif) eq 0) then gengif=0
if(n_elements(genjpg) eq 0) then genjpg=0
if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse
if(n_elements(last) eq 0)then last=0
if(n_elements(namever) eq 0)then namever=2
if(n_elements(selcase) eq 0)then begin
  selcase=1
  printcasenum=0
endif else begin
  printcasenum=1
endelse
if(n_elements(projection) eq 0) then projection=8
if(n_elements(polar) eq 0)then polar=0
if(n_elements(countries) eq 0) then countries=1
if(n_elements(states) eq 0) then states=1
if(n_elements(origin) eq 0) then begin
  automap=1
  origin=[0,0]
endif else if (n_elements(origin) ne 2) then begin
  print,'origin keyword must be specified as origin=[lat,long]'
  return
endif else begin 
  automap=0
  if (origin(0) lt -90) OR (origin(0) gt 90) then begin
    print,'origin latitude must be between -90 and +90'
    return
  endif else if (origin(1) lt -180) OR (origin(1) gt 360) then begin
    print, 'origin longitue must be between -180 and +360'
    return
  endif
endelse

;-----------------------------------------------------------------------
; graphics initialisation

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

xysize=1.0
axiscol=0
!p.font=0
!x.thick=1
!y.thick=1

nxplot=multiplot(1)
nyplot=multiplot(2)


if(plotheader ne 0)then begin
  position=[0.02,0.15,0.98,0.87]
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
  fieldfiles=findfile(datadir+delim+'Fields_'+selgrid+'_*')
endelse

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

readmethead,fieldfiles(0),modtitle,fieldhead,namever=namever

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
    timelist(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'.txt'))-12,12) 
  endelse
endfor

if(last ne 0)then begin
  timelist=timelist(numfiles-1)
  fieldfiles=fieldfiles(numfiles-1)
  numfiles=1
endif

; all fields

if(n_elements(selfield) eq 0) then begin
  if(namever eq 3) then begin
    ; NAME III
    fieldlist=reform(fieldhead(*,2))
  endif else begin
    ; NAME
    fieldlist=reform(fieldhead(*,3))
  endelse
  fieldlist=fieldlist(where(fieldlist))
  numfields=n_elements(fieldlist)
endif else begin
  numfields=n_elements(selfield)
  fieldlist=selfield
endelse 

;-----------------------------------------------------------------------
; plot time sequences

case group of 

'time': begin


;-----------------------------------------------------------------------
; loop through fields and species


for ifield=0,numfields-1 do begin 
 
  field=fieldlist(ifield)

;-----------------------------------------------------------------------     
; open postscript file for printing

  temp=strcompress(field,/remove_all)
  temp=strmid(temp,0,strpos(temp,'('))
  filetem='Plot_'+filetext+$
   temp+'_'+selgrid

  file=datadir+delim+filetem+'.ps'
  set_plot,'ps'
  device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
  device,/color
  device,filename=file
  loadct,5  
  stretch,32
  tek_color
 
     
;-----------------------------------------------------------------------        
; loop through times

  for it=0,numfiles-1 do begin

    filename=fieldfiles(it)
    filetime=timelist(it) 
    
    if (namever eq 3) then begin
      ; NAME III
      selgrid1 = STRING(selgrid,'_C',casenum,'_T',it+1)
      selgrid1 = STRCOMPRESS(selgrid1, /Remove_All)
    endif else begin
      ; NAME
      selgrid1=selgrid
    endelse
 
    readmet,datadir,filetime,selgrid1,field,lon1,lat1,dat1,modtitle,head1,namever=namever

    if(plotpmsl ne 0)then begin
      readmet,datadir,filetime,selgrid1,'Pmsl',lon2,lat2,pmsl,modtitle2,head2,namever=namever
    endif
    
    if(plotuv ne 0)then begin
      readmet,datadir,filetime,selgrid1,'U',lon2,lat2,u,modtitle2,head2,namever=namever
      readmet,datadir,filetime,selgrid1,'V',lon2,lat2,v,modtitle2,head2,namever=namever
    endif else begin
      u=0
      v=0
    endelse
    
    if(plottopog ne 0)then begin
      readmet,datadir,timelist(it),selgrid1,'Topog',lont,latt,topog,titt,headt,namever=namever
    endif else begin
      topog=0
    endelse
    
    pos=(it) mod (nxplot*nyplot)
    multiplot=[pos,nxplot,nyplot]

 ;-----------------------------------------------------------------------     
; plot fields

    metcontours,field,contours1,contourcolors,fieldunits
     
    if(n_elements(fieldunits) eq 0)then begin
      if (namever eq 3) then begin
        ; NAME III
	fieldunits=fieldhead(where(fieldhead(*,2) eq field),4)
      endif else begin
        ; NAME
        fieldunits=fieldhead(where(fieldhead(*,3) eq field),4)
      endelse 
    endif
     
    plotfieldpos,lon1,lat1,dat1, position=position,$
     zoom=zoom,contourlevels=contours1,contourcolors=contourcolors,$
     xyrel=xyrel,fieldunits=fieldunits,multiplot=multiplot,$
     modtitle=field,plotlegend=plotlegend,automap=automap,$
     pmsl=pmsl,topog=topog,u=u,v=v,exact=exact,polar=polar,$
     projection=projection,origin=origin,countries=countries,states=states

;-----------------------------------------------------------------------
; annotation    
   
    
    if(plotheader ne 0 and pos eq 0)then begin

      xyouts,0.5,0.98,modtitle(0),$
        /normal,alignment=0.5,charsize=1.3,color=axiscol 
      xyouts,0.5,0.94,modtitle(1),$
        /normal,alignment=0.5,charsize=1.3,color=axiscol
      if (namever eq 3)then begin
        ; NAME III
  	xyouts,0.5,0.90,'Valid at '+head1(13),$
	/normal,charsize=xysize*1.2,alignment=0.5,color=axiscol
      endif else begin
        ; NAME
        xyouts,0.5,0.90,'Valid at '+head1(6),$
	/normal,charsize=xysize*1.2,alignment=0.5,color=axiscol         
      endelse

      xyouts,0.05,0.12,modtitle(3),/normal,charsize=0.9,color=axiscol  
	xyouts,0.05,0.10,modtitle(2),/normal,charsize=0.9,color=axiscol
	xyouts,0.05,0.08,modtitle(4),/normal,charsize=0.9,color=axiscol
	xyouts,0.05,0.06,modtitle(5),/normal,charsize=0.9,color=axiscol
	
	xyouts,0.6,0.12,modtitle(7),/normal,charsize=0.9,color=axiscol
	xyouts,0.6,0.10,modtitle(8),/normal,charsize=0.9,color=axiscol
	xyouts,0.6,0.08,modtitle(6),/normal,charsize=0.9,color=axiscol


      xyouts,0.5,0.00,'Met Office Crown copyright',$
 	 /normal,alignment=0.5,charsize=0.8,color=axiscol 

      image=read_image('MO_Master_B.jpg')
      loadct,0
      tv,image,0.82,0.91,true=1,xsize=0.15,/normal
      loadct,5
      tek_color

    endif 


;-----------------------------------------------------------------------
    if(pos eq ((nxplot*nyplot)-1))then begin
      erase
    endif      

  endfor

  device,/close_file


;----------------------------------------------------------------------- 
; generate gifs


  if((gengif eq 1 or genjpg eq 1 or genanim eq 1) and pc eq 0)then begin
    genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
        genanim=genanim
  endif

endfor

end   

'field':  begin

numpages=((numfields-1)/(nxplot*nyplot))+1

;-----------------------------------------------------------------------
; loop through fields and species

; determine number of pages (if number of fields > number of plots per page)

for ipage=0,numpages-1 do begin

  filenum=string(ipage,format='(i2.2)') 
  startfield=nxplot*nyplot*ipage
  endfield=(nxplot*nyplot*(ipage+1))
  endfield=min([endfield,numfields])
 

;-----------------------------------------------------------------------     
; open postscript file for printing

  filetem='Plot_'+filetext+'met'+filenum+'_'+selgrid 
  file=datadir+'/'+filetem+'.ps'
  set_plot,'ps'
  device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
  device,/color
  device,filename=file
  loadct,5  
  stretch,32
  tek_color


;-----------------------------------------------------------------------     
; loop through times

  for it=0,numfiles-1 do begin

    filename=fieldfiles(it)
    filetime=timelist(it)
    if (namever eq 3) then begin
      ; NAME III
      selgrid1 = STRING(selgrid,'_C',casenum,'_T',it+1)
      selgrid1 = STRCOMPRESS(selgrid1, /Remove_All)
    endif else begin
      ; NAME
      selgrid1=selgrid
    endelse
   
;-----------------------------------------------------------------------        
; loop through fields and species

    for ifield=startfield,endfield-1 do begin 

      field=fieldlist(ifield)

      readmet,datadir,filetime,selgrid1,field,lon1,lat1,dat1,modtitle,head1,namever=namever

      if(plotpmsl ne 0)then begin
        readmet,datadir,filetime,selgrid1,'Pmsl',lon2,lat2,pmsl,modtitle2,head2,namever=namever
      endif else begin
        pmsl=0
      endelse
      
      if(plotuv ne 0)then begin
        readmet,datadir,filetime,selgrid1,'U',lon2,lat2,u,modtitle2,head2,namever=namever
        readmet,datadir,filetime,selgrid1,'V',lon2,lat2,v,modtitle2,head2,namever=namever
      endif else begin
        u=0
	v=0
      endelse
      
      if(plottopog ne 0)then begin
        readmet,datadir,timelist(it),selgrid1,'Topog',lont,latt,topog,titt,headt,namever=namever
      endif else begin
	topog=0
      endelse
       
      pos=(ifield) mod (nxplot*nyplot)
      multiplot=[pos,nxplot,nyplot]

;-----------------------------------------------------------------------     
; plot fields

      metcontours,field,contours1,contourcolors,fieldunits
      
      if(n_elements(fieldunits) eq 0)then begin
        if (namever eq 3) then begin
          ; NAME III
	  fieldunits=fieldhead(where(fieldhead(*,2) eq field),4)
        endif else begin
          ; NAME
          fieldunits=fieldhead(where(fieldhead(*,3) eq field),4)
        endelse 
      endif

      plotfieldpos,lon1,lat1,dat1, position=position,$
	zoom=zoom,contourlevels=contours1,contourcolors=contourcolors,$
	xyrel=xyrel,fieldunits=fieldunits,multiplot=multiplot,$
	modtitle=field,plotlegend=plotlegend,pmsl=pmsl,u=u,v=v,$
	topog=topog,exact=exact,polar=polar,automap=automap,$
	projection=projection,origin=origin,countries=countries,states=states


;-----------------------------------------------------------------------
; annotation    

      if(plotheader ne 0 and pos eq 0)then begin

        xyouts,0.5,0.98,modtitle(0),$
	 /normal,alignment=0.5,charsize=1.3,color=axiscol 
	xyouts,0.5,0.94,modtitle(1),$
	 /normal,alignment=0.5,charsize=1.3,color=axiscol
	if (namever eq 3)then begin
          ; NAME III
  	   xyouts,0.5,0.90,'Valid at '+head1(13),$
	   /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol
        endif else begin
          ; NAME
  	  xyouts,0.5,0.90,'Valid at '+head1(6),$
	  /normal,charsize=xysize*1.2,alignment=0.5,color=axiscol         
        endelse

	xyouts,0.05,0.12,modtitle(3),/normal,charsize=0.9,color=axiscol  
	xyouts,0.05,0.10,modtitle(2),/normal,charsize=0.9,color=axiscol
	xyouts,0.05,0.08,modtitle(4),/normal,charsize=0.9,color=axiscol
	xyouts,0.05,0.06,modtitle(5),/normal,charsize=0.9,color=axiscol
	
	xyouts,0.6,0.12,modtitle(7),/normal,charsize=0.9,color=axiscol
	xyouts,0.6,0.10,modtitle(8),/normal,charsize=0.9,color=axiscol
	xyouts,0.6,0.08,modtitle(6),/normal,charsize=0.9,color=axiscol


	xyouts,0.5,0.00,'Met Office Crown copyright',$
	    /normal,alignment=0.5,charsize=0.8,color=axiscol 

       image=read_image('MO_Master_B.jpg')
       loadct,0
       tv,image,0.82,0.91,true=1,xsize=0.15,/normal
       loadct,5
       tek_color

      endif 


;-----------------------------------------------------------------------
; new page

      if(pos eq ((nxplot*nyplot)-1) or ifield eq endfield-1)then begin
         erase
      endif      

     endfor  
   endfor    

   device,/close_file
   set_plot,'X'


;----------------------------------------------------------------------- 
; generate gifs
 
   if((gengif eq 1 or genjpg eq 1 or genanim eq 1) and pc eq 0)then begin
     genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
        genanim=genanim
   endif
  
 endfor


end   

endcase

end
