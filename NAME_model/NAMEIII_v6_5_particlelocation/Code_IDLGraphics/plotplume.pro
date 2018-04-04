pro plotplume,datadir,zoom=zoom,filetext=filetext,$
  zrange=zrange,boundarylayer=boundarylayer,$
  multiplot=multiplot,plotheader=plotheader,plotlegend=plotlegend,$
  genanim=genanim,plumecolour=plumecolour,cross_section=cross_section,$
  plumeres=plumeres,plottopog=plottopog,namever=namever,labelgrid=labelgrid,$
  plotuv=plotuv,plotpmsl=plotpmsl,selgrid=selgrid,projection=projection,$
  gengif=gengif,genjpg=genjpg,xpixel=xpixel,plotqfe=plotqfe,polar=polar,$
  origin=origin,countries=countries,coasts=coasts,states=states,$
  automap=automap,plotloc=plotloc,placenames=placenames,$
  usecontours=usecontours,plotfl=plotfl,ncolors=ncolors

;--------------------------------------------------------------------
; Procedure to generate sequence of plume plots
; 
;
; DBR March 2002
;
;
; Arguments
;
; required:
;  datadir      : (string)  run directory
;
; optional:
;  zoom         : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  zrange       : fltarr(2) height range (eg zrange=[0,10000]
;                  plots only particles within this height band and
;                  scales crss sections to this height
;  boundarylayer: integer (0 or 1) =1 plots only boundary layer particles
;                 note zrange still apllies 
;  filetext     : (string)  text added to filenames
;  multiplot    : intarr(3) number of plots, 
;                 (eg multiplot=[0,2,2] for 2 by 2 plots)
;                 first element can be any value
;  plotheader   : integer (0 or 1) =1 adds run information and title to plots
;  plotlegend   : integer (0 or 1) =1 plots legend below each plot
;  genanim      : integer (0 or 1) =1 generates gifs and gif animation
;  plumecolour  : string ('age','height','bl') selects colouring 
;                 of particles
;                 'bl' colours boundary layer particles in red 
;  cross_section: integer (0,1,2 or 3) =1,2,3 plots a cross section of the plume
;                 particles selected according to zoom, zrange and
;                 boundary layer. Note legend not plotted if coloured 
;                 by 'age' or 'height'
;                 =1 Choses X or Y slice depending on X and Y ranges
;                 =2 X (longitude) slice
;                 =3 Y (latitude) slice
;  plumeres     : integer  selects every plumeres particle for plotting
;               : eq if =10 plots one in every 10 particles 
;  plotpmsl     : integer (0 or 1) =1 plots pressure contours
;  plotuv       : integer (0 or 1) =1 plots wind vectors
;  selgrid      : dtring - grid for met data (pmsl,u,v)
;  plottopog    : integer (0 or 1) =1 plots topography
;  gengif       : integer (0 or 1) =1 generates gif images (on unix)
;  genjpg       : integer (0 or 1) =1 generates jpg images (on unix)
;  xpixel       : size of gif/jpg images in call to genanim
;  plotqfe      : integer (0 or 1) =1 to use zqfe, otherwise uses zqnh
;  labelgrid    : integer (0 or 1) =1 to label gridlines (0 default) 
;  projection   : integer (1 - 16) - projection for mapping (for options see
;                 IDL map option) #8 (cylindrical) is standard
;  polar        : integer, 0 for non-polar projection, -1 for south pole, 
;                 1 for north pole, projection modes 1 (stereographic)
;                 and 2 (orthographic) are recommended, but it will work with other
;                 projection modes.  The polar keyword overrides the origin keyword.
;  origin       : 2-element array [lat,long] to be used for origin of map projection.  
;                 If this keyword is not set, the automap=-1 option will be enabled
;                 (unless 'polar' is also set).
;  automap      : (default 0) set to 1 to get map projection with origin in centre of grid
;                 or -1 for map centred by longitude only.  If an orgin has been set, 
;                 automap will override it, unless 'polar' is set to -1 or +1.
;  countries    : set to 0 to not plot national borders
;  coasts       : set to 1 to plot lakes and extra islands.  N.B. do not use this
;                 with the 'hires' keyword on continent-sized maps, and do not use
;                 it without 'hires' set on anything smaller (however it will plot
;                 every lake, however small).  You will need to use it without 'hires' 
;                 for getting the Great Lakes on maps of the whole of N. America, and 
;                 the Caspian Sea on large-scale maps of Asia, and with 'hires' to
;                 show the smaller islands around the UK.
;  states       : set to 0 to not plot individual states on maps of the USA
;  plotloc      : string filename conatining locations to overplot,
;                 see plotloc.pro for file format
;  placenames   : set to 1 to plot names of locations in plotloc file
;  usecontours  : fltarr(number of species,*) array of contour values for each field
;  plotfl       : set to 1 if plotting flight levels  

; example calls:
;  plotplume,datadir,multiplot=[0,3,4],plotheader=1,plotlegend=1,$
;   /genanim,plumeres=10,zoom=[-8,8,48,52],zrange=[0,2000],/cross_section,$
;   plumecolour='bl'
;
;  plotplume,datadir,multiplot=[0,1,1],plotheader=1,plotlegend=1,$
;   /genanim,zrange=[0,10000],plumecolour='age'
;
;
; subroutine calls:
;  readtraj
;  maplimits
;  genposition
;  genanim
;  readmet
;
;
; Changes
;  29/05/2002 DBR  Projection and stretch changes
;----------------------1-------------------------------------------------



;-----------------------------------------------------------------------
; arguments

if(n_elements(plottopog) eq 0)then plottopog=0
if(n_elements(plotlegend) eq 0)then plotlegend=1
if(n_elements(plotheader) eq 0)then plotheader=1
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(genanim) eq 0)then genanim=1
if(n_elements(gengif) eq 0) then gengif=0
if(n_elements(genjpg) eq 0) then genjpg=0
if(n_elements(xpixel) eq 0) then xpixel=595

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse
if(n_elements(plumecolour) eq 0)then plumecolour='age'
if(n_elements(boundarylayer) eq 0)then boundarylayer=0
if(n_elements(cross_section) eq 0)then cross_section=0
if(n_elements(plumeres) eq 0)then plumeres=1
if(n_elements(plotpmsl) eq 0)then plotpmsl=0
if(n_elements(plotuv) eq 0)then plotuv=0
if(n_elements(selgrid) eq 0)then selgrid='grid1'
if(n_elements(plotqfe) eq 0)then plotqfe=1
if(n_elements(labelgrid) eq 0) then labelgrid=0
if(n_elements(namever) eq 0)then namever=2
if(n_elements(polar) eq 0) then polar=0
if(n_elements(countries) eq 0) then countries=1
if(n_elements(coasts) eq 0) then coasts=0
if(n_elements(states) eq 0) then states=1
if(n_elements(origin) eq 0) then begin
  origin=[0,0]
  if (n_elements(automap) eq 0) then automap=-1
endif else if (n_elements(origin) ne 2) then begin
  print,'origin keyword must be specified as origin=[lat,long]'
  return
endif else begin
  if (origin(0) lt -90) OR (origin(0) gt 90) then begin
    print,'origin latitude must be between -90 and +90'
    return
  endif else if (origin(1) lt -180) OR (origin(1) gt 360) then begin
    print, 'origin longitude must be between -180 and +360'
    return
  endif
endelse
if(n_elements(automap) eq 0) then automap=0
if(n_elements(plotloc) eq 0)then plotloc=''
if(n_elements(placenames) eq 0) then placenames=0
if(n_elements(plotfl) eq 0) then plotfl=0

  ;if(n_elements(ncolors) ne 0) then ncolors=ncolors+2
  if(n_elements(ncolors) eq 0) then ncolors=6

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

fieldfiles=findfile(datadir+delim+'Plume_'+'*')

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

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif


; all times

timelist=strarr(numfiles)
for it=0,numfiles-1 do begin
  timelist(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Plume'))+6,12) 
endfor

;-----------------------------------------------------------------------   
; map limits

; read first and last

readtraj,fieldfiles(0),modtitle,date1,xlon1,ylat1,temp,press,ptemp,bldepth,$
         cloud,rh,wspd,wdir,zqfe1,zqnh1,zp1,zfl1,age1,namever=namever
readtraj,fieldfiles(numfiles-1),modtitle,date2,xlon2,ylat2,temp,press,ptemp,$
         bldepth,cloud,rh,wspd,wdir,zqfe2,zqnh2,zp2,zfl2,age2,namever=namever	     

if(n_elements(zoom) eq 0)then begin

  lon1=[xlon1(0:n_elements(xlon1)-1),xlon2(0:n_elements(xlon2)-1)]
  lat1=[ylat1(0:n_elements(ylat1)-1),ylat2(0:n_elements(ylat2)-1)]
  tot=lon1*0.0+1.0
  zoom=fltarr(4)
  maplimits,lon1,lat1,tot,xmin,xmax,ymin,ymax
  zoom(0)=xmin
  zoom(1)=xmax
  zoom(2)=ymin
  zoom(3)=ymax
endif else begin
  xmin=zoom(0)
  xmax=zoom(1)
  ymin=zoom(2)
  ymax=zoom(3)
endelse

; vertical range

if(n_elements(zrange) eq 0)then begin
  zqfe1=zqfe1(0:n_elements(zqfe1)-1)
  zqfe2=zqfe2(0:n_elements(zqfe2)-1)
  zrange=[0,max([zqfe1,zqfe2])]
  if (plotfl ne 0) then zrange=[0,max(zfl1,zfl2)]  
endif


; age range

if(n_elements(agerange) eq 0)then begin
  age1=age1(0:n_elements(age1)-1)
  age2=age2(0:n_elements(age2)-1)
  agerange=[min([age1,age2]),max([age1,age2])]
endif

;-----------------------------------------------------------------------   
; grid lines
  
dxy=ymax-ymin

case 1 of 
  (dxy ge 90): gridxy=30.0
  (dxy ge 50) and (dxy lt 90): gridxy=20.0
  (dxy ge 20) and (dxy lt 50): gridxy=10.0
  (dxy ge 10) and (dxy lt 20): gridxy=5.0
  (dxy ge 5) and (dxy lt 10): gridxy=2.0
  (dxy ge 1) and (dxy lt 5): gridxy=1.0
  (dxy ge 0.5) and (dxy lt 1.0): gridxy=0.2
  (dxy lt 0.5) : gridxy=0.1
endcase
res=gridxy

; character sizes

maxplot=multiplot(1)
case 1 of 
  (maxplot ge 5): lcharsize=0.5
  (maxplot ge 3) and (maxplot lt 5): lcharsize=0.6
  (maxplot ge 2) and (maxplot lt 3): lcharsize=0.9
  (maxplot ge 1) and (maxplot lt 2): lcharsize=1.2
endcase


; set projection and stretch

if(n_elements(projection) eq 0)then begin
  projection=8
endif
 
stretch=0



if(cross_section ne 0)then begin
  plotlegend=0
  margin=[0.08,0.03,0.08,0.03]
endif else begin
  margin=[0.0,0.0,0.0,0.0]
endelse


;-----------------------------------------------------------------------   
; open postscript file for printing

filetem='PlotPlume'+filetext
file=datadir+delim+filetem+'.ps'
set_plot,'ps'
device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
device,/color
device,filename=file
loadct,5  
stretch,32
tek_color

col1=!d.n_colors-30
col2=35

;-----------------------------------------------------------------------   
for it=0,numfiles-1 do begin

  pos=(it) mod (nxplot*nyplot)
  multiplot=[pos,nxplot,nyplot]
  
;----------------------------------------------------------------------  
; read plume

  readtraj,fieldfiles(it),modtitle,dates,xlon,ylat,temp,press,ptemp,bldepth,$
           cloud,rh,wspd,wdir,zqfe,zqnh,zp,zfl,age,namever=namever
	     
   if(plotpmsl ne 0)then begin
     readmet,datadir,timelist(it),selgrid,'Pmsl',lonp,latp,pmsl,titp,headp
   endif else begin
     pmsl=0
   endelse

   if(plotuv ne 0)then begin
     readmet,datadir,timelist(it),selgrid,'U',lonp,latp,u,titp,headp
     readmet,datadir,timelist(it),selgrid,'V',lonp,latp,v,titp,headp
   endif else begin
     u=0
     v=0
   endelse
   
   if(plottopog ne 0)then begin
     readmet,datadir,timelist(it),selgrid,'Topog',lont,latt,topog,titt,headt
   endif else begin
     topog=0
   endelse

   if(plotqfe ne 0)then begin
     zuse=zqfe
   endif else begin
     zuse=zqnh
   endelse
   
   if(plotfl ne 0)then zuse=zfl
  
;---------------------------------------------------------------------- 
; select subset of particles

  if(plumeres gt 1)then begin
    nplume=n_elements(xlon)
    plumenum=findgen(nplume)
    index=where(plumenum mod plumeres eq 0,nindex)
    xlon=xlon(index)
    ylat=ylat(index)
    zuse=zuse(index)
    age=age(index)
    blflag=bldepth(index)-zuse(index)
  endif

   index=reverse(sort(age))
   xlon=xlon(index)
   ylat=ylat(index)
   zuse=zuse(index)
   age=age(index)
   blflag=bldepth(index)-zuse(index)

;---------------------------------------------------------------------- 
; check plot area

  if(xmax gt 180)then begin
    index=where(xlon lt 0,nindex)
    if(nindex gt 0)then xlon(index)=xlon(index)+360.0
  endif


;---------------------------------------------------------------------- 
; select particles to plot and generate titles

; blflag=1 when the particle in the boundary layer, 0 otherwise
; zqfe : Height above Ground (zqnh = Height above sea level)

  if(boundarylayer ne 0) then begin
    index=where(blflag ge 0 and zuse gt zrange(0) and zuse lt zrange(1) and $
       xlon gt xmin and xlon lt xmax and ylat gt ymin and ylat lt ymax,nindex)    
    plottitle='Boundary layer particles between '+$
     strcompress(string(fix(zrange(0)),format='(i5)'),/remove_all)+' to '+$
     strcompress(string(fix(zrange(1)),format='(i5)'),/remove_all)+'m(agl)'
  endif else begin
    index=where(zuse gt zrange(0) and zuse lt zrange(1) and $
       xlon gt xmin and xlon lt xmax and ylat gt ymin and ylat lt ymax,nindex)
     if (plotfl ne 0) then begin 
      plottitle='Particles from FL '+$
       strcompress(string(fix(zrange(0)),format='(i5)'),/remove_all)+' to FL '+$
       strcompress(string(fix(zrange(1)),format='(i5)'),/remove_all)
     endif else begin  
      plottitle='Particles between '+$
       strcompress(string(fix(zrange(0)),format='(i5)'),/remove_all)+' to '+$
       strcompress(string(fix(zrange(1)),format='(i5)'),/remove_all)+'m(agl)'
     endelse
  endelse

  if(cross_section ne 0)then begin
     plottitle2='Longitude: '+$
     strcompress(string(xmin,format='(f7.1)'),/remove_all)+' to '+$
     strcompress(string(xmax,format='(f7.1)'),/remove_all)+'  Latitude: '+$
     strcompress(string(ymin,format='(f7.1)'),/remove_all)+' to '+$
     strcompress(string(ymax,format='(f7.1)'),/remove_all)
  endif else begin
     plottitle2=''
  endelse


;---------------------------------------------------------------------- 
; colours
  
   
  if(plumecolour eq 'age')then begin
    pcolor=col1+((age-agerange(0))/(agerange(1)-agerange(0)))*(col2-col1)
    contourlevels=agerange(0)+findgen(ncolors+1)/ncolors*(agerange(1)-agerange(0))
    contourcolors=col1+(contourlevels/agerange(1)-agerange(0))*(col2-col1)  
    legendtitle='Age (days)'
  endif
  if(plumecolour eq 'height')then begin
    pcolor=col1+(zuse/zrange(1))*(col2-col1)
    contourlevels=zrange(0)+findgen(7)/6*(zrange(1)-zrange(0))  
    contourcolors=col1+(contourlevels/zrange(1))*(col2-col1)   
    legendtitle='Height (m agl)'
    if (plotfl ne 0) then legentitle='Height (FL)'
  endif
  if(plumecolour eq 'bl')then begin
    pcolor=zuse*0.0+4.0
    blindex=where(blflag gt 0,blindexn)
    if(blindexn gt 0)then pcolor(blindex)=2
    plotlegend=0
  endif
  if (n_elements(usecontours) ne 0) then begin
    contourlevels=usecontours
    ncol=n_elements(contourlevels)-1
    dcol=(!d.n_colors-50)/(max([ncol-1,1]))
    contourcolors=(!d.n_colors-10)-(indgen(ncol)*dcol)
    contourcolors(0)=15
    contourcolors(1)=6
    contourcolors(2)=30
    contourcolors(3)=20
    For i=0,n_elements(contourlevels)-2 do begin
      colIndex=Where((age ge contourlevels(i)) and $
       (age lt contourlevels(i+1)),NData) 
      If (NData gt 0) Then Begin
        pcolor(colIndex)=ContourColors(i) 
      EndIf
    EndFor
  endif
  
;----------------------------------------------------------------------  
; plot position
  
  genposition,multiplot,plotposition,position,$
    plotlegend,plottitle,margin=margin

  
;----------------------------------------------------------------------
; If necessary determine projection name from number

projection_name,projection=projection,proj_name=proj_name

;-----------------------------------------------------------------------
; set projection origin

if (automap eq 1) then begin
  origin(0)=(ymin+ymax)/2.0
  origin(1)=(xmin+xmax)/2.0
endif else if (xmax-origin(1) gt 180) OR (xmin-origin(1) lt -180) then begin
  origin(1)=(xmin+xmax)/2.0
endif else if (automap eq -1) then begin
  origin[0] = 0
  origin(1)=(xmin+xmax)/2.0
endif 

; set up polar projection if requested (overrides other mapping specifications)
if (polar eq 1) then begin
  xmin=0
  xmax=360
  if (ymin lt 0) then ymin=0
  ymax=90
  origin(0)=90
  origin(1)=0
endif else if (polar eq -1) then begin
  xmin=0
  xmax=360
  ymin=-90
  if (ymax gt 0) then ymax=0
  origin(0)=-90
  origin(1)=0
endif

;-----------------------------------------------------------------------
; plot map area

  if(cross_section eq 0)then begin

    map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=plotposition,$
    color=0,/noerase,name=proj_name,/hires,/isotropic

    if(n_elements(topog) gt 1)then begin
      tcontourlevels=[-100.,0.001,200,500,1000,2000,5000,10000]
      tcontourcolors=[27,26,25,24,17,16,7]
      pixelplot,topog,lont,latt,tcontourlevels,tcontourcolors,$
      xmin,xmax,ymin,ymax
    endif

;-----------------------------------------------------------------------
; plot release point

    if(nindex gt 0)then plots,xlon(index),ylat(index),psym=3,$
     color=pcolor(index)

;-----------------------------------------------------------------------
; pressure

    if(n_elements(pmsl) gt 1)then begin
      contour,pmsl,reform(lonp(*,0)),reform(latp(0,*)),$
       levels=952+(findgen(25)*4.0),c_colors=intarr(25)+15
    endif

;-----------------------------------------------------------------------
; plot wind vectors
; map onto 50x50 grid

    if(n_elements(u) gt 1 and n_elements(v) gt 1)then begin
      u=congrid(u,50,50)
      v=congrid(v,50,50)
      lonp=congrid(lonp,50,50)
      latp=congrid(latp,50,50)
      velovect,u,v,reform(lonp(*,0)),reform(latp(0,*)),color=15,length=2
    endif

;-----------------------------------------------------------------------
; plot release point

    if (plotloc ne '')then begin
      plotloc,plotloc,zoom=[xmin,xmax,ymin,ymax],placenames=placenames
    endif 
    
;-----------------------------------------------------------------------
; plot continents etc.

    map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=plotposition,$
    color=0,/noerase,name=proj_name,/hires,/isotropic
    
    map_grid,color=18,linestyle=1,label=labelgrid

    if (countries eq 1) then map_continents,/countries,/hires,color=0
    if (coasts eq 1) then begin 
      if (res le 1) then begin
        map_continents,/coasts,/hires,color=0
      endif else begin
        map_continents,/coasts,color=0
      endelse
    endif
    if (states eq 1) then map_continents,/usa,/hires,color=3
    map_lakes,thick=2
   
;----------------------------------------------------------------------  
; cross section

  endif else begin
    
    if cross_section eq 1 then begin
; Decide on which slice to do based on X-Y ranges
      if(xmax-xmin gt ymax-ymin) then begin
; Longitude slice
        cross_section=2
      endif else begin
; Latitude slice
        cross_section=3
      endelse
    endif
    
; Longitude slice
    if cross_section eq 2 then begin  
    
      plot,xlon(index),zuse(index),/nodata,/xstyle,/ystyle,position=plotposition,/noerase,$
       xrange=[xmin,xmax],yrange=zrange,ytitle='Height',xtitle='Longitude',$
       charsize=lcharsize*0.75
        
      plots,xlon(index),zuse(index),color=pcolor(index),psym=3
      
; Latitude slice
    endif else if cross_section eq 3 then begin
      
      plot,ylat(index),zuse(index),/nodata,/xstyle,/ystyle,position=plotposition,$
       /noerase, $
       xrange=[ymin,ymax],yrange=zrange,ytitle='Height (m)',xtitle='Latitude',$
       charsize=lcharsize*0.75

      plots,ylat(index),zuse(index),color=pcolor(index),psym=3
         
    endif else begin
      print,'Incorrect cross_section flag - only 0-3 are valid'
      stop
    endelse

  endelse

;---------------------------------------------------------------------- 
; title
  
  day=strtrim(dates(0).day,2)
  month=strtrim(dates(0).month,2)
  year=strtrim(dates(0).year,2)
  hour=strtrim(dates(0).hour,2)
  minute=strtrim(dates(0).minute,2)
  
  dt_to_str,dates(0),date,time,date_fmt=2,time_fmt=-2
  valtime='valid at '+date+'  '+time+' UTC'

  xyouts,(plotposition(2)+plotposition(0))/2.0,$
        (plotposition(3)+(plotposition(3)-plotposition(1))/20.0),$
         valtime,/normal,charsize=lcharsize,Alignment=0.5    

      
;-----------------------------------------------------------------------
; legend

  if(plotlegend ne 0) then begin

    zoomleg=plotposition
  
; colours
  
    XRange=Zoomleg(2)-Zoomleg(0)
    XMinc=Zoomleg(0)+XRange*0.025
    Xinc=0.95*xrange/(N_Elements(Contourlevels)-1)

    YRange=Zoomleg(3)-Zoomleg(1)
    YMinp=Zoomleg(1)-YRange*0.15
    Yinc=Yrange*0.05

    For Colour=0,N_Elements(ContourLevels)-2 Do Begin
      Xminp=xminc+Colour*Xinc
      PolyFill,[XMinp,XMinp+xinc,XMinp+xinc,XMinp],$ 
       [YMinp,YMinp,YMinp+yinc,Yminp+yinc],$
       Color=ContourColors(Colour),/normal
    Endfor

; text

    For Colour=0,N_Elements(ContourLevels)-1 Do Begin
      Xminp=xminc+Colour*Xinc
      if (n_elements(usecontours) eq 0) then begin
        XYOuts,XMinp,YMinp-0.75*yinc,Alignment=0.5,$
         StrTrim(String(ContourLevels(Colour),Format='(e12.4)'),2),$
         charsize=lcharsize*0.75,color=axiscol,/normal
      endif else begin
        XYOuts,XMinp,YMinp-0.75*yinc,Alignment=0.5,$
         StrTrim(String(ContourLevels(Colour),Format='(f4.1)'),2),$
         charsize=lcharsize*0.75,color=axiscol,/normal
      endelse
    EndFor

    XYOuts,XMinc+0.48*Xrange,YMinp+1.3*yinc,Alignment=0.5,legendtitle,$
    charsize=lcharsize*0.75,color=axiscol,/normal

  Endif 

;---------------------------------------------------------------------- 
; annotation    
       
  if(plotheader ne 0 and pos eq 0)then begin

    xyouts,0.5,0.98,modtitle(0),$
     /normal,alignment=0.5,charsize=1.3,color=axiscol 

    xyouts,0.5,0.95,modtitle(1),$
     /normal,alignment=0.5,charsize=1.3,color=axiscol

    xyouts,0.5,0.92,plottitle,$
     /normal,alignment=0.5,charsize=1.3,color=axiscol
     
    if(plottitle2 ne '')then begin
      xyouts,0.5,0.895,plottitle2,$
       /normal,alignment=0.5,charsize=1.3,color=axiscol
    endif

    npart=n_elements(xlon)
    partstr='Number of particles: '+strtrim(npart,2)    
    xyouts,0.05,0.12,partstr,/normal,charsize=0.9,color=axiscol 
    relloc='Release location:      '+strtrim(xlon1(0),2)+'  '+strtrim(ylat1(0),2)
    xyouts,0.05,0.10,relloc,/normal,charsize=0.9,color=axiscol
    dt_to_str,date1(0),date,time,date_fmt=2,time_fmt=-2
    reltime='Release time:           '+date+'  '+time
    xyouts,0.05,0.08,reltime,/normal,charsize=0.9,color=axiscol

    xyouts,0.6,0.12,modtitle(3),/normal,charsize=0.9,color=axiscol
    xyouts,0.6,0.10,modtitle(2),/normal,charsize=0.9,color=axiscol
    xyouts,0.5,0.00,'Met Office Crown copyright',$
       /normal,alignment=0.5,charsize=0.8,color=axiscol 

    image=read_image('MO_Master_B.jpg')
    loadct,0
    tv,image,0.82,0.91,true=1,xsize=0.15,/normal
    loadct,5
    tek_color

  endif 
      
;-----------------------------------------------------------------------
; printer

  if(pos eq ((nxplot*nyplot)-1))then begin
    erase
  endif
     

endfor
  
device,/close_file


;----------------------------------------------------------------------- 
; generate gifs

if((gengif eq 1 or genjpg eq 1 or genanim eq 1) and pc eq 0)then begin
  genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
    genanim=genanim,xpixel=xpixel
endif


end
