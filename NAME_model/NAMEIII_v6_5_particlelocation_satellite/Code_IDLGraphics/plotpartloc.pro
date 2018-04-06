pro plotpartloc,datadir,sources=sources,clabel=clabel,zoom=zoom,$
   labint=labint,filetext=filetext,gengif=gengif,fit=fit,numtraj=numtraj,$
   projection=projection,cycle=cycle,genjpg=genjpg,ht_contours=ht_contours,$
   namever=namever,xloc=xloc,yloc=yloc,PlotByTime=PlotByTime,$
   plotmet=plotmet,origin=origin,countries=countries,coasts=coasts,$
   states=states,automap=automap

;-----------------------------------------------------------------------
;
; procedure to display NAMEIII particle location and trajectory output
;
; CSW 03/2005
;
; Based on plottraj.pro
;
; Changes:
; 22/06/2007  CSW Updated to handle NAME III particle output files
; 2/10/2007   CSW Added plot locations capability
; 5/10/2007   CSW Added plot by time capability
; 30/5/2008   LB Added plot met variables capability for NAMEIII
; 11/6/2009   CSW modified plot met to only work when requested

;-----------------------------------------------------------------------
; Input variables:
; numtraj       - the number of trajectories to be plotted
; ht_contours   - if set to 1 then the trajectories are plotted with colours
;                that correspond to their final height (useful for back trajectories)
; xloc and yloc - specify lon lat coords as arrays for points to be plotted 
;                 on the map
; PlotByTime    - plots the trajectories by colour in 3 hour blocks
; plotmet       - if set to 1 then met variables along trajectories are plotted
;                 (not set up for NAME)
;  projection      : integer (1 - 16) - projection for mapping (for options see
;                    IDL map option) #8 (cylindrical) is standard
;  polar           : integer, 0 for non-polar projection, -1 for south pole, 
;                    1 for north pole, projection modes 1 (stereographic)
;                    and 2 (orthographic) are recommended, but it will work with other
;                    projection modes.  The polar keyword overrides the origin keyword.
;  origin       : 2-element array [lat,long] to be used for origin of map projection.  
;                 If this keyword is not set, the automap=-1 option will be enabled
;                 (unless 'polar' is also set).
;  automap      : (default 0) set to 1 to get map projection with origin in centre of grid
;                 or -1 for map centred by longitude only.  If an orgin has been set, 
;                 automap will override it, unless 'polar' is set to -1 or +1.
;  countries       : set to 0 to not plot national borders
;  coasts          : set to 1 to plot lakes and extra islands.  N.B. do not use this
;                    with the 'hires' keyword on continent-sized maps, and do not use
;                    it without 'hires' set on anything smaller (however it will plot
;                    every lake, however small).  You will need to use it without 'hires' 
;                    for getting the Great Lakes on maps of the whole of N. America, and 
;                    the Caspian Sea on large-scale maps of Asia, and with 'hires' to
;                    show the smaller islands around the UK.
;  states          : set to 0 to not plot individual states on maps of the USA
;-----------------------------------------------------------------------


; check arguments

if(n_elements(namever) eq 0) then namever=3
if(n_elements(cycle) eq 0) then cycle=8
if(n_elements(fit) eq 0)then fit=0
if(n_elements(numtraj) eq 0) then numtraj=10
if(n_elements(gengif) eq 0)then gengif=1
if(n_elements(genjpg) eq 0)then genjpg=0
if(n_elements(labint) eq 0)then labint=12
if(n_elements(clabel) eq 0) then clabel=1
if(n_elements(ht_contours) eq 0) then ht_contours=0
if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext='_'+filetext
endelse
if((n_elements(xloc) eq 0) or (n_elements(yloc) eq 0)) then begin
   plot_loc=0
endif else begin
  plot_loc=1
endelse
if(n_elements(PlotByTime) eq 0) then PlotByTime=0
if(n_elements(plotmet) eq 0) then plotmet=0
labelgrid=1
if(n_elements(projection) eq 0)then begin
  projection=8
  margin=[0.0,0.0,0.0,0.0]
endif else begin
  margin=[0.03,0.03,0.0,0.0]
endelse
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

;-----------------------------------------------------------------------
; printer setup

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse 

lcharsize=0.8
!p.charsize=lcharsize
set_plot,'ps'
tek_color
colors=[3,18,5,10]   
device,/portrait,ysize=25,yoffset=2
device,/times,/color
file=datadir+delim+'PartLoc_plot'+filetext+'.ps'
device,filename=file
!p.font=0
backcol=!d.n_colors
axiscol=0


;-----------------------------------------------------------------------
; find all files and determine number if times and species

fieldfiles=findfile(datadir+delim+'Data_*.txt')
numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)
  endfor
endif

;-----------------------------------------------------------------------
; determine map limits    


if(n_elements(zoom) eq 0)then begin 

  for i=0,numtraj-1 do begin

    filename=fieldfiles(i)
    readtraj,filename,modtitle,date1,xt,yt,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
             zqfet,zqnht,zpt,zflt,namever=namever
   
    if(i eq 0)then begin
        x=xt
        y=yt
    endif else begin
        x=[x,xt]
        y=[y,yt]
    endelse
   
  endfor

;-----------------------------------------------------------------------
; map limits

  tot=x*0.0+1.0
  zoom=fltarr(4)
  maplimits,x,y,tot,xmin,xmax,ymin,ymax
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

; grid lines
  
dxy=ymax-ymin
case 1 of 
  (dxy ge 90): gridxy=20.0
  (dxy ge 40) and (dxy lt 90): gridxy=10.0
  (dxy ge 20) and (dxy lt 40): gridxy=5.0
  (dxy ge 10) and (dxy lt 20): gridxy=2.0
  (dxy ge 5) and (dxy lt 10): gridxy=1.0
  (dxy ge 2) and (dxy lt 5): gridxy=0.5
  (dxy ge 1)   and (dxy lt 2): gridxy=0.2
  (dxy ge 0.5) and (dxy lt 1): gridxy=0.1
  (dxy ge 0)   and (dxy lt 0.5): gridxy=0.05
endcase

res=fix(max([gridxy,1])*2)

;-----------------------------------------------------------------------
; map position

position=[0.0,0.50,1.0,0.9]
 
;-----------------------------------------------------------------------
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

map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=position,$
color=0,/noerase,name=proj_name,/hires,/isotropic

if (countries eq 1) then map_continents,/countries,/hires
if (coasts eq 1) then begin 
  if (res le 1) then begin
    map_continents,/coasts,/hires
  endif else begin
    map_continents,/coasts
  endelse
endif
if (states eq 1) then map_continents,/usa,/hires
map_lakes,thick=2

;-----------------------------------------------------------------------
; plot lon/lat

map_grid,color=18
map_continents,/hires

if (polar ne 0) then map_grid,color=15,/label,/nogrid

for iy=-90.,90.,gridxy do begin

  if(labelgrid gt 0 and (~strmatch(strlowcase(proj_name),'stereographic')$
      or ~strmatch(strlowcase(proj_name),'orthographic')) and $
     iy gt ymin+0.05*(ymax-ymin) and $
     iy lt ymax-0.05*(ymax-ymin) )then begin
       
     if(gridxy lt 1)then begin
       txt=strcompress(string(abs(iy),format='(f5.2)'),/remove_all)
     endif else begin
       txt=strcompress(string(fix(abs(iy)),format='(i2)'),/remove_all)
     endelse

     xyouts,/data, xmin+0.1*dxy,iy-(0.05*gridxy),txt,$
       charsize=0.6*lcharsize,alignment=0.5,color=15
  endif

endfor

for ix=-180.,360.,gridxy do begin

  if(labelgrid gt 0 and (~strmatch(strlowcase(proj_name),'stereographic')$
      or ~strmatch(strlowcase(proj_name),'orthographic')) and $
     ix gt xmin+0.5*gridxy and $
     ix lt xmax-0.5*gridxy )then begin

     if(gridxy lt 1)then begin
       txt=strcompress(string(abs(ix),format='(f6.2)'),/remove_all)
     endif else begin
       txt=strcompress(string(fix(abs(ix)),format='(i3)'),/remove_all)
     endelse

     xyouts, /Data, ix,ymin+0.1*dxy,txt,$
      charsize=0.6*lcharsize,alignment=0.5,color=15
  endif

endfor

;-----------------------------------------------------------------------
; loop through all particle location trajectories
maxht=0
ht_range=500.
;ht_col=indgen(9)
ht_col=[0,2,3,4,5,6,7,8,9]

numfiles=numtraj

for s=0,numfiles-1 do begin

; read data
  filename=fieldfiles(s)

  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

; date and time strings

  dt_to_str,dates,date,time,date_fmt=2,time_fmt=-2

; colors

  line=0
  col = 2*(s mod cycle)+2
 
; check plot area

  if(xmax gt 180.)then begin
    index=where(x1 lt 0,nindex)
    if(nindex gt 0)then  x1(index)=x1(index)+360.0
  endif

; determine maximum height for fitted profile

  if (n_elements(zqfe) gt n_elements(zqnh)) then begin
    zqnh=zqfe
  endif

  maxa=max(zqnh)
  if (maxa gt maxht) then maxht=maxa
  
  if (ht_contours eq 1) then begin
    ht_level = NINT([maxa/ht_range])
    col = ht_col(ht_level)
    print,'ht_level',ht_level
  endif

; plot particle trajectory   
   
  for n=0,n_elements(x1)-2 do begin
  
   if(abs(x1(n+1)-x1(n)) lt 10 and $
      x1(n) gt xmin and x1(n+1) gt xmin and $
      x1(n) lt xmax and x1(n+1) lt xmax and $
      y1(n) gt ymin and y1(n+1) gt ymin and $
      y1(n) lt ymax and y1(n+1) lt ymax)then begin 

      if (PlotByTime eq 1) then begin
        firsttime=time(0)
        if ((fix(time(n)) ge 0)    and (fix(time(n)) lt 300))  then col=2
        if ((fix(time(n)) ge 300)  and (fix(time(n)) lt 600))  then col=3
        if ((fix(time(n)) ge 600)  and (fix(time(n)) lt 900))  then col=4
        if ((fix(time(n)) ge 900)  and (fix(time(n)) lt 1200)) then col=5
        if ((fix(time(n)) ge 1200) and (fix(time(n)) lt 1500)) then col=6
        if ((fix(time(n)) ge 1500) and (fix(time(n)) lt 1800)) then col=9
        if ((fix(time(n)) ge 1800) and (fix(time(n)) lt 2100)) then col=8
        if ((fix(time(n)) ge 2100) and (fix(time(n)) le 2359)) then col=10
      endif

      oplot,[x1(n),x1(n+1)],[y1(n),y1(n+1)],$
            color=col,thick=1,linestyle=line      

   endif
  endfor


; symbol at start and end

  plots,[x1(0)],[y1(0)],psym=2,color=0

  if(x1(n_elements(x1)-1) gt xmin and x1(n_elements(x1)-1) lt xmax and $
     y1(n_elements(x1)-1) gt ymin and y1(n_elements(x1)-1) lt ymax)then begin
      
       plots,[x1(n_elements(x1)-1)],[y1(n_elements(x1)-1)],psym=4,$
        color=col,symsize=0.6
   endif 

; Plot additional location(s) on map
   if plot_loc eq 1 then begin
      nloc=n_elements(xloc)      
      for l=0,nloc-1 do begin
         plots,xloc(l),yloc(l),psym=2,color=0
      endfor 
   endif

; label with times

  if(clabel ne 0 and numfiles le 12)then begin
    xyouts,x1(0),y1(0),time(0)+'Z '+date(0),$
      charsize=lcharsize*0.7,color=axiscol
  endif

  label=where(dates.hour mod labint eq 0 and dates.minute eq 0,npoints)
  label=[label,n_elements(dates)-1]
  
  for n=0,npoints-1 do begin

   if(x1(label(n)) gt xmin and x1(label(n)) lt xmax and $
      y1(label(n)) gt ymin and y1(label(n)) lt ymax)then begin

    if(clabel ne 0 and numfiles le 12)then begin

      plots,[x1(label(n))],[y1(label(n))],psym=4,color=col,symsize=0.6

    endif
   endif
  endfor

endfor

;-----------------------------------------------------------------------     
; set x-axis (time) 
     
juldat=dates.julian
range=abs(double(max(juldat)-min(juldat)))
timeaxis,range
        
;-----------------------------------------------------------------------
; plot vertical particle trajectories

if(fit ne 0)then begin
  yrange=[0,maxht*1.1]
endif else begin
  yrange=[0.0,10000]
endelse

ylab='Height (m amsl)' 
position=[0.09,0.25,0.92,0.45]
  

for s=0,numfiles-1 do begin

  filename=fieldfiles(s)

  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  if (n_elements(zqfe) gt n_elements(zqnh)) then begin
    zqnh=zqfe
    ylab='Height (m agl)'
  endif

  line=0
  col = 2*(s mod cycle)+2

  if (ht_contours eq 1) then begin
    maxa=max(zqnh)
    ht_level = NINT([maxa/ht_range])
    col = ht_col(ht_level)
  endif

  if(s eq 0)then begin
    plot,juldat,zqnh,/nodata,$
     yrange=yrange,ytitle=ylab,$
     /xstyle,/ystyle,position=position,/noerase,$
     color=axiscol,background=backcol
  endif
  
  oplot,juldat,zqnh,color=col,thick=2,linestyle=line

endfor

;-----------------------------------------------------------------------
; add legend(s)

if (ht_contours eq 1) then begin

  label=['0-500','500-1000','1000-1500','1500-2000','2000-2500',$
         '2500-3000','3500-4000','4000-4500','4500-5000']

  
  XMinc=0.07
  Xinc=0.80/(N_Elements(ht_col)-1)

  YMinp=0.18
  
  XYOuts,0.07,0.193,Alignment=0.1,'Max height:',$
      charsize=lcharsize*0.8,color=0,/normal

  For Colour=0,N_Elements(ht_col)-1 Do Begin
     Xminp=xminc+Colour*Xinc
     XYOuts,XMinp,YMinp,Alignment=0.1,$
      StrTrim(String(label(Colour),Format='(A10)'),2),$
      charsize=lcharsize*0.75,color=ht_col(Colour),/normal
  EndFor
Endif

if (PlotByTime eq 1) then begin

  times=['0000','0300','0600','0900','1200','1500','1800','2100','0000','0300','0600','0900','1200','1500','1800','2100']
  Colours=[2,3,4,5,6,9,8,10,2,3,4,5,6,9,8,10]

  for t=0,7 do begin
    if (firsttime eq times(t)) then begin
      tc=t
    endif
  endfor    

  label=[times(tc)+'-'+times(tc+1),times(tc+1)+'-'+times(tc+2),times(tc+2)+'-'+times(tc+3),times(tc+3)+'-'+times(tc+4),$
         times(tc+4)+'-'+times(tc+5),times(tc+5)+'-'+times(tc+6),times(tc+6)+'-'+times(tc+7),times(tc+7)+'-'+times(tc+8)] 

  Colourarr=[Colours(tc),Colours(tc+1),Colours(tc+2),Colours(tc+3),Colours(tc+4),Colours(tc+5),Colours(tc+6),Colours(tc+7)]

  XMinc=0.07
  Xinc=0.80/(N_Elements(label)-1)

  YMinp=0.47

  For c=0,N_Elements(label)-1 Do Begin
     Xminp=xminc+c*Xinc
     XYOuts,XMinp,YMinp,Alignment=0.1,$
      StrTrim(String(label(c),Format='(A10)'),2),$
      charsize=lcharsize*0.75,color=Colourarr(c),/normal  
  EndFor 

Endif
;-----------------------------------------------------------------------
; add text

xyouts,0.5,0.97,'NAME',/normal,charsize=lcharsize*1.6,alignment=0.5,color=axiscol 
xyouts,0.5,0.94,'Particle Locations',/normal,charsize=lcharsize*1.6,alignment=0.5,color=axiscol

if(clabel ne 0 and numfiles le 12)then begin
  xyouts,0.5,0.91,'Markers every '+string(labint)+' hours',/normal,charsize=lcharsize*1.2, $
           alignment=0.5,color=axiscol
endif

xyouts,0.1,0.16,'Simulation Description',charsize=lcharsize*1.2,/normal,color=axiscol

xyouts,0.1,0.12,'Number of particles:',/normal,charsize=lcharsize,color=axiscol
xyouts,0.3,0.12,numfiles,/normal,charsize=lcharsize,color=axiscol

xyouts,0.1,0.1,'Release location:',/normal,charsize=lcharsize,color=axiscol
xyouts,0.3,0.1,x1(0),/normal,charsize=lcharsize,color=axiscol
xyouts,0.4,0.1,y1(0),/normal,charsize=lcharsize,color=axiscol
xyouts,0.1,0.08,'Release time:',/normal,charsize=lcharsize,color=axiscol
xyouts,0.3,0.08,date(0),/normal,charsize=lcharsize,color=axiscol
xyouts,0.4,0.08,time(0),/normal,charsize=lcharsize,color=axiscol

xyouts,0.1,0.06,'Name version:',/normal,charsize=lcharsize,color=axiscol
xyouts,0.3,0.06,modtitle(0),/normal,charsize=lcharsize,color=axiscol

xyouts,0.1,0.04,modtitle(1),charsize=lcharsize,/normal,color=axiscol


xyouts,0.5,0.01,'Met Office Crown copyright',$
    /normal,alignment=0.5,charsize=lcharsize*.8,color=axiscol
   
image=read_image('MO_Master_B.jpg')
loadct,0
tv,image,0.82,0.91,true=1,xsize=0.15,/normal

; printer
device,/close_file

;-----------------------------------------------------------------------

; plotting met variables (first page) 

if (plotmet eq 1) then begin

  ; making new post script doc for plotting met variables
  lcharsize=0.8
  !p.charsize=lcharsize
  set_plot,'ps'
  tek_color
  colors=[3,18,5,10]   
  device,/portrait,ysize=25,yoffset=2
  device,/times,/color
  file=datadir+delim+'PartLoc_plot_met1'+filetext+'.ps'
  device,filename=file
  !p.font=0
  backcol=!d.n_colors
  axiscol=0

 
  ; plot temperature along particle trajectories
    
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
  
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
   
    temp=temp-273.15
    yrange=[-100,50]
    ylab='Temperature (degC)'

  if (dates(1).julian lt dates(0).julian) then begin
    dt_range=[dates(n_elements(dates)-1).julian,dates(0).julian]
  endif else begin
    dt_range=[dates(0).julian,dates(n_elements(dates)-1).julian]
  endelse       
    
       if (s eq 0) then begin
       
  	 plot,juldat,temp,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.75,0.92,0.9],/noerase,$
  	   color=axiscol,background=backcol
       endif 
    oplot,juldat,temp,color=col,thick=2,linestyle=line
       
 endfor
       

; plot relative humidity along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
  
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
    
    yrange=[0,100]
    ylab='Relative Humidity (%)'
   
   
  	if (s eq 0) then begin

  	  plot,juldat,rh,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.55,0.92,0.7],/noerase,$
  	   color=axiscol,background=backcol
       endif 
    oplot,juldat,rh,color=col,thick=2,linestyle=line
   
 endfor
     
     
; plot wind speed along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
   
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
   
    yrange=[0,50]
    ylab='Wind Speed (m/s)' 
   
  
       if(s eq 0)then begin
  	  plot,juldat,wspd,/nodata,$
  	  yrange=yrange,ytitle=ylab,$
  	  /xstyle,/ystyle,position=[0.09,0.35,0.92,0.5],/noerase,$
  	  color=axiscol,background=backcol
       endif
     oplot,juldat,wspd,color=col,thick=2,linestyle=line
  
 endfor
   
       
; plot wind direction along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
   
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever
       
  line=0
  col = 2*(s mod cycle)+2
   

    yrange=[0,360]
    ylab='Wind Direction (deg)'
     
   
       if(s eq 0)then begin
  	   plot,juldat,wdir,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.15,0.92,0.3],/noerase,$
  	   color=axiscol,background=backcol
       endif
      oplot,juldat,wdir,color=col,thick=2,linestyle=line
 endfor


  xyouts,0.5,0.00,'Met Office Crown copyright',$
        /normal,alignment=0.5,charsize=lcharsize*.8,color=axiscol
   
  image=read_image('MO_Master_B.jpg')
  loadct,0
  tv,image,0.82,0.91,true=1,xsize=0.15,/normal


  device,/close_file
;-----------------------------------------------------------------------

; plotting met variables (second page) 


  ; making new post script doc for plotting met variables
  lcharsize=0.8
  !p.charsize=lcharsize
  set_plot,'ps'
  tek_color
  colors=[3,18,5,10]   
  device,/portrait,ysize=25,yoffset=2
  device,/times,/color
  file=datadir+delim+'PartLoc_plot_met2'+filetext+'.ps'
  device,filename=file
  !p.font=0
  backcol=!d.n_colors
  axiscol=0

 
; plot pressure along particle trajectories
    
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
   
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
   
     press=press/100
     yrange=[0,1100]
     ylab='Pressure (hPa)' 
       
    
       if (s eq 0) then begin
  	   plot,juldat,press,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.75,0.92,0.9],/noerase,$
  	   color=axiscol,background=backcol
       endif 
    oplot,juldat,press,color=col,thick=2,linestyle=line
       
 endfor
       

; plot cloud fraction along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
  
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
    
    yrange=[0,8]
    ylab='Cloud Fraction (Oktas)'
   
   
  	if (s eq 0) then begin
  	   plot,juldat,cloud,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.55,0.92,0.7],/noerase,$
  	   color=axiscol,background=backcol
       endif 
    oplot,juldat,cloud,color=col,thick=2,linestyle=line
   
 endfor
     
     
; plot potential temperature along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
    
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
   
    
     yrange=[200,400]
     ylab='Potential Temperature (K)' 

  
       if(s eq 0)then begin
  	  plot,juldat,ptemp,/nodata,$
  	  yrange=yrange,ytitle=ylab,$
  	  /xstyle,/ystyle,position=[0.09,0.35,0.92,0.5],/noerase,$
  	  color=axiscol,background=backcol
       endif
     oplot,juldat,ptemp,color=col,thick=2,linestyle=line
  
 endfor
   
       
; plot boundary layer depth along particle trajectories
 for s=0,numfiles-1 do begin

  filename=fieldfiles(s)
  position=[0.09,0.25,0.92,0.45]
  
  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2
   

    yrange=[0,4000]
    ylab='Boundary Layer Depth (m)'
     
   
       if(s eq 0)then begin
  	   plot,juldat,bldepth,/nodata,$
  	   yrange=yrange,ytitle=ylab,$
  	   /xstyle,/ystyle,position=[0.09,0.15,0.92,0.3],/noerase,$
  	   color=axiscol,background=backcol
       endif
      oplot,juldat,bldepth,color=col,thick=2,linestyle=line
 endfor


  xyouts,0.5,0.00,'Met Office Crown copyright',$
        /normal,alignment=0.5,charsize=lcharsize*.8,color=axiscol
   
  image=read_image('MO_Master_B.jpg')
  loadct,0
  tv,image,0.82,0.91,true=1,xsize=0.15,/normal

  device,/close_file

endif

;----------------------------------------------------------------------- 
; generate gifs

filetem='none'
if((gengif eq 1 or genjpg eq 1) and pc eq 0)then begin
  genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
    genanim=0
endif
  
End

