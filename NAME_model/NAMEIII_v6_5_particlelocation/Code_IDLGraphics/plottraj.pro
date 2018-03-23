pro plottraj,datadir,sources=sources,clabel=clabel,zoom=zoom,$
   labint=labint,filetext=filetext,gengif=gengif,fit=fit,$
   projection=projection,cycle=cycle,genjpg=genjpg,namever=namever,$
   origin=origin,countries=countries,coasts=coasts,states=states,$
   automap=automap,polar=polar
   
;-----------------------------------------------------------------------

; procedure to display NAMEIII trajectory model output
;
; DBR may 1996
;
; Changes
; 29/05/2002 DBR  Remove labelling, revised projection and stretch options
; 10/07/2002 DBR  plot only if in map area, add projection argument,
;                 Number of colors to cycle through
; 31/01/2006 ARJ  Updated to handle NAME III particle/puff output files
; 01/07/2009 CSW  Updated to use new readtraj.pro

;-----------------------------------------------------------------------
; argument list

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
;  states          : set to 0 to not plot individual states on maps of
;                    the USA
;  clabel          : set to 0 to not plot time labels on trajectories
;  labint          : interval (in hours) between labels on  trajectories
;  fit             : set to 1 to yrange on height versus time plot to
;                    highest height reached by trajectories
;  zoom            : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  gengif          : integer (0 or 1) =1 generates gif images (on unix)
;  genjpg          : integer (0 or 1) =1 generates jpg images (on unix)
;  namever         : integer indicating version of model e.g. 3 - NAME III(PPM)
;  filetext        : (string)  text added to filenames
;  cycle           : number of different colours to be used for
;                    trajectories (default is 8)

;-----------------------------------------------------------------------
; check arguments

if(n_elements(namever) eq 0) then namever=3
if(n_elements(cycle) eq 0) then cycle=8
if(n_elements(fit) eq 0)then fit=0
if(n_elements(gengif) eq 0)then gengif=1
if(n_elements(genjpg) eq 0)then genjpg=0
if(n_elements(labint) eq 0)then labint=12
if(n_elements(clabel) eq 0) then clabel=1
if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext='_'+filetext
endelse
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
file=datadir+delim+'Trajectory_plot'+filetext+'.ps'
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

maxzqnh = 0
if(n_elements(zoom) eq 0)then begin 

  for i=0,numfiles-1 do begin

    filename=fieldfiles(i)
    readtraj,filename,modtitle,date1,xt,yt,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
             zqfet,zqnht,zpt,zflt,namever=namever

    if (max(zqnht) ge maxzqnh) then begin
        maxzqnh = max(zqnht)
     endif
   
    if(i eq 0)then begin
        x=xt
        y=yt
    endif else begin
        x=[x,xt]
        y=[y,yt]
    endelse
   
  endfor

 ; If only one trajectory loaded dulplicate its
 ; latitude/longitude arrays for maplimits.pro
  xdims = size(x, /dimensions)
  if (xdims(0) eq 1) then begin
     x = [x,x]
     y = [y,y]
  endif

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
;If necessary convert projection number to a projection name

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

map_grid,color=Col18
if (polar ne 0) then begin
  map_grid,color=Col15,/label,/nogrid,charsize=0.9*lcharsize
endif else begin
  latlab=xmin+(xmax-xmin)/20.0
  lonlab=ymin+(ymax-ymin)/20.0
  map_grid,/nogrid,/label,color=15,latlab=latlab,lonlab=lonlab,charsize=0.9*lcharsize
endelse

;-----------------------------------------------------------------------
; loop through all trajectories

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

  if(xmax gt 180)then begin
    index=where(x1 lt 0,nindex)
    if(nindex gt 0)then  x1(index)=x1(index)+360.0
  endif
  
  
; plot trajectory   
   
  for n=0,n_elements(x1)-2 do begin
  
   if(abs(x1(n+1)-x1(n)) lt 10 and $
      x1(n) gt xmin and x1(n+1) gt xmin and $
      x1(n) lt xmax and x1(n+1) lt xmax and $
      y1(n) gt ymin and y1(n+1) gt ymin and $
      y1(n) lt ymax and y1(n+1) lt ymax)then begin 

      plots,[x1(n),x1(n+1)],[y1(n),y1(n+1)],$
            color=col,thick=2,linestyle=line      

   endif
  endfor


; symbol at start and end

  plots,[x1(0)],[y1(0)],psym=2,color=col

  if(x1(n_elements(x1)-1) gt xmin and x1(n_elements(x1)-1) lt xmax and $
     y1(n_elements(x1)-1) gt ymin and y1(n_elements(x1)-1) lt ymax)then begin
      
       plots,[x1(n_elements(x1)-1)],[y1(n_elements(x1)-1)],psym=4,$
        color=col,symsize=0.6
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

    plots,[x1(label(n))],[y1(label(n))],psym=4,color=col,symsize=0.6

    if(clabel ne 0 and numfiles le 12)then begin

      xyouts,x1(label(n)),y1(label(n)), strmid(time(label(n)),0,2)+'Z'+$
        strmid(date(label(n)),0,2),charsize=lcharsize*0.7,color=axiscol

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
; plot vertical trajectories

if(fit ne 0)then begin
    if ( maxzqnh ne 0 ) then begin
        yrange=[0,maxzqnh*1.1]
    endif else begin
        yrange=[0,max(zqnh)*1.1]
    endelse
endif else begin
  yrange=[0.0,16000]
endelse

ylab='Height (m amsl)' 
position=[0.08,0.25,0.92,0.45]
  

for s=0,numfiles-1 do begin

  filename=fieldfiles(s)

  readtraj,filename,modtitle,dates,x1,y1,temp,press,ptemp,bldepth,cloud,rh,wspd,wdir,$
           zqfe,zqnh,zp,zfl,namever=namever

  line=0
  col = 2*(s mod cycle)+2

  if (n_elements(zqfe) gt n_elements(zqnh)) then begin
    zqnh=zqfe
    ylab='Height (m agl)'
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
; add text

xyouts,0.5,0.97,'NAME',/normal,charsize=lcharsize*1.6,alignment=0.5,color=axiscol 
xyouts,0.5,0.94,modtitle(5),/normal,charsize=lcharsize*1.6,alignment=0.5,color=axiscol

xyouts,0.1,0.16,'Simulation Description',charsize=lcharsize*1.2,/normal,color=axiscol

xyouts,0.1,0.12,'Number of trajectories:',/normal,charsize=lcharsize,color=axiscol
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
; generate gifs

filetem='none'
if((gengif eq 1 or genjpg eq 1) and pc eq 0)then begin
  genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
    genanim=0
endif
  
End

