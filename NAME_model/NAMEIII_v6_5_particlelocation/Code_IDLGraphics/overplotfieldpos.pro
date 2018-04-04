pro overplotfieldpos,lont,latt,plotdata, $
 zoom=zoom,contourlevels=contourlevels,contourcolors=contourcolors,$
 fieldunits=fieldunits,plotlegend=plotlegend,$
 position=position,xyrel=xyrel,projection=projection,$
 multiplot=multiplot,modtitle=modtitle,labelgrid=labelgrid,$
 pmsl=pmsl,u=u,v=v,topog=topog,exact=exact,highres=highres,$
 plotloc=plotloc,type=type,tlon=tlon,tlat=tlat,$
 plotlincomtopog=plotlincomtopog,polar=polar,origin=origin,automap=automap,$
 countries=countries,coasts=coasts,states=states,placenames=placenames,$
 plotvaacobs=plotvaacobs,imagename=imagename

;-----------------------------------------------------------------------
; procedure to draw contour plot overlaid by map
; DBR  March 2002
;
; 29/05/2002  DBR Apply projection=1, with stretch=1 only for large areas,
;  modify map colours
;
; 11/07/2002  DBR Add filled contour option - through
;  type option ('pixel' ot 'contour') 
;
;-----------------------------------------------------------------------
; Arguments
; required:
;  lont            : (float) array of longitudes
;  latt            : (float) array of latitudes
;  plotdata        : (float) array of values to plot
;
; optional:
;  zoom            : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  contourlevels   : 
;  contourcolors   :
;  fieldunits      :
;  plotlegend      : integer (0 or 1) =1 plots legend below each plot
;  position        : 
;  xyrel           : source position
;  projection      : integer (1 - 16) - projection for mapping (for options see
;                    IDL map option) #8 (cylindrical) is standard
;  multiplot       : intarr(3) number of plots, 
;                    (eg multiplot=[0,2,2] for 2 by 2 plots)
;                    first element can be any value
;  modtitle        : 
;  labelgrid       :
;  pmsl            : pressure contours
;  u               : horizontal windspeed
;  v               : horizontal windspeed
;  topog           : topography
;  exact           : integer (0 or 1) =1 map area based on lon/lat limits
;  highres         : integer (0 or 1) =1 plots highest resolution map possible
;  plotloc         : string filename conatining locations to overplot,
;                    assumes same format as tslocs.txt file
;  type            : string ('pixel' or 'contour') choose pixels, 
;                    filled contours or unfilled contours (default)
;  tlon            :
;  tlat            : 
;  plotlincomtopog :.
;  polar           : integer, 0 for non-polar projection, -1 for south pole, 
;                    1 for north pole, projection modes 1 (stereographic)
;                    and 2 (orthographic) are recommended, but it will work with other
;                    projection modes.  The polar keyword overrides the origin keyword.
;  origin          : 2-element array [lat,long] to be used for origin of map projection.  
;                    If this keyword is not set, the automap option will be enabled and the
;                    map origin will be the centre of the map itself, (unless the polar 
;                    keyword is set).
;  automap         : integer (0, 1 or -1) =1 map projection origin determined automatically
;                  : =-1 only longitude determined automatically and latitude set to 0
;  countries       : set to 0 to not plot national borders
;  coasts          : set to 1 to plot lakes and extra islands.  N.B. do not use this
;                    with the 'hires' keyword on continent-sized maps, and do not use
;                    it without 'hires' set on anything smaller (however it will plot
;                    every lake, however small).  You will need to use it without 'hires' 
;                    for getting the Great Lakes on maps of the whole of N. America, and 
;                    the Caspian Sea on large-scale maps of Asia, and with 'hires' to
;                    show the smaller islands around the UK.
;  states          : set to 0 to not plot individual states on maps of the USA
;  placenames      : set to 0 to not plot names of locations in  plotloc file 
;  plotvaacobs         : string filename conatining observations to overplot


;-----------------------------------------------------------------------

if(n_elements(plotlincomtopog) eq 0)then plotlincomtopog=0

;-----------------------------------------------------------------------
; lon/lat

lon=lont
lat=latt
dlon=lon(1,0)-lon(0,0)
dlat=lat(0,1)-lat(0,0)


; determine array size

s=size(plotdata)
itx=s(1)
ity=s(2)

index=where((plotdata gt 0 and lon gt 170) or $
            (plotdata gt 0 and lon lt -170),nindex)
	    
if(nindex gt 0)then begin
  type='pixel'
endif

;-----------------------------------------------------------------------
; check optional arguments

if(n_elements(type) eq 0)then type='nofill'
if(n_elements(exact) eq 0)then exact=0
if(n_elements(zoom) eq 0)then begin
  zoom=fltarr(4)
  maplimits,lont,latt,plotdata,xmin,xmax,ymin,ymax,exact=exact
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

if(n_elements(plotloc) eq 0)then plotloc=''
if(n_elements(plotvaacobs) eq 0)then plotvaacobs=''
if(n_elements(labelgrid) eq 0)then labelgrid=1
if(n_elements(position) eq 0)then position=[0.1,0.1,0.9,0.9]
if(n_elements(modtitle) eq 0)then modtitle=''
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(xyrel) eq 0)then xyrel=[-999,-999]
if(n_elements(lcharsize) eq 0)then lcharsize=1.0
if(n_elements(plotlegend) eq 0)then plotlegend=1
if(n_elements(fieldunits) eq 0)then fieldunits='Units'

if(n_elements(contourlevels) eq 0)then begin
  maxval=fix(alog10(max(plotdata)))
  contourlevels=[10.^(maxval-5),10.^(maxval-4),10.^(maxval-3),10.^(maxval-2),$
        10.^(maxval-1),10.^(maxval),10.^(maxval+1)]
endif

if(n_elements(contourcolors) eq 0)then begin
  ncol=n_elements(contourlevels)-1
  dcol=(!d.n_colors-50)/(max([ncol-1,1]))
  contourcolors=(!d.n_colors-10)-(indgen(ncol)*dcol)
endif

; set projection and stretch

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
if(n_elements(automap) eq 0) then automap=0
if(n_elements(origin) eq 0) then origin=[0,0]

;-----------------------------------------------------------------------
; graphics initialisation

Col0=0


; if plot crosses dateline then ensure data is 0-360

if (xmax gt 180)then begin
  index=where(lon lt 0,nindex)

  if(nindex gt 0)then begin
    lon(index)=lon(index)+360.0
    index=sort(lon(*,0))
    nshift=index(0)-2
    plotdata=Shift(plotdata,nshift,0)
    lat=shift(lat,nshift,0)
    lon=shift(lon,nshift,0)  
   endif
endif 


; grid lines
  
dxy=ymax-ymin
case 1 of 
  (dxy ge 90): gridxy=20
  (dxy ge 40) and (dxy lt 90): gridxy=10
  (dxy ge 20) and (dxy lt 40): gridxy=5
  (dxy ge 10) and (dxy lt 20): gridxy=2
  (dxy ge 5) and (dxy lt 10): gridxy=1
  (dxy ge 2) and (dxy lt 5): gridxy=0.5
  (dxy ge 1)   and (dxy lt 2): gridxy=0.2
  (dxy ge 0.5) and (dxy lt 1): gridxy=0.1
  (dxy ge 0)   and (dxy lt 0.5): gridxy=0.05
endcase


; determine plot position

genposition,multiplot,plotposition,position,$
  plotlegend,modtitle,margin=margin

; character sizes

maxplot=multiplot(1)
case 1 of 
  (maxplot ge 5): lcharsize=0.5
  (maxplot ge 3) and (maxplot lt 5): lcharsize=0.6
  (maxplot ge 2) and (maxplot lt 3): lcharsize=0.9
  (maxplot ge 1) and (maxplot lt 2): lcharsize=1.2
endcase

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

map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=plotposition,$
color=0,/noerase,name=proj_name,/hires,/isotropic

if (countries eq 1) then map_continents,/countries,/hires
if (coasts eq 1) then begin 
  if (highres eq 1) then begin
    map_continents,/coasts,/hires
  endif else begin
    map_continents,/coasts
  endelse
endif
if (states eq 1) then map_continents,/usa,/hires
map_lakes,thick=2
	
;-----------------------------------------------------------------------
; topography

if(n_elements(topog) gt 1)then begin
 tcontourlevels=[-100.,0.001,200,500,1000,2000,5000,10000]
 tcontourcolors=[27,26,25,24,17,16,7]
 pixelplot,topog,lont,latt,tcontourlevels,tcontourcolors,$
 xmin,xmax,ymin,ymax
endif

;-----------------------------------------------------------------------
; LINCOM topography

if(plotlincomtopog ne 0)then begin 
  if(n_elements(topog) gt 1)then begin
   tcontourlevels=[-0.001,50,500,1000,1500,2000,2500,3000]
   tcontourcolors=[27,26,25,24,17,16,7]
   pixelplot,topog,tlon,tlat,tcontourlevels,tcontourcolors,$
   xmin,xmax,ymin,ymax
  endif
endif

;-----------------------------------------------------------------------
; Add image

image=read_image(imagename)
tv,image,plotposition(0),plotposition(1),xsize=0.85,true=1,/normal

;-----------------------------------------------------------------------
; overplot data 

if(type eq 'pixel') then begin
  pixelplot,plotdata,lon,lat,contourlevels,contourcolors,$
   xmin,xmax,ymin,ymax
endif else if (type eq 'nofill') then begin
  contour,plotdata,reform(lon(*,0)),reform(lat(0,*)),$
    /overplot,levels=contourlevels, c_thick=3,c_colors=contourcolors
endif else begin
  contour,plotdata,reform(lon(*,0)),reform(lat(0,*)),/fill,$
    /overplot,levels=contourlevels,c_colors=contourcolors
  contour,plotdata,reform(lon(*,0)),reform(lat(0,*)),$
    /overplot,levels=contourlevels, c_thick=2
endelse 

;-----------------------------------------------------------------------

; pressure

if(n_elements(pmsl) gt 1)then begin
  contour,pmsl/100.,reform(lon(*,0)),reform(lat(0,*)),/overplot,$
   levels=952+(findgen(25)*4.0),c_colors=intarr(25)+14
endif


;-----------------------------------------------------------------------
; plot wind vectors
; map onto 50x50 grid

if(n_elements(u) gt 1 and n_elements(v) gt 1)then begin
  u=congrid(u,50,50)
  v=congrid(v,50,50)
  lon=congrid(lon,50,50)
  lat=congrid(lat,50,50)
  velovect,u,v,reform(lon(*,0)),reform(lat(0,*)),color=15,length=2,$
    /overplot
endif


;-----------------------------------------------------------------------
; plot release point

if(min(xyrel) ge -180)then begin
  plots,[xyrel(0)],[xyrel(1)],psym=7,symsize=0.7
endif  

;-----------------------------------------------------------------------
; plot lon/lat

map_grid,color=18, lats=ymin, latdel=gridxy, lons=xmin, londel=gridxy
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
; plot locations

if (plotloc ne '')then begin
  plotloc,plotloc,zoom=[xmin,xmax,ymin,ymax],placenames=placenames
endif


;-----------------------------------------------------------------------
; title

if(n_elements(modtitle) ne 0)then begin
  xyouts,(plotposition(2)+plotposition(0))/2.0,$
         (plotposition(3)+(plotposition(3)-plotposition(1))/20.0),$
	 modtitle,/normal,charsize=lcharsize,Alignment=0.5
endif

 
;-----------------------------------------------------------------------
; legend

if(plotlegend ne 0) then begin

  zoomleg=plotposition

; colours
  
  XRange=Zoomleg(2)-Zoomleg(0)
  XMinc=Zoomleg(0)+XRange*0.10
  Xinc=0.80*xrange/(N_Elements(Contourlevels)-1)

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
  
  nc=n_elements(contourlevels) 
  case 1 of 
    (nc le 5): legint=1
    (nc gt 5) and (nc le 11): legint=2
    (nc gt 11) and (nc le 15): legint=3
    (nc gt 15) and (nc le 20): legint=4
    (nc gt 20) and (nc le 25): legint=5
    (nc gt 25) : legint=6
  endcase

  if(max(contourlevels) lt 99.99 and $
     min(contourlevels) gt 0.01)then begin
     form='(f6.2)'
  endif else begin
     form='(e9.2)'
  endelse
     
  For Colour=0,N_Elements(ContourLevels)-1,legint Do Begin
     Xminp=xminc+Colour*Xinc
     XYOuts,XMinp,YMinp-0.85*yinc,Alignment=0.5,$
      StrTrim(String(ContourLevels(Colour),Format=form),2),$
      charsize=lcharsize*0.75,color=col0,/normal
  EndFor

  XYOuts,XMinc+0.48*(Xrange*0.80),YMinp+1.3*yinc,$
   Alignment=0.5,'Maximum value = '+$
   StrTrim(String(Max(PlotData),Format=form),2)+' '+$
   FieldUnits,charsize=lcharsize*0.75,color=col0,/normal

Endif 

;-----------------------------------------------------------------------
; plot observations

if (plotvaacobs ne '')then begin
  plotvaacobs,plotvaacobs,zoom=[xmin,xmax,ymin,ymax]
endif


End
