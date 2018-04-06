pro plotvertpos,lont,latt,plotdata, $
 zoom=zoom,contourlevels=contourlevels,contourcolors=contourcolors,$
 fieldunits=fieldunits,plotlegend=plotlegend,$
 position=position,hrel=hrel,zrel=zrel,$
 multiplot=multiplot,modtitle=modtitle,labelgrid=labelgrid,$
 exact=exact,highres=highres,$
 plotloc=plotloc,type=type,$
 placenames=placenames

;-----------------------------------------------------------------------
; procedure to draw contour plot of a vertical slice of data
;
; SJL August 2009
; (modified plotfieldpos.pro to create)
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
;  hrel           : source position (latitude or longitude)
;  zrel            : source position (height)
;  projection      : integer (1 - 16) - projection for mapping (for options see
;                    IDL map option) #8 (cylindrical) is standard
;  multiplot       : intarr(3) number of plots, 
;                    (eg multiplot=[0,2,2] for 2 by 2 plots)
;                    first element can be any value
;  modtitle        : 
;  labelgrid       :
;  exact           : integer (0 or 1) =1 map area based on lon/lat limits
;  plotloc         : string filename conatining locations to overplot,
;                    assumes same format as tslocs.txt file
;  type            : string ('pixel' or 'contour') choose pixels or 
;                    filled contours
;  placenames   : set to 0 to not plot names of locations in plotloc file 


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

if(n_elements(type) eq 0)then type='pixel'
if(n_elements(exact) eq 0)then exact=0
if(n_elements(zoom) eq 0)then begin
  zoom=fltarr(4)
;Rewrite
;  maplimits,lont,latt,plotdata,xmin,xmax,ymin,ymax,exact=exact
  xmin = min(lon)
  xmax = max(lon)
  ymin = min(lat)
  ymax = max(lat)
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
if(n_elements(labelgrid) eq 0)then labelgrid=1
if(n_elements(position) eq 0)then position=[0.1,0.1,0.9,0.9]
if(n_elements(modtitle) eq 0)then modtitle=''
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(hrel) eq 0)then hrel=-999
if(n_elements(zrel) eq 0)then zrel=-999
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
  dcol=(!d.n_colors-50)/(ncol-1)
  contourcolors=(!d.n_colors-10)-(indgen(ncol)*dcol)
endif

;-----------------------------------------------------------------------
; graphics initialisation

Col0=0


; if plot crosses dateline then ensure data is 0-360
; this should work since when lon is actually latitude it shouldn't
; exceed 180

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
;Allows for extra space required by axis labelling
margin=[0.1,0.1,0.0,0.05]

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
;
;-----------------------------------------------------------------------
; set up plot

!P.MULTI = multiplot[1]
plot,[xmin,xmax], [ymin,ymax], /nodata, position=plotposition

;----------------------------------------------------------------------
; overplot data 

if(type eq 'pixel') then begin
  pixelplot,plotdata,lon,lat,contourlevels,contourcolors,$
   xmin,xmax,ymin,ymax
endif else begin
  contourplot,plotdata,lon,lat,contourlevels,contourcolors,$
   xmin,xmax,ymin,ymax
endelse 

;-----------------------------------------------------------------------
; plot release point
oplot,[hrel],[zrel],psym=7,symsize=0.7

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

;NEED TO MOVE LEGEND to allow for axis labelling

if(plotlegend ne 0) then begin

    zoomleg = plotposition
; y-shifted
    legshift = 0.01
    zoomleg[3] = zoomleg[3]-legshift
    zoomleg[1] = zoomleg[1]-legshift

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


End
