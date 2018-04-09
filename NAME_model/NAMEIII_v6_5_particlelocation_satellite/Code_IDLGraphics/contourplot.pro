pro contourplot,plotdata,lon,lat,contourlevels,contourcolors,$
 xmin,xmax,ymin,ymax

;-----------------------------------------------------------------------
; procedure to plot filled contours on map
; DBR  July 2002
;-----------------------------------------------------------------------


; ensure closed contours

lon1d=reform(lon(*,0))
lat1d=reform(lat(0,*))
dat=reform(plotdata)
sz=size(plotdata)
datc=fltarr(sz(1)+2,sz(2)+2)
datc(1,1)=plotdata
lon1d=[lon1d(0),lon1d,lon1d(sz(1)-1)]
lat1d=[lat1d(0),lat1d,lat1d(sz(2)-1)]


; filled contours

map_contour,datc,lon1d,lat1d,$
 Levels=contourlevels,c_colors=contourcolors,$
 /filled


; line contours

map_contour,datc,lon1d,lat1d,$
 Levels=contourlevels,c_colors=contourcolors*0,$
 c_thick=[2,2,2,2]     


end








