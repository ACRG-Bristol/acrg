pro pixelplot,plotdata,lon,lat,contourlevels,contourcolors,$
 xmin,xmax,ymin,ymax

;-----------------------------------------------------------------------
; procedure to plot pixels on map
; DBR  April 2002
;-----------------------------------------------------------------------

; resolution

dlon=lon(1,0)-lon(0,0)
dlat=lat(0,1)-lat(0,0)


; Set up colour array

ArraySize=Size(PlotData)

itx=ArraySize(1)
ity=ArraySize(2)

ColourArray=intarr(itx,ity)

For i=0,n_elements(contourlevels)-2 do begin

  Index=Where((PlotData ge contourlevels(i)) and $
              (PlotData lt contourlevels(i+1)),NData) 
  If (NData gt 0) Then Begin
    ColourArray(Index)=ContourColors(i) 
  EndIf
EndFor


;-----------------------------------------------------------------------
; overplot pixels 

For ix=0,Itx-1 Do Begin
  For iy=0,Ity-1 Do Begin

    If(ColourArray(ix,iy) gt 0)then begin

       Lon1=Lon(ix,iy)-dlon/2.
       Lon2=Lon(ix,iy)+dlon/2.

       Lat1=Lat(ix,iy)-dlat/2.
       Lat2=Lat(ix,iy)+dlat/2.

       lon1=max([lon1,xmin])
       lon2=max([lon2,xmin]) 
       lon1=min([lon1,xmax])
       lon2=min([lon2,xmax])
       lat1=max([lat1,ymin])
       lat2=max([lat2,ymin]) 
       lat1=min([lat1,ymax])
       lat2=min([lat2,ymax])

       polyfill, /Data, [lon1,lon2,lon2,lon1,lon1],$
        [lat1,lat1,lat2,lat2,lat1],$
        Color=ColourArray(ix,iy)
    Endif 
  EndFor
EndFor       

end
