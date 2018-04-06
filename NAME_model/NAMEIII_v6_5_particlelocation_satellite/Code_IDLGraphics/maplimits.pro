pro maplimits,lon,lat,tot,xmin,xmax,ymin,ymax,$
  meridcross=meridcross,datecross=datecross,exact=exact
  
; procedure to determine plotting limits
; dbr 26/9/00
;
; returns limits within range 0 to 360 (centred on dateline) 
; or -180 to 180 (centred on meridian)


; arguments

if(n_elements(exact) eq 0)then exact=0

dx=lon(1,0)-lon(0,0)

 
; determine if plume crosses meridian or dateline

index=where(tot gt 0.0 and (lon gt 175.0 or lon lt -175.0),nindex)
if(nindex gt 0)then begin
  print,'maplimits: crosses dateline'
  datecross=1
endif else begin
  datecross=0
endelse

index=where(tot gt 0.0 and lon gt -5.0 and lon lt 5.0,nindex)
if(nindex gt 0)then begin
  meridcross=1
  print,'maplimits: crosses meridian'    
endif else begin
  meridcross=0
endelse


; plot from -180 to +180 unless crosses dateline only

if(datecross eq 1 and meridcross eq 0)then begin
  index=where(lon lt 0.0,nindex)
  if (nindex NE 0) then lon(index)=lon(index)+360.0
endif


; limits

if(max(tot) ne 0 and exact eq 0)then begin
  xmin=min(lon(where(tot gt 0.0)))
  xmax=max(lon(where(tot gt 0.0)))
  ymin=min(lat(where(tot gt 0.0)))
  ymax=max(lat(where(tot gt 0.0)))
endif else begin
  xmin=min(lon)
  xmax=max(lon)
  ymin=min(lat)
  ymax=max(lat)
endelse


; balance x and y ranges

if(exact eq 0)then begin

  xrange=xmax-xmin
  yrange=ymax-ymin
  xmid=xmin+xrange/2
  ymid=ymin+yrange/2

  xfactor=1.6
  
  xyrange=max([25.0*dx,xrange/xfactor,yrange,2.0])*1.1

  ymin=min([ymin,ymid-xyrange*0.5]) 
  ymax=max([ymax,ymid+xyrange*0.5])
  xmin=min([xmin,xmid-xyrange*xfactor*0.5])
  xmax=max([xmax,xmid+xyrange*xfactor*0.5])
  
endif

; shift whole area to be within global limits

if(datecross eq 1 and meridcross eq 0)then begin

  If(xmax gt 360)then begin
     xshift=xmax-360.0
     xmax=xmax-xshift
     xmin=xmin-xshift
  endif 
  If(xmin lt 0.0)then begin
     xshift=xmin
     xmax=xmax-xshift
     xmin=xmin-xshift
  endif 

  xmax=min([xmax,360.0])  
  xmin=max([xmin,0.0])  

endif else begin

  If(xmax gt 180)then begin
     xshift=xmax-180.0
     xmax=xmax-xshift
     xmin=xmin-xshift
  endif 
  If(xmin lt -180.0)then begin
     xshift=xmin+180
     xmax=xmax-xshift
     xmin=xmin-xshift
  endif 
   xmax=min([xmax,180.0])  
  xmin=max([xmin,-180.0])
endelse

If(ymax gt 90)then begin
   yshift=ymax-90.0
   ymax=ymax-yshift
   ymin=ymin-yshift
endif 
If(ymin lt -90.0)then begin
   yshift=ymin+90
   ymax=ymax-yshift
   ymin=ymin-yshift
endif 

ymax=min([ymax,90.0])  
ymin=max([ymin,-90.0])

print,'maplimits: ',xmin,xmax,ymin,ymax
end
