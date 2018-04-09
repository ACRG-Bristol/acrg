pro mapcalc,datr,sdat,xyrel,nrecs,xdim,ydim,margin,pacific,ll,tr,ie
    
;  procedure to calculate the lower left and top right corners of the
;  map for area at risk plots from NAME gridpoints
;
;  Lois Huggett 24/03/09

;  datr    : array of gridpoints, longitude, latitude, value and distance
;  sdat    : indices of datr sorted by distance
;  xyrel   : source coords in lat/long
;  nrecs   : numbers of records in datr
;  xdim    : x dimension of map in km
;  ydim    : y dimension of map in km
;  margin  : fractional border for map
;  pacific : flag for swapping to 0-360 longitudes
;  ll,tr   : (output) arrays with lower left and top right gridpoints of map to output
;  ie      : (output) index of highest priority point excluded from map


ie=long(1000000000)
ll=dblarr(2)
tr=dblarr(2)

if xyrel(0) lt -900 then begin
  xyrel(0)=datr(sdat(1),2)
  xyrel(1)=datr(sdat(1),3)
endif

if pacific eq 1 AND xyrel(0) lt 0 then xyrel(0)=xyrel(0)+360.0

ll(0)=xyrel(0) ; min x point
ll(1)=xyrel(1) ; min y point
tr(0)=xyrel(0) ; max x point
tr(1)=xyrel(1) ; max y point

for i=long(1),nrecs-1 do begin

  x=datr(sdat(i),2)
  y=datr(sdat(i),3)
  
  if pacific eq 1 AND x lt 0 then x=x+360.0
  
  xmin=min([x,ll(0)])
  xmax=max([x,tr(0)])
  ymin=min([y,ll(1)])
  ymax=max([y,tr(1)])

  distance,xsize,lon1=xmin,lat1=xyrel(1),lon2=xmax,lat2=xyrel(1)
  distance,ysize,lon1=xyrel(0),lat1=ymin,lon2=xyrel(0),lat2=ymax

  if (xsize lt xdim) AND (ysize lt ydim) then begin
    ll(0)=xmin
    ll(1)=ymin
    tr(0)=xmax
    tr(1)=ymax
  endif else begin
    ie=min([ie,i])
    i=nrecs-1
  endelse  
  
endfor

;-----------------------------------------------------------------------
; increase map size if map is smaller than desired

distance,dx,lon1=xyrel(0)-0.5,lat1=xyrel(1),lon2=xyrel(0)+0.5,lat2=xyrel(1)
distance,dy,lon1=xyrel(0),lat1=xyrel(1)-0.5,lon2=xyrel(0),lat2=xyrel(1)+0.5

ll(0)=ll(0)-(xdim-xsize)/2.0/dx
ll(1)=ll(1)-(ydim-ysize)/2.0/dy
tr(0)=tr(0)+(xdim-xsize)/2.0/dx
tr(1)=tr(1)+(ydim-ysize)/2.0/dy

;----------------------------------------------------------------------- 
; add margins

ll(0)=ll(0)-margin*xdim/dx
ll(1)=ll(1)-margin*ydim/dy
tr(0)=tr(0)+margin*xdim/dx
tr(1)=tr(1)+margin*ydim/dy

;-----------------------------------------------------------------------
; check latitudes between -90 and +90

if ll(1) lt -90 then begin
  ll(1)=-90.0
  print, 'mapcalc warning: minimum latitude set to South Pole'
endif 
if tr(1) gt 90 then begin
  tr(1)=90.0
  print, 'mapcalc warning: maximum latitude set to North Pole'
endif

;-----------------------------------------------------------------------
; check longitudes between -180 and +360

if ll(0) lt -180 then ll(0)=ll(0)+360
if ll(0) gt 360 then ll(0)=ll(0)-360
if tr(0) lt -180 then tr(0)=tr(0)+360
if tr(0) gt 360 then tr(0)=tr(0)-360

;-----------------------------------------------------------------------

end
