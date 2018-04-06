pro plotvaacobs,file,zoom=zoom

; DBR 06/05/2002
; procedure to read location file and overplot locations on map
; AJM 19/6/2002 - Swapped lon and lat around in tslocs file
; LH 23/06/2009 format for source locations file is now the same as for 
; a NAME source locations file, e.g.:
;
;Locations: Source Locations
;Conc       ,    H-Coord,   X,    Y
;Value1, Lat-Long, -3.2000,55.9500
;Value2, Lat-Long, 2.3500,48.8666
;Value3, Lat-Long, -3.7000,40.4000
;Value4, Lat-Long, 9.1666,45.4666

; Example:
; 18/05/10 1200Z 56.5 2.9 8 d 1 !Central_estimate
;  - This would plot a square (8) centered at 56.2 lat, 2.9 lon, in colour 1 (dark blue)

; Column1 - date, not used
; Column2 - hour, not used
; Column3 - lat - eg 56.5 -  set 'd' in column 5 to use decimal
;               - eg 52°8.5N - set 'm' in column 5 to use minutes - need N/S and '°' or 'o'
; Column4 - lon - defined in same way as lat
; Column5 - 'd' for decimal, 'm' for minutes
; Column6 - psym to use - 8=square, 7=cross, 3=dot
;                       - 9 - gives a cross at current location and line to next location
;                       - 93 - gives a dot at current location and line to next location
; Column7 - colour to use to fill in square:
;         - 0 - white
;         - 1 - dark blue 
;         - 2 - mid blue 
;         - 3 - light blue 
;         - 4 - green
;         - 5 - lime
;         - 6 - yellow 
;         - 7 - light orange 
;         - 8 - deep orange 
;         - 9 - light red 
;         - 10- deep red 
;         - 11- black 


;file='/home/h03/apaj/volc_graph/aircraft_flight_04052010_12Z.txt'

readdata,file,data,nrecs,nskip=0,delim=' '
x=(data(3,*))
y=(data(2,*))
sym=(data(4,*))

sizedata=size(data,/dimensions)
units=strarr(nrecs)
units(*)='m'
if (sizedata[0]) ge 6 then units=(data(5,*))

colour=intarr(nrecs)
colour=1
if (sizedata[0]) ge 7 then colour=(data(6,*))

xx=fltarr(nrecs)
yy=fltarr(nrecs)

for i=0L,nrecs-1 do begin
  if units(i) eq 'd' then begin
    xx(i)=float(x(i))
    yy(i)=float(y(i))
  endif else begin ;units='m' or default
    
    ;x
    degpos=-1
    degpos=strpos(x(i),'o')
    if degpos eq -1 then degpos=strpos(x(i),'°')
    if degpos eq -1 then print,'degree symbol not found'
    xdeg=float(strmid(x(i),0,degpos))
    xmin=float(strmid(x(i),degpos+1,strlen(x(i)-1)))
    xx(i)=xdeg+xmin/60.
    if strmid(x(i),strlen(x(i))-1) eq 'W' then xx(i)=-xx(i)

    ;y
    degpos=-1
    degpos=strpos(y(i),'o')
    if degpos eq -1 then degpos=strpos(y(i),'°')
    if degpos eq -1 then print,'degree symbol not found'
    ydeg=float(strmid(y(i),0,degpos))
    ymin=float(strmid(y(i),degpos+1,strlen(x(i)-1)))
    yy(i)=ydeg+ymin/60.
    if strmid(y(i),strlen(y(i))-1) eq 'S' then yy(i)=-yy(i)

  endelse
    
endfor 

if(n_elements(zoom) eq 0)then begin
  dx=0.0
endif else begin
  dx=(zoom(3)-zoom(2))/200.
endelse
 
X = [-1,-1, 1, 1, -1]
Y = [-1, 1, 1,-1, -1]
USERSYM, X, Y, /fill

color_blue2red,red,green,blue
red[1]=255
green[1]=255
blue[1]=255
tvlct,red,green,blue

colourlist=[1,   $ ;white (<min)
            246, $ ;dark blue (1e-17 - 3.16e-17)
            224, $ ;mid blue (3.16e-17 - 1e-16)
            202, $ ;light blue (1e-16 - 3.16e-16)
            180, $ ;green (3.16e-16 - 1e-15)
            158, $ ;lime (1e-15 - 3.16e-15)
            136, $ ;yellow (3.16e-15 - 1e-14)
            114, $ ;light orange (1e-14, 3.16e-14)
            92,  $ ;deep orange (3.16e-14, 1e-13)
            70,  $ ;light red (1e-13, 3.16e-13)
            48,  $ ;deep red (3.16e-13 - 1e-12)
            0 ]    ;black (> max)
            
for i=0L,nrecs-1 do begin
  if(xx(i) gt zoom(0) and xx(i) lt zoom(1) and $
     yy(i) gt zoom(2) and yy(i) lt zoom(3))then begin
    ;plots,[xx(i)],[yy(i)],psym=sym(i),symsize=1.2,color=1
    ;plots,[xx(i)],[yy(i)],psym=sym(i),symsize=0.9,color=0

    if sym(i) eq 9 then begin
      if i ne nrecs-1 then begin
        plots,[xx(i),xx(i+1)],[yy(i),yy(i+1)],color=colourlist[colour(i)],thick=2
      endif
      plots,[xx(i)],[yy(i)],psym=7,color=colourlist[colour(i)]
    endif else if sym(i) eq 93 then begin
      if i ne nrecs-1 then begin
        plots,[xx(i),xx(i+1)],[yy(i),yy(i+1)],color=colourlist[colour(i)],thick=2
      endif
      plots,[xx(i)],[yy(i)],psym=3,color=colourlist[colour(i)]
    endif else begin
      plots,[xx(i)],[yy(i)],psym=sym(i),symsize=1.1,color=48
      plots,[xx(i)],[yy(i)],psym=sym(i),symsize=0.9,color=colourlist[colour(i)]
    endelse
    
  endif
endfor

end
