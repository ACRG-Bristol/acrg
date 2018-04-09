pro plotloc,file,zoom=zoom,placenames=placenames

; DBR 06/05/2002
; procedure to read location file and overplot locations on map
; AJM 19/6/2002 - Swapped lon and lat around in tslocs file
; LH 23/06/2009 format for source locations file is now the same as for 
; a NAME source locations file, e.g.:
;
;Locations: Source Locations
;Name       ,    H-Coord,   X,    Y
;EDINBURGH, Lat-Long, -3.2000,55.9500
;PARIS, Lat-Long, 2.3500,48.8666
;MADRID, Lat-Long, -3.7000,40.4000
;MILAN, Lat-Long, 9.1666,45.4666

IF (NOT(KEYWORD_SET(placenames))) THEN placenames=0

readdata,file,data,nrecs,nfields=4,nskip=2

data=strtrim(data,2)
place=data(0,*)
x=float(data(2,*))
y=float(data(3,*))

if(n_elements(zoom) eq 0)then begin
  dx=0.0
endif else begin
  dx=(zoom(3)-zoom(2))/200.
endelse
 
for i=0,nrecs-1 do begin
  if(x(i) gt zoom(0) and x(i) lt zoom(1) and $
     y(i) gt zoom(2) and y(i) lt zoom(3))then begin
    plots,[x(i)],[y(i)],psym=1,symsize=1
    if placenames eq 1 then xyouts,x(i)+dx,y(i)-dx,place(i),charsize=1,color=14
  endif
endfor

end
