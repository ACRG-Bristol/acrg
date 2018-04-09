pro distance, d,lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2

; Calculates distance between two points on the earth's surface.  
; Assumes spherical earth
;
; Lois Huggett
; 20/03/09

r=double(6371.0)
pi=double(3.1415926535)

; Check for appropriate angles

if (lat1 gt 90) OR (lat1 lt -90) OR (lat2 gt 90) OR (lat2 lt -90) OR $
(lon1 gt 360) OR (lon1 lt -180) OR (lon2 gt 360) OR (lon2 lt -180) then begin
  print, 'Distance fatal error: coordinate out of range'
  retall
endif 

phi1=pi/180.0*(90.0-double(lat1))
phi2=pi/180.0*(90.0-double(lat2))
theta1=pi/180.0*double(lon1)
theta2=pi/180.0*double(lon2)

d=r*acos(sin(phi1)*cos(theta1)*sin(phi2)*cos(theta2)$
  +sin(phi1)*sin(theta1)*sin(phi2)*sin(theta2)+cos(phi1)*cos(phi2))  

end 
