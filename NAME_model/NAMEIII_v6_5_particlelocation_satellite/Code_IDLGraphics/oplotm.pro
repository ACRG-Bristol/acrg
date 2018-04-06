pro oplotm,x,y, $
 linestyle=linestyle,color=color,psym=psym,missing=missing,maxdx=maxdx,$
 thick=thick

; procedure to oplot data with missing data
; works with datetime, real or integer data
; dbr 09/96


s=size(x)

if(n_elements(linestyle) eq 0)then linestyle=!p.linestyle
if(n_elements(color) eq 0)then color=!p.color
if(n_elements(psym) eq 0)then psym=!p.psym
if(n_elements(missing) eq 0)then missing=0.0
if(n_elements(maxdx) eq 0)then maxdx=1.0e99
if(n_elements(thick) eq 0)then thick=1

if(s(n_elements(s)-2) lt 8)then begin

  for n=0l,n_elements(y)-2 do begin
    if(y(n) ne missing and y(n+1) ne missing and $
       x(n+1)-x(n) le maxdx)then begin

      oplot,[x(n),x(n+1)],[y(n),y(n+1)], $
            psym=psym,linestyle=linestyle,color=color,thick=thick
    endif
  endfor

endif else begin

  for n=0l,n_elements(y)-2 do begin
    duration=x(n+1).julian-x(n).julian
    if(y(n) ne missing and y(n+1) ne missing and $
       duration le maxdx/24.0)then begin

      oplot,[x(n),x(n+1)],[y(n),y(n+1)], $
            psym=psym,linestyle=linestyle,color=color,thick=thick
    endif
  endfor

endelse



end
