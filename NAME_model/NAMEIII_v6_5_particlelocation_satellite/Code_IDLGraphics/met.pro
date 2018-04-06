pro met,input,output

for i=0,453 do begin
  for j=0,414 do begin
    ;output(i,j)=input(0,i,j)*65536L+input(1,i,j)*256L+input(2,i,j)
    if input(0,i,j)gt 250 then begin
      output(i,j)=16777215
    endif else if input(0,i,j)lt 10 then begin
      output(i,j)=0
    endif else begin
      output(i,j)=3407837
    endelse
  endfor
endfor

end
