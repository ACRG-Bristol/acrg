pro metlogo,image

image=read_image('MO_Master_B.jpg')
help,image
 for j=0,453 do begin
  for k=0,414 do begin
   if image(0,j,k) gt 180 and image(0,j,k) lt 230 then begin 
     print,'1',image(0,j,k),image(1,j,k),image(2,j,k)
     ;image(0,j,k)=14
     ;image(1,j,k)=219
     image(2,j,k)=185
     print,'2',image(0,j,k),image(1,j,k),image(2,j,k)
   endif
  endfor
 endfor
end
