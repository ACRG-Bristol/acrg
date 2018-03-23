pro genposition,multiplot,plotposition,position,$
 plotlegend,modtitle,margin=margin

;-----------------------------------------------------------------------
;
; generate  position for multiple plots
; DBR March 2002
;
;-----------------------------------------------------------------------

if(n_elements(margin) eq 0)then margin=[0.,0.,0.,0.]
plotnum=multiplot(0)
nx=multiplot(1)
ny=multiplot(2)
plotposition=fltarr(4)


; space above and below plot for annotation

if(plotlegend ne 0)then begin
  ybotmar=0.19+margin(2)
endif else begin
  ybotmar=0.07+margin(2)
endelse  

if(modtitle ne '')then begin
  ytopmar=0.08+margin(3)
 endif else begin
  ytopmar=0.05+margin(3)
endelse  

if(max([nx,ny]) gt 1)then begin
  xleftmar=0.02+margin(0)
  xrightmar=0.02+margin(1)
 endif else begin
  xleftmar=0.02+margin(0)
  xrightmar=0.02+margin(1)
endelse 


; overall area for plots

if(n_elements(position) eq 0)then begin
  xleft=0.0
  xright=1.0
  ytop=1.0
  ybot=0.0
endif else begin
  xleft=position(0)
  xright=position(2)
  ytop=position(3)
  ybot=position(1)
endelse
    

; position of plot

iy=fix((plotnum)/nx)
ix=(plotnum)-(iy*nx)
  
dx=(xright-xleft)/nx
dy=(ytop-ybot)/ny

plotposition(0)=xleft+(dx*(ix))+(xleftmar*dx)
plotposition(2)=xleft+(dx*(ix+1))-(xrightmar*dx)
plotposition(1)=ytop-(dy*(iy+1))+(ybotmar*dy)
plotposition(3)=ytop-(dy*(iy))-(ytopmar*dy)     
 
end

