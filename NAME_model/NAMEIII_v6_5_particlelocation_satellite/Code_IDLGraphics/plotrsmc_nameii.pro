pro plotrsmc,datadir,gridname,zoom=zoom,singleps=singleps,gengif=gengif,$
 filetext=filetext,exact=exact,genjpg=genjpg
 

;--------------------------------------------------------------------
; Procedure to generate rsmc plots
; 
; Reads all fields in a given timestep, extracts air con, dosage, wet
; and total deposition. Generates annotated 4 up plot
; first and last fields read first to generate map limits and contours
;
; DBR Feb 2002
;
; Arguments:
; 
; Datatdir    :(string) run directory
; gridname    :(string) 'grid1' or 'grid2' - grid to use ie 
; singleps    :(integer) =1 to generate one ps file per species, 
;                        =0 for separate ps files per species and timestep
; zoom        :(fltarr(4) plot area
; filetext    :(string) text to add to filenames
; gengif      :(integer) =1 to generate gif images
;                        =0 for postscript output only
; exact       :(integer) =1 to fit map area to grid
; genjpg      :(integer) =1 to generate jpg images
;
;
;
;
;
; Changes
;  29/05/2002 DBR set highres=0 (in plotfieldpos call), 
;   removed projection block, increased gif image size, labelgrid=0
;  11/07/2002 DBR add contour plot option, fix colours
;  25/07/2002 DBR changed gif resolution
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
; arguments

if(n_elements(exact) eq 0)then exact=0
if(n_elements(gengif) eq 0)then gengif=1
if(n_elements(axiscol) eq 0)then axiscol=0
if(n_elements(singleps) eq 0)then singleps=0
if(n_elements(genjpg) eq 0)then genjpg=0

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse


;-----------------------------------------------------------------------
; graphics

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

xysize=1.0
axiscol=0
colors=[3,18,5,10] 
!p.font=0
!x.thick=1
!y.thick=1
position=[0.02,0.25,0.98,0.84]


;-----------------------------------------------------------------------
; required fields

selgrid=[gridname,gridname,gridname,gridname]
selfield=['Dosage','Dosage','Dosage','Total deposition']
sellevel=['From     0 -   500m agl',$
          'From     0 -   500m agl',$
	  'From     0 -   500m agl',$
	  'Boundary layer']
psfile=['Dosage1','Dosage2','Dosage3','Total_deposition']
nplots=4


;-----------------------------------------------------------------------
; find all files and determine number if times and species

fieldfiles=findfile(datadir+delim+'Fields_'+gridname+'*')
numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

seldate={datetime,year:0,month:0,day:0,hour:0,minute:0,second:0}

seltime=strarr(numfiles)
seldate=replicate({datetime},numfiles)
for it=0,numfiles-1 do begin
  seltime(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12) 
  seldate(it).year=strmid(seltime(it),0,4)
  seldate(it).month=strmid(seltime(it),4,2)
  seldate(it).day=strmid(seltime(it),6,2)
  seldate(it).hour=strmid(seltime(it),8,2)
  seldate(it).minute=strmid(seltime(it),10,2)
endfor

; select midday and midnight output

index=where(seldate.hour eq 12 or seldate.hour eq 0)
seltimeindex=[index(1),index(3),index(5),index(5)]


;------------------------------------------------------------------------
; determined species


readfieldhead,fieldfiles(0),modtitle,fieldhead

array=fieldhead(where(fieldhead(*,1)),1)
specieslist=array[uniq(array,sort(array))]

rate=float(strmid(modtitle(6),21,9))


;-----------------------------------------------------------------------
; loop through all fields to find map area

for ispec=0,n_elements(specieslist)-1 do begin

  species=specieslist(ispec)

  for it=0,nplots-1 do begin
 
    readfield,datadir,seltime(seltimeindex(it)),selgrid(it),$
      selfield(it),sellevel(it),species,lon,lat,dat,head,xyrel=xyrel
   
    if(it eq 0 and ispec eq 0)then begin
      tot=dat
    endif else begin
      tot=tot+dat
    endelse
    
  endfor
  
endfor


;-----------------------------------------------------------------------
; map limits

if(n_elements(zoom) eq 0)then begin
  zoom=fltarr(4)
  maplimits,lon,lat,dat,xmin,xmax,ymin,ymax,exact=exact
  zoom(0)=xmin
  zoom(1)=xmax
  zoom(2)=ymin
  zoom(3)=ymax
endif else begin
  xmin=zoom(0)
  xmax=zoom(1)
  ymin=zoom(2)
  ymax=zoom(3)
endelse  


;-----------------------------------------------------------------------
; loop through all species

for ispec=0,n_elements(specieslist)-1 do begin

  species=specieslist(ispec)
  allfiles=[species]

; open postscript file for printing

  if(singleps eq 1)then begin
    file=datadir+delim+'RSMC_'+filetext+gridname+'_'+species+'.ps'
    allfiles=[allfiles,file]
    set_plot,'ps'
    device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
    device,/color
    device,filename=file
    loadct,0 
    stretch,32
    tek_color
  endif


; loop through times

  for it=0,nplots-1 do begin
  
    readfield,datadir,seltime(seltimeindex(it)),selgrid(it),$
    selfield(it),sellevel(it),species,lon,lat,dat,head
 
 
; subtract previous dosage

    if(it eq 1 or it eq 2)then begin

      readfield,datadir,seltime(seltimeindex(it-1)),selgrid(it-1),$
       selfield(it-1),sellevel(it-1),species,lonold,latold,datold,headold

      dat=dat-datold

    endif   


; open postscript file for printing

    if(singleps eq 0)then begin
      file=datadir+delim+'RSMC_'+filetext+gridname+'_'+species+'_'+psfile(it)+'.ps'
      allfiles=[allfiles,file]
      set_plot,'ps'
      device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
      device,/color
      device,filename=file
      loadct,0  
      stretch,32
      tek_color
    endif

    backcol=!d.n_colors
    axiscol=0


; ------------------------------------------------------------------
; contours

    maxdat=max(dat)*1.2
    maxval=fix(alog10(maxdat))
    if (maxval ge 0)then maxval=maxval+1
      
    ncontours=5
    contourstep=100
      
    contours1=fltarr(ncontours)
    contours1(0)=10.0^maxval
    for ic=1,ncontours-1 do begin
     contours1(ic)=contours1(ic-1)/contourstep
    endfor 
    contours1=contours1(sort(contours1))

    multiplot=[0,1,1]
    plotfieldpos,lon,lat,dat, position=position,type='contour',$
     zoom=zoom,contourlevels=contours1,labelgrid=0,$
     xyrel=xyrel,fieldunits=head(4),multiplot=multiplot,$
     plotlegend=plotlegend,contourcolors=[240,220,180,100]


;-----------------------------------------------------------------------
; annotation    

    xyouts,0.5,0.96,head(3)+' '+head(5),$
     /normal,charsize=xysize*1.4,alignment=0.5,color=axiscol
 
    xyouts,0.5,0.92,'Valid at '+head(6),$
     /normal,charsize=xysize*1.4,alignment=0.5,color=axiscol  

    if(rate gt 4.00e-5 and rate lt 5.00e-5)then begin
      xyouts, 0.5,0.88,'Existence or strength of event is unknown',$
      /normal,alignment=0.5,color=axiscol,charsize=xysize*1.3 
      xyouts, 0.5,0.85,'Results based on default initial values',$
      /normal,alignment=0.5,color=axiscol,charsize=xysize*1.3 

    endif

    xyouts,0.05,0.22,'Simulation Description',charsize=1.2,/normal,color=axiscol

    xyouts,0.05,0.18,modtitle(4),/normal,charsize=xysize,color=axiscol  
    xyouts,0.05,0.16,modtitle(5),/normal,charsize=xysize,color=axiscol
    xyouts,0.05,0.14,modtitle(6),/normal,charsize=xysize,color=axiscol
    xyouts,0.05,0.12,modtitle(7),/normal,charsize=xysize,color=axiscol
    xyouts,0.05,0.10,modtitle(8),/normal,charsize=xysize,color=axiscol
    xyouts,0.55,0.18,'Pollutant: '+head(1),/normal,charsize=xysize,color=axiscol
    xyouts,0.55,0.16,modtitle(3),/normal,charsize=xysize,color=axiscol
    xyouts,0.55,0.14,modtitle(2),/normal,charsize=xysize,color=axiscol

    xyouts,0.5,-0.02,'Met Office Crown copyright',$
	/normal,alignment=0.5,charsize=xysize*0.8,color=axiscol 

    image=read_image('MO_Master_B.jpg')
    loadct,0
    tv,image,0.82,0.91,true=1,xsize=0.15,/normal
    tek_color

    if(singleps eq 0)then begin     
      device,/close_file
    endif else begin
      erase
    endelse     


  endfor

  if(singleps eq 1)then begin
    device,/close_file
  endif


;----------------------------------------------------------------------- 
; generate gifs

  filetem='none'
  if((gengif eq 1 or genjpg eq 1) and pc eq 0)then begin
    genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,$
      genanim=0
  endif


endfor

end
