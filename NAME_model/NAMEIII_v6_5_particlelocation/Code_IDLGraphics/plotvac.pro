pro plotvac,datadir,gridname,zoom=zoom,singleps=singleps,filetext=filetext,$
 colht=colht,gengif=gengif,exact=exact,sumht=sumht,projection=projection,$
 genjpg=genjpg,polar=polar,origin=origin,automap=automap

;--------------------------------------------------------------------
; Procedure to generate VAC plots
; 
; Reads all fields in a given timestep, extracts air concentrations 
; at each required level. 
; first and last fields read first to generate map limits and contours
;
; DBR Feb 2002
; CSW Mar 2006 Generates txt file for GIS with critical concentration value 
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
; genjpg      :(integer) =1 to generate jpg images
; exact       :(integer) =1 to fit map area to grid
; polar       :(integer) =1 for north pole-centred map, -1 for south
; origin      :(2-element floating-point array, [lat,long]) set to specify
;               origin point for map projection
; automap     :(integer) =1 for grid-centred origin of map projection
;
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
; arguments

if(n_elements(exact) eq 0)then exact=0
if(n_elements(gengif) eq 0)then gengif=1
if(n_elements(genjpg) eq 0)then genjpg=0
if(n_elements(singleps) eq 0)then singleps=0

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse

if(n_elements(colht) eq 0)then colht=0.0
if(n_elements(projection) eq 0)then projection=8
if(n_elements(polar) eq 0)then polar=0
if(n_elements(placenames) eq 0) then placenames=1
if(n_elements(origin) eq 0) then begin
  origin=[0,0]
  if (n_elements(automap) eq 0) then automap=-1
endif else if (n_elements(origin) ne 2) then begin
  print,'origin keyword must be specified as origin=[lat,long]'
  return
endif else begin
  if (origin(0) lt -90) OR (origin(0) gt 90) then begin
    print,'origin latitude must be between -90 and +90'
    return
  endif else if (origin(1) lt -180) OR (origin(1) gt 360) then begin
    print, 'origin longitude must be between -180 and +360'
    return
  endif
endelse
if(n_elements(automap) eq 0) then automap=0

;-----------------------------------------------------------------------
; threshold concentrations
; taken from VAFTAD (NOAA)

a=fltarr(23,10)

a(22,*)=[1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20]
a(21,*)=[1e-19,1e-19,1e-19,1e-19,1e-19,1e-19,1e-19,1e-19,1e-19,1e-19]
a(20,*)=[1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18]
a(19,*)=[1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18]
a(18,*)=[1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18]
a(17,*)=[1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-18]
a(16,*)=[1e-18,1e-18,1e-18,1e-18,1e-18,1e-18,1e-17,1e-17,1e-17,1e-17]
a(15,*)=[1e-18,1e-18,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17]
a(14,*)=[1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-16]
a(13,*)=[1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-17,1e-16,1e-16,1e-16]
a(12,*)=[1e-17,1e-17,1e-17,1e-17,1e-17,1e-16,1e-16,1e-16,1e-16,1e-16]
a(11,*)=[1e-17,1e-17,1e-17,1e-16,1e-16,1e-16,1e-16,1e-16,1e-16,1e-15]
a(10,*)=[1e-17,1e-16,1e-16,1e-16,1e-16,1e-16,1e-16,1e-15,1e-15,1e-15]
a(9,*)=[1e-16,1e-16,1e-16,1e-16,1e-16,1e-16,1e-15,1e-15,1e-15,1e-15]
a(8,*)=[1e-16,1e-16,1e-16,1e-16,1e-16,1e-15,1e-15,1e-15,1e-15,1e-15]
a(7,*)=[1e-16,1e-16,1e-16,1e-16,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(6,*)=[1e-16,1e-16,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(5,*)=[1e-16,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(4,*)=[1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(3,*)=[1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(2,*)=[1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(1,*)=[1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]
a(0,*)=[1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15]

colhtrange=[0,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,$
  22000,24000,26000,28000,30000,32000,34000,36000,38000,40000,$
  50000,70000,200000]*0.3048

sumhtrange=[0,2000,4000,6000,8000,10000,12000,14000,16000,18000,200000]*0.3048

colindex=where(colhtrange gt colht)
sumindex=where(sumhtrange gt sumht)

critconc=a(colindex(0)-1,sumindex(0)-1)
print,'Summit height (m): ',sumht
print,'Column height (m): ',colht
print,'Critical concentration (g/m3): ',critconc


;-----------------------------------------------------------------------
; delimeters

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse


;-----------------------------------------------------------------------
; create output text file with critical concentration used 

outfile=datadir+delim+'ash_threshold.txt'

openw,fout,outfile,/get_lun  

printf,fout,'Critical concentration for visual ash (g/m3):'
printf,fout,critconc

close,fout
free_lun,fout


;-----------------------------------------------------------------------
; graphics initialisation

xysize=1.0
axiscol=0
colors=[3,18,5,10] 
!p.font=0
!x.thick=1
!y.thick=1
stretch=0

a=findgen(8)*(!pi*2/8)
usersym,cos(a),sin(a),/fill
ssize=0.5
spsym=8


;-----------------------------------------------------------------------
; required fields

selgrid=[gridname,gridname,gridname,gridname,gridname,gridname]
selfield=['Air concentration','Air concentration','Air concentration',$
          'Air concentration','Air concentration','Air concentration']
sellevel=['From FL350 - FL550','From FL200 - FL350','From FL100 - FL200',$
          'From FL030 - FL100','From FL015 - FL030','From FL000 - FL015']

;-----------------------------------------------------------------------
; find all files and determine number if times and species

fieldfiles=findfile(datadir+delim+'Fields_'+gridname+'*')
numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

seltime=strarr(numfiles)
for it=0,numfiles-1 do begin
  seltime(it)=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12) 
endfor

readfieldhead,fieldfiles(0),modtitle,fieldhead,xyrel
array=fieldhead(where(fieldhead(*,1)),1)
specieslist=array[uniq(array,sort(array))]


;-----------------------------------------------------------------------
; loop through species

for ispec=0,n_elements(specieslist)-1 do begin

  species=specieslist(ispec)
  allfiles=[species]

;-----------------------------------------------------------------------
; read first and last for contours and map limits
 
  if(n_elements(zoom) eq 0)then begin

    readfield,datadir,seltime(0),selgrid(0),selfield(0),sellevel(0),$
      species,lon11,lat11,dat11,head11
    readfield,datadir,seltime(0),selgrid(1),selfield(1),sellevel(1),$
      species,lon12,lat12,dat12,head12
    readfield,datadir,seltime(0),selgrid(2),selfield(2),sellevel(2),$
      species,lon13,lat13,dat13,head13
    readfield,datadir,seltime(0),selgrid(3),selfield(3),sellevel(3),$
      species,lon14,lat14,dat14,head14
    readfield,datadir,seltime(0),selgrid(4),selfield(4),sellevel(4),$
      species,lon15,lat15,dat15,head15
    readfield,datadir,seltime(0),selgrid(5),selfield(5),sellevel(5),$
      species,lon16,lat16,dat16,head16


    readfield,datadir,seltime(numfiles-1),selgrid(0),selfield(0),$
      sellevel(0),species,lon21,lat21,dat21,head21  
    readfield,datadir,seltime(numfiles-1),selgrid(1),selfield(1),$
      sellevel(1),species,lon22,lat22,dat22,head22  
    readfield,datadir,seltime(numfiles-1),selgrid(2),selfield(2),$
      sellevel(2),species,lon23,lat23,dat23,head23  
    readfield,datadir,seltime(numfiles-1),selgrid(3),selfield(3),$
      sellevel(3),species,lon24,lat24,dat24,head24  
    readfield,datadir,seltime(numfiles-1),selgrid(4),selfield(4),$
      sellevel(4),species,lon25,lat25,dat25,head25  
    readfield,datadir,seltime(numfiles-1),selgrid(5),selfield(5),$
      sellevel(5),species,lon26,lat26,dat26,head26  


    dat11(where(dat11 le critconc,n1))=0.0
    dat12(where(dat12 le critconc,n2))=0.0
    dat13(where(dat13 le critconc,n3))=0.0
    dat14(where(dat14 le critconc,n4))=0.0
    dat15(where(dat15 le critconc,n5))=0.0
    dat16(where(dat16 le critconc,n6))=0.0
    
    dat21(where(dat21 le critconc,n1))=0.0
    dat22(where(dat22 le critconc,n2))=0.0
    dat23(where(dat23 le critconc,n3))=0.0
    dat24(where(dat24 le critconc,n4))=0.0
    dat25(where(dat25 le critconc,n5))=0.0
    dat26(where(dat26 le critconc,n6))=0.0

    tot=dat11+dat12+dat13+dat14+dat15+dat16+$
        dat21+dat22+dat23+dat24+dat25+dat26
 

;-----------------------------------------------------------------------
; map limits

    zoom=fltarr(4)
    maplimits,lon11,lat11,tot,xmin,xmax,ymin,ymax,exact=exact
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

  dxy=ymax-ymin
  case 1 of 
    (dxy ge 90): gridxy=30
    (dxy ge 50) and (dxy lt 90): gridxy=20
    (dxy ge 20) and (dxy lt 50): gridxy=10
    (dxy ge 10) and (dxy lt 20): gridxy=5
    (dxy ge 5) and (dxy lt 10): gridxy=2
    (dxy ge 0) and (dxy lt 5): gridxy=1
  endcase
  res=gridxy
  
  

stretch=0


;-----------------------------------------------------------------------     
; open postscript file for printing

  if(singleps eq 1)then begin
    file=datadir+delim+'VAC_'+filetext+gridname+'_'+species+'.ps'
    allfiles=[allfiles,file]
    set_plot,'ps'
    device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
    device,/color
    device,filename=file
    loadct,0  
    stretch,20
    tek_color
  endif
  
  
;-----------------------------------------------------------------------        
; loop through times

  for it=0,numfiles-1 do begin

    filename=fieldfiles(it)
    filetime=strmid(fieldfiles(it),(strpos(fieldfiles(it),'Field'))+13,12) 


; read all fields for timestep

    readfield,datadir,seltime(it),selgrid(0),selfield(0),sellevel(0),$
      species,lon1,lat1,dat1,head1
    readfield,datadir,seltime(it),selgrid(1),selfield(1),sellevel(1),$
      species,lon2,lat2,dat2,head2
    readfield,datadir,seltime(it),selgrid(2),selfield(2),sellevel(2),$
      species,lon3,lat3,dat3,head3
    readfield,datadir,seltime(it),selgrid(3),selfield(3),sellevel(3),$
      species,lon4,lat4,dat4,head4
    readfield,datadir,seltime(it),selgrid(4),selfield(4),sellevel(4),$
      species,lon5,lat5,dat5,head5
    readfield,datadir,seltime(it),selgrid(5),selfield(5),sellevel(5),$
      species,lon6,lat6,dat6,head6

    
    dat1(where(dat1 le critconc,n1))=0.0
    dat2(where(dat2 le critconc,n1))=0.0
    dat3(where(dat3 le critconc,n1))=0.0
    dat4(where(dat4 le critconc,n1))=0.0
    dat5(where(dat5 le critconc,n1))=0.0
    dat6(where(dat6 le critconc,n1))=0.0

    lon=lon1
    lat=lat1
    
    
;-----------------------------------------------------------------------     
; open postscript file for printing

    if(singleps eq 0 and it mod 2 eq 0)then begin
      file=datadir+delim+'VAC_'+filetext+gridname+'_'+species+'_'+filetime+'.ps'
      allfiles=[allfiles,file]      
      set_plot,'ps'
      device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
      device,/color
      device,filename=file
      loadct,0  
      stretch,20
      tek_color
    endif
    
    contourlevels=[critconc*0.1,1.0e21]
    contourcolors=[14]

;-----------------------------------------------------------------------     
; shift if plotting centred on dateline

  if(xmax gt 180)then begin
    index=where(lon lt 0.0)
    lon(index)=lon(index)+360.0
  endif



; plot location

  position=fltarr(4,4)

  if(it mod 2  eq 0 )then begin
    position(0,*)=[0.05,0.75,0.4,0.95]
    position(1,*)=[0.05,0.55,0.4,0.75]
    position(2,*)=[0.05,0.35,0.4,0.55]
    position(3,*)=[0.05,0.1 ,0.4,0.3 ]
    erase
    noerase=0
  endif else begin
    position(0,*)=[0.6,0.75,0.95,0.95]
    position(1,*)=[0.6,0.55,0.95,0.75]
    position(2,*)=[0.6,0.35,0.95,0.55]
    position(3,*)=[0.6,0.1,0.95,0.3]
    noerase=1
  endelse

;-----------------------------------------------------------------------
;If necessary convert projection number to a projection name

projection_name,projection=projection,proj_name=proj_name

;-----------------------------------------------------------------------
; set projection origin

if (automap eq 1) then begin
  origin(0)=(ymin+ymax)/2.0
  origin(1)=(xmin+xmax)/2.0
endif else if (xmax-origin(1) gt 180) OR (xmin-origin(1) lt -180) then begin
  origin(1)=(xmin+xmax)/2.0
endif else if (automap eq -1) then begin
  origin[0] = 0
  origin(1)=(xmin+xmax)/2.0
endif 

; set up polar projection if requested (overrides other mapping specifications)
  if (polar eq 1) then begin
    xmin=0
    xmax=360
    if (ymin lt 0) then ymin=0
    ymax=90
    origin(0)=90
    origin(1)=0
  endif else if (polar eq -1) then begin
    xmin=0
    xmax=360
    ymin=-90
    if (ymax gt 0) then ymax=0
    origin(0)=-90
    origin(1)=0
  endif
  
;-----------------------------------------------------------------------
; plot map 
; i=0: FL350 to FL550
; i=1: FL200 to FL350
; i=2: surface to FL200
; i=3: surface to FL550

  for i=0,3 do begin

    map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=position(i,*),$
           /continents,/noerase,name=proj_name,/isotropic,/grid,$
	   color=3,/hires,latdel=gridxy,londel=gridxy	 
    map_continents,/countries,color=3,/hires
      
    if i eq 0 then pixelplot,dat1,lon,lat,contourlevels,contourcolors,$
                   xmin,xmax,ymin,ymax
    if i eq 1 then pixelplot,dat2,lon,lat,contourlevels,contourcolors,$
                   xmin,xmax,ymin,ymax
    if i eq 2 then pixelplot,dat3+dat4+dat5+dat6,lon,lat,contourlevels,$
                   contourcolors,xmin,xmax,ymin,ymax
    if i eq 3 then pixelplot,dat1+dat2+dat3+dat4+dat5+dat6,$
                   lon,lat,contourlevels,contourcolors,$
                   xmin,xmax,ymin,ymax
    
    map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=position(i,*),$
           /continents,/noerase,name=proj_name,/isotropic,/grid,$
  	   color=3,/hires,latdel=gridxy,londel=gridxy	 
    map_continents,/countries,color=3,/hires
    map_lakes,color=3,thick=2
     
    map_set,origin(0),origin(1),limit=[ymin,xmin,ymax,xmax],position=position(i,*),$
           /noerase,name=proj_name,/isotropic,color=0	 
     
    plots,[xyrel(0)],[xyrel(1)],psym=2,symsize=2,color=axiscol     

  endfor



;-----------------------------------------------------------------------
; annotation    

  if(it mod 2 eq 0)then begin
    xyouts,0.5,0.85,'FL550!cFL350',/normal,alignment=0.5,charsize=xysize*1.2,color=axiscol 
    xyouts,0.5,0.65,'FL350!cFL200',/normal,alignment=0.5,charsize=xysize*1.2,color=axiscol 
    xyouts,0.5,0.45,'FL200!cSURFACE',/normal,alignment=0.5,charsize=xysize*1.2,color=axiscol 
    xyouts,0.5,0.2,'FL550!cSURFACE!c(composite)',/normal,alignment=0.5,$
       charsize=xysize*1.2,color=axiscol 
    xyouts,0.25,0.325,'Valid '+head1(6),/normal,alignment=0.5,charsize=xysize*1.1,color=axiscol 

  endif else begin
    xyouts,0.75,0.325,'Valid '+head1(6),/normal,alignment=0.5,charsize=xysize*1.1,color=axiscol 
  endelse
  
  xyouts,0.5,0.97,'VOLCANIC ASH GRAPHICS',$
    /normal,alignment=0.5,charsize=xysize*1.5,color=axiscol 
  
  xyouts,0.05,0.08,modtitle(4),/normal,charsize=xysize*0.9,color=axiscol  
  xyouts,0.05,0.06,modtitle(5),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.05,0.04,modtitle(6),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.05,0.02,modtitle(7),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.05,0.00,modtitle(8),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.6,0.08,'Pollutant: '+head1(1),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.6,0.06,modtitle(3),/normal,charsize=xysize*0.9,color=axiscol
  xyouts,0.6,0.04,modtitle(2),/normal,charsize=xysize*0.9,color=axiscol

  xyouts,0.6,0.02,'Threshold conc: '+$
        string(critconc,format='(e10.2)')+'g/m!e3!n',/normal,$
	charsize=xysize*0.9,color=axiscol 

  xyouts,0.6,0.0,'Summit height: '+$
        string(fix(sumht),format='(i6)')+'m',/normal,$
	charsize=xysize*0.9,color=axiscol 

  xyouts,0.5,-0.02,'Met Office Crown copyright',$
	/normal,alignment=0.5,charsize=xysize*0.7,color=axiscol 

  image=read_image('MO_Master_B.jpg')
  loadct,0
  tv,image,0.82,0.93,true=1,xsize=0.12,/normal
  tek_color

;-----------------------------------------------------------------------
    if(singleps eq 0 and it mod 2 eq 1)then begin     
      device,/close_file
    endif    

    if (singleps eq 1 and it mod 2 eq 1)then begin
      erase
    endif
  
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

