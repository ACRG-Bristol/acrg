pro plotrsmc,datadir,gridname,zoom=zoom,singleps=singleps,gengif=gengif,$
 filetext=filetext,exact=exact,genjpg=genjpg,starttime=starttime
 

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
; starttime   :(integer) != 0 to show beginning of integration time on plots
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
if(n_elements(starttime) eq 0)then starttime=0

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

selgrid=[gridname,gridname,gridname,gridname,gridname,gridname]
selfield=['Dosage','Dosage','Dosage','Total deposition','Total deposition','Total deposition']
sellevel=['From     0 -   500m agl',$
          'From     0 -   500m agl',$
	  'From     0 -   500m agl',$
	  'Boundary layer',$
	  'Boundary layer',$
	  'Boundary layer']
psfile=['Dosage1','Dosage2','Dosage3','Total_deposition1','Total_deposition2','Total_deposition3']
nplots=6


;-----------------------------------------------------------------------
; find all files and determine number if times and species

fieldfiles=findfile(datadir+delim+'Fields_'+gridname+'*')

numfiles=n_elements(fieldfiles)

if(pc eq 1)then begin
  for i=0,numfiles-1 do begin
   fieldfiles(i)=datadir+delim+fieldfiles(i)    
  endfor
endif

seldate={datetime,year:0,month:0,day:0,hour:0,minute:0,second:0,julian:0.0D}

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
seltimeindex=[index(0),index(1),index(2),index(0),index(1),index(2)]

print,'seltimeindex',seltimeindex
print,'index',index

;------------------------------------------------------------------------
; determined species


readfieldhead,fieldfiles(0),modtitle,fieldhead

array=fieldhead(where(fieldhead(*,1)),1)
specieslist=array[uniq(array,sort(array))]

rate=float(strmid(modtitle(6),21,12))

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
  maplimits,lon,lat,tot,xmin,xmax,ymin,ymax,exact=exact
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


; open postscript file for printing

    if(singleps eq 0)then begin
      file=datadir+delim+'RSMC_'+filetext+gridname+'_'+species+'_'+psfile(it)+'.ps'
      allfiles=[allfiles,file]
      set_plot,'ps'
      ;device,papername='A4'
      device,/portrait,ysize=26.7,xsize=18,yoffset=1.5,xoffset=1.5
      ;device,Path_Points=1500,/color
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
    
    ; set minimum contour at 10E-20
    if (maxval lt -12) then maxval=-12
      
    ncontours=5
    contourstep=100
      
    contours1=fltarr(ncontours)
    contours1(0)=10.0^maxval
    for ic=1,ncontours-1 do begin
     contours1(ic)=contours1(ic-1)/contourstep
    endfor 
    contours1=contours1(sort(contours1))

    multiplot=[0,1,1]
    plotfieldpos,lon,lat,dat, position=position,type='pixel',$
     zoom=zoom,contourlevels=contours1,labelgrid=1,$
     xyrel=xyrel,fieldunits=head(4),multiplot=multiplot,$
     plotlegend=plotlegend,contourcolors=[240,220,180,100]

    map_grid, color=0, glinethick=2.0


;-----------------------------------------------------------------------
; annotation    
    
    if (strmatch(strcompress(head(3),/remove_all),'Dosage')) then begin
        althead3=head(2)+' Air Concentration'
    endif else begin
        althead3=head(3)
    endelse
    xyouts,0.5,0.96,althead3,$
     /normal,charsize=xysize*1.4,alignment=0.5,color=axiscol


    if (strmatch(strcompress(head(3),/remove_all),'Totaldeposition')) then begin
        althead5='Surface'
    endif else begin
        althead5=head(5)
    endelse
    xyouts,0.5,0.92,althead5,$
     /normal,charsize=xysize*1.4,alignment=0.5,color=axiscol
 
    
    if (starttime ne 0) then begin
    
      pattern='hr'
      timestring=strsplit(head(2),pattern,/regex,/extract)
      timestring=strsplit(timestring(0),' ',/extract)
      timestring=timestring(n_elements(timestring)-1)  
      
      timeint=double(timestring)
  
      date1=replicate({datetime},1)
      date2=replicate({datetime},1)
  
      endtime=strtrim(head(6),2)
      date1.hour=strmid(endtime,0,2)
      date1.minute=strmid(endtime,2,2)
      date1.day=strmid(endtime,8,2)
      date1.month=strmid(endtime,11,2)
      date1.year=strmid(endtime,14,4)
  
      date1.julian=julday(date1.month,date1.day,date1.year,date1.hour,date1.minute,date1.second)
      date2.julian=date1.julian-timeint/24.0D
      caldat,date2.julian,month,day,year,hour,minute
      date2.month=month
      date2.day=day
      date2.year=year
      date2.hour=hour
      date2.minute=minute
 
      DT_TO_STR,date2,date_string,time_string,Date_Fmt=2,Time_Fmt=-2
 
      postitle = 'Time integrated from '+time_string+'UTC '+date_string+' to '+endtime
      
    endif else begin
    
      postitle='Valid at '+head(6)
      
    endelse
 
    xyouts,0.5,0.88,postitle,$
     /normal,charsize=xysize*1.4,alignment=0.5,color=axiscol  

    if(rate gt 4.00e-5 and rate lt 5.00e-5)then begin
      xyouts, 0.5,0.84,'Existence or strength of event is unknown',$
      /normal,alignment=0.5,color=axiscol,charsize=xysize*1.3 
      xyouts, 0.5,0.81,'Results based on default initial values',$
      /normal,alignment=0.5,color=axiscol,charsize=xysize*1.3 

    endif

    xyouts,0.05,0.26,'Contour values may change from chart to chart',$
      /normal,charsize=xysize,color=axiscol  
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
