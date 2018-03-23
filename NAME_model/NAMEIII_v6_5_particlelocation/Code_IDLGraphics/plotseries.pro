pro plotseries,datadir,selgrid,selspecies=selspecies,$
 sellocation=sellocation,selfield=selfield,sellevel=sellevel,$
 filetext=filetext,multiplot=multiplot,plotheader=plotheader,$
 plotlegend=plotlegend,genanim=genanim,tssplit=tssplit,$
 rimnet=rimnet,gengif=gengif,genjpg=genjpg,keepgif=keepgif,log=log,$
 namever=namever,timeframe=timeframe,selcase=selcase,maxy=maxy,$
 fixdate=fixdate,yrange=yrange
;--------------------------------------------------------------------
; Procedure to generate time series plots
; 
;
; DBR Feb 2002
;
; Arguments
; required:
;  datadir      : (string)  run directory
;  selgrid      : (string)  grid ('grid1' or 'grid2')
;
; optional:
;  sellocation  : (string or string array) location(s) to plot  
;  selspecies   : (string or string array) species(s) to plot  
;  selfield     : (string or string array) field(s) to plot
;  sellevel     : (string or string array) level(s) to plot
;  selcase      : (string) NAME III case to plot
;
;  zoom         : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  filetext     : (string)  text added to filenames
;  multiplot    : intarr(3) number of plots, 
;                 (eg multiplot=[0,2,2] for 2 by 2 plots)
;                 first element can be any value
;  plotheader   : integer (0 or 1) =1 adds run information and title to plots
;  plotlegend   : integer (0 or 1) =1 plots legend below each plot
;  genanim      : integer (0 or 1) =1 generates gifs and gif animation
;  gengif       : integer (0 or 1) =1 generates gif images (on unix)
;  genjpg       : integer (0 or 1) =1 generates jpg images (on unix)
;  keepgif      : integer not used - for backward compatibility
;  log          : integer (0 or 1) =1 to use log scale
;                 note if rimnet=1 log scale will be used
;  namever      : integer indicating version of model e.g. 3 - NAME III
;  timeframe    : string - 'absolute' or 'relative' time frame for NAME III
;  maxy         : maximum value for y axis - sets yrange=[0,maxy]
;  yrange       : sets range of y-axis = [ymin, ymax] - note that yrange
;                 and maxy cannot be used together if both are
;                 specified maxy will be used
;  fixdate      : integer (0 or 1) =1 to correct placement of month string
;    
;
;-----------------------------------------------------------------------
;
; For a more full description of NAME III output format and the appropriate 
; NAME III options for standard idl plotting please refer to the NAME III 
; document: md2_4_v2(Output).
;
; NAME III output options that produce output that can be processed by
; plot series are:
; Separate File  = 'XY'  or blank 
; Across         = 'XYZ' (D also needed for quantities depending on a data grid - must not be a 'floating' D-Grid; can also include D if no D dependence)
; Output Format  = 'AZ'  2 can be used I think. I can not be used with NAME III format
;
; Example file name = Example_TS_C1.txt
; Generated with the following output options:
; Across        = XYZ
; Separate File = left blank
; Output Format = AZ
;
; Can be plotted using:
;
;selfield=['Air concentration']
;sellevel=['Boundary layer average']	
;selgrid='Example_TS'
;plotseries,'${WORKDIR}',selgrid=selgrid,selfield=selfield,sellevel=sellevel,$
;multiplot=[0,1,1],namever=3



;-----------------------------------------------------------------------
; arguments

if(n_elements(plotlegend) eq 0)then plotlegend=1
if(n_elements(plotheader) eq 0)then plotheader=1
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(gengif) eq 0) then gengif=1
if(n_elements(genjpg) eq 0) then genjpg=0
genanim=0

if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext='_'+filetext
endelse
if(n_elements(tssplit) eq 0)then tssplit='NONE'
if(n_elements(rimnet) eq 0)then rimnet=0
if(n_elements(log) eq 0)then log=1
if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(selcase) eq 0)then selcase='1'
if(n_elements(fixdate) eq 0) then fixdate=0
if(n_elements(yrange) gt 0)then yrange1=yrange

if(rimnet ne 0)then log=1
;-----------------------------------------------------------------------
; graphics initialisation

if(strpos(!version.os,'NT') ge 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

if(strpos(!version.os,'Linux') ge 0)then begin
  linux=1
endif else begin
  linux=0
endelse

xysize=1.0
axiscol=0
!p.font=0
!x.thick=1
!y.thick=1

nxplot=multiplot(1)
nyplot=multiplot(2)


if(plotheader ne 0)then begin
  position=[0.02,0.15,0.98,0.87]
endif else begin
  position=[0.0,0.0,0.98,0.98]
endelse

;-----------------------------------------------------------------------
; find all files and readheader 


; define species

readtshead,datadir,modtitle,fieldhead,selgrid=selgrid,$
           namever=namever,timeframe=timeframe,selcase=selcase

; all species

if(n_elements(selspecies) eq 0)then begin

  if (namever eq 3)then begin
    ;NAME III
    array=fieldhead(where(fieldhead(*,3)),3)
    specieslist=array[uniq(array,sort(array))] 
  endif else begin
    ;NAME
    array=fieldhead(where(fieldhead(*,4)),4)
    specieslist=array[uniq(array,sort(array))]
  endelse

  numspecies=n_elements(specieslist)
endif else begin
  numspecies=n_elements(selspecies)
  specieslist=selspecies
endelse


; locations

if(n_elements(sellocation) eq 0)then begin
  if (namever eq 3)then begin
    ;NAME III
    array=fieldhead(where(fieldhead(*,13)),13)
    locationlist=array[uniq(array,sort(array))]
  endif else begin
    ;NAME
    array=fieldhead(where(fieldhead(*,2)),2)
    locationlist=array[uniq(array,sort(array))]
  endelse

  numlocation=n_elements(locationlist)
endif else begin
  numlocation=n_elements(sellocation)
  locationlist=sellocation
endelse


; fields

if(n_elements(selfield) eq 0 or $
   n_elements(sellevel) eq 0) then begin

  if (namever eq 3)then begin
    ;NAME III
    numfields=n_elements(fieldhead(*,0))
    fieldlist=reform(fieldhead(*,2))
    levellist=reform(fieldhead(*,16))
  endif else begin
    ;NAME
    numfields=n_elements(fieldhead(*,0))
    fieldlist=reform(fieldhead(*,5))
    levellist=reform(fieldhead(*,6))
  endelse

endif else begin
  numfields=n_elements(selfield)
  fieldlist=selfield
  levellist=sellevel
endelse 


;-----------------------------------------------------------------------     
; open postscript file for printing

filetem='Plot_time_series_'+selgrid+filetext
file=datadir+delim+filetem+'.ps'
set_plot,'ps'
device,/portrait,ysize=26.7,yoffset=1.5,xsize=18,xoffset=1.5  
device,/color
device,filename=file 
stretch,20
loadct,5
tek_color
 
;-----------------------------------------------------------------------
; loop through fields and species

it=-1

for iloc=0,numlocation-1 do begin

  location=locationlist(iloc)

  for ifield=0,numfields-1 do begin 
 
    field=fieldlist(ifield)
    level=levellist(ifield)
  
    for ispec=0,numspecies-1 do begin

      it=it+1
      
      species=specieslist(ispec)

      readts,datadir,tssplit,location,$
       field,level,species,selgrid, $
       dates,dat,datunits,namever=namever,timeframe=timeframe,$
       selcase=selcase


; generate accumulations
 
      if(strpos(field,'deposition') ge 0) then begin
        dat_acc=dat
        for ih=1,n_elements(dat)-1 do begin
          dat_acc(ih)=dat_acc(ih-1)+dat(ih)
        endfor
	dat=dat_acc
      endif     


; set range

      if(rimnet gt 0)then begin

	if(strpos(field,'deposition') ge 0) then begin
          yrange=[10,1e9]
	endif     
	if(strpos(field,'concentration') ge 0) then begin
          yrange=[0.1,1e7]
	endif     
	if(strpos(field,'Dose') ge 0) then begin
          yrange=[1e-8,100]
	endif     
	
      endif else begin

          if (n_elements(yrange1) ne 0 && n_elements(maxy) eq 0) then begin
              yrange=yrange1
          endif else begin
              if(log eq 1)then begin
                  maxdat=max(dat)
                  maxdat=max([maxdat,1.e-21])
                  maxval=fix(alog10(maxdat))+1      
                  yrange=[10.0^(maxval-8),10.0^maxval]
                  if (n_elements(maxy) ne 0 ) then yrange=[0,maxy]
              endif else begin
                  maxdat=max(dat)*1.2
                  mindat=min(dat)
                  if mindat > 0.0 then begin
                     yrange=[0,maxdat]
                  endif else begin
                     yrange=[mindat,maxdat/1.2]
                  endelse
                  if (n_elements(maxy) ne 0 ) then yrange=[0,maxy]
              endelse
          endelse
	
      endelse
      
      indexlow=where(dat lt yrange(0),nlow)
      if(nlow gt 0)then dat(indexlow)=yrange(0)   


; position
     
      pos=(it) mod (nxplot*nyplot)
      multiplot=[pos,nxplot,nyplot]
      
      genposition,multiplot,plotposition,position,$
       plotlegend,location,margin=[0.1,0.0,0.06,0.0]

       
      maxplot=multiplot(1)
      case 1 of 
	(maxplot ge 5): lcharsize=0.5
	(maxplot ge 3) and (maxplot lt 5): lcharsize=0.6
	(maxplot ge 2) and (maxplot lt 3): lcharsize=0.8
	(maxplot ge 1) and (maxplot lt 2): lcharsize=1.0
      endcase
       

;-----------------------------------------------------------------------     
; set x-axis (time) 
     
     juldat=julday(dates.month,dates.day,dates.year,dates.hour,dates.minute,dates.second)
     range=abs(double(max(juldat)-min(juldat)))
     timeaxis,range
          
;-----------------------------------------------------------------------
; get round an IDL bug in plot routine when using !X.ticklayout=2 which puts the
; label in the wrong place
 
     if (fixdate eq 1) then begin
        COMMON label_date_com, cFormatArray, cMonths, cOffset, $
          cRoundup, cAmpm, cDaysWeek
        ; reformat the label to move it left by 60 characters
        cFormatArray[1]='(C(CMoA," ",CYI-60))'  
     endif   
                
;-----------------------------------------------------------------------     
; plot     

     if(log ne 0)then begin
        
       plot_io,juldat,dat,/nodata,/xstyle,/ystyle,yrange=yrange,$
        position=plotposition,ytitle=datunits,/noerase,$
        title=location+'!c'+level+' '+field,$
        charsize=lcharsize*0.7
    
       oplot,juldat,dat,color=2
       
     endif else begin
     
       plot,juldat,dat,/nodata,/xstyle,/ystyle,yrange=yrange,$
        position=plotposition,ytitle=datunits,/noerase,$
        title=location+'!c'+level+' '+field,$
        charsize=lcharsize*0.7
    
       oplot,juldat,dat,color=2
       
     endelse 
     
;-----------------------------------------------------------------------
; annotation    
   
    
      if(plotheader ne 0 and pos eq 0)then begin

       xyouts,0.5,0.98,modtitle(0),$
        /normal,alignment=0.5,charsize=1.3,color=axiscol 
       xyouts,0.5,0.94,modtitle(1),$
        /normal,alignment=0.5,charsize=1.3,color=axiscol
       xyouts,0.5,0.90,'Time series',$
 	/normal,charsize=xysize*1.2,alignment=0.5,color=axiscol

       xyouts,0.05,0.12,modtitle(4),/normal,charsize=0.9,color=axiscol  
       xyouts,0.05,0.10,modtitle(5),/normal,charsize=0.9,color=axiscol
       xyouts,0.05,0.08,modtitle(6),/normal,charsize=0.9,color=axiscol
       xyouts,0.05,0.06,modtitle(7),/normal,charsize=0.9,color=axiscol
       xyouts,0.05,0.04,modtitle(8),/normal,charsize=0.9,color=axiscol
       xyouts,0.60,0.12,'Pollutant: '+species,$
	/normal,charsize=0.9,color=axiscol
       xyouts,0.6,0.10,modtitle(3),/normal,charsize=0.9,color=axiscol
       xyouts,0.6,0.08,modtitle(2),/normal,charsize=0.9,color=axiscol


       xyouts,0.5,0.00,'Met Office Crown copyright',$
 	  /normal,alignment=0.5,charsize=0.8,color=axiscol 

       image=read_image('MO_Master_B.jpg')
       loadct,0
       tv,image,0.82,0.91,true=1,xsize=0.15,/normal
       loadct,5
       tek_color

      endif 


;-----------------------------------------------------------------------

      if(pos eq ((nxplot*nyplot)-1))then begin
        erase
      endif      

    endfor
  endfor
endfor


device,/close_file


;----------------------------------------------------------------------- 
; generate gifs


if((gengif gt 0 or genjpg gt 0) and pc eq 0)then begin
 genanim,file,datadir,filetem,gengif=gengif,genjpg=genjpg,genanim=genanim
endif

end
