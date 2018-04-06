pro plotattrib,datadir,selgrid,sellocation=sellocation,selspecies=selspecies,$
 selfield=selfield,sellevel=sellevel,zoom=zoom,filetext=filetext,$
 multiplot=multiplot,plotheader=plotheader,plotlegend=plotlegend,$
 genanim=genanim,contours=contours,smooth=smooth,gengif=gengif,genjpg=genjpg,$
 projection=projectionpolar=polar,origin=origin,countries=countries,$
 coasts=coasts,states=states,automap=automap

  
;--------------------------------------------------------------------
; Procedure to plot attribution data
;
; DBR March 2002
;
; 12/07/2002 - added linear contours for travel time
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
;
;  zoom         : fltarr(4) area to plot (eg zoom=[-20,30,35,65])
;  filetext     : (string)  text added to filenames
;  multiplot    : intarr(3) number of plots, 
;                 (eg multiplot=[0,2,2] for 2 by 2 plots)
;                 first element can be any value
;  plotheader   : integer (0 or 1) =1 adds run information and title to plots
;  plotlegend   : integer (0 or 1) =1 plots legend below each plot
;  genanim      : integer (0 or 1) =1 generates gifs and gif animation
;  contours     : fltarr  array of contour values
;  smooth       : integer size of smoothing window (pixels)
;  gengif       : integer (0 or 1) =1 generates gif images (on unix)
;  genjpg       : integer (0 or 1) =1 generates jpg images (on unix)
;  projection   : integer (1 - 16) - projection for mapping (for options see
;                 IDL map option) #8 (cylindrical) is standard
;  polar        : integer, 0 for non-polar projection, -1 for south pole, 
;                 1 for north pole, projection modes 1 (stereographic)
;                 and 2 (orthographic) are recommended, but it will work with other
;                 projection modes.  The polar keyword overrides the origin keyword.
;  origin       : 2-element array [lat,long] to be used for origin of map projection.  
;                 If this keyword is not set, the automap option will be enabled and the
;                 map origin will be the centre of the map itself, (unless the polar 
;                 keyword is set).
;  automap      : integer (0 or 1) =1 map projection origin determined automatically
;  countries    : set to 0 to not plot national borders
;  coasts       : set to 1 to plot lakes and extra islands.  N.B. do not use this
;                 with the 'hires' keyword on continent-sized maps, and do not use
;                 it without 'hires' set on anything smaller (however it will plot
;                 every lake, however small).  You will need to use it without 'hires' 
;                 for getting the Great Lakes on maps of the whole of N. America, and 
;                 the Caspian Sea on large-scale maps of Asia, and with 'hires' to
;                 show the smaller islands around the UK.
;  states       : set to 0 to not plot individual states on maps of the USA
 
;  
; note:
; (1) selspecies,selfield and sellevel must contain same number of requests,
;     otherwise all attribution fields plotted
;
; example calls:
;  plotattrib,'/home/fr1100/apnm/namev66/testv66',$
;   'grid2',sellocation='ASCOT',selspecies='TRACER',$
;    selfield='Air concentration',sellevel='Boundary layer',$
;    multiplot=[0,2,2],plotheader=1,plotlegend=1,/genanim
;
;  plotattrib,'/home/fr1100/apnm/namev66/testv66',$
;  'grid2',multiplot=[0,1,1],plotheader=1,plotlegend=1,/genanim
;
;
; subroutine calls:
;  readattribhead
;  readattrib
;  maplimits
;  genposition
;  plotfieldpos
;  genanim
;
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
; argument checks

if(n_elements(plotlegend) eq 0)then plotlegend=1
if(n_elements(plotheader) eq 0)then plotheader=1
if(n_elements(multiplot) eq 0)then multiplot=[0,1,1]
if(n_elements(genanim) eq 0)then genanim=1
if(n_elements(singleps) eq 0)then singleps=0
if(n_elements(smooth) eq 0)then smooth=0
if(n_elements(gengif) eq 0) then gengif=0
if(n_elements(genjpg) eq 0) then genjpg=0
if(n_elements(projection) eq 0) then projection=8
if(n_elements(polar) eq 0) then polar=0
if(n_elements(countries) eq 0) then countries=1
if(n_elements(coasts) eq 0) then coasts=0
if(n_elements(states) eq 0) then states=1
if(n_elements(automap) eq 0) then begin
  if(n_elements(origin) eq 0) then begin
    automap=1
    origin=[0,0]
  endif else begin
    automap=0
  endelse
endif


if(n_elements(filetext) eq 0)then begin
  filetext=''
endif else begin
  filetext=filetext+'_'
endelse


;-----------------------------------------------------------------------
; graphics initialisation

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

axiscol=0
colors=[3,18,5,10] 
!p.font=0
!x.thick=1
!y.thick=1
nxplot=multiplot(1)
nyplot=multiplot(2)


; plot limits

if(plotheader ne 0)then begin
  position=[0.02,0.15,0.98,0.87]
endif else begin
  position=[0.0,0.0,1.0,1.0]
endelse


;-----------------------------------------------------------------------
; determine number of locations and fields
; if not passed determine from files in datadir and attribution header


; location

if (n_elements(sellocation) eq 0)then begin

  fieldfiles=findfile(datadir+delim+'Attrib_'+selgrid+'*')
  numlocs=n_elements(fieldfiles)
  
  if(pc eq 1)then begin
    for i=0,numlocs-1 do begin
     fieldfiles(i)=datadir+delim+fieldfiles(i)    
    endfor
  endif
  
  locationlist=strarr(numlocs)
  for it=0,numlocs-1 do begin
    pos1=strpos(fieldfiles(it),'Attrib')
    pos2=strpos(fieldfiles(it),'.')
    pos1=pos1+13
    locationlist(it)=strmid(fieldfiles(it),pos1,pos2-pos1) 
   print,locationlist(it)
  endfor

endif else begin

  numlocs=n_elements(sellocation)
  locationlist=sellocation
endelse


; field

if(n_elements(selfield) eq 0 or $
   n_elements(sellevel) eq 0 or $
   n_elements(selspecies) eq 0) then begin

  readattribhead,fieldfiles(0),modtitle,attribhead

  numfields=n_elements(attribhead(*,0))
  fieldlist=reform(attribhead(*,3))
  levellist=reform(attribhead(*,5))
  specieslist=reform(attribhead(*,1))

endif else begin

  numfields=n_elements(selfield)
  fieldlist=selfield
  levellist=sellevel
  specieslist=selspecies
  
endelse 


;-----------------------------------------------------------------------
; loop through locations and fields

for iloc=0,numlocs-1 do begin

  sellocation=locationlist(iloc)


  for ifield=0,numfields-1 do begin

    selfield=fieldlist(ifield)
    sellevel=levellist(ifield)
    selspecies=specieslist(ifield)


;-----------------------------------------------------------------------
; read attribution matrix

    readattrib,datadir,selgrid,sellocation,selspecies,selfield,sellevel, $
      date,lon1,lat1,matrixarr,head1,modtitle,ntime

    tot=reform(total(matrixarr,d=[0]))

;-----------------------------------------------------------------------
; map limits

    if(n_elements(zoom) eq 0)then begin
      zoom=fltarr(4)
      maplimits,lon1,lat1,tot,xmin,xmax,ymin,ymax
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

    allfiles=['file']


;-----------------------------------------------------------------------
; determine contours

  if(n_elements(contours) eq 0)then begin
    if(strpos(selfield,'travel') ge 0)then begin
      maxdat=fix(max([matrixarr]))+1.0
      contours1=(indgen(8)+1)*(maxdat)/8.0
    endif else begin
      maxdat=max([matrixarr])
      maxdat=max([maxdat,1.e-12])
      maxval=fix(alog10(maxdat))
      contours1=10.0^(maxval-indgen(8)+1)
      contours1=contours1(sort(contours1))
    endelse
  endif else begin
    contours1=contours
  endelse
  

;-----------------------------------------------------------------------     
; open postscript file for printing

    file=datadir+delim+'Attribplot_'+filetext+selgrid+'_'+sellocation+'_'+$
     strcompress(selfield,/remove_all)+'_'+$
     strcompress(sellevel,/remove_all)+'_'+$     
     selspecies+'.ps'
    allfiles=[allfiles,file] 
    set_plot,'ps'
    device,/portrait,ysize=25,yoffset=2  
    device,/color
    device,filename=file
    loadct,5  
    stretch,25
    tek_color


;-----------------------------------------------------------------------        
; loop through times

    for it=0,ntime-1 do begin
 
      dat1=reform(matrixarr(it,*,*))
 
      yy=date(it).year
      mm=date(it).month
      dd=date(it).day
      hh=date(it).hour
      mn=date(it).minute
      filetime=string(yy,mm,dd,hh,mn,format='(i4.4,4i2.2)')
      timetext=string(hh,mn,dd,mm,yy,$
       format='(i2.2,'':'',I2.2,'' '',i2.2,''/'',i2.2,''/'',i4.4)')


;-----------------------------------------------------------------------     
; plot fields

      pos=(it) mod (nxplot*nyplot)
      multiplot=[pos,nxplot,nyplot]
 
      if(smooth ne 0)then begin
        dat1=smooth(dat1,smooth)
      endif
      
      plotfieldpos,lon1,lat1,dat1,position=position,$
        zoom=zoom,contourlevels=contours1,contourcolors=contourcolors,$
        fieldunits=head1(4),multiplot=multiplot,modtitle=timetext,$
        plotlegend=plotlegend,projection=projection,polar=polar,$
	countries=countries,coasts=coasts,states=states,$
	origin=origin,automap=automap


;-----------------------------------------------------------------------   
; annotation

      if(plotheader ne 0 and pos eq 0)then begin
  
        xyouts,0.5,0.98,modtitle(0),$
         /normal,alignment=0.5,charsize=1.3,color=axiscol 
        xyouts,0.5,0.94,'Attribution for '+head1(0),$
         /normal,charsize=1.2,alignment=0.5,color=axiscol
        xyouts,0.5,0.90,head1(5)+' '+head1(2)+' '+head1(3),$
         /normal,charsize=1.2,alignment=0.5,color=axiscol

        xyouts,0.05,0.12,modtitle(4),/normal,charsize=0.9,color=axiscol  
        xyouts,0.05,0.10,modtitle(5),/normal,charsize=0.9,color=axiscol
        xyouts,0.05,0.08,modtitle(6),/normal,charsize=0.9,color=axiscol
        xyouts,0.05,0.06,modtitle(7),/normal,charsize=0.9,color=axiscol
        xyouts,0.05,0.04,modtitle(8),/normal,charsize=0.9,color=axiscol
        xyouts,0.60,0.12,'Pollutant: '+head1(1),$
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
; printer

      if(pos eq ((nxplot*nyplot)-1))then begin
        erase
      endif
     

    endfor
  
    device,/close_file


;----------------------------------------------------------------------- 
; generate gifs

    if((gengif eq 1 or genjpg eq 1 or genanim eq 1) and pc eq 0)then begin
      filetem='Attribplot_'+filetext+selgrid+'_'+sellocation+'_'+$
       strcompress(selfield,/remove_all)+'_'+$
       strcompress(sellevel,/remove_all)+'_'+$      
       selspecies
      genanim,allfiles(1:*),datadir,filetem,gengif=gengif,genjpg=genjpg,$
       genanim=genanim
    endif

  endfor
endfor

end




