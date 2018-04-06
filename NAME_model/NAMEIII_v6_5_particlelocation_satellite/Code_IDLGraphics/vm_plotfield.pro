PRO vm_plotfield, optfile
  ;Program is designed to adapt plotfield for use with the IDL Virtual Machine
  ;It obtains all parameters for input into plotfield from an options file.
  ;Keyword and positional parameters may be listed in any order. The name of the file
  ;containing the parameters is entered as an argument to the program call
  ;N.B.1 Quotation marks are not required in optfile
  ;N.B.2 Array values should be separated by commas but not surrounded by brackets
  ;
  ;-----------------------------------------------------------------------------
  ;04 Sep 2009   sjl   created
  ;22 Feb 2010   sjl   altered 'IF' loop to 'CASE'
  ;-----------------------------------------------------------------------------
  ;
  ;Specify number of command line arguments and read into program
  count=strarr(1)
  optfile=(COMMAND_LINE_ARGS(count=count))[0]
  ;
  if (count eq 0) then begin
    optfile=dialog_pickfile(/read,filter='*.txt')
  endif
  ;
  ;Open file containing parameters
  openr,1,optfile
  ;
  ;Loop to read header
  header = 'mmLLLLLLLLLLLLL'
  mm = strmatch(strcompress(header,/remove_all),';;PLOTFIELDINPUTS')
  while (mm eq 0) do begin
    readf, 1, header
    mm = strmatch(strcompress(header,/remove_all),';;PLOTFIELDINPUTS')
    if EOF(1) then begin
      print, 'FATAL ERROR: Inputs not found. First line of inputs should read:'
      print, '             ;;PLOTFIELD INPUTS'
      return
    endif
  endwhile
  ;
  ;Set variables to type string
  keyw=''
  ;Read remainder of lines to end of file checking format as reading
  while ~EOF(1) do begin
    readf,1, keyw
    result=strmatch(keyw,'*=*')
    if result eq 0 then begin
      ;Check for blank lines and exit loop if found
      if strmatch(strcompress(keyw,/remove_all),'') then break
      print, 'WARNING: Keyword parameters should be in the format "keyword=A"'
      print, '         ignoring "', keyw, '"'
      continue
    endif
    ;
    ;Split input
    splitkey = strsplit(keyw,'=',/extract)
    splitkey0 = splitkey[0]
    ;
    ;apply to relevant keyword parameter
    case splitkey0 of
      'datadir': datadir=splitkey[1]
      'selgrid': selgrid=splitkey[1]
      'selfield': begin
        selfield1=splitkey[1]
        selfield=strsplit(selfield1,',',/extract)
      end
      'selspecies': begin
        selspecies1=splitkey[1]
        selspecies=strsplit(selspecies1,',',/extract)
      end
      'sellevel': begin
        sellevel1=strmid(keyw,9,200)
        sellevel=strsplit(sellevel1,',',/extract)
      end
      'seltimeaverage': begin
        seltimeaverage1=splitkey[1]
        seltimeaverage=strsplit(seltimeaverage1,',',/extract)
      end
      'seldatavalue': begin
        seldatavalue1=splitkey[1]
        seldatavalue=strsplit(seldatavalue1,',',/extract)
      end
      'selname': begin
        selname1=splitkey[1]
        selname=strsplit(selname1,',',/extract)
      end
      'selcase' :  selcase=splitkey[1]
      'back' : back=splitkey[1]
      'type' :  type=splitkey[1]
      'group' :  group=splitkey[1]
      'plotloc' :  plotloc=splitkey[1]
      'titletext' :  titletext=splitkey[1]
      'rimnet' :  rimnet=splitkey[1]
      'hires' :  hires=splitkey[1]
      'exact' :  exact=splitkey[1]
      'plottopog' :  plottopog=splitkey[1]
      'plotpmsl' :  plotpmsl=splitkey[1]
      'plotuv' :  plotuv=splitkey[1]
      'time' :  time=splitkey[1]
      'plotlegend' :  plotlegend=splitkey[1]
      'plotheader' :  plotheader=splitkey[1]
      'multiplot' :  begin
        multiplot=[0,1,1]
        multiplot[*]= fix(strsplit(splitkey[1],',',/extract))
      end
      'genanim' :  genanim=splitkey[1]
      'gengif' :  gengif=splitkey[1]
      'genjpg' :  genjpg=splitkey[1]
      'projection' :  projection=fix(splitkey[1])
      'ncontours' :  ncontours=splitkey[1]
      'filetext' :  filetext=splitkey[1]
      'template' :  template=splitkey[1]
      'colourtable' :  colourtable=splitkey[1]
      'ct_blue2red' :  ct_blue2red=splitkey[1]
      'percent' :  percent=splitkey[1]
      'last' :  last=splitkey[1]
      'xpixel' :  xpixel=splitkey[1]
      'contourstep' :  contourstep=splitkey[1]
      'scalefactor' :  scalefactor=splitkey[1]
      'namever' :  namever=splitkey[1]
      'timeframe' :  timeframe=splitkey[1]
      'plotlincomtopog' :  plotlincomtopog=splitkey[1]
      'relduration' :  relduration=splitkey[1]
      'metoffice' :  metoffice=splitkey[1]
      'polar' :  polar=splitkey[1]
      'aar' :  aar=splitkey[1]
      'countries' :  countries=splitkey[1]
      'coasts' :  coasts=splitkey[1]
      'states' :  states=splitkey[1]
      'placenames' :  placenames=splitkey[1]
      'ContourGrid' :  ContourGrid=splitkey[1]
      'automap' : automap=splitkey[1]
      'zoom' : begin
        zoom = fltarr(4)
        zoom[*] = float(strsplit(splitkey[1],',',/extract))
      end
      'origin' : begin
        origin = fltarr(2)
        origin[*] = float(strsplit(splitkey[1],',',/extract))
      end
      'contours' : begin
        nfields = max([1,n_elements(selfield)])
        nspecies = max([1,n_elements(selspecies)])
        clevs = float(strsplit(splitkey[1],',',/extract))
        contours = fltarr(nfields,nspecies,n_elements(clevs))
        for i = 1,nspecies do begin
          for j = 1,nfields do begin
            contours[j-1,i-1,*] = clevs
          endfor
        endfor
      end
      else: print, 'WARNING: ', splitkey0, ' is not a recognised input keyword'
    endcase
  endwhile
  ;
  ;Check that datadir and selgrid have been included
  if (n_elements(datadir) eq 0) or (n_elements(selgrid) eq 0) then begin
    print, 'FATAL ERROR: DATADIR and SELGRID are required for plotfield to run'
    return
  endif
  ;
  ;Use parameter information to call plotfield
  print, 'Calling plotfield'
  plotfield,datadir,selgrid,selspecies=selspecies,$
    selfield=selfield,sellevel=sellevel,zoom=zoom,filetext=filetext,$
    multiplot=multiplot,plotheader=plotheader,plotlegend=plotlegend,$
    genanim=genanim,group=group,ct_blue2red=ct_blue2red,$
    plotuv=plotuv,plotpmsl=plotpmsl,plottopog=plottopog,$
    exact=exact,highres=highres,rimnet=rimnet,template=template,$
    plotloc=plotloc,type=type,gengif=gengif,genjpg=genjpg,$
    projection=projection,ncontours=ncontours,back=back,$
    percent=percent,last=last,xpixel=xpixel,contourstep=contourstep,$
    scalefactor=scalefactor,namever=namever,timeframe=timeframe,$
    titletext=titletext,SelTimeAverage=SelTimeAverage,$
    plotlincomtopog=plotlincomtopog,relduration=relduration,$
    seldatavalue=seldatavalue,selname=selname,usecontours=contours,$
    metoffice=metoffice,colourtable=colourtable,selcase=selcase,$
    polar=polar,aar=aar,countries=countries,coasts=coasts,states=states,$
    placenames=placenames,origin=origin,ContourGrid=ContourGrid,automap=automap
    
  ;base = WIDGET_BASE(XOFFSET=300,YOFFSET=300)
  ;button1 = WIDGET_BUTTON(base, VALUE='Program VM_PLOTFIELD Complete',$
  ;  UVALUE='Done', XSIZE=400, YSIZE=100)
  ;WIDGET_CONTROL, base, /REALIZE
  ;wait, 3
  ;WIDGET_CONTROL, base, /DESTROY
  
END



