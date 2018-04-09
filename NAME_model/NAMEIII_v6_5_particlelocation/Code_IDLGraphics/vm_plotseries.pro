PRO vm_plotseries, optfile
  ;Program is designed to adapt plotfield for use with the IDL Virtual Machine
  ;It obtains all parameters for input into plotfield from an options file.
  ;Keyword and positional parameters may be listed in any order. The name of the file
  ;containing the parameters is either entered as an argument to the program call or
  ;may be selected from a drop-down menu
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
  mm = strmatch(strcompress(header,/remove_all),';;PLOTSERIESINPUTS')
  while (mm eq 0) do begin
    readf, 1, header
    mm = strmatch(strcompress(header,/remove_all),';;PLOTSERIESINPUTS')
    if EOF(1) then begin
      print, 'FATAL ERROR: Inputs not found. First line of inputs should read:'
      print, '             ;;PLOTSERIES INPUTS'
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
      'datadir' : datadir=splitkey[1]
      'selgrid' : selgrid=splitkey[1]
      'selfield' : begin
        selfield1=splitkey[1]
        selfield=strsplit(selfield1,',',/extract)
      end
      'selspecies' : begin
        selspecies1=splitkey[1]
        selspecies=strsplit(selspecies1,',',/extract)
      end
      'sellevel' : begin
        sellevel1=strmid(keyw,9,200)
        sellevel=strsplit(sellevel1,',',/extract)
      end
      'sellocation' : begin
        sellocation1=splitkey[1]
        sellocation=strsplit(sellocation1,',',/extract)
      end
      'selcase' :  selcase=splitkey[1]
      'rimnet' :  rimnet=splitkey[1]
      'timeframe' :  timeframe=splitkey[1]
      'plotlegend' :  plotlegend=splitkey[1]
      'plotheader' :  plotheader=splitkey[1]
      'multiplot' :  begin
        multiplot=[0,1,1]
        multiplot[*]= fix(strsplit(splitkey[1],',',/extract))
      end
      'genanim' :  genanim=splitkey[1]
      'gengif' :  gengif=splitkey[1]
      'genjpg' :  genjpg=splitkey[1]
      'filetext' :  filetext=splitkey[1]
      'namever' :  namever=splitkey[1]
      'log' :  log=splitkey[1]
      'tssplit' :  tssplit=splitkey[1]
      'keepgif' :  keepgif=splitkey[1]
      else:  print, 'WARNING: ', splitkey0, ' is not a recognised input keyword'
    endcase
  endwhile
  ;
  ;Check that datadir and selgrid have been included
  if (n_elements(datadir) eq 0) or (n_elements(selgrid) eq 0) then begin
    print, 'FATAL ERROR: DATADIR and SELGRID are required for plotseries to run'
    return
  endif
  ;
  ;Use parameter information to call plotfield
  print, 'Calling plotseries'
  plotseries,datadir,selgrid,selspecies=selspecies,$
    sellocation=sellocation,selfield=selfield,sellevel=sellevel,$
    filetext=filetext,multiplot=multiplot,plotheader=plotheader,$
    plotlegend=plotlegend,genanim=genanim,tssplit=tssplit,$
    rimnet=rimnet,gengif=gengif,genjpg=genjpg,keepgif=keepgif,log=log,$
    namever=namever,timeframe=timeframe,selcase=selcase
    
END



