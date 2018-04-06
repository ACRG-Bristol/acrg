PRO vm_plottraj, optfile
  ;Program is designed to adapt plottraj for use with the IDL Virtual Machine
  ;It obtains all parameters for input into plotfield from an options file.
  ;Keyword and positional parameters may be listed in any order. The name of the file
  ;containing the parameters is either entered as an argument to the program call or
  ;may be selected from a drop-down menu
  ;N.B.1 Quotation marks are not required in optfile
  ;N.B.2 Array values should be separated by commas but not surrounded by brackets
  ;
  ;-----------------------------------------------------------------------------
  ;22 Oct 2009   sjl   created
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
  ;Error messages
  
  err2 = 'Keyword parameters should be in the format "type=pixel"'
  ;
  ;Open file containing parameters
  openr,1,optfile
  ;
  ;Loop to read header
  header = 'mmLLLLLLLLLLLLL'
  mm = strmatch(strcompress(header,/remove_all),';;PLOTTRAJINPUTS')
  while (mm eq 0) do begin
    readf, 1, header
    mm = strmatch(strcompress(header,/remove_all),';;PLOTTRAJINPUTS')
    if EOF(1) then begin
      print, 'FATAL ERROR: Inputs not found. First line of inputs should read:'
      print, '             ;;PLOTTRAJ INPUTS'
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
      'sources' : sources=splitkey[1]
      'clabel' : clabel=splitkey[1]
      'labint' : labint=splitkey[1]
      'cycle' : cycle=splitkey[1]
      'gengif' :  gengif=splitkey[1]
      'genjpg' :  genjpg=splitkey[1]
      'fit' : fit=splitkey[1]
      'projection' :  projection=fix(splitkey[1])
      'filetext' :  filetext=splitkey[1]
      'namever' :  namever=splitkey[1]
      'polar' :  polar=splitkey[1]
      'countries' :  countries=splitkey[1]
      'coasts' :  coasts=splitkey[1]
      'states' :  states=splitkey[1]
      'automap' :  automap=splitkey[1]
      'zoom' : begin
        zoom = fltarr(4)
        zoom[*] = float(strsplit(splitkey[1],',',/extract))
      end
      'origin' : begin
        origin = fltarr(2)
        origin[*] = float(strsplit(splitkey[1],',',/extract))
      end
      else:  print, 'WARNING: ', splitkey0, ' is not a recognised input keyword'
    endcase
  endwhile
  ;
  ;Check that datadir and selgrid have been included
  if (n_elements(datadir) eq 0) then begin
    print, 'FATAL ERROR: DATADIR is required for plotfield to run'
    return
  endif
  ;
  ;Use parameter information to call plotfield
  print, 'Calling plottraj'
  plottraj,datadir,sources=sources,clabel=clabel,zoom=zoom,$
    labint=labint,filetext=filetext,gengif=gengif,fit=fit,$
    projection=projection,cycle=cycle,genjpg=genjpg,namever=namever,$
    origin=origin,countries=countries,coasts=coasts,states=states,$
    automap=automap,polar=polar
    
END



