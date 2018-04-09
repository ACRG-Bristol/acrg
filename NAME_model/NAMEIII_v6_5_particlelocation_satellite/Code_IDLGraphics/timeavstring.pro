pro timeavstring,endtime,avstring,outstring,namever=namever

; Procedure to make a validity time string for time averaged output
; Lois Huggett
; 19/05/10

if (n_elements(namever) eq 0) then namever=2

if (((strpos(avstring,'no') lt 0) AND (strpos(avstring,'No') lt 0 )AND (strpos(avstring,'NO') lt 0)) $
   AND (((strpos(avstring,'aver') ge 0) OR (strpos(avstring,'AVER') ge 0)))) then begin
   
   
  if (strpos(avstring,'hr') ge 0) then begin
    pattern='hr' 
  endif else if (strpos(avstring,'HR') ge 0) then begin
    pattern='HR'
  endif else if (strpos(avstring,'hour') ge 0) then begin
    pattern='hour'
  endif else if (strpos(avstring,'HOUR') ge 0) then begin
    pattern='HOUR'
  endif else begin
    outstring='Valid at '+strtrim(endtime,2)
    return
  endelse

  timestring=strsplit(avstring,pattern,/regex,/extract)
  timestring=strsplit(timestring(0),' ',/extract)
  timestring=timestring(n_elements(timestring)-1)  
  
  timeav=double(timestring)
  
  datetimestruct
  
  date1=replicate({datetime},1)
  date2=replicate({datetime},1)
  
  endtime=strtrim(endtime,2)
  date1.hour=strmid(endtime,0,2)
  date1.minute=strmid(endtime,2,2)
  date1.day=strmid(endtime,8,2)
  date1.month=strmid(endtime,11,2)
  date1.year=strmid(endtime,14,4)
  
  date1.julian=julday(date1.month,date1.day,date1.year,date1.hour,date1.minute,date1.second)
  date2.julian=date1.julian-timeav/24.0D
  caldat,date2.julian,month,day,year,hour,minute
  date2.month=month
  date2.day=day
  date2.year=year
  date2.hour=hour
  date2.minute=minute
 
  DT_TO_STR,date2,date_string,time_string,Date_Fmt=2,Time_Fmt=-2
  outstring='Valid from '+time_string+'UTC '+date_string+' to '+endtime

endif else begin
  outstring='Valid at '+strtrim(endtime,2)
endelse


end 
