PRO DT_TO_STR,dt_var,date_string,time_string,Date_Fmt=date_fmt,Time_Fmt=time_fmt

; dt_var is a date/time variable
; Returns date string in format 'dd/mm/yyyy' if date_fmt is undefined or 2
; Returns time string in format 'hhmm' if time_fmt is undefined or -2

IF (NOT(KEYWORD_SET(date_fmt))) THEN date_fmt=2
IF (NOT(KEYWORD_SET(time_fmt))) THEN time_fmt=-2

IF (date_fmt eq 2) THEN BEGIN

  yystr=STRING(dt_var.year,format='(i4.4)')
  mmstr=STRING(dt_var.month,format='(i2.2)')
  ddstr=STRING(dt_var.day,format='(i2.2)')

  date_string=ddstr+'/'+mmstr+'/'+yystr
  
ENDIF ELSE BEGIN 

  PRINT,'Need to define dt_to_str option:',date_format
  STOP
  
ENDELSE

IF (time_fmt eq -2) THEN BEGIN

  hstr=STRING(dt_var.hour,format='(i2.2)')
  mstr=STRING(dt_var.minute,format='(i2.2)')
  
  time_string=hstr+mstr
  
ENDIF ELSE BEGIN 

  PRINT,'Need to define dt_to_str option:',time_format
  STOP
  
ENDELSE
  
END
