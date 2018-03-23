FUNCTION VAR_TO_DT,Year,Month,Day,Hour,Minute,Second

; Returns a date/time variable with date and time set
; according to other paramenters (all integers except second)
; Minimum requirement is Year, Month and Day
; When other parameters are not set they are set to zero. 
; If array of year, month, day etc (they must have same number of elements)
; then the return will be an array of date/time variables.

sz=size(year)

if sz(0) eq 0 then begin
  nelem=1
  year=reform(year,1)
  month=reform(month,1)
  day=reform(day,1)
endif else if sz(0) eq 1 then begin
  nelem=sz(1)
  sz1=size(month)
  sz2=size(day)
  if (nelem ne sz1(1)) or (nelem ne sz2(1)) then begin
    print,'ERROR var_to_dt: Arguments must have same number of elements'
    print,nelem,sz1(1),sz2(1),' These should be equal'
    stop
  endif
endif else begin
  print,'Incorrect number of dimensions in dt variable!'
  print,'ERROR - var_to_dt
  stop
endelse
 
CASE N_PARAMS() OF

; Only year, month, day
    3: BEGIN
      Hour=replicate(0,nelem)
      Minute=replicate(0,nelem)
      Second=replicate(0.,nelem)
     END

; Only year, month, day, hour
    4: BEGIN
      if sz(0) eq 0 then hour=reform(hour,1)
      Minute=replicate(0,nelem)
      Second=replicate(0.,nelem)
     END

; Only year, month, day, hour, minute
    5: BEGIN
      if sz(0) eq 0 then begin
        hour=reform(hour,1)
        minute=reform(minute,1)
      endif
      Second=replicate(0.,nelem)
     END

; All available
    6: BEGIN
      if sz(0) eq 0 then begin
        hour=reform(hour,1)
        minute=reform(minute,1)
        second=reform(second,1)
      endif
     END

    ELSE: BEGIN
       print,'Number parameters:',n_params()
       print,'This is not allowed in var_to_dt!
       stop
     END
  ENDCASE

date=REPLICATE({dt},nelem)

FOR i=0L,nelem-1L do begin

  if (month(i) lt 1) or (month(i) gt 12) or $
     (day(i) lt 1) or (day(i) gt 31) or $
     (hour(i) lt 0) or (hour(i) gt 24) or $
     (minute(i) lt 0) or (minute(i) gt 60) or $
     (second(i) lt 0.) or (second(i) gt 60.) then begin
    print,'ERROR: var_to_dt: date/time variables out of range!'
    print,year(i),month(i),day(i),hour(i),minute(i),second(i)
    stop   
  endif
     
  date(i).year=year(i)
  date(i).month=month(i)
  date(i).day=day(i)

  add_hour=0
  if (hour(i) eq 24) then begin
    add_hour=1
    hour(i)=23
  endif  
  date(i).hour=hour(i)
  
  add_minute=0
  if (minute(i) eq 60) then begin
    add_minute=1
    minute(i)=59
  endif  
  date(i).minute=minute(i)
  
  add_sec=0
  if (second(i) eq 60.) then begin
    add_sec=1
    second(i)=59.
  endif  
  date(i).second=FLOAT(second(i))
  
  date(i).julian=JULDAY(date(i).month,date(i).day,date(i).year,$
                   date(i).hour,date(i).minute,FIX(date(i).second))
  
  if add_sec eq 1 then date(i)=dt_add(date(i),second=1.)
  if add_minute eq 1 then date(i)=dt_add(date(i),minute=1)
  if add_hour eq 1 then date(i)=dt_add(date(i),hour=1)
                    
ENDFOR

if sz(0) eq 0 then begin
  year=year(0)
  month=month(0)
  day=day(0)
  hour=hour(0)
  minute=minute(0)
  second=second(0)
endif  

RETURN, date
END                        
