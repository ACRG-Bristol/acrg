FUNCTION DT_ADD,dt_time,Day=nday,Month=nmonth,Year=nyear,Hour=nhour,$
                        Minute=nminute,Second=nsecond

; dt_time is a date/time variable
; returns a date/time variable

jultime=dt_time.julian
mon=dt_time.month
dd=dt_time.day
yy=dt_time.year
hr=dt_time.hour
minu=dt_time.minute
sec=NINT(dt_time.second)

IF KEYWORD_SET(nday) THEN BEGIN
  newjultime=jultime+DOUBLE(nday)
ENDIF

IF KEYWORD_SET(nmonth) THEN BEGIN
  mon=mon+nmonth
  IF ((mon lt 1) or (mon gt 12)) THEN BEGIN
    IF (mon lt 1) THEN yy_to_add=(mon/12)-1
    IF (mon gt 12) THEN yy_to_add=((mon-1)/12)
    nmm=mon-(yy_to_add*12)
    temptime=var_to_dt(yy,nmm,dd,hr,minu,sec)
    newtime=dt_add(temptime,Year=yy_to_add)
    newjultime=newtime.julian
  ENDIF ELSE BEGIN
    newjultime=JULDAY(mon,dd,yy,hr,minu,sec)
  ENDELSE
ENDIF

IF KEYWORD_SET(nyear) THEN BEGIN
  nyy=yy+nyear
  newjultime=JULDAY(mon,dd,nyy,hr,minu,sec)
ENDIF

IF KEYWORD_SET(nhour) THEN BEGIN
  hr=hr+nhour
  IF ((hr lt 0) or (hr gt 23)) THEN BEGIN
    IF (hr lt 0) THEN dd_to_add=(hr/24)-1
    IF (hr gt 23) THEN dd_to_add=(hr/24)
    nhr=hr-(dd_to_add*24)
    temptime=var_to_dt(yy,mon,dd,nhr,minu,sec)
    newtime=dt_add(temptime,Day=dd_to_add)
    newjultime=newtime.julian
  ENDIF ELSE BEGIN
    newjultime=JULDAY(mon,dd,yy,hr,minu,sec)
  ENDELSE
ENDIF

IF KEYWORD_SET(nminute) THEN BEGIN
  minu=minu+nminute
  IF ((minu lt 0) or (minu gt 59)) THEN BEGIN
    IF (minu lt 0) THEN hr_to_add=(minu/60)-1
    IF (minu gt 59) THEN hr_to_add=(minu/60)
    nminu=minu-(hr_to_add*60)
    temptime=var_to_dt(yy,mon,dd,hr,nminu,sec)
    newtime=dt_add(temptime,Hour=hr_to_add)
    newjultime=newtime.julian
  ENDIF ELSE BEGIN
    newjultime=JULDAY(mon,dd,yy,hr,minu,sec)
  ENDELSE
ENDIF

IF KEYWORD_SET(nsecond) THEN BEGIN
  sec=sec+FIX(nsecond)
  IF ((sec lt 0) or (sec gt 59)) THEN BEGIN
    IF (sec lt 0) THEN minu_to_add=(sec/60)-1
    IF (sec gt 59) THEN minu_to_add=(sec/60)
    nsec=sec-(minu_to_add*60)
    temptime=var_to_dt(yy,mon,dd,hr,minu,nsec)
    newtime=dt_add(temptime,Minute=minu_to_add)
    newjultime=newtime.julian
  ENDIF ELSE BEGIN
    newjultime=JULDAY(mon,dd,yy,hr,minu,sec)
  ENDELSE
ENDIF

date=JUL_TO_DT(newjultime)

RETURN, date
END                        
