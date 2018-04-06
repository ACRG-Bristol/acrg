FUNCTION NINT,float_number,Long=nlong

; Returns nearest integer to floating point number
; If Long keyword is present returns longword integer

IF KEYWORD_SET(nlong) THEN BEGIN
  int_number=LONG(float_number+0.5)
ENDIF ELSE BEGIN
  int_number=FIX(float_number+0.5)
ENDELSE


RETURN, int_number
END                        
