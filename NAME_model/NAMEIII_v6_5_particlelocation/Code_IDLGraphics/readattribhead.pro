pro readattribhead,filename,modtitle,fieldhead

;--------------------------------------------------------------
; procedure to read attribution header
; reads model header, time series headers and all series
;
; DBR 20/02/2002
;--------------------------------------------------------------


;-----------------------------------------------------------------------  
; read title text

readtitle,filename,modtitle,modhead1,modhead2

;-----------------------------------------------------------------------  
; extract grid details

nfields=modhead2(where(modhead1 eq 'Number of series:'))
nfields=fix(nfields(0))


;-----------------------------------------------------------------------  
; series headers

readdata,filename,fieldhead,nrows=6,nskip=17

fieldhead=strtrim(fieldhead,2)
fieldhead=fieldhead(5:*,*)

end
