pro readmethead,filename,modtitle,fieldhead,namever=namever

;--------------------------------------------------------------
; procedure to read met datafile header
;
; DBR 09/04/2002
;--------------------------------------------------------------

if(n_elements(namever) eq 0)then namever=2

;-----------------------------------------------------------------------  
; read title text

readtitle,filename,modtitle,modhead1,modhead2,namever=namever

if (namever eq 3)then begin
  ; NAME III
  nfields=modhead2(where(modhead1 eq 'Number of field cols:'))
endif else begin
  ;NAME
  nfields=modhead2(where(modhead1 eq 'Number of fields:'))
endelse

nfields=fix(nfields(0))

;-----------------------------------------------------------------------  
; field headers

if (namever eq 3) then begin
  ; NAME III
  readdata,filename,fieldhead,nfields=nfields+5,nrows=17,nskip=19
endif else begin
  ; NAME
  readdata,filename,fieldhead,nfields=nfields+5,nrows=7,nskip=17
endelse

fieldhead=strtrim(fieldhead,2)
fieldhead=fieldhead(4:*,*)

end
