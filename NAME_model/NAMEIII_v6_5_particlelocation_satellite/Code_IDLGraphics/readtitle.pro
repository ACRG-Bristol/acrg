PRO readtitle, filename, modtitle, modhead1, modhead2, namever=namever, $
 nrecs=nrecs, nskip=nskip

;-------------------------------------------------------------------------------
; Procedure to read field title
;
; Lois Huggett 10/02/2009
;------------------------------------------------------------------------------- 
status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'readtitle: cannot open file:',filename
  return
endif

if N_ELEMENTS(nskip) EQ 0 THEN nskip=0
if N_ELEMENTS(namever) EQ 0 THEN namever=2
  
if (namever eq 3)then begin
  ; NAME III
  if N_ELEMENTS(nrecs) EQ 0 THEN nrecs=18
  restore, 'fieldheadtemplate.sav'
  ; fieldheadtemplate reads 2 strings per line, first string width 28, skips 0 rows automatically
  template=fieldheadtemplate 
endif else begin
  ;NAME
  if N_ELEMENTS(nrecs) EQ 0 THEN nrecs=17
  restore, 'fieldheadtemplateV2.sav'
  ; fieldheadtemplateV2 reads 2 strings per line, first string width 21, skips 0 rows automatically
  template=fieldheadtemplateV2
endelse

modhead=read_ascii(filename,num_records=nrecs,record_start=nskip,template=template,header=title)

modtitle=modhead.field1+modhead.field2
modtitle=strtrim(modtitle,2)

; NAME III
if namever eq 3 then modtitle=[title,modtitle]

modhead1=strtrim(modhead.field1,2)
modhead2=strtrim(modhead.field2,2)

END
