pro readfieldhead,filename,modtitle,fieldhead,xyrel,$
  grid,namever=namever,timeframe=timeframe


;--------------------------------------------------------------
; procedure to read field header
;
; DBR 20/02/2002
;--------------------------------------------------------------

; inputs:
;  filename   :
;  namever    : 
;  timeframe  : fixed or relative time

; outputs     
;  modtitle   : (string) array of title information from file
;  fieldhead  : (string) array of header information from file
;  xyrel      : source location
;  grid       : grid information

;--------------------------------------------------------------
if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'

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


;-----------------------------------------------------------------------  
; release location

cloc=strsplit(modhead2(where(modhead1 eq 'Release location:')),' ',/extract)
len=strlen(cloc)

dir=strarr(2)
dir(0)=strmid(cloc(0),len(0)-1)
dir(1)=strmid(cloc(1),len(1)-1)

if(dir(0) eq 'W' OR dir(0) eq 'E')then begin

  xyrel=dblarr(2)
  xyrel(0)=double(strmid(cloc(0),0,len(0)-1))
  xyrel(1)=double(strmid(cloc(1),0,len(1)-1))

  if(dir(0) eq 'W')then begin
    xyrel(0)=-xyrel(0)
  endif

  if(dir(1) eq 'S')then begin
    xyrel(1)=-xyrel(1)
  endif

endif else begin
  xyrel=[-999.0,-999.0]
endelse

;-----------------------------------------------------------------------  
; map grid

grid=dblarr(6)
grid(0)=float(modhead2(where(modhead1 eq 'X grid resolution:')))
grid(1)=float(modhead2(where(modhead1 eq 'Y grid resolution:')))
grid(2)=float(modhead2(where(modhead1 eq 'X grid size:')))
grid(3)=float(modhead2(where(modhead1 eq 'Y grid size:')))
grid(4)=float(modhead2(where(modhead1 eq 'X grid origin:')))
grid(5)=float(modhead2(where(modhead1 eq 'Y grid origin:')))

end

