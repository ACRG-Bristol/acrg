pro timeaxis,range,ticklayout=ticklayout

; Sets up x-axis to display dates and times for ranges of between 1 second and 10 years 

; Lois Huggett
; 17/07/2009
; ----------------------------------------------------------------------------------------            

if (n_elements(ticklayout) eq 0) then ticklayout=2

; ----------------------------------------------------------------------------------------

!X.TICKLAYOUT = ticklayout

!X.TICKINTERVAL = 1     
if range gt 365*5+10 then begin
; 5 years +
  date_label=LABEL_DATE(DATE_FORMAT=['%Y'])
  !X.TICKFORMAT=['LABEL_DATE']
  !X.TICKUNITS=['Year']
endif else if range gt 95 then begin
; 3 months - 5 years
  date_label=LABEL_DATE(DATE_FORMAT=['%M','%Y'])
  !X.TICKFORMAT=['LABEL_DATE','LABEL_DATE'] 
  !X.TICKUNITS=['Month','Year']
  if range gt 240 then !X.TICKINTERVAL=2 
  if range gt 370 then !X.TICKINTERVAL=3
  if range gt 740 then !X.TICKINTERVAL=6
endif else if range gt 4 then begin
; 4 days - 3 months
  date_label=LABEL_DATE(DATE_FORMAT=['%D','%M %Y'])
  !X.TICKFORMAT=['LABEL_DATE','LABEL_DATE']
  !X.TICKUNITS=['Day','Month']
  if range gt 8 then !X.TICKINTERVAL=2
  if range gt 20 then !X.TICKINTERVAL=5
  if range gt 31 then !X.TICKINTERVAL=10
endif else if range gt 0.125 then begin
; 3 hours - 3 days
  date_label=LABEL_DATE(DATE_FORMAT=['%H:00','%D %M %Y'])
  !X.TICKFORMAT=['LABEL_DATE','LABEL_DATE']
  !X.TICKUNITS=['Hour','Day']
  if range gt 0.25 then !X.TICKINTERVAL=2
  if range gt 0.5 then !X.TICKINTERVAL=3
  if range gt 1 then !X.TICKINTERVAL=6
  if range gt 2 then !X.TICKINTERVAL=12
endif else if range gt 2.1E-3 then begin
; 3 minutes - 3 hours
  date_label=LABEL_DATE(DATE_FORMAT=['%H:%I','%D %M %Y'])
  !X.TICKFORMAT=['LABEL_DATE','LABEL_DATE']
  !X.TICKUNITS=['Time','Day']
  if range gt 3.5E-3 then !X.TICKINTERVAL=2
  if range gt 0.011 then !X.TICKINTERVAL=5
  if range gt 0.021 then !X.TICKINTERVAL=10
  if range gt 0.042 then !X.TICKINTERVAL=15
  if range gt 0.084 then !X.TICKINTERVAL=30
endif else begin
; less than 3 minutes
  date_label=LABEL_DATE(DATE_FORMAT=['%H:%I:%S','%D %M %Y'])
  !X.TICKFORMAT=['LABEL_DATE','LABEL_DATE']
  !X.TICKUNITS=['Time','Day']  
  if range lt 3.47E-5 then !X.TICKINTERVAL=0.5
  if range lt 2.31E-5 then !X.TICKINTERVAL=0.2
  if range gt 5.79E-4 then !X.TICKINTERVAL=2
  if range gt 1.74E-4 then !X.TICKINTERVAL=5
  if range gt 3.47E-4 then !X.TICKINTERVAL=10
  if range gt 6.94E-4 then !X.TICKINTERVAL=20
  if range gt 1.39E-3 then !X.TICKINTERVAL=30
endelse

End
